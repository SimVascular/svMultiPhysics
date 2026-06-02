#include "LagrangeBasisPyramid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "Basis/BasisExceptions.h"
#include "BasisTolerance.h"
#include "Math/DenseLinearAlgebra.h"
#include "Math/DenseTransformKernels.h"
#include "LagrangeBasisUtility.h"
#include "PyramidModalBasis.h"

namespace svmp {
namespace FE {
namespace basis {
namespace detail {

class PyramidLagrangeCache {
public:
    using ModalTerm = pyramid_modal::Term;

    struct UvPolynomial {
        using Power = std::pair<int, int>;
        std::vector<std::pair<Power, Real>> coeffs;

        void add_term(int pu, int pv, Real coeff, Real tol = Real(1e-14)) {
            if (std::abs(coeff) <= tol) {
                return;
            }
            const auto key = std::make_pair(pu, pv);
            const auto found = std::lower_bound(
                coeffs.begin(),
                coeffs.end(),
                key,
                [](const auto& entry, const Power& value) { return entry.first < value; });
            if (found == coeffs.end() || found->first != key) {
                coeffs.insert(found, {key, coeff});
                return;
            }

            found->second += coeff;
            if (std::abs(found->second) <= tol) {
                coeffs.erase(found);
            }
        }

        void add_scaled(const UvPolynomial& other, Real scale, Real tol = Real(1e-14)) {
            if (std::abs(scale) <= tol) {
                return;
            }
            for (const auto& [powers, coeff] : other.coeffs) {
                add_term(powers.first, powers.second, scale * coeff, tol);
            }
        }

        bool empty(Real tol = Real(1e-12)) const {
            for (const auto& [powers, coeff] : coeffs) {
                (void)powers;
                if (std::abs(coeff) > tol) {
                    return false;
                }
            }
            return true;
        }

        bool is_constant(Real tol = Real(1e-12)) const {
            for (const auto& [powers, coeff] : coeffs) {
                if (std::abs(coeff) <= tol) {
                    continue;
                }
                if (powers.first != 0 || powers.second != 0) {
                    return false;
                }
            }
            return true;
        }

        Real constant_value(Real tol = Real(1e-12)) const {
            Real value = Real(0);
            for (const auto& [powers, coeff] : coeffs) {
                if (std::abs(coeff) <= tol) {
                    continue;
                }
                if (powers.first == 0 && powers.second == 0) {
                    value += coeff;
                }
            }
            return value;
        }
    };

    struct ApexSeries {
        std::vector<std::pair<int, UvPolynomial>> by_power;

        void add_term(int beta, int pu, int pv, Real coeff, Real tol = Real(1e-14)) {
            const auto found = find_or_insert(beta);
            found->second.add_term(pu, pv, coeff, tol);
            if (found->second.empty(tol)) {
                by_power.erase(found);
            }
        }

        void add_scaled(const ApexSeries& other, Real scale, Real tol = Real(1e-14)) {
            if (std::abs(scale) <= tol) {
                return;
            }
            for (const auto& [beta, poly] : other.by_power) {
                const auto found = find_or_insert(beta);
                found->second.add_scaled(poly, scale, tol);
                if (found->second.empty(tol)) {
                    by_power.erase(found);
                }
            }
        }

    private:
        std::vector<std::pair<int, UvPolynomial>>::iterator find_or_insert(int beta) {
            const auto found = std::lower_bound(
                by_power.begin(),
                by_power.end(),
                beta,
                [](const auto& entry, int value) { return entry.first < value; });
            if (found != by_power.end() && found->first == beta) {
                return found;
            }
            return by_power.insert(found, {beta, UvPolynomial{}});
        }
    };

    using GradientSeries = std::array<ApexSeries, 3>;
    using HessianSeries = std::array<std::array<ApexSeries, 3>, 3>;

    enum class ApexLimitKind {
        Constant,
        DirectionDependent,
        Singular,
    };

    enum class ApexRankStatus {
        Exact,
        DirectionDependent,
        Singular,
    };

    struct ApexClassification {
        ApexLimitKind kind{ApexLimitKind::Constant};
        Real constant_value{0};
        int leading_power{1};
    };

    struct ApexData {
        std::vector<Real> values;
        std::vector<Gradient> gradients;
        std::vector<Hessian> hessians;
        ApexRankStatus gradient_status{ApexRankStatus::Exact};
        ApexRankStatus hessian_status{ApexRankStatus::Exact};
    };

    struct OrderData {
        int order{0};
        std::vector<math::Vector<Real, 3>> nodes;
        std::vector<ModalTerm> modal_terms;
        std::vector<Real> modal_to_nodal;
        ApexData apex;
    };

    struct EvaluationScratch {
        std::vector<Real> modal_values;
        std::vector<Real> modal_gradient_components;
        std::vector<Real> modal_hessian_components;
        std::vector<Gradient> modal_gradients;
        std::vector<Hessian> modal_hessians;
        pyramid_modal::EvaluationPoint modal_point;

        void prewarm(std::size_t max_size, std::size_t max_qpts) {
            const std::size_t batched_size = max_size * std::max<std::size_t>(max_qpts, 1u);
            modal_values.reserve(batched_size);
            modal_gradient_components.reserve(batched_size * 3u);
            modal_hessian_components.reserve(batched_size * 9u);
            modal_gradients.reserve(max_size);
            modal_hessians.reserve(max_size);
        }
    };

    static EvaluationScratch& evaluation_scratch() {
        // Scratch is intentionally thread-local: production assembly uses a
        // persistent worker-thread team, so buffers stay warm on each worker.
        static thread_local EvaluationScratch scratch;
        return scratch;
    }

    static void prewarm_scratch(std::size_t max_size, std::size_t max_qpts) {
        evaluation_scratch().prewarm(max_size, max_qpts);
    }

    static bool is_apex_point(const math::Vector<Real, 3>& xi) {
        const Real tol = apex_coord_tolerance();
        return std::abs(xi[0]) <= tol &&
               std::abs(xi[1]) <= tol &&
               std::abs(Real(1) - xi[2]) <= tol;
    }

    static bool on_degenerate_top_plane(const math::Vector<Real, 3>& xi) {
        return basis_near_zero(Real(1) - xi[2]);
    }

    static void validate_top_plane_query(const math::Vector<Real, 3>& xi) {
        if (on_degenerate_top_plane(xi) && !is_apex_point(xi)) [[unlikely]] {
            throw BasisEvaluationException(
                "Pyramid reference evaluation on the degenerate z=1 plane is only defined at the apex",
                __FILE__, __LINE__, __func__);
        }
    }

    static OrderData build_order_data(int order) {
        OrderData data;
        data.order = order;

        data.nodes = build_public_nodes(order);
        data.modal_terms = pyramid_modal::build_terms(order);

        const std::size_t n = data.nodes.size();
        if (data.modal_terms.size() != n) {
            throw BasisConstructionException("LagrangeBasis pyramid modal basis size mismatch",
                                             __FILE__, __LINE__, __func__);
        }

        std::vector<Real> vandermonde(n * n, Real(0));
        for (std::size_t row = 0; row < n; ++row) {
            pyramid_modal::EvaluationPoint modal_point;
            pyramid_modal::prepare_evaluation_point(
                data.modal_terms, data.nodes[row], modal_point);
            for (std::size_t col = 0; col < n; ++col) {
                Real value = Real(0);
                pyramid_modal::evaluate_term(data.modal_terms[col], modal_point, value);
                vandermonde[row * n + col] = value;
            }
        }

        const auto inverse_result = math::invert_dense_matrix_with_diagnostics(
            std::move(vandermonde),
            n,
            "LagrangeBasis pyramid Vandermonde");
        math::validate_dense_inverse_diagnostics(
            inverse_result,
            n,
            "LagrangeBasis pyramid Vandermonde");
        const std::vector<Real>& inverse = inverse_result.inverse;

        data.modal_to_nodal.assign(n * n, Real(0));
        for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
            for (std::size_t modal_j = 0; modal_j < n; ++modal_j) {
                data.modal_to_nodal[basis_i * n + modal_j] =
                    inverse[modal_j * n + basis_i];
            }
        }
        data.apex = build_apex_data(data);
        return data;
    }

    static bool has_low_order_fast_modal_to_nodal(const OrderData& data) noexcept {
        return data.order == 1 || data.order == 2;
    }

    static const OrderData& get(int order) {
        constexpr int kMaxOnceCachedOrder = 12;
        if (order >= 0 && order <= kMaxOnceCachedOrder) {
            static std::array<std::once_flag, kMaxOnceCachedOrder + 1> flags;
            static std::array<std::unique_ptr<OrderData>, kMaxOnceCachedOrder + 1> cache;
            const auto idx = static_cast<std::size_t>(order);
            std::call_once(flags[idx], [idx, order]() {
                cache[idx] = std::make_unique<OrderData>(build_order_data(order));
            });
            return *cache[idx];
        }

        static std::mutex fallback_mutex;
        static std::map<int, std::unique_ptr<OrderData>> fallback_cache;

        std::lock_guard<std::mutex> lock(fallback_mutex);
        const auto found = fallback_cache.find(order);
        if (found != fallback_cache.end()) {
            return *found->second;
        }

        auto data = std::make_unique<OrderData>(build_order_data(order));
        const auto [it, inserted] = fallback_cache.emplace(order, std::move(data));
        (void)inserted;
        return *it->second;
    }

    static void evaluate_values(const OrderData& data,
                                const math::Vector<Real, 3>& xi,
                                std::vector<Real>& values) {
        validate_top_plane_query(xi);
        if (is_apex_point(xi)) {
            values = data.apex.values;
            return;
        }

        auto& scratch = evaluation_scratch();
        auto& modal = scratch.modal_values;
        auto& modal_point = scratch.modal_point;
        modal.resize(data.modal_terms.size());
        pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);
        for (std::size_t m = 0; m < data.modal_terms.size(); ++m) {
            pyramid_modal::evaluate_term(data.modal_terms[m], modal_point, modal[m]);
        }
        if (has_low_order_fast_modal_to_nodal(data)) {
            apply_sparse_basis_to_nodal(data, modal, values);
        } else {
            apply_modal_to_nodal(data, modal, values);
        }
    }

    static void evaluate_gradients(const OrderData& data,
                                   const math::Vector<Real, 3>& xi,
                                   std::vector<Gradient>& gradients) {
        validate_top_plane_query(xi);
        if (is_apex_point(xi)) {
            if (data.apex.gradient_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("gradient", data.apex.gradient_status),
                    __FILE__, __LINE__, __func__);
            }
            gradients = data.apex.gradients;
            return;
        }

        auto& scratch = evaluation_scratch();
        auto& modal_gradients = scratch.modal_gradients;
        auto& modal_point = scratch.modal_point;
        modal_gradients.resize(data.modal_terms.size());
        pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);
        for (std::size_t m = 0; m < data.modal_terms.size(); ++m) {
            Real value = Real(0);
            pyramid_modal::evaluate_term(data.modal_terms[m], modal_point, value, &modal_gradients[m]);
        }
        if (has_low_order_fast_modal_to_nodal(data)) {
            apply_sparse_basis_to_nodal(data, modal_gradients, gradients);
        } else {
            apply_modal_to_nodal(data, modal_gradients, gradients);
        }
    }

    static void evaluate_hessians(const OrderData& data,
                                  const math::Vector<Real, 3>& xi,
                                  std::vector<Hessian>& hessians) {
        validate_top_plane_query(xi);
        if (is_apex_point(xi)) {
            if (data.apex.hessian_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("Hessian", data.apex.hessian_status),
                    __FILE__, __LINE__, __func__);
            }
            hessians = data.apex.hessians;
            return;
        }

        auto& scratch = evaluation_scratch();
        auto& modal_hessians = scratch.modal_hessians;
        auto& modal_point = scratch.modal_point;
        modal_hessians.resize(data.modal_terms.size());
        pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);
        for (std::size_t m = 0; m < data.modal_terms.size(); ++m) {
            Real value = Real(0);
            pyramid_modal::evaluate_term(data.modal_terms[m], modal_point, value, nullptr, &modal_hessians[m]);
        }
        if (has_low_order_fast_modal_to_nodal(data)) {
            apply_sparse_basis_to_nodal(data, modal_hessians, hessians);
        } else {
            apply_modal_to_nodal(data, modal_hessians, hessians);
        }
    }

    static void evaluate_all(const OrderData& data,
                             const math::Vector<Real, 3>& xi,
                             std::vector<Real>& values,
                             std::vector<Gradient>& gradients,
                             std::vector<Hessian>& hessians) {
        validate_top_plane_query(xi);
        if (is_apex_point(xi)) {
            if (data.apex.gradient_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("gradient", data.apex.gradient_status),
                    __FILE__, __LINE__, __func__);
            }
            if (data.apex.hessian_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("Hessian", data.apex.hessian_status),
                    __FILE__, __LINE__, __func__);
            }
            values = data.apex.values;
            gradients = data.apex.gradients;
            hessians = data.apex.hessians;
            return;
        }

        const std::size_t n = data.modal_terms.size();
        auto& scratch = evaluation_scratch();
        auto& modal_values = scratch.modal_values;
        auto& modal_gradients = scratch.modal_gradients;
        auto& modal_hessians = scratch.modal_hessians;
        auto& modal_point = scratch.modal_point;
        modal_values.resize(n);
        modal_gradients.resize(n);
        modal_hessians.resize(n);
        pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);

        for (std::size_t m = 0; m < n; ++m) {
            pyramid_modal::evaluate_term(
                data.modal_terms[m], modal_point, modal_values[m], &modal_gradients[m], &modal_hessians[m]);
        }

        if (has_low_order_fast_modal_to_nodal(data)) {
            apply_sparse_basis_to_nodal_all(
                data, modal_values, modal_gradients, modal_hessians, values, gradients, hessians);
            return;
        }

        values.resize(n);
        gradients.resize(n);
        hessians.resize(n);
        for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
            const Real* row = data.modal_to_nodal.data() + basis_i * n;
            Gradient gradient{};
            Hessian hessian{};
            Real value = Real(0);
            for (std::size_t modal_j = 0; modal_j < n; ++modal_j) {
                const Real coeff = row[modal_j];
                value += coeff * modal_values[modal_j];

                const Real* modal_gradient = modal_gradients[modal_j].data();
                gradient[0] += coeff * modal_gradient[0];
                gradient[1] += coeff * modal_gradient[1];
                gradient[2] += coeff * modal_gradient[2];

                const Real* modal_hessian = modal_hessians[modal_j].data();
                Real* hessian_data = hessian.data();
                hessian_data[0] += coeff * modal_hessian[0];
                hessian_data[1] += coeff * modal_hessian[1];
                hessian_data[2] += coeff * modal_hessian[2];
                hessian_data[4] += coeff * modal_hessian[4];
                hessian_data[5] += coeff * modal_hessian[5];
                hessian_data[8] += coeff * modal_hessian[8];
            }
            values[basis_i] = value;
            gradients[basis_i] = gradient;
            Real* hessian_data = hessian.data();
            hessian_data[3] = hessian_data[1];
            hessian_data[6] = hessian_data[2];
            hessian_data[7] = hessian_data[5];
            hessians[basis_i] = hessian;
        }
    }

    static void evaluate_values_to(const OrderData& data,
                                   const math::Vector<Real, 3>& xi,
                                   Real* SVMP_RESTRICT values_out) {
        validate_top_plane_query(xi);
        if (is_apex_point(xi)) {
            std::copy(data.apex.values.begin(), data.apex.values.end(), values_out);
            return;
        }

        auto& scratch = evaluation_scratch();
        auto& modal = scratch.modal_values;
        auto& modal_point = scratch.modal_point;
        modal.resize(data.modal_terms.size());
        pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);
        for (std::size_t m = 0; m < data.modal_terms.size(); ++m) {
            pyramid_modal::evaluate_term(data.modal_terms[m], modal_point, modal[m]);
        }
        if (has_low_order_fast_modal_to_nodal(data)) {
            apply_sparse_basis_to_nodal_to(data, modal, values_out);
        } else {
            apply_modal_to_nodal_to(data, modal, values_out);
        }
    }

    static void evaluate_gradients_to(const OrderData& data,
                                      const math::Vector<Real, 3>& xi,
                                      Real* SVMP_RESTRICT gradients_out) {
        validate_top_plane_query(xi);
        if (is_apex_point(xi)) {
            if (data.apex.gradient_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("gradient", data.apex.gradient_status),
                    __FILE__, __LINE__, __func__);
            }
            for (std::size_t i = 0; i < data.apex.gradients.size(); ++i) {
                gradients_out[i * 3u + 0u] = data.apex.gradients[i][0];
                gradients_out[i * 3u + 1u] = data.apex.gradients[i][1];
                gradients_out[i * 3u + 2u] = data.apex.gradients[i][2];
            }
            return;
        }

        auto& scratch = evaluation_scratch();
        auto& modal_gradients = scratch.modal_gradients;
        auto& modal_point = scratch.modal_point;
        modal_gradients.resize(data.modal_terms.size());
        pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);
        for (std::size_t m = 0; m < data.modal_terms.size(); ++m) {
            Real value = Real(0);
            pyramid_modal::evaluate_term(data.modal_terms[m], modal_point, value, &modal_gradients[m]);
        }
        if (has_low_order_fast_modal_to_nodal(data)) {
            apply_sparse_basis_to_nodal_to(data, modal_gradients, gradients_out);
        } else {
            apply_modal_to_nodal_to(data, modal_gradients, gradients_out);
        }
    }

    static void evaluate_hessians_to(const OrderData& data,
                                     const math::Vector<Real, 3>& xi,
                                     Real* SVMP_RESTRICT hessians_out) {
        validate_top_plane_query(xi);
        if (is_apex_point(xi)) {
            if (data.apex.hessian_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("Hessian", data.apex.hessian_status),
                    __FILE__, __LINE__, __func__);
            }
            for (std::size_t i = 0; i < data.apex.hessians.size(); ++i) {
                store_hessian(data.apex.hessians[i], hessians_out + i * 9u);
            }
            return;
        }

        auto& scratch = evaluation_scratch();
        auto& modal_hessians = scratch.modal_hessians;
        auto& modal_point = scratch.modal_point;
        modal_hessians.resize(data.modal_terms.size());
        pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);
        for (std::size_t m = 0; m < data.modal_terms.size(); ++m) {
            Real value = Real(0);
            pyramid_modal::evaluate_term(data.modal_terms[m], modal_point, value, nullptr, &modal_hessians[m]);
        }
        if (has_low_order_fast_modal_to_nodal(data)) {
            apply_sparse_basis_to_nodal_to(data, modal_hessians, hessians_out);
        } else {
            apply_modal_to_nodal_to(data, modal_hessians, hessians_out);
        }
    }

    static void evaluate_all_to(const OrderData& data,
                                const math::Vector<Real, 3>& xi,
                                Real* SVMP_RESTRICT values_out,
                                Real* SVMP_RESTRICT gradients_out,
                                Real* SVMP_RESTRICT hessians_out) {
        validate_top_plane_query(xi);
        if (is_apex_point(xi)) {
            if (data.apex.gradient_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("gradient", data.apex.gradient_status),
                    __FILE__, __LINE__, __func__);
            }
            if (data.apex.hessian_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("Hessian", data.apex.hessian_status),
                    __FILE__, __LINE__, __func__);
            }
            std::copy(data.apex.values.begin(), data.apex.values.end(), values_out);
            for (std::size_t i = 0; i < data.apex.gradients.size(); ++i) {
                gradients_out[i * 3u + 0u] = data.apex.gradients[i][0];
                gradients_out[i * 3u + 1u] = data.apex.gradients[i][1];
                gradients_out[i * 3u + 2u] = data.apex.gradients[i][2];
            }
            for (std::size_t i = 0; i < data.apex.hessians.size(); ++i) {
                const Real* hessian = data.apex.hessians[i].data();
                std::copy(hessian, hessian + 9u, hessians_out + i * 9u);
            }
            return;
        }

        const std::size_t n = data.modal_terms.size();
        auto& scratch = evaluation_scratch();
        auto& modal_values = scratch.modal_values;
        auto& modal_gradients = scratch.modal_gradients;
        auto& modal_hessians = scratch.modal_hessians;
        auto& modal_point = scratch.modal_point;
        modal_values.resize(n);
        modal_gradients.resize(n);
        modal_hessians.resize(n);
        pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);

        for (std::size_t m = 0; m < n; ++m) {
            pyramid_modal::evaluate_term(
                data.modal_terms[m], modal_point, modal_values[m], &modal_gradients[m], &modal_hessians[m]);
        }

        if (has_low_order_fast_modal_to_nodal(data)) {
            apply_sparse_basis_to_nodal_all_to(
                data, modal_values, modal_gradients, modal_hessians, values_out, gradients_out, hessians_out);
            return;
        }

        for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
            const Real* row = data.modal_to_nodal.data() + basis_i * n;
            Real value = Real(0);
            Real gradient[3] = {Real(0), Real(0), Real(0)};
            Real hessian[9] = {};
            for (std::size_t modal_j = 0; modal_j < n; ++modal_j) {
                const Real coeff = row[modal_j];
                value += coeff * modal_values[modal_j];

                const Real* modal_gradient = modal_gradients[modal_j].data();
                gradient[0] += coeff * modal_gradient[0];
                gradient[1] += coeff * modal_gradient[1];
                gradient[2] += coeff * modal_gradient[2];

                const Real* modal_hessian = modal_hessians[modal_j].data();
                hessian[0] += coeff * modal_hessian[0];
                hessian[1] += coeff * modal_hessian[1];
                hessian[2] += coeff * modal_hessian[2];
                hessian[4] += coeff * modal_hessian[4];
                hessian[5] += coeff * modal_hessian[5];
                hessian[8] += coeff * modal_hessian[8];
            }

            values_out[basis_i] = value;
            Real* gradient_out = gradients_out + basis_i * 3u;
            gradient_out[0] = gradient[0];
            gradient_out[1] = gradient[1];
            gradient_out[2] = gradient[2];

            Real* hessian_out = hessians_out + basis_i * 9u;
            hessian_out[0] = hessian[0];
            hessian_out[1] = hessian[1];
            hessian_out[2] = hessian[2];
            hessian_out[3] = hessian[1];
            hessian_out[4] = hessian[4];
            hessian_out[5] = hessian[5];
            hessian_out[6] = hessian[2];
            hessian_out[7] = hessian[5];
            hessian_out[8] = hessian[8];
        }
    }

    static void evaluate_at_quadrature_points_strided(
        const OrderData& data,
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT gradients_out,
        Real* SVMP_RESTRICT hessians_out) {
        const unsigned mask = (values_out != nullptr ? 1u : 0u) |
                              (gradients_out != nullptr ? 2u : 0u) |
                              (hessians_out != nullptr ? 4u : 0u);
        switch (mask) {
            case 0u:
                validate_strided_points(points);
                return;
            case 1u:
                evaluate_at_quadrature_points_strided_impl<true, false, false>(
                    data, points, output_stride, values_out, gradients_out, hessians_out);
                return;
            case 2u:
                evaluate_at_quadrature_points_strided_impl<false, true, false>(
                    data, points, output_stride, values_out, gradients_out, hessians_out);
                return;
            case 3u:
                evaluate_at_quadrature_points_strided_impl<true, true, false>(
                    data, points, output_stride, values_out, gradients_out, hessians_out);
                return;
            case 4u:
                evaluate_at_quadrature_points_strided_impl<false, false, true>(
                    data, points, output_stride, values_out, gradients_out, hessians_out);
                return;
            case 5u:
                evaluate_at_quadrature_points_strided_impl<true, false, true>(
                    data, points, output_stride, values_out, gradients_out, hessians_out);
                return;
            case 6u:
                evaluate_at_quadrature_points_strided_impl<false, true, true>(
                    data, points, output_stride, values_out, gradients_out, hessians_out);
                return;
            case 7u:
                evaluate_at_quadrature_points_strided_impl<true, true, true>(
                    data, points, output_stride, values_out, gradients_out, hessians_out);
                return;
            default:
                return;
        }
    }

private:
    static void validate_strided_points(const std::vector<math::Vector<Real, 3>>& points) {
        for (const auto& xi : points) {
            validate_top_plane_query(xi);
        }
    }

    template <bool NeedValues, bool NeedGradients, bool NeedHessians>
    static void write_apex_strided(const OrderData& data,
                                   std::size_t q,
                                   std::size_t output_stride,
                                   Real* SVMP_RESTRICT values_out,
                                   Real* SVMP_RESTRICT gradients_out,
                                   Real* SVMP_RESTRICT hessians_out) {
        const std::size_t n = data.modal_terms.size();
        if constexpr (NeedValues) {
            for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
                values_out[basis_i * output_stride + q] = data.apex.values[basis_i];
            }
        }
        if constexpr (NeedGradients) {
            if (data.apex.gradient_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("gradient", data.apex.gradient_status),
                    __FILE__, __LINE__, __func__);
            }
            for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
                Real* g = gradients_out + basis_i * 3u * output_stride;
                g[0u * output_stride + q] = data.apex.gradients[basis_i][0];
                g[1u * output_stride + q] = data.apex.gradients[basis_i][1];
                g[2u * output_stride + q] = data.apex.gradients[basis_i][2];
            }
        }
        if constexpr (NeedHessians) {
            if (data.apex.hessian_status != ApexRankStatus::Exact) {
                throw BasisEvaluationException(
                    apex_status_message("Hessian", data.apex.hessian_status),
                    __FILE__, __LINE__, __func__);
            }
            for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
                const Real* hessian = data.apex.hessians[basis_i].data();
                Real* H = hessians_out + basis_i * 9u * output_stride;
                for (std::size_t component = 0; component < 9u; ++component) {
                    H[component * output_stride + q] = hessian[component];
                }
            }
        }
    }

    template <int Px,
              int Py,
              int Pz,
              int DenomPower,
              bool NeedValues,
              bool NeedGradients,
              bool NeedHessians>
    static void fill_low_order_modal_jet(std::size_t modal_i,
                                         const Real* SVMP_RESTRICT xp,
                                         const Real* SVMP_RESTRICT yp,
                                         const Real* SVMP_RESTRICT zp,
                                         const Real* SVMP_RESTRICT inv_tp,
                                         Real* SVMP_RESTRICT modal_values,
                                         Real (*SVMP_RESTRICT modal_gradients)[3],
                                         Real (*SVMP_RESTRICT modal_hessians)[9]) {
        const Real xy_base = xp[Px] * yp[Py];
        const Real base = xy_base * zp[Pz];
        const Real inv_denom = inv_tp[DenomPower];
        const Real value = base * inv_denom;

        if constexpr (NeedValues) {
            modal_values[modal_i] = value;
        }
        if constexpr (NeedGradients) {
            Real* g = modal_gradients[modal_i];
            if constexpr (Px > 0) {
                g[0] = static_cast<Real>(Px) * xp[Px - 1] * yp[Py] * zp[Pz] * inv_denom;
            } else {
                g[0] = Real(0);
            }
            if constexpr (Py > 0) {
                g[1] = static_cast<Real>(Py) * xp[Px] * yp[Py - 1] * zp[Pz] * inv_denom;
            } else {
                g[1] = Real(0);
            }
            Real gz = Real(0);
            if constexpr (Pz > 0) {
                gz += static_cast<Real>(Pz) * xy_base * zp[Pz - 1] * inv_denom;
            }
            if constexpr (DenomPower > 0) {
                gz += static_cast<Real>(DenomPower) * base * inv_tp[DenomPower + 1];
            }
            g[2] = gz;
        }
        if constexpr (NeedHessians) {
            Real* H = modal_hessians[modal_i];
            if constexpr (Px > 1) {
                H[0] = static_cast<Real>(Px * (Px - 1)) *
                       xp[Px - 2] * yp[Py] * zp[Pz] * inv_denom;
            } else {
                H[0] = Real(0);
            }
            if constexpr (Py > 1) {
                H[4] = static_cast<Real>(Py * (Py - 1)) *
                       xp[Px] * yp[Py - 2] * zp[Pz] * inv_denom;
            } else {
                H[4] = Real(0);
            }
            Real hxy = Real(0);
            if constexpr (Px > 0 && Py > 0) {
                hxy = static_cast<Real>(Px * Py) *
                      xp[Px - 1] * yp[Py - 1] * zp[Pz] * inv_denom;
            }
            H[1] = hxy;
            H[3] = hxy;

            Real hxz = Real(0);
            if constexpr (Px > 0) {
                constexpr Real px_real = static_cast<Real>(Px);
                const Real x_deriv_y = px_real * xp[Px - 1] * yp[Py];
                if constexpr (Pz > 0) {
                    hxz += x_deriv_y * static_cast<Real>(Pz) *
                           zp[Pz - 1] * inv_denom;
                }
                if constexpr (DenomPower > 0) {
                    hxz += x_deriv_y * static_cast<Real>(DenomPower) *
                           zp[Pz] * inv_tp[DenomPower + 1];
                }
            }
            H[2] = hxz;
            H[6] = hxz;

            Real hyz = Real(0);
            if constexpr (Py > 0) {
                constexpr Real py_real = static_cast<Real>(Py);
                const Real x_y_deriv = py_real * xp[Px] * yp[Py - 1];
                if constexpr (Pz > 0) {
                    hyz += x_y_deriv * static_cast<Real>(Pz) *
                           zp[Pz - 1] * inv_denom;
                }
                if constexpr (DenomPower > 0) {
                    hyz += x_y_deriv * static_cast<Real>(DenomPower) *
                           zp[Pz] * inv_tp[DenomPower + 1];
                }
            }
            H[5] = hyz;
            H[7] = hyz;

            Real hzz = Real(0);
            if constexpr (Pz > 1) {
                hzz += static_cast<Real>(Pz * (Pz - 1)) *
                       xy_base * zp[Pz - 2] * inv_denom;
            }
            if constexpr (Pz > 0 && DenomPower > 0) {
                hzz += static_cast<Real>(2 * Pz * DenomPower) * xy_base *
                       zp[Pz - 1] * inv_tp[DenomPower + 1];
            }
            if constexpr (DenomPower > 0) {
                hzz += static_cast<Real>(DenomPower * (DenomPower + 1)) *
                       base * inv_tp[DenomPower + 2];
            }
            H[8] = hzz;
        }
    }

    template <bool NeedValues, bool NeedGradients, bool NeedHessians>
    static void evaluate_low_order_modal_jets(const OrderData& data,
                                              const math::Vector<Real, 3>& xi,
                                              Real* SVMP_RESTRICT modal_values,
                                              Real (*SVMP_RESTRICT modal_gradients)[3],
                                              Real (*SVMP_RESTRICT modal_hessians)[9]) {
        const Real x = xi[0];
        const Real y = xi[1];
        const Real z = xi[2];
        const Real inv_t = Real(1) / (Real(1) - z);
        const Real xp[3] = {Real(1), x, x * x};
        const Real yp[3] = {Real(1), y, y * y};
        const Real zp[3] = {Real(1), z, z * z};
        Real inv_tp[5] = {Real(1), inv_t, Real(0), Real(0), Real(0)};
        inv_tp[2] = inv_tp[1] * inv_t;
        inv_tp[3] = inv_tp[2] * inv_t;
        inv_tp[4] = inv_tp[3] * inv_t;

        fill_low_order_modal_jet<0, 0, 0, 0, NeedValues, NeedGradients, NeedHessians>(
            0u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<1, 0, 0, 0, NeedValues, NeedGradients, NeedHessians>(
            1u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        if (data.order == 1) {
            fill_low_order_modal_jet<0, 1, 0, 0, NeedValues, NeedGradients, NeedHessians>(
                2u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
            fill_low_order_modal_jet<1, 1, 0, 1, NeedValues, NeedGradients, NeedHessians>(
                3u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
            fill_low_order_modal_jet<0, 0, 1, 0, NeedValues, NeedGradients, NeedHessians>(
                4u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
            return;
        }

        fill_low_order_modal_jet<2, 0, 0, 0, NeedValues, NeedGradients, NeedHessians>(
            2u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<0, 1, 0, 0, NeedValues, NeedGradients, NeedHessians>(
            3u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<1, 1, 0, 1, NeedValues, NeedGradients, NeedHessians>(
            4u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<2, 1, 0, 1, NeedValues, NeedGradients, NeedHessians>(
            5u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<0, 2, 0, 0, NeedValues, NeedGradients, NeedHessians>(
            6u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<1, 2, 0, 1, NeedValues, NeedGradients, NeedHessians>(
            7u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<2, 2, 0, 2, NeedValues, NeedGradients, NeedHessians>(
            8u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<0, 0, 1, 0, NeedValues, NeedGradients, NeedHessians>(
            9u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<1, 0, 1, 0, NeedValues, NeedGradients, NeedHessians>(
            10u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<0, 1, 1, 0, NeedValues, NeedGradients, NeedHessians>(
            11u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<1, 1, 1, 1, NeedValues, NeedGradients, NeedHessians>(
            12u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
        fill_low_order_modal_jet<0, 0, 2, 0, NeedValues, NeedGradients, NeedHessians>(
            13u, xp, yp, zp, inv_tp, modal_values, modal_gradients, modal_hessians);
    }

    template <bool NeedValues, bool NeedGradients, bool NeedHessians>
    static bool try_evaluate_low_order_strided(
        const OrderData& data,
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT gradients_out,
        Real* SVMP_RESTRICT hessians_out) {
        if (!has_low_order_fast_modal_to_nodal(data)) {
            return false;
        }
        for (const auto& xi : points) {
            validate_top_plane_query(xi);
            if (is_apex_point(xi)) {
                return false;
            }
        }

        Real modal_values[14];
        Real modal_gradients[14][3];
        Real modal_hessians[14][9];
        for (std::size_t q = 0; q < points.size(); ++q) {
            evaluate_low_order_modal_jets<NeedValues, NeedGradients, NeedHessians>(
                data, points[q], modal_values, modal_gradients, modal_hessians);
            if constexpr (NeedValues) {
                apply_low_order_combination(
                    data,
                    1u,
                    [&](std::size_t modal_i, std::size_t) {
                        return modal_values[modal_i];
                    },
                    [&](std::size_t basis_i, std::size_t, Real value) {
                        values_out[basis_i * output_stride + q] = value;
                    });
            }
            if constexpr (NeedGradients) {
                apply_low_order_combination(
                    data,
                    3u,
                    [&](std::size_t modal_i, std::size_t component) {
                        return modal_gradients[modal_i][component];
                    },
                    [&](std::size_t basis_i, std::size_t component, Real value) {
                        gradients_out[basis_i * 3u * output_stride +
                                      component * output_stride + q] = value;
                    });
            }
            if constexpr (NeedHessians) {
                apply_low_order_combination(
                    data,
                    9u,
                    [&](std::size_t modal_i, std::size_t component) {
                        return modal_hessians[modal_i][component];
                    },
                    [&](std::size_t basis_i, std::size_t component, Real value) {
                        hessians_out[basis_i * 9u * output_stride +
                                     component * output_stride + q] = value;
                    });
            }
        }
        return true;
    }

    template <bool NeedValues, bool NeedGradients, bool NeedHessians>
    static void evaluate_at_quadrature_points_strided_impl(
        const OrderData& data,
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT gradients_out,
        Real* SVMP_RESTRICT hessians_out) {
        const std::size_t n = data.modal_terms.size();
        if (points.empty() || n == 0u) {
            return;
        }
        if (try_evaluate_low_order_strided<NeedValues, NeedGradients, NeedHessians>(
                data, points, output_stride, values_out, gradients_out, hessians_out)) {
            return;
        }

        auto& scratch = evaluation_scratch();
        auto& modal_values = scratch.modal_values;
        auto& modal_gradients = scratch.modal_gradients;
        auto& modal_hessians = scratch.modal_hessians;
        auto& modal_point = scratch.modal_point;
        if constexpr (NeedValues) {
            modal_values.resize(n);
        }
        if constexpr (NeedGradients) {
            modal_gradients.resize(n);
        }
        if constexpr (NeedHessians) {
            modal_hessians.resize(n);
        }
        const bool use_fast_modal_to_nodal = has_low_order_fast_modal_to_nodal(data);

        if (!use_fast_modal_to_nodal) {
            bool has_apex_query = false;
            for (const auto& xi : points) {
                validate_top_plane_query(xi);
                has_apex_query = has_apex_query || is_apex_point(xi);
            }

            if (!has_apex_query) {
                const std::size_t num_qpts = points.size();
                if constexpr (NeedValues) {
                    modal_values.resize(n * num_qpts);
                }
                if constexpr (NeedGradients) {
                    scratch.modal_gradient_components.resize(n * 3u * num_qpts);
                }
                if constexpr (NeedHessians) {
                    scratch.modal_hessian_components.resize(n * 9u * num_qpts);
                }

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const auto& xi = points[q];
                    pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);
                    for (std::size_t modal_j = 0; modal_j < n; ++modal_j) {
                        Real modal_value = Real(0);
                        Gradient modal_gradient{};
                        Hessian modal_hessian{};
                        pyramid_modal::evaluate_term(
                            data.modal_terms[modal_j],
                            modal_point,
                            modal_value,
                            NeedGradients ? &modal_gradient : nullptr,
                            NeedHessians ? &modal_hessian : nullptr);
                        if constexpr (NeedValues) {
                            modal_values[modal_j * num_qpts + q] = modal_value;
                        }
                        if constexpr (NeedGradients) {
                            for (std::size_t component = 0; component < 3u; ++component) {
                                scratch.modal_gradient_components[
                                    (modal_j * 3u + component) * num_qpts + q] =
                                    modal_gradient[component];
                            }
                        }
                        if constexpr (NeedHessians) {
                            for (std::size_t component = 0; component < 9u; ++component) {
                                scratch.modal_hessian_components[
                                    (modal_j * 9u + component) * num_qpts + q] =
                                    modal_hessian.data()[component];
                            }
                        }
                    }
                }

                const Real* transform = data.modal_to_nodal.data();
                if constexpr (NeedValues) {
                    math::dense_transform_batched_row_major(
                        transform,
                        n,
                        n,
                        modal_values.data(),
                        num_qpts,
                        values_out,
                        output_stride,
                        num_qpts);
                }
                if constexpr (NeedGradients) {
                    for (std::size_t component = 0; component < 3u; ++component) {
                        math::dense_transform_batched_row_major(
                            transform,
                            n,
                            n,
                            scratch.modal_gradient_components.data() + component * num_qpts,
                            3u * num_qpts,
                            gradients_out + component * output_stride,
                            3u * output_stride,
                            num_qpts);
                    }
                }
                if constexpr (NeedHessians) {
                    for (std::size_t component = 0; component < 9u; ++component) {
                        math::dense_transform_batched_row_major(
                            transform,
                            n,
                            n,
                            scratch.modal_hessian_components.data() + component * num_qpts,
                            9u * num_qpts,
                            hessians_out + component * output_stride,
                            9u * output_stride,
                            num_qpts);
                    }
                }
                return;
            }
        }

        for (std::size_t q = 0; q < points.size(); ++q) {
            const auto& xi = points[q];
            validate_top_plane_query(xi);

            if (is_apex_point(xi)) {
                write_apex_strided<NeedValues, NeedGradients, NeedHessians>(
                    data, q, output_stride, values_out, gradients_out, hessians_out);
                continue;
            }

            pyramid_modal::prepare_evaluation_point(data.modal_terms, xi, modal_point);
            for (std::size_t modal_j = 0; modal_j < n; ++modal_j) {
                Gradient* gradient_out = nullptr;
                Hessian* hessian_out = nullptr;
                if constexpr (NeedGradients) {
                    gradient_out = &modal_gradients[modal_j];
                }
                if constexpr (NeedHessians) {
                    hessian_out = &modal_hessians[modal_j];
                }
                if constexpr (NeedValues) {
                    pyramid_modal::evaluate_term(
                        data.modal_terms[modal_j],
                        modal_point,
                        modal_values[modal_j],
                        gradient_out,
                        hessian_out);
                } else {
                    Real value = Real(0);
                    pyramid_modal::evaluate_term(
                        data.modal_terms[modal_j],
                        modal_point,
                        value,
                        gradient_out,
                        hessian_out);
                }
            }

            if (use_fast_modal_to_nodal) {
                if constexpr (NeedValues) {
                    apply_low_order_combination(
                        data,
                        1u,
                        [&](std::size_t modal_i, std::size_t) {
                            return modal_values[modal_i];
                        },
                        [&](std::size_t basis_i, std::size_t, Real value) {
                            values_out[basis_i * output_stride + q] = value;
                        });
                }
                if constexpr (NeedGradients) {
                    apply_low_order_combination(
                        data,
                        3u,
                        [&](std::size_t modal_i, std::size_t component) {
                            return modal_gradients[modal_i][component];
                        },
                        [&](std::size_t basis_i, std::size_t component, Real value) {
                            gradients_out[basis_i * 3u * output_stride +
                                          component * output_stride + q] = value;
                        });
                }
                if constexpr (NeedHessians) {
                    apply_low_order_combination(
                        data,
                        9u,
                        [&](std::size_t modal_i, std::size_t component) {
                            return modal_hessians[modal_i].data()[component];
                        },
                        [&](std::size_t basis_i, std::size_t component, Real value) {
                            hessians_out[basis_i * 9u * output_stride +
                                         component * output_stride + q] = value;
                        });
                }
                continue;
            }

            for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
                const Real* matrix_row = data.modal_to_nodal.data() + basis_i * n;
                [[maybe_unused]] Real value = Real(0);
                [[maybe_unused]] std::array<Real, 3> gradient{};
                [[maybe_unused]] std::array<Real, 9> hessian{};

                for (std::size_t modal_j = 0; modal_j < n; ++modal_j) {
                    const Real coeff = matrix_row[modal_j];
                    if constexpr (NeedValues) {
                        value += coeff * modal_values[modal_j];
                    }
                    if constexpr (NeedGradients) {
                        const Real* modal_gradient = modal_gradients[modal_j].data();
                        gradient[0] += coeff * modal_gradient[0];
                        gradient[1] += coeff * modal_gradient[1];
                        gradient[2] += coeff * modal_gradient[2];
                    }
                    if constexpr (NeedHessians) {
                        const Real* modal_hessian = modal_hessians[modal_j].data();
                        for (std::size_t component = 0; component < 9u; ++component) {
                            hessian[component] += coeff * modal_hessian[component];
                        }
                    }
                }

                if constexpr (NeedValues) {
                    values_out[basis_i * output_stride + q] = value;
                }
                if constexpr (NeedGradients) {
                    Real* g = gradients_out + basis_i * 3u * output_stride;
                    g[0u * output_stride + q] = gradient[0];
                    g[1u * output_stride + q] = gradient[1];
                    g[2u * output_stride + q] = gradient[2];
                }
                if constexpr (NeedHessians) {
                    Real* H = hessians_out + basis_i * 9u * output_stride;
                    for (std::size_t component = 0; component < 9u; ++component) {
                        H[component * output_stride + q] = hessian[component];
                    }
                }
            }
        }
    }

    static Real apex_coord_tolerance() noexcept {
        return basis_scaled_tolerance();
    }

    // Coefficient pruning for symbolic apex series, not a reference-coordinate
    // roundoff test. Keep this strict and separate from BasisTolerance.
    static constexpr Real kSeriesTolerance = Real(1e-12);

    static Real binomial_coeff(int n, int k) {
        if (k < 0 || k > n) {
            return Real(0);
        }
        if (k == 0 || k == n) {
            return Real(1);
        }
        k = std::min(k, n - k);
        Real coeff = Real(1);
        for (int i = 1; i <= k; ++i) {
            coeff *= static_cast<Real>(n - (k - i));
            coeff /= static_cast<Real>(i);
        }
        return coeff;
    }

    static void add_z_expansion(ApexSeries& series,
                                int z_power,
                                int beta0,
                                int pu,
                                int pv,
                                Real coeff) {
        for (int q = 0; q <= z_power; ++q) {
            const Real z_coeff = coeff * binomial_coeff(z_power, q) *
                                 ((q % 2 == 0) ? Real(1) : Real(-1));
            series.add_term(beta0 + q, pu, pv, z_coeff, kSeriesTolerance);
        }
    }

    static ApexSeries modal_value_asymptotic(const ModalTerm& term) {
        ApexSeries series;
        add_z_expansion(series,
                        term.pz,
                        term.px + term.py - term.denom_power,
                        term.px,
                        term.py,
                        Real(1));
        return series;
    }

    static GradientSeries modal_gradient_asymptotic(const ModalTerm& term) {
        GradientSeries gradient_series{};

        if (term.px > 0) {
            add_z_expansion(gradient_series[0],
                            term.pz,
                            term.px - 1 + term.py - term.denom_power,
                            term.px - 1,
                            term.py,
                            static_cast<Real>(term.px));
        }

        if (term.py > 0) {
            add_z_expansion(gradient_series[1],
                            term.pz,
                            term.px + term.py - 1 - term.denom_power,
                            term.px,
                            term.py - 1,
                            static_cast<Real>(term.py));
        }

        if (term.pz > 0) {
            add_z_expansion(gradient_series[2],
                            term.pz - 1,
                            term.px + term.py - term.denom_power,
                            term.px,
                            term.py,
                            static_cast<Real>(term.pz));
        }
        if (term.denom_power > 0) {
            add_z_expansion(gradient_series[2],
                            term.pz,
                            term.px + term.py - term.denom_power - 1,
                            term.px,
                            term.py,
                            static_cast<Real>(term.denom_power));
        }

        return gradient_series;
    }

    static HessianSeries modal_hessian_asymptotic(const ModalTerm& term) {
        HessianSeries hessian_series{};

        if (term.px > 1) {
            add_z_expansion(hessian_series[0][0],
                            term.pz,
                            term.px - 2 + term.py - term.denom_power,
                            term.px - 2,
                            term.py,
                            static_cast<Real>(term.px * (term.px - 1)));
        }

        if (term.py > 1) {
            add_z_expansion(hessian_series[1][1],
                            term.pz,
                            term.px + term.py - 2 - term.denom_power,
                            term.px,
                            term.py - 2,
                            static_cast<Real>(term.py * (term.py - 1)));
        }

        if (term.px > 0 && term.py > 0) {
            add_z_expansion(hessian_series[0][1],
                            term.pz,
                            term.px + term.py - 2 - term.denom_power,
                            term.px - 1,
                            term.py - 1,
                            static_cast<Real>(term.px * term.py));
            hessian_series[1][0] = hessian_series[0][1];
        }

        if (term.px > 0 && term.pz > 0) {
            add_z_expansion(hessian_series[0][2],
                            term.pz - 1,
                            term.px - 1 + term.py - term.denom_power,
                            term.px - 1,
                            term.py,
                            static_cast<Real>(term.px * term.pz));
        }
        if (term.px > 0 && term.denom_power > 0) {
            add_z_expansion(hessian_series[0][2],
                            term.pz,
                            term.px - 1 + term.py - term.denom_power - 1,
                            term.px - 1,
                            term.py,
                            static_cast<Real>(term.px * term.denom_power));
        }
        hessian_series[2][0] = hessian_series[0][2];

        if (term.py > 0 && term.pz > 0) {
            add_z_expansion(hessian_series[1][2],
                            term.pz - 1,
                            term.px + term.py - 1 - term.denom_power,
                            term.px,
                            term.py - 1,
                            static_cast<Real>(term.py * term.pz));
        }
        if (term.py > 0 && term.denom_power > 0) {
            add_z_expansion(hessian_series[1][2],
                            term.pz,
                            term.px + term.py - 1 - term.denom_power - 1,
                            term.px,
                            term.py - 1,
                            static_cast<Real>(term.py * term.denom_power));
        }
        hessian_series[2][1] = hessian_series[1][2];

        if (term.pz > 1) {
            add_z_expansion(hessian_series[2][2],
                            term.pz - 2,
                            term.px + term.py - term.denom_power,
                            term.px,
                            term.py,
                            static_cast<Real>(term.pz * (term.pz - 1)));
        }
        if (term.pz > 0 && term.denom_power > 0) {
            add_z_expansion(hessian_series[2][2],
                            term.pz - 1,
                            term.px + term.py - term.denom_power - 1,
                            term.px,
                            term.py,
                            static_cast<Real>(2 * term.pz * term.denom_power));
        }
        if (term.denom_power > 0) {
            add_z_expansion(hessian_series[2][2],
                            term.pz,
                            term.px + term.py - term.denom_power - 2,
                            term.px,
                            term.py,
                            static_cast<Real>(term.denom_power * (term.denom_power + 1)));
        }

        return hessian_series;
    }

    static ApexClassification classify_series(const ApexSeries& series) {
        for (const auto& [beta, poly] : series.by_power) {
            if (poly.empty(kSeriesTolerance)) {
                continue;
            }
            if (beta < 0) {
                return {ApexLimitKind::Singular, Real(0), beta};
            }
            if (beta > 0) {
                return {ApexLimitKind::Constant, Real(0), beta};
            }
            if (poly.is_constant(kSeriesTolerance)) {
                return {ApexLimitKind::Constant, poly.constant_value(kSeriesTolerance), beta};
            }
            return {ApexLimitKind::DirectionDependent, Real(0), beta};
        }
        return {ApexLimitKind::Constant, Real(0), 1};
    }

    static void accumulate_rank_status(ApexRankStatus& status,
                                       const ApexClassification& classification) {
        if (classification.kind == ApexLimitKind::Singular) {
            status = ApexRankStatus::Singular;
            return;
        }
        if (classification.kind == ApexLimitKind::DirectionDependent &&
            status != ApexRankStatus::Singular) {
            status = ApexRankStatus::DirectionDependent;
        }
    }

    static std::string apex_status_message(const char* rank,
                                           ApexRankStatus status) {
        switch (status) {
            case ApexRankStatus::DirectionDependent:
                return std::string("Pyramid rational nodal ") + rank +
                       " at the exact apex is not uniquely defined under admissible interior approaches";
            case ApexRankStatus::Singular:
                return std::string("Pyramid rational nodal ") + rank +
                       " at the exact apex is singular for this basis family";
            case ApexRankStatus::Exact:
                return std::string("Pyramid rational nodal ") + rank +
                       " apex evaluation unexpectedly reported non-exact status";
        }
        return std::string("Pyramid rational nodal ") + rank +
               " apex evaluation is not available";
    }

    static ApexData build_apex_data(const OrderData& data) {
        const std::size_t n = data.modal_terms.size();

        std::vector<ApexSeries> modal_values(n);
        std::vector<GradientSeries> modal_gradients(n);
        std::vector<HessianSeries> modal_hessians(n);
        for (std::size_t m = 0; m < n; ++m) {
            modal_values[m] = modal_value_asymptotic(data.modal_terms[m]);
            modal_gradients[m] = modal_gradient_asymptotic(data.modal_terms[m]);
            modal_hessians[m] = modal_hessian_asymptotic(data.modal_terms[m]);
        }

        std::vector<ApexSeries> nodal_values(n);
        std::vector<GradientSeries> nodal_gradients(n);
        std::vector<HessianSeries> nodal_hessians(n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t m = 0; m < n; ++m) {
                const Real coeff = data.modal_to_nodal[i * n + m];
                nodal_values[i].add_scaled(modal_values[m], coeff, kSeriesTolerance);
                for (int d = 0; d < 3; ++d) {
                    nodal_gradients[i][static_cast<std::size_t>(d)].add_scaled(
                        modal_gradients[m][static_cast<std::size_t>(d)], coeff, kSeriesTolerance);
                }
                for (int r = 0; r < 3; ++r) {
                    for (int c = 0; c < 3; ++c) {
                        nodal_hessians[i][static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]
                            .add_scaled(
                                modal_hessians[m][static_cast<std::size_t>(r)][static_cast<std::size_t>(c)],
                                coeff,
                                kSeriesTolerance);
                    }
                }
            }
        }

        ApexData apex;
        apex.values.assign(n, Real(0));
        apex.gradients.assign(n, Gradient{});
        apex.hessians.assign(n, Hessian{});

        for (std::size_t i = 0; i < n; ++i) {
            const ApexClassification value_class = classify_series(nodal_values[i]);
            if (value_class.kind != ApexLimitKind::Constant) {
                throw BasisConstructionException(
                    "Pyramid nodal value at apex is not uniquely defined for basis index " +
                    std::to_string(i),
                    __FILE__, __LINE__, __func__);
            }
            apex.values[i] = value_class.constant_value;

            for (int d = 0; d < 3; ++d) {
                const ApexClassification grad_class = classify_series(
                    nodal_gradients[i][static_cast<std::size_t>(d)]);
                accumulate_rank_status(apex.gradient_status, grad_class);
                if (grad_class.kind == ApexLimitKind::Constant) {
                    apex.gradients[i][static_cast<std::size_t>(d)] = grad_class.constant_value;
                }
            }

            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    const ApexClassification hess_class = classify_series(
                        nodal_hessians[i][static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]);
                    accumulate_rank_status(apex.hessian_status, hess_class);
                    if (hess_class.kind == ApexLimitKind::Constant) {
                        apex.hessians[i](static_cast<std::size_t>(r),
                                         static_cast<std::size_t>(c)) = hess_class.constant_value;
                    }
                }
            }
        }

        if (apex.gradient_status != ApexRankStatus::Exact) {
            apex.gradients.clear();
        }
        if (apex.hessian_status != ApexRankStatus::Exact) {
            apex.hessians.clear();
        }

        return apex;
    }

    static std::vector<math::Vector<Real, 3>> build_public_nodes(int order) {
        if (order == 0) {
            return {math::Vector<Real, 3>{Real(0), Real(0), Real(0.25)}};
        }

        std::vector<math::Vector<Real, 3>> nodes;
        nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 2) * (2 * order + 3) / 6));

        nodes.push_back(math::Vector<Real, 3>{Real(-1), Real(-1), Real(0)});
        nodes.push_back(math::Vector<Real, 3>{Real(1), Real(-1), Real(0)});
        nodes.push_back(math::Vector<Real, 3>{Real(1), Real(1), Real(0)});
        nodes.push_back(math::Vector<Real, 3>{Real(-1), Real(1), Real(0)});
        nodes.push_back(math::Vector<Real, 3>{Real(0), Real(0), Real(1)});

        for (int m = 1; m < order; ++m) {
            nodes.push_back(math::Vector<Real, 3>{equispaced_pm_one_coord(m, order), Real(-1), Real(0)});
        }
        for (int m = 1; m < order; ++m) {
            nodes.push_back(math::Vector<Real, 3>{Real(1), equispaced_pm_one_coord(m, order), Real(0)});
        }
        for (int m = order - 1; m >= 1; --m) {
            nodes.push_back(math::Vector<Real, 3>{equispaced_pm_one_coord(m, order), Real(1), Real(0)});
        }
        for (int m = order - 1; m >= 1; --m) {
            nodes.push_back(math::Vector<Real, 3>{Real(-1), equispaced_pm_one_coord(m, order), Real(0)});
        }

        for (int level = 1; level < order; ++level) {
            const Real z = static_cast<Real>(level) / static_cast<Real>(order);
            const Real scale = Real(1) - z;
            nodes.push_back(math::Vector<Real, 3>{-scale, -scale, z});
            nodes.push_back(math::Vector<Real, 3>{scale, -scale, z});
            nodes.push_back(math::Vector<Real, 3>{scale, scale, z});
            nodes.push_back(math::Vector<Real, 3>{-scale, scale, z});
        }

        for (int j = 1; j < order; ++j) {
            for (int i = 1; i < order; ++i) {
                nodes.push_back(math::Vector<Real, 3>{equispaced_pm_one_coord(i, order),
                                                      equispaced_pm_one_coord(j, order),
                                                      Real(0)});
            }
        }

        for (int level = 1; level < order - 1; ++level) {
            const int n = order - level;
            const Real z = static_cast<Real>(level) / static_cast<Real>(order);
            const Real scale = Real(1) - z;

            for (int m = 1; m < n; ++m) {
                const Real s = equispaced_pm_one_coord(m, n) * scale;
                nodes.push_back(math::Vector<Real, 3>{s, -scale, z});
            }
            for (int m = 1; m < n; ++m) {
                const Real s = equispaced_pm_one_coord(m, n) * scale;
                nodes.push_back(math::Vector<Real, 3>{scale, s, z});
            }
            for (int m = n - 1; m >= 1; --m) {
                const Real s = equispaced_pm_one_coord(m, n) * scale;
                nodes.push_back(math::Vector<Real, 3>{s, scale, z});
            }
            for (int m = n - 1; m >= 1; --m) {
                const Real s = equispaced_pm_one_coord(m, n) * scale;
                nodes.push_back(math::Vector<Real, 3>{-scale, s, z});
            }
        }

        for (int level = 1; level < order - 1; ++level) {
            const int n = order - level;
            const Real z = static_cast<Real>(level) / static_cast<Real>(order);
            const Real scale = Real(1) - z;
            for (int j = 1; j < n; ++j) {
                for (int i = 1; i < n; ++i) {
                    nodes.push_back(math::Vector<Real, 3>{equispaced_pm_one_coord(i, n) * scale,
                                                          equispaced_pm_one_coord(j, n) * scale,
                                                          z});
                }
            }
        }

        return nodes;
    }

    struct VectorValueSink {
        std::vector<Real>& output;
        void resize(std::size_t n) const { output.resize(n); }
        void write(std::size_t i, Real value) const { output[i] = value; }
    };

    struct RawValueSink {
        Real* output;
        void resize(std::size_t) const {}
        void write(std::size_t i, Real value) const { output[i] = value; }
    };

    struct VectorGradientSink {
        std::vector<Gradient>& output;
        void resize(std::size_t n) const { output.resize(n); }
        void write(std::size_t i, const Gradient& value) const { output[i] = value; }
    };

    struct RawGradientSink {
        Real* output;
        void resize(std::size_t) const {}
        void write(std::size_t i, const Gradient& value) const {
            Real* dst = output + i * 3u;
            dst[0] = value[0];
            dst[1] = value[1];
            dst[2] = value[2];
        }
    };

    struct VectorHessianSink {
        std::vector<Hessian>& output;
        void resize(std::size_t n) const { output.resize(n); }
        void write(std::size_t i, const Hessian& value) const { output[i] = value; }
    };

    struct RawHessianSink {
        Real* output;
        void resize(std::size_t) const {}
        void write(std::size_t i, const Hessian& value) const {
            store_hessian(value, output + i * 9u);
        }
    };

    template <typename Get, typename Set>
    static void apply_order1_combination(std::size_t components,
                                         const Get& get,
                                         const Set& set) {
        for (std::size_t c = 0; c < components; ++c) {
            const Real m0 = get(0u, c);
            const Real m1 = get(1u, c);
            const Real m2 = get(2u, c);
            const Real m3 = get(3u, c);
            const Real m4 = get(4u, c);
            set(0u, c, Real(0.25) * (m0 - m1 - m2 + m3 - m4));
            set(1u, c, Real(0.25) * (m0 + m1 - m2 - m3 - m4));
            set(2u, c, Real(0.25) * (m0 + m1 + m2 + m3 - m4));
            set(3u, c, Real(0.25) * (m0 - m1 + m2 - m3 - m4));
            set(4u, c, m4);
        }
    }

    template <typename Get, typename Set>
    static void apply_order2_combination(std::size_t components,
                                         const Get& get,
                                         const Set& set) {
        for (std::size_t c = 0; c < components; ++c) {
            const Real m0 = get(0u, c);
            const Real m1 = get(1u, c);
            const Real m2 = get(2u, c);
            const Real m3 = get(3u, c);
            const Real m4 = get(4u, c);
            const Real m5 = get(5u, c);
            const Real m6 = get(6u, c);
            const Real m7 = get(7u, c);
            const Real m8 = get(8u, c);
            const Real m9 = get(9u, c);
            const Real m10 = get(10u, c);
            const Real m11 = get(11u, c);
            const Real m12 = get(12u, c);
            const Real m13 = get(13u, c);
            set(0u, c, Real(0.25) * (m4 - m5 - m7 + m8 - m9 + m10 + m11 - Real(2) * m12 + m13));
            set(1u, c, Real(0.25) * (-m4 - m5 + m7 + m8 - m9 - m10 + m11 + Real(2) * m12 + m13));
            set(2u, c, Real(0.25) * (m4 + m5 + m7 + m8 - m9 - m10 - m11 - Real(2) * m12 + m13));
            set(3u, c, Real(0.25) * (-m4 + m5 - m7 + m8 - m9 + m10 - m11 + Real(2) * m12 + m13));
            set(4u, c, -m9 + Real(2) * m13);
            set(5u, c, Real(0.5) * (-m3 + m5 + m6 - m8 + m11));
            set(6u, c, Real(0.5) * (m1 + m2 - m7 - m8 - m10));
            set(7u, c, Real(0.5) * (m3 - m5 + m6 - m8 - m11));
            set(8u, c, Real(0.5) * (-m1 + m2 + m7 - m8 + m10));
            set(9u, c, m9 - m10 - m11 + m12 - m13);
            set(10u, c, m9 + m10 - m11 - m12 - m13);
            set(11u, c, m9 + m10 + m11 + m12 - m13);
            set(12u, c, m9 - m10 + m11 - m12 - m13);
            set(13u, c, m0 - m2 - m6 + m8 - Real(2) * m9 + m13);
        }
    }

    template <typename Get, typename Set>
    static void apply_low_order_combination(const OrderData& data,
                                            std::size_t components,
                                            const Get& get,
                                            const Set& set) {
        if (data.order == 1) {
            apply_order1_combination(components, get, set);
            return;
        }
        apply_order2_combination(components, get, set);
    }

    static void apply_sparse_basis_to_nodal(const OrderData& data,
                                            const std::vector<Real>& modal_values,
                                            std::vector<Real>& nodal_values) {
        const std::size_t n = modal_values.size();
        nodal_values.resize(n);
        apply_low_order_combination(
            data,
            1u,
            [&](std::size_t modal_i, std::size_t) { return modal_values[modal_i]; },
            [&](std::size_t basis_i, std::size_t, Real value) { nodal_values[basis_i] = value; });
    }

    static void apply_sparse_basis_to_nodal_to(const OrderData& data,
                                               const std::vector<Real>& modal_values,
                                               Real* SVMP_RESTRICT nodal_values) {
        apply_low_order_combination(
            data,
            1u,
            [&](std::size_t modal_i, std::size_t) { return modal_values[modal_i]; },
            [&](std::size_t basis_i, std::size_t, Real value) { nodal_values[basis_i] = value; });
    }

    static void apply_sparse_basis_to_nodal(const OrderData& data,
                                            const std::vector<Gradient>& modal_gradients,
                                            std::vector<Gradient>& nodal_gradients) {
        const std::size_t n = modal_gradients.size();
        nodal_gradients.resize(n);
        apply_low_order_combination(
            data,
            3u,
            [&](std::size_t modal_i, std::size_t component) {
                return modal_gradients[modal_i][component];
            },
            [&](std::size_t basis_i, std::size_t component, Real value) {
                nodal_gradients[basis_i][component] = value;
            });
    }

    static void apply_sparse_basis_to_nodal_to(const OrderData& data,
                                               const std::vector<Gradient>& modal_gradients,
                                               Real* SVMP_RESTRICT nodal_gradients) {
        apply_low_order_combination(
            data,
            3u,
            [&](std::size_t modal_i, std::size_t component) {
                return modal_gradients[modal_i][component];
            },
            [&](std::size_t basis_i, std::size_t component, Real value) {
                nodal_gradients[basis_i * 3u + component] = value;
            });
    }

    static void apply_sparse_basis_to_nodal(const OrderData& data,
                                            const std::vector<Hessian>& modal_hessians,
                                            std::vector<Hessian>& nodal_hessians) {
        const std::size_t n = modal_hessians.size();
        nodal_hessians.resize(n);
        apply_low_order_combination(
            data,
            9u,
            [&](std::size_t modal_i, std::size_t component) {
                return modal_hessians[modal_i].data()[component];
            },
            [&](std::size_t basis_i, std::size_t component, Real value) {
                nodal_hessians[basis_i].data()[component] = value;
            });
    }

    static void apply_sparse_basis_to_nodal_to(const OrderData& data,
                                               const std::vector<Hessian>& modal_hessians,
                                               Real* SVMP_RESTRICT nodal_hessians) {
        apply_low_order_combination(
            data,
            9u,
            [&](std::size_t modal_i, std::size_t component) {
                return modal_hessians[modal_i].data()[component];
            },
            [&](std::size_t basis_i, std::size_t component, Real value) {
                nodal_hessians[basis_i * 9u + component] = value;
            });
    }

    static void apply_sparse_basis_to_nodal_all(
        const OrderData& data,
        const std::vector<Real>& modal_values,
        const std::vector<Gradient>& modal_gradients,
        const std::vector<Hessian>& modal_hessians,
        std::vector<Real>& nodal_values,
        std::vector<Gradient>& nodal_gradients,
        std::vector<Hessian>& nodal_hessians) {
        const std::size_t n = modal_values.size();
        nodal_values.resize(n);
        nodal_gradients.resize(n);
        nodal_hessians.resize(n);
        apply_low_order_combination(
            data,
            1u,
            [&](std::size_t modal_i, std::size_t) { return modal_values[modal_i]; },
            [&](std::size_t basis_i, std::size_t, Real value) { nodal_values[basis_i] = value; });
        apply_low_order_combination(
            data,
            3u,
            [&](std::size_t modal_i, std::size_t component) {
                return modal_gradients[modal_i][component];
            },
            [&](std::size_t basis_i, std::size_t component, Real value) {
                nodal_gradients[basis_i][component] = value;
            });
        apply_low_order_combination(
            data,
            9u,
            [&](std::size_t modal_i, std::size_t component) {
                return modal_hessians[modal_i].data()[component];
            },
            [&](std::size_t basis_i, std::size_t component, Real value) {
                nodal_hessians[basis_i].data()[component] = value;
            });
    }

    static void apply_sparse_basis_to_nodal_all_to(
        const OrderData& data,
        const std::vector<Real>& modal_values,
        const std::vector<Gradient>& modal_gradients,
        const std::vector<Hessian>& modal_hessians,
        Real* SVMP_RESTRICT nodal_values,
        Real* SVMP_RESTRICT nodal_gradients,
        Real* SVMP_RESTRICT nodal_hessians) {
        apply_low_order_combination(
            data,
            1u,
            [&](std::size_t modal_i, std::size_t) { return modal_values[modal_i]; },
            [&](std::size_t basis_i, std::size_t, Real value) { nodal_values[basis_i] = value; });
        apply_low_order_combination(
            data,
            3u,
            [&](std::size_t modal_i, std::size_t component) {
                return modal_gradients[modal_i][component];
            },
            [&](std::size_t basis_i, std::size_t component, Real value) {
                nodal_gradients[basis_i * 3u + component] = value;
            });
        apply_low_order_combination(
            data,
            9u,
            [&](std::size_t modal_i, std::size_t component) {
                return modal_hessians[modal_i].data()[component];
            },
            [&](std::size_t basis_i, std::size_t component, Real value) {
                nodal_hessians[basis_i * 9u + component] = value;
            });
    }

    template <typename Sink>
    // Keep modal transform helpers free of forced-inline attributes unless
    // compiler-versioned benchmarks and LLVM IR checks show a stable benefit.
    static void apply_modal_values_to_nodal(const OrderData& data,
                                            const std::vector<Real>& modal_values,
                                            const Sink& sink) {
        const std::size_t n = modal_values.size();
        sink.resize(n);
        for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
            const Real* row = data.modal_to_nodal.data() + basis_i * n;
            Real value = Real(0);
            for (std::size_t modal_j = 0; modal_j < n; ++modal_j) {
                value += row[modal_j] * modal_values[modal_j];
            }
            sink.write(basis_i, value);
        }
    }

    template <typename Sink>
    static void apply_modal_gradients_to_nodal(const OrderData& data,
                                               const std::vector<Gradient>& modal_gradients,
                                               const Sink& sink) {
        const std::size_t n = modal_gradients.size();
        sink.resize(n);
        for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
            const Real* row = data.modal_to_nodal.data() + basis_i * n;
            Gradient gradient{};
            for (std::size_t modal_j = 0; modal_j < n; ++modal_j) {
                const Real coeff = row[modal_j];
                for (std::size_t component = 0; component < 3u; ++component) {
                    gradient[component] += coeff * modal_gradients[modal_j][component];
                }
            }
            sink.write(basis_i, gradient);
        }
    }

    template <typename Sink>
    static void apply_modal_hessians_to_nodal(const OrderData& data,
                                              const std::vector<Hessian>& modal_hessians,
                                              const Sink& sink) {
        const std::size_t n = modal_hessians.size();
        sink.resize(n);
        for (std::size_t basis_i = 0; basis_i < n; ++basis_i) {
            const Real* matrix_row = data.modal_to_nodal.data() + basis_i * n;
            Hessian hessian{};
            for (std::size_t modal_j = 0; modal_j < n; ++modal_j) {
                const Real coeff = matrix_row[modal_j];
                for (std::size_t row = 0; row < 3u; ++row) {
                    for (std::size_t col = 0; col < 3u; ++col) {
                        hessian(row, col) += coeff * modal_hessians[modal_j](row, col);
                    }
                }
            }
            sink.write(basis_i, hessian);
        }
    }

    static void apply_modal_to_nodal(const OrderData& data,
                                     const std::vector<Real>& modal_values,
                                     std::vector<Real>& nodal_values) {
        apply_modal_values_to_nodal(data, modal_values, VectorValueSink{nodal_values});
    }

    static void apply_modal_to_nodal(const OrderData& data,
                                     const std::vector<Gradient>& modal_gradients,
                                     std::vector<Gradient>& nodal_gradients) {
        apply_modal_gradients_to_nodal(data, modal_gradients, VectorGradientSink{nodal_gradients});
    }

    static void apply_modal_to_nodal(const OrderData& data,
                                     const std::vector<Hessian>& modal_hessians,
                                     std::vector<Hessian>& nodal_hessians) {
        apply_modal_hessians_to_nodal(data, modal_hessians, VectorHessianSink{nodal_hessians});
    }

    static void apply_modal_to_nodal_to(const OrderData& data,
                                        const std::vector<Real>& modal_values,
                                        Real* nodal_values) {
        apply_modal_values_to_nodal(data, modal_values, RawValueSink{nodal_values});
    }

    static void apply_modal_to_nodal_to(const OrderData& data,
                                        const std::vector<Gradient>& modal_gradients,
                                        Real* nodal_gradients) {
        apply_modal_gradients_to_nodal(data, modal_gradients, RawGradientSink{nodal_gradients});
    }

    static void apply_modal_to_nodal_to(const OrderData& data,
                                        const std::vector<Hessian>& modal_hessians,
                                        Real* nodal_hessians) {
        apply_modal_hessians_to_nodal(data, modal_hessians, RawHessianSink{nodal_hessians});
    }
};

namespace lagrange_pyramid {

const std::vector<math::Vector<Real, 3>>& nodes(int order) {
    return PyramidLagrangeCache::get(order).nodes;
}

void prewarm_scratch(int order, std::size_t max_qpts) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::prewarm_scratch(data.modal_terms.size(), max_qpts);
}

void evaluate_values(int order,
                     const math::Vector<Real, 3>& xi,
                     std::vector<Real>& values) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::evaluate_values(data, xi, values);
}

void evaluate_gradients(int order,
                        const math::Vector<Real, 3>& xi,
                        std::vector<Gradient>& gradients) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::evaluate_gradients(data, xi, gradients);
}

void evaluate_hessians(int order,
                       const math::Vector<Real, 3>& xi,
                       std::vector<Hessian>& hessians) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::evaluate_hessians(data, xi, hessians);
}

void evaluate_all(int order,
                  const math::Vector<Real, 3>& xi,
                  std::vector<Real>& values,
                  std::vector<Gradient>& gradients,
                  std::vector<Hessian>& hessians) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::evaluate_all(data, xi, values, gradients, hessians);
}

void evaluate_values_to(int order,
                        const math::Vector<Real, 3>& xi,
                        Real* SVMP_RESTRICT values_out) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::evaluate_values_to(data, xi, values_out);
}

void evaluate_gradients_to(int order,
                           const math::Vector<Real, 3>& xi,
                           Real* SVMP_RESTRICT gradients_out) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::evaluate_gradients_to(data, xi, gradients_out);
}

void evaluate_hessians_to(int order,
                          const math::Vector<Real, 3>& xi,
                          Real* SVMP_RESTRICT hessians_out) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::evaluate_hessians_to(data, xi, hessians_out);
}

void evaluate_all_to(int order,
                     const math::Vector<Real, 3>& xi,
                     Real* SVMP_RESTRICT values_out,
                     Real* SVMP_RESTRICT gradients_out,
                     Real* SVMP_RESTRICT hessians_out) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::evaluate_all_to(data, xi, values_out, gradients_out, hessians_out);
}

void evaluate_at_quadrature_points_strided(
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    const auto& data = PyramidLagrangeCache::get(order);
    PyramidLagrangeCache::evaluate_at_quadrature_points_strided(
        data, points, output_stride, values_out, gradients_out, hessians_out);
}

} // namespace lagrange_pyramid

} // namespace detail
} // namespace basis
} // namespace FE
} // namespace svmp
