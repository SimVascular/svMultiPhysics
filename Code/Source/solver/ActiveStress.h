// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_H
#define ACTIVE_STRESS_H

#include "Array.h"
#include "Parameters.h"
#include "Vector.h"
#include "consts.h"
#include "factory.h"
#include "mat_fun.h"

#include "CmMod.h"

#include <memory>

/**
 * @brief Return whether a certain equation type can be solved with active
 * stress.
 */
bool supports_active_stress(const consts::EquationType eq_type);

/**
 * @brief Abstract active stress class.
 *
 * ## Table of contents
 *
 * - @ref activestress-overview
 * - @ref activestress-directions
 * - @ref activestress-implementing
 * - @ref activestress-coupling
 *
 * ## Overview {#activestress-overview}
 *
 * This class provides an interface for defining active stress models, i.e.
 * models that, in the context of structural mechanics of muscular tissue,
 * compute an active tension representing the contribution of muscular
 * contraction to the constitiutive law.
 *
 * The class assumes that the active tension can be expressed as
 * @f[
 *   \Tact = \Tact(t, \calcium, \fiberstretch, \fiberstretchrate,
 *                 \astressstate),
 * @f]
 * where @f$\calcium@f$ is the intracellular calcium concentration,
 * @f$\fiberstretch@f$ is the fiber stretch, @f$\fiberstretchrate@f$ is the
 * fiber stretch rate, and @f$\astressstate@f$ is a vector of internal state
 * variables, representing the state of contraction.
 *
 * The expression assumed above implies that the active tension is a local
 * function of the variables it depends on, that is the active tension at a
 * given point only depends on the value of other variables at that same point.
 * Accordingly, this class works nodally, storing the state of contraction at
 * every mesh node and evaluating the active tension from it, on demand,
 * through @ref ActiveStress::compute_tension.
 *
 * ## Directional distribution of active stress {#activestress-directions}
 *
 * In muscular mechanics models, active stress normally acts only along the
 * direction of fibers @f$\fiberdirection@f$, reflecting the fact that
 * contractile units are aligned with fibers. However, one might want to account
 * for fiber dispersion, i.e. the fact that fibers are not perfectly and
 * regularly aligned, but rather have a certain distribution of orientations
 * centered around the principal direction @f$\fiberdirection@f$.
 *
 * This can be surrogated by defining the active stress tensor as
 * @f[
 *   S_\text{act} = \Tact \left(
 *     \eta_f \fiberdirection \otimes \fiberdirection +
 *     \eta_s \sheetdirection \otimes \sheetdirection +
 *     \eta_n \sheetnormaldirection \otimes \sheetnormaldirection
 *   \right),
 * @f]
 * where @f$\sheetdirection@f$ and @f$\sheetnormaldirection@f$ are the sheet and
 * sheet-normal directions, respectively, and @f$\eta_f@f$, @f$\eta_s@f$, and
 * @f$\eta_n@f$ are coefficients that define the distribution the active tension
 * along the three principal directions. The coefficients must be such that
 * @f$\eta_f + \eta_s + \eta_n = 1@f$.
 *
 * This class stores the values of @f$\eta_f@f$, @f$\eta_s@f$ and @f$\eta_n@f$,
 * and applies them in @ref ActiveStress::compute_tension, which returns
 * @f$\eta_f \Tact@f$, @f$\eta_s \Tact@f$ and @f$\eta_n \Tact@f$, bundled in an
 * @ref ActiveTension.
 *
 * ## Implementing concrete active stress models {#activestress-implementing}
 *
 * To implement a new active stress model, the following steps need to be taken:
 *
 * 1. Create a new class derived from @ref ActiveStress.
 * 2. Override the methods @ref init_local, @ref advance_time_step_local,
 *    @ref compute_active_tension_local and
 *    @ref compute_active_tension_derivative_local, defining the initial
 *    condition, time evolution, active tension and its partial derivative with
 *    respect to the fiber stretch, respectively, for a single node.
 * 3. Create a new class derived from @ref ActiveStressModelParameters to store
 *    the parameters specific to the new active stress model.
 * 4. Override the methods @ref get_parameters,
 *    @ref read_model_specific_parameters and
 *    @ref distribute_model_specific_parameters to manage the parameters of the
 *    new active stress model.
 * 5. Register the new class into the active stress model factory by using the
 *    macro @ref REGISTER_ACTIVE_STRESS_MODEL. The macro should be called in a
 *    `.cpp` file, not in a header file.
 *
 * Notice that if the model is expressed in terms of a system of ODEs, it can
 * be implemented by deriving from @ref ActiveStressODE, which already addresses
 * some of the points above.
 *
 * ## Coupling with the structural mechanics problem {#activestress-coupling}
 *
 * The active tension depends on the fiber stretch @f$\fiberstretch@f$ both
 * directly, through the expression of @f$\Tact@f$, and indirectly, through the
 * state @f$\astressstate@f$, which is itself driven by the fiber stretch. The
 * mechanics problem, in turn, depends on @f$\Tact@f$.
 *
 * Explicit time discretization for the direct dependence was observed to lead
 * to instabilities. Accordingly, the direct dependence is discretized
 * implicitly, that is the active tension is recomputed within the nonlinear
 * iterations for the structure problem that uses ActiveStress. To facilitate
 * the convergence of nonlinear iterations, this class also allows to compute
 * the derivative @f$\frac{\partial\Tact}{\partial\fiberstretch}@f$, which is
 * used to assemble tangent terms associated to this in the structural system.
 *
 * Indirect dependence was not observed to give rise to instabilities.
 * Accordingly, it is discretized explicitly by default, meaning that the state
 * is updated once per time step evaluating the fiber stretch
 * @f$\fiberstretch@f$ using the displacement from the previous time step. The
 * user can change this behavior by setting the parameter @c
 * Implicit_state_coupling to @c true in the XML file. This will make the state
 * update every nonlinear iteration. No tangent terms are computed for this
 * contribution, so nonlinear iterations can be expected to converge more slowly
 * when this is enabled.
 */
class ActiveStress {
public:
  /**
   * @brief Active tension information at a point.
   *
   * This struct bundles the active tension along the three principal
   * directions (fibers @f$\fiberdirection@f$, sheets @f$\sheetdirection@f$ and
   * sheet normals @f$\sheetnormaldirection@f$) and their partial derivatives
   * with respect to the fiber stretch.
   *
   * It is a convenience data structure used to pass this information to
   * functions that consume active tension information (e.g. the structural
   * mechanics assembly functions).
   */
  struct ActiveTension {
    /// Tension along the fiber direction, @f$\eta_f \Tact@f$.
    double fibers = 0.0;

    /// Tension along the sheet direction, @f$\eta_s \Tact@f$.
    double sheets = 0.0;

    /// Tension along the sheet-normal direction, @f$\eta_n \Tact@f$.
    double sheet_normals = 0.0;

    /// Derivative of @ref fibers with respect to the fiber stretch, at fixed
    /// state, @f$\eta_f \pdv*{\Tact}{\fiberstretch}@f$.
    double d_fibers = 0.0;

    /// Derivative of @ref sheets with respect to the fiber stretch, at fixed
    /// state, @f$\eta_s \pdv*{\Tact}{\fiberstretch}@f$.
    double d_sheets = 0.0;

    /// Derivative of @ref sheet_normals with respect to the fiber stretch, at
    /// fixed state, @f$\eta_n \pdv*{\Tact}{\fiberstretch}@f$.
    double d_sheet_normals = 0.0;
  };

  /**
   * @brief Helper to evaluate active tension at quadrature points.
   *
   * The active stress model stores the state at the degrees of freedom, but the
   * active tension needs to be evaluated at quadrature points. This class
   * allows to do that by interpolating the state to quadrature points and then
   * evaluating the active tension there.
   */
  class Evaluator {
  public:
    /**
     * @brief Update the evaluator for a given element.
     *
     * @param[in] active_stress Active stress model of the domain the element
     *   belongs to.
     * @param[in] nodes Indices of the mesh nodes of the element.
     */
    void update(const ActiveStress &active_stress, const Vector<int> &nodes);

    /**
     * @brief Reset this evaluator so that @ref evaluate returns zero tension,
     * until the next call to @ref update.
     *
     * Used for elements whose domain has no active stress model.
     */
    void clear() { active_stress_ = nullptr; }

    /**
     * @brief Compute the active tension at a quadrature point.
     *
     * @param[in] N Shape functions at the quadrature point, of the same nodes
     *   the state was gathered at by @ref update.
     * @param[in] F Deformation gradient at the quadrature point.
     * @param[in] fN Fiber directions of the element, one per column, the first
     *   column being the fiber direction itself.
     *
     * @return Active tension along fibers, sheets and sheet normals, and their
     *   derivatives with respect to the fiber stretch, bundled in an object of
     *   type @ref ActiveTension.
     *
     * @tparam nsd Number of spatial dimensions.
     */
    template <int nsd>
    ActiveTension
    evaluate(const Eigen::Ref<const Eigen::VectorXd> &N,
             const mat_fun::Matrix<nsd> &F,
             const Eigen::Ref<const Eigen::Matrix<double, nsd, Eigen::Dynamic>>
                 &fN) const;

  private:
    /// Active stress model of the domain the element belongs to, or null if
    /// @ref clear was called last, or if this evaluator was never updated.
    const ActiveStress *active_stress_ = nullptr;

    /// State variables at the element nodes, of size (n_states, element
    /// nodes).
    Eigen::MatrixXd state_;

    /// State variables interpolated at the quadrature point, of size
    /// n_states. Scratch storage for @ref evaluate, mutable so that it can
    /// stay const.
    mutable Vector<double> interpolated_state;
  };

  /// Grants @ref Evaluator direct access to @ref states.
  friend class Evaluator;

  /**
   * @brief Constructor.
   *
   * @param n_states_ Number of state variables for this model.
   * @param needs_fiber_stretch Whether this model uses the fiber stretch
   *   passed to @ref update. This flag can be used to determine whether fiber
   *   stretch computation can be skipped for efficiency.
   * @param needs_fiber_stretch_rate Whether this model uses the fiber stretch
   *   rate passed to @ref update. This flag can be used to determine whether
   *   fiber stretch rate computation can be skipped for efficiency.
   */
  ActiveStress(const unsigned int n_states_, const bool needs_fiber_stretch,
               const bool needs_fiber_stretch_rate)
      : n_states(n_states_), needs_fiber_stretch_(needs_fiber_stretch),
        needs_fiber_stretch_rate_(needs_fiber_stretch_rate) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~ActiveStress() = default;

  /**
   * @brief Construct an instance of model parameters for this model.
   */
  virtual std::unique_ptr<ActiveStressModelParameters>
  get_parameters() const = 0;

  /**
   * @brief Read model parameters from a parameter object.
   */
  void read_parameters(const ActiveStressParameters &params);

  /**
   * @brief Distribute model parameters to all parallel processes.
   */
  void distribute_parameters(const CmMod &cm_mod, const cmType &cm);

  /**
   * @brief Get the state vector at a given mesh node.
   *
   * @param[in] Ac Index of the mesh node.
   */
  Vector<double> get_state(const unsigned int Ac) const {
    return states.col(Ac);
  }

  /**
   * @brief Compute the active tension given the state vector and fiber stretch.
   *
   * @param[in] state State vector at the point.
   * @param[in] fiber_stretch Fiber stretch at the point.
   *
   * @return Active tension along fibers, sheets and sheet normals, and their
   *   derivatives with respect to the fiber stretch, bundled in an object of
   *   type @ref ActiveTension.
   */
  ActiveTension compute_tension(const Vector<double> &state,
                                const double fiber_stretch) const {
    const double tension = compute_active_tension_local(state, fiber_stretch);
    const double derivative =
        compute_active_tension_derivative_local(state, fiber_stretch);

    return {eta_f * tension,    eta_s * tension,    eta_n * tension,
            eta_f * derivative, eta_s * derivative, eta_n * derivative};
  }

  /**
   * @brief Initialize the model.
   *
   * Allocates the internal state vector and initializes it with the model's
   * initial conditions.
   *
   * @param[in] tnNo Total number of mesh nodes for the current rank.
   */
  virtual void init(const unsigned int tnNo);

  /**
   * @brief Begin a new time step.
   *
   * Stores the current state as the initial condition of the time step. Must be
   * called once per time step, before any call to @ref update.
   */
  virtual void time_advance();

  /**
   * @brief Update the state over the current time step.
   *
   * Advances the state stored by @ref time_advance over one time step, using
   * the given calcium, fiber stretch and fiber stretch rate.
   *
   * This function may be called more than once per time step: every call
   * restarts from the state stored by @ref time_advance, so the resulting state
   * depends only on the arguments of the last call.
   *
   * @param[in] t Current time (i.e. the time instant being advanced to).
   * @param[in] dt Time step size.
   * @param[in] calcium Calcium concentration at every node.
   * @param[in] fiber_stretch Fiber stretch at every node. This is usually
   *   computed with post::fib_stretch.
   * @param[in] fiber_stretch_rate Fiber stretch rate at every node. This is
   *   usually computed with post::fib_stretch_rate.
   */
  virtual void update(const double t, const double dt,
                      const Vector<double> &calcium,
                      const Vector<double> &fiber_stretch,
                      const Vector<double> &fiber_stretch_rate);

  /**
   * @brief Whether this model uses the fiber stretch passed to @ref update.
   * This flag can be used to determine whether fiber stretch computation can be
   * skipped for efficiency.
   */
  bool needs_fiber_stretch() const { return needs_fiber_stretch_; }

  /**
   * @brief Whether this model uses the fiber stretch rate passed to
   * @ref update. This flag can be used to determine whether fiber stretch rate
   * computation can be skipped for efficiency.
   */
  bool needs_fiber_stretch_rate() const { return needs_fiber_stretch_rate_; }

  /**
   * @brief Whether the state of this model is updated within the nonlinear
   * iterations of the mechanics problem, i.e. whether the indirect dependence
   * of the active tension on the fiber stretch is treated implicitly.
   */
  bool implicit_state_coupling() const { return implicit_state_coupling_; }

  /// Number of state variables for this model.
  const unsigned int n_states;

protected:
  /**
   * @brief Read model parameters from a parameter object.
   *
   * This method needs to be overridden by derived classes to read the
   * parameters specific to the concrete model they implement.
   */
  virtual void
  read_model_specific_parameters(const ActiveStressModelParameters &params) = 0;

  /**
   * @brief Distribute model parameters to all parallel processes.
   *
   * This method needs to be overridden by derived classes to distribute the
   * parameters specific to the concrete model they implement to all parallel
   * processes.
   */
  virtual void distribute_model_specific_parameters(const CmMod &cm_mod,
                                                    const cmType &cm) = 0;

  /**
   * @brief Initialize the state vector for a single node.
   *
   * @param[out] state State vector for a single node, to be initialized by
   *   this function.
   */
  virtual void init_local(Vector<double> &state) const = 0;

  /**
   * @brief Advance in time for a single node.
   *
   * @param[in] t Current time (i.e. the time instant being advanced to).
   * @param[in] dt Time step size.
   * @param[in] calcium Calcium concentration at the current node.
   * @param[in] fiber_stretch Fiber stretch at the current node.
   * @param[in] fiber_stretch_rate Fiber stretch rate at the current node.
   * @param[in,out] state State vector for a single node, to be updated by
   *   this function.
   */
  virtual void advance_time_step_local(const double t, const double dt,
                                       const double calcium,
                                       const double fiber_stretch,
                                       const double fiber_stretch_rate,
                                       Vector<double> &state) const = 0;

  /**
   * @brief Compute the active tension for a single node.
   *
   * @param[in] state State vector for a single node.
   * @param[in] fiber_stretch Fiber stretch at the current node.
   */
  virtual double
  compute_active_tension_local(const Vector<double> &state,
                               const double fiber_stretch) const = 0;

  /**
   * @brief Compute the partial derivative of the active tension with respect to
   * the fiber stretch, at fixed state, for a single node.
   *
   * @param[in] state State vector for a single node.
   * @param[in] fiber_stretch Fiber stretch at the current node.
   */
  virtual double
  compute_active_tension_derivative_local(const Vector<double> &state,
                                          const double fiber_stretch) const = 0;

  /**
   * @brief Backing store for @ref needs_fiber_stretch.
   */
  bool needs_fiber_stretch_;

  /**
   * @brief Backing store for @ref needs_fiber_stretch_rate.
   */
  bool needs_fiber_stretch_rate_;

  /// Time instant being advanced to. Set by @ref update.
  double time = 0.0;

  /// State variables for the model.
  Array<double> states;

  /**
   * @brief State variables at the beginning of the current time step.
   *
   * Set by @ref time_advance and used by @ref update as the initial condition
   * of every call within the time step.
   */
  Array<double> states_at_time_step_start;

  /**
   * @brief Whether the state of this model is updated within the nonlinear
   * iterations of the mechanics problem.
   */
  bool implicit_state_coupling_;

  /// Active tension coefficient along the fiber direction.
  double eta_f;

  /// Active tension coefficient along the sheet direction.
  double eta_s;

  /// Active tension coefficient along the sheet-normal direction.
  double eta_n;
};

/**
 * @brief Alias for the active stress model factory.
 *
 * See the documentation for @ref Factory for more details on how this works.
 */
using ActiveStressFactory = Factory<ActiveStress>;

/**
 * @brief Macro to register an active stress model in the factory.
 */
#define REGISTER_ACTIVE_STRESS_MODEL(name, type)                               \
  REGISTER_IN_FACTORY(ActiveStress, type, name)

#endif