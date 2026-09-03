#ifndef DARCY_H
#define DARCY_H

#include "ComMod.h"
#include "SolutionStates.h"

/**
 * @brief Pressure-based Darcy flow in porous media.
 *
 * This namespace implements the Darcy equation for perfusion of porous media
 * with intrinsic dimensions 2 and 3. Material coefficients are homogeneous
 * within each solver domain, permeability is isotropic, and Stokes-flow
 * assumptions apply. The assembled pressure strong form is
 * \f[
 *   \rho \beta \frac{\partial p}{\partial t}
 *   - \nabla \cdot \left(\frac{\rho K}{\mu}\nabla p\right)
 *   = \rho s.
 * \f]
 *
 * This discretizes the pressure-only strong form.
 * Velocity is not an independent unknown.
 * After pressure is solved, Darcy velocity is
 * evaluated as the derived field
 * \f[
 *   \boldsymbol{q} = -\frac{K}{\mu}\nabla p.
 * \f]
 *
 * The model quantities and their admissible ranges are:
 * - \f$p\f$: pressure unknown.
 * - \f$\boldsymbol{q}\f$: derived Darcy velocity.
 * - \f$K\f$: configured intrinsic scalar permeability. `Darcy_permeability`
 *   defaults to \f$10^{-15}\f$ and must satisfy \f$K > 0\f$.
 * - \f$\mu\f$: configured dynamic viscosity. `Darcy_fluid_viscosity` defaults
 *   to 1 and must satisfy \f$\mu > 0\f$.
 * - \f$\rho\f$: configured reference fluid density. `Fluid_density` defaults
 *   to 0.5 and must satisfy \f$\rho > 0\f$.
 * - \f$\beta\f$: configured storage/compressibility.
 *   `Darcy_media_compressibility` defaults to 0 and must satisfy
 *   \f$\beta \ge 0\f$.
 * - \f$s\f$: configured volumetric source provided by `Source_term`; it
 *   defaults to 0 and is constant within each configured domain.
 *
 * @par Darcy flux output
 * On supported two- and three-dimensional meshes, the pressure gradient is
 * reconstructed in the mesh coordinates and the derived Darcy flux is
 * \f[
 *   \boldsymbol{q} = -\frac{K}{\mu}\nabla p.
 * \f]
 *
 * @par Cardiovascular porous-flow context
 * The following works describe future multi-compartment and microcirculation
 * model extensions than the single-compartment formulation implemented here:
 * - C. Michler et al., "A computationally efficient framework for the
 *   simulation of cardiac perfusion using a multi-compartment Darcy
 *   porous-media flow model," DOI
 *   <a href="https://doi.org/10.1002/cnm.2520">10.1002/cnm.2520</a>.
 * - G. Montino Pelagi et al., "Modeling cardiac microcirculation for the
 *   simulation of coronary flow and 3D myocardial perfusion," DOI
 *   <a href="https://doi.org/10.1007/s10237-024-01873-z">10.1007/s10237-024-01873-z</a>.
 */
namespace darcy {

    /// Validate the configured Darcy material coefficients for a domain.
    /// @param[in] domain Solver domain whose material properties are checked.
    void validate_material_properties(const dmnType& domain);

    /// Reject mesh types that Darcy assembly and flux output do not implement.
    /// @param[in] mesh Mesh whose element type is checked.
    /// @note `mshType::lFib` denotes a one-dimensional mesh embedded in the
    /// ambient geometry, not a myocardial material fiber direction. Darcy is
    /// currently limited to intrinsic dimensions 2 and 3.
    void validate_element_support(const mshType& mesh);

    /// Assemble a Darcy boundary contribution into the element residual.
    /// @param[in] com_mod Common solver state retained for the common assembly interface.
    /// @param[in] eNoN Number of element nodes.
    /// @param[in] w Weighted boundary quadrature measure.
    /// @param[in] N Shape-function values at the quadrature point.
    /// @param[in] h Prescribed boundary flux contribution.
    /// @param[in,out] lR Element residual.
    void b_darcy(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const double h, Array<double>& lR);

    /// Assemble Darcy volume contributions for all supported elements in a mesh.
    /// @param[in,out] com_mod Common solver state and assembly interface.
    /// @param[in] lM Mesh whose Darcy elements are assembled.
    /// @param[in] solutions Solution states used for element-local fields.
    void construct_darcy(ComMod& com_mod, const mshType& lM, const SolutionStates& solutions);

    /// Assemble the residual and tangent for an intrinsic two-dimensional element.
    /// @param[in] com_mod Common solver state.
    /// @param[in] eNoN Number of element nodes.
    /// @param[in] w Weighted volume quadrature measure.
    /// @param[in] N Shape-function values at the quadrature point.
    /// @param[in] Nx Mapped spatial shape-function derivatives.
    /// @param[in] al Element-local pressure rates.
    /// @param[in] yl Element-local pressure state.
    /// @param[in,out] lR Element residual.
    /// @param[in,out] lK Element tangent matrix.
    void darcy_2d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
                  const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK);

    /// Assemble the residual and tangent for an intrinsic three-dimensional element.
    /// @param[in] com_mod Common solver state.
    /// @param[in] eNoN Number of element nodes.
    /// @param[in] w Weighted volume quadrature measure.
    /// @param[in] N Shape-function values at the quadrature point.
    /// @param[in] Nx Mapped spatial shape-function derivatives.
    /// @param[in] al Element-local pressure rates.
    /// @param[in] yl Element-local pressure state.
    /// @param[in,out] lR Element residual.
    /// @param[in,out] lK Element tangent matrix.
    void darcy_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
                  const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK);
}

#endif //DARCY_H
