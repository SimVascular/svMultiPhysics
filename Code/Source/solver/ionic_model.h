// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_MODEL_H
#define IONIC_MODEL_H

#include "Array.h"
#include "Vector.h"

#include <map>
#include <string>
#include <utility>

/**
 * @brief Abstract ionic model class.
 *
 * ### Mathematical model
 *
 * This class represents an abstract ionic model, i.e. an ODE system in the form
 * @f[ \begin{aligned}
 *   \frac{\partial v}{\partial t} + I_\text{ion}(v, \mathbf{w}) &=
 *      I_\text{ext}(t) \\
 *   \frac{\partial \mathbf{w}}{\partial t} &= \mathbf{F}(v, \mathbf{w})
 * \end{aligned} @f]
 * where @f$v@f$ is the transmembrane potential, @f$\mathbf{w}@f$ is a vector of
 * ionic model variables (which may be gating variables or ionic
 * concentrations), @f$I_\text{ion}@f$ is the ionic current, and
 * @f$I_\text{ext}@f$ is an externally applied current.
 *
 * Individual models differ in the number of state variables @f$\mathbf{w}@f$
 * and in the expressions of @f$I_\text{ion}@f$ and @f$\mathbf{F}@f$. These are
 * specified by classes derived from this.
 *
 * **References**:
 *  - Colli Franzone, Pavarino, Scacchi. Mathematical Cardiac Electrophysiology.
 *    Springer, 2014.
 *
 * @todo Document stretch-activated current.
 *
 * ### Numerical methods
 *
 * @todo
 *
 * ### Implementation details
 *
 * @todo
 *
 * @todo The number of state variables in the model should be moved to this
 * class.
 */
class IonicModel {
public:
  /// State variable information. Bundles a string label for the variable, to be
  /// used for parameter files, and the variable's initial value.
  class StateVariable {
  public:
    std::string label;
    double initial_value;
  };

  /// Alias for state variables vector.
  using StateVector = std::vector<StateVariable>;

  /// Class for managing initial values for a generic ionic model.
  class InitialValues : public ParameterLists {};

  /// Constructor.
  IonicModel(const StateVector &states_X_, const StateVector &states_Xg_,
             const double Vrest_)
      : states_X(states_X_), states_Xg(states_Xg_), Vrest(Vrest_), Vscale(1.0),
        Tscale(1.0), Voffset(0.0) {}

  /// Constructor with scaling factors.
  IonicModel(const StateVector &states_X_, const StateVector &states_Xg_,
             const double Vrest_, const double Vscale_, const double Tscale_,
             const double Voffset_)
      : states_X(states_X_), states_Xg(states_Xg_), Vrest(Vrest_),
        Vscale(Vscale_), Tscale(Tscale_), Voffset(Voffset_) {}

  /// Virtual destructor.
  virtual ~IonicModel() = default;

  /**
   * @brief Setup model initial conditions.
   *
   * @param[in] nX Number of state variables.
   * @param[in] nG Number of gating variables.
   * @param[out] X Vector of state variables to be initialized.
   * @param[out] Xg Vector of gating variables to be initialized.
   */
  void init(const int nX, const int nG, Vector<double> &X,
            Vector<double> &Xg) const;

  /**
   * @name Integration methods.
   * @{
   */

  /**
   * @brief Integrate the model with the Crank-Nicolson method.
   *
   * @todo Document numerical formulation.
   *
   * @todo IPAR and RPAR can probably be made const here. Also, the meaning of
   *       their entries needs to be documented.
   */
  void integ_cn2(const unsigned int zone_id, const int nX, const int nG,
                 Vector<double> &X, Vector<double> &Xg, const double Ts,
                 const double Ti, const double Istim, const double Ksac,
                 Vector<int> &IPAR, Vector<double> &RPAR) const;

  /**
   * @brief Integrate the model with the forward Euler method.
   *
   * @todo Document numerical formulation.
   */
  void integ_fe(const unsigned int zone_id, const int nX, const int nG,
                Vector<double> &X, Vector<double> &Xg, const double Ts,
                const double Ti, const double Istim, const double Ksac) const;

  /**
   * @brief Integrate the model with the Runge-Kutta method.
   *
   * @todo Document numerical formulation.
   */
  void integ_rk(const unsigned int zone_id, const int nX, const int nG,
                Vector<double> &X, Vector<double> &Xg, const double Ts,
                const double Ti, const double Istim, const double Ksac) const;

  /**
   * @}
   */

protected:
  /**
   * @brief Update variables with analytical solution.
   *
   * @todo Extend documentation.
   */
  virtual void update_g(const unsigned int zone_id, const double dt,
                        const int nX, const int nG, const Vector<double> &X,
                        Vector<double> &Xg) const = 0;

  /**
   * @brief Model right hand side.
   *
   * Defines the right-hand side function for the potential and ionic equations.
   * Must be ovverridden in derived classes.
   *
   * @todo Document the meaning of the individual parameters.
   */
  virtual void getf(const unsigned int zone_id, const int nX, const int nG,
                    const Vector<double> &X, const Vector<double> &Xg,
                    Vector<double> &f, const double I_stim,
                    const double I_sac) const = 0;

  /**
   * @brief Model jacobian.
   *
   * Defines the jacobian matrix of the model equations, that is the matirx of
   * derivatives of the function evaulated by getf.
   *
   * @todo Document the meaning of the individual parameters.
   */
  virtual void getj(const unsigned int zone_id, const int nX, const int nG,
                    const Vector<double> &X, const Vector<double> &Xg,
                    Array<double> &Jac, const double Ksac) const = 0;

  /// Initial states.
  StateVector states_X;

  /// Initial gating variables.
  StateVector states_Xg;

  /// Resting transmembrane potential. It is used to define the
  /// stretch-activated current.
  const double Vrest;

  /**
   * @name Scaling factors.
   *
   * Individual ionic models may need to rescale the time or voltage variable,
   * e.g. to bring them into dimensionless form. These are the factors used for
   * that purpose. They are assigned in the constructor of this class.
   *
   * @todo Document units of measure.
   *
   * @{
   */

  /// Voltage scaling.
  const double Vscale;

  /// Time scaling.
  const double Tscale;

  /// Voltage offset parameter.
  const double Voffset;

  /**
   * @}
   */
};

#endif