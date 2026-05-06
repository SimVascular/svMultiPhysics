// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_MODEL_H
#define IONIC_MODEL_H

#include "Array.h"
#include "Parameters.h"
#include "Vector.h"

#include "FE/Common/FEException.h"

#include <string>
#include <utility>
#include <vector>

#include "CmMod.h"

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
  /// Alias for initial states vector. Each initial state is a pair of a label
  /// for that state variable and its initial value.
  /// @todo This would work better with a struct, due to the fields having
  /// meaningful names instead of first and second.
  using InitialStates = std::vector<std::pair<std::string, double>>;

  /// Constructor.
  IonicModel(const InitialStates &initial_X_, const InitialStates &initial_Xg_,
             const double Vrest_)
      : initial_X(initial_X_), initial_Xg(initial_Xg_), Vrest(Vrest_),
        Vscale(1.0), Tscale(1.0), Voffset(0.0) {}

  /// Constructor with scaling factors.
  IonicModel(const InitialStates &initial_X_, const InitialStates &initial_Xg_,
             const double Vrest_, const double Vscale_, const double Tscale_,
             const double Voffset_)
      : initial_X(initial_X_), initial_Xg(initial_Xg_), Vrest(Vrest_),
        Vscale(Vscale_), Tscale(Tscale_), Voffset(Voffset_) {}

  /// Virtual destructor.
  virtual ~IonicModel() = default;

  /**
   * @brief Construct an instance of model parameters for this model.
   */
  virtual std::unique_ptr<IonicModelParameters> get_parameters() const {
    return nullptr;
  };

  /**
   * @brief Read model parameters from a parameter object.
   */
  virtual void read_parameters(const IonicModelParameters &params);

  /**
   * @brief Distribute model parameters to all parallel processes.
   */
  virtual void distribute_parameters(const CmMod &cm_mod, const cmType &cm);

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

  /**
   * @brief Get initial conditions for the model.
   */
  const InitialStates &get_initial_X() const { return initial_X; }

  /**
   * @brief Get initial gating variables for the model.
   */
  const InitialStates &get_initial_Xg() const { return initial_Xg; }

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
  InitialStates initial_X;

  /// Initial gating variables.
  InitialStates initial_Xg;

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

/**
 * @brief Self-registering factory for ionic models.
 *
 * This class gives a way to register ionic models when they are defined, and
 * then instantiate concrete ionic models, derived from IonicModel, by name. To
 * be compatible with this factory, classes derived from IonicModel must be
 * default constructible.
 *
 * It combines the
 * [factory](https://en.wikipedia.org/wiki/Abstract_factory_pattern) and
 * [singleton](https://en.wikipedia.org/wiki/Singleton_pattern) patterns. There
 * should always exist only one instance of this class, which cannot be accessed
 * directly but only manipulated through the static methods of this class.
 *
 * To register a new ionic model into the factory, you can call the
 * register_model static method, passing a class derived from Ionic as template
 * argument and a label for the model as argument. A shortcut for this is to
 * use the macro REGISTER_IONIC_MODEL.
 */
class IonicModelFactory {
public:
  /**
   * @brief Register a child model.
   */
  template <class Model> static bool register_model(const std::string &name) {
    auto &factory_instance = get_instance();

    if (factory_instance.children.find(name) !=
        factory_instance.children.end()) {
      svmp::raise<svmp::FE::InvalidArgumentException>(
          SVMP_HERE,
          "A model with name '" + name +
              "' was already registered in the ionic model factory.");
    }

    factory_instance.children[name] = []() {
      return std::make_unique<Model>();
    };

    return true;
  }

  /**
   * @brief Instantiate a model from its name.
   */
  static std::unique_ptr<IonicModel> create_model(const std::string &name);

  /**
   * @brief Iterate through registered ionic models.
   *
   * For every registered ionic model, creates a dummy instance of it, and then
   * calls the provided function on that model. All the dummy model instances
   * are destroyed after the function call.
   */
  static void
  visit(const std::function<void(const std::string &, const IonicModel &)> &f);

protected:
  /**
   * @brief Default constructor.
   */
  IonicModelFactory() = default;

  /**
   * @brief Access the singleton instance.
   */
  static IonicModelFactory &get_instance() {
    static IonicModelFactory instance;
    return instance;
  }

  /**
   * @brief Registered ionic models.
   *
   * Each ionic model is represented by a function that takes no argument and
   * returns a unique_ptr<IonicModel> constructing an instance of that model.
   * This requires classes derived from IonicModel to be default
   * constructible.
   */
  std::map<std::string, std::function<std::unique_ptr<IonicModel>()>> children;
};

/**
 * @brief Macro to register a ionic model in the factory.
 */
#define REGISTER_IONIC_MODEL(name, type)                                       \
  namespace IonicModelFactoryInternals {                                       \
  static inline volatile const bool ionic_model_factory_registered_##type =    \
      IonicModelFactory::register_model<type>(name);                           \
  }

#endif