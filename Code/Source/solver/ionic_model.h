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

// Forward declarations.
class outputType;

/**
 * @brief Enumeration of time integration types for ODEs.
 */
enum class TimeIntegrationType {
  NA = 200,  ///< Undefined integration method.
  FE = 201,  ///< Forward Euler.
  RK4 = 202, ///< 4th order explicit Runge-Kutta
  CN2 = 203  ///< Crank-Nicolson.
};

/// Map from a string representation of the integration type to the actual type.
extern const std::map<std::string, TimeIntegrationType> cep_time_int_to_type;

/// Write a TimeIntegrationType to an output stream.
static std::ostream &operator<<(std::ostream &strm, TimeIntegrationType type) {
  const std::map<TimeIntegrationType, std::string> names = {
      {TimeIntegrationType::NA, "NA"},
      {TimeIntegrationType::FE, "FE"},
      {TimeIntegrationType::RK4, "RK4"},
      {TimeIntegrationType::CN2, "CN2"},
  };

  return strm << names.at(type);
}

/// @brief Time integration scheme and related parameters
class odeType {
public:
  odeType() {};

  /// @brief Time integration method type
  TimeIntegrationType tIntType = TimeIntegrationType::NA;

  /// @brief Max. iterations for Newton-Raphson method
  int maxItr = 5;

  /// @brief Absolute tolerance
  double absTol = 1.E-8;

  /// @brief Relative tolerance
  double relTol = 1.E-4;
};

/**
 * @brief Abstract ionic model class.
 *
 * ### Mathematical model
 *
 * This class represents an abstract ionic model, i.e. an ODE system in the form
 * @f[ \begin{aligned}
 *   \frac{\partial v}{\partial t} + I_\text{ion}(v, \mathbf{w}) &=
 *      I_\text{ext}(t) \\
 *   \frac{\partial \mathbf{w}}{\partial t} &= \mathbf{f}(v, \mathbf{w})
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
 *  - Goktepe, Kuhl. Computational Modeling of Cardiac Electrophysiology: a
 *    Novel Finite Element Approach. International Journal for Numerical Methods
 *    in Biomedical Engineering, 2009.
 *
 * ### Numerical methods
 *
 * The vector of state variables @f$\mathbf{w}@f$ is partitioned into two
 * groups, the state variables @f$\mathbf{x}@f$ (including the transmembrane
 * potential @f$v@f$) and the gating variables @f$\mathbf{x}_g@f$. The ODE
 * system is advanced through a Rush-Larsen scheme, in which the gating
 * variables are updated through an analytical expression and the state
 * variables through a time-stepping scheme.
 *
 * This class currently supports forward Euler (FE), fourth-order explicit
 * Runge-Kutta (RK) and Crank-Nicolson (CN) for advancing the state variables.
 * Refer to the documentation of the functions integ_fe, integ_rk and integ_cn2
 * for details on the individual methods. Beware that the CN method is only
 * supported for ionic models that implement the getj method, and an exception
 * is raised otherwise.
 *
 * ### Implementing concrete ionic models
 *
 * To implement a new ionic model, the following steps need to be taken:
 *
 * 1. Create a class derived from @ref IonicModel.
 * 2. Override the methods @ref init, @ref update_g, @ref getf for that ionic
 *    model. If implicit solvers are desired, also override @ref getj.
 * 3. Create a class derived from @ref IonicModelParameters to manage the
 *    parameters specific to the new ionic model.
 * 4. Override the methods @ref get_parameters, @ref read_parameters and @ref
 *    distribute_parameters to manage the parameters of the new ionic model.
 * 5. Override the @ref get_calcium_index method to return the index of the
 *    calcium proxy variable.
 * 6. Register the new class into the ionic model factory by using the macro
 *    @ref REGISTER_IONIC_MODEL. The macro should be called in a `.cpp`
 *    file, not in a header file.
 * 7. Edit the files `CepMod.h` and `CepMod.cpp` to add the label for the new
 *    ionic model type.
 */
class IonicModel {
public:
  /// Alias for initial states vector. Each initial state is a pair of a label
  /// for that state variable and its initial value.
  /// @todo[michelebucelli] This would work better with a struct, due to the
  /// fields having meaningful names instead of first and second.
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
   *
   * By default, this method only takes care of the parameters related to
   * initial states. If a derived ionic model has other parameters, then this
   * method should be overridden to read those as well.
   */
  virtual void read_parameters(const IonicModelParameters &params);

  /**
   * @brief Distribute model parameters to all parallel processes.
   *
   * By default, this method only takes care of the parameters related to
   * initial states. If a derived ionic model has other parameters, then this
   * method should be overridden to distribute those as well.
   */
  virtual void distribute_parameters(const CmMod &cm_mod, const cmType &cm);

  /**
   * @brief Setup model initial conditions.
   *
   * @param[out] X Vector of state variables to be initialized.
   * @param[out] Xg Vector of gating variables to be initialized.
   */
  void init(Vector<double> &X, Vector<double> &Xg) const;

  /**
   * @brief Integrate over one time step.
   *
   * @param[in] ode_solver_params ODE solver parameters structure, including
   *   the solver method and the stopping criterion information.
   * @param[in] zone_id Identifier for the transmural zone (epicardium,
   *   endocardium, myocardium).
   * @param[in] t Current time.
   * @param[in] dt Integration time step.
   * @param[in] Istim Amplitude of the stimulus current at the current time.
   * @param[in] Ksac Amplitude of the stretch-activated current.
   * @param[in,out] X Vector of state variables to be updated.
   * @param[in,out] Xg Vector of gating variables to be updated.
   */
  void integ(const odeType &ode_solver_params, const int zone_id,
             const double t, const double dt, const double Istim,
             const double Ksac, Vector<double> &X, Vector<double> &Xg) const;

  /**
   * @brief Get the number of state variables.
   */
  unsigned int nX() const { return initial_X.size(); }

  /**
   * @brief Get the number of gating variables.
   */
  unsigned int nG() const { return initial_Xg.size(); }

  /**
   * @brief Get the index of the intracellular calcium concentration in the
   * state vector.
   *
   * This is the index of the variable to be used for electromechanics coupling.
   */
  virtual unsigned int get_calcium_index() const = 0;

  /**
   * @brief Get a list of state variables to export to VTU.
   *
   * By default, this function returns the calcium index as the only exported
   * state variable. Derived classes can override this to export additional
   * states if needed. Beware that, for phenomenological models, the calcium
   * variable might actually be a proxy for calcium, rather than the actual
   * concentration (and in particular it might be non-dimensional or have
   * different units than a molar concentration).
   *
   * @return A vector of pairs {variable_name, state_vector_index}. The vector
   * need not include the transmembrane potential V, which is exported by other
   * means.
   */
  virtual std::vector<std::pair<std::string, int>>
  get_output_variables() const {
    return {{"Calcium", get_calcium_index()}};
  }

  /**
   * @brief Get output variable information for output registration.
   */
  std::vector<outputType> get_registered_outputs() const;

protected:
  /**
   * @name Integration methods.
   * @{
   */

  /**
   * @brief Integrate the model with the forward Euler method.
   *
   * @f[ \begin{aligned}
   *   \mathbf{x}^{n+1} &=
   *     \mathbf{x}^n + \Delta t \mathbf{f}(\mathbf{x}^n, \mathbf{x}_g^n) \\
   *   \mathbf{x}_g^{n+1} &=
   *     \texttt{update_g}(\Delta t, \mathbf{x}^n, \mathbf{x}_g^n)
   * \end{aligned} @f]
   *
   * @param[in] zone_id Identifier for the transmural zone (epicardium,
   *   endocardium, myocardium).
   * @param[in,out] X Vector of state variables to be updated.
   * @param[in,out] Xg Vector of gating variables to be updated.
   * @param[in] Ts Current time.
   * @param[in] Ti Time step.
   * @param[in] Istim Applied current.
   * @param[in] Ksac Stretch-activated current coefficient.
   */
  void integ_fe(const unsigned int zone_id, Vector<double> &X,
                Vector<double> &Xg, const double Ts, const double Ti,
                const double Istim, const double Ksac) const;

  /**
   * @brief Integrate the model with the fourth-order explicit Runge-Kutta
   * method.
   *
   * @f[ \begin{aligned}
   *   \mathbf{x}^{(1)} &= \mathbf{x}^n \\
   *   \mathbf{f}^{(1)} &= \mathbf{f}(\mathbf{x}^{(1)}, \mathbf{x}_g^n) \\
   *   \mathbf{x}_g^{(1)} &=
   *     \texttt{update_g}(\Delta t/2, \mathbf{x}^{(1)}, \mathbf{x}_g^n) \\
   *   \\
   *   \mathbf{x}^{(2)} &= \mathbf{x}^n + \frac{\Delta t}{2}\mathbf{f}^{(1)} \\
   *   \mathbf{f}^{(2)} &= \mathbf{f}(\mathbf{x}^{(2)}, \mathbf{x}_g^{(1)}) \\
   *   \\
   *   \mathbf{x}^{(3)} &= \mathbf{x}^n + \frac{\Delta t}{2}\mathbf{f}^{(2)} \\
   *   \mathbf{f}^{(3)} &= \mathbf{f}(\mathbf{x}^{(3)}, \mathbf{x}_g^{(1)}) \\
   *   \\
   *   \mathbf{x}_g^{(4)} &=
   *     \texttt{update_g}(\Delta t, \mathbf{x}^{(1)}, \mathbf{x}_g^{(1)}) \\
   *   \mathbf{x}^{(4)} &= \mathbf{x} + \Delta t \mathbf{f}^{(3)} \\
   *   \mathbf{f}^{(4)} &= \mathbf{f}(\mathbf{x}^{(4)}, \mathbf{x}_g^{(4)}) \\
   *   \\
   *   \mathbf{x}^{n+1} &= \mathbf{x} + \frac{\Delta t}{6} \left(
   *     \mathbf{f}^{(1)} + 2\mathbf{f}^{(2)} +
   *     2\mathbf{f}^{(3)} + \mathbf{f}^{(4)} \right) \\
   *   \mathbf{x}_g^{n+1} &= \mathbf{x}_g^{(4)}
   * \end{aligned} @f]
   *
   * @param[in] zone_id Identifier for the transmural zone (epicardium,
   *   endocardium, myocardium).
   * @param[in,out] X Vector of state variables to be updated.
   * @param[in,out] Xg Vector of gating variables to be updated.
   * @param[in] Ts Current time.
   * @param[in] Ti Time step.
   * @param[in] Istim Applied current.
   * @param[in] Ksac Stretch-activated current coefficient.
   */
  void integ_rk(const unsigned int zone_id, Vector<double> &X,
                Vector<double> &Xg, const double Ts, const double Ti,
                const double Istim, const double Ksac) const;

  /**
   * @brief Integrate the model with the Crank-Nicolson method.
   *
   * The state variables @f$\mathbf{x}@f$ are updated by solving the following
   * non-linear problem:
   * @f[
   *   \frac{\mathbf{x}^{n+1} - \mathbf{x}^n}{\Delta t} = \frac{1}{2} \left(
   *     \mathbf{f}(\mathbf{x}^n, \mathbf{x}_g^n) +
   *     \mathbf{f}(\mathbf{x}^{n+1}, \mathbf{x}_g^n) \right)
   * @f]
   * The problem is solved through Newton's method. Subsequently, the gating
   * variables are updated:
   * @f[
   *   \mathbf{x}_g^{n+1} =
   *     \texttt{update_g}(\Delta t, \mathbf{x}^{n+1}, \mathbf{x}_g^n)
   * @f]
   *
   * Since this method involves solving a nonlinear system of equations, it is
   * only available for those derived classes that implement the Jacobian matrix
   * of the system through getj. An exception will be raised otherwise.
   *
   * @param[in] zone_id Identifier for the transmural zone (epicardium,
   *   endocardium, myocardium).
   * @param[in,out] X Vector of state variables to be updated.
   * @param[in,out] Xg Vector of gating variables to be updated.
   * @param[in] Ts Current time.
   * @param[in] Ti Time step.
   * @param[in] Istim Applied current.
   * @param[in] Ksac Stretch-activated current coefficient.
   * @param[in] max_iter Maximum number of Newton iterations. Beware that no
   *   error is raised if the maximum number of iterations is reached, so that
   *   the solution might be affected by the truncation of the iterations.
   * @param[in] rtol Relative tolerance for the Newton method.
   * @param[in] atol Absolute tolerance for the Newton method.
   */
  void integ_cn2(const unsigned int zone_id, Vector<double> &X,
                 Vector<double> &Xg, const double Ts, const double Ti,
                 const double Istim, const double Ksac,
                 const unsigned int max_iter, const double rtol,
                 const double atol) const;

  /**
   * @}
   */

  /**
   * @brief Update variables with analytical solution.
   *
   * @param[in] zone_id Identifier for the transmural zone (epicardium,
   *   endocardium, myocardium).
   * @param[in] dt Time step.
   * @param[in] X Vector of state variables.
   * @param[out] Xg Vector of gating variables.
   */
  virtual void update_g(const unsigned int zone_id, const double dt,
                        const Vector<double> &X, Vector<double> &Xg) const = 0;

  /**
   * @brief Model right hand side.
   *
   * Defines the right-hand side function for the potential and ionic equations.
   * Must be ovverridden in derived classes.
   *
   * @param[in] zone_id Identifier for the transmural zone (epicardium,
   *   endocardium, myocardium).
   * @param[in] X Vector of state variables.
   * @param[in] Xg Vector of gating variables.
   * @param[in] I_stim Applied current.
   * @param[in] I_sac Stretch-activated current.
   *
   * @return A vector containing the right-hand side of the model equations.
   */
  virtual Vector<double> getf(const unsigned int zone_id,
                              const Vector<double> &X, const Vector<double> &Xg,
                              const double I_stim,
                              const double I_sac) const = 0;

  /**
   * @brief Model jacobian.
   *
   * Defines the jacobian matrix of the model equations, that is the matirx of
   * derivatives of the function evaulated by getf.
   *
   * @param[in] zone_id Identifier for the transmural zone (epicardium,
   *   endocardium, myocardium).
   * @param[in] X Vector of state variables.
   * @param[in] Xg Vector of gating variables.
   * @param[in] Ksac Stretch-activated current coefficient.
   *
   * @return A matrix containing the jacobian of the model equations.
   */
  virtual Array<double> getj(const unsigned int zone_id,
                             const Vector<double> &X, const Vector<double> &Xg,
                             const double Ksac) const {
    svmp::raise<svmp::FE::NotImplementedException>(
        SVMP_HERE, "getj method not implemented for this ionic model.");

    // Dummy return statement to avoid compiler warnings.
    Array<double> dummy(X.size(), X.size());
    return dummy;
  }

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
   * @{
   */

  /// Voltage scaling [mV].
  const double Vscale;

  /// Time scaling [ms].
  const double Tscale;

  /// Voltage offset parameter [mV].
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