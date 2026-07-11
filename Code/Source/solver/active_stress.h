// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_H
#define ACTIVE_STRESS_H

#include "Array.h"
#include "Parameters.h"
#include "Vector.h"
#include "consts.h"
#include "factory.h"

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
 * Accordingly, this class works nodally, by evaluating the active tension at
 * every mesh node and storing it in a vector, whose values can be accessed
 * through @ref ActiveStress::get_tension_fibers.
 *
 * ### Directional distribution of active stress
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
 * and provides the functions @ref ActiveStress::get_tension_fibers,
 * @ref ActiveStress::get_tension_sheets and @ref
 * ActiveStress::get_tension_sheet_normals to access @f$\eta_f \Tact@f$,
 * @f$\eta_s \Tact@f$ and @f$\eta_n \Tact@f$, respectively.
 *
 * ### Implementing concrete active stress models
 *
 * To implement a new active stress model, the following steps need to be taken:
 *
 * 1. Create a new class derived from @ref ActiveStress.
 * 2. Override the methods @ref init_local, @ref advance_time_step_local and
 *    @ref compute_active_tension_local, defining the initial condition,
 *    time evolution and active tension computation, respectively, for a single
 *    node.
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
 */
class ActiveStress {
public:
  /**
   * @brief Constructor.
   *
   * @param n_states_ Number of state variables for this model.
   */
  ActiveStress(const unsigned int n_states_) : n_states(n_states_) {}

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
   * @brief Get the tension along fibers @f$\eta_f \Tact@f$ at a given point.
   */
  double get_tension_fibers(const int idx) const {
    return eta_f * active_tension[idx];
  }

  /**
   * @brief Get the tension along sheets @f$\eta_s \Tact@f$ at a given point.
   */
  double get_tension_sheets(const int idx) const {
    return eta_s * active_tension[idx];
  }

  /**
   * @brief Get the tension along sheet normals @f$\eta_n \Tact@f$ at a given
   * point.
   */
  double get_tension_sheet_normals(const int idx) const {
    return eta_n * active_tension[idx];
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
   * @brief Advance in time.
   *
   * @param[in] t Current time (i.e. the time instant being advanced to).
   * @param[in] dt Time step size.
   * @param[in] calcium Calcium concentration at every node.
   * @param[in] fiber_stretch Fiber stretch at every node. This is usually
   *   computed with post::fib_stretch.
   * @param[in] fiber_stretch_rate Fiber stretch rate at every node. This is
   *   usually computed with post::fib_stretch_rate.
   */
  virtual void advance_time_step(const double t, const double dt,
                                 const Vector<double> &calcium,
                                 const Vector<double> &fiber_stretch,
                                 const Vector<double> &fiber_stretch_rate);

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
   */
  virtual double
  compute_active_tension_local(const Vector<double> &state) const = 0;

  /// Current time. Updated whenever calling @ref advance_time_step.
  double time;

  /// State variables for the model.
  Array<double> states;

  /// Active tension at every node.
  Vector<double> active_tension;

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