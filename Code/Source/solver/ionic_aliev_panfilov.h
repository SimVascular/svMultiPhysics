// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_ALIEV_PANFILOV_H
#define IONIC_ALIEV_PANFILOV_H

#include "ionic_model.h"

#include "Vector.h"

/**
 * @brief Aliev-Panfilov ionic model.
 *
 * This class implements the split-parameter Aliev-Panfilov formulation of
 * Goktepe and Kuhl [2], Eqs. 19--20, derived from the original
 * Aliev-Panfilov model [1]. In [1], the single parameter @f$a@f$ appears in
 * both the excitation cubic and the recovery equation. In [2], these roles
 * are split into @f$\alpha@f$ and @f$b@f$, respectively, while @f$k@f$ is
 * renamed @f$c@f$ and @f$\varepsilon_0@f$ is renamed @f$\gamma@f$.
 *
 * The svMultiPhysics parameter mapping to [2] is:
 * @f$\texttt{alpha}\leftrightarrow\alpha@f$,
 * @f$\texttt{b}\leftrightarrow b@f$,
 * @f$\texttt{a}\leftrightarrow\gamma@f$,
 * @f$\texttt{c}\leftrightarrow c@f$, and
 * @f$\texttt{mu1},\texttt{mu2}\leftrightarrow\mu_1,\mu_2@f$.
 * The parameter defaults, including @f$\texttt{alpha}=0.01@f$, and the
 * initial recovery value are svMultiPhysics defaults; they are not presented
 * as the parameter set of either source. Applied and stretch-activated
 * currents are added using the svMultiPhysics sign convention.
 *
 * **References**:
 * 1. [Aliev, Panfilov (1996)](https://doi.org/10.1016/0960-0779(95)00089-5)
 * 2. [Goktepe, Kuhl (2009)](https://doi.org/10.1002/nme.2571)
 */
class AlievPanfilov : public IonicModel {
public:
  /// Model label.
  static inline const std::string label = "AP";

  /// State variables.
  static inline const InitialStates initial_X = {{"V", -80.0}, {"w", 1.0e-3}};

  /// Gating variables.
  static inline const InitialStates initial_Xg = {};

  /// Index of the recovery variable (w), used as calcium proxy for
  /// electromechanical coupling.
  static constexpr unsigned int calcium_index = 1;

  /// Model parameters class.
  class Parameters : public IonicModelParameters {
  public:
    Parameters() : IonicModelParameters(label, initial_X, initial_Xg) {
      constexpr bool required = true;

      add_parameter("alpha", 1.0e-2, required);
      add_parameter("a", 2.0e-3, required);
      add_parameter("b", 0.15, required);
      add_parameter("c", 8.0, required);
      add_parameter("mu1", 0.20, required);
      add_parameter("mu2", 0.30, required);
    }
  };

  /// Constructor.
  AlievPanfilov()
      : IonicModel(initial_X, initial_Xg,
                   /* Vrest_ = */ -80.0, /* Vscale_ = */ 100.0,
                   /* Tscale_ = */ 12.90, /* Voffset_ = */ -80.0) {}

  /// Construct an instance of model parameters.
  virtual std::unique_ptr<IonicModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

  /// Read model parameters from a parameter object.
  virtual void read_parameters(const IonicModelParameters &params) override;

  /// Distribute model parameters to all parallel processes.
  virtual void distribute_parameters(const CmMod &cm_mod,
                                     const cmType &cm) override;

  /// Get the calcium proxy index.
  virtual unsigned int get_calcium_index() const override {
    return calcium_index;
  }

protected:
  /// @name Model parameters
  /// @{

  /// Excitation threshold @f$\alpha@f$ in Goktepe-Kuhl [2]. Together with
  /// @c b, this splits the two roles of @f$a@f$ in Aliev-Panfilov [1].
  double alpha = 1.0e-2;

  /// Baseline recovery rate @f$\gamma@f$ in Goktepe-Kuhl [2], corresponding
  /// to @f$\varepsilon_0@f$ in Aliev-Panfilov [1].
  double a = 2.0e-3;

  /// Recovery-equation shift @f$b@f$ in Goktepe-Kuhl [2]. Together with
  /// @c alpha, this splits the two roles of @f$a@f$ in Aliev-Panfilov [1].
  double b = 0.15;

  /// Kinetic coefficient @f$c@f$ in Goktepe-Kuhl [2], corresponding to
  /// @f$k@f$ in Aliev-Panfilov [1].
  double c = 8.0;

  double mu1 = 0.20; ///< @f$\mu_1@f$ in [1, 2].
  double mu2 = 0.30; ///< @f$\mu_2@f$ in [1, 2].

  /// @}

  /// Update variable with analytical solution. This model has none, so this
  /// method does nothing.
  virtual void update_g(const unsigned int zone_id, const double dt,
                        const Vector<double> &X,
                        Vector<double> &Xg) const override {}

  /// Model right-hand side.
  virtual Vector<double> getf(const unsigned int zone_id,
                              const Vector<double> &X, const Vector<double> &Xg,
                              const double I_stim,
                              const double I_sac) const override;

  /// Model jacobian.
  virtual Array<double> getj(const unsigned int zone_id,
                             const Vector<double> &X, const Vector<double> &Xg,
                             const double Ksac) const override;
};

#endif
