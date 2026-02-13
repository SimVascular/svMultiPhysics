// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SIMULATION_H 
#define SIMULATION_H 

#include "ComMod.h"
#include "Parameters.h"
#include "SimulationLogger.h"
#include "LinearAlgebra.h"

#include <string>
#include <memory>

// Forward declaration
class Integrator;

class Simulation {

  public:
    Simulation();
    ~Simulation();

    const mshType& get_msh(const std::string& name);

    CepMod& get_cep_mod() { return cep_mod; };
    ChnlMod& get_chnl_mod() { return chnl_mod; };
    ComMod& get_com_mod() { return com_mod; };
    Integrator& get_integrator();

    // Initialize the Integrator object after simulation setup is complete
    // Takes ownership of solution states via move semantics
    void initialize_integrator(SolutionStates&& solutions);

    // Read a solver paramerer input XML file.
    void read_parameters(const std::string& fileName);

    // Set simulation and module member data from Parameters.
    void set_module_parameters();

    //----- Fortran subroutines -----//
    //void read_msh();

    //----- Fortran modules -----//
    CepMod cep_mod;
    ChnlMod chnl_mod;
    CmMod cm_mod;
    ComMod com_mod;

    // Solver parameters read in from solver input XML file.
    Parameters parameters;

    // Log solution information.
    SimulationLogger logger;

    // Number of time steps
    int nTs;

    // Simulation initialization file path
    std::string fTmp;

    // Spectral radius of infinite time step; this is later used in equations.
    double roInf;

    // Simulation requires remeshing

    bool isReqd;

    // Name of the history file.
    std::string history_file_name;

    LinearAlgebra* linear_algebra = nullptr;

  private:
    // Time integrator for Newton iteration loop
    std::unique_ptr<Integrator> integrator_;
};

#endif

