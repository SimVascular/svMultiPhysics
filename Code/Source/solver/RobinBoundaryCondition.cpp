#include "RobinBoundaryCondition.h"

#define debug_robin_bc

void RobinBoundaryCondition::init_from_vtp(const std::string& vtp_file_path, const faceType& face) {

    // Set array names specific to Robin BC
    std::vector<std::string> array_names{"Stiffness", "Damping"};

    // Call the base class method
    BoundaryCondition::init_from_vtp(vtp_file_path, array_names, face);
}

void RobinBoundaryCondition::init_uniform(double uniform_stiffness, double uniform_damping, const faceType& face) {

    // Set uniform values with names specific to Robin BC
    std::map<std::string, double> uniform_values = {
        {"Stiffness", uniform_stiffness},
        {"Damping", uniform_damping}
    };

    // Call the base class method
    BoundaryCondition::init_uniform(uniform_values, face);
}
