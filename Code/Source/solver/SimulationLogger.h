// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SIMULATION_LOGGER_H 
#define SIMULATION_LOGGER_H 

#include <fstream>
#include <string>

/// @brief The SimulationLogger class is used to write information 
/// to a text file and optionally to cout.
//
class SimulationLogger {

  public:
    SimulationLogger() { }

    SimulationLogger(const std::string& file_name, bool cout_write=false)
    { 
      this->initialize(file_name, cout_write);
    }

    void initialize(const std::string& file_name, bool cout_write=false) 
    {
      log_file_.open(file_name);
      if (log_file_.fail()) {
        throw std::runtime_error("[SimulationLogger] Unable to open the file '" + file_name + "' for writing.");
      }

      cout_write_ = cout_write;
      file_name_ = file_name;
    }

    ~SimulationLogger() 
    {
      log_file_.close();
    }

  template <class T> SimulationLogger& operator<< (const T& value)
  {
    if (file_name_ == "") { 
      return *this;
    }

    log_file_ << value;

    if (cout_write_) {
      std::cout << value;
    }

    return *this;
  }

  SimulationLogger& operator<<(std::ostream&(*value)(std::ostream& o))
  {
    if (file_name_ == "") { 
      return *this;
    }

    log_file_ << value;

    if (cout_write_) {
      std::cout << value;
    }

    return *this;
  };

  private:
    bool cout_write_ = false;
    std::string file_name_;
    std::ofstream log_file_;
};

#endif


