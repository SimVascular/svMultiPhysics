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

    void initialize(const std::string& file_name, bool cout_write=false) const
    {
      log_file_.open(file_name);
      if (log_file_.fail()) {
        throw std::runtime_error("[SimulationLogger] Unable to open the file '" + file_name + "' for writing.");
      }

      cout_write_ = cout_write;
      file_name_ = file_name;
    }

    bool is_initialized() const { return log_file_.is_open(); }

    ~SimulationLogger() 
    {
      log_file_.close();
    }

  template <class T> const SimulationLogger& operator<< (const T& value) const
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

  // Special handling for vector<string>
  const SimulationLogger& operator<< (const std::vector<std::string>& values) const
  {
    if (file_name_ == "") { 
      return *this;
    }

    bool first = true;
    for (const auto& value : values) {
      if (!first) {
        log_file_ << ", ";
        if (cout_write_) std::cout << ", ";
      }
      log_file_ << value;
      if (cout_write_) std::cout << value;
      first = false;
    }

    return *this;
  }

  const SimulationLogger& operator<<(std::ostream&(*value)(std::ostream& o)) const
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

  /// @brief Log a message with automatic space separation
  /// @param args Arguments to log
  template<typename... Args>
  const SimulationLogger& log_message(const Args&... args) const {
    if (file_name_ == "") return *this;

    bool first = true;
    ((first ? (first = false, (*this << args)) : (*this << ' ' << args)), ...);
    *this << std::endl;
    return *this;
  }

  private:
    // These members are marked mutable because they can be modified in const member functions.
    // While logging writes to a file (modifying these members), it doesn't change the logical
    // state of the logger - this is exactly what mutable is designed for: implementation
    // details that need to change even when the object's public state remains constant.
    mutable bool cout_write_ = false;
    mutable bool initialized_ = false;
    mutable std::string file_name_;
    mutable std::ofstream log_file_;
};

#endif


