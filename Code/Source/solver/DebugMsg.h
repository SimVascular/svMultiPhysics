// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DEBUG_MSG_H 
#define DEBUG_MSG_H 

#include <iostream>

/// @brief The DebugMsg is class is used to print debugging messages. 
class DebugMsg 
{
  public:
    DebugMsg() = delete;
    DebugMsg(const char* function, const int task_id, bool add_endl=true) 
    {
      function_name_ = function;
      task_id_ = task_id;
      prefix_ = "[" + function_name_ + ":" + std::to_string(task_id) + "] ";
      banner_ = prefix_ + "=============== " + function_name_ + " ===============";
      add_endl_ = add_endl;
    }

    void banner() { std::cout << banner_ << std::endl; };
    std::string prefix() { return prefix_; };

    template <class T> DebugMsg& operator<< (const T& x) 
    { 
      if (start_) {
        std::cout << prefix_; 
      }

      std::cout << x;
      count_ += 1;

      if (add_endl_ && count_ == 2) {
        std::cout << std::endl; 
        start_ = true;
        count_ = 0;
      } else {
        start_ = false;
      }

      return *this; 
    };

    DebugMsg& operator<<(std::ostream& (*f)(std::ostream& o)) 
    {
      std::cout << f;
      start_ = true;
      count_ = 0;
      return *this;
    };

  private:
    bool add_endl_ = true;
    std::string banner_;
    int count_ = 0;
    std::string function_name_;
    std::string prefix_;
    bool start_ = true;
    int task_id_;
};

#endif

