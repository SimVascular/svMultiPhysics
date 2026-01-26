// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Array.h"
#include "utils.h"

// Check if index checking warning should be shown (only once, only on rank 0)
bool show_index_check_warning() {
  static bool shown = false;
  if (shown) {
    return false;
  }
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
      return false;
    }
  }
  shown = true;
  return true;
}

template<>
double Vector<double>::memory_in_use = 0;

template<>
double Vector<double>::memory_returned = 0;

template<>
int Vector<double>::num_allocated = 0;

template<>
int Vector<double>::active = 0;

template<>
bool Vector<double>::write_enabled = false;

template<>
void Vector<double>::memory(const std::string& prefix)
{
  utils::print_mem("Vector<double>", prefix, memory_in_use, memory_returned);
}

template<>
void Vector<double>::stats(const std::string& prefix)
{
  utils::print_stats("Vector<double>", prefix, num_allocated, active);
}

//------//
//  int  //
//------//

template<>
double Vector<int>::memory_in_use = 0;

template<>
double Vector<int>::memory_returned = 0;

template<>
int Vector<int>::num_allocated = 0;

template<>
int Vector<int>::active = 0;

template<>
bool Vector<int>::write_enabled = false;

template<>
void Vector<int>::memory(const std::string& prefix)
{
  utils::print_mem("Vector<int>", prefix, memory_in_use, memory_returned);
}

template<>
void Vector<int>::stats(const std::string& prefix)
{
  utils::print_stats("Vector<int>", prefix, num_allocated, active);
}

// Vector<Vector<double>>

template<>
double Vector<Vector<double>>::memory_in_use = 0;

template<>
double Vector<Vector<double>>::memory_returned = 0;

template<>
int Vector<Vector<double>>::num_allocated = 0;

template<>
int Vector<Vector<double>>::active = 0;

template<>
bool Vector<Vector<double>>::write_enabled = false;

// float //

template<>
double Vector<float>::memory_in_use = 0;

template<>
double Vector<float>::memory_returned = 0;

template<>
int Vector<float>::num_allocated = 0;

template<>
int Vector<float>::active = 0;

template<>
bool Vector<float>::write_enabled = false;

/// @brief Build a prefix for a file name from label which may be from a debugging message.
//
std::string build_file_prefix(const std::string& label)
{
  std::string file_prefix;

  for (auto c : label) {
    if (c == '[') {
      continue;
    }
    if (c == ']') {
      file_prefix.push_back('_');
      continue;
    }
    if (c == ' ') {
      continue;
    }
    if (c == ':') {
      c =  '_';
    }
    file_prefix.push_back(c);
  }

  return file_prefix;
}

