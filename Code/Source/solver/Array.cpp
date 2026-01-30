// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Array.h"

// b o o l //

template<>
int Array<bool>::id = 0;

template<>
double Array<bool>::memory_in_use = 0;

template<>
double Array<bool>::memory_returned = 0;

template<>
int Array<bool>::num_allocated = 0;

template<>
int Array<bool>::active = 0;

template<>
bool Array<bool>::write_enabled = false;

// d o u b l e //

template<>
int Array<double>::id = 0;

template<>
double Array<double>::memory_in_use = 0;

template<>
double Array<double>::memory_returned = 0;

template<>
int Array<double>::num_allocated = 0;

template<>
int Array<double>::active = 0;

template<>
void Array<double>::memory(const std::string& prefix)
{
  utils::print_mem("Array<double>", prefix, memory_in_use, memory_returned);
}

template<>
void Array<double>::stats(const std::string& prefix)
{ 
  utils::print_stats("Array<double>", prefix, num_allocated, active);
}

template<>
bool Array<double>::write_enabled = false;

//  i n t  //

template<>
int Array<int>::id = 0;

template<>
double Array<int>::memory_in_use = 0;

template<>
double Array<int>::memory_returned = 0;

template<>
int Array<int>::num_allocated = 0;

template<>
int Array<int>::active = 0;

template<>
void Array<int>::memory(const std::string& prefix)
{
  utils::print_mem("Array<int>", prefix, memory_in_use, memory_returned);
}

template<>
void Array<int>::stats(const std::string& prefix)
{
  utils::print_stats("Array<int>", prefix, num_allocated, active);
}

template<>
bool Array<int>::write_enabled = false;

