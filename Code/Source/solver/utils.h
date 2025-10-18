// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef UTILS_H 
#define UTILS_H 

#include "Array.h"

namespace utils {

class stackType
{
  public:
    /// Maximum length of the stack
    int maxN = 0;

    /// Current size of stack
    int n = 0;

    /// Values inside stack
    Vector<int> v;
};

class queueType
{
  public:
    int n = 0;
    int maxN = 0;
    Vector<int> v;
};


bool btest(int value, int pos);

int CountBits(int n);

double cput();

Vector<double> cross(const Array<double>& V);

bool dequeue(queueType& que, int& iVal);
void enqueue(queueType& que, int iVal);

int ibclr(int value, int pos);
int ibset(int value, int pos);
bool is_zero(double value1, double value2=0.0);

double mem_usage(const bool print_usage=false, const std::string& prefix="");

double norm(const Vector<double>& U); 
double norm(const Vector<double>& U, const Vector<double>& V);
double norm(const Array<double>& U);

void print_mem(const std::string& type, const std::string& prefix, const double memory_in_use, const double memory_returned);
void print_stats(const std::string& type, const std::string& prefix, const int allocated, const int active);

bool pull_stack(stackType& stk, int& iVal);
void push_stack(stackType& stk, int iVal);
void push_stack(stackType& stk, std::initializer_list<int> values);

int sign(double value);

void swap(int& value1, int& value2);

void find_loc(const Array<int>& array, int value, std::array<int, 2>& ind);

};

#endif

