/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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

