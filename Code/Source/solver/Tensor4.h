// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TENSOR4_H
#define TENSOR4_H

#include <cstring>
#include <iostream>

#include "Vector.h"

#ifdef ENABLE_ARRAY_INDEX_CHECKING
#define Tensor4_check_enabled
#endif

/// @brief The Tensor4 template class implements a simple interface to 4th order tensors.
//
template<typename T>
class Tensor4
{
  public:
    std::string name_ = "";
    int ni_ = 0;
    int nj_ = 0;
    int nk_ = 0;
    int nl_ = 0;
    int p1_ = 0;
    int p2_ = 0;
    int size_ = 0;
    T *data_ = nullptr;

    Tensor4() 
    {
      std::string name_ = "";
      int ni_ = 0;
      int nj_ = 0;
      int nk_ = 0;
      int nl_ = 0;
      int p1_ = 0;
      int p2_ = 0;
      int size_ = 0;
      T *data_ = nullptr;
    };

    Tensor4(const int num_i, const int num_j, const int num_k, const int num_l)
    {
      allocate(num_i, num_j, num_k, num_l);
    }

    ~Tensor4() 
    {
      //std::cout << "- - - - - Tensor4 dtor - - - - - " << std::endl;
      if (data_ != nullptr) {
        //std::cout << "[Tensor4 dtor] delete[] data: " << data_ << std::endl;
        delete [] data_;
        data_ = nullptr;
       }
     }

    // Copy
    Tensor4(const Tensor4& rhs) 
    {
      if (rhs.ni_ <= 0 || rhs.nj_ <= 0 || rhs.nk_ <= 0 || rhs.nl_ <= 0) {
        return;
      }
      allocate(rhs.ni_, rhs.nj_, rhs.nk_, rhs.nl_);
      memcpy(data_, rhs.data_, size_*sizeof(T));
    }

    int num_i() { return ni_; }
    int num_j() { return nj_; }
    int num_k() { return nk_; }
    int num_l() { return nl_; }

    // Tensor4 assignment.
    //
    Tensor4& operator=(const Tensor4& rhs)
    {
      if (this == &rhs) {
        return *this;
      }

      if (rhs.ni_ <= 0 || rhs.nj_ <= 0 || rhs.nk_ <= 0 || rhs.nl_ <= 0) {
        return *this;
      }

      if (ni_ != rhs.ni_ || nj_ != rhs.nj_ || nk_ != rhs.nk_ || nl_ != rhs.nl_ <= 0) {
        clear();
        //delete [] data_;
        //data_ = nullptr;
        allocate(rhs.ni_, rhs.nj_, rhs.nk_, rhs.nl_);
      }

      memcpy(data_, rhs.data_, sizeof(T) * size_);
      return *this;
    }

    Tensor4& operator=(const double value)
    {
      for (int i = 0; i < size_; i++) {
        data_[i] = value;
      }
      return *this;
    }

    void write(const std::string& label, bool memory=true, T offset={}) const
    {
      int n = 0;
      char file_prefix[1000];
      for (auto c : label) {
        if (c == '[') {
          continue;
        }
        if (c == ']') {
          file_prefix[n] = '_';
          n += 1;
          continue;
        }
        if (c == ' ') {
          continue;
        }
        if (c == ':') {
          c =  '_';
        }
        file_prefix[n] = c;
        n += 1;
      }
      file_prefix[n] = '\0';

      FILE* fp;
      std::string file_name;

      // Write binary file.
      file_name = std::string(file_prefix) + "_cm.bin";
      std::cout << "[Tensor4::write] file_name: " << file_name << std::endl;
      fp = fopen(file_name.c_str(), "wb");
      fwrite(&size_, sizeof(int), 1, fp);
      fwrite(data_, sizeof(double), size_, fp);
      fclose(fp);
    }

    friend std::ostream& operator << (std::ostream& out, const Tensor4<T>& lhs)
    {
      for (int i = 0; i < lhs.size_; i++) {
        out << lhs.data_[i];
        if (i != lhs.size_-1) {
          out << ", ";
         }
       }
       return out;
    }

    //----------
    // allocate
    //----------
    //
    void allocate(const int num_i, const int num_j, const int num_k, const int num_l)
    {
      //std::cout << "----- Tensor4::allocate -----" << std::endl;
      ni_ = num_i;
      nj_ = num_j;
      nk_ = num_k;
      nl_ = num_l;
      p1_ = num_i * num_j;
      p2_ = p1_ * num_l;
      size_ =  ni_ * nj_ * nk_ * nl_;
      data_ = new T [size_];
      memset(data_, 0, sizeof(T)*size_);
      //std::cout << "[Tensor4::allocate] data_: " << data_ << std::endl;
    }

    //-------------
    // check_index
    //-------------
    //
    void check_index(const int i, const int j, const int k, const int l) const
    {
      if (data_ == nullptr) {
        throw std::runtime_error(name_+"Accessing null data in Tensor4.");
      } 
      if ((i < 0) || (i >= ni_) || (j < 0) || (j >= nj_) || (k < 0) || (k >= nk_) || (l < 0) || (l >= nl_)) {
        auto i_str = std::to_string(ni_);
        auto j_str = std::to_string(nj_);
        auto k_str = std::to_string(nk_);
        auto l_str = std::to_string(nl_);
        auto dims = i_str + " x " + j_str + " x " + k_str + " x " + l_str;
        auto index_str = " " + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "," + std::to_string(l);
        throw std::runtime_error("Index (i,j,k,l)=" + index_str + " is out of bounds for " + dims + " array.");
      }
    }

    //-------
    // clear 
    //-------
    // Free the array data. 
    //
    // This is to replicate the Fortran DEALLOCATE().
    //
    void clear()
    {
      //std::cout << "----- Tensor4::erase -----" << std::endl;
      if (data_ != nullptr) {
        //std::cout << "[Tensor4::erase] data_: " << data_ << std::endl;
        delete [] data_;
      }

      ni_ = 0;
      nj_ = 0;
      nk_ = 0;
      nl_ = 0;
      p1_ = 0;
      p2_ = 0;
      size_ = 0;
      data_ = nullptr;
    }

    //--------
    // resize
    //--------
    // Resize the array.
    //
    void resize(const int num_i, const int num_j, const int num_k, const int num_l)
    {
      if (data_ != nullptr) {
        //std::cout << "[Tensor4::resize] data_: " << data_ << std::endl;
        delete [] data_;
        data_ = nullptr;
      }

      allocate(num_i, num_j, num_k, num_l);
    }

    /////////////////////////
    //  O p e r a t o r s  //
    /////////////////////////

    const T& operator()(const int i, const int j, const int k, const int l) const
    {
      #ifdef Tensor4_check_enabled
      check_index(i, j, k , l);
      #endif
      return data_[i + j*ni_ + p1_*k + p2_*l ];
    }

    T& operator()(const int i, const int j, const int k, const int l)
    {
      #ifdef Tensor4_check_enabled
      check_index(i, j, k , l);
      #endif
      return data_[i + j*ni_ + p1_*k + p2_*l ];
    }

    ///////////////////////////////////////////////////
    //  M a t h e m a t i c a l   o p e r a t o r s  //
    ///////////////////////////////////////////////////

    // Multiply by a scalar.
    //
    Tensor4<T> operator*(const T value) const
    {
      Tensor4<T> result(ni_, nj_, nk_, nl_);
      for (int i = 0; i < size_; i++) {
        result.data_[i] = value * data_[i];
      }
      return result;
    }

    friend const Tensor4<T> operator * (const T value, const Tensor4& rhs)
    {
      if (rhs.data_ == nullptr) {
        throw std::runtime_error("Null data for rhs Tensor4.");
      }
      Tensor4<T> result(rhs.ni_, rhs.nj_, rhs.nk_, rhs.nl_);
      for (int i = 0; i < rhs.size_; i++) {
        result.data_[i] = value * rhs.data_[i];
      }
      return result;
    }

    // Compound add assignment. 
    //
    Tensor4<T> operator+=(const Tensor4<T>& rhs) const
    {
      for (int i = 0; i < size_; i++) {
        data_[i] += rhs.data_[i];
      }
      return *this;
    }

    // Compound subtract assignment. 
    //
    Tensor4<T> operator-=(const Tensor4<T>& rhs) const
    { 
      for (int i = 0; i < size_; i++) {
        data_[i] -= rhs.data_[i];
      }
      return *this;
    }

    // Add and subtract arrays.
    //
    Tensor4<T> operator+(const Tensor4<T>& rhs) const
    {
      Tensor4<T> result(rhs.ni_, rhs.nj_, rhs.nk_, rhs.nl_);
      for (int i = 0; i < size_; i++) {
        result.data_[i] = data_[i] + rhs.data_[i];
      }
      return result;
    }

    Tensor4<T> operator-(const Tensor4<T>& rhs) const
    {
      Tensor4<T> result(rhs.ni_, rhs.nj_, rhs.nk_, rhs.nl_);
      for (int i = 0; i < size_; i++) {
        result.data_[i] = data_[i] - rhs.data_[i];
      }
      return result;
    }

};

#endif

