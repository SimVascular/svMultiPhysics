// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Array.h"

#ifndef ARRAY3_H 
#define ARRAY3_H 

#include <float.h>
#include <iostream>
#include <string>
#include <cstring>

#ifdef ENABLE_ARRAY_INDEX_CHECKING
#define Array3_check_enabled
#endif

template<typename T>

/*! @brief The Array3 template class implements a simple interface to 3D arrays.
*
* A 3D array is defined by (num_rows * num_cols) * num_slices values.
*/
class Array3 
{
  public:

    static int num_allocated;
    static int active;
    static double memory_in_use;
    static double memory_returned;
    static bool write_enabled;
    static void memory(const std::string& prefix="");
    static void stats(const std::string& prefix="");

    Array3() 
    {
      nrows_ = 0;
      ncols_ = 0;
      nslices_ = 0;
      slice_size_ = 0;
      size_ = 0;
      data_ = nullptr;
      num_allocated += 1;
      active += 1;
    };

    Array3(const int num_rows, const int num_cols, const int num_slices, bool row_major=true) 
    {
      // [NOTE] This is tfu but need to mimic fortran.
      nrows_ = num_rows;
      ncols_ = num_cols;
      nslices_ = num_slices;

      if ((num_rows <= 0) || (num_cols <= 0) || (num_slices <= 0)) { 
        return; 
        //throw std::runtime_error(+"Allocating zero size Array3.");
      }

      allocate(num_rows, num_cols, num_slices, row_major);
      num_allocated += 1;
      active += 1;
    }

    /// @brief Array copy
    Array3(const Array3 &rhs)
    {
      if ((rhs.nrows_ <= 0) || (rhs.ncols_ <= 0) || (rhs.nslices_ <= 0)) {
        return;
      }

      allocate(rhs.nrows_, rhs.ncols_, rhs.nslices_);
      memcpy(data_, rhs.data_, size_*sizeof(T));
      num_allocated += 1;
      active += 1;
    }

    ~Array3() 
    {
      if (data_ != nullptr) {
        memory_in_use -= sizeof(T) * size_;;
        memory_returned += sizeof(T) * size_;;
        active -= 1;
        delete [] data_;
        data_ = nullptr;
       }
     }

    int ncols() const { return ncols_; }
    int nrows() const { return nrows_; }
    int nslices() const { return nslices_; }

    void allocate(const int num_rows, const int num_cols, const int num_slices, const bool row_major=true)
    {
      nrows_ = num_rows;
      ncols_ = num_cols;
      nslices_ = num_slices;
      slice_size_ = ncols_ * nrows_;
      size_ = nrows_ * ncols_ * nslices_;
      data_ = new T [size_];
      memset(data_, 0, sizeof(T)*size_);
      memory_in_use += sizeof(T) * size_;;
    }

    void check_index(const int i, const int j, const int k) const
    {
      if (data_ == nullptr) {
        throw std::runtime_error(+"Accessing null data in Array3.");
      } 

      if ((i < 0) || (i >= nrows_) or (j < 0) || (j >= ncols_) or (k < 0) || (k >= nslices_)) {
        auto i_str = std::to_string(nrows_);
        auto j_str = std::to_string(ncols_);
        auto k_str = std::to_string(nslices_);
        auto dims = i_str + " x " + j_str + " x " + k_str;
        auto index_str = " " + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + " ";
        throw std::runtime_error("Index (i,j,k)=" + index_str + " is out of bounds for " + dims + " array.");
      }
    }

    friend std::ostream& operator << (std::ostream& out, const Array3<T>& lhs)
    {
      if (lhs.data_ == nullptr) {
        throw std::runtime_error("[Array3] Accessing null data in ostream.");
      } 

      for (int i = 0; i < lhs.size(); i++) {
        out << lhs.data_[i];
        if (i != lhs.size()-1) {
          out << ", ";
         }
       }
       return out;
    }

    /// @brief Free the array data 
    ///
    /// This is to replicate the Fortran DEALLOCATE().
    void clear()
    {
      if (data_ != nullptr) {
        delete [] data_;
        memory_in_use -= sizeof(T) * size_;;
        memory_returned += sizeof(T) * size_;;
      }

      nrows_ = 0;
      ncols_ = 0;
      nslices_ = 0;
      slice_size_ = 0;
      size_ = 0;
      data_ = nullptr;
    }

    /// @brief rslice
    ///
    /// Return an Array with data pointing into the Array3 internal data.
    Array<T> rslice(const int slice) const
    { 
      #ifdef Array3_check_enabled
      check_index(0, 0, slice);
      #endif
      
      Array<T> array_slice(nrows_, ncols_, &data_[slice*slice_size_]);
      
      return array_slice;
    }

    T* slice_data(const int slice) { 
      return &data_[slice*slice_size_];
    }

    void print(const std::string& label)
    {
      printf("%s (%d): \n", label.c_str(), size_);
      for (int i = 0; i < size_; i++) {
        if (data_[i] != 0.0) {
          printf("%s %d %g\n", label.c_str(), i+1, data_[i]);
        }
      }
    }

    /// @brief Get a slice.
    Array<T> slice(const int slice) const
    {
      #ifdef Array3_check_enabled
      check_index(0, 0, slice);
      #endif
      Array<T> array_slice(nrows_, ncols_);

      for (int col = 0; col  < ncols_; col++) {
        for (int row = 0; row  < nrows_; row++) {
          array_slice(row, col) = data_[row + col*nrows_ + slice*slice_size_];
        }
      }

      return array_slice;
    }

    void set_row(const int col, const int slice, const std::vector<T>& values) const
    {
      #ifdef Array3_check_enabled
      check_index(0, col, slice);
      #endif
      for (int row = 0; row < values.size(); row++) {
        data_[row + col*nrows_ + slice*slice_size_] = values[row];
      }
    }

    void set_slice(const int slice, const Array<T>& values) const
    {
      #ifdef Array3_check_enabled
      check_index(0, 0, slice);
      #endif
      for (int col = 0; col  < ncols_; col++) {
        for (int row = 0; row  < nrows_; row++) {
          data_[row + col*nrows_ + slice*slice_size_] = values(row,col);
        }
      }
    }

    int size() const
    {
      return size_;
    }

    /// @brief  Test if an array has rows or columns or slices set but no data.
    ///
    /// [NOTE] This is tfu but need to mimic fortran.
    bool allocated()
    {
      if ((nrows_ > 0) || (ncols_ > 0) || nslices_ > 0) {
        return true;
      }

      return false;
    }

    int msize() const
    {
      return size_ * sizeof(T);
    }

    /// @brief Resize the array.
    void resize(const int num_rows, const int num_cols, const int num_slices)
    {
      // [NOTE] This is tfu but need to mimic fortran.
      nrows_ = num_rows;
      ncols_ = num_cols;
      nslices_ = num_slices;

      if ((num_rows <= 0) || (num_cols <= 0) || (num_slices <= 0)) { 
        return; 
      }

      if (data_ != nullptr) {
        delete [] data_;
        memory_in_use -= sizeof(T) * size_;;
        memory_returned += sizeof(T) * size_;;
        data_ = nullptr;
      }

      allocate(num_rows, num_cols, num_slices);
    }

    /// @brief Set the array values from an Array with the equivalent number of values.
    ///
    /// This sort of replicates the Fortan reshape function.
    void set_values(Array<T>& rhs)
    {
      int rhs_size = rhs.size();

      if (size_ != rhs_size) {
        throw std::runtime_error("The RHS size " + std::to_string(rhs_size) + " does not have the same size " +
            std::to_string(size_) + " of this array.");
      }

      auto rhs_data = rhs.data();

      for (int i = 0; i < size_; i++) {
        data_[i] = rhs_data[i];
      }
    }

    void read(const std::string& file_name) 
    { 
      auto fp = fopen(file_name.c_str(), "rb");
      fread(&size_, sizeof(int), 1, fp);
      fread(data_, sizeof(double), size_, fp);
      fclose(fp);
    }

    void write(const std::string& label, bool memory=true)
    {
      if (!write_enabled) {
        return;
      }

      auto file_prefix = build_file_prefix(label);
      auto file_name = file_prefix + "_cm.bin";

     // Write binary file.
      auto fp = fopen(file_name.c_str(), "wb");
      fwrite(&size_, sizeof(int), 1, fp);
      fwrite(data_, sizeof(double), size_, fp);
      fclose(fp);
    }

    void append(const std::string& label, bool memory=true)
    {
      if (!write_enabled) {
        return;
      }

      auto file_prefix = build_file_prefix(label);
      auto file_name = file_prefix + "_cm.bin";

      // Append to binary file.
      auto fp = fopen(file_name.c_str(), "ab");
      fwrite(data_, sizeof(double), size_, fp);
      fclose(fp);
    }

    /////////////////////////
    //  O p e r a t o r s  //
    /////////////////////////

    /// @brief Array assignment 
    ///
    /// Copy data to replicate Fortran behavior.
    Array3& operator = (const Array3& rhs)
    {
      if ((rhs.nrows_ <= 0) || (rhs.ncols_ <= 0) || (rhs.nslices_ <= 0)) { 
        return *this;
      }

      if (rhs.data_ == nullptr) { 
        throw std::runtime_error(+"RHS has null data.");
      }

      if (this == &rhs) {
        return *this;
      }

      if (size_ != rhs.size_) {
        delete [] data_;
        allocate(rhs.nrows_, rhs.ncols_, rhs.nslices_);
      }

      memcpy(data_, rhs.data_, sizeof(T) * size_);

      return *this;
    }

    /// @brief Get the array value at (row,col).
    const T& operator()(const int row, const int col, const int slice) const
    {
      #ifdef Array3_check_enabled
      check_index(row, col, slice);
      #endif
      return data_[row + col*nrows_ + slice*slice_size_];
    }

    /// @brief Set the array value at (row,col).
    T& operator()(const int row, const int col, const int slice)
    {
      #ifdef Array3_check_enabled
      check_index(row, col, slice);
      #endif
      return data_[row + col*nrows_ + slice*slice_size_];
    }

    Array3& operator=(const double value)
    {
      for (int i = 0; i < size_; i++) {
        data_[i] = value;
      }
      return *this;
    }

    ///////////////////////////////////////////////////
    //  M a t h e m a t i c a l   o p e r a t o r s  //
    ///////////////////////////////////////////////////

    /// @brief Multiply by a scalar.
    Array3<T> operator * (const T value) const
    {
      Array3<T> result(nrows_, ncols_, nslices_);
      for (int i = 0; i < size_; i++) {
        result.data_[i] = value * data_[i];
      }
      return result;
    }

    Array3<T>& operator *= (const T value) 
    {
      for (int i = 0; i < size_; i++) {
        data_[i] *= value;
      }
      return *this;
    }

    friend const Array3<T> operator * (const T value, const Array3& rhs)
    {
      if (rhs.data_ == nullptr) { 
        throw std::runtime_error("Null data for rhs Array3.");
      }
      Array3<T> result(rhs.nrows_, rhs.ncols_, rhs.nslices_);
      for (int i = 0; i < rhs.size_; i++) {
        result.data_[i] = value * rhs.data_[i];
      }
      return result;
    }

   T* data() const { 
      return data_;
    }

  private:
    int nrows_ = 0;
    int ncols_ = 0;
    int nslices_ = 0;
    int slice_size_ = 0;
    int size_ = 0;
    T *data_ = nullptr;

};

#endif

