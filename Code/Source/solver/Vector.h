// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef VECTOR_H 
#define VECTOR_H 

#include <algorithm>
#include <float.h>
#include <iostream>
#include <string>
#include <vector>

std::string build_file_prefix(const std::string& label);

#ifdef ENABLE_ARRAY_INDEX_CHECKING
#define Vector_check_enabled
#endif

/// @brief The Vector template class is used for storing int and double data.
//
template<typename T>
class Vector 
{
  public:

    static int num_allocated;
    static int active;
    static double memory_in_use;
    static double memory_returned;
    static bool write_enabled;
    static void memory(const std::string& prefix="");
    static void stats(const std::string& prefix="");

    Vector() 
    {
      check_type();
      size_ = 0;
      data_ = nullptr;
      num_allocated += 1;
      active += 1;
    };

    Vector(const int size)
    {
      // [NOTE] This is tfu but need to mimic Fortran that can allocate 0-sized arrays.
      if (size <= 0) {
        is_allocated_ = true;
        return;
      }
      check_type();
      allocate(size);
      num_allocated += 1;
      active += 1;
    }

    Vector(const int size, T* data)
    {
      reference_data_ = true;
      size_ = size;
      data_ = data;
    }

    Vector(std::initializer_list<T> values)
    {
      if (values.size() == 0) { 
        return;
      }
      check_type();
      allocate(values.size());
      std::copy(values.begin(), values.end(), data_);
      num_allocated += 1;
      active += 1;
    }

    ~Vector()
    {
      if (data_ != nullptr) {
        if (!reference_data_) { 
          delete[] data_; 
        }
        memory_in_use -= sizeof(T)*size_;
        memory_returned += sizeof(T)*size_;
        active -= 1;
      } 

      size_ = 0;
      data_ = nullptr;
    }

    // Vector copy.
    Vector(const Vector& rhs)
    {
      if (rhs.size_ <= 0) {
        return;
      }
      allocate(rhs.size_);
      for (int i = 0; i < rhs.size_; i++) {
        data_[i] = rhs.data_[i];
      }
      num_allocated += 1;
      active += 1;
    }
  
    bool allocated() const
    {
      return is_allocated_ ;
    }

    /// @brief Free the array data. 
    ///
    /// This is to replicate the Fortran DEALLOCATE().
    //
    void clear()
    {
      if (data_ != nullptr) {
        if (reference_data_) { 
          throw std::runtime_error("[Vector] Can't clear a Vector with reference data.");
        }
        delete [] data_;
        memory_in_use -= sizeof(T) * size_;;
        memory_returned += sizeof(T) * size_;;
      }

      is_allocated_  = false;
      size_ = 0;
      data_ = nullptr;
    }

    void print(const std::string& label)
    {
      printf("%s (%d): \n", label.c_str(), size_);
      for (int i = 0; i < size_; i++) {
        printf("%s %d %g\n", label.c_str(), i+1, data_[i]);
      }
    }   

    /// @brief Resize the vector's memory.
    //
    void resize(const int size)
    { 
      if (size <= 0) {
        //throw std::runtime_error(+"Allocating a zero size Vector.");
        return;
      }
      if (data_ != nullptr) {
        if (reference_data_) { 
          throw std::runtime_error("[Vector] Can't resize a Vector with reference data.");
        }
        delete[] data_; 
        memory_in_use -= sizeof(T) * size_;;
        memory_returned += sizeof(T) * size_;;
        size_ = 0;
        data_ = nullptr;
      } 
      allocate(size);
    }
    
    /// @brief Grow the vector.
    //
    void grow(const int size, T value={})
    {
      if (size <= 0) {
        return;
      }

      memory_in_use += sizeof(T) * size;;
      int new_size = size_ + size;
      T* new_data = new T [new_size];
      for (int i = 0; i < size; i++) {
        new_data[i+size_] = value;
      }
      memcpy(new_data, data_, sizeof(T)*size_);
      if (reference_data_) { 
        throw std::runtime_error("[Vector] Can't grow a Vector with reference data.");
      }
      delete[] data_; 
      size_ = new_size;
      data_ = new_data;
    }

    void set_values(std::initializer_list<T> values)
    {
      if (values.size() == 0) {
        return;
      }
      check_type();
      allocate(values.size());
      std::copy(values.begin(), values.end(), data_);
    }

    void set_values(std::vector<T> values)
    {
      if (values.size() == 0) {
        return;
      }
      check_type();
      allocate(values.size());
      std::copy(values.begin(), values.end(), data_);
    }

    void read(const std::string& file_name) 
    { 
      auto fp = fopen(file_name.c_str(), "rb");
      int size;
      fread(&size, sizeof(int), 1, fp);
      fread(data_, sizeof(T), size_, fp);
      fclose(fp);
    }

    void write(const std::string& label, const T offset={}) const
    {
      if (!write_enabled) {
        return;
      }

      auto file_prefix = build_file_prefix(label);
      auto file_name = file_prefix + "_cm.bin";

      // Write binary file.
      //
      auto fp = fopen(file_name.c_str(), "wb");
      fwrite(&size_, sizeof(int), 1, fp);
      fwrite(data_, sizeof(T), size_, fp);
      fclose(fp);
    }

    /////////////////////////
    //  O p e r a t o r s  //
    /////////////////////////

    /// @brief Vector assigment.
    ///
    /// Note: There are two ways to do this:
    /// 
    ///   1) Swap pointers to data_
    ///
    ///   2) Copy data_
    ///
    /// Fortran using 2) I think.
    //
    Vector& operator=(const Vector& rhs)
    {
      if (rhs.size_ <= 0) {
        return *this;
      }

      if (this == &rhs) {
        return *this;
      }

      if (size_ != rhs.size_) {
        clear();
        allocate(rhs.size_);
      }

      memcpy(data_, rhs.data_, sizeof(T) * size_);

      return *this;
    }

    Vector& operator=(const double value)
    {
      for (int i = 0; i < size_; i++) {
        data_[i] = value;
      }
      return *this;
    }

    int size() const {
      return size_;
    }

    int msize() const
    {
      return size_ * sizeof(T);
    }

    // Index operators.
    //
    const T& operator()(const int i) const
    {
      #ifdef Vector_check_enabled
      check_index(i);
      #endif
      return data_[i];
    }

    T& operator()(const int i)
    {
      #ifdef Vector_check_enabled
      check_index(i);
      #endif
      return data_[i];
    }

    const T& operator[](const int i) const
    {
      #ifdef Vector_check_enabled
      check_index(i);
      #endif
      return data_[i];
    }

    T& operator[](const int i)
    {
      #ifdef Vector_check_enabled
      check_index(i);
      #endif
      return data_[i];
    }

    friend std::ostream& operator << (std::ostream& out, const Vector<T>& lhs)
    {
      for (int i = 0; i < lhs.size(); i++) {
        out << lhs[i];
        if (i != lhs.size()-1) {
          out << ", ";
         }
       }
       return out;
    }

    ///////////////////////////////////////////////////
    //  M a t h e m a t i c a l   o p e r a t o r s  //
    ///////////////////////////////////////////////////

    /// @brief Add and subtract vectors.
    //
    Vector<T> operator+(const Vector<T>& vec) const
    {
      if (size_ != vec.size()) {
        throw std::runtime_error("[Vector dot product] Vectors have diffrenct sizes: " + 
            std::to_string(size_) +  " != " + std::to_string(vec.size()) + ".");
      }
      Vector<T> result(size_);
      for (int i = 0; i < size_; i++) {
        result(i) = data_[i] + vec[i];
      }
      return result;
    }

    Vector<T> operator-(const Vector<T>& x) const
    {
      Vector<T> result(size_);
      for (int i = 0; i < size_; i++) {
        result(i) = data_[i] - x(i);
      }
      return result;
    }

    /// @brief Add and subtract a scalar from a vector.
    //
    Vector<T> operator+(const T value) const
    {
      Vector<T> result(size_);
      for (int i = 0; i < size_; i++) {
        result(i) = data_[i] + value;
      }
      return result;
    }

    friend const Vector<T> operator+(const T value, const Vector& rhs) 
    {
      Vector<T> result(rhs.size_);
      for (int i = 0; i < rhs.size_; i++) {
        result(i) = rhs.data_[i] + value;
      }
      return result;
    }

    Vector<T> operator-(const T value) const
    {
      Vector<T> result(size_);
      for (int i = 0; i < size_; i++) {
        result(i) = data_[i] - value;
      }
      return result;
    }

    friend Vector<T> operator-(T value, const Vector& rhs)
    {
      Vector<T> result(rhs.size_);
      for (int i = 0; i < rhs.size_; i++) {
        result(i) = value - rhs.data_[i];
      }
      return result;
    }

    /// @brief Divide by a scalar.
    //
    Vector<T> operator/(const T value) const
    {
      Vector<T> result(size_);
      for (int i = 0; i < size_; i++) {
        result(i) = data_[i] / value;
      }
      return result;
    }

    friend const Vector<T> operator/(const T value, const Vector& rhs)
    {
      Vector<T> result(rhs.size_);
      for (int i = 0; i < rhs.size_; i++) {
        result(i) = rhs.data_[i] / value;
      }
      return result;
    }

    /// @brief Multiply by a scalar.
    //
    Vector<T> operator*(const T value) const
    {
      Vector<T> result(size_);
      for (int i = 0; i < size_; i++) {
        result(i) = value * data_[i];
      }
      return result;
    }

    friend const Vector<T> operator*(const T value, const Vector& rhs)
    {
      Vector<T> result(rhs.size_);
      for (int i = 0; i < rhs.size_; i++) {
        result(i) = value * rhs.data_[i];
      }
      return result;
    }


    /// @brief Negate.
    Vector<T> operator-() const
    {
      Vector<T> result(size_);
      for (int i = 0; i < size_; i++) {
        result(i) = -data_[i];
      }
      return result;
    }

    /// @brief Absolute value.
    Vector<T> abs() const
    { 
      Vector<T> result(size_);
      for (int i = 0; i < size_; i++) {
        result(i) = std::abs(data_[i]);
      }
      return result;
    }

    /// @brief Cross product
    Vector<T> cross(const Vector<T>& v2)
    {
      Vector<T> result(size_);

      result(0) = (*this)(1)*v2(2) - (*this)(2)*v2(1); 
      result(1) = (*this)(2)*v2(0) - (*this)(0)*v2(2); 
      result(2) = (*this)(0)*v2(1) - (*this)(1)*v2(0); 

      return result;
    }

    /// @brief Dot product.
    T dot(const Vector<T>& v2)
    {
      if (size_ != v2.size()) {
        throw std::runtime_error("[Vector dot product] Vectors have diffrenct sizes: " + 
            std::to_string(size_) +  " != " + std::to_string(v2.size()) + ".");
      }
      T sum = 0.0;
      for (int i = 0; i < size_; i++) {
        sum += (*this)(i) * v2(i); 
      }
      return sum;
    }

    friend T operator*(const Vector& v1, const Vector& v2) 
    {
      if (v1.size() != v2.size()) {
        throw std::runtime_error("[Vector dot product] Vectors have diffrenct sizes: " + 
            std::to_string(v1.size()) +  " != " + std::to_string(v2.size()) + ".");
      }
      T sum = 0.0;
      for (int i = 0; i < v1.size(); i++) {
        sum += v1[i] * v2[i]; 
      }
      return sum;
    }

    T min() const
    {
      /*
      if (size_ == 0) {
        return std::numeric_limits<T>::max();
      }
      */
      return *std::min_element((*this).begin(), (*this).end());
    }

    T max() const
    {
      /*
      if (size_ == 0) {
        return -std::numeric_limits<T>::max();
      }
      */
      return *std::max_element((*this).begin(), (*this).end());
    }

    T sum() const
    {
      T sum = {};
      for (int i = 0; i < size_; i++) {
        sum += data_[i];
      }
      return sum;
    }

    /////////////////////////////////////
    //  i t e r a t o r   c l a s s s  //
    /////////////////////////////////////
   
    /// @brief This class provides an interface to access Vector like STL containers.
    //
    class Iterator
    {
      public:
        typedef T value_type;
        typedef T& reference;
        typedef T* pointer;
        typedef int difference_type;
        typedef std::forward_iterator_tag iterator_category;

        Iterator(T* ptr) : ptr_{ptr} {}

        Iterator& operator++() { this->ptr_ ++; return *this; }
        Iterator& operator--() { this->ptr_ --; return *this; }
        Iterator& operator++(int) { this->ptr_ ++; return *this; }
        Iterator& operator--(int) { this->ptr_ --; return *this; }
        T& operator*() { return *this->ptr_; };
        bool operator==(const Iterator& iter) { return this->ptr_ == iter.ptr_; }
        bool operator!=(const Iterator& iter) { return this->ptr_ != iter.ptr_; }
      private:
        T* ptr_;
    };

    Iterator begin() const
    {  
      return Iterator(data_);
    }

    Iterator end() const
    {  
      return Iterator(data_+size_);
    }

    T* data() const
    { 
      return data_;
    }

    void allocate(const int size)
    {
      if (size <= 0) {
        //throw std::runtime_error(+"Allocating a zero size Vector.");
        return;
      }

      size_ = size;
      data_ = new T [size_];
      memset(data_, 0, sizeof(T)*(size_));
      memory_in_use += sizeof(T)*size_;
    }

    void check_index(const int i) const
    {
      if (data_ == nullptr) {
        std::cout << "[Vector] WARNING: Accessing null data in Vector at " << i << std::endl;
        return;
        //throw std::runtime_error(+"Accessing null data in Vector.");
      }

      if ((i < 0) || (i >= size_)) {
        auto index_str = std::to_string(i);
        auto dims = std::to_string(size_);
        throw std::runtime_error( + "Index i=" + index_str + " is out of bounds for " + dims + " vector.");
      }
    }

    /// @brief Check that the Vector template type is int or double.
    //
    void check_type() const
    {
      if (!std::is_same<T, double>::value && !std::is_same<T, int>::value &&
          !std::is_same<T, Vector<double>>::value && !std::is_same<T, float>::value) {
        std::string msg = std::string("Cannot use Vector class template for type '") + typeid(T).name() + "'.";
        throw std::runtime_error(msg);
      }
    }

  private:
    bool is_allocated_ = false;
    int size_ = 0;
    bool reference_data_ = false;
    T *data_ = nullptr;
};

#endif

