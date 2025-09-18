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

#ifndef VTK_DATA_H
#define VTK_DATA_H

#include "Array.h"
#include "Vector.h"

#include <string>

class VtkData {
  public:
    VtkData();
    virtual ~VtkData();

    virtual Array<int> get_connectivity() const = 0;
    virtual Array<double> get_points() const = 0;
    virtual int num_elems() const = 0;
    virtual int elem_type() const = 0;
    virtual int np_elem() const = 0;
    virtual int num_points() const = 0;
    virtual void read_file(const std::string& file_name) = 0;

    virtual void set_element_data(const std::string& data_name, const Array<double>& data) = 0;
    virtual void set_element_data(const std::string& data_name, const Array<int>& data) = 0;

    virtual void set_point_data(const std::string& data_name, const Array<double>& data) = 0;
    virtual void set_point_data(const std::string& data_name, const Array<int>& data) = 0;
    virtual void set_point_data(const std::string& data_name, const Vector<int>& data) = 0;

    virtual void set_points(const Array<double>& points) = 0;
    virtual void set_connectivity(const int nsd, const Array<int>& conn, const int pid = 0) = 0;

    virtual bool has_point_data(const std::string& data_name) = 0;

    virtual void copy_points(Array<double>& points) = 0;
    virtual void copy_point_data(const std::string& data_name, Array<double>& mesh_data) = 0;
    virtual void copy_point_data(const std::string& data_name, Vector<double>& mesh_data) = 0;
    virtual void write() = 0;

    static VtkData* create_reader(const std::string& file_name);
    static VtkData* create_writer(const std::string& file_name);

    std::string file_name;
};

class VtkVtpData : public VtkData {
  public:
    VtkVtpData();
    VtkVtpData(const std::string& file_name, bool reader=true);
    ~VtkVtpData();

    virtual Array<int> get_connectivity() const override;
    virtual Array<double> get_points() const override;
    virtual int elem_type() const override;
    virtual int num_elems() const override;
    virtual int np_elem() const override;
    virtual int num_points() const override;
    virtual void read_file(const std::string& file_name) override;

    void copy_points(Array<double>& points) override;
    void copy_point_data(const std::string& data_name, Array<double>& mesh_data) override;
    void copy_point_data(const std::string& data_name, Vector<double>& mesh_data) override;
    void copy_point_data(const std::string& data_name, Vector<int>& mesh_data);
    Array<double> get_point_data(const std::string& data_name);
    std::vector<std::string> get_point_data_names();
    bool has_point_data(const std::string& data_name) override;
    virtual void set_connectivity(const int nsd, const Array<int>& conn, const int pid = 0) override;

    virtual void set_element_data(const std::string& data_name, const Array<double>& data) override;
    virtual void set_element_data(const std::string& data_name, const Array<int>& data) override;

    virtual void set_point_data(const std::string& data_name, const Array<double>& data) override;
    virtual void set_point_data(const std::string& data_name, const Array<int>& data) override;
    virtual void set_point_data(const std::string& data_name, const Vector<int>& data) override;

    virtual void set_points(const Array<double>& points) override;
    virtual void write() override;

  private:
    class VtkVtpDataImpl;
    VtkVtpDataImpl* impl;
};

class VtkVtuData : public VtkData {
  public:
    VtkVtuData();
    VtkVtuData(const std::string& file_name, bool reader=true);
    ~VtkVtuData();

    virtual Array<int> get_connectivity() const override;
    virtual int elem_type() const override;
    virtual int num_elems() const override;
    virtual int np_elem() const override;
    virtual int num_points() const override;
    virtual void read_file(const std::string& file_name) override;

    void copy_points(Array<double>& points) override;
    void copy_point_data(const std::string& data_name, Array<double>& mesh_data) override;
    void copy_point_data(const std::string& data_name, Vector<double>& mesh_data) override;
    void copy_point_data(const std::string& data_name, Vector<int>& mesh_data);

    Array<double> get_point_data(const std::string& data_name);
    std::vector<std::string> get_point_data_names();
    virtual Array<double> get_points() const override;
    bool has_point_data(const std::string& data_name) override;
    virtual void set_connectivity(const int nsd, const Array<int>& conn, const int pid = 0) override;

    virtual void set_element_data(const std::string& data_name, const Array<double>& data) override;
    virtual void set_element_data(const std::string& data_name, const Array<int>& data) override;

    virtual void set_point_data(const std::string& data_name, const Array<double>& data) override;
    virtual void set_point_data(const std::string& data_name, const Array<int>& data) override;
    virtual void set_point_data(const std::string& data_name, const Vector<int>& data) override;

    virtual void set_points(const Array<double>& points) override;
    virtual void write() override;

  private:
    class VtkVtuDataImpl;
    VtkVtuDataImpl* impl;
};

#endif
