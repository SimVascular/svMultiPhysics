// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef VTK_DATA_H
#define VTK_DATA_H

#include "Array.h"
#include "Vector.h"

#include <string>
#include <utility>
#include <vector>

#include <vtkIdList.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

/**
 * @brief A mesh stored in one of the VTK XML file formats.
 *
 * The mesh is held as a vtkPointSet, and consists of the point coordinates,
 * the element connectivity, and any number of named data arrays associated
 * with the points, with the elements, or with the mesh as a whole (field
 * data). The class provides the operations needed to read those from a file,
 * to build them up before writing a file, and to copy them to and from the
 * Array and Vector types used by the solver.
 *
 * The file format determines how the mesh is represented, which is what the
 * derived classes provide: VtkVtpData holds a surface mesh as a vtkPolyData,
 * VtkVtuData a volume mesh as a vtkUnstructuredGrid. Use create_reader() and
 * create_writer() to obtain an object of the type matching a given file name.
 */
class VtkData {
  public:
    /**
     * @brief Default constructor.
     */
    VtkData() = default;

    /**
     * @brief Virtual destructor.
     */
    virtual ~VtkData() = default;

    /**
     * @brief Read the mesh data from a VTK file.
     *
     * @throws svmp::FileFormatException if the file cannot be read, or if it
     *   contains no points or no elements
     */
    virtual void read_file(const std::string &file_name);

    /**
     * @brief Create an empty grid.
     *
     * Derived classes need to override this function by implementing the
     * appropriate logic to create an empty grid. This will include initializing
     * the correct vtkPointSet object.
     */
    virtual void create_grid() = 0;

    /**
     * @brief Write the mesh data to a VTK file.
     *
     * @throws svmp::CoreException if the file cannot be written
     */
    virtual void write() const = 0;

    /**
     * @brief Get the connectivity of the mesh elements.
     *
     * @return An array of size (num_points_per_elem, num_elems) containing the
     *   connectivity of the mesh elements. Each column corresponds to an
     *   element, and each row corresponds to a point index in that element.
     */
    Array<int> get_connectivity() const;

    /**
     * @brief Get the points of the mesh.
     *
     * @return An array of size (3, num_points) containing the coordinates of
     * the mesh points. Each column corresponds to a point, and each row
     * corresponds to a coordinate (x, y, z).
     */
    Array<double> get_points() const;

    /**
     * @brief Get the number of elements in the mesh.
     */
    int num_elems() const;

    /**
     * @brief Get the element type.
     */
    int elem_type() const;

    /**
     * @brief Get the number of points per element.
     */
    int num_points_per_elem() const;

    /**
     * @brief Get the number of points in the mesh.
     */
    int num_points() const;

    /**
     * @brief Set a double-valued element data array.
     *
     * @param[in] data_name The name of the data array to set.
     * @param[in] data The data array to set, holding one value per element.
     *
     * @throws svmp::FE::InvalidArgumentException if the number of values
     *   differs from the number of elements of the mesh
     */
    void set_element_data(const std::string &data_name,
                          const Array<double> &data);

    /**
     * @brief Set an int-valued element data array.
     *
     * @param[in] data_name The name of the data array to set.
     * @param[in] data The data array to set, holding one value per element.
     *
     * @throws svmp::FE::InvalidArgumentException if the number of values
     *   differs from the number of elements of the mesh
     */
    void set_element_data(const std::string &data_name, const Array<int> &data);

    /**
     * @brief Set an int-valued element data vector.
     *
     * @param[in] data_name The name of the data array to set.
     * @param[in] data The data vector to set, holding one value per element.
     *
     * @throws svmp::FE::InvalidArgumentException if the number of values
     *   differs from the number of elements of the mesh
     */
    void set_element_data(const std::string &data_name,
                          const Vector<int> &data);

    /**
     * @brief Set a double-valued point data array.
     *
     * @param[in] data_name The name of the data array to set.
     * @param[in] data The data array to set, holding one value per point.
     *
     * @throws svmp::FE::InvalidArgumentException if the number of values
     *   differs from the number of points of the mesh
     */
    void set_point_data(const std::string &data_name,
                        const Array<double> &data);

    /**
     * @brief Set an int-valued point data array.
     *
     * @param[in] data_name The name of the data array to set.
     * @param[in] data The data array to set, holding one value per point.
     *
     * @throws svmp::FE::InvalidArgumentException if the number of values
     *   differs from the number of points of the mesh
     */
    void set_point_data(const std::string &data_name, const Array<int> &data);

    /**
     * @brief Set an int-valued point data array.
     *
     * @param[in] data_name The name of the data array to set.
     * @param[in] data The data array to set, holding one value per point.
     *
     * @throws svmp::FE::InvalidArgumentException if the number of values
     *   differs from the number of points of the mesh
     */
    void set_point_data(const std::string &data_name, const Vector<int> &data);

    /**
     * @brief Set the point coordinates of the mesh.
     *
     * @param[in] points The coordinates, of size (3, num_points).
     *
     * @throws svmp::FE::InvalidArgumentException if there are no points or
     *   fewer than three coordinates are given for each point
     */
    void set_points(const Array<double> &points);

    /**
     * @brief Set the mesh connectivity to define the elements.
     *
     * The elements are appended to those already defined, so a mesh made of
     * several parts is built up by calling this once per part.
     *
     * @param[in] nsd The number of spatial dimensions, which together with the
     *   number of points per element determines the element type.
     * @param[in] conn The connectivity, of size (num_points_per_elem,
     *   num_elems). Each column holds the point indices of one element.
     *
     * @throws svmp::FE::InvalidArgumentException if a point index does not
     *   refer to one of the points of the mesh
     */
    void set_connectivity(const int nsd, const Array<int> &conn);

    /**
     * @brief Store a time value as field data.
     *
     * The value is written as a single-tuple Float64 field data array named
     * 'TimeValue'. VTK XML readers such as ParaView turn this array into the
     * pipeline time of the data object.
     *
     * @param[in] time The time value to associate with the data.
     */
    void set_time_value(const double time);

    /**
     * @brief Check if a given cell data array exists.
     */
    bool has_cell_data(const std::string &data_name) const;

    /**
     * @brief Check if a given point data array exists.
     */
    bool has_point_data(const std::string &data_name) const;

    /**
     * @brief Copy the mesh points to an Array.
     *
     * @param[out] points The array to copy the mesh points into. It must be
     *   of size (3, num_points).
     *
     * @throws svmp::FE::InvalidArgumentException if the points do not fit in
     *   the array
     */
    void copy_points(Array<double> &points) const;

    /**
     * @brief Copy an array of point data from the mesh into the given Array.
     *
     * @param[in] data_name The name of the point data array to copy.
     * @param[out] mesh_data The array to copy the point data into. It must be
     *   of size (num_components, num_points).
     *
     * @throws svmp::FE::InvalidArgumentException if the mesh has no
     *   double-valued point data array with the given name, or if its values do
     *   not fit in mesh_data
     */
    void copy_point_data(const std::string &data_name,
                         Array<double> &mesh_data) const;

    /**
     * @brief Copy an array of point data from the mesh into the given Vector.
     *
     * @param[in] data_name The name of the point data array to copy.
     * @param[out] mesh_data The vector to copy the point data into. It must be
     *   of size (num_points).
     *
     * @throws svmp::FE::InvalidArgumentException if the mesh has no
     *   double-valued point data array with the given name, or if its values do
     *   not fit in mesh_data
     */
    void copy_point_data(const std::string &data_name,
                         Vector<double> &mesh_data) const;

    /**
     * @brief Copy an array of int-valued point data from the mesh into the
     * given Vector.
     *
     * @param[in] data_name The name of the point data array to copy.
     * @param[out] mesh_data The vector to copy the point data into. It must be
     *   of size (num_points).
     *
     * @throws svmp::FE::InvalidArgumentException if the mesh has no int-valued
     *   point data array with the given name, or if its values do not fit in
     *   mesh_data
     */
    void copy_point_data(const std::string &data_name,
                         Vector<int> &mesh_data) const;

    /**
     * @brief Copy an array of cell data from the mesh into the given Array.
     *
     * @param[in] data_name The name of the cell data array to copy.
     * @param[out] mesh_data The array to copy the cell data into. It must be
     *   of size (num_components, num_cells).
     *
     * @throws svmp::FE::InvalidArgumentException if the mesh has no
     *   double-valued cell data array with the given name, or if its values do
     *   not fit in mesh_data
     */
    void copy_cell_data(const std::string &data_name,
                        Array<double> &mesh_data) const;

    /**
     * @brief Copy an array of cell data from the mesh into the given Vector.
     *
     * @param[in] data_name The name of the cell data array to copy.
     * @param[out] mesh_data The vector to copy the cell data into. It must be
     *   of size (num_points).
     *
     * @throws svmp::FE::InvalidArgumentException if the mesh has no
     *   double-valued cell data array with the given name, or if its values do
     *   not fit in mesh_data
     */
    void copy_cell_data(const std::string &data_name,
                        Vector<double> &mesh_data) const;

    /**
     * @brief Copy an array of int-valued cell data from the mesh into the given
     * Vector.
     *
     * @param[in] data_name The name of the cell data array to copy.
     * @param[out] mesh_data The vector to copy the cell data into. It must be
     *   of size (num_points).
     *
     * @throws svmp::FE::InvalidArgumentException if the mesh has no int-valued
     *   cell data array with the given name, or if its values do not fit in
     *   mesh_data
     */
    void copy_cell_data(const std::string &data_name,
                        Vector<int> &mesh_data) const;

    /**
     * @brief Get an array of point data from the mesh.
     *
     * @throws svmp::FE::InvalidArgumentException if the mesh has no
     *   double-valued point data array with the given name
     *
     * @todo[michelebucelli] This should fall back onto copy_point_data.
     */
    Array<double> get_point_data(const std::string &data_name) const;

    /**
     * @brief Get a list of point data names.
     *
     * @return A vector of strings containing the names of the point data
     *   arrays.
     */
    std::vector<std::string> get_point_data_names() const;

    /**
     * @brief Get the dimensions of a cell data array.
     *
     * @param[in] data_name The name of the cell data array to get the
     *   dimensions of.
     *
     * @return A pair of integers representing the number of components and the
     *   number of tuples in the array.
     */
    std::pair<int, int>
    get_cell_data_dimensions(const std::string &data_name) const;

    /**
     * @brief Create an object to read a mesh from a VTK file.
     *
     * The concrete type is selected from the file extension, and the mesh is
     * read as part of the construction. The file extension must be 'vtp' or
     * 'vtu'.
     *
     * @param[in] file_name The name of the VTK file to read.
     *
     * @return A pointer to a newly allocated object holding the mesh read from
     *   the file. The caller owns the object and must delete it.
     *
     * @throws svmp::FE::InvalidArgumentException if the file extension is not
     *   'vtp' or 'vtu'
     */
    static VtkData *create_reader(const std::string &file_name);

    /**
     * @brief Create an object to write a mesh to a VTK file.
     *
     * The concrete type is selected from the file extension, and the object is
     * created holding an empty mesh. The file extension must be 'vtp' or 'vtu'.
     * The file itself is written by write(), once the mesh has been defined.
     *
     * @param[in] file_name The name of the VTK file to write.
     *
     * @return A pointer to a newly allocated object holding an empty mesh. The
     *   caller owns the object and must delete it.
     *
     * @throws svmp::FE::InvalidArgumentException if the file extension is not
     *   'vtp' or 'vtu'
     */
    static VtkData *create_writer(const std::string &file_name);

  protected:
    /**
     * @brief Read the mesh data from a file.
     *
     * Derived classes need to override this function by implementing the
     * appropriate reading logic. This will include selecting the correct VTK
     * reader, and initializing vtk_data.
     */
    virtual void read_file_internal(const std::string &file_name) = 0;

    /**
     * @brief Get the VTK cell type for a given number of spatial dimensions and
     * number of points per element.
     *
     * Derived classes must override this to implement the appropriate mapping
     * from the number of spatial dimensions and number of points per element to
     * the VTK cell type.
     *
     * @throws svmp::FE::InvalidArgumentException if there is no cell type that
     *   the file format can hold for the given element
     */
    virtual int cell_type(int nsd, int np_elem) const = 0;

    /**
     * @brief Insert a new cell into the VTK data object.
     *
     * Derived classes must override this to implement the appropriate insertion
     * call.
     *
     * @param[in] vtk_cell_type The VTK cell type of the element.
     * @param[in] elem_nodes The list of point IDs that define the element.
     */
    virtual void insert_cell(int vtk_cell_type,
                             vtkSmartPointer<vtkIdList> elem_nodes) = 0;
    /**
     * Pointer to the underlying VTK data object. The concrete type will be
     * either vtkPolyData, for VTP files, or vtkUnstructuredGrid, for VTU files.
     */
    vtkSmartPointer<vtkPointSet> vtk_data;

    /// Filename that the mesh is read from, or written to.
    std::string file_name_;

    /// Type of elements in the mesh.
    int elem_type_ = -1;

    /// Number of elements.
    int num_elems_ = 0;

    /// Number of points per element.
    int num_points_per_elem_ = 0;

    /// Number of points.
    int num_points_ = 0;
};

/**
 * @brief A mesh stored in the VTK XML polygonal data format ('.vtp' files).
 *
 * The mesh is held as a vtkPolyData, and its elements are polygons: the
 * surface meshes that define the faces of a volume mesh, and the line meshes
 * used for one-dimensional domains.
 */
class VtkVtpData : public VtkData {
public:
  /**
   * @brief Default constructor.
   *
   * Creates an empty VtkVtpData object with an empty vtkPolyData grid.
   */
  VtkVtpData();

  /**
   * @brief Constructor.
   *
   * @param[in] file_name The name of the VTP file to read from or write to.
   * @param[in] reader If true, the constructor reads the mesh data from the
   *   given file. If false, it creates an empty grid.
   */
  VtkVtpData(const std::string &file_name, bool reader = true);

  /**
   * @brief Create an empty grid.
   */
  virtual void create_grid() override;

  /**
   * @brief Write the mesh data to a file.
   */
  virtual void write() const override;

protected:
  /**
   * @brief Read the mesh data from a file.
   */
  virtual void read_file_internal(const std::string &file_name) override;

  /**
   * @brief Get the VTK cell type of a surface element.
   */
  virtual int cell_type(int nsd, int np_elem) const override;

  /**
   * @brief Insert a new cell into the vtkPolyData object.
   */
  virtual void insert_cell(int vtk_cell_type,
                           vtkSmartPointer<vtkIdList> elem_nodes) override;

  /**
   * The mesh, represented as a vtkPolyData object.
   */
  vtkSmartPointer<vtkPolyData> vtk_polydata;
};

/**
 * @brief A mesh stored in the VTK XML unstructured grid format ('.vtu' files).
 *
 * The mesh is held as a vtkUnstructuredGrid, whose elements may be of any VTK
 * cell type: the volume meshes the equations are solved on, and the results
 * written at each saved time step.
 */
class VtkVtuData : public VtkData {
public:
  /**
   * @brief Default constructor.
   *
   * Creates an empty VtkVtuData object with an empty vtkUnstructuredGrid grid.
   */
  VtkVtuData();

  /**
   * @brief Constructor.
   *
   * @param[in] file_name The name of the VTP file to read from or write to.
   * @param[in] reader If true, the constructor reads the mesh data from the
   *   given file. If false, it creates an empty grid.
   */
  VtkVtuData(const std::string &file_name, bool reader = true);

  /**
   * @brief Create an empty grid.
   */
  virtual void create_grid() override;

  /**
   * @brief Write the mesh data to a file.
   */
  virtual void write() const override;

protected:
  /**
   * @brief Read the mesh data from a file.
   */
  virtual void read_file_internal(const std::string &file_name) override;

  /**
   * @brief Get the VTK cell type of a volume element.
   */
  virtual int cell_type(int nsd, int np_elem) const override;

  /**
   * @brief Insert a new cell into the vtkUnstructuredGrid object.
   */
  virtual void insert_cell(int vtk_cell_type,
                           vtkSmartPointer<vtkIdList> elem_nodes) override;

  /**
   * The mesh, represented as a vtkUnstructuredGrid object.
   */
  vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid;
};

#endif
