// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIBERS_H
#define FIBERS_H

#include "Array.h"
#include "Array3.h"

#include "eigen3/Eigen/Core"

class ComMod;
class CmMod;
class MeshParameters;
class cmType;
class mshType;

/**
 * @brief Fiber directions associated with one mesh.
 */
class Fibers {
public:
  /// Location of the fiber directions in the input.
  enum class DataLocation {
    none,
    /// The fiber directions are associated with the mesh nodes.
    point,
    /// The fiber directions are associated with the mesh cells.
    cell,
    /// The fiber directions are constant vectors defined in the input file.
    constant
  };

  /// Fiber representation requested by a physics or output consumer.
  enum class Representation {
    raw,
    orthonormal_basis
  };

  /**
   * @brief Clears the fiber directions.
   */
  void clear();

  /**
   * @brief Reads fiber directions from file paths or constant vectors.
   */
  void read(const MeshParameters &mesh_parameters, mshType &mesh, int nsd);

  /**
   * @brief Parallel distribution of fiber directions.
   */
  void distribute(ComMod &com_mod, const CmMod &cm_mod, const cmType &cm,
                  mshType &mesh);

  /**
   * @brief Builds a completed orthonormal basis into @ref quadrature_fiber_basis.
   *
   * Returns immediately when a basis has already been built. Raw
   * @ref quadrature_fibers are left unchanged. The input must contain
   * nsd - 1 or nsd mutually orthogonal unit directions. For point data, the
   * independently interpolated directions are Gram--Schmidt orthonormalized at
   * every quadrature point.
   */
  void initialize_fiber_bases(int nsd);

  /**
   * @brief Returns the location of the fiber directions in the input.
   */
  DataLocation get_data_location() const { return data_location; }

  /**
   * @brief Returns the number of directions supplied by the user.
   */
  int get_num_input_directions() const { return num_input_fibers; }

  /**
   * @brief Returns the number of directions in the completed basis.
   *
   * Returns zero when no completed basis is available.
   */
  int get_num_basis_directions() const { return num_basis_fibers; }

  /**
   * @brief Returns the number of directions in a requested representation.
   */
  int get_num_directions(Representation representation) const;

  /**
   * @brief Returns the completed basis when available, otherwise raw fibers.
   */
  Representation get_default_representation() const
  {
    return has_orthonormal_basis()
      ? Representation::orthonormal_basis
      : Representation::raw;
  }

  /**
   * @brief Returns whether a completed orthonormal basis is available.
   */
  bool has_orthonormal_basis() const { return num_basis_fibers != 0; }

  /**
   * @brief Returns a requested fiber representation at a quadrature point.
   */
  Array<double> get_quadrature_fibers(int nsd, int g, int e,
                                      Representation representation) const;

  /**
   * @brief Evaluates a requested representation with supplied basis values.
   */
  Array<double> get_quadrature_fibers(
      int nsd, const Vector<double> &basis_functions, int e,
      Representation representation) const;

  /**
   * @brief Gathers the representation used for timestep-zero output.
   *
   * Raw directions are returned when no completed basis exists. Otherwise a
   * completed basis is returned at the original input association.
   */
  Array<double> gather_fibers(const ComMod &com_mod, const CmMod &cm_mod,
                              const mshType &mesh, int nsd) const;

  /**
   * @brief Returns a read-only Eigen view of a requested representation.
   *
   * The map aliases @ref quadrature_fibers or @ref quadrature_fiber_basis
   * and remains valid until those arrays are resized or cleared.
   */
  Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
  get_eigen_fibers(int nsd, int g, int e,
                   Representation representation) const;

  /**
   * @brief Returns true if there are fiber directions.
   */
  bool has_fibers() const { return num_input_fibers != 0; }

  /**
   * @brief Returns an nsd-by-1 zero fiber matrix for meshes with no directions.
   */
  static Array<double> zero_fibers(int nsd);

  /**
   * @brief Returns a read-only Eigen view of a zero nsd-by-1 fiber matrix.
   *
   * Backed by thread-local storage; valid until the next call to
   * zero_eigen_fibers on the same thread.
   */
  static Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
  zero_eigen_fibers(int nsd);

private:
  /**
   * @brief Fiber directions read from the input.
   *
   * Each column stores all supplied directions (nsd*num_input_fibers
   * components). Column count depends on the data location: one for constant
   * directions, one per cell for cell data, and one per point for point data.
   * Not distributed: after @ref distribute only the master keeps it, and only
   * for point data, where it is needed for timestep-zero output.
   */
  Array<double> input_fibers;

  /**
   * @brief Raw fiber directions at the quadrature points.
   */
  Array3<double> quadrature_fibers;

  /**
   * @brief Completed orthonormal fiber basis at the quadrature points.
   *
   * Empty until @ref initialize_fiber_bases. Sized nsd*nsd components per
   * location.
   */
  Array3<double> quadrature_fiber_basis;

  /**
   * @brief Point-associated input fibers gathered at each element's nodes.
   *
   * Keeping this element-local representation allows interpolation with
   * alternative basis functions after the mesh has been partitioned.
   */
  Array3<double> element_nodal_fibers;

  /// Location of the fiber directions in the input.
  DataLocation data_location = DataLocation::none;

  /// Number of fiber directions supplied as input.
  int num_input_fibers = 0;

  /// Number of directions in the completed basis, or zero when absent.
  int num_basis_fibers = 0;

  /**
   * @brief Reads fiber directions from VTK files.
   *
   * expected_location must be point or cell and selects which FIB_DIR
   * association is required.
   */
  void read_from_file(const MeshParameters &mesh_parameters, mshType &mesh,
                      int nsd, DataLocation expected_location);

  /**
   * @brief Reads constant fiber directions as specified in the input file.
   */
  void read_directions(const MeshParameters &mesh_parameters, int nsd);

  /**
   * @brief Checks that every input fiber direction has unit length.
   */
  void check_input_fibers_unit(int nsd) const;

  /**
   * @brief Checks mutual orthogonality of the input directions.
   */
  void check_input_fibers_orthogonal(int nsd) const;

  /**
   * @brief Populates the raw quadrature representation.
   */
  void populate_quadrature_fibers(const mshType &mesh, int nsd);

  /**
   * @brief Completes already-orthonormal input fibers to an nsd-dimensional
   * basis via cross product (3D) or 90-degree rotation (2D).
   */
  void complete_fiber_basis(Array<double> &fibers, int nsd,
                            int num_directions) const;

  /**
   * @brief Interpolates the original point directions for one element.
   *
   * Uses the supplied basis-function values with the cached
   * @ref element_nodal_fibers.
   */
  Array<double>
  interpolate_point_fibers_element(int nsd,
                                   const Vector<double> &basis_functions,
                                   int e) const;

  /**
   * @brief Interpolates nodal fiber directions onto the mesh quadrature
   * points.
   *
   * Builds @ref element_nodal_fibers and fills @ref quadrature_fibers using
   * the mesh shape functions.
   */
  void interpolate_point_fibers_mesh(const mshType &mesh, int nsd);
};

#endif
