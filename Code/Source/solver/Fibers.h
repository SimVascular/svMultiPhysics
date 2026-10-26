// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FIBERS_H
#define FIBERS_H

#include "Array.h"
#include "Array3.h"

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
   * @brief Returns the location of the fiber directions in the input.
   */
  DataLocation get_data_location() const { return data_location; }

  /// Returns the number of raw directions supplied by the user.
  int get_num_directions() const { return num_input_fibers; }

  /**
   * @brief Returns fibers at a quadrature point.
   *
   * Raw directions are returned by default. When @p orthonormalize is true,
   * an owned, completed orthonormal basis is returned without modifying or
   * caching the raw fibers.
   */
  Array<double> get_quadrature_fibers(int nsd, int g, int e,
                                      bool orthonormalize = false) const;

  /**
   * @brief Evaluates fibers with supplied basis-function values.
   *
   * For point data, raw directions are re-interpolated with the supplied
   * values before optional orthonormalization.
   */
  Array<double> get_quadrature_fibers(
      int nsd, const Vector<double> &basis_functions, int e,
      bool orthonormalize = false) const;

  /**
   * @brief Gathers raw fibers for timestep-zero output.
   */
  Array<double> gather_fibers(const ComMod &com_mod, const CmMod &cm_mod,
                              const mshType &mesh, int nsd) const;

  /**
   * @brief Returns true if there are fiber directions.
   */
  bool has_fibers() const { return num_input_fibers != 0; }

  /**
   * @brief Returns an nsd-by-1 zero fiber matrix for meshes with no directions.
   */
  static Array<double> zero_fibers(int nsd);

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

  /// Whether every raw input location contains mutually orthogonal directions.
  bool basis_input_is_orthogonal = true;

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
   * @brief Returns whether input directions are mutually orthogonal.
   */
  bool input_fibers_are_orthogonal(int nsd) const;

  /**
   * @brief Populates the raw quadrature representation.
   */
  void populate_quadrature_fibers(const mshType &mesh, int nsd);

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
