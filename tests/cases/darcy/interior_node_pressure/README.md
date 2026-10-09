# **Problem Description**

Prescribe Darcy pressure at interior nodes of a 2D unit square. The pressure is 10 at an injection node at (0.25, 0.5) and 0 along three producer nodes at x = 0.75. No faces are defined; the outer boundary is therefore impermeable (zero normal flux), and fluid flows from the injector to the producer.

A prescribed pressure fixes the solution value at the selected nodes. The corresponding injection or extraction rate is part of the solution; it is not a source term.

The compressibility is zero and `Spectral_radius_of_infinite_time_step` is 0, so the saved pressure is the steady solution.

Run the example from this directory with `svmultiphysics solver.xml`, or through the test suite with `pytest tests/test_darcy.py -k interior_node_pressure`.

## Node sets

A node set names existing points of a volume mesh. It is defined inside `Add_mesh` and needs no face file, surface elements, areas or normals, so it can select isolated or disconnected interior nodes.

```
<Add_mesh name="msh" >
  <Mesh_file_path> mesh/mesh-complete.mesh.vtu </Mesh_file_path>

  <Add_node_set name="injector">
    <Node_IDs> 216 </Node_IDs>
  </Add_node_set>

  <Add_node_set name="producer">
    <Node_IDs_file_path> producer_nodes.dat </Node_IDs_file_path>
  </Add_node_set>
</Add_mesh>
```

* Node IDs are one-based point indices in the file given by `Mesh_file_path`: ID 1 is the first point in that file. They are not values of a `GlobalNodeID` array, and the same IDs select the same points for any number of MPI processes.
* Give exactly one of `Node_IDs` (whitespace-separated IDs) or `Node_IDs_file_path` (a text file of whitespace-separated IDs).
* A set must contain at least one ID, and IDs may not repeat. Set names must be unique within a mesh.

## Dirichlet conditions on node sets

A Dirichlet condition selects a node set with `Mesh_name` and `Node_set`, which must be given together. The `Add_BC` name is then only a label and is not looked up as a face.

```
<Add_BC name="injection_pressure" >
  <Mesh_name> msh </Mesh_name>
  <Node_set> injector </Node_set>
  <Type> Dirichlet </Type>
  <Time_dependence> Steady </Time_dependence>
  <Value> 10.0 </Value>
</Add_BC>
```

* The equation containing the condition determines the prescribed variable: the scalar unknown of a scalar equation, such as Darcy pressure or temperature, or the velocity (or, with `Impose_on_state_variable_integral`, displacement) components of a vector equation.
* `Effective_direction` selects Cartesian components, for example `1 0 0`. Without it, every component of the equation's Dirichlet variable receives the value. Values are never multiplied by a normal.
* `Time_dependence` may be `Steady` (`Value`), `Unsteady` (`Temporal_values_file_path` or `Fourier_coefficients_file_path`), or `General` (`Temporal_and_spatial_values_file_path`). `Steady` and `Unsteady` apply one value to every selected node; `General` gives each node its own values.

### Per-node values

A `General` values file lists the number of components, the number of time points and the number of nodes, followed by the time points and one record per node of the set:

```
1 2 4
0.0 1.0
205  0.0  0.0
226  0.0  0.0
247  0.0  0.0
216 10.0 10.0
```

This file assigns constant values to the nodes of a set that contains nodes 205, 226, 247 and 216. The first time is zero and times increase strictly; at least two are required. Each record starts with a node ID and then lists, for each time point in turn, one value per component. Every node of the set appears exactly once, in any order. The number of components must equal the number of prescribed components. Values are interpolated linearly in time and repeat with a period equal to the last time, as for face conditions; repeat the same value at two times to hold it constant.

## Restrictions

* Only strongly imposed Dirichlet values are supported. `Weakly_applied`, `Apply_along_normal_direction`, `Impose_flux`, `Zero_out_perimeter`, profiles other than `Flat`, face data files (`Spatial_profile_file_path`, `Spatial_values_file_path`, `Bct_file_path`, `Traction_values_file_path`), `CST_shell_bc_type`, `Coupling_interface`, `Undeforming_neu_face`, `Follower_pressure_load`, and time dependences other than those above are rejected.
* Node sets are not supported in `ustruct` equations, including FSI with a `ustruct` solid, or in simulations that require remeshing.
* A component prescribed by a node set may not be prescribed by any other Dirichlet condition of the same equation, face or node set, even with the same value; the solver stops with an error instead of choosing one by input order. Different conditions may prescribe different components of the same node.
