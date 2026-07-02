# **Problem Description**
Simulation of a 3D pipe flow problem with an unfitted RIS (resistive immersed surfaces) valve model.



The model parameters are specified in the `Add_URIS_mesh` sub-section
```
<Add_URIS_mesh name="MV" > 
  <Add_URIS_face name="LCC" > 
    <Face_file_path> meshes/LCC_mesh.vtu </Face_file_path>
    <Open_motion_file_path> meshes/LCC_motion_open.dat </Open_motion_file_path>
    <Close_motion_file_path> meshes/LCC_motion_close.dat </Close_motion_file_path>
  </Add_URIS_face>
  <Mesh_scale_factor> 1.0 </Mesh_scale_factor>
  <Thickness> 0.2 </Thickness>
  <Closed_thickness> 0.2 </Closed_thickness>
  <Resistance> 1.0e5 </Resistance>
  <Invert_normal> false </Invert_normal>
  <Include_URIS_velocity> true </Include_URIS_velocity>
  <Positive_flow_normal_file_path> meshes/normal.dat </Positive_flow_normal_file_path>
  <Scaffold_file_path> meshes/mv_scaffold.vtu </Scaffold_file_path>
</Add_URIS_mesh>
```

In this test, the valve velocity $u_{\Gamma_i}$ associated with the prescribed valve motion is included in the URIS forcing term by setting `Include_URIS_velocity` as `true`,

$$
f_{\mathrm{URIS}} = \sum_{i=1}^{n} \left( u - u_{\Gamma_i} \right) \text{ .}
$$

The valve velocity $u_{\Gamma_i}$ is computed using finite differences based on the input valve motion data. If the input motion is sampled with relatively large time steps, the resulting velocity may be noisy or inaccurate. To improve the temporal resolution of the prescribed motion, the valve trajectory can be interpolated using the utility provided in

`svMultiPhysics/utilities/interpolate_uris_valve_velocity/`

This script generates a refined valve motion dataset with additional intermediate frames, leading to a smoother and more accurate estimate of $u_{\Gamma_i}$.