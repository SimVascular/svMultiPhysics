# **Problem Description**
Simulation of a 3D pipe FSI problem with an unfitted RIS (resistive immersed surfaces) valve model.



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
  <Resistance> 1.0e5 </Resistance>
  <Positive_flow_normal_file_path> meshes/normal.dat </Positive_flow_normal_file_path>
</Add_URIS_mesh>
```