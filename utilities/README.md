# Utilities

The utilities directory contains pre- and post-processing scripts that can be useful for svMultiphysics. 

The following scripts are provided:
- `fourier_coefficients`: generation of fourier coefficients after providing temporal flow data.
- `fiber_generation`: generate myocardial fiber orientations for biventricular heart models using the Bayer et al. (2012) and Doste et al. (2018) rule-based methods.
- `generate_boundary_condition_data`: generate svMultiPhysics boundary condition inputs, including transient load profiles and spatially varying Robin boundary condition fields defined from user-specified analytical functions.
- `generate_cap_surface`: generate cap surfaces for open VTP surface meshes by filling holes, recomputing normals, and saving both the capped mesh and the extracted cap as separate VTP files.
- `interpolate_uris_valve_velocity`: generate interpolated prescribed valve motion for selected time intervals with a user-specified number of output frames to improve valve velocity estimation.