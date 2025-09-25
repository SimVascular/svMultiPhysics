This test case simulates a spatially variable Robin boundary condition on a slab of material described by the Guccione material model. This case is identical to `struct/spatially_variable_robin`, except it uses `ustruct` physics.

- Primary fibers run along the length of the slab (z-direction) and secondary fibers run across the width of the slab
(x-direction).

- The slab is loaded on the +Y surface with a uniform pressure load. The load profile is a ramp to 10 dynes/cm^2 over 0.5 seconds, then held there until
2 seconds. The load is defined in `load.dat`, which can be generated with
`utilities/generate_boundary_condition_data/generate_load.py`. The load tends to push the slab downward. 

![Load Profile](load.png)

- This is resisted by a spatially varying Robin boundary condition on the -Y surface. The stiffness is 0 at z = 0, and 50 at the far end. This is provided in `Y0_spatially_varying_robin.vtp`, which can be generated with `utilities/generate_boundary_condition_data/generate_spatially_variable_robin.py`.

![Spatially varying Robin BC](Y0_spatially_varying_robin.png)


- The slab is also constrained by Dirichlet boundary conditions on the +-X and +-Z
surfaces, applied in the normal direction. These prevent the slab from moving
in the x and z directions.

- The resulting deformation is shown in the video below:
![Deformation](animation.gif)
The black outline shows the initial configuration. As you can see, the displacement of the slab is greatest at z = 0, where the Robin BC stiffness is zero. At the far end of the slab, the displacement is very little, where the stiffness is greatest. The oscillations are due to the absence of any damping in the Robin BC.

- Note, the deformation is fairly different from `struct/spatially_variable_robin`, like due to the different underlying physics formulation and the coarse mesh. Therefore, the `result_002.vtu` files are not identical.

