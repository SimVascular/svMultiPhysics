# Test Case: Tensile Adventitia with Holzapfel-Ogden Model and Directional Active Stress

This test case demonstrates prescribed active contraction of a slab of material
described by the Holzapfel-Ogden Modified Anisotropy material model. 

## Mesh

The mesh is a simple rectangular slab with dimensions: 4.95792 cm × 1.48517 cm × 0.254564 cm

## Fiber Configuration

Primary fibers run along the length of the slab (z-direction) and secondary fibers 
run across the width of the slab (x-direction), with the sheet-normal direction 
in the y-direction.

## Active Stress Distribution

This test case demonstrates the **directional active stress distribution feature**, which
allows the total active stress to be distributed among the three orthogonal fiber directions:

- **70%** in the fiber direction (z-direction)
- **20%** in the sheet direction (x-direction)
- **10%** in the sheet-normal direction (y-direction)

The fractions are specified in the solver.xml file as:
```xml
<Fiber_reinforcement_stress type="Unsteady" >
   <Temporal_values_file_path> stress.dat </Temporal_values_file_path> 
   <Directional_distribution>
      <Fiber_direction> 0.7 </Fiber_direction>
      <Sheet_direction> 0.2 </Sheet_direction>
      <Sheet_normal_direction> 0.1 </Sheet_normal_direction>
   </Directional_distribution>
</Fiber_reinforcement_stress>
```

These fractions must sum to 1.0.

## Stress Profile

The active stress profile is a linear ramp from 0 to 50 dynes/cm² over 1 second. 
The active stress data is defined in `stress.dat`, which can be regenerated using 
`generate_stress.py`.

![Stress Profile](media/stress.png)

## Material Model

The Holzapfel-Ogden Modified Anisotropy model is used with parameters from 
`LV_HolzapfelOgden_active`, but in CGS units (dynes/cm²).

## Boundary Conditions
- Z0 face is constrained in z-direction
- X0 face is constrained in x-direction
- Y0 face is constrained in y-direction
  
## Purpose

This test case serves to test and demonstrate the directional active stress distribution feature.

## Results

The following animations show the deformation under different active stress distributions.
Note that the displacements have been scaled 1000x in Paraview for easy visualization.

### Combined Stress Distribution (70% fiber, 20% sheet, 10% sheet-normal)
![Combined Animation](media/animation.gif)

### Pure Fiber Direction Stress (100% fiber)
![Fiber Animation](media/animation_fiber_1.gif)

### Pure Sheet Direction Stress (100% sheet)
![Sheet Animation](media/animation_sheet_1.gif)

### Pure Sheet-Normal Direction Stress (100% sheet-normal)
![Sheet-Normal Animation](media/animation_sheet_normal_1.gif)

