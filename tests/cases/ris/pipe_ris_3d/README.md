# **Problem Description**
Simulation of a 3D pipe flow problem with an explicit RIS (resistive immersed surfaces) valve model.



The model parameters are specified in the `Add_RIS_projection` sub-section
```
<Add_RIS_projection name="left_ris" >
  <Project_from_face> right_ris </Project_from_face>
  <Resistance> 1.e6 </Resistance>
  <Projection_tolerance> 1.e-8 </Projection_tolerance>
</Add_RIS_projection>
```

## Reference
Astorino, M., Hamers, J., Shadden, S. C., & Gerbeau, J. F. (2012). A robust and efficient valve model based on resistive immersed surfaces. International Journal for Numerical Methods in Biomedical Engineering, 28(9), 937-959. https://doi.org/10.1002/cnm.2474.