
# **Problem Description**

Simulate dye transport in a cylinder govered by a given velocity field. The dye is passively transported by the velocity field through advection and diffusion.

## Scalar transport equation
The advection-diffusion equation that governs dye transport uses a known velocity field read in from a VTK-format file named **precomputed_velocity.vtu** containing a *PointData DataArray* named **Velocity**.

```
<Precomputed_solution>
  <Use_precomputed_solution> true </Use_precomputed_solution>
  <File_path> precomputed_velocity.vtu </File_path>
  <Field_name> Velocity </Field_name>
</Precomputed_solution>
```

## Inflow boundary condition
The inflow boundary condition for a constant release of dye at the **lumen_inlet** inlet face.
```
<Add_BC name="lumen_inlet" >
  <Type> Dirichlet </Type>
  <Time_dependence> Steady </Time_dependence>
  <Value> 1.0 </Value>
</Add_BC>
```
