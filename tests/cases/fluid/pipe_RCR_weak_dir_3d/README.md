
# **Problem Description**

Simulate unsteady fluid flow in a pipe. Uses weak Dirichlet boundary condition to
impose no-slip at the pipe wall.

# Inlet flow boundary condition

Time-dependent volumetric flow data is read in from the **lumen_inlet.flow** text file for the **lumen_inlet** face. The flow is converted into a **Parabolic** spatial profile.

```
<Add_BC name="lumen_inlet" > 
  <Type> Dir </Type> 
  <Time_dependence> Unsteady </Time_dependence> 
  <Temporal_values_file_path> lumen_inlet.flow</Temporal_values_file_path> 
  <Profile> Parabolic </Profile> 
  <Impose_flux> true </Impose_flux> 
</Add_BC>
```

# Outlet RCR boundary condition

An RCR boundary condition is defined for the **lumen_outlet** outlet face.

```
<Add_BC name="lumen_outlet" > 
  <Type> Neu </Type> 
  <Time_dependence> RCR </Time_dependence> 
  <RCR_values> 
    <Capacitance> 1.5e-5 </Capacitance> 
    <Distal_resistance> 1212 </Distal_resistance> 
    <Proximal_resistance> 121 </Proximal_resistance> 
     <Distal_pressure> 0 </Distal_pressure> 
    <Initial_pressure> 0 </Initial_pressure> 
  </RCR_values> 
</Add_BC>
```
