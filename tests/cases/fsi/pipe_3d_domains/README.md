# Test Case: Tissue engineered vascular graft with two domains
 
This test case demonstrates the capability to read domains from a `.vtu` file. 

## Mesh

The mesh is a cylinder consisting of a lumen mesh and a wall mesh. The wall mesh is divided into two domains (defined in `mesh/domainIDs.vtu`) as shown below (blue=1, red=2),

![Domains](mesh/domains.png)

For comparison, a `domainIDs.dat` file is also included. Running `solver_dat.xml` will use this file instead. Both approaches are equivalent.
  
## Results

The following animations show the velocity and displacement for 300 timesteps.

![Results](vid.gif)