"""
Script to fill holes in a surface mesh, and save the cap meshes as separate
vtp files. In general, there can be multiple holes to cap. The files are called
cap_#.vtp. If there are two caps, then this script will produce two files
- cap_0.vtp
- cap_1.vtp

You should rename these files as appropriate. 

Also, produces the capped surface mesh as endo_capped.vtp. This is useful to
check the geometry and normals of the cap.

USAGE:
    - modify set the full path to the surface mesh that you want to cap below,
    and run the script.
OR 
    - Run the script at the command line with the path to surface mesh that
    you want to cap
        python3 generate_cap.py path/to/mesh_to_cap.vtp
"""


import os
import pyvista as pv  # for VTK mesh manipulations
import vtk
from vtk import VTK_TRIANGLE, VTK_QUADRATIC_TRIANGLE, VTK_QUAD, VTK_QUADRATIC_QUAD, VTK_POLYGON
import numpy as np
import sys

## ------------------------- PARAMETERS TO CHANGE --------------------------- ## 
# Full path to surface mesh that you want to cap. This will be overwritten if 
# the script is run with a file path as a command-line argument.
mesh_to_cap = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./mesh/mesh-surfaces/endo.vtp")
## ------------------------- DO NOT MODIFY BELOW ---------------------------- ##

def save_data(file_name, data):

    """ 
    Write the given VTK object to a file. This function writes the mesh
    in a format that can be read by svFSI.                                                                                                                                                            
                                                                                                                                                                                                             
    Args:                                                                                                                                                                                                    
        file_name (str): The name of the file to write.                                                                                                                                                      
        data (mesh object): Mesh data to write. 
    Returns:
        None, saves data mesh as a file called file_name.                                                                                                                                                                                                                                                                                           
    """

    # Check filename format.
    file_ext = file_name.split(".")[-1]

    data = pv.wrap(data)
    data.points_to_double()

    if file_ext == "":
        raise RuntimeError("The file does not have an extension")

    # Get writer.
    if file_ext == "vtp":
        writer = vtk.vtkXMLPolyDataWriter()
    elif file_ext == "vtu":
        writer = vtk.vtkXMLUnstructuredGridWriter()
    else:
        raise RuntimeError("Unknown file type %s" % file_ext)

    # Set file name and the data to write.
    writer.SetFileName(file_name)
    writer.SetInputData(data)
    # compressor = vtk.vtkZLibDataCompressor()
    # compressor.SetCompressionLevel(9)
    # writer.SetCompressor(compressor)
    # writer.SetDataModeToAscii()
    # writer.EncodeAppendedDataOff()
    writer.Update()

    # Write the data.
    writer.Write()




# Optionally read mesh_to_cap as command line argument
if len(sys.argv) == 2:
    print(sys.argv)
    mesh_to_cap = os.path.abspath(sys.argv[1])
    print("Mesh to cap: ", mesh_to_cap)

# Set folder to output cap meshes
output_folder = os.path.dirname(mesh_to_cap)
print("\nOutput folder: ")
print(output_folder)

# Get base name of mesh to cap and remove file extension
mesh_to_cap_base = os.path.basename(mesh_to_cap).split(".")[0]
print(f"\nMesh to cap base name: {mesh_to_cap_base}")

print("\n## Capping holes ##")

# Load surface mesh that we want to cap
endo = pv.read(f"{mesh_to_cap}")

# Convert to unstructured grid. This is necessary if dealing with higher order
# cells (e.g. Quadratic triangles). We need to convert to polydata with simple
# triangles (linear cells), and to do that we need to first convert to 
# Unstructured Grid
endo = pv.UnstructuredGrid(endo)

# Determine the cell type (e.g. VTK_QUADRATIC_TRIANGLE, VTK_HEXAHEDRON) by
# checking the number of points per cell. If we don't properly set the celltypes
# array, extract_surface() will not produce the simple triangle polydata surface
# we need
#points_per_cell = endo.cell_n_points(0) # Check num points in first cell
points_per_cell = endo.get_cell(0).n_points
if points_per_cell == 3:
    endo.celltypes[:] = VTK_TRIANGLE
elif points_per_cell == 6:
    endo.celltypes[:] = VTK_QUADRATIC_TRIANGLE
elif points_per_cell == 4:
    endo.celltypes[:] = VTK_QUAD
elif points_per_cell == 8:
    endo.celltypes[:] = VTK_QUADRATIC_QUAD
else:
    endo.celltypes[:] = VTK_POLYGON
    print("Could not determine cell type")

# Extract the surface from the UnstructuredGrid. This produces a polydata mesh
# of simple triangles. For a mesh that originally had quadratic triangles, 
# each quadratic triangle is divided into 3 linear triangles.
endo = endo.extract_surface()

# Now we can fill the holes
endo_capped = endo.fill_holes(100)  # 100 is the largest size of hole to fill

# Recompute normals, incase the normals of the cap are opposite
endo_capped.compute_normals(inplace=True)

# Save capped endo polydata surface (to check geometry and normals)
# (Hopefully the normals on the filled cap will be consistent with the normals
# on the rest of the surface, but you should check to make sure.)
#endo_capped.save(f'{output_folder}/endo_capped.vtp')
save_data(os.path.join(output_folder, f"{mesh_to_cap_base}_capped.vtp"), endo_capped)

# Extract just the cap polydata surface. The cap cells were added to the end of
# endo when we capped then, so we can easily extract them by index
print(range(endo.number_of_cells, endo_capped.number_of_cells))
cap = endo_capped.extract_cells(
    range(endo.number_of_cells, endo_capped.number_of_cells)
)

# Extract surface again to get a polydata surface
cap = cap.extract_surface()

# Save cap surface
save_data(os.path.join(output_folder, f"{mesh_to_cap_base}_cap.vtp"), cap)

# Identify connected components of the surface. This will identify individual 
# cap surfaces, if there are more than one, and store them in the list caps.
#caps = cap.split_bodies()



# Save each cap individually as cap_0.vtp, cap_1.vtp, ...
# for k, cap in enumerate(caps):
#     # cap.save(f'{output_folder}/cap.vtp', binary = False)
#     save_data(os.path.join(output_folder,f"cap_{k}.vtp"), cap.extract_surface())




