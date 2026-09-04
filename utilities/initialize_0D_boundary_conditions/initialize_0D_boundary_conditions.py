#!/usr/bin/env python3
"""Update svZeroD coupled initial conditions from an svMultiPhysics VTU restart.

The script reads an svMultiPhysics XML file, finds:

  * Add_mesh/Initial_pressures_file_path
  * Add_mesh/Initial_velocities_file_path
  * svZeroDSolver_interface/Configuration_file
  * coupled Add_BC face names and svZeroDSolver_block names

It then integrates the VTU pressure/velocity on each coupled face and builds
the matching svZeroDSolver `initial_condition` keys.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def text_of(parent: ET.Element, tag: str) -> str | None:
    element = parent.find(tag)
    if element is None or element.text is None:
        return None
    return element.text.strip()


def resolve_path(base_dir: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return base_dir / path


def validate_vtk_xml_file(path: Path) -> None:
    try:
        with path.open("rb") as file:
            header = file.read(128)
    except OSError as error:
        raise RuntimeError(f"Could not open VTK XML file {path}: {error}") from error

    if header.startswith(b"version https://git-lfs.github.com/spec/v1"):
        raise RuntimeError(
            f"{path} is a Git LFS pointer file, not downloaded VTK data. "
            "Fetch the LFS contents before running this utility."
        )


def read_unstructured_grid(path: Path) -> vtk.vtkUnstructuredGrid:
    validate_vtk_xml_file(path)
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    if grid is None or grid.GetNumberOfPoints() == 0:
        raise RuntimeError(f"Could not read VTU grid from {path}")
    return grid


def read_polydata(path: Path) -> vtk.vtkPolyData:
    validate_vtk_xml_file(path)
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    polydata = reader.GetOutput()
    if polydata is None or polydata.GetNumberOfPoints() == 0:
        raise RuntimeError(f"Could not read VTP surface from {path}")
    return polydata


def make_locator(grid: vtk.vtkUnstructuredGrid) -> vtk.vtkStaticPointLocator:
    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(grid)
    locator.BuildLocator()
    return locator


def get_point_array(grid: vtk.vtkUnstructuredGrid, name: str) -> np.ndarray:
    array = grid.GetPointData().GetArray(name)
    if array is None:
        available = [
            grid.GetPointData().GetArrayName(i)
            for i in range(grid.GetPointData().GetNumberOfArrays())
        ]
        raise RuntimeError(f"Point data array {name!r} not found. Available: {available}")
    return vtk_to_numpy(array)


def face_average_pressure_and_flux(
    surface_path: Path,
    grid: vtk.vtkUnstructuredGrid,
    locator: vtk.vtkStaticPointLocator,
    pressure: np.ndarray,
    velocity: np.ndarray,
) -> tuple[float, float, float]:
    surface = read_polydata(surface_path)
    volume_point_ids = [
        locator.FindClosestPoint(surface.GetPoint(i))
        for i in range(surface.GetNumberOfPoints())
    ]

    point_area = np.zeros(len(volume_point_ids))
    oriented_flux = 0.0
    area = 0.0

    for cell_id in range(surface.GetNumberOfCells()):
        cell = surface.GetCell(cell_id)
        points = cell.GetPoints()
        num_points = points.GetNumberOfPoints()
        if num_points < 3:
            continue

        p0 = np.array(points.GetPoint(0))
        local_ids = [cell.GetPointId(i) for i in range(num_points)]

        for i in range(1, num_points - 1):
            p1 = np.array(points.GetPoint(i))
            p2 = np.array(points.GetPoint(i + 1))
            area_vector = 0.5 * np.cross(p1 - p0, p2 - p0)
            tri_area = float(np.linalg.norm(area_vector))
            if tri_area == 0.0:
                continue

            tri_local_ids = [local_ids[0], local_ids[i], local_ids[i + 1]]
            tri_volume_ids = [volume_point_ids[j] for j in tri_local_ids]
            tri_velocity = velocity[tri_volume_ids].mean(axis=0)

            oriented_flux += float(np.dot(tri_velocity, area_vector))
            area += tri_area
            for local_id in tri_local_ids:
                point_area[local_id] += tri_area / 3.0

    if area == 0.0 or point_area.sum() == 0.0:
        raise RuntimeError(f"Surface {surface_path} has no measurable area")

    face_pressures = pressure[volume_point_ids]
    average_pressure = float(np.dot(face_pressures, point_area) / point_area.sum())
    return average_pressure, oriented_flux, area


def parse_xml_config(xml_path: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    tree = ET.parse(xml_path)
    root = tree.getroot()
    base_dir = xml_path.parent

    mesh = root.find("Add_mesh")
    if mesh is None:
        raise RuntimeError("Could not find Add_mesh in XML")

    pressure_vtu = text_of(mesh, "Initial_pressures_file_path")
    velocity_vtu = text_of(mesh, "Initial_velocities_file_path")
    if not pressure_vtu or not velocity_vtu:
        raise RuntimeError("XML must define Initial_pressures_file_path and Initial_velocities_file_path")
    if pressure_vtu != velocity_vtu:
        raise RuntimeError(
            "This script expects pressure and velocity initial conditions in the same VTU. "
            f"Got {pressure_vtu!r} and {velocity_vtu!r}."
        )

    face_paths = {}
    for add_face in mesh.findall("Add_face"):
        name = add_face.attrib.get("name", "").strip()
        face_path = text_of(add_face, "Face_file_path")
        if name and face_path:
            face_paths[name] = resolve_path(base_dir, face_path)

    equation = root.find("Add_equation")
    if equation is None:
        raise RuntimeError("Could not find Add_equation in XML")

    interface = equation.find("svZeroDSolver_interface")
    if interface is None:
        raise RuntimeError("Could not find svZeroDSolver_interface in XML")

    json_file = text_of(interface, "Configuration_file")
    if not json_file:
        raise RuntimeError("Could not find svZeroDSolver_interface/Configuration_file in XML")

    coupled_faces = []
    for bc in equation.findall("Add_BC"):
        if text_of(bc, "Time_dependence") != "Coupled":
            continue
        face_name = bc.attrib.get("name", "").strip()
        coupling = bc.find("Coupling_interface")
        block_name = text_of(coupling, "svZeroDSolver_block") if coupling is not None else None
        if not face_name or not block_name:
            raise RuntimeError(f"Coupled BC is missing face name or svZeroDSolver_block: {ET.tostring(bc)}")
        if face_name not in face_paths:
            raise RuntimeError(f"Coupled BC {face_name!r} does not have a matching Add_face")
        coupled_faces.append(
            {
                "face_name": face_name,
                "surface_path": face_paths[face_name],
                "coupling_block": block_name,
            }
        )

    if not coupled_faces:
        raise RuntimeError("No coupled Add_BC entries found in XML")

    return resolve_path(base_dir, pressure_vtu), resolve_path(base_dir, json_file), coupled_faces


def build_initial_condition(
    json_config: dict,
    coupled_faces: list[dict[str, str]],
    measurements: dict[str, tuple[float, float, float]],
    flow_mode: str,
) -> dict[str, float]:
    coupling_blocks = {
        block["name"]: block["connected_block"]
        for block in json_config.get("external_solver_coupling_blocks", [])
    }
    rp_by_rcr = {
        bc["bc_name"]: float(bc["bc_values"]["Rp"])
        for bc in json_config.get("boundary_conditions", [])
        if bc.get("bc_type") == "RCR"
    }

    initial_condition = {}
    for face in coupled_faces:
        face_name = face["face_name"]
        coupling_block = face["coupling_block"]
        if coupling_block not in coupling_blocks:
            raise RuntimeError(f"Coupling block {coupling_block!r} from XML is not in the JSON")

        rcr_block = coupling_blocks[coupling_block]
        if rcr_block not in rp_by_rcr:
            raise RuntimeError(f"Connected RCR block {rcr_block!r} has no Rp in the JSON")

        pressure_value, oriented_flux, _area = measurements[face_name]
        flow_value = abs(oriented_flux) if flow_mode == "abs" else oriented_flux
        pressure_c = pressure_value - rp_by_rcr[rcr_block] * flow_value

        initial_condition[f"pressure_c:{rcr_block}"] = pressure_c
        initial_condition[f"flow:{coupling_block}:{rcr_block}"] = flow_value
        initial_condition[f"pressure:{coupling_block}:{rcr_block}"] = pressure_value

    return initial_condition


def format_initial_condition(initial_condition: dict[str, float]) -> str:
    return json.dumps({"initial_condition": initial_condition}, indent=2)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate svZeroD coupled initial_condition values from an svMultiPhysics VTU restart."
    )
    parser.add_argument(
        "xml",
        nargs="?",
        default="solver_vtu_and_bcs.xml",
        help="svMultiPhysics XML file. Default: solver_vtu_and_bcs.xml",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Update the JSON Configuration_file in place. Default: print only.",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not create a .bak copy when using --write.",
    )
    parser.add_argument(
        "--flow-mode",
        choices=["abs", "oriented"],
        default="abs",
        help="Use absolute or oriented surface flux for JSON flow values. Default: abs.",
    )
    parser.add_argument(
        "--pressure-array",
        default="Pressure",
        help="VTU point-data pressure array name. Default: Pressure.",
    )
    parser.add_argument(
        "--velocity-array",
        default="Velocity",
        help="VTU point-data velocity array name. Default: Velocity.",
    )
    args = parser.parse_args()

    xml_path = Path(args.xml).resolve()
    vtu_path, json_path, coupled_faces = parse_xml_config(xml_path)

    grid = read_unstructured_grid(vtu_path)
    locator = make_locator(grid)
    pressure = get_point_array(grid, args.pressure_array)
    velocity = get_point_array(grid, args.velocity_array)

    measurements = {}
    for face in coupled_faces:
        measurements[face["face_name"]] = face_average_pressure_and_flux(
            Path(face["surface_path"]),
            grid,
            locator,
            pressure,
            velocity,
        )

    with json_path.open() as json_file:
        json_config = json.load(json_file)

    initial_condition = build_initial_condition(
        json_config,
        coupled_faces,
        measurements,
        args.flow_mode,
    )

    print(f"VTU: {vtu_path}")
    print(f"JSON: {json_path}")
    print()
    print("Computed coupled face values:")
    for face in coupled_faces:
        pressure_value, oriented_flux, area = measurements[face["face_name"]]
        flow_value = abs(oriented_flux) if args.flow_mode == "abs" else oriented_flux
        print(
            f"  {face['face_name']}: pressure={pressure_value:.10e}, "
            f"flow={flow_value:.10e}, area={area:.10e}"
        )
    print()
    print(format_initial_condition(initial_condition))

    if args.write:
        if not args.no_backup:
            backup_path = json_path.with_suffix(json_path.suffix + ".bak")
            shutil.copy2(json_path, backup_path)
            print(f"\nBackup written: {backup_path}")
        json_config["initial_condition"] = initial_condition
        with json_path.open("w") as json_file:
            json.dump(json_config, json_file, indent=4)
            json_file.write("\n")
        print(f"Updated JSON: {json_path}")
    else:
        print("\nRun again with --write to update the JSON file.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
