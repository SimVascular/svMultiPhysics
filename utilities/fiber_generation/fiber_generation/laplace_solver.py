#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
# -*-coding:utf-8 -*-
"""Laplace solver module for biventricular fiber generation.

This module provides the LaplaceSolver class that generates XML configuration
files and runs the svMultiPhysics Laplace-Dirichlet solver for computing
scalar fields needed for fiber direction generation.

Supports both:
    - Bayer et al. (2012): Truncated BiV geometry
    - Doste et al. (2019): BiV geometry with outflow tracts
"""

import os
from xml.etree import ElementTree as ET
from xml.dom import minidom
from .surface_names import SurfaceName


class LaplaceSolver:
    """Laplace solver for biventricular fiber generation.
    
    This class generates XML configuration files for the svMultiPhysics solver
    and executes the Laplace-Dirichlet problems needed for fiber generation.
    
    Attributes:
        mesh_path: Path to the volumetric mesh file (.vtu).
        surface_paths: Dictionary mapping SurfaceName enum values to full file paths.
        exec_svmultiphysics: Command to execute svMultiPhysics solver.
    
    Example:
        >>> from src.SurfaceNames import SurfaceName
        >>> solver = LaplaceSolver(
        ...     mesh_path="/path/to/mesh.vtu",
        ...     surface_paths={
        ...         SurfaceName.EPICARDIUM: '/path/to/epicardium.vtp',
        ...         SurfaceName.ENDOCARDIUM_LV: '/path/to/lv_endocardium.vtp',
        ...         SurfaceName.ENDOCARDIUM_RV: '/path/to/rv_endocardium.vtp',
        ...         SurfaceName.BASE: '/path/to/base.vtp',
        ...         SurfaceName.EPICARDIUM_APEX: '/path/to/apex.vtp',
        ...     },
        ...     exec_svmultiphysics="svmultiphysics "
        ... )
        >>> result_file = solver.run("bayer", "/path/to/output")
    """
    
    # Solver configuration constants
    DEFAULT_TOLERANCE = 1e-6
    DEFAULT_MAX_LS_ITERATIONS = 2000
    DEFAULT_MAX_NL_ITERATIONS = 5
    
    def __init__(self, mesh_path, surface_paths, exec_svmultiphysics = "svmultiphysics "):
        """Initialize the LaplaceSolver.
        
        Args:
            mesh_path: Path to the volumetric mesh file (.vtu).
            surface_paths: Dictionary mapping SurfaceName enum values to full file paths.
                For Bayer method, required: SurfaceName.EPICARDIUM, ENDOCARDIUM_LV, ENDOCARDIUM_RV, BASE, EPICARDIUM_APEX
                For Doste method, required: SurfaceName.EPICARDIUM, ENDOCARDIUM_LV, ENDOCARDIUM_RV, EPICARDIUM_APEX, MITRAL_VALVE, AORTIC_VALVE, TRICUSPID_VALVE, PULMONARY_VALVE
            exec_svmultiphysics: Command to execute svMultiPhysics (e.g., "svmultiphysics ").
            
        Raises:
            TypeError: If surface_paths contains non-enum keys.
        """
        self.mesh_path = mesh_path
        
        # Validate that all keys are SurfaceName enum values
        for key in surface_paths.keys():
            if not isinstance(key, SurfaceName):
                raise TypeError(
                    f"surface_paths must use SurfaceName enum keys, not {type(key).__name__}. "
                    f"Got key: {key}"
                )
        
        self.surface_paths = surface_paths
        self.exec_svmultiphysics = exec_svmultiphysics
    
    def _create_general_params(self, output_dir):
        """Create the GeneralSimulationParameters XML element.
        
        Args:
            output_dir: Directory to save results.
            
        Returns:
            ET.Element: XML element for general simulation parameters.
        """
        params = ET.Element("GeneralSimulationParameters")
        
        elements = [
            ("Continue_previous_simulation", "0"),
            ("Number_of_spatial_dimensions", "3"),
            ("Number_of_time_steps", "1"),
            ("Time_step_size", "1"),
            ("Spectral_radius_of_infinite_time_step", "0."),
            ("Searched_file_name_to_trigger_stop", "STOP_SIM"),
            ("Save_results_to_VTK_format", "1"),
            ("Name_prefix_of_saved_VTK_files", "result"),
            ("Increment_in_saving_VTK_files", "1"),
            ("Save_results_in_folder", output_dir),
            ("Start_saving_after_time_step", "1"),
            ("Increment_in_saving_restart_files", "1"),
            ("Convert_BIN_to_VTK_format", "0"),
            ("Verbose", "1"),
            ("Warning", "1"),
            ("Debug", "0"),
        ]
        
        for name, value in elements:
            elem = ET.SubElement(params, name)
            elem.text = f" {value} "
        
        return params
    
    def _create_mesh_element(self, face_names):
        """Create the Add_mesh XML element with faces.
        
        Args:
            face_names: List of face names to include.
            
        Returns:
            ET.Element: XML element for mesh definition.
        """
        mesh = ET.Element("Add_mesh", name="msh")
        
        mesh_path_elem = ET.SubElement(mesh, "Mesh_file_path")
        mesh_path_elem.text = f" {self.mesh_path} "
        
        for face_name in face_names:
            face = ET.SubElement(mesh, "Add_face", name=face_name)
            face_path = ET.SubElement(face, "Face_file_path")
            face_path.text = f" {self._get_surface_path(face_name)} "
        
        return mesh
    
    def _get_surface_path(self, face_name):
        """Get the surface file path for a given face name.
        
        Args:
            face_name: Internal face name used in XML.
            
        Returns:
            str: Path to the surface file.
            
        Raises:
            KeyError: If the surface is not found in surface_paths.
        """
        # Convert XML face name to SurfaceName enum
        surface_enum = SurfaceName.from_xml_face_name(face_name)
        
        if surface_enum is None:
            raise KeyError(f"Unknown XML face name: {face_name}")
        
        if surface_enum not in self.surface_paths:
            raise KeyError(f"Surface {surface_enum.value} not found in surface_paths")
        
        return self.surface_paths[surface_enum]
    
    def _create_equation(self, output_alias, boundary_conditions):
        """Create a heat equation XML element for Laplace problem.
        
        Args:
            output_alias: Name for the output temperature field.
            boundary_conditions: List of tuples (face_name, value) for Dirichlet BCs.
            
        Returns:
            ET.Element: XML element for the heat equation.
        """
        eq = ET.Element("Add_equation", type="heatS")
        
        # Equation parameters
        eq_params = [
            ("Coupled", "0"),
            ("Min_iterations", "1"),
            ("Max_iterations", str(self.DEFAULT_MAX_NL_ITERATIONS)),
            ("Tolerance", str(self.DEFAULT_TOLERANCE)),
            ("Conductivity", "1.0"),
            ("Source_term", "0.0"),
            ("Density", "0.0"),
        ]
        
        for name, value in eq_params:
            elem = ET.SubElement(eq, name)
            elem.text = f" {value} "
        
        # Output configuration
        output_spatial = ET.SubElement(eq, "Output", type="Spatial")
        temp_spatial = ET.SubElement(output_spatial, "Temperature")
        temp_spatial.text = " 1 "
        
        output_alias_elem = ET.SubElement(eq, "Output", type="Alias")
        temp_alias = ET.SubElement(output_alias_elem, "Temperature")
        temp_alias.text = f" {output_alias} "
        
        # Linear solver
        ls = ET.SubElement(eq, "LS", type="CG")
        la = ET.SubElement(ls, "Linear_algebra", type="fsils")
        precond = ET.SubElement(la, "Preconditioner")
        precond.text = " fsils "
        max_iter = ET.SubElement(ls, "Max_iterations")
        max_iter.text = f" {self.DEFAULT_MAX_LS_ITERATIONS} "
        tol = ET.SubElement(ls, "Tolerance")
        tol.text = f" {self.DEFAULT_TOLERANCE} "
        
        # Boundary conditions
        for face_name, value in boundary_conditions:
            bc = ET.SubElement(eq, "Add_BC", name=face_name)
            bc_type = ET.SubElement(bc, "Type")
            bc_type.text = " Dir "
            bc_value = ET.SubElement(bc, "Value")
            bc_value.text = f" {value} "
            zero_perim = ET.SubElement(bc, "Zero_out_perimeter")
            zero_perim.text = " 0 "
        
        return eq
    
    def _get_bayer_equations(self):
        """Get equation definitions for the Bayer method.
        
        Returns:
            list: List of (alias, boundary_conditions) tuples.
        """
        return [
            # Trans_EPI: Transmural field (epicardium=1, endo=0)
            ("Trans_EPI", [
                ("epicardium", 1.0),
                ("lv_endocardium", 0.0),
                ("rv_endocardium", 0.0),
            ]),
            # Trans_LV: LV field (lv_endo=1, others=0)
            ("Trans_LV", [
                ("lv_endocardium", 1.0),
                ("rv_endocardium", 0.0),
                ("epicardium", 0.0),
            ]),
            # Trans_RV: RV field (rv_endo=1, others=0)
            ("Trans_RV", [
                ("rv_endocardium", 1.0),
                ("lv_endocardium", 0.0),
                ("epicardium", 0.0),
            ]),
            # Long_AB: Apex-to-base field (base=1, apex=0)
            ("Long_AB", [
                ("base", 1.0),
                ("epi_apex", 0.0),
            ]),
        ]
    
    def _get_doste_equations(self):
        """Get equation definitions for the Doste method.
        
        Returns:
            list: List of (alias, boundary_conditions) tuples.
        """
        return [
            # Trans_BiV: Ventricular transmural (LV=-2, RV=1)
            ("Trans_BiV", [
                ("lv_endocardium", -2.0),
                ("rv_endocardium", 1.0),
            ]),
            # Long_AV: LV longitudinal from AV (apex=1, aortic_valve=0)
            ("Long_AV", [
                ("epi_apex", 1.0),
                ("aortic_valve", 0.0),
            ]),
            # Long_MV: LV longitudinal from MV (apex=1, mitral_valve=0)
            ("Long_MV", [
                ("epi_apex", 1.0),
                ("mitral_valve", 0.0),
            ]),
            # Long_PV: RV longitudinal from PV (apex=1, pulmonary_valve=0)
            ("Long_PV", [
                ("epi_apex", 1.0),
                ("pulmonary_valve", 0.0),
            ]),
            # Long_TV: RV longitudinal from TV (apex=1, tricuspid_valve=0)
            ("Long_TV", [
                ("epi_apex", 1.0),
                ("tricuspid_valve", 0.0),
            ]),
            # Weight_LV: LV valve weights (mitral_valve=1, aortic_valve=0)
            ("Weight_LV", [
                ("aortic_valve", 0.0),
                ("mitral_valve", 1.0),
            ]),
            # Weight_RV: RV valve weights (tricuspid_valve=1, pulmonary_valve=0)
            ("Weight_RV", [
                ("pulmonary_valve", 0.0),
                ("tricuspid_valve", 1.0),
            ]),
            # Trans_EPI: Epicardial transmural (epicardium=1, endo=0)
            ("Trans_EPI", [
                ("epicardium", 1.0),
                ("lv_endocardium", 0.0),
                ("rv_endocardium", 0.0),
            ]),
            # Trans_LV: LV transmural (lv_endo=1, others=0)
            ("Trans_LV", [
                ("epicardium", 0.0),
                ("rv_endocardium", 0.0),
                ("lv_endocardium", 1.0),
            ]),
            # Trans_RV: RV transmural (rv_endo=1, others=0)
            ("Trans_RV", [
                ("epicardium", 0.0),
                ("lv_endocardium", 0.0),
                ("rv_endocardium", 1.0),
            ]),
            # Trans: transmural (rv_endo=1, epi=0, lv_endo=-2)
            ("Trans", [
                ("epicardium", 0.0),
                ("lv_endocardium", -2.0),
                ("rv_endocardium", 1.0),
            ]),
        ]
    
    def _get_face_names(self, method):
        """Get the list of face names for a given method.
        
        Args:
            method: Either "bayer" or "doste".
            
        Returns:
            list: List of face names.
        """
        if method == "bayer":
            return ["epicardium", "base", "epi_apex", "lv_endocardium", "rv_endocardium"]
        elif method == "doste":
            return ["epicardium", "mitral_valve", "aortic_valve", "tricuspid_valve", "pulmonary_valve", "epi_apex", "lv_endocardium", "rv_endocardium"]
        else:
            raise ValueError(f"Unknown method: {method}. Use 'bayer' or 'doste'.")
    
    def _prettify_xml(self, elem):
        """Convert an XML element to a prettified string.
        
        Args:
            elem: ET.Element to prettify.
            
        Returns:
            str: Prettified XML string.
        """
        rough_string = ET.tostring(elem, encoding='unicode')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="\t")
    
    def generate_solver_xml(self, method, output_dir):
        """Generate the svMultiPhysics XML configuration file.
        
        Args:
            method: Either "bayer" or "doste".
            output_dir: Directory to save solver results.
            
        Returns:
            str: Path to the generated XML file.
        """
        # Create root element
        root = ET.Element("svMultiPhysicsFile", version="0.1")
        
        # Add general parameters
        root.append(self._create_general_params(output_dir))
        
        # Add mesh with faces
        face_names = self._get_face_names(method)
        root.append(self._create_mesh_element(face_names))
        
        # Add equations based on method
        if method == "bayer":
            equations = self._get_bayer_equations()
        elif method == "doste":
            equations = self._get_doste_equations()
        else:
            raise ValueError(f"Unknown method: {method}. Use 'bayer' or 'doste'.")
        
        for alias, bcs in equations:
            root.append(self._create_equation(alias, bcs))
        
        # Generate XML string
        xml_declaration = '<?xml version="1.0" encoding="UTF-8"?>\n'
        xml_content = self._prettify_xml(root)
        
        # Remove the default XML declaration from prettify and use our own
        if xml_content.startswith('<?xml'):
            xml_content = xml_content.split('?>', 1)[1].strip()
        xml_content = xml_declaration + xml_content
        
        # Determine output path (same directory as mesh file)
        xml_dir = os.path.dirname(self.mesh_path)
        xml_path = os.path.join(xml_dir, "svFibers_BiV.xml")
        
        with open(xml_path, 'w') as f:
            f.write(xml_content)
        
        return xml_path
    
    def run(self, method, output_dir, delete_xml = False):
        """Generate XML and run the svMultiPhysics Laplace solver.
        
        This is the main entry point that generates the XML configuration
        and executes the solver.
        
        Args:
            method: Either "bayer" or "doste".
            output_dir: Directory to save solver results.
            
        Returns:
            str: Path to the result file containing Laplace solutions.
            
        Raises:
            ValueError: If required surfaces are missing for the specified method.
        """
        # Validate surfaces before running
        self.validate_surfaces(method)
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate XML configuration
        xml_path = self.generate_solver_xml(method, output_dir)
        print(f"Generated solver XML at: {xml_path}")
        
        if delete_xml:
            os.remove(xml_path)
        
        # Run solver
        print("   Running svMultiPhysics solver")
        print(f"   {self.exec_svmultiphysics + xml_path}")
        os.system(self.exec_svmultiphysics + xml_path)
        
        return os.path.join(output_dir, 'result_001.vtu')
    
    def validate_surfaces(self, method):
        """Validate that all required surfaces are provided.
        
        Args:
            method: Either "bayer" or "doste".
            
        Raises:
            ValueError: If required surfaces are missing.
        """
        required = SurfaceName.get_required_for_method(method)
        available = set(self.surface_paths.keys())
        
        missing = required - available
        if missing:
            missing_names = [s.value for s in missing]
            raise ValueError(f"Missing required surfaces for {method} method: {missing_names}")
