#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
# -*-coding:utf-8 -*-
"""Surface name enumerations for biventricular heart models.

This module defines enumerations for surface names used in fiber generation
and Laplace solver configuration.
"""

from enum import Enum


class SurfaceName(Enum):
    """Enumeration of surface names for biventricular heart models.
    
    Attributes:
        EPICARDIUM: Epicardial surface
        BASE: Base surface (all valves together)
        EPICARDIUM_APEX: Epicardial apex surface
        ENDOCARDIUM_LV: Left ventricle endocardial surface
        ENDOCARDIUM_RV: Right ventricle endocardial surface
        MITRAL_VALVE: Mitral valve surface
        AORTIC_VALVE: Aortic valve surface
        TRICUSPID_VALVE: Tricuspid valve surface
        PULMONARY_VALVE: Pulmonary valve surface
    """
    EPICARDIUM = "epicardium"
    BASE = "base"
    EPICARDIUM_APEX = "epi_apex"
    ENDOCARDIUM_LV = "lv_endocardium"
    ENDOCARDIUM_RV = "rv_endocardium"
    MITRAL_VALVE = "mitral_valve"
    AORTIC_VALVE = "aortic_valve"
    TRICUSPID_VALVE = "tricuspid_valve"
    PULMONARY_VALVE = "pulmonary_valve"
    
    @classmethod
    def from_xml_face_name(cls, xml_name):
        """Convert XML face name to SurfaceName enum.
        
        Args:
            xml_name: XML face name (e.g., 'epicardium').
            
        Returns:
            SurfaceName: Corresponding enum value.
        """
        # Map XML face names to enum values
        xml_to_enum = {
            'epicardium': cls.EPICARDIUM,
            'base': cls.BASE,
            'epi_apex': cls.EPICARDIUM_APEX,
            'lv_endocardium': cls.ENDOCARDIUM_LV,
            'rv_endocardium': cls.ENDOCARDIUM_RV,
            'mitral_valve': cls.MITRAL_VALVE,
            'aortic_valve': cls.AORTIC_VALVE,
            'tricuspid_valve': cls.TRICUSPID_VALVE,
            'pulmonary_valve': cls.PULMONARY_VALVE,
        }
        return xml_to_enum.get(xml_name, None)
    
    @classmethod
    def get_required_for_method(cls, method):
        """Get required surface names for a given method.
        
        Args:
            method: Either "bayer" or "doste".
            
        Returns:
            set: Set of required SurfaceName enum values.
        """
        if method == "bayer":
            return {cls.EPICARDIUM, cls.ENDOCARDIUM_LV, cls.ENDOCARDIUM_RV, cls.BASE, cls.EPICARDIUM_APEX}
        elif method == "doste":
            return {cls.EPICARDIUM, cls.ENDOCARDIUM_LV, cls.ENDOCARDIUM_RV, cls.EPICARDIUM_APEX, 
                   cls.MITRAL_VALVE, cls.AORTIC_VALVE, cls.TRICUSPID_VALVE, cls.PULMONARY_VALVE}
        else:
            raise ValueError(f"Unknown method: {method}. Use 'bayer' or 'doste'.")
