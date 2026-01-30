#!/usr/bin/env python
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
    EPICARDIUM = "epi"
    BASE = "base"
    EPICARDIUM_APEX = "epi_apex"
    ENDOCARDIUM_LV = "endo_lv"
    ENDOCARDIUM_RV = "endo_rv"
    MITRAL_VALVE = "mv"
    AORTIC_VALVE = "av"
    TRICUSPID_VALVE = "tv"
    PULMONARY_VALVE = "pv"
    
    @classmethod
    def from_xml_face_name(cls, xml_name):
        """Convert XML face name to SurfaceName enum.
        
        Args:
            xml_name: XML face name (e.g., 'epi_top').
            
        Returns:
            SurfaceName: Corresponding enum value.
        """
        # Map XML face names to enum values
        xml_to_enum = {
            'epi': cls.EPICARDIUM,
            'epi_top': cls.BASE,
            'epi_apex': cls.EPICARDIUM_APEX,
            'endo_lv': cls.ENDOCARDIUM_LV,
            'endo_rv': cls.ENDOCARDIUM_RV,
            'mv': cls.MITRAL_VALVE,
            'av': cls.AORTIC_VALVE,
            'tv': cls.TRICUSPID_VALVE,
            'pv': cls.PULMONARY_VALVE,
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
