// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/// @brief Define a map type used to set the bounds of element shape functions.
///
/// Sets   
/// \code {.cpp}
///   mesh.xib(2, com_mod.nsd);   // Bounds on Gauss integration points in parametric space
///   mesh.Nb2, mesh.eNoN);       // Bounds on shape functions
/// \endcode

using SetElementShapeBoundsMapType = std::map<ElementType, std::function<void(int, mshType&)>>;

SetElementShapeBoundsMapType set_element_shape_bounds_data = {

  {ElementType::HEX8, [](int nsd, mshType& mesh) -> void { 
      for (int i = 0; i < nsd; i++) {
        mesh.xib(0,i) = -1.0;
        mesh.xib(1,i) =  1.0;
      }
    }
  },

  {ElementType::LIN1, [](int nsd, mshType& mesh) -> void { 
    std::cout << "[set_element_shape_bounds_data] **************************" << std::endl;
    std::cout << "[set_element_shape_bounds_data] ERROR: LIN1 not supported." << std::endl;
    std::cout << "[set_element_shape_bounds_data] **************************" << std::endl;
    }
  },

  {ElementType::TET4, [](int nsd, mshType& mesh) -> void { 
      for (int i = 0; i < nsd; i++) {
        mesh.xib(0,i) = 0.0;
      }
    }
  },

  {ElementType::TRI3, [](int nsd, mshType& mesh) -> void { 
      for (int i = 0; i < nsd; i++) {
        mesh.xib(0,i) = 0.0;
      }
    }
  },

  {ElementType::WDG, [](int nsd, mshType& mesh) -> void { 
      mesh.xib(0,0) = 0.0;
      mesh.xib(0,1) = 0.0;
    }
  },

};

