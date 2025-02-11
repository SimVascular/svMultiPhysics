/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/// @brief Define a map type used to set element properties.
using SetElementPropsMapType = std::map<int, std::function<void(int, mshType&)>>;

/// @brief Map used to set 3D element properties.
///
/// This replicates the case statement in the Fortran 'SUBROUTINE SELECTELE(lM)' 
/// defined in NN.f. 
//
SetElementPropsMapType set_3d_element_props = {

  {4, [](int insd, mshType& mesh) -> void { 
    mesh.eType = ElementType::TET4; 
    mesh.nG = 4; 
    mesh.vtkType = 10; 
    mesh.nEf = 4; 
    mesh.lShpF = true; 
    }
  },

 {6, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::WDG;
    mesh.nG = 6;
    mesh.vtkType = 13;
    mesh.nEf = 3;
    mesh.lShpF = true;
    }
  },

 {8, [](int insd, mshType& mesh) -> void { 
    mesh.eType = ElementType::HEX8;
    mesh.nG = 8; 
    mesh.vtkType = 12;
    mesh.nEf = 6; 
    mesh.lShpF = false;
    }
  },

 {20, [](int insd, mshType& mesh) -> void { 
    mesh.eType = ElementType::HEX20;
    mesh.nG = 27; 
    mesh.vtkType = 25;
    mesh.nEf = 6; 
    mesh.lShpF = false;
    }
  },

 {27, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::HEX27;
    mesh.nG = 27;
    mesh.vtkType = 29;
    mesh.nEf = 6;
    mesh.lShpF = false;
    }
  },

 {10, [](int insd, mshType& mesh) -> void { 
    mesh.eType = ElementType::TET10;
    mesh.nG = 15; 
    mesh.vtkType = 24;
    mesh.nEf = 4; 
    mesh.lShpF = false;
    }
  },

 {20, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::HEX20;
    mesh.nG = 27;
    mesh.vtkType = 25;
    mesh.nEf = 6;
    mesh.lShpF = false;
    }
  },

 {27, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::HEX27;
    mesh.nG = 27;
    mesh.vtkType = 29;
    mesh.nEf = 6;
    mesh.lShpF = false;
    }
  },
};

//---------------------
// set_2d_element_props
//---------------------
// Map used to set 2D element properties.
//
SetElementPropsMapType set_2d_element_props = {

  {3, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::TRI3;
    mesh.nG = 3;
    mesh.vtkType = 5;
    mesh.nEf = 3;
    mesh.lShpF = true;
    }
  },

  {4, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::QUD4;
    mesh.nG = 4;
    mesh.vtkType = 9;
    mesh.nEf = 4;
    mesh.lShpF = false;
    }
  },

  {6, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::TRI6;
    mesh.nG = 7;
    mesh.vtkType = 22;
    mesh.nEf = 3;
    mesh.lShpF = false;
    }
  },

  {8, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::QUD8;
    mesh.nG = 9;
    mesh.vtkType = 23;
    mesh.nEf = 4;
    mesh.lShpF = false;
    }
  },

  {9, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::QUD9;
    mesh.nG = 9;
    mesh.vtkType = 28;
    mesh.nEf = 4;
    mesh.lShpF = false;
    }
  },
};

/// @brief Map used to set 1D element properties.
//
SetElementPropsMapType set_1d_element_props = {

  {2, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::LIN1;
    mesh.nG = 2;
    mesh.vtkType = 3;
    mesh.nEf = 2;
    mesh.lShpF = true;
    }
  },

  {3, [](int insd, mshType& mesh) -> void {
    mesh.eType = ElementType::LIN2;
    mesh.nG = 3;
    mesh.vtkType = 21;
    mesh.nEf = 2;
    mesh.lShpF = false;
    }
  },
};

/// @brief Define a map type used to set face properties.
using SetFacePropsMapType = std::map<int, std::function<void(int, faceType&)>>;

/// @brief Map used to set face properties.
///
/// set_face_element_props[eNoN].
///
/// Replaces the Fortan 'select' statemnt in 'SELECTELEB'.
//
SetFacePropsMapType set_face_element_props = {

  {1, [](int insd, faceType& face) -> void 
    {
      face.eType = ElementType::PNT;
      face.nG = 1;
    }
  },

  // Two nodes per element.
  //
  {2, [](int insd, faceType& face) -> void 
    {
      if (insd == 1) {
        face.eType = ElementType::LIN1;
        face.nG = 2;
      }
    }
  },

  {3, [](int insd, faceType& face) -> void 
    {
      if (insd == 2) {
         face.eType = ElementType::TRI3;
      } else {
         face.eType = ElementType::LIN2;
      }
      face.nG = 3;
    }
  },

  {4, [](int insd, faceType& face) -> void {
    face.eType = ElementType::QUD4;
    face.nG = 4;
    }
  },

  {6, [](int insd, faceType& face) -> void {
    face.eType = ElementType::TRI6;
    face.nG = 7;
    }
  },

  {8, [](int insd, faceType& face) -> void {
    face.eType = ElementType::QUD8;
    face.nG = 9;
    }
  },

  {9, [](int insd, faceType& face) -> void {
    face.eType = ElementType::QUD9;
    face.nG = 9;
    }
  },


};

