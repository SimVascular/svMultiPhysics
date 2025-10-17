// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/// brief Define a map type used to set element Gauss integration data.
//
using GetElementGausIntMapType = std::map<ElementType, std::function<void(const int, const int,
    Vector<double>&, Array<double>&)>>;

/// @brief Replicates 'SUBROUTINE GETGIP(insd, eType, nG, w, xi)' defined in NN.f.
//
GetElementGausIntMapType get_element_gauss_int_data = {

  {ElementType::HEX8, [](const int insd, const int nG, Vector<double>& w, Array<double>& xi) -> void {
    w = 1.0;
    double s = 1.0 / sqrt(3.0);
    double t = -1.0 / sqrt(3.0);

    xi(0,0) = t; 
    xi(1,0) = t; 
    xi(2,0) = t;

    xi(0,1) = s; 
    xi(1,1) = t; 
    xi(2,1) = t;

    xi(0,2) = s; 
    xi(1,2) = s; 
    xi(2,2) = t;

    xi(0,3) = t; 
    xi(1,3) = s; 
    xi(2,3) = t;

    xi(0,4) = t; 
    xi(1,4) = t; 
    xi(2,4) = s;

    xi(0,5) = s; 
    xi(1,5) = t; 
    xi(2,5) = s;

    xi(0,6) = s; 
    xi(1,6) = s; 
    xi(2,6) = s;

    xi(0,7) = t; 
    xi(1,7) = s; 
    xi(2,7) = s;
    }
  },

  {ElementType::HEX20, [](const int insd, const int nG, Vector<double>& w, Array<double>& xi) -> void {
    for (int i = 0; i < 8; i++) {
      w(i) = 125.0 / 729.0;
    }

    for (int i = 8; i < 20; i++) {
      w(i) = 200.0 / 729.0;
    }

    for (int i = 20; i < 26; i++) {
      w(i) = 320.0 / 729.0;
    }

    w(26) = 512.0 / 729.0;

    double s = sqrt(0.6);
    double  t = 0.0;

    int n = 0;
    xi(0, n) = -s; xi(1, n) = -s; xi(2, n) = -s; n++;
    xi(0, n) =  s; xi(1, n) = -s; xi(2, n) = -s; n++;
    xi(0, n) =  s; xi(1, n) =  s; xi(2, n) = -s; n++;
    xi(0, n) = -s; xi(1, n) =  s; xi(2, n) = -s; n++;
    xi(0, n) = -s; xi(1, n) = -s; xi(2, n) =  s; n++;
    xi(0, n) =  s; xi(1, n) = -s; xi(2, n) =  s; n++;
    xi(0, n) =  s; xi(1, n) =  s; xi(2, n) =  s; n++;
    xi(0, n) = -s; xi(1, n) =  s; xi(2, n) =  s; n++;

    xi(0, n) =  t; xi(1, n) = -s; xi(2, n) = -s; n++;
    xi(0,n) =  s; xi(1,n) =  t; xi(2,n) = -s; n++;
    xi(0,n) =  t; xi(1,n) =  s; xi(2,n) = -s; n++;
    xi(0,n) = -s; xi(1,n) =  t; xi(2,n) = -s; n++;
    xi(0,n) =  t; xi(1,n) = -s; xi(2,n) =  s; n++;
    xi(0,n) =  s; xi(1,n) =  t; xi(2,n) =  s; n++;
    xi(0,n) =  t; xi(1,n) =  s; xi(2,n) =  s; n++;
    xi(0,n) = -s; xi(1,n) =  t; xi(2,n) =  s; n++;
    xi(0,n) = -s; xi(1,n) = -s; xi(2,n) =  t; n++;
    xi(0,n) =  s; xi(1,n) = -s; xi(2,n) =  t; n++;
    xi(0,n) =  s; xi(1,n) =  s; xi(2,n) =  t; n++;
    xi(0,n) = -s; xi(1,n) =  s; xi(2,n) =  t; n++;

    xi(0,n) = -s; xi(1,n) =  t; xi(2,n) =  t; n++;
    xi(0,n) =  s; xi(1,n) =  t; xi(2,n) =  t; n++;
    xi(0,n) =  t; xi(1,n) = -s; xi(2,n) =  t; n++;
    xi(0,n) =  t; xi(1,n) =  s; xi(2,n) =  t; n++;
    xi(0,n) =  t; xi(1,n) =  t; xi(2,n) = -s; n++;
    xi(0,n) =  t; xi(1,n) =  t; xi(2,n) =  s; n++;

    xi(0,n) =  t; xi(1,n) =  t; xi(2,n) =  t;
    }
  },

  {ElementType::HEX27, [](const int insd, const int nG, Vector<double>& w, Array<double>& xi) -> void {

    for (int i = 0; i < 8; i++) {
      w(i) = 125.0 / 729.0;
    }

    for (int i = 8; i < 20; i++) {
      w(i) = 200.0 / 729.0;
    }

    for (int i = 20; i < 26; i++) {
      w(i) = 320.0 / 729.0;
    }

    w(26) = 512.0 / 729.0;

    double s = sqrt(0.6);
    double  t = 0.0;

    int n = 0;
    xi(0, n) = -s; xi(1, n) = -s; xi(2, n) = -s; n++;
    xi(0, n) =  s; xi(1, n) = -s; xi(2, n) = -s; n++;
    xi(0, n) =  s; xi(1, n) =  s; xi(2, n) = -s; n++;
    xi(0, n) = -s; xi(1, n) =  s; xi(2, n) = -s; n++;
    xi(0, n) = -s; xi(1, n) = -s; xi(2, n) =  s; n++;
    xi(0, n) =  s; xi(1, n) = -s; xi(2, n) =  s; n++;
    xi(0, n) =  s; xi(1, n) =  s; xi(2, n) =  s; n++;
    xi(0, n) = -s; xi(1, n) =  s; xi(2, n) =  s; n++;

    xi(0, n) =  t; xi(1, n) = -s; xi(2, n) = -s; n++;
    xi(0,n) =  s; xi(1,n) =  t; xi(2,n) = -s; n++;
    xi(0,n) =  t; xi(1,n) =  s; xi(2,n) = -s; n++;
    xi(0,n) = -s; xi(1,n) =  t; xi(2,n) = -s; n++;
    xi(0,n) =  t; xi(1,n) = -s; xi(2,n) =  s; n++;
    xi(0,n) =  s; xi(1,n) =  t; xi(2,n) =  s; n++;
    xi(0,n) =  t; xi(1,n) =  s; xi(2,n) =  s; n++;
    xi(0,n) = -s; xi(1,n) =  t; xi(2,n) =  s; n++;
    xi(0,n) = -s; xi(1,n) = -s; xi(2,n) =  t; n++;
    xi(0,n) =  s; xi(1,n) = -s; xi(2,n) =  t; n++;
    xi(0,n) =  s; xi(1,n) =  s; xi(2,n) =  t; n++;
    xi(0,n) = -s; xi(1,n) =  s; xi(2,n) =  t; n++;

    xi(0,n) = -s; xi(1,n) =  t; xi(2,n) =  t; n++;
    xi(0,n) =  s; xi(1,n) =  t; xi(2,n) =  t; n++;
    xi(0,n) =  t; xi(1,n) = -s; xi(2,n) =  t; n++;
    xi(0,n) =  t; xi(1,n) =  s; xi(2,n) =  t; n++;
    xi(0,n) =  t; xi(1,n) =  t; xi(2,n) = -s; n++;
    xi(0,n) =  t; xi(1,n) =  t; xi(2,n) =  s; n++;

    xi(0,n) =  t; xi(1,n) =  t; xi(2,n) =  t;
    }
  },

  {ElementType::LIN1, [](const int insd, const int nG, Vector<double>& w, Array<double>& xi) -> void {
    w = 1.0; 
    double s = 1.0 / sqrt(3.0);
    xi(0,0) = -s;
    xi(0,1) =  s;
    }
  },

  {ElementType::QUD4, [](const int insd, const int nG, Vector<double>& w, Array<double>& xi) -> void {
    w = 1.0;
    double s = 1.0 / sqrt(3.0);
    xi(0,0) = -s; xi(1,0) = -s;
    xi(0,1) =  s; xi(1,1) = -s;
    xi(0,2) =  s; xi(1,2) =  s;
    xi(0,3) = -s; xi(1,3) =  s;
    }
  },

  {ElementType::QUD9, [](const int insd, const int nG, Vector<double>& w, Array<double>& xi) -> void {
     for (int i = 0; i < 4; i++) {
        w(i) = 25.0 / 81.0; 
        w(i+4) = 40.0 / 81.0; 
      }
      w(8) = 64.0 / 81.0; 

      double s = sqrt(6.0);
      double t = 0.0;

      xi(0,0) = -s; 
      xi(1,0) = -s;

      xi(0,1) =  s; 
      xi(1,1) = -s;

      xi(0,2) =  s; 
      xi(1,2) =  s;

      xi(0,3) = -s; 
      xi(1,3) =  s;

      xi(0,4) =  t; 
      xi(1,4) = -s;

      xi(0,5) =  s; 
      xi(1,5) =  t;

      xi(0,6) =  t; 
      xi(1,6) =  s;

      xi(0,7) = -s; 
      xi(1,7) =  t;

      xi(0,8) =  t; 
      xi(1,8) =  t;
    }
  },

  {ElementType::TET4, [](const int insd, const int nG, Vector<double>& w, Array<double>& xi) -> void {
     for (int i = 0; i < nG; i++) {
        w(i) = 1.0 / 24.0;
      }
      double s = (5.0 + 3.0*sqrt(5.0)) / 20.0;
      double t = (5.0 - sqrt(5.0)) / 20.0;
      xi(0,0) = s; xi(1,0) = t; xi(2,0) = t;
      xi(0,1) = t; xi(1,1) = s; xi(2,1) = t;
      xi(0,2) = t; xi(1,2) = t; xi(2,2) = s;
      xi(0,3) = t; xi(1,3) = t; xi(2,3) = t;
    }
  },

  {ElementType::TRI3, [](const int insd, const int nG, Vector<double>& w, Array<double>& xi) -> void {
     for (int i = 0; i < nG; i++) {
        w(i) = 1.0 / 6.0;
      }
    double s = 2.0/3.0;
    double t = 1.0/6.0;
    xi(0,0) = t; 
    xi(1,0) = t;
    xi(0,1) = s; 
    xi(1,1) = t;
    xi(0,2) = t; 
    xi(1,2) = s;
    }
  },

  {ElementType::WDG, [](const int insd, const int nG, Vector<double>& w, Array<double>& xi)  -> void { 
     for (int i = 0; i < nG; i++) { 
       w(i) = 1.0 / 6.0;
     }
    double s  =  2.0 / 3.0;
    double t  =  1.0 / 6.0;
    double uz =  1.0 / sqrt(3.0);
    double lz = -1.0 / sqrt(3.0);
    xi(0,0) = s; xi(1,0) = t; xi(2,0) = lz;
    xi(0,1) = t; xi(1,1) = s; xi(2,1) = lz;
    xi(0,2) = t; xi(1,2) = t; xi(2,2) = lz;
    xi(0,3) = s; xi(1,3) = t; xi(2,3) = uz;
    xi(0,4) = t; xi(1,4) = s; xi(2,4) = uz;
    xi(0,5) = t; xi(1,5) = t; xi(2,5) = uz;
    }
  },
  
};


/// @brief Define a map type used to set element Gauss integration data.
using SetElementGausIntMapType = std::map<ElementType, std::function<void(mshType&)>>;

/// @brief Replicates 'SUBROUTINE GETGIP(insd, eType, nG, w, xi)' defined in NN.f.
///
/// \todo [TODO:DaveP] I screwed up here, should just have a single map for mesh and face types.
//
SetElementGausIntMapType set_element_gauss_int_data = {

  {ElementType::HEX8, [](mshType& mesh) -> void { 
      for (int i = 0; i < mesh.nG; i++) { 
        mesh.w(i) =  1.0;
      }
      double s =  1.0 / sqrt(3.0);
      double t = -1.0 / sqrt(3.0);
      mesh.xi(0,0) = t; mesh.xi(1,0) = t; mesh.xi(2,0) = t;
      mesh.xi(0,1) = s; mesh.xi(1,1) = t; mesh.xi(2,1) = t;
      mesh.xi(0,2) = s; mesh.xi(1,2) = s; mesh.xi(2,2) = t;
      mesh.xi(0,3) = t; mesh.xi(1,3) = s; mesh.xi(2,3) = t;
      mesh.xi(0,4) = t; mesh.xi(1,4) = t; mesh.xi(2,4) = s;
      mesh.xi(0,5) = s; mesh.xi(1,5) = t; mesh.xi(2,5) = s;
      mesh.xi(0,6) = s; mesh.xi(1,6) = s; mesh.xi(2,6) = s;
      mesh.xi(0,7) = t; mesh.xi(1,7) = s; mesh.xi(2,7) = s;
    }
  },

  {ElementType::HEX20, [](mshType& mesh) -> void { 

    for (int i = 0; i < 8; i++) {
      mesh.w(i) = 125.0 / 729.0;
    }

    for (int i = 8; i < 20; i++) {
      mesh.w(i) = 200.0 / 729.0;
    }

    for (int i = 20; i < 26; i++) {
      mesh.w(i) = 320.0 / 729.0;
    }

    mesh.w(26) = 512.0 / 729.0;

    double s = sqrt(0.6);
    double  t = 0.0;

    int n = 0;
    mesh.xi(0, n) = -s; mesh.xi(1, n) = -s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0, n) =  s; mesh.xi(1, n) = -s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0, n) =  s; mesh.xi(1, n) =  s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0, n) = -s; mesh.xi(1, n) =  s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0, n) = -s; mesh.xi(1, n) = -s; mesh.xi(2, n) =  s; n++;
    mesh.xi(0, n) =  s; mesh.xi(1, n) = -s; mesh.xi(2, n) =  s; n++;
    mesh.xi(0, n) =  s; mesh.xi(1, n) =  s; mesh.xi(2, n) =  s; n++;
    mesh.xi(0, n) = -s; mesh.xi(1, n) =  s; mesh.xi(2, n) =  s; n++;

    mesh.xi(0, n) =  t; mesh.xi(1, n) = -s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) =  t; mesh.xi(2,n) = -s; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  s; mesh.xi(2,n) = -s; n++;
    mesh.xi(0,n) = -s; mesh.xi(1,n) =  t; mesh.xi(2,n) = -s; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) = -s; mesh.xi(2,n) =  s; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) =  t; mesh.xi(2,n) =  s; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  s; mesh.xi(2,n) =  s; n++;
    mesh.xi(0,n) = -s; mesh.xi(1,n) =  t; mesh.xi(2,n) =  s; n++;
    mesh.xi(0,n) = -s; mesh.xi(1,n) = -s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) = -s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) =  s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) = -s; mesh.xi(1,n) =  s; mesh.xi(2,n) =  t; n++;

    mesh.xi(0,n) = -s; mesh.xi(1,n) =  t; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) =  t; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) = -s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  t; mesh.xi(2,n) = -s; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  t; mesh.xi(2,n) =  s; n++;

    mesh.xi(0,n) =  t; mesh.xi(1,n) =  t; mesh.xi(2,n) =  t;
    }
  },

  {ElementType::HEX27, [](mshType& mesh) -> void { 

    for (int i = 0; i < 8; i++) {
      mesh.w(i) = 125.0 / 729.0;
    }

    for (int i = 8; i < 20; i++) {
      mesh.w(i) = 200.0 / 729.0;
    }

    for (int i = 20; i < 26; i++) {
      mesh.w(i) = 320.0 / 729.0;
    }

    mesh.w(26) = 512.0 / 729.0;

    double s = sqrt(0.6);
    double  t = 0.0;

    int n = 0;
    mesh.xi(0, n) = -s; mesh.xi(1, n) = -s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0, n) =  s; mesh.xi(1, n) = -s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0, n) =  s; mesh.xi(1, n) =  s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0, n) = -s; mesh.xi(1, n) =  s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0, n) = -s; mesh.xi(1, n) = -s; mesh.xi(2, n) =  s; n++;
    mesh.xi(0, n) =  s; mesh.xi(1, n) = -s; mesh.xi(2, n) =  s; n++;
    mesh.xi(0, n) =  s; mesh.xi(1, n) =  s; mesh.xi(2, n) =  s; n++;
    mesh.xi(0, n) = -s; mesh.xi(1, n) =  s; mesh.xi(2, n) =  s; n++;

    mesh.xi(0, n) =  t; mesh.xi(1, n) = -s; mesh.xi(2, n) = -s; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) =  t; mesh.xi(2,n) = -s; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  s; mesh.xi(2,n) = -s; n++;
    mesh.xi(0,n) = -s; mesh.xi(1,n) =  t; mesh.xi(2,n) = -s; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) = -s; mesh.xi(2,n) =  s; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) =  t; mesh.xi(2,n) =  s; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  s; mesh.xi(2,n) =  s; n++;
    mesh.xi(0,n) = -s; mesh.xi(1,n) =  t; mesh.xi(2,n) =  s; n++;
    mesh.xi(0,n) = -s; mesh.xi(1,n) = -s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) = -s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) =  s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) = -s; mesh.xi(1,n) =  s; mesh.xi(2,n) =  t; n++;

    mesh.xi(0,n) = -s; mesh.xi(1,n) =  t; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  s; mesh.xi(1,n) =  t; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) = -s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  s; mesh.xi(2,n) =  t; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  t; mesh.xi(2,n) = -s; n++;
    mesh.xi(0,n) =  t; mesh.xi(1,n) =  t; mesh.xi(2,n) =  s; n++;

    mesh.xi(0,n) =  t; mesh.xi(1,n) =  t; mesh.xi(2,n) =  t;
    }
  },

  {ElementType::LIN1, [](mshType& mesh) -> void { 
      for (int i = 0; i < mesh.nG; i++) { 
        mesh.w(i) =  1.0;
      }
      double s = 1.0 / sqrt(3.0);
      mesh.xi(0,0) = -s;
      mesh.xi(0,1) =  s;
    }
  },

  {ElementType::LIN2, [](mshType& mesh) -> void {
      mesh.w[0] = 5.0 / 9.0;
      mesh.w[1] = 5.0 / 9.0;
      mesh.w[2] = 8.0 / 9.0;

      double s = sqrt(0.6);

      mesh.xi(0,0) = -s;
      mesh.xi(0,1) =  s;
      mesh.xi(0,2) = 0.0;
    }
  },

  {ElementType::QUD4, [](mshType& mesh) -> void {
     for (int i = 0; i < mesh.nG; i++) {
       mesh.w[i] = 1.0;
     }
    double s = 1.0 / sqrt(3.0);
    mesh.xi(0,0) = -s; mesh.xi(1,0) = -s;
    mesh.xi(0,1) =  s; mesh.xi(1,1) = -s;
    mesh.xi(0,2) =  s; mesh.xi(1,2) =  s;
    mesh.xi(0,3) = -s; mesh.xi(1,3) =  s;
    }
  },

  {ElementType::QUD9, [](mshType& mesh) -> void {
    mesh.w(0) = 25.0 / 81.0;
    mesh.w(1) = 25.0 / 81.0;
    mesh.w(2) = 25.0 / 81.0;
    mesh.w(3) = 25.0 / 81.0;

    mesh.w(4) = 40.0 / 81.0;
    mesh.w(5) = 40.0 / 81.0;
    mesh.w(6) = 40.0 / 81.0;
    mesh.w(7) = 40.0 / 81.0;
    mesh.w(8) = 64.0 / 81.0;

    double s = sqrt(0.6);
    double t = 0.0;

    mesh.xi(0,0) = -s; 
    mesh.xi(1,0) = -s;

    mesh.xi(0,1) =  s; 
    mesh.xi(1,1) = -s;

    mesh.xi(0,2) =  s; 
    mesh.xi(1,2) =  s;

    mesh.xi(0,3) = -s; 
    mesh.xi(1,3) =  s;

    mesh.xi(0,4) =  t; 
    mesh.xi(1,4) = -s;

    mesh.xi(0,5) =  s; 
    mesh.xi(1,5) =  t;

    mesh.xi(0,6) =  t; 
    mesh.xi(1,6) =  s;

    mesh.xi(0,7) = -s; 
    mesh.xi(1,7) =  t;

    mesh.xi(0,8) =  t; 
    mesh.xi(1,8) =  t;
    }
  },


  {ElementType::TET4, [](mshType& mesh) -> void { 
      for (int i = 0; i < mesh.nG; i++) { 
        mesh.w(i) = 1.0 / 24.0;
      }
      // s=0.25: central quadrature rule
      // s=(5.0+3.0*sqrt(5.0))/20.0: original gaussian quadrature rule
      // s=1.0: nodal quadrature
      double s = mesh.qmTET4;
      // t=0.25: central quadrature rule
      // t=(5.0-sqrt(5.0))/20.0: original gaussian quadrature rule
      // t=0.0: nodal quadrature
      double t = (1-s)/3.0;
      mesh.xi(0,0) = s; mesh.xi(1,0) = t; mesh.xi(2,0) = t;
      mesh.xi(0,1) = t; mesh.xi(1,1) = s; mesh.xi(2,1) = t;
      mesh.xi(0,2) = t; mesh.xi(1,2) = t; mesh.xi(2,2) = s;
      mesh.xi(0,3) = t; mesh.xi(1,3) = t; mesh.xi(2,3) = t;
    }
  },

  {ElementType::TET10, [](mshType& mesh) -> void {
      mesh.w(0) = 0.0302836780970890;

      mesh.w(1) = 0.0060267857142860;
      mesh.w(2) = 0.0060267857142860;
      mesh.w(3) = 0.0060267857142860;
      mesh.w(4) = 0.0060267857142860;

      mesh.w(5) = 0.0116452490860290;
      mesh.w(6) = 0.0116452490860290;
      mesh.w(7) = 0.0116452490860290;
      mesh.w(8) = 0.0116452490860290;

      mesh.w(9) = 0.0109491415613860;
      mesh.w(10) = 0.0109491415613860;
      mesh.w(11) = 0.0109491415613860;
      mesh.w(12) = 0.0109491415613860;
      mesh.w(13) = 0.0109491415613860;
      mesh.w(14) = 0.0109491415613860;

      double s = 0.250;
      auto& xi = mesh.xi;
      xi(0,0) = s; xi(1,0) = s; xi(2,0) = s;

      s = 0.3333333333333330;
      double t = 0.0;
      xi(0,1) = t; xi(1,1) = s; xi(2,1) = s;
      xi(0,2) = s; xi(1,2) = t; xi(2,2) = s;
      xi(0,3) = s; xi(1,3) = s; xi(2,3) = t;
      xi(0,4) = s; xi(1,4) = s; xi(2,4) = s;

      s = 0.0909090909090910;
      t = 0.7272727272727270;
      xi(0,5) = t; xi(1,5) = s; xi(2,5) = s;
      xi(0,6) = s; xi(1,6) = t; xi(2,6) = s;
      xi(0,7) = s; xi(1,7) = s; xi(2,7) = t;
      xi(0,8) = s; xi(1,8) = s; xi(2,8) = s;

      s = 0.0665501535736640;
      t = 0.4334498464263360;
      xi(0, 9) = s; xi(1, 9) = s; xi(2, 9) = t;
      xi(0,10) = s; xi(1,10) = t; xi(2,10) = s;
      xi(0,11) = s; xi(1,11) = t; xi(2,11) = t;
      xi(0,12) = t; xi(1,12) = t; xi(2,12) = s;
      xi(0,13) = t; xi(1,13) = s; xi(2,13) = t;
      xi(0,14) = t; xi(1,14) = s; xi(2,14) = s;
    }
  },

  {ElementType::TRI3, [](mshType& mesh) -> void { 
     for (int i = 0; i < mesh.nG; i++) { 
       mesh.w(i) = 1.0 / 6.0;
     }
    double s = 2.0 / 3.0;
    double t = 1.0 / 6.0;
    mesh.xi(0,0) = t; mesh.xi(1,0) = t;
    mesh.xi(0,1) = s; mesh.xi(1,1) = t;
    mesh.xi(0,2) = t; mesh.xi(1,2) = s;
    }
  },

  {ElementType::TRI6, [](mshType& mesh) -> void {
    mesh.w[0] = 0.225000000000000 * 5.0e-1 ;
    mesh.w[1] = 0.125939180544827 * 5.0e-1 ;
    mesh.w[2] = 0.125939180544827 * 5.0e-1 ;
    mesh.w[3] = 0.125939180544827 * 5.0e-1 ;
    mesh.w(4) = 0.132394152788506 * 5.0e-1 ;
    mesh.w(5) = 0.132394152788506 * 5.0e-1 ;
    mesh.w(6) = 0.132394152788506 * 5.0e-1 ;

    double s = 0.333333333333333;
    mesh.xi(0,0) = s; mesh.xi(1,0) = s;

    s = 0.797426985353087;
    double t = 0.101286507323456;
    mesh.xi(0,1) = s; mesh.xi(1,1) = t;
    mesh.xi(0,2) = t; mesh.xi(1,2) = s;
    mesh.xi(0,3) = t; mesh.xi(1,3) = t;

    s = 0.059715871789770;
    t = 0.470142064105115;
    mesh.xi(0,4) = s; mesh.xi(1,4) = t;
    mesh.xi(0,5) = t; mesh.xi(1,5) = s;
    mesh.xi(0,6) = t; mesh.xi(1,6) = t;
    }
  },

  {ElementType::WDG, [](mshType& mesh) -> void { 
     for (int i = 0; i < mesh.nG; i++) { 
       mesh.w(i) = 1.0 / 6.0;
     }
    double s  =  2.0 / 3.0;
    double t  =  1.0 / 6.0;
    double uz =  1.0 / sqrt(3.0);
    double lz = -1.0 / sqrt(3.0);
    mesh.xi(0,0) = s; mesh.xi(1,0) = t; mesh.xi(2,0) = lz;
    mesh.xi(0,1) = t; mesh.xi(1,1) = s; mesh.xi(2,1) = lz;
    mesh.xi(0,2) = t; mesh.xi(1,2) = t; mesh.xi(2,2) = lz;
    mesh.xi(0,3) = s; mesh.xi(1,3) = t; mesh.xi(2,3) = uz;
    mesh.xi(0,4) = t; mesh.xi(1,4) = s; mesh.xi(2,4) = uz;
    mesh.xi(0,5) = t; mesh.xi(1,5) = t; mesh.xi(2,5) = uz;
    }
  },

};

/// @brief Define a map type used to set face Gauss integration data.
//
using SetFaceGausIntMapType = std::map<ElementType, std::function<void(faceType&)>>;

/// \todo [TODO:DaveP] I screwed up here, should just have a single map for mesh and face types.
//
SetFaceGausIntMapType set_face_gauss_int_data = {

  {ElementType::PNT, [](faceType& face) -> void {
      for (int i = 0; i < face.nG; i++) {
        face.w(i) =  1.0;
      }
    }
  },

  {ElementType::LIN1, [](faceType& face) -> void {
      for (int i = 0; i < face.nG; i++) {
        face.w(i) =  1.0;
      }
      double s = 1.0 / sqrt(3.0);
      face.xi(0,0) = -s;
      face.xi(0,1) =  s;
    }
  },

  {ElementType::LIN2, [](faceType& face) -> void {
      face.w[0] = 5.0 / 9.0;
      face.w[1] = 5.0 / 9.0;
      face.w[2] = 8.0 / 9.0;

      double s = sqrt(0.6);

      face.xi(0,0) = -s;
      face.xi(0,1) =  s;
      face.xi(0,2) = 0.0;
    }
  },

  {ElementType::QUD4, [](faceType& face) -> void {
      face.w = 1.0;
      double s =  1.0 / sqrt(3.0);
      face.xi(0,0) = -s; face.xi(1,0) = -s;
      face.xi(0,1) =  s; face.xi(1,1) = -s;
      face.xi(0,2) =  s; face.xi(1,2) =  s;
      face.xi(0,3) = -s; face.xi(1,3) =  s;
    }
  },

  {ElementType::QUD8, [](faceType& face) -> void {
      face.w(0) = 25.0/81.0;
      face.w(1) = 25.0/81.0;
      face.w(2) = 25.0/81.0;
      face.w(3) = 25.0/81.0;

      face.w(4) = 40.0/81.0;
      face.w(5) = 40.0/81.0;
      face.w(6) = 40.0/81.0;
      face.w(7) = 40.0/81.0;

      face.w(8) = 64.0/81.0;

      double s = sqrt(0.6);
      double t = 0.0;

      face.xi(0,0) = -s; face.xi(1,0) = -s;
      face.xi(0,1) =  s; face.xi(1,1) = -s;
      face.xi(0,2) =  s; face.xi(1,2) =  s;
      face.xi(0,3) = -s; face.xi(1,3) =  s;
      face.xi(0,4) =  t; face.xi(1,4) = -s;
      face.xi(0,5) =  s; face.xi(1,5) =  t;
      face.xi(0,6) =  t; face.xi(1,6) =  s;
      face.xi(0,7) = -s; face.xi(1,7) =  t;
      face.xi(0,8) =  t; face.xi(1,8) =  t;
    }
  },

  {ElementType::QUD9, [](faceType& face) -> void {
    face.w(0) = 25.0 / 81.0;
    face.w(1) = 25.0 / 81.0;
    face.w(2) = 25.0 / 81.0;
    face.w(3) = 25.0 / 81.0;

    face.w(4) = 40.0 / 81.0;
    face.w(5) = 40.0 / 81.0;
    face.w(6) = 40.0 / 81.0;
    face.w(7) = 40.0 / 81.0;
    face.w(8) = 64.0 / 81.0;

    double s = sqrt(0.6);
    double t = 0.0;

    face.xi(0,0) = -s; 
    face.xi(1,0) = -s;

    face.xi(0,1) =  s; 
    face.xi(1,1) = -s;

    face.xi(0,2) =  s; 
    face.xi(1,2) =  s;

    face.xi(0,3) = -s; 
    face.xi(1,3) =  s;

    face.xi(0,4) =  t; 
    face.xi(1,4) = -s;

    face.xi(0,5) =  s; 
    face.xi(1,5) =  t;

    face.xi(0,6) =  t; 
    face.xi(1,6) =  s;

    face.xi(0,7) = -s; 
    face.xi(1,7) =  t;

    face.xi(0,8) =  t; 
    face.xi(1,8) =  t;
    }
  },


  {ElementType::TRI3, [](faceType& face) -> void {
      //for (int i = 0; i < face.nG; i++) {
      //  face.w(i) =  1.0 / 6.0;
      //}
      face.w =  1.0 / 6.0;

      // s=1.0/3.0: central quadrature rule
      // s=2.0/3.0: original gaussian quadrature rule
      // s=1.0: nodal quadrature
      double s = face.qmTRI3;
      // t=1.0/3.0: central quadrature rule
      // t=1.0/6.0: original gaussian quadrature rule
      // t=0.0: nodal quadrature
      double t = -0.5*s+0.5;
      face.xi(0,0) = t; face.xi(1,0) = t;
      face.xi(0,1) = s; face.xi(1,1) = t;
      face.xi(0,2) = t; face.xi(1,2) = s;
    }
  },

  {ElementType::TRI6, [](faceType& face) -> void {
    face.w[0] = 0.225000000000000 * 5.0e-1;
    face.w[1] = 0.125939180544827 * 5.0e-1;
    face.w[2] = 0.125939180544827 * 5.0e-1;
    face.w[3] = 0.125939180544827 * 5.0e-1 ;
    face.w(4) = 0.132394152788506 * 5.0e-1 ;
    face.w(5) = 0.132394152788506 * 5.0e-1 ;
    face.w(6) = 0.132394152788506 * 5.0e-1 ;

    double s = 0.333333333333333;
    face.xi(0,0) = s; face.xi(1,0) = s;

    s = 0.797426985353087;
    double t = 0.101286507323456;
    face.xi(0,1) = s; face.xi(1,1) = t;
    face.xi(0,2) = t; face.xi(1,2) = s;
    face.xi(0,3) = t; face.xi(1,3) = t;

    s = 0.059715871789770;
    t = 0.470142064105115;
    face.xi(0,4) = s; face.xi(1,4) = t;
    face.xi(0,5) = t; face.xi(1,5) = s;
    face.xi(0,6) = t; face.xi(1,6) = t;
    }
  },




};
