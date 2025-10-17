// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/// @brief Define a map type used to compute 2nd direivatives of element shape function data.
///
/// Replicates 'SUBROUTINE GETGNNxx(insd, ind2, eType, eNoN, xi, Nxx)'
//
static double fp = 4.0;
static double fn = -4.0;
static double en = -8.0;
static double ze =  0.0;

using GetElement2ndDerivMapType = std::map<ElementType, std::function<void(const int, const int, const int, 
    const int, const Array<double>&, Array3<double>&)>>;

GetElement2ndDerivMapType get_element_2nd_derivs = {

  {ElementType::QUD8, [](const int insd, const int ind2, const int eNoN, const int g, const Array<double>& xi, 
       Array3<double>& Nxx) -> void {

    double lx = 1.0 - xi(0);
    double ly = 1.0 - xi(1);
    double ux = 1.0 + xi(0);
    double uy = 1.0 + xi(1);
    double mx = xi(0);
    double my = xi(1);

    Nxx(0,0,g) =  ly*0.50;
    Nxx(1,0,g) =  lx*0.50;
    Nxx(2,0,g) =  (lx+lx+ly+ly-3.0)/4.0;

    Nxx(0,1,g) =  ly*0.50;
    Nxx(1,1,g) =  ux*0.50;
    Nxx(2,1,g) = -(ux+ux+ly+ly-3.0)/4.0;

    Nxx(0,2,g) =  uy*0.50;
    Nxx(1,2,g) =  ux*0.50;
    Nxx(2,3,g) =  (ux+ux+uy+uy-3.0)/4.0;

    Nxx(0,3,g) =  uy*0.50;
    Nxx(1,3,g) =  lx*0.50;
    Nxx(2,3,g) = -(lx+lx+uy+uy-3.0)/4.0;

    Nxx(0,4,g) = -ly;
    Nxx(1,4,g) =  0.0;
    Nxx(2,4,g) =  mx;

    Nxx(0,5,g) =  0.0;
    Nxx(1,5,g) = -ux;
    Nxx(2,5,g) = -my;

    Nxx(0,6,g) = -uy;
    Nxx(1,6,g) =  0.0;
    Nxx(2,6,g) = -mx;

    Nxx(0,7,g) =  0.0;
    Nxx(1,7,g) = -lx;
    Nxx(2,7,g) =  my;
    }
  },

  {ElementType::QUD9, [](const int insd, const int ind2, const int eNoN, const int g, const Array<double>& xi, 
       Array3<double>& Nxx) -> void { 

    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double mx = xi(0,g);
    double my = xi(1,g);

    Nxx(0,0,g) = -ly*my*0.5;
    Nxx(1,0,g) = -lx*mx*0.5;
    Nxx(2,0,g) =  (lx-mx)*(ly-my)/4.0;

    Nxx(0,1,g) = -ly*my*0.5;
    Nxx(1,1,g) =  ux*mx*0.5;
    Nxx(2,1,g) = -(ux+mx)*(ly-my)/4.0;

    Nxx(0,2,g) =  uy*my*0.5;
    Nxx(1,2,g) =  ux*mx*0.5;
    Nxx(2,2,g) =  (ux+mx)*(uy+my)/4.0;

    Nxx(0,3,g) =  uy*my*0.5;
    Nxx(1,3,g) = -lx*mx*0.5;
    Nxx(2,3,g) = -(lx-mx)*(uy+my)/4.0;

    Nxx(0,4,g) =  ly*my;
    Nxx(1,4,g) =  lx*ux;
    Nxx(2,4,g) =  mx*(ly-my);

    Nxx(0,5,g) =  ly*uy;
    Nxx(1,5,g) = -ux*mx;
    Nxx(2,5,g) = -(ux+mx)*my;

    Nxx(0,6,g) = -uy*my;
    Nxx(1,6,g) =  lx*ux;
    Nxx(2,6,g) = -mx*(uy+my);

    Nxx(0,7,g) =  ly*uy;
    Nxx(1,7,g) =  lx*mx;
    Nxx(2,7,g) =  (lx-mx)*my;

    Nxx(0,8,g) = -ly*uy*2.0;
    Nxx(1,8,g) = -lx*ux*2.0;
    Nxx(2,8,g) =  mx*my*4.0;
    }
  },

  {ElementType::TET10, [](const int insd, const int ind2, const int eNoN, const int g, const Array<double>& xi, 
       Array3<double>& Nxx) -> void { 
    Nxx.set_row(0, g, {fp, ze, ze, ze, ze, ze});
    Nxx.set_row(1, g, {ze, fp, ze, ze, ze, ze});
    Nxx.set_row(2, g, {ze, ze, fp, ze, ze, ze});
    Nxx.set_row(3, g, {fp, fp, fp, fp, fp, fp});
    Nxx.set_row(4, g, {ze, ze, ze, fp, ze, ze});
    Nxx.set_row(5, g, {ze, ze, ze, ze, fp, ze});
    Nxx.set_row(6, g, {ze, ze, ze, ze, ze, fp});
    Nxx.set_row(7, g, {en, ze, ze, fn, ze, fn});
    Nxx.set_row(8, g, {ze, en, ze, fn, fn, ze});
    Nxx.set_row(9, g, {ze, ze, en, ze, fn, fn});
    }
  },

  {ElementType::TRI6, [](const int insd, const int ind2, const int eNoN, const int g, const Array<double>& xi,
       Array3<double>& Nxx) -> void {

    Nxx.set_row(0, g, {fp, ze, ze});
    Nxx.set_row(1, g, {ze, fp, ze});
    Nxx.set_row(2, g, {fp, fp, fp});
    Nxx.set_row(3, g, {ze, ze, fp});
    Nxx.set_row(4, g, {ze, en, fn});
    Nxx.set_row(5, g, {en, ze, fn});
    }
  },

};


