// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/// @brief Define a map type used to set element shape function data.
///
/// Reproduces the Fortran 'GETGNN' subroutine.
//
using GetElementShapeMapType = std::map<ElementType, std::function<void(const int, const int, const int, 
    Array<double>&, Array<double>&, Array3<double>&)>>;

GetElementShapeMapType get_element_shape_data = {

  {ElementType::HEX8, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double lz = 1.0 - xi(2,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double uz = 1.0 + xi(2,g);

    N(0,g) = lx*ly*lz/8.0;
    N(1,g) = ux*ly*lz/8.0;
    N(2,g) = ux*uy*lz/8.0;
    N(3,g) = lx*uy*lz/8.0;
    N(4,g) = lx*ly*uz/8.0;
    N(5,g) = ux*ly*uz/8.0;
    N(6,g) = ux*uy*uz/8.0;
    N(7,g) = lx*uy*uz/8.0;

    Nx(0,0,g) = -ly*lz/8.0;
    Nx(1,0,g) = -lx*lz/8.0;
    Nx(2,0,g) = -lx*ly/8.0;

    Nx(0,1,g) =  ly*lz/8.0;
    Nx(1,1,g) = -ux*lz/8.0;
    Nx(2,1,g) = -ux*ly/8.0;

    Nx(0,2,g) =  uy*lz/8.0;
    Nx(1,2,g) =  ux*lz/8.0;
    Nx(2,2,g) = -ux*uy/8.0;

    Nx(0,3,g) = -uy*lz/8.0;
    Nx(1,3,g) =  lx*lz/8.0;
    Nx(2,3,g) = -lx*uy/8.0;

    Nx(0,4,g) = -ly*uz/8.0;
    Nx(1,4,g) = -lx*uz/8.0;
    Nx(2,4,g) =  lx*ly/8.0;

    Nx(0,5,g) =  ly*uz/8.0;
    Nx(1,5,g) = -ux*uz/8.0;
    Nx(2,5,g) =  ux*ly/8.0;

    Nx(0,6,g) =  uy*uz/8.0;
    Nx(1,6,g) =  ux*uz/8.0;
    Nx(2,6,g) =  ux*uy/8.0;

    Nx(0,7,g) = -uy*uz/8.0;
    Nx(1,7,g) =  lx*uz/8.0;
    Nx(2,7,g) =  lx*uy/8.0;
    }
  },

  {ElementType::HEX20, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N,
      Array3<double>& Nx) -> void
    {
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double lz = 1.0 - xi(2,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double uz = 1.0 + xi(2,g);

    double mx = lx*ux;
    double my = ly*uy;
    double mz = lz*uz;

    N(0, g) = lx*ly*lz*(lx+ly+lz-5.0)/8.0;
    N(1, g) = ux*ly*lz*(ux+ly+lz-5.0)/8.0;
    N(2, g) = ux*uy*lz*(ux+uy+lz-5.0)/8.0;
    N(3, g) = lx*uy*lz*(lx+uy+lz-5.0)/8.0;
    N(4, g) = lx*ly*uz*(lx+ly+uz-5.0)/8.0;
    N(5, g) = ux*ly*uz*(ux+ly+uz-5.0)/8.0;
    N(6, g) = ux*uy*uz*(ux+uy+uz-5.0)/8.0;
    N(7, g) = lx*uy*uz*(lx+uy+uz-5.0)/8.0;
    N(8, g) = mx*ly*lz/4.0;
    N(9, g) = ux*my*lz/4.0;
    N(10, g) = mx*uy*lz/4.0;
    N(11, g) = lx*my*lz/4.0;
    N(12, g) = mx*ly*uz/4.0;
    N(13, g) = ux*my*uz/4.0;
    N(14, g) = mx*uy*uz/4.0;
    N(15, g) = lx*my*uz/4.0;
    N(16, g) = lx*ly*mz/4.0;
    N(17, g) = ux*ly*mz/4.0;
    N(18, g) = ux*uy*mz/4.0;
    N(19, g) = lx*uy*mz/4.0;

    // N(1)  = lx*ly*lz*(lx+ly+lz-5.0)/8.0;
    int n = 0;
    Nx(0,n,g) = -ly*lz*(lx+ly+lz-5.0+lx)/8.0;
    Nx(1,n,g) = -lx*lz*(lx+ly+lz-5.0+ly)/8.0;
    Nx(2,n,g) = -lx*ly*(lx+ly+lz-5.0+lz)/8.0;

//c   N(n,g) = ux*ly*lz*(ux+ly+lz-5.0)/8.0;
    n += 1;
    Nx(0,n,g) =  ly*lz*(ux+ly+lz-5.0+ux)/8.0;
    Nx(1,n,g) = -ux*lz*(ux+ly+lz-5.0+ly)/8.0;
    Nx(2,n,g) = -ux*ly*(ux+ly+lz-5.0+lz)/8.0;

//c   N(n,g) = ux*uy*lz*(ux+uy+lz-5.0)/8.0
    n += 1;
    Nx(0,n,g) =  uy*lz*(ux+uy+lz-5.0+ux)/8.0;
    Nx(1,n,g) =  ux*lz*(ux+uy+lz-5.0+uy)/8.0;
    Nx(2,n,g) = -ux*uy*(ux+uy+lz-5.0+lz)/8.0;

//c   N(n,g) = lx*uy*lz*(lx+uy+lz-5.0)/8.0
    n += 1;
    Nx(0,n,g) = -uy*lz*(lx+uy+lz-5.0+lx)/8.0;
    Nx(1,n,g) =  lx*lz*(lx+uy+lz-5.0+uy)/8.0;
    Nx(2,n,g) = -lx*uy*(lx+uy+lz-5.0+lz)/8.0;

//c   N(n,g) = lx*ly*uz*(lx+ly+uz-5.0)/8.0
    n += 1;
    Nx(0,n,g) = -ly*uz*(lx+ly+uz-5.0+lx)/8.0;
    Nx(1,n,g) = -lx*uz*(lx+ly+uz-5.0+ly)/8.0;
    Nx(2,n,g) =  lx*ly*(lx+ly+uz-5.0+uz)/8.0;

//c   N(n,g) = ux*ly*uz*(ux+ly+uz-5.0)/8.0
    n += 1;
    Nx(0,n,g) =  ly*uz*(ux+ly+uz-5.0+ux)/8.0;
    Nx(1,n,g) = -ux*uz*(ux+ly+uz-5.0+ly)/8.0;
    Nx(2,n,g) =  ux*ly*(ux+ly+uz-5.0+uz)/8.0;

//c   N(n,g) = ux*uy*uz*(ux+uy+uz-5.0)/8.0
    n += 1;
    Nx(0,n,g) =  uy*uz*(ux+uy+uz-5.0+ux)/8.0;
    Nx(1,n,g) =  ux*uz*(ux+uy+uz-5.0+uy)/8.0;
    Nx(2,n,g) =  ux*uy*(ux+uy+uz-5.0+uz)/8.0;

//c   N(n,g) = lx*uy*uz*(lx+uy+uz-5.0)/8.0
    n += 1;
    Nx(0,n,g) = -uy*uz*(lx+uy+uz-5.0+lx)/8.0;
    Nx(1,n,g) =  lx*uz*(lx+uy+uz-5.0+uy)/8.0;
    Nx(2,n,g) =  lx*uy*(lx+uy+uz-5.0+uz)/8.0;

//c   N(n,g) = mx*ly*lz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*ly*lz/4.0;
    Nx(1,n,g) = -mx*lz/4.0;
    Nx(2,n,g) = -mx*ly/4.0;

//c   N(0n,g) = ux*my*lz/4.0
    n += 1;
    Nx(0,n,g) =  my*lz/4.0;
    Nx(1,n,g) =  (ly - uy)*ux*lz/4.0;
    Nx(2,n,g) = -ux*my/4.0;

//c   N(0n,g) = mx*uy*lz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*uy*lz/4.0;
    Nx(1,n,g) =  mx*lz/4.0;
    Nx(2,n,g) = -mx*uy/4.0;

//c   N(0n,g) = lx*my*lz/4.0
    n += 1;
    Nx(0,n,g) = -my*lz/4.0;
    Nx(1,n,g) =  (ly - uy)*lx*lz/4.0;
    Nx(2,n,g) = -lx*my/4.0;

//c   N(0n,g) = mx*ly*uz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*ly*uz/4.0;
    Nx(1,n,g) = -mx*uz/4.0;
    Nx(2,n,g) =  mx*ly/4.0;

//c   N(0n,g) = ux*my*uz/4.0
    n += 1;
    Nx(0,n,g) =  my*uz/4.0;
    Nx(1,n,g) =  (ly - uy)*ux*uz/4.0;
    Nx(2,n,g) =  ux*my/4.0;

//c   N(0n,g) = mx*uy*uz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*uy*uz/4.0;
    Nx(1,n,g) =  mx*uz/4.0;
    Nx(2,n,g) =  mx*uy/4.0;

//c   N(0n,g) = lx*my*uz/4.0
    n += 1;
    Nx(0,n,g) = -my*uz/4.0;
    Nx(1,n,g) =  (ly - uy)*lx*uz/4.0;
    Nx(2,n,g) =  lx*my/4.0;

//c   N(0n,g) = lx*ly*mz/4.0
    n += 1;
    Nx(0,n,g) = -ly*mz/4.0;
    Nx(1,n,g) = -lx*mz/4.0;
    Nx(2,n,g) =  (lz - uz)*lx*ly/4.0;

//c   N(0n,g) = ux*ly*mz/4.0
    n += 1;
    Nx(0,n,g) =  ly*mz/4.0;
    Nx(1,n,g) = -ux*mz/4.0;
    Nx(2,n,g) =  (lz - uz)*ux*ly/4.0;

//c   N(0n,g) = ux*uy*mz/4.0
    n += 1;
    Nx(0,n,g) =  uy*mz/4.0;
    Nx(1,n,g) =  ux*mz/4.0;
    Nx(2,n,g) =  (lz - uz)*ux*uy/4.0;

//c   N(n,g) = lx*uy*mz/4.0
    n += 1;
    Nx(0,n,g) = -uy*mz/4.0;
    Nx(1,n,g) =  lx*mz/4.0;
    Nx(2,n,g) =  (lz - uz)*lx*uy/4.0;
    }
  },

  {ElementType::HEX27, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N,
      Array3<double>& Nx) -> void
    {
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double lz = 1.0 - xi(2,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double uz = 1.0 + xi(2,g);

    double mx = xi(0,g);
    double my = xi(1,g);
    double mz = xi(2,g);

    N(0,g)  = -mx*lx*my*ly*mz*lz/8.0;
    N(1,g)  =  mx*ux*my*ly*mz*lz/8.0;
    N(2,g)  = -mx*ux*my*uy*mz*lz/8.0;
    N(3,g)  =  mx*lx*my*uy*mz*lz/8.0;
    N(4,g)  =  mx*lx*my*ly*mz*uz/8.0;
    N(5,g)  = -mx*ux*my*ly*mz*uz/8.0;
    N(6,g)  =  mx*ux*my*uy*mz*uz/8.0;
    N(7,g)  = -mx*lx*my*uy*mz*uz/8.0;
    N(8,g)  =  lx*ux*my*ly*mz*lz/4.0;
    N(9,g)  = -mx*ux*ly*uy*mz*lz/4.0;
    N(10,g) = -lx*ux*my*uy*mz*lz/4.0;
    N(11,g) =  mx*lx*ly*uy*mz*lz/4.0;
    N(12,g) = -lx*ux*my*ly*mz*uz/4.0;
    N(13,g) =  mx*ux*ly*uy*mz*uz/4.0;
    N(14,g) =  lx*ux*my*uy*mz*uz/4.0;
    N(15,g) = -mx*lx*ly*uy*mz*uz/4.0;
    N(16,g) =  mx*lx*my*ly*lz*uz/4.0;
    N(17,g) = -mx*ux*my*ly*lz*uz/4.0;
    N(18,g) =  mx*ux*my*uy*lz*uz/4.0;
    N(19,g) = -mx*lx*my*uy*lz*uz/4.0;

    N(20,g) = -mx*lx*ly*uy*lz*uz/2.0;
    N(21,g) =  mx*ux*ly*uy*lz*uz/2.0;
    N(22,g) = -lx*ux*my*ly*lz*uz/2.0;
    N(23,g) =  lx*ux*my*uy*lz*uz/2.0;
    N(24,g) = -lx*ux*ly*uy*mz*lz/2.0;
    N(25,g) =  lx*ux*ly*uy*mz*uz/2.0;

    N(26,g) =  lx*ux*ly*uy*lz*uz;

    // N(0)  = -mx*lx*my*ly*mz*lz/8.0
    int n = 0;
    Nx(0,n,g)  = -(lx - mx)*my*ly*mz*lz/8.0;
    Nx(1,n,g)  = -(ly - my)*mx*lx*mz*lz/8.0;
    Nx(2,n,g)  = -(lz - mz)*mx*lx*my*ly/8.0;

    // N(n,g)  =  mx*ux*my*ly*mz*lz/8.0
    n += 1;
    Nx(0,n,g)  =  (mx + ux)*my*ly*mz*lz/8.0;
    Nx(1,n,g)  =  (ly - my)*mx*ux*mz*lz/8.0;
    Nx(2,n,g)  =  (lz - mz)*mx*ux*my*ly/8.0;

    // N(n,g)  = -mx*ux*my*uy*mz*lz/8.0
    n += 1;
    Nx(0,n,g)  = -(mx + ux)*my*uy*mz*lz/8.0;
    Nx(1,n,g)  = -(my + uy)*mx*ux*mz*lz/8.0;
    Nx(2,n,g)  = -(lz - mz)*mx*ux*my*uy/8.0;

    // N(n,g)  =  mx*lx*my*uy*mz*lz/8.0
    n += 1;
    Nx(0,n,g)  =  (lx - mx)*my*uy*mz*lz/8.0;
    Nx(1,n,g)  =  (my + uy)*mx*lx*mz*lz/8.0;
    Nx(2,n,g)  =  (lz - mz)*mx*lx*my*uy/8.0;

    // N(n,g)  =  mx*lx*my*ly*mz*uz/8.0
    n += 1;
    Nx(0,n,g)  =  (lx - mx)*my*ly*mz*uz/8.0;
    Nx(1,n,g)  =  (ly - my)*mx*lx*mz*uz/8.0;
    Nx(2,n,g)  =  (mz + uz)*mx*lx*my*ly/8.0;

    // N(n,g)  = -mx*ux*my*ly*mz*uz/8.0
    n += 1;
    Nx(0,n,g)  = -(mx + ux)*my*ly*mz*uz/8.0;
    Nx(1,n,g)  = -(ly - my)*mx*ux*mz*uz/8.0;
    Nx(2,n,g)  = -(mz + uz)*mx*ux*my*ly/8.0;

    // N(n,g)  =  mx*ux*my*uy*mz*uz/8.0
    n += 1;
    Nx(0,n,g)  =  (mx + ux)*my*uy*mz*uz/8.0;
    Nx(1,n,g)  =  (my + uy)*mx*ux*mz*uz/8.0;
    Nx(2,n,g)  =  (mz + uz)*mx*ux*my*uy/8.0;

    // N(n,g)  = -mx*lx*my*uy*mz*uz/8.0
    n += 1;
    Nx(0,n,g)  = -(lx - mx)*my*uy*mz*uz/8.0;
    Nx(1,n,g)  = -(my + uy)*mx*lx*mz*uz/8.0;
    Nx(2,n,g)  = -(mz + uz)*mx*lx*my*uy/8.0;

    // N(n,g)  =  lx*ux*my*ly*mz*lz/4.0
    n += 1;
    Nx(0,n,g)  =  (lx - ux)*my*ly*mz*lz/4.0;
    Nx(1,n,g)  =  (ly - my)*lx*ux*mz*lz/4.0;
    Nx(2,n,g)  =  (lz - mz)*lx*ux*my*ly/4.0;

    // N(n,g) = -mx*ux*ly*uy*mz*lz/4.0
    n += 1;
    Nx(0,n,g) = -(mx + ux)*ly*uy*mz*lz/4.0;
    Nx(1,n,g) = -(ly - uy)*mx*ux*mz*lz/4.0;
    Nx(2,n,g) = -(lz - mz)*mx*ux*ly*uy/4.0;

    //   N(n,g) = -lx*ux*my*uy*mz*lz/4.0
    n += 1;
    Nx(0,n,g) = -(lx - ux)*my*uy*mz*lz/4.0;
    Nx(1,n,g) = -(my + uy)*lx*ux*mz*lz/4.0;
    Nx(2,n,g) = -(lz - mz)*lx*ux*my*uy/4.0;

    //   N(n,g) =  mx*lx*ly*uy*mz*lz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - mx)*ly*uy*mz*lz/4.0;
    Nx(1,n,g) =  (ly - uy)*mx*lx*mz*lz/4.0;
    Nx(2,n,g) =  (lz - mz)*mx*lx*ly*uy/4.0;

    //   N(n,g) = -lx*ux*my*ly*mz*uz/4.0
    n += 1;
    Nx(0,n,g) = -(lx - ux)*my*ly*mz*uz/4.0;
    Nx(1,n,g) = -(ly - my)*lx*ux*mz*uz/4.0;
    Nx(2,n,g) = -(mz + uz)*lx*ux*my*ly/4.0;

    //   N(n,g) =  mx*ux*ly*uy*mz*uz/4.0
    n += 1;
    Nx(0,n,g) =  (mx + ux)*ly*uy*mz*uz/4.0;
    Nx(1,n,g) =  (ly - uy)*mx*ux*mz*uz/4.0;
    Nx(2,n,g) =  (mz + uz)*mx*ux*ly*uy/4.0;

    //   N(n,g) =  lx*ux*my*uy*mz*uz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*my*uy*mz*uz/4.0;
    Nx(1,n,g) =  (my + uy)*lx*ux*mz*uz/4.0;
    Nx(2,n,g) =  (mz + uz)*lx*ux*my*uy/4.0;

    //   N(n,g) = -mx*lx*ly*uy*mz*uz/4.0
    n += 1;
    Nx(0,n,g) = -(lx - mx)*ly*uy*mz*uz/4.0;
    Nx(1,n,g) = -(ly - uy)*mx*lx*mz*uz/4.0;
    Nx(2,n,g) = -(mz + uz)*mx*lx*ly*uy/4.0;

    //   N(n,g) =  mx*lx*my*ly*lz*uz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - mx)*my*ly*lz*uz/4.0;
    Nx(1,n,g) =  (ly - my)*mx*lx*lz*uz/4.0;
    Nx(2,n,g) =  (lz - uz)*mx*lx*my*ly/4.0;

    //   N(n,g) = -mx*ux*my*ly*lz*uz/4.0
    n += 1;
    Nx(0,n,g) = -(mx + ux)*my*ly*lz*uz/4.0;
    Nx(1,n,g) = -(ly - my)*mx*ux*lz*uz/4.0;
    Nx(2,n,g) = -(lz - uz)*mx*ux*my*ly/4.0;

    //   N(n,g) =  mx*ux*my*uy*lz*uz/4.0
    n += 1;
    Nx(0,n,g) =  (mx + ux)*my*uy*lz*uz/4.0;
    Nx(1,n,g) =  (my + uy)*mx*ux*lz*uz/4.0;
    Nx(2,n,g) =  (lz - uz)*mx*ux*my*uy/4.0;

    //   N(n,g) = -mx*lx*my*uy*lz*uz/4.0
    n += 1;
    Nx(0,n,g) = -(lx - mx)*my*uy*lz*uz/4.0;
    Nx(1,n,g) = -(my + uy)*mx*lx*lz*uz/4.0;
    Nx(2,n,g) = -(lz - uz)*mx*lx*my*uy/4.0;

    //   N(n,g) = -mx*lx*ly*uy*lz*uz/2.0
    n += 1;
    Nx(0,n,g) = -(lx - mx)*ly*uy*lz*uz/2.0;
    Nx(1,n,g) = -(ly - uy)*mx*lx*lz*uz/2.0;
    Nx(2,n,g) = -(lz - uz)*mx*lx*ly*uy/2.0;

    //   N(n,g) =  mx*ux*ly*uy*lz*uz/2.0
    n += 1;
    Nx(0,n,g) =  (mx + ux)*ly*uy*lz*uz/2.0;
    Nx(1,n,g) =  (ly - uy)*mx*ux*lz*uz/2.0;
    Nx(2,n,g) =  (lz - uz)*mx*ux*ly*uy/2.0;

    //   N(n,g) = -lx*ux*my*ly*lz*uz/2.0
    n += 1;
    Nx(0,n,g) = -(lx - ux)*my*ly*lz*uz/2.0;
    Nx(1,n,g) = -(ly - my)*lx*ux*lz*uz/2.0;
    Nx(2,n,g) = -(lz - uz)*lx*ux*my*ly/2.0;

    //   N(n,g) =  lx*ux*my*uy*lz*uz/2.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*my*uy*lz*uz/2.0;
    Nx(1,n,g) =  (my + uy)*lx*ux*lz*uz/2.0;
    Nx(2,n,g) =  (lz - uz)*lx*ux*my*uy/2.0;

    //   N(n,g) = -lx*ux*ly*uy*mz*lz/2.0
    n += 1;
    Nx(0,n,g) = -(lx - ux)*ly*uy*mz*lz/2.0;
    Nx(1,n,g) = -(ly - uy)*lx*ux*mz*lz/2.0;
    Nx(2,n,g) = -(lz - mz)*lx*ux*ly*uy/2.0;

    //   N(n,g) =  lx*ux*ly*uy*mz*uz/2.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*ly*uy*mz*uz/2.0;
    Nx(1,n,g) =  (ly - uy)*lx*ux*mz*uz/2.0;
    Nx(2,n,g) =  (mz + uz)*lx*ux*ly*uy/2.0;

    //   N(n,g) =  lx*ux*ly*uy*lz*uz
    n += 1;
    Nx(0,n,g) =  (lx - ux)*ly*uy*lz*uz;
    Nx(1,n,g) =  (ly - uy)*lx*ux*lz*uz;
    Nx(2,n,g) =  (lz - uz)*lx*ux*ly*uy;
    }
  },

  {ElementType::LIN1, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    N(0,g) = (1.0 - xi(0,g))*0.5;
    N(1,g) = (1.0 + xi(0,g))*0.5;

    Nx(0,0,g) = -0.5;
    Nx(0,1,g) =  0.5;
    }
  },

  {ElementType::LIN2, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    N(0,g) = -xi(0,g)*(1.0 - xi(0,g))*0.50;
    N(1,g) =  xi(0,g)*(1.0 + xi(0,g))*0.50;
    N(2,g) = (1.0 - xi(0,g))*(1.0 + xi(0,g));

    Nx(0,0,g) = -0.50 + xi(0,g);
    Nx(0,1,g) =  0.50 + xi(0,g);
    Nx(0,2,g) = -2.0*xi(0,g);
    }
  },

  {ElementType::QUD4, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);

    N(0,g) = lx*ly / 4.0;
    N(1,g) = ux*ly / 4.0;
    N(2,g) = ux*uy / 4.0;
    N(3,g) = lx*uy / 4.0;

    Nx(0,0,g) = -ly / 4.0;
    Nx(1,0,g) = -lx / 4.0;
    Nx(0,1,g) =  ly / 4.0;
    Nx(1,1,g) = -ux / 4.0;
    Nx(0,2,g) =  uy / 4.0;
    Nx(1,2,g) =  ux / 4.0;
    Nx(0,3,g) = -uy / 4.0;
    Nx(1,3,g) =  lx / 4.0;
    }
  },

  {ElementType::QUD9, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double mx = xi(0,g);
    double my = xi(1,g);

    N(0,g) =  mx*lx*my*ly/4.0;
    N(1,g) = -mx*ux*my*ly/4.0;
    N(2,g) =  mx*ux*my*uy/4.0;
    N(3,g) = -mx*lx*my*uy/4.0;
    N(4,g) = -lx*ux*my*ly*0.50;
    N(5,g) =  mx*ux*ly*uy*0.50;
    N(6,g) =  lx*ux*my*uy*0.50;
    N(7,g) = -mx*lx*ly*uy*0.50;
    N(8,g) =  lx*ux*ly*uy;

    Nx(0,0,g) =  (lx - mx)*my*ly/4.0;
    Nx(1,0,g) =  (ly - my)*mx*lx/4.0;
    Nx(0,1,g) = -(ux + mx)*my*ly/4.0;
    Nx(1,1,g) = -(ly - my)*mx*ux/4.0;
    Nx(0,2,g) =  (ux + mx)*my*uy/4.0;
    Nx(1,2,g) =  (uy + my)*mx*ux/4.0;
    Nx(0,3,g) = -(lx - mx)*my*uy/4.0;
    Nx(1,3,g) = -(uy + my)*mx*lx/4.0;
    Nx(0,4,g) = -(lx - ux)*my*ly*0.50;
    Nx(1,4,g) = -(ly - my)*lx*ux*0.50;
    Nx(0,5,g) =  (ux + mx)*ly*uy*0.50;
    Nx(1,5,g) =  (ly - uy)*mx*ux*0.50;
    Nx(0,6,g) =  (lx - ux)*my*uy*0.50;
    Nx(1,6,g) =  (uy + my)*lx*ux*0.50;
    Nx(0,7,g) = -(lx - mx)*ly*uy*0.50;
    Nx(1,7,g) = -(ly - uy)*mx*lx*0.50;
    Nx(0,8,g) =  (lx - ux)*ly*uy;
    Nx(1,8,g) =  (ly - uy)*lx*ux;
    }
  },

  {ElementType::TET4, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    //std::cout << "[get_element_shape_data] TET4 " << std::endl;

    N(0,g) = xi(0,g);
    N(1,g) = xi(1,g);
    N(2,g) = xi(2,g);
    N(3,g) = 1.0 - xi(0,g) - xi(1,g) - xi(2,g);

    Nx(0,0,g) =  1.0;
    Nx(1,0,g) =  0.0;
    Nx(2,0,g) =  0.0;
    Nx(0,1,g) =  0.0;
    Nx(1,1,g) =  1.0;
    Nx(2,1,g) =  0.0;
    Nx(0,2,g) =  0.0;
    Nx(1,2,g) =  0.0;
    Nx(2,2,g) =  1.0;
    Nx(0,3,g) = -1.0;
    Nx(1,3,g) = -1.0;
    Nx(2,3,g) = -1.0;
    }
  },

  {ElementType::TET10, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    double s = 1.0 - xi(0,g) - xi(1,g) - xi(2,g);
    N(0,g) = xi(0,g)*(2.0*xi(0,g) - 1.0);
    N(1,g) = xi(1,g)*(2.0*xi(1,g) - 1.0);
    N(2,g) = xi(2,g)*(2.0*xi(2,g) - 1.0);
    N(3,g) = s * (2.0*s - 1.0);
    N(4,g) = 4.0*xi(0,g)*xi(1,g);
    N(5,g) = 4.0*xi(1,g)*xi(2,g);
    N(6,g) = 4.0*xi(0,g)*xi(2,g);
    N(7,g) = 4.0*xi(0,g)*s;
    N(8,g) = 4.0*xi(1,g)*s;
    N(9,g) = 4.0*xi(2,g)*s;

    Nx(0,0,g)  =  4.0*xi(0,g) - 1.0;
    Nx(1,0,g)  =  0.0;
    Nx(2,0,g)  =  0.0;

    Nx(0,1,g)  =  0.0;
    Nx(1,1,g)  =  4.0*xi(1,g) - 1.0;
    Nx(2,1,g)  =  0.0;

    Nx(0,2,g)  =  0.0;
    Nx(1,2,g)  =  0.0;
    Nx(2,2,g)  =  4.0*xi(2,g) - 1.0;

    Nx(0,3,g)  =  1.0 - 4.0*s;
    Nx(1,3,g)  =  1.0 - 4.0*s;
    Nx(2,3,g)  =  1.0 - 4.0*s;

    Nx(0,4,g)  =  4.0*xi(1,g);
    Nx(1,4,g)  =  4.0*xi(0,g);
    Nx(2,4,g)  =  0.0;

    Nx(0,5,g)  =  0.0;
    Nx(1,5,g)  =  4.0*xi(2,g);
    Nx(2,5,g)  =  4.0*xi(1,g);

    Nx(0,6,g)  =  4.0*xi(2,g);
    Nx(1,6,g)  =  0.0;
    Nx(2,6,g)  =  4.0*xi(0,g);

    Nx(0,7,g)  =  4.0*( s - xi(0,g));
    Nx(1,7,g)  = -4.0*xi(0,g);
    Nx(2,7,g)  = -4.0*xi(0,g);

    Nx(0,8,g)  = -4.0*xi(1,g);
    Nx(1,8,g)  =  4.0*( s - xi(1,g));
    Nx(2,8,g)  = -4.0*xi(1,g);

    Nx(0,9,g) = -4.0*xi(2,g);
    Nx(1,9,g) = -4.0*xi(2,g);
    Nx(2,9,g) =  4.0*( s - xi(2,g));
    }
  },

  {ElementType::TRI3, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    //std::cout << "[get_element_shape_data] TRI3 " << std::endl;
    N(0,g) = xi(0,g);
    N(1,g) = xi(1,g);
    N(2,g) = 1.0 - xi(0,g) - xi(1,g);

    Nx(0,0,g) =  1.0;
    Nx(1,0,g) =  0.0;
    Nx(0,1,g) =  0.0;
    Nx(1,1,g) =  1.0;
    Nx(0,2,g) = -1.0;
    Nx(1,2,g) = -1.0;
    }
  },

  {ElementType::TRI6, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    double s = 1.0 - xi(0,g) - xi(1,g);
    N(0,g) = xi(0,g) * (2.0*xi(0,g) - 1.0);
    N(1,g) = xi(1,g) * (2.0*xi(1,g) - 1.0);
    N(2,g) = s * (2.0*s - 1.0);
    N(3,g) = 4.0*xi(0,g)*xi(1,g);
    N(4,g) = 4.0*xi(1,g)*s;
    N(5,g) = 4.0*xi(0,g)*s;

    Nx(0,0,g) =  4.0*xi(0,g) - 1.0;
    Nx(1,0,g) =  0.0;

    Nx(0,1,g) =  0.0;
    Nx(1,1,g) =  4.0*xi(1,g) - 1.0;

    Nx(0,2,g) =  1.0 - 4.0*s;
    Nx(1,2,g) =  1.0 - 4.0*s;

    Nx(0,3,g) =  4.0*xi(1,g);
    Nx(1,3,g) =  4.0*xi(0,g);

    Nx(0,4,g) = -4.0*xi(1,g);
    Nx(1,4,g) =  4.0*( s - xi(1,g) );

    Nx(0,5,g) =  4.0*( s - xi(0,g) );
    Nx(1,5,g) = -4.0*xi(0,g);
    }
  },

  {ElementType::WDG, [](const int insd, const int eNoN, const int g, Array<double>& xi, Array<double>& N, 
      Array3<double>& Nx) -> void 
    { 
    double ux = xi(0,g);
    double uy = xi(1,g);
    double uz = 1.0 - ux - uy;
    double s = (1.0 + xi(2,g))*0.5;
    double t = (1.0 - xi(2,g))*0.5;
    N(0,g) = ux*t;
    N(1,g) = uy*t;
    N(2,g) = uz*t;
    N(3,g) = ux*s;
    N(4,g) = uy*s;
    N(5,g) = uz*s;

    Nx(0,0,g) =  t;
    Nx(1,0,g) =  0.0;
    Nx(2,0,g) = -ux*0.50;

    Nx(0,1,g) =  0.0;
    Nx(1,1,g) =  t;
    Nx(2,1,g) = -uy*0.50;

    Nx(0,2,g) = -t;
    Nx(1,2,g) = -t;
    Nx(2,2,g) = -uz*0.50;

    Nx(0,3,g) =  s;
    Nx(1,3,g) =  0.0;
    Nx(2,3,g) =  ux*0.50;

    Nx(0,4,g) =  0.0;
    Nx(1,4,g) =  s;
    Nx(2,4,g) =  uy*0.50;

    Nx(0,5,g) = -s;
    Nx(1,5,g) = -s;
    Nx(2,5,g) =  uz*0.50;
    }
  },



};


//------------------------
// set_element_shape_data 
//------------------------
// Replicates 'SUBROUTINE GETGNN(insd, eType, eNoN, xi, N, Nxi)' defined in NN.f.
//
using SetElementShapeMapType = std::map<ElementType, std::function<void(int, mshType&)>>;

SetElementShapeMapType set_element_shape_data = {

  {ElementType::HEX8, [](int g, mshType& mesh) -> void { 
    auto& xi = mesh.xi;
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double lz = 1.0 - xi(2,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double uz = 1.0 + xi(2,g);

    auto& N = mesh.N;
    N(0,g) = lx*ly*lz/8.0;
    N(1,g) = ux*ly*lz/8.0;
    N(2,g) = ux*uy*lz/8.0;
    N(3,g) = lx*uy*lz/8.0;
    N(4,g) = lx*ly*uz/8.0;
    N(5,g) = ux*ly*uz/8.0;
    N(6,g) = ux*uy*uz/8.0;
    N(7,g) = lx*uy*uz/8.0;

    auto& Nx = mesh.Nx;
    Nx(0,0,g) = -ly*lz/8.0;
    Nx(1,0,g) = -lx*lz/8.0;
    Nx(2,0,g) = -lx*ly/8.0;

    Nx(0,1,g) =  ly*lz/8.0;
    Nx(1,1,g) = -ux*lz/8.0;
    Nx(2,1,g) = -ux*ly/8.0;

    Nx(0,2,g) =  uy*lz/8.0;
    Nx(1,2,g) =  ux*lz/8.0;
    Nx(2,2,g) = -ux*uy/8.0;

    Nx(0,3,g) = -uy*lz/8.0;
    Nx(1,3,g) =  lx*lz/8.0;
    Nx(2,3,g) = -lx*uy/8.0;

    Nx(0,4,g) = -ly*uz/8.0;
    Nx(1,4,g) = -lx*uz/8.0;
    Nx(2,4,g) =  lx*ly/8.0;

    Nx(0,5,g) =  ly*uz/8.0;
    Nx(1,5,g) = -ux*uz/8.0;
    Nx(2,5,g) =  ux*ly/8.0;

    Nx(0,6,g) =  uy*uz/8.0;
    Nx(1,6,g) =  ux*uz/8.0;
    Nx(2,6,g) =  ux*uy/8.0;

    Nx(0,7,g) = -uy*uz/8.0;
    Nx(1,7,g) =  lx*uz/8.0;
    Nx(2,7,g) =  lx*uy/8.0;
    }
  },

  {ElementType::HEX20, [](int g, mshType& mesh) -> void {

    auto& xi = mesh.xi;
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double lz = 1.0 - xi(2,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double uz = 1.0 + xi(2,g);

    double mx = lx*ux;
    double my = ly*uy;
    double mz = lz*uz;

    auto& N = mesh.N;
    N(0, g) = lx*ly*lz*(lx+ly+lz-5.0)/8.0;
    N(1, g) = ux*ly*lz*(ux+ly+lz-5.0)/8.0;
    N(2, g) = ux*uy*lz*(ux+uy+lz-5.0)/8.0;
    N(3, g) = lx*uy*lz*(lx+uy+lz-5.0)/8.0;
    N(4, g) = lx*ly*uz*(lx+ly+uz-5.0)/8.0;
    N(5, g) = ux*ly*uz*(ux+ly+uz-5.0)/8.0;
    N(6, g) = ux*uy*uz*(ux+uy+uz-5.0)/8.0;
    N(7, g) = lx*uy*uz*(lx+uy+uz-5.0)/8.0;
    N(8, g) = mx*ly*lz/4.0;
    N(9, g) = ux*my*lz/4.0;
    N(10, g) = mx*uy*lz/4.0;
    N(11, g) = lx*my*lz/4.0;
    N(12, g) = mx*ly*uz/4.0;
    N(13, g) = ux*my*uz/4.0;
    N(14, g) = mx*uy*uz/4.0;
    N(15, g) = lx*my*uz/4.0;
    N(16, g) = lx*ly*mz/4.0;
    N(17, g) = ux*ly*mz/4.0;
    N(18, g) = ux*uy*mz/4.0;
    N(19, g) = lx*uy*mz/4.0;

    // N(1)  = lx*ly*lz*(lx+ly+lz-5.0)/8.0;
    auto& Nx = mesh.Nx;
    int n = 0;
    Nx(0,n,g) = -ly*lz*(lx+ly+lz-5.0+lx)/8.0;
    Nx(1,n,g) = -lx*lz*(lx+ly+lz-5.0+ly)/8.0;
    Nx(2,n,g) = -lx*ly*(lx+ly+lz-5.0+lz)/8.0;

//c   N(n,g) = ux*ly*lz*(ux+ly+lz-5.0)/8.0;
    n += 1;
    Nx(0,n,g) =  ly*lz*(ux+ly+lz-5.0+ux)/8.0;
    Nx(1,n,g) = -ux*lz*(ux+ly+lz-5.0+ly)/8.0;
    Nx(2,n,g) = -ux*ly*(ux+ly+lz-5.0+lz)/8.0;

//c   N(n,g) = ux*uy*lz*(ux+uy+lz-5.0)/8.0
    n += 1;
    Nx(0,n,g) =  uy*lz*(ux+uy+lz-5.0+ux)/8.0;
    Nx(1,n,g) =  ux*lz*(ux+uy+lz-5.0+uy)/8.0;
    Nx(2,n,g) = -ux*uy*(ux+uy+lz-5.0+lz)/8.0;

//c   N(n,g) = lx*uy*lz*(lx+uy+lz-5.0)/8.0
    n += 1;
    Nx(0,n,g) = -uy*lz*(lx+uy+lz-5.0+lx)/8.0;
    Nx(1,n,g) =  lx*lz*(lx+uy+lz-5.0+uy)/8.0;
    Nx(2,n,g) = -lx*uy*(lx+uy+lz-5.0+lz)/8.0;

//c   N(n,g) = lx*ly*uz*(lx+ly+uz-5.0)/8.0
    n += 1;
    Nx(0,n,g) = -ly*uz*(lx+ly+uz-5.0+lx)/8.0;
    Nx(1,n,g) = -lx*uz*(lx+ly+uz-5.0+ly)/8.0;
    Nx(2,n,g) =  lx*ly*(lx+ly+uz-5.0+uz)/8.0;

//c   N(n,g) = ux*ly*uz*(ux+ly+uz-5.0)/8.0
    n += 1;
    Nx(0,n,g) =  ly*uz*(ux+ly+uz-5.0+ux)/8.0;
    Nx(1,n,g) = -ux*uz*(ux+ly+uz-5.0+ly)/8.0;
    Nx(2,n,g) =  ux*ly*(ux+ly+uz-5.0+uz)/8.0;

//c   N(n,g) = ux*uy*uz*(ux+uy+uz-5.0)/8.0
    n += 1;
    Nx(0,n,g) =  uy*uz*(ux+uy+uz-5.0+ux)/8.0;
    Nx(1,n,g) =  ux*uz*(ux+uy+uz-5.0+uy)/8.0;
    Nx(2,n,g) =  ux*uy*(ux+uy+uz-5.0+uz)/8.0;

//c   N(n,g) = lx*uy*uz*(lx+uy+uz-5.0)/8.0
    n += 1;
    Nx(0,n,g) = -uy*uz*(lx+uy+uz-5.0+lx)/8.0;
    Nx(1,n,g) =  lx*uz*(lx+uy+uz-5.0+uy)/8.0;
    Nx(2,n,g) =  lx*uy*(lx+uy+uz-5.0+uz)/8.0;

//c   N(n,g) = mx*ly*lz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*ly*lz/4.0;
    Nx(1,n,g) = -mx*lz/4.0;
    Nx(2,n,g) = -mx*ly/4.0;

//c   N(0n,g) = ux*my*lz/4.0
    n += 1;
    Nx(0,n,g) =  my*lz/4.0;
    Nx(1,n,g) =  (ly - uy)*ux*lz/4.0;
    Nx(2,n,g) = -ux*my/4.0;

//c   N(0n,g) = mx*uy*lz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*uy*lz/4.0;
    Nx(1,n,g) =  mx*lz/4.0;
    Nx(2,n,g) = -mx*uy/4.0;

//c   N(0n,g) = lx*my*lz/4.0
    n += 1;
    Nx(0,n,g) = -my*lz/4.0;
    Nx(1,n,g) =  (ly - uy)*lx*lz/4.0;
    Nx(2,n,g) = -lx*my/4.0;

//c   N(0n,g) = mx*ly*uz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*ly*uz/4.0;
    Nx(1,n,g) = -mx*uz/4.0;
    Nx(2,n,g) =  mx*ly/4.0;

//c   N(0n,g) = ux*my*uz/4.0
    n += 1;
    Nx(0,n,g) =  my*uz/4.0;
    Nx(1,n,g) =  (ly - uy)*ux*uz/4.0;
    Nx(2,n,g) =  ux*my/4.0;

//c   N(0n,g) = mx*uy*uz/4.0
    n += 1;
    Nx(0,n,g) =  (lx - ux)*uy*uz/4.0;
    Nx(1,n,g) =  mx*uz/4.0;
    Nx(2,n,g) =  mx*uy/4.0;

//c   N(0n,g) = lx*my*uz/4.0
    n += 1;
    Nx(0,n,g) = -my*uz/4.0;
    Nx(1,n,g) =  (ly - uy)*lx*uz/4.0;
    Nx(2,n,g) =  lx*my/4.0;

//c   N(0n,g) = lx*ly*mz/4.0
    n += 1;
    Nx(0,n,g) = -ly*mz/4.0;
    Nx(1,n,g) = -lx*mz/4.0;
    Nx(2,n,g) =  (lz - uz)*lx*ly/4.0;

//c   N(0n,g) = ux*ly*mz/4.0
    n += 1;
    Nx(0,n,g) =  ly*mz/4.0;
    Nx(1,n,g) = -ux*mz/4.0;
    Nx(2,n,g) =  (lz - uz)*ux*ly/4.0;

//c   N(0n,g) = ux*uy*mz/4.0
    n += 1;
    Nx(0,n,g) =  uy*mz/4.0;
    Nx(1,n,g) =  ux*mz/4.0;
    Nx(2,n,g) =  (lz - uz)*ux*uy/4.0;

//c   N(n,g) = lx*uy*mz/4.0
    n += 1;
    Nx(0,n,g) = -uy*mz/4.0;
    Nx(1,n,g) =  lx*mz/4.0;
    Nx(2,n,g) =  (lz - uz)*lx*uy/4.0;
    }
  },

  {ElementType::HEX27, [](int g, mshType& mesh) -> void {

    auto& xi = mesh.xi;
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double lz = 1.0 - xi(2,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double uz = 1.0 + xi(2,g);

    double mx = xi(0,g);
    double my = xi(1,g);
    double mz = xi(2,g);

    auto& N = mesh.N;
    N(0,g)  = -mx*lx*my*ly*mz*lz/8.0;
    N(1,g)  =  mx*ux*my*ly*mz*lz/8.0;
    N(2,g)  = -mx*ux*my*uy*mz*lz/8.0;
    N(3,g)  =  mx*lx*my*uy*mz*lz/8.0;
    N(4,g)  =  mx*lx*my*ly*mz*uz/8.0;
    N(5,g)  = -mx*ux*my*ly*mz*uz/8.0;
    N(6,g)  =  mx*ux*my*uy*mz*uz/8.0;
    N(7,g)  = -mx*lx*my*uy*mz*uz/8.0;
    N(8,g)  =  lx*ux*my*ly*mz*lz/4.0;
    N(9,g) = -mx*ux*ly*uy*mz*lz/4.0;
    N(10,g) = -lx*ux*my*uy*mz*lz/4.0;
    N(11,g) =  mx*lx*ly*uy*mz*lz/4.0;
    N(12,g) = -lx*ux*my*ly*mz*uz/4.0;
    N(13,g) =  mx*ux*ly*uy*mz*uz/4.0;
    N(14,g) =  lx*ux*my*uy*mz*uz/4.0;
    N(15,g) = -mx*lx*ly*uy*mz*uz/4.0;
    N(16,g) =  mx*lx*my*ly*lz*uz/4.0;
    N(17,g) = -mx*ux*my*ly*lz*uz/4.0;
    N(18,g) =  mx*ux*my*uy*lz*uz/4.0;
    N(19,g) = -mx*lx*my*uy*lz*uz/4.0;

    N(20,g) = -mx*lx*ly*uy*lz*uz/2.0;
    N(21,g) =  mx*ux*ly*uy*lz*uz/2.0;
    N(22,g) = -lx*ux*my*ly*lz*uz/2.0;
    N(23,g) =  lx*ux*my*uy*lz*uz/2.0;
    N(24,g) = -lx*ux*ly*uy*mz*lz/2.0;
    N(25,g) =  lx*ux*ly*uy*mz*uz/2.0;

    N(26,g) =  lx*ux*ly*uy*lz*uz;

    auto& Nxi = mesh.Nx;
    int n = 0;
    Nxi(0,n,g)  = -(lx - mx)*my*ly*mz*lz/8.0;
    Nxi(1,n,g)  = -(ly - my)*mx*lx*mz*lz/8.0;
    Nxi(2,n,g)  = -(lz - mz)*mx*lx*my*ly/8.0;

    n += 1;
    Nxi(0,n,g)  =  (mx + ux)*my*ly*mz*lz/8.0;
    Nxi(1,n,g)  =  (ly - my)*mx*ux*mz*lz/8.0;
    Nxi(2,n,g)  =  (lz - mz)*mx*ux*my*ly/8.0;

    n += 1;
    Nxi(0,n,g)  = -(mx + ux)*my*uy*mz*lz/8.0;
    Nxi(1,n,g)  = -(my + uy)*mx*ux*mz*lz/8.0;
    Nxi(2,n,g)  = -(lz - mz)*mx*ux*my*uy/8.0;

    n += 1;
    Nxi(0,n,g)  =  (lx - mx)*my*uy*mz*lz/8.0;
    Nxi(1,n,g)  =  (my + uy)*mx*lx*mz*lz/8.0;
    Nxi(2,n,g)  =  (lz - mz)*mx*lx*my*uy/8.0;

    n += 1;
    Nxi(0,n,g)  =  (lx - mx)*my*ly*mz*uz/8.0;
    Nxi(1,n,g)  =  (ly - my)*mx*lx*mz*uz/8.0;
    Nxi(2,n,g)  =  (mz + uz)*mx*lx*my*ly/8.0;

    n += 1;
    Nxi(0,n,g)  = -(mx + ux)*my*ly*mz*uz/8.0;
    Nxi(1,n,g)  = -(ly - my)*mx*ux*mz*uz/8.0;
    Nxi(2,n,g)  = -(mz + uz)*mx*ux*my*ly/8.0;

    n += 1;
    Nxi(0,n,g)  =  (mx + ux)*my*uy*mz*uz/8.0;
    Nxi(1,n,g)  =  (my + uy)*mx*ux*mz*uz/8.0;
    Nxi(2,n,g)  =  (mz + uz)*mx*ux*my*uy/8.0;

    n += 1;
    Nxi(0,n,g)  = -(lx - mx)*my*uy*mz*uz/8.0;
    Nxi(1,n,g)  = -(my + uy)*mx*lx*mz*uz/8.0;
    Nxi(2,n,g)  = -(mz + uz)*mx*lx*my*uy/8.0;

    n += 1;
    Nxi(0,n,g)  =  (lx - ux)*my*ly*mz*lz/4.0;
    Nxi(1,n,g)  =  (ly - my)*lx*ux*mz*lz/4.0;
    Nxi(2,n,g)  =  (lz - mz)*lx*ux*my*ly/4.0;

    n += 1;
    Nxi(0,n,g) = -(mx + ux)*ly*uy*mz*lz/4.0;
    Nxi(1,n,g) = -(ly - uy)*mx*ux*mz*lz/4.0;
    Nxi(2,n,g) = -(lz - mz)*mx*ux*ly*uy/4.0;

    n += 1;
    Nxi(0,n,g) = -(lx - ux)*my*uy*mz*lz/4.0;
    Nxi(1,n,g) = -(my + uy)*lx*ux*mz*lz/4.0;
    Nxi(2,n,g) = -(lz - mz)*lx*ux*my*uy/4.0;

    n += 1;
    Nxi(0,n,g) =  (lx - mx)*ly*uy*mz*lz/4.0;
    Nxi(1,n,g) =  (ly - uy)*mx*lx*mz*lz/4.0;
    Nxi(2,n,g) =  (lz - mz)*mx*lx*ly*uy/4.0;

    n += 1;
    Nxi(0,n,g) = -(lx - ux)*my*ly*mz*uz/4.0;
    Nxi(1,n,g) = -(ly - my)*lx*ux*mz*uz/4.0;
    Nxi(2,n,g) = -(mz + uz)*lx*ux*my*ly/4.0;

    n += 1;
    Nxi(0,n,g) =  (mx + ux)*ly*uy*mz*uz/4.0;
    Nxi(1,n,g) =  (ly - uy)*mx*ux*mz*uz/4.0;
    Nxi(2,n,g) =  (mz + uz)*mx*ux*ly*uy/4.0;

    n += 1;
    Nxi(0,n,g) =  (lx - ux)*my*uy*mz*uz/4.0;
    Nxi(1,n,g) =  (my + uy)*lx*ux*mz*uz/4.0;
    Nxi(2,n,g) =  (mz + uz)*lx*ux*my*uy/4.0;

    n += 1;
    Nxi(0,n,g) = -(lx - mx)*ly*uy*mz*uz/4.0;
    Nxi(1,n,g) = -(ly - uy)*mx*lx*mz*uz/4.0;
    Nxi(2,n,g) = -(mz + uz)*mx*lx*ly*uy/4.0;

    n += 1;
    Nxi(0,n,g) =  (lx - mx)*my*ly*lz*uz/4.0;
    Nxi(1,n,g) =  (ly - my)*mx*lx*lz*uz/4.0;
    Nxi(2,n,g) =  (lz - uz)*mx*lx*my*ly/4.0;

    n += 1;
    Nxi(0,n,g) = -(mx + ux)*my*ly*lz*uz/4.0;
    Nxi(1,n,g) = -(ly - my)*mx*ux*lz*uz/4.0;
    Nxi(2,n,g) = -(lz - uz)*mx*ux*my*ly/4.0;

    n += 1;
    Nxi(0,n,g) =  (mx + ux)*my*uy*lz*uz/4.0;
    Nxi(1,n,g) =  (my + uy)*mx*ux*lz*uz/4.0;
    Nxi(2,n,g) =  (lz - uz)*mx*ux*my*uy/4.0;

    n += 1;
    Nxi(0,n,g) = -(lx - mx)*my*uy*lz*uz/4.0;
    Nxi(1,n,g) = -(my + uy)*mx*lx*lz*uz/4.0;
    Nxi(2,n,g) = -(lz - uz)*mx*lx*my*uy/4.0;

    n += 1;
    Nxi(0,n,g) = -(lx - mx)*ly*uy*lz*uz/2.0;
    Nxi(1,n,g) = -(ly - uy)*mx*lx*lz*uz/2.0;
    Nxi(2,n,g) = -(lz - uz)*mx*lx*ly*uy/2.0;

    n += 1;
    Nxi(0,n,g) =  (mx + ux)*ly*uy*lz*uz/2.0;
    Nxi(1,n,g) =  (ly - uy)*mx*ux*lz*uz/2.0;
    Nxi(2,n,g) =  (lz - uz)*mx*ux*ly*uy/2.0;

    n += 1;
    Nxi(0,n,g) = -(lx - ux)*my*ly*lz*uz/2.0;
    Nxi(1,n,g) = -(ly - my)*lx*ux*lz*uz/2.0;
    Nxi(2,n,g) = -(lz - uz)*lx*ux*my*ly/2.0;

    n += 1;
    Nxi(0,n,g) =  (lx - ux)*my*uy*lz*uz/2.0;
    Nxi(1,n,g) =  (my + uy)*lx*ux*lz*uz/2.0;
    Nxi(2,n,g) =  (lz - uz)*lx*ux*my*uy/2.0;

    n += 1;
    Nxi(0,n,g) = -(lx - ux)*ly*uy*mz*lz/2.0;
    Nxi(1,n,g) = -(ly - uy)*lx*ux*mz*lz/2.0;
    Nxi(2,n,g) = -(lz - mz)*lx*ux*ly*uy/2.0;

    n += 1;
    Nxi(0,n,g) =  (lx - ux)*ly*uy*mz*uz/2.0;
    Nxi(1,n,g) =  (ly - uy)*lx*ux*mz*uz/2.0;
    Nxi(2,n,g) =  (mz + uz)*lx*ux*ly*uy/2.0;

    n += 1;
    Nxi(0,n,g) =  (lx - ux)*ly*uy*lz*uz;
    Nxi(1,n,g) =  (ly - uy)*lx*ux*lz*uz;
    Nxi(2,n,g) =  (lz - uz)*lx*ux*ly*uy;
    }
  },

  {ElementType::LIN1, [](int g, mshType& mesh) -> void { 
    //std::cout << "[set_element_shape_data] **************************" << std::endl;
    //std::cout << "[set_element_shape_data] ERROR: LIN1 not supported." << std::endl;
    //std::cout << "[set_element_shape_data] **************************" << std::endl;
    auto& xi = mesh.xi;
    auto& N = mesh.N;
    N(0,g) = (1.0 - xi(0,g))*0.5;
    N(1,g) = (1.0 + xi(0,g))*0.5;

    auto& Nx = mesh.Nx;
    Nx(0,0,g) = -0.5;
    Nx(0,1,g) =  0.5;
    }
  },

  {ElementType::LIN2, [](int g, mshType& mesh) -> void {
    auto& xi = mesh.xi;
    auto& N = mesh.N;
    N(0,g) = -xi(0,g)*(1.0 - xi(0,g))*0.50;
    N(1,g) =  xi(0,g)*(1.0 + xi(0,g))*0.50;
    N(2,g) = (1.0 - xi(0,g))*(1.0 + xi(0,g));

    auto& Nx = mesh.Nx;
    Nx(0,0,g) = -0.50 + xi(0,g);
    Nx(0,1,g) =  0.50 + xi(0,g);
    Nx(0,2,g) = -2.0*xi(0,g);
    }
  },

  {ElementType::QUD4, [](int g, mshType& mesh) -> void {
    auto& xi = mesh.xi;
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);

    auto& N = mesh.N;
    N(0,g) = lx*ly / 4.0;    
    N(1,g) = ux*ly / 4.0;    
    N(2,g) = ux*uy / 4.0;    
    N(3,g) = lx*uy / 4.0;    

    auto& Nx = mesh.Nx;
    Nx(0,0,g) = -ly / 4.0;
    Nx(1,0,g) = -lx / 4.0;
    Nx(0,1,g) =  ly / 4.0;
    Nx(1,1,g) = -ux / 4.0;
    Nx(0,2,g) =  uy / 4.0;
    Nx(1,2,g) =  ux / 4.0;
    Nx(0,3,g) = -uy / 4.0;
    Nx(1,3,g) =  lx / 4.0;
    }
  },

 {ElementType::QUD9, [](int g, mshType& mesh) -> void {
    auto& xi = mesh.xi;
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double mx = xi(0,g);
    double my = xi(1,g);

    auto& N = mesh.N;
    N(0,g) =  mx*lx*my*ly/4.0;
    N(1,g) = -mx*ux*my*ly/4.0;
    N(2,g) =  mx*ux*my*uy/4.0;
    N(3,g) = -mx*lx*my*uy/4.0;
    N(4,g) = -lx*ux*my*ly*0.50;
    N(5,g) =  mx*ux*ly*uy*0.50;
    N(6,g) =  lx*ux*my*uy*0.50;
    N(7,g) = -mx*lx*ly*uy*0.50;
    N(8,g) =  lx*ux*ly*uy;

    auto& Nx = mesh.Nx;
    Nx(0,0,g) =  (lx - mx)*my*ly/4.0;
    Nx(1,0,g) =  (ly - my)*mx*lx/4.0;

    Nx(0,1,g) = -(ux + mx)*my*ly/4.0;
    Nx(1,1,g) = -(ly - my)*mx*ux/4.0;

    Nx(0,2,g) =  (ux + mx)*my*uy/4.0;
    Nx(1,2,g) =  (uy + my)*mx*ux/4.0;

    Nx(0,3,g) = -(lx - mx)*my*uy/4.0;
    Nx(1,3,g) = -(uy + my)*mx*lx/4.0;

    Nx(0,4,g) = -(lx - ux)*my*ly*0.50;
    Nx(1,4,g) = -(ly - my)*lx*ux*0.50;

    Nx(0,5,g) =  (ux + mx)*ly*uy*0.50;
    Nx(1,5,g) =  (ly - uy)*mx*ux*0.50;

    Nx(0,6,g) =  (lx - ux)*my*uy*0.50;
    Nx(1,6,g) =  (uy + my)*lx*ux*0.50;

    Nx(0,7,g) = -(lx - mx)*ly*uy*0.50;
    Nx(1,7,g) = -(ly - uy)*mx*lx*0.50;

    Nx(0,8,g) =  (lx - ux)*ly*uy;
    Nx(1,8,g) =  (ly - uy)*lx*ux;
    }
  },

  {ElementType::TET4, [](int g, mshType& mesh) -> void { 
    auto& xi = mesh.xi;
    auto& N = mesh.N;
    N(0,g) = xi(0,g);
    N(1,g) = xi(1,g);
    N(2,g) = xi(2,g);
    N(3,g) = 1.0 - xi(0,g) - xi(1,g) - xi(2,g);

    auto& Nx = mesh.Nx;
    Nx(0,0,g) =  1.0;
    Nx(1,0,g) =  0.0;
    Nx(2,0,g) =  0.0;
    Nx(0,1,g) =  0.0;
    Nx(1,1,g) =  1.0;
    Nx(2,1,g) =  0.0;
    Nx(0,2,g) =  0.0;
    Nx(1,2,g) =  0.0;
    Nx(2,2,g) =  1.0;
    Nx(0,3,g) = -1.0;
    Nx(1,3,g) = -1.0;
    Nx(2,3,g) = -1.0;
    }
  },

  {ElementType::TET10, [](int g, mshType& mesh) -> void {
    auto& xi = mesh.xi;
    auto& N = mesh.N;
    double s = 1.0 - xi(0,g) - xi(1,g) - xi(2,g);
    N(0,g)  = xi(0,g)*(2.0*xi(0,g) - 1.0);
    N(1,g)  = xi(1,g)*(2.0*xi(1,g) - 1.0);
    N(2,g)  = xi(2,g)*(2.0*xi(2,g) - 1.0);
    N(3,g)  = s    *(2.0*s     - 1.0);
    N(4,g)  = 4.0*xi(0,g)*xi(1,g);
    N(5,g)  = 4.0*xi(1,g)*xi(2,g);
    N(6,g)  = 4.0*xi(0,g)*xi(2,g);
    N(7,g)  = 4.0*xi(0,g)*s;
    N(8,g)  = 4.0*xi(1,g)*s;
    N(9,g) = 4.0*xi(2,g)*s;

    auto& Nx = mesh.Nx;
    Nx(0,0,g)  =  4.0*xi(0,g) - 1.0;
    Nx(1,0,g)  =  0.0;
    Nx(2,0,g)  =  0.0;

    Nx(0,1,g)  =  0.0;
    Nx(1,1,g)  =  4.0*xi(1,g) - 1.0;
    Nx(2,1,g)  =  0.0;

    Nx(0,2,g)  =  0.0;
    Nx(1,2,g)  =  0.0;
    Nx(2,2,g)  =  4.0*xi(2,g) - 1.0;

    Nx(0,3,g)  =  1.0 - 4.0*s;
    Nx(1,3,g)  =  1.0 - 4.0*s;
    Nx(2,3,g)  =  1.0 - 4.0*s;

    Nx(0,4,g)  =  4.0*xi(1,g);
    Nx(1,4,g)  =  4.0*xi(0,g);
    Nx(2,4,g)  =  0.0;

    Nx(0,5,g)  =  0.0;
    Nx(1,5,g)  =  4.0*xi(2,g);
    Nx(2,5,g)  =  4.0*xi(1,g);

    Nx(0,6,g)  =  4.0*xi(2,g);
    Nx(1,6,g)  =  0.0;
    Nx(2,6,g)  =  4.0*xi(0,g);

    Nx(0,7,g)  =  4.0*( s - xi(0,g));
    Nx(1,7,g)  = -4.0*xi(0,g);
    Nx(2,7,g)  = -4.0*xi(0,g);

    Nx(0,8,g)  = -4.0*xi(1,g);
    Nx(1,8,g)  =  4.0*( s - xi(1,g));
    Nx(2,8,g)  = -4.0*xi(1,g);

    Nx(0,9,g) = -4.0*xi(2,g);
    Nx(1,9,g) = -4.0*xi(2,g);
    Nx(2,9,g) =  4.0*( s - xi(2,g));
    }
  },

  {ElementType::TRI3, [](int g, mshType& mesh) -> void { 
    auto& xi = mesh.xi;
    auto& N = mesh.N;
    N(0,g) = xi(0,g);
    N(1,g) = xi(1,g);
    N(2,g) = 1.0 - xi(0,g) - xi(1,g);

    auto& Nxi = mesh.Nx;
    Nxi(0,0,g) =  1.0;
    Nxi(1,0,g) =  0.0;
    Nxi(0,1,g) =  0.0;
    Nxi(1,1,g) =  1.0;
    Nxi(0,2,g) = -1.0;
    Nxi(1,2,g) = -1.0;
    }
  },

  {ElementType::TRI6, [](int g, mshType& mesh) -> void {
    auto& xi = mesh.xi;
    auto& N = mesh.N;

    double s = 1.0 - xi(0,g) - xi(1,g);
    N(0,g) = xi(0,g)*( 2.0*xi(0,g) - 1.0 );
    N(1,g) = xi(1,g)*( 2.0*xi(1,g) - 1.0 );
    N(2,g) = s    *( 2.0*s     - 1.0 );
    N(3,g) = 4.0*xi(0,g)*xi(1,g);
    N(4,g) = 4.0*xi(1,g)*s;
    N(5,g) = 4.0*xi(0,g)*s;

    auto& Nxi = mesh.Nx;
    Nxi(0,0,g) =  4.0*xi(0,g) - 1.0;
    Nxi(1,0,g) =  0.0;
    Nxi(0,1,g) =  0.0;
    Nxi(1,1,g) =  4.0*xi(1,g) - 1.0;
    Nxi(0,2,g) =  1.0 - 4.0*s;
    Nxi(1,2,g) =  1.0 - 4.0*s;
    Nxi(0,3,g) =  4.0*xi(1,g);
    Nxi(1,3,g) =  4.0*xi(0,g);
    Nxi(0,4,g) = -4.0*xi(1,g);
    Nxi(1,4,g) =  4.0*( s - xi(1,g) );
    Nxi(0,5,g) =  4.0*( s - xi(0,g) );
    Nxi(1,5,g) = -4.0*xi(0,g);
    }
  },

  {ElementType::WDG, [](int g, mshType& mesh) -> void 
    { 
    auto& xi = mesh.xi;
    auto& N = mesh.N;
    double ux = xi(0,g);
    double uy = xi(1,g);
    double uz = 1.0 - ux - uy;
    double s = (1.0 + xi(2,g))*0.5;
    double t = (1.0 - xi(2,g))*0.5;
    N(0,g) = ux*t;
    N(1,g) = uy*t;
    N(2,g) = uz*t;
    N(3,g) = ux*s;
    N(4,g) = uy*s;
    N(5,g) = uz*s;

    auto& Nxi = mesh.Nx;
    Nxi(0,0,g) =  t;
    Nxi(1,0,g) =  0.0;
    Nxi(2,0,g) = -ux*0.50;

    Nxi(0,1,g) =  0.0;
    Nxi(1,1,g) =  t;
    Nxi(2,1,g) = -uy*0.50;

    Nxi(0,2,g) = -t;
    Nxi(1,2,g) = -t;
    Nxi(2,2,g) = -uz*0.50;

    Nxi(0,3,g) =  s;
    Nxi(1,3,g) =  0.0;
    Nxi(2,3,g) =  ux*0.50;

    Nxi(0,4,g) =  0.0;
    Nxi(1,4,g) =  s;
    Nxi(2,4,g) =  uy*0.50;

    Nxi(0,5,g) = -s;
    Nxi(1,5,g) = -s;
    Nxi(2,5,g) =  uz*0.50;
    }
  },

};

//---------------------
// set_face_shape_data
//---------------------
// Define a map type used to face element shape function data.
//
// This reproduces 'SUBROUTINE GETGNN(insd, eType, eNoN, xi, N, Nxi)' in NN.f.
//
using SetFaceShapeMapType = std::map<ElementType, std::function<void(int, faceType&)>>;

SetFaceShapeMapType set_face_shape_data = {

  {ElementType::PNT, [](int g, faceType& face) -> void 
    {
    face.N(0,g) = 1.0;
    }
  },

  {ElementType::QUD8, [](int g, faceType& face) -> void 
    {
    auto& xi = face.xi;
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double mx = lx*ux;
    double my = ly*uy;

    auto& N = face.N;
    N(0,g) = lx*ly*(lx+ly-3.0)/4.0;
    N(1,g) = ux*ly*(ux+ly-3.0)/4.0;
    N(2,g) = ux*uy*(ux+uy-3.0)/4.0;
    N(3,g) = lx*uy*(lx+uy-3.0)/4.0;
    N(4,g) = mx*ly*0.50;
    N(5,g) = ux*my*0.50;
    N(6,g) = mx*uy*0.50;
    N(7,g) = lx*my*0.50;

    auto& Nxi = face.Nx;
    Nxi(0,0,g) = -ly*(lx+ly-3.0+lx)/4.0;
    Nxi(1,0,g) = -lx*(lx+ly-3.0+ly)/4.0;

    Nxi(0,1,g) =  ly*(ux+ly-3.0+ux)/4.0;
    Nxi(1,1,g) = -ux*(ux+ly-3.0+ly)/4.0;

    Nxi(0,2,g) =  uy*(ux+uy-3.0+ux)/4.0;
    Nxi(1,2,g) =  ux*(ux+uy-3.0+uy)/4.0;

    Nxi(0,3,g) = -uy*(lx+uy-3.0+lx)/4.0;
    Nxi(1,3,g) =  lx*(lx+uy-3.0+uy)/4.0;

    Nxi(0,4,g) =  (lx - ux)*ly*0.50;
    Nxi(1,4,g) = -mx*0.50;

    Nxi(0,5,g) =  my*0.50;
    Nxi(1,5,g) =  (ly - uy)*ux*0.50;

    Nxi(0,6,g) =  (lx - ux)*uy*0.50;
    Nxi(1,6,g) =  mx*0.50;

    Nxi(0,7,g) = -my*0.50;
    Nxi(1,7,g) =  (ly - uy)*lx*0.50;
    }
  },

  {ElementType::QUD9, [](int g, faceType& face) -> void 
    {
    auto& xi = face.xi;
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);
    double mx = xi(0,g);
    double my = xi(1,g);

    auto& N = face.N;
    N(0,g) =  mx*lx*my*ly/4.0;
    N(1,g) = -mx*ux*my*ly/4.0;
    N(2,g) =  mx*ux*my*uy/4.0;
    N(3,g) = -mx*lx*my*uy/4.0;
    N(4,g) = -lx*ux*my*ly*0.50;
    N(5,g) =  mx*ux*ly*uy*0.50;
    N(6,g) =  lx*ux*my*uy*0.50;
    N(7,g) = -mx*lx*ly*uy*0.50;
    N(8,g) =  lx*ux*ly*uy;

    auto& Nx = face.Nx;
    Nx(0,0,g) =  (lx - mx)*my*ly/4.0;
    Nx(1,0,g) =  (ly - my)*mx*lx/4.0;
    Nx(0,1,g) = -(ux + mx)*my*ly/4.0;
    Nx(1,1,g) = -(ly - my)*mx*ux/4.0;
    Nx(0,2,g) =  (ux + mx)*my*uy/4.0;
    Nx(1,2,g) =  (uy + my)*mx*ux/4.0;
    Nx(0,3,g) = -(lx - mx)*my*uy/4.0;
    Nx(1,3,g) = -(uy + my)*mx*lx/4.0;
    Nx(0,4,g) = -(lx - ux)*my*ly*0.50;
    Nx(1,4,g) = -(ly - my)*lx*ux*0.50;
    Nx(0,5,g) =  (ux + mx)*ly*uy*0.50;
    Nx(1,5,g) =  (ly - uy)*mx*ux*0.50;
    Nx(0,6,g) =  (lx - ux)*my*uy*0.50;
    Nx(1,6,g) =  (uy + my)*lx*ux*0.50;
    Nx(0,7,g) = -(lx - mx)*ly*uy*0.50;
    Nx(1,7,g) = -(ly - uy)*mx*lx*0.50;
    Nx(0,8,g) =  (lx - ux)*ly*uy;
    Nx(1,8,g) =  (ly - uy)*lx*ux;
    }
  },

  {ElementType::LIN1, [](int g, faceType& face) -> void 
    {
    face.N(0,g) = 0.5 * (1.0 - face.xi(0,g));
    face.N(1,g) = 0.5 * (1.0 + face.xi(0,g));

    face.Nx(0,0,g) = -0.5;
    face.Nx(0,1,g) =  0.5;
    }
  },

  {ElementType::LIN2, [](int g, faceType& face) -> void
    {
    auto& xi = face.xi;
    auto& N = face.N;
    N(0,g) = -xi(0,g)*(1.0 - xi(0,g))*0.50;
    N(1,g) =  xi(0,g)*(1.0 + xi(0,g))*0.50;
    N(2,g) = (1.0 - xi(0,g))*(1.0 + xi(0,g));

    auto& Nx = face.Nx;
    Nx(0,0,g) = -0.50 + xi(0,g);
    Nx(0,1,g) =  0.50 + xi(0,g);
    Nx(0,2,g) = -2.0*xi(0,g);
    }
  },

  {ElementType::QUD4, [](int g, faceType& face) -> void {
    auto& xi = face.xi;
    double lx = 1.0 - xi(0,g);
    double ly = 1.0 - xi(1,g);
    double ux = 1.0 + xi(0,g);
    double uy = 1.0 + xi(1,g);

    auto& N =face.N;
    N(0,g) = lx*ly / 4.0;
    N(1,g) = ux*ly / 4.0;
    N(2,g) = ux*uy / 4.0;
    N(3,g) = lx*uy / 4.0;

    auto& Nx = face.Nx;
    Nx(0,0,g) = -ly / 4.0;
    Nx(1,0,g) = -lx / 4.0;
    Nx(0,1,g) =  ly / 4.0;
    Nx(1,1,g) = -ux / 4.0;
    Nx(0,2,g) =  uy / 4.0;
    Nx(1,2,g) =  ux / 4.0;
    Nx(0,3,g) = -uy / 4.0;
    Nx(1,3,g) =  lx / 4.0;
    }
  },

  {ElementType::TRI3, [](int g, faceType& face) -> void 
    {
    face.N(0,g) = face.xi(0,g);
    face.N(1,g) = face.xi(1,g);
    face.N(2,g) = 1.0 - face.xi(0,g) - face.xi(1,g);

    face.Nx(0,0,g) = 1.0;
    face.Nx(1,0,g) = 0.0;

    face.Nx(0,1,g) = 0.0;
    face.Nx(1,1,g) = 1.0;

    face.Nx(0,2,g) = -1.0;
    face.Nx(1,2,g) = -1.0;
    }
  },

  {ElementType::TRI6, [](int g, faceType& face) -> void
    {
    auto& xi = face.xi;
    auto& N = face.N;

    double s = 1.0 - xi(0,g) - xi(1,g);
    N(0,g) = xi(0,g)*( 2.0*xi(0,g) - 1.0 );
    N(1,g) = xi(1,g)*( 2.0*xi(1,g) - 1.0 );
    N(2,g) = s    *( 2.0*s     - 1.0 );
    N(3,g) = 4.0*xi(0,g)*xi(1,g);
    N(4,g) = 4.0*xi(1,g)*s;
    N(5,g) = 4.0*xi(0,g)*s;

    auto& Nxi = face.Nx;
    Nxi(0,0,g) =  4.0*xi(0,g) - 1.0;
    Nxi(1,0,g) =  0.0;

    Nxi(0,1,g) =  0.0;
    Nxi(1,1,g) =  4.0*xi(1,g) - 1.0;

    Nxi(0,2,g) =  1.0 - 4.0*s;
    Nxi(1,2,g) =  1.0 - 4.0*s;

    Nxi(0,3,g) =  4.0*xi(1,g);
    Nxi(1,3,g) =  4.0*xi(0,g);

    Nxi(0,4,g) = -4.0*xi(1,g);
    Nxi(1,4,g) =  4.0*( s - xi(1,g) );

    Nxi(0,5,g) =  4.0*( s - xi(0,g) );
    Nxi(1,5,g) = -4.0*xi(0,g);
    }
  },


};
