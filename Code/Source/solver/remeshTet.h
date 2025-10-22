// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

void remesh3d_tetgen(const int nPoints, const int nFacets, const double* pointList,   
    const int* facetList, const std::array<double,3>& params, int* pOK);

