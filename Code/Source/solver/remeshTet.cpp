// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "tetgen.h"
#include <array>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

/// @brief Interface to Tetgen for remeshing purposes.
class tetOptions {
  public:
    double maxRadRatio;
    double minDihedAng;
    double maxEdgeSize;
    int optimLevel;
    int optimScheme;
    tetOptions();
    double maxTetVol(double r) {return (pow(r,3)/(6.0*sqrt(2)));}
  };

tetOptions::tetOptions () {
  maxRadRatio = 1.15;
  minDihedAng = 10.0;
  maxEdgeSize = 0.15;
  optimLevel  = 2;
  optimScheme = 7;
}

void remesh3d_tetgen(const int nPoints, const int nFacets, const double* pointList, 
                     const int* facetList, const std::array<double,3>& params, int* pOK)
{
   //std::cout << "========== remesh3d_tetgen ==========" << std::endl;
   tetgenio in, out;
   tetgenio::facet *f;
   tetgenio::polygon *p;
   char fname [250];
   char switches [250];
   tetOptions options;

   *pOK = 0;
   in.firstnumber = 0;
   //in.firstnumber = 1;
   in.numberofpoints = nPoints;
   in.pointlist = new REAL [in.numberofpoints * 3];

   for (int i=0; i < in.numberofpoints; i++)
   {
      for (int j=0; j < 3; j++) {
         in.pointlist[3*i+j] = *pointList;
         ++pointList;
      }
   }

   in.numberoffacets = nFacets;
   in.facetlist = new tetgenio::facet[in.numberoffacets];
   in.facetmarkerlist = new int[in.numberoffacets];

   for (int i=0; i < in.numberoffacets; i++)
   {
      f = &in.facetlist[i];
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
      f->numberofholes = 0;
      f->holelist = NULL;
      p = &f->polygonlist[0];
      p->numberofvertices = 3;
      p->vertexlist = new int [p->numberofvertices];
      for (int j=0; j < 3; j++)
      {
         p->vertexlist[j] = *facetList;
         ++facetList;
      }
      in.facetmarkerlist[i] = 0;
   }

   options.maxRadRatio = params[0];
   options.minDihedAng = params[1];
   options.maxEdgeSize = params[2];

   //std::cout << "[remesh3d_tetgen] Using parameter <maxRadRatio> " << options.maxRadRatio << "\n";
   //std::cout << "[remesh3d_tetgen] Using parameter <minDihedAng> " << options.minDihedAng << "\n";
   //std::cout << "[remesh3d_tetgen] Using parameter <maxEdgeSize> " << options.maxEdgeSize << "\n\n";

   int len;
   len = sprintf(switches,"pYq%.2f/%.1fa%8.3eO%d/%d",    \
         options.maxRadRatio, options.minDihedAng,   \
         options.maxTetVol(options.maxEdgeSize),      \
         options.optimLevel,   options.optimScheme);

   if ( len > 250 )
   {
      std::cout << "    ERROR: Length of switch exceeded limit (250 char)\n";
      *pOK = -1;
      return;
   }
   tetrahedralize(switches, &in, &out);

   strcpy(fname, "new-vol-mesh-cpp");
   out.save_nodes(fname);
   out.save_elements(fname);

   return;
}
