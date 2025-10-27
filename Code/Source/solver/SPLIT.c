// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

//--------------------------------------------------------------------
//
// Interface to Metis for partitioning the mesh.
//
//--------------------------------------------------------------------

#ifndef SEQ

#include"parmetislib.h"

int split_(int *nElptr, int *eNoNptr, int *eNoNbptr, int *IEN,
   int *nPartsPtr, idx_t *iElmdist, float *iWgt, idx_t *part)
{

   int i, e, a, nEl=*nElptr, eNoN=*eNoNptr, eNoNb=*eNoNbptr,
      nparts, nTasks=*nPartsPtr, wgtflag, numflag, ncon, task,
      ncommonnodes, options[10], *exRanks, nExRanks, *map, edgecut;

   float ubvec[MAXNCON], *wgt;
   idx_t *eptr, *eind, *elmdist;

   map     = (int *)malloc(nTasks*sizeof(int));
   exRanks = (int *)malloc(nTasks*sizeof(int));
   wgt     = (float *)malloc(nTasks*sizeof(float));
   elmdist = (idx_t *)malloc((nTasks+1)*sizeof(idx_t));
   MPI_Group newGrp, tmpGrp;
   MPI_Comm comm;

   MPI_Comm_rank(MPI_COMM_WORLD, &task);

// This is for the case one of the processors doesn't posses any
// part of this mesh
   nExRanks   = 0;
   nparts     = 0;
   elmdist[0] = 0;
   for (i=0; i<nTasks ; i++ ) {
      if ((iElmdist[i+1] - iElmdist[i]) == 0) {
         exRanks[nExRanks] = i;
         nExRanks++;
      } else {
         map[nparts] = i;
         wgt[nparts] = iWgt[i];
         elmdist[nparts+1] = iElmdist[i+1];
         nparts++;
      }
   }
   if (nExRanks == 0) {
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);
   } else {
      MPI_Comm_group(MPI_COMM_WORLD, &tmpGrp);
      MPI_Group_excl(tmpGrp, nExRanks, exRanks, &newGrp);
      MPI_Comm_create(MPI_COMM_WORLD, newGrp, &comm);
      MPI_Group_free(&tmpGrp);
      MPI_Group_free(&newGrp);
      if (nEl == 0) return 0;
   }
// If there is just one processor left, we give all the element to that
// one
   if ( nparts == 1 ) {
      for (e=0; e<nEl; e++) part[e] = task;
      return -1;
   }

   eptr = (idx_t *)malloc((nEl+1)*sizeof(idx_t));
   eind = (idx_t *)malloc(nEl*eNoN*sizeof(idx_t));

   for (e=0; e<=nEl; e++) {
      eptr[e] = e*eNoN;
   }
   for (a=0; a<nEl*eNoN; a++) {
      eind[a] = IEN[a] - 1;
   }
   wgtflag = 0;
   numflag = 0;
   ncon = 1;
   ncommonnodes = eNoNb;

   for (i=0; i<ncon; i++) ubvec[i] = UNBALANCE_FRACTION;

   options[0] = 1;
   options[PMV3_OPTION_DBGLVL] = 0;
   options[PMV3_OPTION_SEED] = 10;

   ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL, &wgtflag,
      &numflag, &ncon, &ncommonnodes, &nparts, wgt, ubvec,
      options, &edgecut, part, &comm);

   MPI_Comm_free(&comm);

// Mapping proc ID to the global numbering
   if (edgecut == 0) {
// In this case ParMETIS has failed. So I assume each
      for (e=0; e<nEl; e++) {
         part[e] = task;
      }
   } else {
      for (e=0; e<nEl; e++) {
         i = part[e];
         i = map[i];
         part[e] = i;
      }
   }

   // free arrays
   free (eind);
   free (eptr);
   free (wgt);
   free (elmdist);
   free (map);
   free (exRanks);

   return edgecut;
}
#else
int split_(int *nElptr, int *eNoNptr, int *eNoNbptr, int *IEN,
   int *nPartsptr, int *iElmdist, float *iWgt, int *part)  {
   return 0;
}
#endif

//###################################################################

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

double cput_()
{
   struct timespec time;
   double res;

   #ifdef __MACH__    // OSX does not have clock_gettime
      clock_serv_t cclock;
      mach_timespec_t mts;
      host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
      clock_get_time(cclock, &mts);
      mach_port_deallocate(mach_task_self(), cclock);
      time.tv_sec = mts.tv_sec;
      time.tv_nsec = mts.tv_nsec;
   #else
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time);
   #endif

   res = ((double)time.tv_nsec)/1000000000.0 + (double)(time.tv_sec%1000000);

   return res;
}

