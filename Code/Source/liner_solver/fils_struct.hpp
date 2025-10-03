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

#ifndef FSI_LINEAR_SOLVER_H 
#define FSI_LINEAR_SOLVER_H 

#include "CmMod.h"

#include "mpi.h"

#include <map>

/// SELECTED_REAL_KIND(P,R) returns the kind value of a real data type with 
///
///  1) decimal precision of at least P digits, 
///
///  2) exponent range of at least R
///
/// \code {.f}
/// INTEGER, PARAMETER :: LSRP4  = SELECTED_REAL_KIND(6, 37)
/// INTEGER, PARAMETER :: LSRP8  = SELECTED_REAL_KIND(15, 307)      // C++ double
/// INTEGER, PARAMETER :: LSRP16 = SELECTED_REAL_KIND(33, 4931)     // C++ long double
/// INTEGER, PARAMETER :: LSRP = LSRP8
/// \endcode
///
/// The enums here replicate the PARAMETERs defined
/// in FSILS_STRUCT.h.
//
namespace fsi_linear_solver {

enum class BcType
{
  BC_TYPE_Dir = 0,
  BC_TYPE_Neu = 1
};

enum class BcopType
{
  BCOP_TYPE_ADD = 0,
  BCOP_TYPE_PRE = 1
};

enum LinearSolverType
{
  LS_TYPE_CG = 798,
  LS_TYPE_GMRES = 797, 
  LS_TYPE_NS = 796, 
  LS_TYPE_BICGS = 795
};

class FSILS_commuType 
{
  public:
    /// Free of created          (USE)
    int foC;         

    /// If this the master       (USE)    
    int masF;            
    
    /// Master ID                (USE)
    int master;          

    /// ID of this proc.         (USE)
    int task;            

    /// Task in FORTRAN indexing (USE)
    int tF;              

    /// Total number of tasks    (USE)
    int nTasks;          

    /// MPI communicator         (IN)
    MPI_Comm comm;       
};

class FSILS_cSType
{
  public:
    /// The processor to communicate with
    int iP = -1;

    /// Number of data to be commu
    int n = 0;

    /// Pointer to the data for commu
    Vector<int> ptr;
};

class FSILS_faceType
{
  public:
    /// Free or created                (USE)
    bool foC = false;

    /// Neu: P/Q coupling              (USE)
    bool coupledFlag = false;

    /// Neu: shared between proces     (USE)
    bool sharedFlag = false;

    /// Included in the computations   (IN)
    bool incFlag = false;

    /// Number of nodes                (IN)
    int nNo = 0;

    /// Degrees of freedom for val     (IN)
    int dof = 0;

    /// Dir/Neu                        (IN)
    BcType bGrp = BcType::BC_TYPE_Dir;

    /// Only for data alignment
    int reserved = 0;

    /// Global node number             (IN)
    Vector<int> glob;

    /// ||Sai||**2._LSRP                   (USE)
    double nS = 0.0;

    /// Neu: P = res*Q                 (IN)
    double res = 0.0;

    /// nodal Sai for Neu              (IN)
    Array<double> val;

    /// Neu W*Sai                      (TMP)
    Array<double> valM;

    /// Flag for virtual face          (USE)
    bool vrtual = false;

    /// Index of mesh in msh associated with this face
    int iM = -1;

    /// Index of face in msh[iM].fa that caps this face
    int iFa = -1;

    /// Index of face in lhs.face that caps this face
    int faInCap = -1;

    /// Flag indicating if this face is capped by another face
    bool isCapped = false;
};

/// @brief Modified in:
///
///  fsils_lhs_create()
//
class FSILS_lhsType 
{
  public:
    /// Free of created                     (USE)
    bool foC = false; 

    bool debug_active = false; 

    /// Global number of nodes              (IN)
    int gnNo = 0;

    /// Number of nodes                     (IN)
    int nNo = 0;

    /// Number of non-zero in lhs           (IN)
    int nnz = 0;

    ///  Number of faces                     (IN)
    int nFaces = 0;

    /// nNo of this proc                    (USE)
    int mynNo = 0;

    /// nNo of shared with lower proc       (USE)
    int shnNo = 0;

    /// Number of communication requests    (USE)
    int nReq = 0;

    /// Column pointer                      (USE)
    Vector<int> colPtr;

    /// Row pointer                         (USE)
    Array<int> rowPtr;

    /// Diagonal pointer                    (USE)
    Vector<int> diagPtr;

    /// Mapping of nodes                    (USE)
    Vector<int> map;

    FSILS_commuType commu;

    std::vector<FSILS_cSType> cS;

    std::vector<FSILS_faceType> face;
};

class FSILS_subLsType 
{
  public:
    /// Successful solving          (OUT)
    bool suc;       

    //int suc;       // Successful solving          (OUT)

    /// Maximum iteration           (IN)
    int mItr;      

    /// Space dimension             (IN)
    int sD;        

    /// Number of iteration         (OUT)
    int itr;       

    /// Number of Ax multiply       (OUT)
    int cM;        

    /// Number of |x| norms         (OUT)
    int cN;        

    /// Number of <x.y> dot products(OUT)
    int cD;        

    /// Only for data alignment     (-)
    int reserve;   

    /// Absolute tolerance          (IN)
    double absTol; 

    /// Relative tolerance          (IN)
    double relTol; 

    /// Initial norm of residual    (OUT)
    double iNorm;  

    /// Final norm of residual      (OUT)
    double fNorm;  

    /// Res. rduction in last itr.  (OUT)
    double dB;     

    /// Calling duration            (OUT)
    double callD;  
};

class FSILS_lsType 
{
  public:
    /// Free of created             (USE)
    int foC;                     

    /// Which one of LS             (IN)
    LinearSolverType LS_type;    

    /// Contribution of mom. res.   (OUT)
    int Resm;                    

    /// Contribution of cont. res.  (OUT)
    int Resc;                    
    
    FSILS_subLsType GM;
    FSILS_subLsType CG;
    FSILS_subLsType RI;
};


};

#endif

