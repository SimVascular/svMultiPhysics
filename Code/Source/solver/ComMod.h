// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The classes defined here duplicate the data structures in the Fortran COMMOD module
// defined in MOD.f. 
//
// All of the data structures used for the mesh, boundarsy conditions and solver parameters, etc. 
// are defined here.

#ifndef COMMOD_H 
#define COMMOD_H 

#include "Array.h"
#include "Array3.h"
#include "CepMod.h"
#include "ChnlMod.h"
#include "CmMod.h"
#include "Parameters.h"
#include "RobinBoundaryCondition.h"
#include "Timer.h"
#include "Vector.h"

#include "DebugMsg.h"

#include "consts.h"

#include "fils_struct.hpp"

#include "Parameters.h"

#include "ArtificialNeuralNetMaterial.h"

#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

class LinearAlgebra;

/// @brief Fourier coefficients that are used to specify unsteady BCs
//
class fcType
{
  public:

    bool defined() { return n != 0; };

    // If this is a ramp function
    bool lrmp = false;

    // Number of Fourier coefficient
    int n = 0;

    // No. of dimensions (scalar or vector)
    int d = 0;
   
    // Initial value
    Vector<double> qi;

    // Time derivative of linear part
    Vector<double> qs;

    // Period
    double T = 0.0;

    // Initial time
    double ti = 0.0;

    // Imaginary part of coefficint
    Array<double> i;

    // Real part of coefficint
    Array<double> r;
};

/// @brief Moving boundary data structure (used for general BC)
//
class MBType
{
  public:

    bool defined() { return dof != 0; };

    // Degrees of freedom of d(:,.,.)
    int dof = 0;

    // Number of time points to be read
    int nTP = 0;

    // The period of data
    double period = 0.0;

    // Time points
    Vector<double> t;

    // Displacements at each direction, location, and time point
    Array3<double> d;
};

class rcrType
{
  public:

    // Proximal resistance
    double Rp = 0.0;

    // Capacitance
    double C = 0.0;

    // Distance resistance
    double Rd = 0.0;

    // Distal pressure
    double Pd = 0.0;

    // Initial value
    double Xo = 0.0;
};

/// @brief Boundary condition data type
//
class bcType
{
  public:
    // Strong/Weak application of Dirichlet BC
    bool weakDir = false;

    // Whether load vector changes with deformation
    // (Neu - struct/ustruct only)
    bool flwP = false;

    // Strong/Weak application of Dirichlet BC
    int clsFlgRis = 0;

    // Pre/Res/Flat/Para... boundary types
    //
    // This stores differnt BCs as bitwise values. 
    //
    int bType = 0;

    // Pointer to coupledBC%face
    int cplBCptr = -1;

    // The face index that corresponds to this BC
    int iFa = -1;

    // The mesh index that corresponds to this BC
    int iM = -1;

    // Pointer to FSILS%bc
    int lsPtr = -1;

    // Undeforming Neu BC master-slave node parameters.
    int masN = 0;

    // Defined steady value
    double g = 0.0;

    // Neu: defined resistance
    double r = 0.0;

    // Robin: VTP file path for per-node stiffness and damping
    std::string robin_vtp_file = "";

    // RIS0D: resistance
    double resistance = 0.0;

    // Penalty parameters for weakly applied Dir BC
    Vector<double> tauB{0.0, 0.0};
    //double tauB[2];

    // Direction vector for imposing the BC
    Vector<int> eDrn;

    // Defined steady vector (traction)
    Vector<double> h;

    // Spatial dependant BC (profile data)
    Vector<double> gx;

    // General BC (unsteady and UD combination)
    //
    // This is declare ALLOCATABLE in MOD.f. 
    //
    MBType gm;

    // Time dependant BC (Unsteady imposed value);
    //
    // This is declare ALLOCATABLE in MOD.f. 
    //
    fcType gt;

    // Neu: RCR
    rcrType RCR;

    // Robin BC class
    RobinBoundaryCondition robin_bc;
};

/// @brief Class storing data for B-Splines.
//
class bsType 
{
  public:

    // Number of knots (p + nNo + 1)
    int n = 0;

    // Number of Gauss points for integration
    int nG = 0;

    // Number of knot spans (element)
    int nEl = 0;

    // Number of control points (nodes)
    int nNo = 0;

    // Number of sample points in each element (for output)
    int nSl = 0;

    // The order
    int p = 0;

    // Knot vector.
    Vector<double> xi;
};

/// @brief Function spaces (basis) type.
//
class fsType {

  public:

    fsType();

    void destroy();

    // Whether the basis function is linear
    bool lShpF = false;

    // Element type
    consts::ElementType eType = consts::ElementType::NA;

    // Number of basis functions, typically equals msh%eNoN
    int eNoN = 0;

    // Number of Gauss points for integration
    int nG = 0;

    // Gauss weights
    Vector<double> w;

    // Gauss integration points in parametric space
    Array<double> xi;

    // Bounds on Gauss integration points in parametric space
    Array<double> xib;

    // Parent shape function
    Array<double> N;

    // Bounds on shape functions
    Array<double> Nb;

    // Parent shape functions gradient
    Array3<double> Nx;

    // Second derivatives of shape functions - used for shells & IGA
    Array3<double> Nxx;
};

//--------
// bfType
//--------
// Body force data structure type
//
class bfType
{
  public:

    std::string file_name;
    std::string mesh_name;

    // Type of body force applied
    int bType = 0;

    // No. of dimensions (1 or nsd)
    int dof = 0;

    // Mesh index corresponding to this body force
    int iM = -1;

    // Steady value
    Vector<double> b;

    // Steady but spatially dependant
    Array<double> bx;

    // Time dependant (unsteady imposed value)
    fcType bt;

    // General (unsteady and spatially dependent combination)
    MBType bm;
};

// Fiber stress type
class fibStrsType
{
  public:

    // Type of fiber stress
    int fType = 0;

    // Constant steady value
    double g = 0.0;

    // Directional stress distribution parameters
    // Fraction of active stress in fiber direction (default: 1.0)
    double eta_f = 1.0;
    
    // Fraction of active stress in sheet direction (default: 0.0)
    double eta_s = 0.0;
    
    // Fraction of active stress in sheet-normal direction (default: 0.0)
    double eta_n = 0.0;

    // Unsteady time-dependent values
    fcType gt;
};

/// @brief Structural domain type
//
class stModelType
{
  public:
    // Type of constitutive model (volumetric) for struct/FSI
    consts::ConstitutiveModelType volType = consts::ConstitutiveModelType::stIso_NA;

    // Penalty parameter
    double Kpen = 0.0;

    // Type of constitutive model (isochoric) for struct/FSI
    consts::ConstitutiveModelType isoType = consts::ConstitutiveModelType::stIso_NA;

    // Parameters specific to the constitutive model (isochoric)
    // NeoHookean model (C10 = mu/2)
    double C10 = 0.0;

    // Mooney-Rivlin model (C10, C01)
    double C01 = 0.0;

    // Holzapfel model(a, b, aff, bff, ass, bss, afs, bfs, kap)
    double a = 0.0;
    double b = 0.0;
    double aff = 0.0;
    double bff = 0.0;
    double ass = 0.0;
    double bss = 0.0;
    double afs = 0.0;
    double bfs = 0.0;

    // Collagen fiber dispersion parameter (Holzapfel model)
    double kap = 0.0;

    // Heaviside function parameter (Holzapfel-Ogden model)
    double khs = 100.0;

    // Lee-Sacks model
    double a0 = 0.0;
    double b1 = 0.0;
    double b2 = 0.0;
    double mu0 = 0.0;

    // Fiber reinforcement stress
    fibStrsType Tf;

    // CANN Model/UAnisoHyper_inv
    ArtificialNeuralNetMaterial paramTable;
  
    stModelType();
};


/// @brief Fluid viscosity model type
//
class fluidViscModelType
{
  public:

    // Type of constitutive model for fluid viscosity
    consts::FluidViscosityModelType viscType = consts::FluidViscosityModelType::viscType_NA;

    // Limiting zero shear-rate viscosity value
    double mu_o = 0.0;

    // Limiting high shear-rate viscosity (asymptotic) value (at infinity)
    double mu_i = 0.0;

    // Strain-rate tensor multiplier
    double lam = 0.0;

    // Strain-rate tensor exponent
    double a = 0.0;

    // Power-law exponent
    double n = 0.0;
};

/// @brief Fluid viscosity model type
//
class solidViscModelType
{
  public:

    // Type of constitutive model for fluid viscosity
    consts::SolidViscosityModelType viscType = consts::SolidViscosityModelType::viscType_NA;

    // Viscosity value
    double mu = 0.0;
};

/// @brief Domain type is to keep track with element belong to which domain
/// and also different physical quantities
//
class dmnType
{
  public:
    dmnType();
    ~dmnType();

    // The domain ID. Default includes entire domain
    int Id = -1;

    // Which physics must be solved in this domain
    consts::EquationType phys = consts::EquationType::phys_NA;

    // The volume of this domain
    double v = 0.0;

    // General physical properties such as density, elastic modulus...
    // FIX davep double prop[maxNProp] ;
    std::map<consts::PhysicalProperyType,double> prop;
    //double prop[consts::maxNProp];

    // Electrophysiology model
    cepModelType cep;

    // Structure material model
    stModelType stM;

    // Viscosity model for fluids
    fluidViscModelType fluid_visc;

    // Viscosity model for solids
    solidViscModelType solid_visc;
};

/// @brief Mesh adjacency (neighboring element for each element)
//
class adjType
{
  public:
    void destroy();

    // No of non-zeros
    int nnz = 0;

    // Column pointer
    Vector<int> pcol;

    // Row pointer
    Vector<int> prow;

};

/// @brief Tracer type used for immersed boundaries. Identifies traces of
/// nodes and integration points on background mesh elements
//
class traceType
{
  public:
    // No. of non-zero nodal traces
    int n = 0;

    // No. of non-zero integration point traces
    int nG = 0;

    // Local to global nodes maping nNo --> tnNo
    Vector<int> gN;

    // Self pointer of each trace to the IB integration point and
    // element ID
    Array<int> gE;

    // Nodal trace pointer array stores two values for each trace.
    // (1) background mesh element to which the trace points to,
    // (2) mesh ID
    Array<int> nptr;

    // Integration point tracer array stores two values for each trace
    // (1) background mesh element to which the trace points to,
    // (2) mesh ID
    Array<int> gptr;

    // Parametric coordinate for each nodal trace
    Array<double> xi;

    // Parametric coordinate for each Gauss point trace
    Array<double> xiG;
};

/// @brief The face type containing mesh at boundary
//
class faceType
{
  public:
    faceType();
    ~faceType();

    void destroy();

    //faceType& operator=(const faceType& rhs);

    // Parametric direction normal to this face (NURBS)
    int d = 0;

    // Number of nodes (control points) in a single element
    int eNoN = 0;

    // Element type
    consts::ElementType eType = consts::ElementType::NA;

    // The mesh index that this face belongs to
    int iM = 0;

    // Number of elements
    int nEl = 0;

    // Global number of elements
    int gnEl = 0;

    // Number of function spaces
    int nFs = 0;

    // Number of Gauss points for integration
    int nG = 0;

    // Number of nodes
    int nNo = 0;

    // Global element Ids
    Vector<int> gE;

    // Global node Ids
    Vector<int> gN;

    // Global to local maping tnNo --> nNo
    Vector<int> lN;

    // Connectivity array
    Array<int> IEN;

    // EBC array (gE + gIEN)
    Array<int> gebc;

    // Surface area
    double area = 0.0;

    // Gauss point weights
    Vector<double> w;

    // Position coordinates
    Array<double> x;

    // Gauss points in parametric space
    Array<double> xi;

    // Shape functions at Gauss points
    Array<double> N;

    // Normal vector to each nodal point
    Array<double> nV;

    // Shape functions derivative at Gauss points
    // double Nx(:,:,:);
    Array3<double> Nx;

    // Second derivatives of shape functions - for shells & IGA
    // double Nxx(:,:,:);
    Array3<double> Nxx;

    // Face name for flux files
    std::string name;

    // Face nodal adjacency
    adjType nAdj;

    // Face element adjacency
    adjType eAdj;

    // Function spaces (basis)
    std::vector<fsType> fs;

    // TRI3 quadrature modifier
    double qmTRI3 = 2.0/3.0;
};

/// @brief Store options for output types.
//
struct OutputOptions {
  bool boundary_integral = false;
  bool spatial = false;
  bool volume_integral = false;

  bool no_options_set() {
    return !(boundary_integral | spatial | volume_integral);
  }

  void set_option(consts::OutputType type, bool value) 
  {
    if (type == consts::OutputType::boundary_integral) {
      boundary_integral = value;
    } else if (type == consts::OutputType::spatial) {
      spatial = value;
    } else if (type == consts::OutputType::volume_integral) {
      volume_integral = value;
    }
  }
};
    
/// @brief Declared type for outputed variables
//
class outputType
{
  public:

    // Options to write various output types.
    OutputOptions options;

    // The group that this belong to (one of outType_*)
    consts::OutputNameType grp = consts::OutputNameType::outGrp_NA;

    // Length of the outputed variable
    int l = 0;

    // Offset from the first index
    int o = 0;

    // The name to be used for the output and also in input file
    std::string name;
};

/// @brief Linear system of equations solver type
//
class lsType
{
  public:

    /// @brief LS solver                     (IN)
    consts::SolverType LS_type = consts::SolverType::lSolver_NA;

    /// @brief Successful solving            (OUT)
    bool suc = false;

    /// @brief Maximum iterations            (IN)
    int mItr = 1000;

    /// @brief Space dimension               (IN)
    int sD = 0;

    /// @brief Number of iteration           (OUT)
    int itr = 0;

    /// @brief Number of Ax multiple         (OUT)
    int cM = 0;

    /// @brief Number of |x| norms           (OUT)
    int cN = 0;

    /// @brief Number of <x.y> dot products  (OUT)
    int cD = 0;

    /// @brief Only for data alignment       (-)
    int reserve = 0;

    /// @brief Absolute tolerance            (IN)
    double absTol = 1e-08;

    /// @brief Relative tolerance            (IN)
    double relTol = 0.0;

    /// @brief Initial norm of residual      (OUT)
    double iNorm = 0.0;

    /// @brief Final norm of residual        (OUT)
    double fNorm = 0.0;

    /// @brief Res. rduction in last itr.    (OUT)
    double dB = 0.0;

    /// @brief Calling duration              (OUT)
    double callD = 0.0;

    //@brief Configuration file for linear solvers (Trilinos, PETSc)
    std::string config;
};


/// @brief Contact model type
//
class cntctModelType
{
  public:
    // Contact model
    consts::ContactModelType cType = consts::ContactModelType::cntctM_NA;

    // Penalty parameter
    double k = 0.0;

    // Min depth of penetration
    double h = 0.0;

    // Max depth of penetration
    double c = 0.0;

    // Min norm of face normals in contact
    double al = 0.0;

    // Tolerance
    double tol = 0.0; 
};

class cplFaceType
{
  public:
    // GenBC_Dir/GenBC_Neu
    consts::CplBCType bGrp = consts::CplBCType::cplBC_NA;

    // Pointer to X
    int Xptr = -1;

    // Internal genBC use
    int eqv = 0;

    // Flow rates at t
    double Qo = 0.0;

    // Flow rates at t+dt
    double Qn = 0.0;

    // Pressures at t
    double Po = 0.0;

    // Pressures at t+dt
    double Pn = 0.0;

    // Imposed flow/pressure
    double y = 0.0;

    // Name of the face
    std::string name;

    // RCR type BC
    rcrType RCR;
};

//----------------------------
// svZeroDSolverInterfaceType
//----------------------------
// This class stores information used to interface to
// the svZeroDSolver.
//
class svZeroDSolverInterfaceType
{
  public:

    // The path/name of the 0D solver shared library.
    std::string solver_library;

    // Maps a 0D block name with a 3D face name representing the 
    // coupling of a 0D block with a 3D face.
    std::map<std::string,std::string> block_surface_map;

    // The path/name of the 0D solver JSON file.
    std::string configuration_file;

    // How often the 0D contribution to solver tangent matrix is added.
    std::string coupling_type;

    // Initialize the 0D flows.
    bool have_initial_flows = false;
    double initial_flows;

    // Initialize the 0D pressures.
    bool have_initial_pressures = false;
    double initial_pressures;

    // If the data has been set for the interface. This is only true if 
    // the svZeroDSolver_interface XML parameter has been defined.
    bool has_data = false;

    void set_data(const svZeroDSolverInterfaceParameters& params);
    void add_block_face(const std::string& block_name, const std::string& face_name);
};

/// @brief For coupled 0D-3D problems
//
class cplBCType
{
  public:
    cplBCType();
    /// @brief Is multi-domain active
    bool coupled = false;

    /// @brief Whether to use genBC
    bool useGenBC = false;

    //  Whether to use svZeroD
    bool useSvZeroD = false;

    //  Whether to initialize RCR from flow data
    bool initRCR = false;

    /// @brief Number of coupled faces
    int nFa = 0;

    /// @brief Number of unknowns in the 0D domain
    int nX = 0;

    /// @brief Number of output variables addition to nX
    int nXp = 0;

    /// @brief Implicit/Explicit/Semi-implicit schemes
    consts::CplBCType schm = consts::CplBCType::cplBC_NA;
    //int schm = cplBC_NA;

    /// @brief Path to the 0D code binary file
    std::string binPath;

    /// @brief File name for communication between 0D and 3D
    std::string commuName;
    //std::string commuName = ".CPLBC_0D_3D.tmp";

    svZeroDSolverInterfaceType svzerod_solver_interface;

    /// @brief The name of history file containing "X"
    std::string saveName;
    //std::string(LEN=stdL) :: saveName = "LPN.dat";

    /// @brief New time step unknowns in the 0D domain
    Vector<double> xn;

    /// @brief Old time step unknowns in the 0D domain
    Vector<double> xo;

    /// @brief Output variables to be printed
    Vector<double> xp;

    /// @brief Data structure used for communicating with 0D code
    std::vector<cplFaceType> fa;
};

/// @brief This is the container for a mesh or NURBS patch, those specific
/// to NURBS are noted
//
class mshType
{
  public:
    mshType();
    std::string dname = "";

/*
    mshType(const mshType &other) 
    { 
      std::cout << "c c c c c mshType copy c c c c c" << std::endl;
    } 
      
    mshType& operator = (const mshType &other)
    {
      std::cout << "= = = = = mshType assignment = = = = =" << std::endl;
      return *this;
    } 
*/

    ~mshType() 
    { 
      //std::cout << "- - - - -  mshType dtor - - - - -   dname: " << dname << std::endl;
    };

    /// @brief Whether the shape function is linear
    bool lShpF = false;

    /// @brief Whether the mesh is shell
    bool lShl = false;

    /// @brief Whether the mesh is fibers (Purkinje)
    bool lFib = false;

    /// @brief Element type
    consts::ElementType eType = consts::ElementType::NA;
    //int eType = eType_NA

    /// @brief Number of nodes (control points) in a single element
    int eNoN = 0;

    /// @brief Global number of elements (knot spans)
    int gnEl = 0;

    /// @brief Global number of nodes (control points) on a single mesh
    int gnNo = 0;

    /// @brief Number of element face. Used for reading Gambit mesh files
    int nEf = 0;

    /// @brief Number of elements (knot spans)
    int nEl = 0;

    /// @brief Number of faces
    int nFa = 0;

    /// @brief Number of function spaces
    int nFs = 0;

    /// @brief Number of Gauss points for integration
    int nG = 0;

    /// @brief Number of nodes (control points) for 2D elements?
    int nNo = 0;

    /// @brief Number of elements sample points to be outputs (NURBS)
    int nSl = 0;

    /// @brief The element type recognized by VTK format
    int vtkType = 0;

    /// @brief Number of fiber directions
    int nFn = 0;

    /// @brief Mesh scale factor
    double scF = 0.0;

    /// @brief IB: Mesh size parameter
    double dx = 0.0;

    /// @brief RIS resistance value
    double res = 0.0;

    /// @brief RIS projection tolerance 
    double tol = 0.0;

    /// @brief The volume of this mesh
    double v = 0.0;

    /// @breif ordering: node ordering for boundaries
    std::vector<std::vector<int>> ordering;

    /// @brief Element distribution between processors
    Vector<int> eDist;

    /// @brief Element domain ID number
    Vector<int> eId;

    /// @brief Global nodes maping nNo --> tnNo
    Vector<int> gN;

    /// @brief GLobal projected nodes mapping projected -> unprojected mapping
    Vector<int> gpN;

    /// @brief Global connectivity array mappig eNoN,nEl --> gnNo
    Array<int> gIEN;

    /// @brief The connectivity array mapping eNoN,nEl --> nNo
    Array<int> IEN;

    /// @brief gIEN mapper from old to new
    Vector<int> otnIEN;

    /// @brief Local knot pointer (NURBS)
    Array<int> INN;

    /// @brief Global to local maping tnNo --> nNo
    Vector<int> lN;

    /// @brief Shells: extended IEN array with neighboring nodes
    Array<int> eIEN;

    /// @brief Shells: boundary condition variable
    Array<int> sbc;

    /// @brief IB: Whether a cell is a ghost cell or not
    Vector<int> iGC;

    /// @brief Control points weights (NURBS)
    Vector<double> nW;

    /// @brief Gauss weights
    Vector<double> w;

    /// @brief Gauss integration points in parametric space
    Array<double> xi;

    /// @brief Bounds on parameteric coordinates
    Array<double> xib;

    /// @brief Position coordinates (not always, however, as they get overwritten by read_vtu_pdata())
    Array<double> x;

    /// @brief Parent shape function
    Array<double> N;

    /// @brief Shape function bounds
    Array<double> Nb;

    /// @brief Normal vector to each nodal point (for Shells)
    Array<double> nV;

    /// @brief Fiber orientations stored at the element level - used for
    /// electrophysiology and solid mechanics
    Array<double> fN;

    /// @brief Parent shape functions gradient
    /// double Nx(:,:,:)
    Array3<double> Nx;

    /// @brief Second derivatives of shape functions - used for shells & IGA
    /// davep double Nxx(:,:,:)
    Array3<double> Nxx;

    /// @brief Solution field (displacement, velocity, pressure, etc.) for a known, potentially
    /// time-varying, quantity of interest across a mesh
    Array3<double> Ys;

    /// @brief Mesh Name
    std::string name;

    /// @brief Mesh nodal adjacency
    adjType nAdj;

    /// @brief Mesh element adjacency
    adjType eAdj;

    /// @brief Function spaces (basis)
    std::vector<fsType> fs;

    /// @brief BSpline in different directions (NURBS)
    std::vector<bsType> bs;

    /// @brief Faces are stored in this variable
    std::vector<faceType> fa;

    /// @brief IB: tracers
    traceType trc;

    /// @brief RIS: flags of whether elemets are adjacent to RIS projections
    // std::vector<bool> eRIS;
    Vector<int> eRIS;

    /// @brief RIS: processor ids to change element partitions to
    Vector<int> partRIS;

    /// @brief TET4 quadrature modifier
    double qmTET4 = (5.0+3.0*sqrt(5.0))/20.0;

  private:
    //mshType(const mshType&);
    //mshType& operator=(const mshType&);

};

/// @brief Equation type
//
class eqType
{
  public:
    eqType();
    ~eqType();

    /// @brief Should be satisfied in a coupled/uncoupled fashion
    bool coupled = false;
    //bool coupled = .TRUE.

    /// @brief Satisfied/not satisfied
    bool ok = false;

    /// @brief Use C++ Trilinos framework for the linear solvers
    bool useTLS = false;

    /// @brief Use C++ Trilinos framework for assembly and for linear solvers
    bool assmTLS = false;

    /// @brief Degrees of freedom
    int dof = 0;

    /// @brief Pointer to end of unknown Yo(:,s:e)
    int e = -1;

    /// @brief Pointer to start of unknown Yo(:,s:e)
    int s = -1;

    /// @brief Number of performed iterations
    int itr = 0;

    /// @brief Maximum iteration for this eq.
    int maxItr = 5;

    /// @brief Minimum iteration for this eq.
    int minItr = 1;

    /// @brief Number of possible outputs
    int nOutput = 0;

    /// @brief IB: Number of possible outputs
    int nOutIB = 0;

    /// @brief URIS: Number of possible outputs
    int nOutURIS = 0;

    /// @brief Number of domains
    int nDmn = 0;

    /// @brief IB: Number of immersed domains
    int nDmnIB = 0;

    /// @brief Number of BCs
    int nBc = 0;

    /// @brief Number of BCs on immersed surfaces
    int nBcIB = 0;

    /// @brief Number of BFs
    int nBf = 0;

    /// @brief Type of equation fluid/heatF/heatS/lElas/FSI
    consts::EquationType phys = consts::EquationType::phys_NA;

    // Parameters used for the Generalized α− Method.
    //
    /// @brief \f$\alpha_f\f$
    double af = 0.0;

    /// @brief \f$\alpha_m\f$
    ///
    /// For second order equations: am = (2.0 - roInf) / (1.0 + roInf)
    /// First order equations: am = 0.5 * (3.0 - roInf) / (1.0 + roInf)
    //
    double am = 0.0;

    /// @brief \f$\beta\f$
    double beta = 0.0;

    /// @brief \f$\gamma\f$
    double gam = 0.0;

    /// @brief Initial norm of residual
    double iNorm = 0.0;

    /// @brief First iteration norm
    double pNorm = 0.0;

    /// @brief \f$\rho_{infinity}\f$
    double roInf = 0.0;

    /// @brief Accepted relative tolerance
    double tol = 0.0;

    /// @brief Equation symbol
    std::string sym;
    //std::string(LEN=2) :: sym = "NA";

    /// @brief type of linear solver
    lsType ls;

    /// @brief The type of interface to a numerical linear algebra library.
    consts::LinearAlgebraType linear_algebra_type;

    /// @brief The type of assembly interface to a numerical linear algebra library.
    consts::LinearAlgebraType linear_algebra_assembly_type;

    /// @brief The type of preconditioner used by the interface to a numerical linear algebra library.
    consts::PreconditionerType linear_algebra_preconditioner = consts::PreconditionerType::PREC_FSILS;

    /// @brief Interface to a numerical linear algebra library.
    LinearAlgebra* linear_algebra = nullptr;

    /// @brief FSILS type of linear solver
    fsi_linear_solver::FSILS_lsType FSILS;

    /// @brief BCs associated with this equation;
    std::vector<bcType> bc;

    /// @brief IB: BCs associated with this equation on immersed surfaces
    std::vector<bcType> bcIB;

    /// @brief domains that this equation must be solved
    std::vector<dmnType> dmn;

    /// @brief IB: immersed domains that this equation must be solved
    std::vector<dmnType> dmnIB;

    /// @brief Outputs
    std::vector<outputType> output;

    /// @brief IB: Outputs
    std::vector<outputType> outIB;

    /// @brief URIS: Outputs
    std::vector<outputType> outURIS;

    /// @brief Body force associated with this equation
    std::vector<bfType> bf;
};

/// @brief This type will be used to write data in the VTK files.
//
class dataType
{
  public:

    //  Element number of nodes
    int eNoN = 0;

    //  Number of elements
    int nEl = 0;

    //  Number of nodes
    int nNo = 0;

    //  vtk type
    int vtkType = 0;

    //  Connectivity array
    Array<int> IEN;

    //  Element based variables to be written
    Array<double> xe;

    //  All the variables after transformation to global format
    Array<double> gx;

    //  All the variables to be written (including position)
    Array<double> x;
};


class rmshType
{
  public:

    rmshType();

    /// @brief Whether remesh is required for problem or not
    bool isReqd = false;

    /// @brief Method for remeshing: 1-TetGen, 2-MeshSim
    consts::MeshGeneratorType method = consts::MeshGeneratorType::RMSH_TETGEN;

    /// @brief Counter to track number of remesh done
    int cntr = 0;

    /// @brief Time step from which remeshing is done
    int rTS = 0;

    /// @brief Time step freq for saving data
    int cpVar = 0;

    /// @brief Time step at which forced remeshing is done
    int fTS = 1000;

    /// @brief Time step frequency for forced remeshing
    int freq = 1000;

    /// @brief Time where remeshing starts
    double time = 0.0;

    /// @brief Mesh quality parameters
    double minDihedAng = 0.0;
    double maxRadRatio = 0.0;

    /// @brief Edge size of mesh
    Vector<double> maxEdgeSize;

    /// @brief Initial norm of an equation
    Vector<double> iNorm;

    /// @brief Copy of solution variables where remeshing starts
    Array<double> A0;
    Array<double> Y0;
    Array<double> D0;

    /// @brief Flag is set if remeshing is required for each mesh
    std::vector<bool> flag;
};

class ibCommType
{
  public:
    /// @brief Num traces (nodes) local to each process
    Vector<int> n;

    /// @brief Pointer to global trace (node num) stacked contiguously
    Vector<int> gN;

    /// @brief Num traces (Gauss points) local to each process
    Vector<int> nG;

    /// @brief Pointer to global trace (Gauss point) stacked contiguously
    Vector<int> gE;
};


/// @brief Immersed Boundary (IB) data type
//
class ibType
{
  public:

    /// @brief  Whether any file being saved
    bool savedOnce = false;
    //bool savedOnce = .FALSE.

    /// @brief IB method
    int mthd = 0;
    //int mthd = ibMthd_NA;

    /// @brief IB coupling
    int cpld = 0;
    //int cpld = ibCpld_NA;

    /// @brief IB interpolation method
    int intrp = 0;
    //int intrp = ibIntrp_NA;

    /// @brief Current IB domain ID
    int cDmn = 0;

    /// @brief Current equation
    int cEq = 0;

    /// @brief Total number of IB nodes
    int tnNo = 0;

    /// @brief Number of IB meshes
    int nMsh = 0;

    /// @brief IB call duration (1: total time; 2: update; 3,4: communication)
    double callD[4] = {0.0, 0.0, 0.0, 0.0};

    /// @brief IB Domain ID
    Vector<int> dmnID;

    /// @brief Row pointer (for sparse LHS matrix storage)
    Vector<int> rowPtr;

    /// @brief Column pointer (for sparse LHS matrix storage)
    Vector<int> colPtr;

    /// @brief IB position coordinates
    Array<double> x;

    /// @brief Velocity (new)
    Array<double> Yb;

    /// @brief Time derivative of displacement (old)
    Array<double> Auo;

    /// @brief Time derivative of displacement (new)
    Array<double> Aun;

    /// @brief Time derivative of displacement (n+am)
    Array<double> Auk;

    /// @brief Displacement (old)
    Array<double> Ubo;

    /// @brief Displacement (new)
    Array<double> Ubn;

    /// @brief Displacement (n+af)
    Array<double> Ubk;

    /// @brief Displacement (projected on background mesh, old)
    Array<double> Uo;

    /// @brief Displacement (projected on background mesh, new, n+af)
    Array<double> Un;

    /// @brief Residual (FSI force)
    Array<double> R;

    /// @brief Residual (displacement, background mesh)
    Array<double> Ru;

    /// @brief Residual (displacement, IB mesh)
    Array<double> Rub;

    /// @brief LHS tangent matrix for displacement
    Array<double> Ku;

    /// @brief DERIVED class VARIABLES IB meshes;
    std::vector<mshType> msh;

    /// @brief IB communicator
    ibCommType cm;
};

/// @brief Data type for Resistive Immersed Surface
//
class risFaceType 
{
  public:

    /// @brief Number of RIS surface
    int nbrRIS = 0;

    /// @brief Count time steps where no check is needed
    Vector<int> nbrIter;

    /// @brief List of meshes, and faces connected. The first face is the 
    // proximal pressure's face, while the second is the distal one
    Array3<int> lst;

    /// @brief Resistance value 
    Vector<double> Res;

    /// @brief Flag closed surface active, the valve is considerd open initially 
    std::vector<bool> clsFlg;

    /// @brief Mean distal and proximal pressure (1: distal, 2: proximal)
    Array<double> meanP;

    /// @brief Mean flux on the RIS surface 
    Vector<double> meanFl;

    /// @brief Status RIS interface
    std::vector<bool> status;
};

/// @brief Unfitted Resistive Immersed surface data type
//
class urisType 
{
  public:

    // Name of the URIS instance.
    std::string name;

    // Whether any file has been saved.
    bool savedOnce = false;

    // Total number of IB nodes.
    int tnNo = 0;

    // Number of IB meshes.
    int nFa = 0;

    // Position coordinates (2D array: rows x columns).
    Array<double> x;

    // Displacement (new) (2D array).
    Array<double> Yd;

    // Default signed distance value away from the valve.
    double sdf_default = 1000.0;

    // Default distance value of the valve boundary (valve thickness).
    double sdf_deps = 0.25;

    // Default distance value of the valve boundary when the valve is closed.
    double sdf_deps_close = 0.25;

    // Displacements of the valve when it opens (3D array).
    Array3<double> DxOpen;

    // Displacements of the valve when it closes (3D array).
    Array3<double> DxClose;

    // Normal vector pointing in the positive flow direction (1D array).
    Vector<double> nrm;

    // Close flag.
    bool clsFlg = true;

    // Iteration count.
    int cnt = 1000000;

    // URIS: signed distance function of each node to the uris (1D array).
    Vector<double> sdf;

    // Mesh scale factor.
    double scF = 1.0;

    // Mean pressure upstream.
    double meanPU = 0.0;

    // Mean pressure downstream.
    double meanPD = 0.0;

    // Relaxation factor to compute weighted averages of pressure values.
    double relax_factor = 0.5;

    // Array to store the fluid mesh elements that the uris node is in (2D array).
    Array<int> elemId;

    // Array to count how many times a uris node is found in the fluid mesh of a processor (1D array).
    Vector<int> elemCounter;

    // Derived type variables
    // IB meshes
    std::vector<mshType> msh;

};

/// @brief The ComMod class duplicates the data structures in the Fortran COMMOD module
/// defined in MOD.f. 
///
/// The data members here are the global variables exposed by the COMMOD module.
//
class ComMod {

  public:
    ComMod();
    ~ComMod();

    //----- bool members -----//

    /// @brief Whether there is a requirement to update mesh and Dn-Do variables
    bool dFlag = false;

    /// @brief Whether mesh is moving
    bool mvMsh = false;

    /// @brief Whether to averaged results
    bool saveAve = false;

    /// @brief Whether to save to VTK files
    bool saveVTK = false;

    /// @brief Whether any file being saved
    bool savedOnce = false;

    /// @brief Whether to use separator in output
    bool sepOutput = false;

    /// @brief Whether start from beginning or from simulations
    bool stFileFlag = false;

    /// @brief Whether to overwrite restart file or not
    bool stFileRepl = false;

    /// @brief Restart simulation after remeshing
    bool resetSim = false;

    /// @brief Check IEN array for initial mesh
    bool ichckIEN = false;

    /// @brief Reset averaging variables from zero
    bool zeroAve = false;

    /// @brief Whether CMM equation is initialized
    bool cmmInit = false;

    /// @brief Whether variable wall properties are used for CMM
    bool cmmVarWall = false;

    /// @brief Whether shell equation is being solved
    bool shlEq = false;

    /// @brief Whether PRESTRESS is being solved
    bool pstEq = false;

    /// @brief Whether velocity-pressure based structural dynamics solver is used
    bool sstEq = false;

    /// @brief Whether to detect and apply any contact model
    bool iCntct = false;

    /// @brief Whether any Immersed Boundary (IB) treatment is required
    bool ibFlag = false;

    /// @brief Postprocess step - convert bin to vtk
    bool bin2VTK = false;

    /// @brief Whether any RIS surface is considered 
    bool risFlag = false;

    /// @brief Whether any one-sided RIS surface with 0D coupling is considered 
    bool ris0DFlag = false;

    /// @brief Whether any URIS surface is considered
    bool urisFlag = false;

    /// @brief Whether the URIS surface is active
    bool urisActFlag = false;

    /// @brief Number of URIS surfaces (uninitialized, to be set later)
    int nUris;

    /// @brief URIS resistance
    double urisRes;

    /// @brief URIS resistance when the valve is closed
    double urisResClose;

    /// @brief Whether to use precomputed state-variable solutions
    bool usePrecomp = false;
    //----- int members -----//

    /// @brief Current domain
    int cDmn = 0;

    /// @brief Current equation
    int cEq = 0;

    /// @brief Current time step
    int cTS = 0;

    std::array<double,3> timeP;

    /// @brief Starting time step
    int startTS = 0;

    /// @brief Current equation degrees of freedom
    int dof = 0;

    /// @brief Global total number of nodes, across all meshes (total) and all 
    /// procs (global)
    int gtnNo = 0;

    /// @brief Number of equations
    int nEq = 0;

    /// @brief Number of faces in the LHS passed to FSILS
    int nFacesLS = 0;

    /// @brief Number of meshes
    int nMsh = 0;

    /// @brief Number of spatial dimensions
    int nsd = 0;

    /// @brief Number of time steps
    int nTS = 0;

    /// @brief Number of initialization time steps
    int nITs = 0;

    /// @brief stFiles record length
    int recLn = 0;

    /// @brief Start saving after this number of time step
    int saveATS = 0;

    /// @brief Increment in saving solutions
    int saveIncr = 0;

    /// @brief Stamp ID to make sure simulation is compatible with stFiles
    std::array<int,7> stamp;

    /// @brief Increment in saving restart file
    int stFileIncr = 0;

    /// @brief Total number of degrees of freedom per node
    int tDof = 0;

    /// @brief Total number of nodes (number of nodes on current proc across
    /// all meshes)
    int tnNo = 0;

    /// @brief Restart Time Step
    int rsTS = 0;

    /// @brief Number of stress values to be stored
    int nsymd = 0;

    /// @brief Nbr of iterations 
    int RisnbrIter = 0;


    //----- double members -----//

    /// @brief Time step size
    double dt = 0.0;

    /// @brief Time step size of the precomputed state-variables
    double precompDt = 0.0;

    /// @brief Time
    double time = 0.0;


    //----- string members -----//

    /// @brief Initialization file path
    std::string iniFilePath;

    /// @brief Saved output file name
    std::string saveName;

    /// @brief Restart file name
    std::string stFileName;

    /// @brief Stop_trigger file name
    std::string stopTrigName;

    /// @brief Precomputed state-variable file name
    std::string precompFileName;

    /// @brief Precomputed state-variable field name
    std::string precompFieldName;
    // ALLOCATABLE DATA

    /// @brief Column pointer (for sparse LHS matrix structure)
    /// Modified in: lhsa()
    Vector<int> colPtr;

    /// @brief Domain ID
    Vector<int>  dmnId;

    /// @brief Local to global pointer tnNo --> gtnNo
    Vector<int> ltg;

    /// @brief Row pointer (for sparse LHS matrix structure)
    /// Modified in: lhsa()
    Vector<int> rowPtr;

    /// @brief Array that maps global node id to rowN in the matrix
    /// Modified in: lhsa()
    Vector<int> idMap;

    /// @brief Boundary nodes set for CMM initialization and for zeroing-out
    /// non-wall nodal displacements
    Vector<int> cmmBdry;

    /// @brief IB: iblank used for immersed boundaries (1 => solid, 0 => fluid)
    Vector<int> iblank;

    /// @brief TODO: for now, better to organize these within a class      
    struct Array2D {
        // std::vector<std::vector<int>> map;
        Array<int> map;
    };

    /// @brief RIS mapping array, with local (mesh) enumeration
     std::vector<Array2D> risMapList;

    /// @brief RIS mapping array, with global (total) enumeration
     std::vector<Array2D> grisMapList;

    /// @brief Old time derivative of variables (acceleration); known result at current time step
    Array<double>  Ao;

    /// @brief New time derivative of variables (acceleration); unknown result at next time step
    Array<double>  An;

    /// @brief Old integrated variables (displacement)
    Array<double>  Do;

    /// @brief New integrated variables (displacement)
    Array<double>  Dn;

    /// @brief Residual vector
    Array<double>  R;

    /// @brief LHS matrix
    Array<double>  Val;

    /// @brief Position vector of mesh nodes (in ref config)
    Array<double>  x;

    /// @brief Old variables (velocity); known result at current time step
    Array<double>  Yo;

    /// @brief New variables (velocity); unknown result at next time step
    Array<double>  Yn;

    /// @brief Body force
    Array<double>  Bf;

    //-----------------------------------------------------
    // Additional arrays for velocity-based formulation of 
    // nonlinear solid mechanics.
    //-----------------------------------------------------

    /// @brief Time derivative of displacement
    Array<double>  Ad;

    /// @brief Residual of the displacement equation
    Array<double>  Rd;

    /// @brief LHS matrix for displacement equation
    Array<double>  Kd;

    /// @brief Variables for prestress calculations
    Array<double>  pS0;
    Array<double>  pSn;
    Vector<double>  pSa;

    /// @brief Temporary storage for initializing state variables
    Vector<double> Pinit;
    Array<double>  Vinit;
    Array<double>  Dinit;

    /// @brief CMM-variable wall properties: 1-thickness, 2-Elasticity modulus
    Array<double>  varWallProps;

    //------------------------
    // DERIVED TYPE VARIABLES
    //------------------------

    /// @brief Coupled BCs structures used for multidomain simulations
    cplBCType cplBC;

    /// @brief All data related to equations are stored in this container
    std::vector<eqType> eq;

    /// @brief FSILS data structure to produce LHS sparse matrix
    fsi_linear_solver::FSILS_lhsType lhs;

    /// @brief All the meshes are stored in this variable
    std::vector<mshType> msh;

    /// @brief Input/output to the screen is handled by this structure
    chnlType std, err, wrn, dbg;

    /// @brief To group above channels
    ioType io;

    /// @brief The general communicator
    cmType cm;

    /// @brief Remesher type
    rmshType rmsh;

    /// @brief Contact model type
    cntctModelType cntctM;

    /// @brief IB: Immersed boundary data structure
    ibType ib;

    /// @brief risFace object
    risFaceType ris;

    /// @brief unfitted RIS object
    std::vector<urisType> uris;

    bool debug_active = false;

    Timer timer;
};

#endif

