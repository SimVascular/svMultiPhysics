// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The classes defined here duplicate the data structures in the Fortran CMMOD module
// defined in COMU.f. 

#ifndef CMMOD_H 
#define CMMOD_H 

#include "Array.h"
#include "Vector.h"

#include "mpi.h"
#include "consts.h"

#include <string>

namespace cm_mod {

using MpiCommWorldType = MPI_Comm; 
//using MpiCommWorldType = decltype(MPI_COMM_WORLD);

// Set MPI data type names.
const decltype(MPI_CXX_BOOL) mplog = MPI_CXX_BOOL;
//const decltype(MPI_LOGICAL) mplog  = MPI_LOGICAL;
const decltype(MPI_INTEGER) mpint = MPI_INTEGER;
const decltype(MPI_DOUBLE_PRECISION) mpreal = MPI_DOUBLE_PRECISION;
const decltype(MPI_CHARACTER) mpchar = MPI_CHARACTER;
};

/// @brief The CmMod class duplicates the data structures in the Fortran CMMOD module
/// defined in COMU.f. 
///
/// The data members here are the global variables exposed by the CMMOD module.
class CmMod {

  public:
    CmMod();
    ~CmMod();

    // Size of blocks for openMP communications.
    int mpBs = 1000;

    // master is assumed to have zero ID
    int master = 0;

    //  Abstracted MPI names.
    decltype(MPI_LOGICAL) mplog  = MPI_LOGICAL;
    decltype(MPI_INTEGER) mpint = MPI_INTEGER;
    decltype(MPI_DOUBLE_PRECISION) mpreal = MPI_DOUBLE_PRECISION;
    decltype(MPI_CHARACTER) mpchar = MPI_CHARACTER;
};

/// @brief The cmType class stores data and defines methods used for mpi communication.
class cmType {
  public:
    cmType();
    ~cmType();

    // Communicator handle.
    decltype(MPI_COMM_WORLD) cHndl;

    // Processors ID.
    int taskId = 0;

    // Number of openMP threads in this cm
    int nThreads = 0;

    // Number of processors
    int nProcs = 0;

    //----- M e t h o d s -----//

    void bcast(const CmMod& cm_mod, bool* data) const;
    void bcast(const CmMod& cm_mod, std::vector<bool>& data) const;

    void bcast(const CmMod& cm_mod, std::string& data) const;

    void bcast(const CmMod& cm_mod, double* data) const;
    void bcast(const CmMod& cm_mod, Vector<double>& data, const std::string& name="") const;
    void bcast(const CmMod& cm_mod, Array<double>& data, const std::string& name="") const;

    void bcast(const CmMod& cm_mod, int* data) const;
    void bcast(const CmMod& cm_mod, Vector<int>& data) const;

    void bcast(const CmMod& cm_mod, Array<int>& data, const std::string& name="") const;

    // Gather operations
    void gather(const CmMod& cm_mod, const int* send_data, int send_count, int* recv_data, int recv_count, int root) const;
    void gather(const CmMod& cm_mod, const double* send_data, int send_count, double* recv_data, int recv_count, int root) const;
    void gather(const CmMod& cm_mod, const Vector<int>& send_data, Vector<int>& recv_data, int root) const;
    void gather(const CmMod& cm_mod, const Vector<double>& send_data, Vector<double>& recv_data, int root) const;

    // Gatherv operations
    void gatherv(const CmMod& cm_mod, const Vector<int>& send_data, Vector<int>& recv_data, const Vector<int>& recv_counts, const Vector<int>& displs, int root) const;
    void gatherv(const CmMod& cm_mod, const Vector<double>& send_data, Vector<double>& recv_data, const Vector<int>& recv_counts, const Vector<int>& displs, int root) const;

    // Scatterv operations
    void scatterv(const CmMod& cm_mod, const Vector<int>& send_data, const Vector<int>& send_counts, const Vector<int>& displs, Vector<int>& recv_data, int root) const;
    void scatterv(const CmMod& cm_mod, const Vector<double>& send_data, const Vector<int>& send_counts, const Vector<int>& displs, Vector<double>& recv_data, int root) const;

    // Scatter operations
    void scatter(const CmMod& cm_mod, const int* send_data, int send_count, int* recv_data, int recv_count, int root) const;
    void scatter(const CmMod& cm_mod, const double* send_data, int send_count, double* recv_data, int recv_count, int root) const;
    void scatter(const CmMod& cm_mod, const Vector<int>& send_data, Vector<int>& recv_data, int root) const;
    void scatter(const CmMod& cm_mod, const Vector<double>& send_data, Vector<double>& recv_data, int root) const;

    //------------
    // bcast_enum
    //------------
    //
    template <typename T>
    void bcast_enum(const CmMod& cm_mod, T* data) const
    {
      int idata = static_cast<int>(*data);
      //std::cout << "[bcast_enum] idata in: " << idata << std::endl;
      MPI_Bcast(&idata, 1, cm_mod::mpint, cm_mod.master, com());
      //std::cout << "[bcast_enum] idata out: " << idata << std::endl;
      *data = static_cast<T>(idata);
    }

    //------------
    // bcast_prop
    //------------
    //
    template <typename T>
    void bcast_prop(const CmMod& cm_mod, std::map<T,double>& props) const
    {
      static const int MAX_SIZE = 100;

      if (2*props.size() > MAX_SIZE) {
        throw std::runtime_error("bcast prop is larger than " + std::to_string(MAX_SIZE) + ".");
      }

      double prop_array[MAX_SIZE];
      std::fill_n(prop_array, MAX_SIZE, -1.0);

      int n = 0;
      for (auto& entry : props) {
        prop_array[n++] = static_cast<int>(entry.first);
        prop_array[n++] = entry.second;
      }

      MPI_Bcast(prop_array, MAX_SIZE, cm_mod::mpreal, cm_mod.master, com());

      props.clear();
      int num_props = MAX_SIZE / 2;;

      for (int i = 0; i < num_props; i++) {
        int iprop = static_cast<int>(prop_array[2*i]);
        if (iprop == -1) { 
          break;
        }
        auto prop = static_cast<T>(iprop);
        props[prop] = prop_array[2*i+1];
      }
    }

    // Returns commu handle
    cm_mod::MpiCommWorldType com() const;
    //decltype(MPI_COMM_WORLD) com() const;

    int idcm() const { return taskId; };
    int id() { return taskId; };
    bool mas(const CmMod& cm_mod) const { return (taskId == cm_mod.master); }; 

    // Create a new Communicator
    void new_cm(decltype(MPI_COMM_WORLD) comHandle);

    int np() const { return nProcs; }; 

    int nT() { return nThreads; };


    //--------
    // reduce
    //--------
    // For an int or double scalar.
    //
    template<typename T>
    T reduce(const CmMod& cm_mod, T u, MPI_Op op = MPI_SUM) const
    {
      T gU{};

      MPI_Datatype data_type;
      if (typeid(T) == typeid(double)) {
        data_type = MPI_DOUBLE_PRECISION;
      } else if (typeid(T) == typeid(int)) {
        data_type = MPI_INTEGER;
      } else {
        throw std::runtime_error("[cm_mod::reduce called with unknown data type.");
      }

      if (seq()) {
        gU = u;
      } else {
        MPI_Allreduce(&u, &gU, 1, data_type, op, com());
      }

      return gU;
    }

    //--------
    // reduce
    //--------
    // For an int or double Vector.
    //
    template<typename T>
    Vector<T> reduce(const CmMod& cm_mod, Vector<T>& u, MPI_Op op = MPI_SUM) const
    {
      int size = u.size();
      Vector<T> gU(size);

      MPI_Datatype data_type;
      if (typeid(T) == typeid(double)) {
        data_type = MPI_DOUBLE_PRECISION;
      } if (typeid(T) == typeid(int)) {
        data_type = MPI_INTEGER;
      }

      if (seq()) {
        gU = u;
      } else {
        MPI_Allreduce(u.data(), gU.data(), size, data_type, op, com());
      }

      return gU;
    }

    bool seq() const { return (nProcs == 1); };

    bool slv(const CmMod& cm_mod) const { return (taskId != cm_mod.master); };

    // Returns processor ID in fortran indexing
    int tF(const CmMod& cm_mod) const { return taskId + 1; };

};


#endif

