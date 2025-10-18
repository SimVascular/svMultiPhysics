// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "CmMod.h"

#include "mpi.h"

#include <iostream>

CmMod::CmMod()
{
}

CmMod::~CmMod()
{
}

//////////////////////////////////////////////////////////////////////////
//                c m T y p e    i m p l e m e n t a t i o n            //
//////////////////////////////////////////////////////////////////////////

cmType::cmType()
{
}

cmType::~cmType()
{
}

/// @brief Returns commu handle
decltype(MPI_COMM_WORLD) cmType::com() const
{  
  return cHndl; 
}

void cmType::new_cm(decltype(MPI_COMM_WORLD) comHandle)
{
  cHndl  = comHandle;
  taskId = 0;
  nProcs = 1;
  int ierr;
  ierr = MPI_Comm_rank(comHandle, &taskId);
  ierr = MPI_Comm_size(comHandle, &nProcs);
  nThreads = 1;

  if (nProcs == 1) {
    return;
  }
}

/// @brief bcast bool.
void cmType::bcast(const CmMod& cm_mod, bool* data) const
{
  MPI_Bcast(data, 1, cm_mod::mplog, cm_mod.master, com());
}

/// @brief bcast bool array
void cmType::bcast(const CmMod& cm_mod, std::vector<bool>& data) const
{
  static const int MAX_ARRAY_SIZE = 100;

  if (data.size() == 0) {
    throw std::runtime_error("bcast bool array has size 0.");
  }

  if (data.size() > MAX_ARRAY_SIZE) {
    throw std::runtime_error("bcast bool array is larger than " + std::to_string(MAX_ARRAY_SIZE) + ".");
  }

  bool bool_array[MAX_ARRAY_SIZE];
  int n = data.size();
  for (int i = 0; i < n; i++) {
    bool_array[i] = data[i];
  }

  MPI_Bcast(bool_array, n, cm_mod::mplog, cm_mod.master, com());

  data.clear();
  for (int i = 0; i < n; i++) {
    data.push_back(bool_array[i]);
  }
}

/// @brief bcast char*
void cmType::bcast(const CmMod& cm_mod, std::string& data) const
{
  static const int MAX_ARRAY_SIZE = 400;

  if (data.size() > MAX_ARRAY_SIZE) {
    throw std::runtime_error("bcast string is larger than " + std::to_string(MAX_ARRAY_SIZE) + ".");
  }
  static char buffer[MAX_ARRAY_SIZE]{};

  strcpy(buffer, data.c_str());

  MPI_Bcast(buffer, MAX_ARRAY_SIZE, cm_mod::mpchar, cm_mod.master, com());

  data = buffer;
}

/// @brief bcast double 
void cmType::bcast(const CmMod& cm_mod, double* data) const
{
  MPI_Bcast(data, 1, cm_mod::mpreal, cm_mod.master, com());
}

/// @brief bcast double array
void cmType::bcast(const CmMod& cm_mod, Array<double>& data, const std::string& name) const
{
  MPI_Bcast(data.data(), data.size(), cm_mod::mpreal, cm_mod.master, com());
}

/// @brief bcast int array
void cmType::bcast(const CmMod& cm_mod, Array<int>& data, const std::string& name) const
{
  MPI_Bcast(data.data(), data.size(), cm_mod::mpint, cm_mod.master, com());
}

/// @brief bcast double Vector 
void cmType::bcast(const CmMod& cm_mod, Vector<double>& data, const std::string& name) const
{
  if (data.size() == 0) {
    throw std::runtime_error(name + ":bcast double vector has size 0.");
  }
  MPI_Bcast(data.data(), data.size(), cm_mod::mpreal, cm_mod.master, com());
}

/// @brief bcast int 
void cmType::bcast(const CmMod& cm_mod, int* data) const
{
  MPI_Bcast(data, 1, cm_mod::mpint, cm_mod.master, com());
}

/// @brief bcast int Vector
void cmType::bcast(const CmMod& cm_mod, Vector<int>& data) const
{
  MPI_Bcast(data.data(), data.size(), cm_mod::mpint, cm_mod.master, com());
}

/// @brief gather int array
void cmType::gather(const CmMod& cm_mod, const int* send_data, int send_count, int* recv_data, int recv_count, int root) const
{
  MPI_Gather(send_data, send_count, cm_mod::mpint, recv_data, recv_count, cm_mod::mpint, root, com());
}

/// @brief gather double array
void cmType::gather(const CmMod& cm_mod, const double* send_data, int send_count, double* recv_data, int recv_count, int root) const
{
  MPI_Gather(send_data, send_count, cm_mod::mpreal, recv_data, recv_count, cm_mod::mpreal, root, com());
}

/// @brief gather int Vector
void cmType::gather(const CmMod& cm_mod, const Vector<int>& send_data, Vector<int>& recv_data, int root) const
{
  MPI_Gather(send_data.data(), send_data.size(), cm_mod::mpint, 
             recv_data.data(), send_data.size(), cm_mod::mpint, root, com());
}

/// @brief gather double Vector
void cmType::gather(const CmMod& cm_mod, const Vector<double>& send_data, Vector<double>& recv_data, int root) const
{
  MPI_Gather(send_data.data(), send_data.size(), cm_mod::mpreal, 
             recv_data.data(), send_data.size(), cm_mod::mpreal, root, com());
}

/// @brief scatter int array
void cmType::scatter(const CmMod& cm_mod, const int* send_data, int send_count, int* recv_data, int recv_count, int root) const
{
  MPI_Scatter(send_data, send_count, cm_mod::mpint, recv_data, recv_count, cm_mod::mpint, root, com());
}

/// @brief scatter double array
void cmType::scatter(const CmMod& cm_mod, const double* send_data, int send_count, double* recv_data, int recv_count, int root) const
{
  MPI_Scatter(send_data, send_count, cm_mod::mpreal, recv_data, recv_count, cm_mod::mpreal, root, com());
}

/// @brief scatter int Vector
void cmType::scatter(const CmMod& cm_mod, const Vector<int>& send_data, Vector<int>& recv_data, int root) const
{
  MPI_Scatter(send_data.data(), send_data.size() / nProcs, cm_mod::mpint, 
              recv_data.data(), send_data.size() / nProcs, cm_mod::mpint, root, com());
}

/// @brief scatter double Vector
void cmType::scatter(const CmMod& cm_mod, const Vector<double>& send_data, Vector<double>& recv_data, int root) const
{
  MPI_Scatter(send_data.data(), send_data.size() / nProcs, cm_mod::mpreal, 
              recv_data.data(), send_data.size() / nProcs, cm_mod::mpreal, root, com());
}

/// @brief gatherv int Vector
void cmType::gatherv(const CmMod& cm_mod, const Vector<int>& send_data, Vector<int>& recv_data, const Vector<int>& recv_counts, const Vector<int>& displs, int root) const
{
  MPI_Gatherv(send_data.data(), send_data.size(), cm_mod::mpint,
              recv_data.data(), recv_counts.data(), displs.data(), cm_mod::mpint, root, com());
}

/// @brief gatherv double Vector
void cmType::gatherv(const CmMod& cm_mod, const Vector<double>& send_data, Vector<double>& recv_data, const Vector<int>& recv_counts, const Vector<int>& displs, int root) const
{
  MPI_Gatherv(send_data.data(), send_data.size(), cm_mod::mpreal,
              recv_data.data(), recv_counts.data(), displs.data(), cm_mod::mpreal, root, com());
}

/// @brief scatterv int Vector
void cmType::scatterv(const CmMod& cm_mod, const Vector<int>& send_data, const Vector<int>& send_counts, const Vector<int>& displs, Vector<int>& recv_data, int root) const
{
  MPI_Scatterv(send_data.data(), send_counts.data(), displs.data(), cm_mod::mpint,
               recv_data.data(), recv_data.size(), cm_mod::mpint, root, com());
}

/// @brief scatterv double Vector
void cmType::scatterv(const CmMod& cm_mod, const Vector<double>& send_data, const Vector<int>& send_counts, const Vector<int>& displs, Vector<double>& recv_data, int root) const
{
  MPI_Scatterv(send_data.data(), send_counts.data(), displs.data(), cm_mod::mpreal,
               recv_data.data(), recv_data.size(), cm_mod::mpreal, root, com());
}