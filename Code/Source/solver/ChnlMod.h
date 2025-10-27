// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The classes defined here duplicate the data structures in the Fortran CHNLMOD module
// defined in CHNL.f. 

// This module defines data structures for 

#include "Array.h"
#include "Vector.h"

#include <string>

#ifndef CHNL_MOD_H 
#define CHNL_MOD_H 

/// @brief Channel type, used in I/O
class chnlType
{
  // Whether it is open to the screen
  bool oTS = false;

  // Whether it is open to the file
  bool oTF = false; 

  // Channel identifier
  int id;

  // File ID
  int fId;

  // File address
  std::string fName = "histor";

  // Channel tag
  std::string tag = "";

  //         Creates a new channel
  //       PROCEDURE :: new => newChnl
  //         Closes the channel
  //  PROCEDURE :: close => closeChnl
  //         To send a string to channel
   //        PROCEDURE chnlAssign
      //     GENERIC :: ASSIGNMENT(=) => chnlAssign
};

/// @brief Only to group four channels, in case I rather have them as one
/// variable
class ioType
{
  // Standard output
  chnlType o;

  // Error
  chnlType e;

  // Warning
  chnlType w;

  // Debugging
  chnlType d;

  // Status file
  chnlType s;

  //      CONTAINS
  //!        Opens all as standard channels
           //PROCEDURE :: new => newIO
  //!        Closes the channel
           //PROCEDURE :: close => closeIO
};

class ChnlMod
{
  public:
    ChnlMod();

    // Channels cases: standard output, error output, warning output,
    // debugging output
    int CHNL_O = 601;
    int CHNL_E = 602;
    int CHNL_W = 603;
    int CHNL_D = 604;

    // Whether to use color in printing outputs
    bool pClr = true;

    // A general counter for file ID
    int gFID = 314;

    // Appended path to all files that are going to be saved
    std::string appPath = "";

};

#endif

