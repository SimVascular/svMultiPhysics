#!/bin/bash

c_solver=/Users/parkerda/software/ktbolt/svMultiPhysics/build/svMultiPhysics-build/bin/svmultiphysics

solver=$c_solver

file=$1

#mpiexec -n 6 $solver $file 

$solver  $file  FILE=linear.xml


