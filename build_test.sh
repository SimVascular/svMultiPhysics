#!/bin/bash

# Exit immediately if a command fails
set -e

# Remove the "build" directory if it exists
rm -rf build

# Create a new "build" directory and navigate into it
mkdir build
cd build

# clean the build directory
#make clean

# Run cmake and make
cmake ../ 
#-DCMAKE_BUILD_TYPE=Debug -DENABLE_ARRAY_INDEX_CHECKING=ON
make -j8

# Return to the parent directory
cd ..

# Compile the svZeroDSolver
#cd svZeroDSolver
#rm -rf build
#mkdir build
#cd build
#cmake ..
#make -j8
#cd ..

# Navigate to LV_NeoHookean_passive_sv0D_capped directory
cd tests/cases/struct/LV_NeoHookean_passive_genBC_capped
#cd tests/cases/struct/104-final-sim-9-25
rm -rf GenBC.int
rm -rf InitialData
rm -rf AllData

cd genBC_svMultiPhysics
make clean
make
cd ..

# Run the test
mpirun -np 6 ../../../../build/svMultiPhysics-build/bin/svmultiphysics solver.xml > log6.txt

echo "6-proc test done"

#rm -rf GenBC.int
#rm -rf InitialData
#rm -rf AllData

#cd genBC_svMultiPhysics
#make clean
#make
#cd ..

#mpirun -np 6 ../../../../build/svMultiPhysics-build/bin/svmultiphysics solver.xml > log6.txt

#echo "6-proc test done"

# Return to the parent directory
#cd ../../../../
