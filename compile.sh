cd /home/amin/Desktop/UWOpenSees/OpenSees
rm -rf build
mkdir build
cd build
rm /home/amin/Desktop/UWOpenSees/OpenSees/build/CMakeFiles/OPS_Analysis.dir/SRC/analysis/integrator/Newmark.cpp.o
rm  /home/amin/Desktop/UWOpenSees/OpenSees/build/CMakeFiles/OPS_Element.dir/SRC/element/PML/PML3D.cpp.o
$HOME/.local/bin/conan install .. --build missing
cmake .. -DMUMPS_DIR=$PWD/../../mumps/build -DOPENMPI=TRUE -DSCALAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/libmkl_blacs_openmpi_lp64.so;/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so.2.1"
cmake --build . --config Release --target OpenSees --parallel 4
# cmake --build . --config Release --target OpenSeesMP
cmake --build . --config Release --target OpenSeesSP
cd /home/amin/Projects/Github/OpenSeesProjects/Shakermaker/SP/
mpirun -np 2 /home/amin/Desktop/UWOpenSees/OpenSees/build/bin/OpenSeesSP PML3Dbox.tcl YES

cd /home/amin/Desktop/UWOpenSees/OpenSees
