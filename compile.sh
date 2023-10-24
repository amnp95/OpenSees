cd /home/amnp95/Projects/OpenSees
# rm -rf build
# mkdir build
cd build
rm ./build/CMakeFiles/OPS_Analysis.dir/SRC/analysis/integrator/Newmark.cpp.o
rm  /home/amnp95/Projects/OpenSees/build/CMakeFiles/OPS_Element.dir/SRC/element/PML/PML3D.cpp.o
rm  /home/amnp95/Projects/OpenSees/build/CMakeFiles/OPS_Element.dir/SRC/element/PML/PML3DVISCOUS.cpp.o
rm /home/amnp95/Projects/OpenSees/build/CMakeFiles/OPS_Domain.dir/SRC/domain/pattern/drm/H5DRM.cpp.o
# $HOME/.local/bin/conan install .. --build missing
cmake .. -DMUMPS_DIR=$PWD/../../mumps/build -DOPENMPI=TRUE -DSCALAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/libmkl_blacs_openmpi_lp64.so;/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so.2.1"
cmake --build . --config Release --target OpenSees --parallel 4
# # cmake --build . --config Release --target OpenSeesMP
cmake --build . --config Release --target OpenSeesSP --parallel 4
# cd /home/amnp95/Projects/Github/OpenSeesProjects/Shakermaker/SP
# mpirun -np 6 OpenSeesSP PML3Dbox.tcl YES

