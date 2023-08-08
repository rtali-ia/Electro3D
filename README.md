# Electro3D

Instructions: \

1. First build the exe under /EllipsoidHole3D \
   a. cd EllipsoidHole3D \
   b. mkdir build-anvil \
   c. cd build-anvil \
   d. cmake .. -DCMAKE_BUILD_TYPE=Release -Dtalyfem_DIR=/anvil/scratch/x-rtali/packages/taly_fem/build -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
   e. make -j8 \

2. In config.txt make sure the location of msh file is set to ./mesh/Ellipsoid3D.msh \

3. Now run the job on Anvil \
   a. cd Electro3D \
   b. In "j1.job" change {EXECUTABLE} to point to ../Ellipsoid3D/build-anvil/capacitance3d \
   c. sbatch j1.job \
