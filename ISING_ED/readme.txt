This is the fortran code for the exact diagonalization(ED) of transverse field Ising model Hamiltonian.
You need to install gfortran and LAPACK to build the executable binary file.
On your UBUNTU linux box, you may install gfortran and LAPACK using "apt-get".

apt-get install gfortran
apt-get install libblas-dev liblapack-dev

compiling:

gfortran -o ed.x ed.f90 -llapack

running:

./ed.x

You can compare with the lowest energy level from ED and the ground state energy from NQS-VMC.
You should input the value of hfield and number of sites as inputs.
