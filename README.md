------------------------------------------------
Description:
------------------------------------------------
The Octopus example uses parallel HDF5 file format for a MD based model and analysis.

miniMD is a parallel molecular dynamics (MD) simulation package written
in C++ and intended for use on parallel supercomputers and new 
architechtures for testing purposes. The software package is meant to  
be simple, lightweight, and easily adapted to new hardware. It is 
designed following many of the same algorithm concepts as our LAMMPS 
(http://lammps.sandia.gov) parallel MD code, but is much simpler.


This simple code is a self-contained piece of C++ software 
that performs parallel molecular dynamics simulation 


The sub-directories contain different variants of miniMD:

ref:          model (S) supports MPI+OpenMP hybrid mode.
ext:          analysis (A) supports MPI + OpenMP hybrid code

A sample worfflow is shown in the md_hdf5.sh file.

------------------------------------------------
Compiling the code:
------------------------------------------------
module load cray-hdf5-parallel/1.8.13

cd ref
make cray

cd ext
make cray

Refer to the README files in each folder for more details.

Requiremnts for Octopus:

TODO

