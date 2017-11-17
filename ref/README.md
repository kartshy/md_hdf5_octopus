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


------------------------------------------------
Compiling the code:
------------------------------------------------

There is a simple Makefile that should be easily modified for most 
Unix-like environments.  There are also one or more Makefiles with 
extensions that indicate the target machine and compilers. Read the 
Makefile for further instructions.  If you generate a Makefile for 
your platform and care to share it, please send it to Paul Crozier:
pscrozi@sandia.gov . By default the code compiles with MPI support 
and can be run on one or more processors. There is also a 
Makefile.default which should NOT require a GNU Make compatible 
make. 

==Compiling: for Cray machine user platform as "cray". On Lustre this fails as some locking is still not implemented. Refer to bug
http://bugzilla.us.cray.com/show_bug.cgi?id=847639

The example works fine on th /home directories.

  module load cray-hdf5-parallel/1.8.13

  make cray
  
  Get info on all options, and targets
  
  make -f Makefile.default
  
  Build with simplified Makefile, using defaults for a CPU system
  
  make <platform>

  make clean_<platform>

  or 

  make clean


Usage:


mpiexec -n numproc miniMD_cray (MPI mode)

Example:

mpirun -np 16 ./miniMD_cray 

MiniMD understands a number of command-line options. To get the options 
for each particular variant of miniMD please use "-h" as an argument.

You will also need to provide a simple input script, which you can model
after the ones included in this directory (e.g. in.lj.miniMD). The format and
parameter description is as follows:

Sample input file contents found in "lj.in":
------------------------------------------------

Lennard-Jones input file for MD benchmark

lj             units (lj or metal)
none           data file (none or filename)       
lj             force style (lj or eam)
1.0 1.0        force parameters for LJ (epsilon, sigma)
32 32 32       size of problem
100            timesteps
0.005          timestep size 
1.44           initial temperature 
0.8442         density 
20             reneighboring every this many steps
2.5 0.30       force cutoff and neighbor skin 
100            thermo calculation every this many steps (0 = start,end)


------------------------------------------------

Sample output file contents found in "out.lj.miniMD":
------------------------------------------------

# Create System:
# Done .... 
# miniMD-Reference 1.2 (MPI+OpenMP) output ...
# Run Settings: 
	# MPI processes: 2
	# OpenMP threads: 16
	# Inputfile: in.lj.miniMD
	# Datafile: None
# Physics Settings: 
	# ForceStyle: LJ
	# Force Parameters: 1.00 1.00
	# Units: LJ
	# Atoms: 864000
	# System size: 100.78 100.78 100.78 (unit cells: 60 60 60)
	# Density: 0.844200
	# Force cutoff: 2.500000
	# Timestep size: 0.005000
# Technical Settings: 
	# Neigh cutoff: 2.800000
	# Half neighborlists: 0
	# Neighbor bins: 50 50 50
	# Neighbor frequency: 20
	# Sorting frequency: 20
	# Thermo frequency: 100
	# Ghost Newton: 1
	# Use intrinsics: 0
	# Do safe exchange: 0
	# Size of float: 8

# Starting dynamics ...
# Timestep T U P Time
0 1.440000e+00 -6.773368e+00 -5.019671e+00  0.000
100 7.310629e-01 -5.712170e+00 1.204577e+00  3.650


# Performance Summary:
# MPI_proc OMP_threads nsteps natoms t_total t_force t_neigh t_comm t_other performance perf/thread grep_string t_extra
2 16 100 864000 3.649762 2.584821 0.735003 0.145945 0.183993 23672777.021430 739774.281920 PERF_SUMMARY 0.035863

