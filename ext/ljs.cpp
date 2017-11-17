/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National
   Laboratories ( http://www.mantevo.org ). The primary
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov), Paul Crozier
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation;
   either version 3 of the License, or (at your option) any later
   version.

   This code is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov) or
   Christian Trott (crtrott@sandia.gov).

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"

#include "variant.h"
#include "ljs.h"
#include "atom.h"
#include "neighbor.h"
#include "integrate.h"
#include "thermo.h"
#include "comm.h"
#include "timer.h"
#include "threadData.h"
#include "string.h"
#include "openmp.h"
#include "force_eam.h"
#include "force.h"
#include "force_lj.h"
#include "hdf5.h"
#include "hdf5_hl.h"


#define MAXLINE 256

int input(In &, const char*);
void create_box(Atom &, int, int, int, double);
int create_atoms(Atom &, int, int, int, double);
void create_velocity(double, Atom &, Thermo &);
void output(In &, Atom &, Force*, Neighbor &, Comm &,
            Thermo &, Integrate &, Timer &, int);
int read_lammps_data(Atom &atom, Comm &comm, Neighbor &neighbor, Integrate &integrate, Thermo &thermo, char* file, int units);

int main(int argc, char** argv)
{
  In in;
  in.datafile = NULL;
  int me = 0;                   //local MPI rank
  int nprocs = 1;               //number of MPI ranks
  int num_threads = 1;		//number of OpenMP threads
  int num_steps = -1;           //number of timesteps (if -1 use value from lj.in)
  int system_size = -1;         //size of the system (if -1 use value from lj.in)
  int nx = -1;
  int ny = -1;
  int nz = -1;
  int check_safeexchange = 0;   //if 1 complain if atom moves further than 1 subdomain length between exchanges
  int do_safeexchange = 0;      //if 1 use safe exchange mode [allows exchange over multiple subdomains]
  int use_sse = 0;              //setting for SSE variant of miniMD only
  int screen_yaml = 0;          //print yaml output to screen also
  int yaml_output = 0;          //print yaml output
  int halfneigh = 1;            //1: use half neighborlist; 0: use full neighborlist; -1: use original miniMD version half neighborlist force
  int teams = 1;
  int device = 0;
  int neighbor_size = -1;
  char* input_file = NULL;
  int ghost_newton = 1;
  int sort = -1;
  int ntypes = 4;

  for(int i = 0; i < argc; i++) {
    if((strcmp(argv[i], "-i") == 0) || (strcmp(argv[i], "--input_file") == 0)) {
      input_file = argv[++i];
      continue;
    }
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int error = 0;

  if(input_file == NULL)
    error = input(in, "in.lj.miniMD");
  else
    error = input(in, input_file);

  if(error) {
    MPI_Finalize();
    exit(0);
  }

  srand(5413);

  for(int i = 0; i < argc; i++) {
    if((strcmp(argv[i], "-t") == 0) || (strcmp(argv[i], "--num_threads") == 0)) {
      num_threads = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "--teams") == 0)) {
      teams = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-n") == 0) || (strcmp(argv[i], "--nsteps") == 0))  {
      num_steps = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-s") == 0) || (strcmp(argv[i], "--size") == 0)) {
      system_size = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-nx") == 0)) {
      nx = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-ny") == 0)) {
      ny = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-nz") == 0)) {
      nz = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "--ntypes") == 0)) {
      ntypes = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-b") == 0) || (strcmp(argv[i], "--neigh_bins") == 0))  {
      neighbor_size = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "--half_neigh") == 0))  {
      halfneigh = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-sse") == 0))  {
      use_sse = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "--check_exchange") == 0))  {
      check_safeexchange = 1;
      continue;
    }

    if((strcmp(argv[i], "--sort") == 0))  {
      sort = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-o") == 0) || (strcmp(argv[i], "--yaml_output") == 0))  {
      yaml_output = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "--yaml_screen") == 0))  {
      screen_yaml = 1;
      continue;
    }

    if((strcmp(argv[i], "-f") == 0) || (strcmp(argv[i], "--data_file") == 0)) {
      if(in.datafile == NULL) in.datafile = new char[1000];

      strcpy(in.datafile, argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-u") == 0) || (strcmp(argv[i], "--units") == 0)) {
      in.units = strcmp(argv[++i], "metal") == 0 ? 1 : 0;
      continue;
    }

    if((strcmp(argv[i], "-p") == 0) || (strcmp(argv[i], "--force") == 0)) {
      in.forcetype = strcmp(argv[++i], "eam") == 0 ? FORCEEAM : FORCELJ;
      continue;
    }

    if((strcmp(argv[i], "-gn") == 0) || (strcmp(argv[i], "--ghost_newton") == 0)) {
      ghost_newton = atoi(argv[++i]);
      continue;
    }

    if((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
      printf("\n-----------------------------------------------------------------------------------------------------------\n");
      printf("-------------" VARIANT_STRING "--------------------\n");
      printf("-------------------------------------------------------------------------------------------------------------\n\n");

      printf("miniAnalyser is a simple, parallel molecular dynamics (MD) code,\n"
             "which is part of the Mantevo project at Sandia National\n"
             "Laboratories ( http://www.mantevo.org ).\n"
             "The original authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) ,\n"
             "Paul Crozier (pscrozi@sandia.gov) with current\n"
             "versions written by Christian Trott (crtrott@sandia.gov).\n\n");
      printf("Commandline Options:\n");
      printf("\n  Execution configuration:\n");
      printf("\t--teams <nteams>:             set number of thread-teams used per MPI rank (default 1)\n");
      printf("\t-t / --num_threads <threads>: set number of threads per thread-team (default 1)\n");
      printf("\t--half_neigh <int>:           use half neighborlists (default 1)\n"
             "\t                                0: full neighborlist\n"
             "\t                                1: half neighborlist\n"
             "\t                               -1: original miniMD half neighborlist force (not OpenMP safe)\n");
      printf("\t-d / --device <int>:          choose device to use (only applicable for GPU execution)\n");
      printf("\t-dm / --device_map:           map devices to MPI ranks\n");
      printf("\t-ng / --num_gpus <int>:       give number of GPUs per Node (used in conjuction with -dm\n"
             "\t                              to determine device id: 'id=mpi_rank%%ng' (default 2)\n");
      printf("\t--skip_gpu <int>:             skip the specified gpu when assigning devices to MPI ranks\n"
             "\t                              used in conjunction with -dm (but must come first in arg list)\n");
      printf("\t-sse <sse_version>:           use explicit sse intrinsics (use miniMD-SSE variant)\n");
      printf("\t-gn / --ghost_newton <int>:   set usage of newtons third law for ghost atoms\n"
             "\t                                (only applicable with half neighborlists)\n");
      printf("\n  Simulation setup:\n");
      printf("\t-i / --input_file <string>:   set input file to be used (default: in.lj.miniMD)\n");
      printf("\t--ntypes <int>:               set number of atom types for simulation (default: 4)\n");
      printf("\t-n / --nsteps <int>:          set number of timesteps for simulation\n");
      printf("\t-s / --size <int>:            set linear dimension of systembox\n");
      printf("\t-nx/-ny/-nz <int>:            set linear dimension of systembox in x/y/z direction\n");
      printf("\t-b / --neigh_bins <int>:      set linear dimension of neighbor bin grid\n");
      printf("\t-u / --units <string>:        set units (lj or metal), see LAMMPS documentation\n");
      printf("\t-p / --force <string>:        set interaction model (lj or eam)\n");
      printf("\t-f / --data_file <string>:    read configuration from LAMMPS data file\n");

      printf("\n  Miscelaneous:\n");
      printf("\t--check_exchange:             check whether atoms moved further than subdomain width\n");
      printf("\t--safe_exchange:              perform exchange communication with all MPI processes\n"
             "\t                                within rcut_neighbor (outer force cutoff)\n");
      printf("\t--sort <n>:                   resort atoms (simple bins) every <n> steps (default: use reneigh frequency; never=0)");
      printf("\t-o / --yaml_output <int>:     level of yaml output (default 1)\n");
      printf("\t--yaml_screen:                write yaml output also to screen\n");
      printf("\t-h / --help:                  display this help message\n\n");
      printf("---------------------------------------------------------\n\n");

      exit(0);
    }
  }


  Atom atom(ntypes);
  Neighbor neighbor(ntypes);
  Integrate integrate;
  Thermo thermo;
  Comm comm;
  Timer timer;
  ThreadData threads;

  Force* force;

  if(in.forcetype == FORCEEAM) {
    force = (Force*) new ForceEAM(ntypes);

    if(ghost_newton == 1) {
      if(me == 0)
        printf("# EAM currently requires '--ghost_newton 0'; Changing setting now.\n");

      ghost_newton = 0;
    }
  }

  if(in.forcetype == FORCELJ) force = (Force*) new ForceLJ(ntypes);

  threads.mpi_me = me;
  threads.mpi_num_threads = nprocs;
  threads.omp_me = 0;
  threads.omp_num_threads = num_threads;

  atom.threads = &threads;
  comm.threads = &threads;
  force->threads = &threads;
  integrate.threads = &threads;
  neighbor.threads = &threads;
  thermo.threads = &threads;

  if(in.forcetype == FORCELJ) {
    for(int i=0; i<ntypes*ntypes; i++) {
      force->epsilon[i] = in.epsilon;
      force->sigma[i] = in.sigma;
      force->sigma6[i] = in.sigma*in.sigma*in.sigma*in.sigma*in.sigma*in.sigma;
    }
  }

  neighbor.ghost_newton = ghost_newton;

  omp_set_num_threads(num_threads);

  neighbor.timer = &timer;
  force->timer = &timer;
  comm.check_safeexchange = check_safeexchange;
  comm.do_safeexchange = do_safeexchange;
  force->use_sse = use_sse;
  neighbor.halfneigh = halfneigh;

  if(halfneigh < 0) force->use_oldcompute = 1;

  if(use_sse) {
#ifdef VARIANT_REFERENCE

    if(me == 0) printf("ERROR: Trying to run with -sse with miniMD reference version. Use SSE variant instead. Exiting.\n");

    MPI_Finalize();
    exit(0);
#endif
  }

  if(num_steps > 0) in.ntimes = num_steps;

  if(system_size > 0) {
    in.nx = system_size;
    in.ny = system_size;
    in.nz = system_size;
  }

  if(nx > 0) {
    in.nx = nx;
    if(ny > 0)
      in.ny = ny;
    else if(system_size < 0)
      in.ny = nx;

    if(nz > 0)
      in.nz = nz;
    else if(system_size < 0)
      in.nz = nx;
  }

  if(neighbor_size > 0) {
    neighbor.nbinx = neighbor_size;
    neighbor.nbiny = neighbor_size;
    neighbor.nbinz = neighbor_size;
  }

  if(neighbor_size < 0 && in.datafile == NULL) {
    MMD_float neighscale = 5.0 / 6.0;
    neighbor.nbinx = neighscale * in.nx;
    neighbor.nbiny = neighscale * in.ny;
    neighbor.nbinz = neighscale * in.nz;
  }

  if(neighbor_size < 0 && in.datafile)
    neighbor.nbinx = -1;

  if(neighbor.nbinx == 0) neighbor.nbinx = 1;

  if(neighbor.nbiny == 0) neighbor.nbiny = 1;

  if(neighbor.nbinz == 0) neighbor.nbinz = 1;

  integrate.ntimes = in.ntimes;
  integrate.dt = in.dt;
  integrate.sort_every = sort>0?sort:(sort<0?in.neigh_every:0);
  neighbor.every = in.neigh_every;
  neighbor.cutneigh = in.neigh_cut;
  force->cutforce = in.force_cut;
  thermo.nstat = in.thermo_nstat;


  if(me == 0)
    printf("# Create System:\n");

  if(in.datafile) {
    read_lammps_data(atom, comm, neighbor, integrate, thermo, in.datafile, in.units);
    MMD_float volume = atom.box.xprd * atom.box.yprd * atom.box.zprd;
    in.rho = 1.0 * atom.natoms / volume;
    force->setup();

    if(in.forcetype == FORCEEAM) atom.mass = force->mass;
  } else {
    create_box(atom, in.nx, in.ny, in.nz, in.rho);

    comm.setup(neighbor.cutneigh, atom);

    neighbor.setup(atom);

    integrate.setup();

    force->setup();

    if(in.forcetype == FORCEEAM) atom.mass = force->mass;

    create_atoms(atom, in.nx, in.ny, in.nz, in.rho);
    thermo.setup(in.rho, integrate, atom, in.units);

    create_velocity(in.t_request, atom, thermo);

  }

  if(me == 0)
    printf("# Done .... \n");

  if(me == 0) {
    fprintf(stdout, "# " VARIANT_STRING " output ...\n");
    fprintf(stdout, "# Run Settings: \n");
    fprintf(stdout, "\t# MPI processes: %i\n", neighbor.threads->mpi_num_threads);
    fprintf(stdout, "\t# OpenMP threads: %i\n", neighbor.threads->omp_num_threads);
    fprintf(stdout, "\t# Inputfile: %s\n", input_file == 0 ? "in.lj.miniMD" : input_file);
    fprintf(stdout, "\t# Datafile: %s\n", in.datafile ? in.datafile : "None");
    fprintf(stdout, "# Physics Settings: \n");
    fprintf(stdout, "\t# ForceStyle: %s\n", in.forcetype == FORCELJ ? "LJ" : "EAM");
    fprintf(stdout, "\t# Force Parameters: %2.2lf %2.2lf\n",in.epsilon,in.sigma);
    fprintf(stdout, "\t# Units: %s\n", in.units == 0 ? "LJ" : "METAL");
    fprintf(stdout, "\t# Atoms: %i\n", atom.natoms);
    fprintf(stdout, "\t# Atom types: %i\n", atom.ntypes);
    fprintf(stdout, "\t# System size: %2.2lf %2.2lf %2.2lf (unit cells: %i %i %i)\n", atom.box.xprd, atom.box.yprd, atom.box.zprd, in.nx, in.ny, in.nz);
    fprintf(stdout, "\t# Density: %lf\n", in.rho);
    fprintf(stdout, "\t# Force cutoff: %lf\n", force->cutforce);
    fprintf(stdout, "\t# Timestep size: %lf\n", integrate.dt);
    fprintf(stdout, "# Technical Settings: \n");
    fprintf(stdout, "\t# Neigh cutoff: %lf\n", neighbor.cutneigh);
    fprintf(stdout, "\t# Half neighborlists: %i\n", neighbor.halfneigh);
    fprintf(stdout, "\t# Neighbor bins: %i %i %i\n", neighbor.nbinx, neighbor.nbiny, neighbor.nbinz);
    fprintf(stdout, "\t# Neighbor frequency: %i\n", neighbor.every);
    fprintf(stdout, "\t# Sorting frequency: %i\n", integrate.sort_every);
    fprintf(stdout, "\t# Thermo frequency: %i\n", thermo.nstat);
    fprintf(stdout, "\t# Ghost Newton: %i\n", ghost_newton);
    fprintf(stdout, "\t# Use intrinsics: %i\n", force->use_sse);
    fprintf(stdout, "\t# Do safe exchange: %i\n", comm.do_safeexchange);
    fprintf(stdout, "\t# Size of float: %i\n\n", (int) sizeof(MMD_float));
  }

  comm.exchange(atom);
  if(sort>0)
    atom.sort(neighbor);
  comm.borders(atom);

  int natoms;
  MPI_Allreduce(&atom.nlocal, &natoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
  //printf("%d : %d %d\n", me,  atom.nlocal, natoms);

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,  MPI_COMM_WORLD, MPI_INFO_NULL);

  hid_t file = H5Fopen ("../file.h5", H5F_ACC_RDONLY, plist_id);
  hid_t dataset_x = H5Dopen (file, "/atom.x", H5P_DEFAULT);
  hid_t dataset_v = H5Dopen (file, "/atom.v", H5P_DEFAULT);
  hid_t dataset_f = H5Dopen (file, "/atom.f", H5P_DEFAULT);
  hid_t dataspace = H5Dget_space (dataset_x);    /* dataspace handle */
  int rank      = H5Sget_simple_extent_ndims (dataspace);
  hsize_t     dims_out[rank]; 
  herr_t status  = H5Sget_simple_extent_dims (dataspace, dims_out, NULL);

  hsize_t offset = me * atom.nlocal*3  ;
  hsize_t count = atom.nlocal*3 ;

  //hsize_t offset = me * 16*3  ;
  //hsize_t count = 16*3 ;
  //printf("\nRank: %d\nDimensions: %lu  %lu %lu\n", rank,
  //         (unsigned long)(dims_out[0] ), (unsigned long)(count), (unsigned long)(offset));
  status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, &offset, NULL, 
                                  &count, NULL);

  hid_t memspace = H5Screate_simple (1, &count , NULL);  

  status = H5Dread (dataset_x, H5T_NATIVE_DOUBLE, memspace, dataspace,
                      H5P_DEFAULT, atom.x);

  status = H5Dread (dataset_v, H5T_NATIVE_DOUBLE, memspace, dataspace,
                      H5P_DEFAULT, atom.v);

  status = H5Dread (dataset_f, H5T_NATIVE_DOUBLE, memspace, dataspace,
                      H5P_DEFAULT, atom.f);

  double density;
  H5LTget_attribute_double(file, "/atoms", "density", &density );

  H5Pclose(plist_id);
  H5Dclose(dataset_f);
  H5Dclose(dataset_v);
  H5Dclose(dataset_x);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Fclose(file);
  /*printf("X %d : %f %f\n ", me,  atom.x[me+1], atom.x[me+2]);
  printf("V %d : %f %f\n ", me,  atom.v[me+1], atom.v[me+2]);
  printf("F %d : %f %f\n ", me,  atom.f[me+1], atom.f[me+2]);
  printf("density %f\n", density);
  */
 //This analysis is not science and does not mean anything

 double norm = 0.0;
 for(int i = 0; i < atom.nlocal; i++) {

     norm += atom.v[i*PAD  + 0] * atom.v[i*PAD  + 0] * atom.f[i*PAD  + 0] * atom.x[i*PAD  + 0]*  atom.x[i*PAD  + 0] ;
     norm += atom.v[i*PAD  + 1] * atom.v[i*PAD  + 1] * atom.f[i*PAD  + 1] * atom.x[i*PAD  + 0]*  atom.x[i*PAD  + 1] ;
     norm += atom.v[i*PAD  + 2] * atom.v[i*PAD  + 2] * atom.f[i*PAD  + 2] * atom.x[i*PAD  + 0]*  atom.x[i*PAD  + 2] ;
 }

 norm = norm/density;
 double totalnorm;
 MPI_Allreduce( &norm, &totalnorm, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD );

 if(me==0)
   printf(" Norm %f\n", totalnorm);

  delete force;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

