/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)
------------------------------------------------------------------------- */

#include "fix_ttm.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "modify.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "compute.h"
#include "iostream"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixTTM::FixTTM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  random(NULL), fp(NULL), fpr(NULL), nsum(NULL), nsum_all(NULL),
  T_initial_set(NULL), gfactor1(NULL), gfactor2(NULL), ratio(NULL),
  flangevin(NULL), T_electron(NULL), T_electron_old(NULL), sum_vsq(NULL),
  sum_mass_vsq(NULL), sum_vsq_all(NULL), sum_mass_vsq_all(NULL),
  net_energy_transfer(NULL), net_energy_transfer_all(NULL), u_node(NULL), v_node(NULL), w_node(NULL), nvel(NULL),
  u_node_all(NULL), v_node_all(NULL), w_node_all(NULL), nvel_all(NULL)
{
  if (narg < 15) error->all(FLERR,"Illegal fix ttm command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 1;
  nevery = 1;
  restart_peratom = 1;
  restart_global = 1;

  id_temp = NULL;
  temperature = NULL;
  seed = force->inumeric(FLERR,arg[3]);
  electronic_specific_heat = force->numeric(FLERR,arg[4]);
  electronic_density = force->numeric(FLERR,arg[5]);
  electronic_thermal_conductivity = force->numeric(FLERR,arg[6]);
  gamma_p = force->numeric(FLERR,arg[7]);
  gamma_s = force->numeric(FLERR,arg[8]);
  v_0 = force->numeric(FLERR,arg[9]);
  nxnodes = force->inumeric(FLERR,arg[10]);
  nynodes = force->inumeric(FLERR,arg[11]);
  nznodes = force->inumeric(FLERR,arg[12]);
  nodes_xyz[0] = nxnodes;
  nodes_xyz[1] = nynodes;
  nodes_xyz[2] = nznodes;

  fpr = fopen(arg[13],"r");
  if (fpr == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open file %s",arg[13]);
    error->one(FLERR,str);
  }

  nfileevery = force->inumeric(FLERR,arg[14]);

  if (nfileevery) {
    if (narg != 17) error->all(FLERR,"Illegal fix ttm command");
    MPI_Comm_rank(world,&me);
    if (me == 0) {
      fp = fopen(arg[15],"w");
      if (fp == NULL) {
        char str[128];
        snprintf(str,128,"Cannot open fix ttm file %s",arg[15]);
        error->one(FLERR,str);
      }
    }
  }

  // error check

  if (seed <= 0)
    error->all(FLERR,"Invalid random number seed in fix ttm command");
  if (electronic_specific_heat <= 0.0)
    error->all(FLERR,"Fix ttm electronic_specific_heat must be > 0.0");
  if (electronic_density <= 0.0)
    error->all(FLERR,"Fix ttm electronic_density must be > 0.0");
  if (electronic_thermal_conductivity < 0.0)
    error->all(FLERR,"Fix ttm electronic_thermal_conductivity must be >= 0.0");
  if (gamma_p <= 0.0) error->all(FLERR,"Fix ttm gamma_p must be > 0.0");
  if (gamma_s < 0.0) error->all(FLERR,"Fix ttm gamma_s must be >= 0.0");
  if (v_0 < 0.0) error->all(FLERR,"Fix ttm v_0 must be >= 0.0");
  if (nxnodes <= 0 || nynodes <= 0 || nznodes <= 0)
    error->all(FLERR,"Fix ttm number of nodes must be > 0");

  for (whichfix = 0; whichfix < modify->nfix; whichfix++) {

     if (strcmp(arg[narg-1],modify->fix[whichfix]->id) == 0){

         break;

     }  

     if (whichfix == modify->nfix){

         error->universe_all(FLERR,"ttm fix ID is not defined");

     }
  } 

  v_0_sq = v_0*v_0;

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];

  // allocate 3d grid variables

  total_nnodes = nxnodes*nynodes*nznodes;

  memory->create(nsum,nxnodes,nynodes,nznodes,"ttm:nsum");
  memory->create(nsum_all,nxnodes,nynodes,nznodes,"ttm:nsum_all");
  memory->create(T_initial_set,nxnodes,nynodes,nznodes,"ttm:T_initial_set");
  memory->create(sum_vsq,nxnodes,nynodes,nznodes,"ttm:sum_vsq");
  memory->create(sum_mass_vsq,nxnodes,nynodes,nznodes,"ttm:sum_mass_vsq");
  memory->create(sum_vsq_all,nxnodes,nynodes,nznodes,"ttm:sum_vsq_all");
  memory->create(sum_mass_vsq_all,nxnodes,nynodes,nznodes,
                 "ttm:sum_mass_vsq_all");
  memory->create(T_electron_old,nxnodes,nynodes,nznodes,"ttm:T_electron_old");
  memory->create(T_electron,nxnodes,nynodes,nznodes,"ttm:T_electron");
  memory->create(u_node,nxnodes,nynodes,nznodes,"ttm:u_node");
  memory->create(v_node,nxnodes,nynodes,nznodes,"ttm:v_node");
  memory->create(w_node,nxnodes,nynodes,nznodes,"ttm:w_node");
  memory->create(nvel,nxnodes,nynodes,nznodes,"ttm:nvel"); 
  memory->create(u_node_all,nxnodes,nynodes,nznodes,"ttm:u_node_all");
  memory->create(v_node_all,nxnodes,nynodes,nznodes,"ttm:v_node_all");
  memory->create(w_node_all,nxnodes,nynodes,nznodes,"ttm:w_node_all");
  memory->create(nvel_all,nxnodes,nynodes,nznodes,"ttm:nvel_all");
  memory->create(net_energy_transfer,nxnodes,nynodes,nznodes,
                 "TTM:net_energy_transfer");
  memory->create(net_energy_transfer_all,nxnodes,nynodes,nznodes,
                 "TTM:net_energy_transfer_all");

  flangevin = NULL;
  grow_arrays(atom->nmax);

  // zero out the flangevin array

  for (int i = 0; i < atom->nmax; i++) {
    flangevin[i][0] = 0;
    flangevin[i][1] = 0;
    flangevin[i][2] = 0;
  }

  atom->add_callback(0);
  atom->add_callback(1);

  // set initial electron temperatures from user input file

  if (me == 0) read_initial_electron_temperatures();
  MPI_Bcast(&T_electron[0][0][0],total_nnodes,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

FixTTM::~FixTTM()
{
  if (nfileevery && me == 0) fclose(fp);

  delete random;

  delete [] gfactor1;
  delete [] gfactor2;

  memory->destroy(nsum);
  memory->destroy(nsum_all);
  memory->destroy(T_initial_set);
  memory->destroy(sum_vsq);
  memory->destroy(sum_mass_vsq);
  memory->destroy(sum_vsq_all);
  memory->destroy(sum_mass_vsq_all);
  memory->destroy(T_electron_old);
  memory->destroy(T_electron);
  memory->destroy(u_node);
  memory->destroy(v_node);
  memory->destroy(w_node);
  memory->destroy(nvel);
  memory->destroy(u_node_all);
  memory->destroy(v_node_all);
  memory->destroy(w_node_all);
  memory->destroy(nvel_all);
  memory->destroy(flangevin);
  memory->destroy(net_energy_transfer);
  memory->destroy(net_energy_transfer_all);
}

/* ---------------------------------------------------------------------- */

int FixTTM::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTTM::init()
{
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix ttm with 2d simulation");
  if (domain->nonperiodic != 0)
    error->all(FLERR,"Cannot use non-periodic boundares with fix ttm");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix ttm with triclinic box");

  // set force prefactors

  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor1[i] = - gamma_p / force->ftm2v;
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
  }

  for (int ixnode = 0; ixnode < nxnodes; ixnode++){
    for (int iynode = 0; iynode < nynodes; iynode++){
      for (int iznode = 0; iznode < nznodes; iznode++){
         net_energy_transfer_all[ixnode][iynode][iznode] = 0;
         u_node_all[ixnode][iynode][iznode] = 0;
         v_node_all[ixnode][iynode][iznode] = 0;
         w_node_all[ixnode][iynode][iznode] = 0;
         nvel_all[ixnode][iynode][iznode] = 0;
      }
    }
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixTTM::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force_setup(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa_setup(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double gamma1,gamma2;
//  compute_target();
  
  if (temperature!=0){
     temperature->compute_scalar();
  }

  // apply damping and thermostat to all atoms in fix group

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      int ixnode = static_cast<int>(xscale*nxnodes);
      int iynode = static_cast<int>(yscale*nynodes);
      int iznode = static_cast<int>(zscale*nznodes);
      while (ixnode > nxnodes-1) ixnode -= nxnodes;
      while (iynode > nynodes-1) iynode -= nynodes;
      while (iznode > nznodes-1) iznode -= nznodes;
      while (ixnode < 0) ixnode += nxnodes;
      while (iynode < 0) iynode += nynodes;
      while (iznode < 0) iznode += nznodes;

      if (T_electron[ixnode][iynode][iznode] < 0)
        error->all(FLERR,"Electronic temperature dropped below zero");
      
      if (temperature!=0){
         temperature->remove_bias(i,v[i]);
      } 

      double tsqrt = sqrt(T_electron[ixnode][iynode][iznode]);

      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > v_0_sq) gamma1 *= (gamma_p + gamma_s)/gamma_p;
      gamma2 = gfactor2[type[i]] * tsqrt;

//      temperature->remove_bias(i,v[i]);
      flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
      flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
      flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);

      if (temperature!=0){
         temperature->restore_bias(i,v[i]);
      }

      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];

    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_setup(int /*vflag*/)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // apply langevin forces that have been stored from previous run

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_respa_setup(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force_setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::reset_dt()
{
  for (int i = 1; i <= atom->ntypes; i++)
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
}

/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only called by proc 0
------------------------------------------------------------------------- */

void FixTTM::read_initial_electron_temperatures()
{
  char line[MAXLINE];

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_initial_set[ixnode][iynode][iznode] = 0;

  // read initial electron temperature values from file

  int ixnode,iynode,iznode;
  double T_tmp;
  while (1) {
    if (fgets(line,MAXLINE,fpr) == NULL) break;
    sscanf(line,"%d %d %d %lg",&ixnode,&iynode,&iznode,&T_tmp);
    if (T_tmp < 0.0)
      error->one(FLERR,"Fix ttm electron temperatures must be > 0.0");
    T_electron[ixnode][iynode][iznode] = T_tmp;
    T_initial_set[ixnode][iynode][iznode] = 1;
  }

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        if (T_initial_set[ixnode][iynode][iznode] == 0)
          error->one(FLERR,"Initial temperatures not all set in fix ttm");

  // close file

  fclose(fpr);
}

/* ---------------------------------------------------------------------- */

void FixTTM::end_of_step()
{
  if (temperature!=0){
     temperature->compute_scalar();
  }

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
//  std::cout <<"\n\nnumber: " << nlocal;
  for (int ixnode = 0; ixnode < nxnodes; ixnode++){
    for (int iynode = 0; iynode < nynodes; iynode++){
      for (int iznode = 0; iznode < nznodes; iznode++){

        net_energy_transfer[ixnode][iynode][iznode] = 0;
        u_node[ixnode][iynode][iznode] = 0;
        v_node[ixnode][iynode][iznode] = 0;
        w_node[ixnode][iynode][iznode] = 0;
        nvel[ixnode][iynode][iznode] = 0;

       }
     }
   }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      int ixnode = static_cast<int>(xscale*nxnodes);
      int iynode = static_cast<int>(yscale*nynodes);
      int iznode = static_cast<int>(zscale*nznodes);
      while (ixnode > nxnodes-1) ixnode -= nxnodes;
      while (iynode > nynodes-1) iynode -= nynodes;
      while (iznode > nznodes-1) iznode -= nznodes;
      while (ixnode < 0) ixnode += nxnodes;
      while (iynode < 0) iynode += nynodes;
      while (iznode < 0) iznode += nznodes;

      if (temperature!=0){
         temperature->remove_bias(i,v[i]);
      }

      u_node[ixnode][iynode][iznode] += v[i][0];
      v_node[ixnode][iynode][iznode] += v[i][1];
      w_node[ixnode][iynode][iznode] += v[i][2];
      nvel[ixnode][iynode][iznode] += 1;
      net_energy_transfer[ixnode][iynode][iznode] +=
        (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
         flangevin[i][2]*v[i][2]);
    }

  MPI_Allreduce(&net_energy_transfer[0][0][0],&net_energy_transfer_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&u_node[0][0][0],&u_node_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&v_node[0][0][0],&v_node_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&w_node[0][0][0],&w_node_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&nvel[0][0][0],&nvel_all[0][0][0],total_nnodes,MPI_DOUBLE,MPI_SUM,world);
  MPI_Comm_rank(world,&me);

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;
  double stability_criterion = 1.0 -
    2.0*inner_dt/(electronic_specific_heat*electronic_density) *
    (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
  if (stability_criterion < 0.0) {
    inner_dt = 0.5*(electronic_specific_heat*electronic_density) /
      (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm",0);
  }

  for (int ith_inner_timestep = 0; ith_inner_timestep < num_inner_timesteps;
       ith_inner_timestep++) {

    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++)
          T_electron_old[ixnode][iynode][iznode] =
            T_electron[ixnode][iynode][iznode];

    // compute new electron T profile

    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          int right_xnode = ixnode + 1;
          int right_ynode = iynode + 1;
          int right_znode = iznode + 1;
          if (right_xnode == nxnodes) right_xnode = 0;
          if (right_ynode == nynodes) right_ynode = 0;
          if (right_znode == nznodes) right_znode = 0;
          int left_xnode = ixnode - 1;
          int left_ynode = iynode - 1;
          int left_znode = iznode - 1;
          if (left_xnode == -1) left_xnode = nxnodes - 1;
          if (left_ynode == -1) left_ynode = nynodes - 1;
          if (left_znode == -1) left_znode = nznodes - 1;


// Boundary conditions for test case

          if(iynode == 0){

		T_electron[ixnode][iynode][iznode] = 0.1;

          }

	  else if (iynode == nynodes - 1){

		T_electron[ixnode][iynode][iznode] = 0.8;

          }

/////////////////////////////////////

	  else {

          T_electron[ixnode][iynode][iznode] =
            T_electron_old[ixnode][iynode][iznode] +
            inner_dt/(electronic_specific_heat*electronic_density) *
            (electronic_thermal_conductivity *
             ((T_electron_old[right_xnode][iynode][iznode] +
               T_electron_old[left_xnode][iynode][iznode] -
               2*T_electron_old[ixnode][iynode][iznode])/dx/dx +
              (T_electron_old[ixnode][right_ynode][iznode] +
               T_electron_old[ixnode][left_ynode][iznode] -
               2*T_electron_old[ixnode][iynode][iznode])/dy/dy +
              (T_electron_old[ixnode][iynode][right_znode] +
               T_electron_old[ixnode][iynode][left_znode] -
               2*T_electron_old[ixnode][iynode][iznode])/dz/dz) -
              (net_energy_transfer_all[ixnode][iynode][iznode])/del_vol) -
              u_node_all[ixnode][iynode][iznode]*T_electron_old[right_xnode][iynode][iznode]*(0.5*(inner_dt/dx)/nvel_all[ixnode][iynode][iznode]) +
              u_node_all[ixnode][iynode][iznode]*T_electron_old[left_xnode][iynode][iznode]*(0.5*(inner_dt/dx)/nvel_all[ixnode][iynode][iznode]) -
              v_node_all[ixnode][iynode][iznode]*T_electron_old[ixnode][right_ynode][iznode]*(0.5*(inner_dt/dy)/nvel_all[ixnode][iynode][iznode]) +
              v_node_all[ixnode][iynode][iznode]*T_electron_old[ixnode][left_ynode][iznode]*(0.5*(inner_dt/dy)/nvel_all[ixnode][iynode][iznode]) -
              w_node_all[ixnode][iynode][iznode]*T_electron_old[ixnode][iynode][right_znode]*(0.5*(inner_dt/dz)/nvel_all[ixnode][iynode][iznode]) +
              w_node_all[ixnode][iynode][iznode]*T_electron_old[ixnode][iynode][left_znode]*(0.5*(inner_dt/dz)/nvel_all[ixnode][iynode][iznode]);

            }
        }
  }

  // output nodal temperatures for current timestep

  if ((nfileevery) && !(update->ntimestep % nfileevery)) {

    // compute atomic Ta for each grid point

    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
        for (int iznode = 0; iznode < nznodes; iznode++) {
          nsum[ixnode][iynode][iznode] = 0;
          nsum_all[ixnode][iynode][iznode] = 0;
          sum_vsq[ixnode][iynode][iznode] = 0.0;
          sum_mass_vsq[ixnode][iynode][iznode] = 0.0;
          sum_vsq_all[ixnode][iynode][iznode] = 0.0;
          sum_mass_vsq_all[ixnode][iynode][iznode] = 0.0;
        }

    double massone;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
        double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
        double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
        int ixnode = static_cast<int>(xscale*nxnodes);
        int iynode = static_cast<int>(yscale*nynodes);
        int iznode = static_cast<int>(zscale*nznodes);
        while (ixnode > nxnodes-1) ixnode -= nxnodes;
        while (iynode > nynodes-1) iynode -= nynodes;
        while (iznode > nznodes-1) iznode -= nznodes;
        while (ixnode < 0) ixnode += nxnodes;
        while (iynode < 0) iynode += nynodes;
        while (iznode < 0) iznode += nznodes;
        double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        nsum[ixnode][iynode][iznode] += 1;
        sum_vsq[ixnode][iynode][iznode] += vsq;
        sum_mass_vsq[ixnode][iynode][iznode] += massone*vsq;
        if (temperature!=0){
           temperature->restore_bias(i,v[i]);
        }
      }

    MPI_Allreduce(&nsum[0][0][0],&nsum_all[0][0][0],total_nnodes,
                  MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&sum_vsq[0][0][0],&sum_vsq_all[0][0][0],total_nnodes,
                  MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&sum_mass_vsq[0][0][0],&sum_mass_vsq_all[0][0][0],
                  total_nnodes,MPI_DOUBLE,MPI_SUM,world);

    if (me == 0) {
      fprintf(fp,BIGINT_FORMAT,update->ntimestep);

      double T_a;
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++) {
            T_a = 0;
            if (nsum_all[ixnode][iynode][iznode] > 0)
              T_a = sum_mass_vsq_all[ixnode][iynode][iznode]/
                (3.0*force->boltz*nsum_all[ixnode][iynode][iznode]/force->mvv2e);
            fprintf(fp," %f",T_a);
          }

      fprintf(fp,"\t");
      for (int ixnode = 0; ixnode < nxnodes; ixnode++)
        for (int iynode = 0; iynode < nynodes; iynode++)
          for (int iznode = 0; iznode < nznodes; iznode++)
            fprintf(fp,"%f ",T_electron[ixnode][iynode][iznode]);
      fprintf(fp,"\n");
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of 3d grid
------------------------------------------------------------------------- */

double FixTTM::memory_usage()
{
  double bytes = 0.0;
  bytes += 5*total_nnodes * sizeof(int);
  bytes += 14*total_nnodes * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixTTM::grow_arrays(int ngrow)
{

 memory->grow(flangevin,ngrow,3,"FixLangevinTTM:flangevin");

}

/* ----------------------------------------------------------------------
  return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixTTM::compute_vector(int n)
{
  double e_energy = 0.0;
  double transfer_energy = 0.0;

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        e_energy +=
          T_electron[ixnode][iynode][iznode]*electronic_specific_heat*
          electronic_density*del_vol;
        transfer_energy +=
          net_energy_transfer_all[ixnode][iynode][iznode]*update->dt;
  }

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTTM::write_restart(FILE *fp)
{
  double *rlist;
  memory->create(rlist,nxnodes*nynodes*nznodes+1,"TTM:rlist");

  int n = 0;
  rlist[n++] = seed;

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        rlist[n++] =  T_electron[ixnode][iynode][iznode];

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rlist,sizeof(double),n,fp);
  }

  memory->destroy(rlist);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTTM::restart(char *buf)
{
  int n = 0;
  double *rlist = (double *) buf;

  // the seed must be changed from the initial seed

  seed = static_cast<int> (0.5*rlist[n++]);

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_electron[ixnode][iynode][iznode] = rlist[n++];

  delete random;
  random = new RanMars(lmp,seed+comm->me);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTTM::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = flangevin[i][0];
  buf[2] = flangevin[i][1];
  buf[3] = flangevin[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixTTM::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  flangevin[nlocal][0] = extra[nlocal][m++];
  flangevin[nlocal][1] = extra[nlocal][m++];
  flangevin[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTTM::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTTM::size_restart(int /*nlocal*/)
{
  return 4;
}

/* ---------------------------------------------------------------------- */

int FixTTM::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
 *    extract thermostat properties
 *    ------------------------------------------------------------------------- */

void *FixTTM::extract(const char *str, int &dim)
{
//  dim = 0;
  if (strcmp(str,"electron_temperature") == 0) {
    dim = 3;
    return T_electron;
  }
  if (strcmp(str,"nodes") == 0) {
    dim = 1;
    return nodes_xyz;
  }
  if (strcmp(str,"seed") == 0) {
    dim = 1;
    return &seed;
  }
  return NULL;
}
