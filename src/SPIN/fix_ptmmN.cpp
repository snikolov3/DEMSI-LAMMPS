// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_ptmmN.h"
#include <iostream>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
//#include "compute_ptm_atom.h"
#include "compute_chunk_atom.h"
#include "fix_ave_chunk.h"
#include "error.h"
#include "random_mars.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair_spin_exchange.h"
#include "pair_hybrid.h"
#include "pair_spin.h"
#include "force.h"
#include "comm.h"
#include "math.h"
#include "math_const.h"
#include "fix_neigh_history.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixPTMMN::FixPTMMN(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_ptm(nullptr), id_voro(nullptr), id_stress(nullptr), 
  spin_pair(nullptr), pair(nullptr), output(nullptr), fast(nullptr), random(nullptr)
{
  if (narg != 12)
     error->all(FLERR,"Illegal fix ptm command");

  vector_flag = 1;
  size_vector = 9;
  extvector = 1;

  // store phase/atom, vol/atom, stress/atom IDs used by fix ptmm computation
  // insure they are valid for these computations
  comm_forward = 18;
  id_ptm = utils::strdup(arg[3]);
  id_voro = utils::strdup(arg[4]);
  id_stress = utils::strdup(arg[5]);
  ave_steps = utils::inumeric(FLERR,arg[6],false,lmp);
  rneigh = utils::inumeric(FLERR,arg[7],false,lmp);
  fast = utils::strdup(arg[8]);
  fcc_trans = utils::numeric(FLERR,arg[9],false,lmp);
  hcp_trans = utils::numeric(FLERR,arg[10],false,lmp);
  seed = utils::inumeric(FLERR,arg[11],false,lmp);

  int iptm = modify->find_compute(id_ptm);
  int ivoro = modify->find_compute(id_voro);
  int istress = modify->find_compute(id_stress);

  if (iptm < 0 || ivoro < 0 || istress < 0)
     error->all(FLERR,"Could not find compute stress/atom ID");
  if (strcmp(modify->compute[iptm]->style,"ptm/atom") != 0)
     error->all(FLERR,"Fix ptm compute ID does not compute ptm/atom");
  if (strcmp(modify->compute[ivoro]->style,"voronoi/atom") != 0)
     error->all(FLERR,"Fix ptm compute ID does not compute voronoi/atom");
  if (modify->compute[istress]->pressatomflag != 1)
     error->all(FLERR,"Fix ptm compute ID does not compute stress/atom");

  FixPTMMN::grow_arrays(atom->nmax);
  nmax = 0;

  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

FixPTMMN::~FixPTMMN()
{
  delete [] id_ptm;
  delete [] id_voro;
  delete [] id_stress;
  delete [] spin_pair;
  delete [] output;
  delete [] fast;
  delete random;
}

/* ---------------------------------------------------------------------- */

void FixPTMMN::grow_arrays(int ngrow)
{
  memory->grow(output,ngrow,18,"ptmmN:output");
}

/* ---------------------------------------------------------------------- */

void FixPTMMN::init()
{
  // error checks

  int nlocal = atom->nlocal;
  int nmax_ini = atom->nmax;
  int iptm = modify->find_compute(id_ptm);
  int ivoro = modify->find_compute(id_voro);
  int istress = modify->find_compute(id_stress);
  if (iptm < 0 || istress < 0 || ivoro < 0)
    error->all(FLERR,"Could not find one or more of provided compute IDs");

  c_voro = modify->compute[ivoro];
  c_ptm = modify->compute[iptm];
  c_stress = modify->compute[istress];
  npairspin = 0;
  PairHybrid *hybrid = (PairHybrid *)force->pair_match("^hybrid",0);
  if (force->pair_match("^spin",0,0)) {        // only one Pair/Spin style
    pair = force->pair_match("^spin",0,0);
    if (hybrid == nullptr) npairs = 1;
    else npairs = hybrid->nstyles;
    npairspin = 1;
  } else if (force->pair_match("^spin",0,1)) { // more than one Pair/Spin style
    pair = force->pair_match("^spin",0,1);
    if (hybrid == nullptr) npairs = 1;
    else npairs = hybrid->nstyles;
    for (int i = 0; i<npairs; i++) {
      if (force->pair_match("^spin",0,i)) {
        npairspin ++;
      }
    }
  }

  if (npairspin > 0) {
    spin_pair = new PairSpin*[npairspin];
  }

  int count1 = 0;
  if (npairspin == 1) {
    count1 = 1;
    spin_pair[0] = (PairSpin *) force->pair_match("^spin",0,0);
  } else if (npairspin > 1) {
    for (int i = 0; i<npairs; i++) {
        if (force->pair_match("^spin",0,i)) {
           spin_pair[count1] = (PairSpin *) force->pair_match("^spin",0,i);
           count1++;
        }
    }
  }

  for (int ind = 0; ind < nmax_ini; ind++) {
       output[ind][0] = 0.0; // local summation of stresses
       output[ind][1] = 0.0; // local summation of volumes
       output[ind][2] = 0.0; // ind atom voronoi volume
       output[ind][3] = 0.0; // ind atom voronoi neighbor count
       output[ind][4] = 0.0; // ind atom ptm type
       output[ind][5] = 0.0; // ind atom ptm rmsd
       output[ind][6] = 0.0; // ind atom ptm interatomic distance
       output[ind][7] = 0.0; // ind atom ptm qw
       output[ind][8] = 0.0; // ind atom ptm qx
       output[ind][9] = 0.0; // ind atom ptm qy
       output[ind][10] = 0.0; // ind atom ptm qz
       output[ind][11] = 0.0; // ind atom ptm liquid yes/no
       output[ind][12] = 0.0; // ind atom xx stress
       output[ind][13] = 0.0; // ind atom yy stress
       output[ind][14] = 0.0; // ind atom zz stress
       output[ind][15] = 0.0; // ind atom xy stress
       output[ind][16] = 0.0; // ind atom xz stress
       output[ind][17] = 0.0; // ind atom yz stress
  }

  modify->addstep_compute_all(update->ntimestep);
}

/* ---------------------------------------------------------------------- */

int FixPTMMN::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
//  mask |= INITIAL_INTEGRATE;
//  mask |= END_OF_STEP;
  return mask;
}

//void FixPTMMN::setup(int /*vflag*/)
//{ 
//  end_of_step();
//}

/* ---------------------------------------------------------------------- */

//void FixPTMMN::pre_force(int /*vflag*/)
void FixPTMMN::post_force(int /*vflag*/)
//void FixPTMMN::initial_integrate(int /*vflag*/)
//void FixPTMMN::end_of_step()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int nlocal = atom->nlocal;
  int nnmax = atom->nmax;
  int *mask = atom->mask;
  int *type = atom->type;
  double **sp = atom->sp;
  double **x = atom->x;
  tagint *tag = atom->tag;
  int *ilist,*jlist,*numneigh,**firstneigh;
  inum = spin_pair[0]->list->inum;
  ilist = spin_pair[0]->list->ilist;
  numneigh = spin_pair[0]->list->numneigh;
  firstneigh = spin_pair[0]->list->firstneigh;

  if (atom->nmax > nmax) {
     memory->destroy(output);
     nmax = atom->nmax;
     memory->create(output,nmax,18,"ptmmN/output");
     array_atom = output;
  }

  if (strcmp(fast,"no") == 0 || (update->ntimestep % ave_steps == 0)){
     if (!(c_ptm->invoked_flag & Compute::INVOKED_PERATOM)) {
        c_ptm->compute_peratom();
        c_ptm->invoked_flag |= Compute::INVOKED_PERATOM;
     }
     if (!(c_voro->invoked_flag & Compute::INVOKED_PERATOM)) {
        c_voro->compute_peratom();
        c_voro->invoked_flag |= Compute::INVOKED_PERATOM;
     }
  }  // (strcmp(fast,"no") == 0

  if (!(c_stress->invoked_flag & Compute::INVOKED_PERATOM)) {
     int temp = update->vflag_atom;
     update->vflag_atom = update->ntimestep;
     c_stress->compute_peratom();
     c_stress->invoked_flag |= Compute::INVOKED_PERATOM;
     update->vflag_atom = temp;
  }

  if (strcmp(fast,"no") == 0 || (update->ntimestep % ave_steps == 0)){
     voro = c_voro->array_atom; // output array from compute_voronoi
     ptm = c_ptm->array_atom; // output array from compute_ptm
     stress = c_stress->array_atom; // output array from compute_stress
  }

  for (int ind = 0; ind < inum; ind++) {
     i = ilist[ind];
     itype = type[i];
     output[i][0] = 0.0; // local summation of stresses
     output[i][1] = 0.0; // local summation of volumes
     output[i][2] = voro[i][0]; // ind atom voronoi volume
     output[i][3] = voro[i][1]; // ind atom voronoi neighbor count
     output[i][4] = ptm[i][0]; // ind atom ptm type
     output[i][5] = ptm[i][1]; // ind atom ptm rmsd
     output[i][6] = ptm[i][2]; // ind atom ptm interatomic distance
     output[i][7] = ptm[i][3]; // ind atom ptm qw
     output[i][8] = ptm[i][4]; // ind atom ptm qx
     output[i][9] = ptm[i][5]; // ind atom ptm qy
     output[i][10] = ptm[i][6]; // ind atom ptm qz
     output[i][11] = ptm[i][7]; // ind atom ptm liquid yes/no
     output[i][12] = stress[i][0]; // ind atom xx stress
     output[i][13] = stress[i][1]; // ind atom yy stress
     output[i][14] = stress[i][2]; // ind atom zz stress
     output[i][15] = stress[i][3]; // ind atom xy stress
     output[i][16] = stress[i][4]; // ind atom xz stress
     output[i][17] = stress[i][5]; // ind atom yz stress
  }

  comm->forward_comm(this);
  double magn = 0;
  double loc_press = 0;
  double rsq = rneigh*rneigh;
  double dist = 0;
  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    output[i][0] += (output[i][12] + output[i][13] + output[i][14])/3;
    if (strcmp(fast,"no") == 0){
       output[i][1] += output[i][2]; // volume from voronoi compute
    }
    else if (update->ntimestep % ave_steps == 0 && strcmp(fast,"yes") == 0){
       output[i][1] += output[i][2]*ave_steps;
    }

    // loop on neighbors
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      dist = pow(x[i][0] - x[j][0],2.0) + pow(x[i][1] - x[j][1],2.0) + pow(x[i][2] - x[j][2],2.0);
      if (dist <= rsq){
         output[i][0] += (output[j][12] + output[j][13] + output[j][14])/3;
         if (strcmp(fast,"no") == 0){
            output[i][1] += output[j][2];
         }
         else if (update->ntimestep % ave_steps == 0 && strcmp(fast,"yes") == 0){
            output[i][1] += output[j][2]*ave_steps;
         }
      } // if (dist <= rsq)
    } // for (jj = 0; jj < jnum; jj++)
  } // for (ii = 0; ii < inum; ii++)

  if (update->ntimestep % ave_steps == 0) {
     for (ii = 0; ii < inum; ii++) {
	 i = ilist[ii];
	 itype = type[i];
         loc_press = (-output[i][0] / (output[i][1] * 10000) ); // in GPa
         double oldmag = sp[i][3];
         double mag_zero_flag = 0;
         if (oldmag == 0) 
            mag_zero_flag == 1;

         if (output[i][4] == 3){ //bcc
	    if (loc_press > 530){
	       magn = -0.00275014*loc_press + 2.02263392;
	    }
	    else{
	       magn = (0.000000000049643876*pow(loc_press,4.0) - 0.000000066656241*pow(loc_press,3.0) + 0.000030478248*pow(loc_press,2.0) - 0.0081920216*loc_press + 2.3346597)*(2.2/2.3346597);
	    }
	 }
         else if (output[i][4] == -1){ //liquid 
            //magn = 0;
            magn = (0.000000000049643876*pow(loc_press,4.0) - 0.000000066656241*pow(loc_press,3.0) + 0.000030478248*pow(loc_press,2.0) - 0.0081920216*loc_press + 2.3346597)*(2.2/2.3346597);
         }

         else if (output[i][4] == 1){ //fcc
            magn = IronBCCMag(loc_press) - IronBCCMag(loc_press)*(1/(1+exp(-loc_press+fcc_trans)));
	 }
         else if (output[i][4] == 2){ //hcp 
            magn = IronBCCMag(loc_press) - IronBCCMag(loc_press)*(1/(1+exp(-loc_press+hcp_trans)));
	 }
	 else {
            magn = (0.000000000049643876*pow(loc_press,4.0) - 0.000000066656241*pow(loc_press,3.0) + 0.000030478248*pow(loc_press,2.0) - 0.0081920216*loc_press + 2.3346597)*(2.2/2.3346597);
         }

	 if (magn < 0.0){
	    magn = 0.0;
	 }

         double sunit = sp[i][0] + sp[i][1] + sp[i][2];
         if (magn > 0 && output[i][4] == 2){
            double erase1 = 5;
         }
         
/*         if (sunit == 0 && magn > 0){
            double mx = random->gaussian();
            double my = random->gaussian();
            double mz = random->gaussian();
            double msq = mx*mx + my*my + mz*mz;
     jlist = firstneigh[i];
     jnum = numneigh[i];
            double scale = 1.0/sqrt(msq);
            sp[indexI][0] = mx*scale;
            sp[indexI][1] = my*scale;
            sp[indexI][2] = mz*scale;
         } */

         if (sunit > 0 && magn == 0){
            sp[i][0] = 0;
            sp[i][1] = 0;
            sp[i][2] = 0;
         }
         sp[i][3] = magn;
         output[i][0] = 0;
         output[i][1] = 0;
     }
  }
}

/* ---------------------------------------------------------------------- */

double FixPTMMN::IronBCCMag(double locPress){
            double magnorm = 0;
            if (locPress > 530){
               magnorm =  -0.00275014*locPress + 2.02263392;
            }
            else{
               magnorm = (0.000000000049643876*pow(locPress,4.0) - 0.000000066656241*pow(locPress,3.0) + 0.000030478248*pow(locPress,2.0) - 0.0081920216*locPress + 2.3346597)*(2.2/2.3346597);
            }
            if (magnorm < 0){
               magnorm = 0;
            }
            return magnorm;
}

/* ---------------------------------------------------------------------- */

int FixPTMMN::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{ 
  int i,j,m;
  
  m = 0; 
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = output[j][0];
    buf[m++] = output[j][1];
    buf[m++] = output[j][2];
    buf[m++] = output[j][3];
    buf[m++] = output[j][4];
    buf[m++] = output[j][5];
    buf[m++] = output[j][6];
    buf[m++] = output[j][7];
    buf[m++] = output[j][8];
    buf[m++] = output[j][9];
    buf[m++] = output[j][10];
    buf[m++] = output[j][11];
    buf[m++] = output[j][12];
    buf[m++] = output[j][13];
    buf[m++] = output[j][14];
    buf[m++] = output[j][15];
    buf[m++] = output[j][16];
    buf[m++] = output[j][17];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixPTMMN::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    output[i][0] = buf[m++];
    output[i][1] = buf[m++];
    output[i][2] = buf[m++];
    output[i][3] = buf[m++];
    output[i][4] = buf[m++];
    output[i][5] = buf[m++];
    output[i][6] = buf[m++];
    output[i][7] = buf[m++];
    output[i][8] = buf[m++];
    output[i][9] = buf[m++];
    output[i][10] = buf[m++];
    output[i][11] = buf[m++];
    output[i][12] = buf[m++];
    output[i][13] = buf[m++];
    output[i][14] = buf[m++];
    output[i][15] = buf[m++];
    output[i][16] = buf[m++];
    output[i][17] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double FixPTMMN::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) atom->nmax * 2 * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */
