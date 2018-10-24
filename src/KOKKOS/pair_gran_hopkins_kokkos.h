/* ----------------------------------------------------------------------
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Dan S. Bolintineanu (SNL), Adrian K. Turner (LANL)
   ------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pair_gran_hopkins_kokkos.h"
#include "kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "update.h"
#include "force.h"
#include "fix.h"
#include "fix_neigh_history_kokkos.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 1e-10

#define DEBUGID_1 35
#define DEBUGID_2 55
#define DEBUG_TIMESTEP 32696
/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairGranHopkinsKokkos<DeviceType>::PairGranHopkinsKokkos(LAMMPS *lmp) : PairGranHopkins(lmp)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | OMEGA_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | MASK_MASK | 
                  ENERGY_MASK | VIRIAL_MASK | RMASS_MASK | RADIUS_MASK | THICKNESS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
PairGranHopkinsKokkos<DeviceType>::~PairGranHopkinsKokkos()
{
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairGranHopkinsKokkos<DeviceType>::init_style()
{
  if (history && fix_history == NULL) {
    char dnumstr[16];
    sprintf(dnumstr,"%d",3);
    char **fixarg = new char*[4];
    fixarg[0] = (char *) "NEIGH_HISTORY";
    fixarg[1] = (char *) "all";
    if (execution_space == Device)
      fixarg[2] = (char *) "NEIGH_HISTORY/KK/DEVICE";
    else
      fixarg[2] = (char *) "NEIGH_HISTORY/KK/HOST";
    fixarg[3] = dnumstr;
    modify->add_fix(4,fixarg,1);
    delete [] fixarg;
    fix_history = (FixNeighHistory *) modify->fix[modify->nfix-1];
    fix_history->pair = this;
    fix_historyKK = (FixNeighHistoryKokkos<DeviceType> *)fix_history;
  }

  // PairGranHopkins uses PairGranHookeHistory::init_style
  PairGranHookeHistory::init_style();

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with gran/hopkins/kk");
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairGranHopkinsKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  int historyupdate = 1;
  if (update->setupflag) historyupdate = 0;

  // need this?
  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }
 
  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  d_firsttouch = fix_historyKK->d_firstflag;
  d_firsthistory = fix_historyKK->d_firstvalue;

  EV_FLOAT ev; // need this?

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  min_thick = atomKK->k_min_thickness.view<DeviceType>();
  mean_thick = atomKK->k_mean_thickness.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  if (lmp->kokkos->neighflag == HALF) {
     if (force->newton_pair) {
        if (historyupdate) {
            Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,1,1>>(0,inum),*this);
        } else {
            Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,1,0>>(0,inum),*this);
        }
     } else {
        if (historyupdate) {
            Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,0,1>>(0,inum),*this);
        } else {
            Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALF,0,0>>(0,inum),*this);
        }
     }
  }
  else if (lmp->kokkos->neighflag == HALFTHREAD) {
     if (force->newton_pair) {
        if (historyupdate) {
            Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,1,1>>(0,inum),*this);
        } else {
            Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,1,0>>(0,inum),*this);
        }
     } else {
        if (historyupdate) {
            Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,0,1>>(0,inum),*this);
        } else {
            Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHopkinsCompute<HALFTHREAD,0,0>>(0,inum),*this);
        }
     }
  }

}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int HISTORYUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::operator()(TagPairGranHopkinsCompute<NEIGHFLAG,NEWTON_PAIR,HISTORYUPDATE>, const int ii) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  int i = d_ilist[ii];
  int itype = type[i];
  int jnum = d_numneigh[i];

  //F_FLOAT fx = 0.0;
  //F_FLOAT fy = 0.0;

  //F_FLOAT torque_i = 0.0;
  //F_FLOAT torque_j = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    int jtype = type[j];

    F_FLOAT fx = 0.0;
    F_FLOAT fy = 0.0;
    F_FLOAT torque_i = 0.0;
    F_FLOAT torque_j = 0.0;

    F_FLOAT chi1 = d_firsthistory(i,size_history*jj+8);
    F_FLOAT chi2 = d_firsthistory(i,size_history*jj+9);
    if (chi1 >= chi2){ // Un-bonded, chi1 >= chi2
      compute_nonbonded_kokkos<HISTORYUPDATE>(i,j,jj,fx,fy,torque_i,torque_j);
    }
    else { //Bonded
      compute_bonded_kokkos<HISTORYUPDATE,NEWTON_PAIR>(i,j,jj,fx,fy,torque_i,torque_j);
    }

    a_f(i,0) += fx;
    a_f(i,1) += fy;

    // torque induced by tangential force
    a_torque(i,2) += torque_i;

    if (NEWTON_PAIR || j < nlocal){
       a_f(j,0) -= fx;
       a_f(j,1) -= fy;
       a_torque(j,2) += torque_j;
     }
  }
}

template<class DeviceType>
template<int HISTORYUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::compute_nonbonded_kokkos(const int i, const int j, const int jj,
                     F_FLOAT& fx, F_FLOAT& fy, F_FLOAT& torque_i, F_FLOAT& torque_j) const
{
   F_FLOAT r, rinv, nx, ny, radmin;
   F_FLOAT vnnr, vnx, vny;
   F_FLOAT wrz, vtrx, vtry, vtx, vty, vrel;
   F_FLOAT delta, delta_dot;
 
   F_FLOAT fnx, fny;
   F_FLOAT sig_c, hmin;
   F_FLOAT hprime, kp, kr, ke, L, kt0;
   F_FLOAT num, denom, fnmag_plastic, fnmag_elastic, fnmag;
   F_FLOAT ncrossF;
 
   F_FLOAT ndisp, dispmag, scalefac;
   F_FLOAT ftx, fty, ftmag, ftcrit;
   F_FLOAT var1,var2;

   F_FLOAT hstar = 0.3;

   X_FLOAT delx = x(i,0) - x(j,0);
   X_FLOAT dely = x(i,1) - x(j,1);
   X_FLOAT rsq = delx*delx + dely*dely;
   F_FLOAT radsum = radius[i] + radius[j];

   if (rsq >= radsum*radsum){
     d_firsttouch(i,jj) = 0;
     for (int k = 0; k < size_history; k++) {
        d_firsthistory(i,size_history*jj+k) = 0;
     }
   }
   else{
     if (!d_firsttouch(i,jj)){ //If this is first contact
          d_firsttouch(i,jj) = 1;
          d_firsthistory(i,size_history*jj) = 0;
          d_firsthistory(i,size_history*jj+1) = 0;
          d_firsthistory(i,size_history*jj+2) = 0;
          d_firsthistory(i,size_history*jj+3) = 0;
          d_firsthistory(i,size_history*jj+4) = hprime_0;
          d_firsthistory(i,size_history*jj+5) = 0;
      }
   }

   r = sqrt(rsq);
   rinv = 1.0/r;
   nx = delx/r;
   ny = dely/r;

   radmin = MIN(radius[i],radius[j]); 

   L = 2*radmin*(1+(abs(radius[i] - radius[j])/r));

   // relative translational velocity
   V_FLOAT vrx = v(i,0) - v(j,0);
   V_FLOAT vry = v(i,1) - v(j,1);

   delta = radsum - r - d_firsthistory(i,size_history*jj+5);

   if (delta < 0) return ; //delta = 0;

   // Compute tangential force
   // normal component of relative translational velocity
   vnnr = vrx*nx + vry*ny;
   vnx = nx*vnnr;
   vny = ny*vnnr;

   // subtract to compute tangential component of relative translational velocity
   vtrx = vrx - vnx;
   vtry = vry - vny;

   // total relative tangential velocities at contact
   wrz = radius[i]*omega(i,2) + radius[j]*omega(j,2);
   vtx = vtrx + ny*wrz;
   vty = vtry - nx*wrz;

   vrel = vtx*vtx + vty*vty;
   vrel = sqrt(vrel);

   delta_dot = -vnnr;

   // Compute plastic normal force
   hprime = d_firsthistory(i,size_history*jj+4);
   ke = Emod/L*(1/(1/mean_thick(i) + 1/mean_thick(j)));
   hmin = MIN(min_thick(i), min_thick(j));
   if (hprime < hstar){
      kr = 26126*hprime;
      kp = 928*hprime*hprime;
   }
   else{
      kr = kp = hprime*sig_c;
   }

   num = d_firsthistory(i,size_history*jj)/(kp*update->dt) + delta_dot*L + delta*L*ke/damp_normal + ke*kr/(damp_normal*kp);
   denom = 1/(kp*update->dt) + 1/damp_normal*(1+ke/kp);
   fnmag_plastic = num/denom;

   // Elastic normal force
   fnmag_elastic = ke*delta*L + damp_normal*delta_dot*L;

   if (fabs(fnmag_elastic) < fabs(fnmag_plastic))
      fnmag = fnmag_elastic;
   else
      fnmag = fnmag_plastic;

   // should this be here?
   fnmag = fnmag_elastic;

   fnx = fnmag*nx;
   fny = fnmag*ny;

    // update tangential displacement, rotate if needed
    if (HISTORYUPDATE){
       d_firsthistory(i,size_history*jj) = fnx;
       d_firsthistory(i,size_history*jj+1) = fny;
       ndisp = nx*d_firsthistory(i,size_history*jj+2) + ny*d_firsthistory(i,size_history*jj+3);
       var1 = d_firsthistory(i,size_history*jj+2);
       var2 = d_firsthistory(i,size_history*jj+3);
       dispmag =sqrt(var1*var1 + var2*var2); 
       denom = dispmag - ndisp;
       if (ndisp > EPSILON && denom != 0){
          scalefac = dispmag/denom;
          d_firsthistory(i,size_history*jj+2) -= ndisp*nx;
          d_firsthistory(i,size_history*jj+3) -= ndisp*ny;
          d_firsthistory(i,size_history*jj+2) *= scalefac;
          d_firsthistory(i,size_history*jj+3) *= scalefac;
        }
        d_firsthistory(i,size_history*jj+2) += vtx*update->dt;
        d_firsthistory(i,size_history*jj+3) += vty*update->dt;
    }

    var1 = d_firsthistory(i,size_history*jj+2);
    var2 = d_firsthistory(i,size_history*jj+3);
    dispmag =sqrt(var1*var1 + var2*var2); 

    // total tangential force
    kt0 = Gmod/L*(1/(1/mean_thick(i) + 1/mean_thick(j)));
    ftx = - (kt0*var1 + damp_tangential*vtx);
    fty = - (kt0*var2 + damp_tangential*vty);

    ftmag = sqrt(ftx*ftx + fty*fty);
    ftcrit = friction_tangential*fabs(fnmag);
    if (ftmag > ftcrit){
      if (dispmag != 0){
        ftx *= ftcrit/ftmag;
        fty *= ftcrit/ftmag;
        d_firsthistory(i,size_history*jj+2) = -(ftcrit + damp_tangential*vtx)/kt0;
        d_firsthistory(i,size_history*jj+3) = -(ftcrit + damp_tangential*vty)/kt0;
      }
      else ftx = fty = 0;
    }

    //Apply forces
    fx = fnx + ftx;
    fy = fny + fty;

    //Torque
    ncrossF = nx*fty - ny*ftx;
    torque_i = -radius[i]*ncrossF;
    torque_j = -radius[j]*ncrossF;
}

        
template<class DeviceType>
template<int HISTORYUPDATE, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::compute_bonded_kokkos(const int i, const int j, const int jj,
                     F_FLOAT& fx, F_FLOAT& fy, F_FLOAT& torque_i, F_FLOAT& torque_j) const
{

  //See design document for definitions of these variables
  F_FLOAT s1x, s1y, s2x, s2y, mx, my, mmag, mex, mey, bex, bey, rx, ry, rxj, ryj;
  F_FLOAT An, Bn, Cn, Dn, Bt, Ct, Dt, Bnj, Cnj, Dnj, Btj;
  F_FLOAT Fnmag, Ftmag, Nn, Nt, Nnj, Ntj;

  F_FLOAT chi1 = d_firsthistory(i,size_history*jj+8);
  F_FLOAT chi2 = d_firsthistory(i,size_history*jj+9);
  F_FLOAT chidiff, chidiff2, chidiff3;

  F_FLOAT sig_n1, sig_s1;
  F_FLOAT chi_c, chi_t, chi_s1, chi_s2;
  F_FLOAT nprefac, sprefac;
  F_FLOAT damp_prefac, area_bond, fdampx, fdampy, torquedamp;
  F_FLOAT hmin;

  F_FLOAT dvx, dvy;

  if (HISTORYUPDATE){

    // Update bond end points based on particle translations
    d_firsthistory(i,size_history*jj)   += dt*v(i,0);
    d_firsthistory(i,size_history*jj+1) += dt*v(i,1);
    d_firsthistory(i,size_history*jj+2) += dt*v(i,0);
    d_firsthistory(i,size_history*jj+3) += dt*v(i,1);

    d_firsthistory(i,size_history*jj+4) += dt*v(j,0);
    d_firsthistory(i,size_history*jj+5) += dt*v(j,1);
    d_firsthistory(i,size_history*jj+6) += dt*v(j,0);
    d_firsthistory(i,size_history*jj+7) += dt*v(j,1);

    // Update bond end points based on particle rotations
    d_firsthistory(i,size_history*jj)   += -dt*omega(i,2)*(d_firsthistory(i,size_history*jj+1)-x(i,1));
    d_firsthistory(i,size_history*jj+1) +=  dt*omega(i,2)*(d_firsthistory(i,size_history*jj)  -x(i,0));
    d_firsthistory(i,size_history*jj+2) += -dt*omega(i,2)*(d_firsthistory(i,size_history*jj+3)-x(i,1));
    d_firsthistory(i,size_history*jj+3) +=  dt*omega(i,2)*(d_firsthistory(i,size_history*jj+2)-x(i,0));

    d_firsthistory(i,size_history*jj+4) += -dt*omega(j,2)*(d_firsthistory(i,size_history*jj+5)-x(j,1));
    d_firsthistory(i,size_history*jj+5) +=  dt*omega(j,2)*(d_firsthistory(i,size_history*jj+4)-x(j,0));
    d_firsthistory(i,size_history*jj+6) += -dt*omega(j,2)*(d_firsthistory(i,size_history*jj+7)-x(j,1));
    d_firsthistory(i,size_history*jj+7) +=  dt*omega(j,2)*(d_firsthistory(i,size_history*jj+6)-x(j,0));
  }

  //Compute s_1, s_2, m, m_e, b_e
  s1x = d_firsthistory(i,size_history*jj+4) - d_firsthistory(i,size_history*jj);
  s1y = d_firsthistory(i,size_history*jj+5) - d_firsthistory(i,size_history*jj+1);
  s2x = d_firsthistory(i,size_history*jj+6) - d_firsthistory(i,size_history*jj+2);
  s2y = d_firsthistory(i,size_history*jj+7) - d_firsthistory(i,size_history*jj+3);

  mx = d_firsthistory(i,size_history*jj+2) + 0.5*s2x - d_firsthistory(i,size_history*jj) - 0.5*s1x;
  my = d_firsthistory(i,size_history*jj+3) + 0.5*s2y - d_firsthistory(i,size_history*jj+1) - 0.5*s1y;
  mmag = sqrt(mx*mx + my*my);
  mex = mx/mmag;
  mey = my/mmag;

  bex = mey;
  bey = -mex;

  rx = d_firsthistory(i,size_history*jj)   + 0.5*s1x - x(i,0);
  ry = d_firsthistory(i,size_history*jj+1) + 0.5*s1y - x(i,1);

  //Compute forces and torques
  Dn = s1x*bex + s1y*bey;
  Cn = s2x*bex + s2y*bey - Dn;
  Bn = rx*bey - ry*bex;
  An = mx*bey - my*bex;

  Dt = s1x*mex + s1y*mey;
  Ct = s2x*mex + s2y*mey - Dt;
  Bt = rx*mey - ry*mex;

  chidiff = chi2-chi1;
  chidiff2 = 0.5*(chi2*chi2 - chi1*chi1);
  chidiff3 = MathConst::THIRD*(chi2*chi2*chi2 - chi1*chi1*chi1);
  F_FLOAT kn0 = d_firsthistory(i,size_history*jj+11)*Emod/d_firsthistory(i,size_history*jj+10);
  F_FLOAT kt0 = d_firsthistory(i,size_history*jj+11)*Gmod/d_firsthistory(i,size_history*jj+10);
  nprefac = d_firsthistory(i,size_history*jj+10)*kn0;
  sprefac = d_firsthistory(i,size_history*jj+10)*kt0;

  Fnmag = nprefac*(Dn*chidiff + Cn*chidiff2);
  Ftmag = sprefac*(Dt*chidiff + Ct*chidiff2);
  Nn = nprefac*(An*Cn*chidiff3 + (Bn*Cn+An*Dn)*chidiff2 + Bn*Dn*chidiff);
  Nt = Bt*Ftmag;

  //Damping force
  area_bond = d_firsthistory(i,size_history*jj+10)*d_firsthistory(i,size_history*jj+11)*chidiff;
  damp_prefac = damp_bonded*area_bond;

  dvx = v(j,0) - v(i,0);
  dvy = v(j,1) - v(i,1);
  fdampx = damp_prefac*dvx;
  fdampy = damp_prefac*dvy;

  //Damping torque
  torquedamp = -damp_bonded*area_bond*area_bond*(omega(i,2)-omega(j,2));

  //Update forces and torque
  fx = Fnmag*bex + Ftmag*mex;
  fy = Fnmag*bey + Ftmag*mey;

  //Ensure that damping does not cause force to change sign
 /* if (fx < 0 && fdampx > 0){
    fx = MIN(fx+fdampx, 0);
  }
  else if (fx > 0 && fdampx < 0){
    fx = MAX(fx+fdampx, 0);
  }
  else */
  fx = fx+fdampx;

  /*if (fy < 0 && fdampy > 0){
    fy = MIN(fy+fdampy, 0);
  }
  else if (fy > 0 && fdampy < 0){
    fy = MAX(fy+fdampy, 0);
  }
  else*/
  fy = fy+fdampy;

  torque_i += Nn + Nt + torquedamp;

  if (NEWTON_PAIR || j < nlocal){

    rxj = d_firsthistory(i,size_history*jj) + 0.5*s1x - x(j,0);
    ryj = d_firsthistory(i,size_history*jj+1) + 0.5*s1y - x(j,1);
    Bnj = rxj*bey - ryj*bex;
    Cnj = -Cn;
    Dnj = -Dn;
    Btj = rxj*mey - ryj*mex;

    Nnj = nprefac*(An*Cnj*chidiff3 + (Bnj*Cnj+An*Dnj)*chidiff2 + Bnj*Dnj*chidiff);
    Ntj = -Btj*Ftmag;
    torque_j += Nnj + Ntj - torquedamp;
  }

  //Update chi1, chi2
  hmin = MIN(min_thick(i), min_thick(j));
  if (HISTORYUPDATE){
    F_FLOAT c1, c2;
    c1 = d_firsthistory(i,size_history*jj+8);
    c2 = d_firsthistory(i,size_history*jj+9);
    //update_chi(kn0, kt0, Dn, Cn, Dt, Ct, hmin, d_firsthistory(i,size_history*jj+8), d_firsthistory(i,size_history*jj+9));
    update_chi(kn0, kt0, Dn, Cn, Dt, Ct, hmin, c1, c2);
    d_firsthistory(i,size_history*jj+8) = c1;
    d_firsthistory(i,size_history*jj+9) = c2;
    d_firsttouch(i,jj) = 1;
    if (d_firsthistory(i,size_history*jj+8) >= d_firsthistory(i,size_history*jj+9)){ //Bond just broke
        F_FLOAT dx = x(i,0) - x(j,0);
        F_FLOAT dy = x(i,1) - x(j,1);
        F_FLOAT rij = sqrt(dx*dx + dy*dy);
        F_FLOAT delta_0 = radius(i) + radius(j) - rij;
        if (delta_0 < 0) delta_0 = 0;
        for (int k = 0; k < size_history; ++k)
          d_firsthistory(i,size_history*jj+k) = 0;
        d_firsthistory(i,size_history*jj+4) = hprime_0;
        d_firsthistory(i,size_history*jj+5) = delta_0;
    }
  }
}


template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairGranHopkinsKokkos<DeviceType>::update_chi(F_FLOAT kn0, F_FLOAT kt0, F_FLOAT Dn, F_FLOAT Cn, 
                                   F_FLOAT Dt, F_FLOAT Ct, F_FLOAT hmin, F_FLOAT chi1, F_FLOAT chi2) const
{

  F_FLOAT sig_n1 = kn0*(Dn + Cn*chi1);
  F_FLOAT sig_s1 = kt0*(Dt + Ct*chi1);
  F_FLOAT sig_n2 = kn0*(Dn + Cn*chi2);
  F_FLOAT sig_s2 = kt0*(Dt + Ct*chi2);

  // function pointers for greater efficiency?
  F_FLOAT sig_c;
  if (strcmp(sig_c0_type,"constant") == 0) {
    sig_c = sig_c0;
  } else if (strcmp(sig_c0_type,"KovacsSodhi") == 0) {
    sig_c = sig_c0*pow(hmin,(2.0/3.0)) * 1000.0;
  } else {
    error->all(FLERR,"Unknown sig_c0_type");
  }
  F_FLOAT sig_t;
  if (strcmp(sig_t0_type,"constant") == 0) {
    sig_t = sig_t0;
  } else if (strcmp(sig_t0_type,"multiply_sig_c0") == 0) {
    sig_t = sig_t0 * sig_c;
  } else {
    error->all(FLERR,"Unknown sig_t0_type");
  }

  F_FLOAT denom;
  sig_t = -sig_t;

  F_FLOAT c1, c2;
  c1 = chi1;
  c2 = chi2;
  //Check for purely tensile/compressive failure
  if (sig_n1 < sig_t || sig_n2 < sig_t || sig_n1 > sig_c || sig_n2 > sig_c){
    if (sig_n1 < sig_t){
      if (Cn != 0) chi1 = (sig_t/kn0 - Dn)/Cn;
      else {
        chi1 = chi2;
        return;
      }
    }
    if (sig_n2 < sig_t) chi2 = (sig_t/kn0 - Dn)/Cn;
    //if Cn==0, then sig_n1 = sig_n2, and function would've returned above

    if (sig_n1 > sig_c){
      if (Cn != 0) chi1 = (sig_c/kn0 - Dn)/Cn;
      else{
        chi1 = chi2;
        return;
      }
    }
    if (sig_n2 > sig_c) chi2 = (sig_c/kn0 - Dn)/Cn;
    //if Cn==0, function would've returned above
  }

  //Check for 'cohesion' shear failure at chi1
  if (sig_s1 > tanphi*sig_n1 - tanphi*sig_t){
    denom = tanphi*kn0*Cn - kt0*Ct;
    if (denom != 0) chi1 = (tanphi*(sig_t-kn0*Dn)+kt0*Dt)/denom;
    else if (Cn != 0 && Ct != 0){
      //No need to treat case where Cn = Ct = 0, since sig_s = 0 for that case,
      // and it would've been picked up above
      //For the rare case of kt0*Ct = -kn0*tanphi*Cn, yield criterion is independent of chi,
      // therefore bond breaks.
      chi1 = chi2;
      return;
    }
  }

  if (sig_s1 < -tanphi*sig_n1 + tanphi*sig_t){
    denom = -kn0*tanphi*Cn - kt0*Ct;
    if (denom != 0) chi1 = (kt0*Dt+tanphi*(kn0*Dn-sig_t))/denom;
    else if (Cn != 0 && Ct != 0){
      chi1 = chi2;
      return;
    }
  }

  //Check for 'cohesion' shear failure at chi2
  if (sig_s2 > tanphi*sig_n2 - tanphi*sig_t){
    denom = tanphi*kn0*Cn - kt0*Ct;
    if (denom != 0) chi2 = (tanphi*(sig_t-kn0*Dn)+kt0*Dt)/denom;
    else if (Cn != 0 && Ct != 0){
      //No need to treat case where Cn = Ct = 0, since sig_s = 0 for that case,
      // and it would've been picked up above
      //For the rare case of kt0*Ct = -kn0*tanphi*Cn, yield criterion is independent of chi,
      // therefore bond breaks.
      chi2 = chi1;
      return;
    }
  }
  if (sig_s2 < -tanphi*sig_n2 + tanphi*sig_t){
    denom = -kn0*tanphi*Cn - kt0*Ct;
    if (denom != 0) chi2 = (kt0*Dt+tanphi*(kn0*Dn-sig_t))/denom;
    else if (Cn != 0 && Ct != 0){
      chi2 = chi1;
      return;
    }
  }  
  if (chi1 < 0 || chi1 > 1 || chi2 < 0 || chi2 > 1)
    chi1 = chi2 = 0;

}


namespace LAMMPS_NS {
template class PairGranHopkinsKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class PairGranHopkinsKokkos<LMPHostType>;
#endif
}

