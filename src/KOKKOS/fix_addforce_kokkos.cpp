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

#include <cstring>
#include <cstdlib>
#include "fix_addforce_kokkos.h"
#include "atom_kokkos.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "domain_kokkos.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory_kokkos.h"
#include "error.h"
#include "force.h"
#include "atom_masks.h"
#include "kokkos_base.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixAddForceKokkos<DeviceType>::FixAddForceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixAddForce(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  domainKK = (DomainKokkos *) domain;

  memory->destroy(sforce);
  memoryKK->create_kokkos(k_sforce,sforce,maxatom,4,"addforce:sforce");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixAddForceKokkos<DeviceType>::~FixAddForceKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_sforce,sforce);
  sforce = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixAddForceKokkos<DeviceType>::init()
{
  FixAddForce::init();

  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixAddForceKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space, X_MASK | F_MASK | MASK_MASK | IMAGE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  image = atomKK->k_image.view<DeviceType>();

  int nlocal = atom->nlocal;
  
  //if (vflag)
  //  error->all(FLERR,"vflag not enabled in fix_addforce with Kokkos");
  if (vflag) v_setup(vflag);
  else evflag = 0;

  // update region if necessary

  region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
    DAT::tdual_int_1d k_match = DAT::tdual_int_1d("addforce:k_match",nlocal);
    KokkosBase* regionKKBase = dynamic_cast<KokkosBase*>(region);
    regionKKBase->match_all_kokkos(groupbit,k_match);
    k_match.template sync<DeviceType>();
    d_match = k_match.template view<DeviceType>();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memoryKK->destroy_kokkos(k_sforce,sforce);
    memoryKK->create_kokkos(k_sforce,sforce,maxatom,4,"addforce:sforce");
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  double_4 foriginal_kk;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
    copymode = 1;
    domainKK_prd = Few<double,3>(domain->prd);
    domainKK_h = Few<double,6>(domain->h);
    domainKK_triclinic  = domain->triclinic;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixAddForceConstant>(0,nlocal),*this,foriginal_kk);
    copymode = 0;

  // variable force, wrap with clear/add

  } else {

    atomKK->sync(Host,ALL_MASK); // this can be removed when variable class is ported to Kokkos

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],4,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],4,0);
    if (estyle == ATOM)
      input->variable->compute_atom(evar,igroup,&sforce[0][3],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    if (varflag == ATOM) {  // this can be removed when variable class is ported to Kokkos
      k_sforce.modify<LMPHostType>();
      k_sforce.sync<DeviceType>();
      d_sforce = k_sforce.template view<DeviceType>();
    }

    copymode = 1;
    domainKK_prd = Few<double,3>(domain->prd);
    domainKK_h = Few<double,6>(domain->h);
    domainKK_triclinic  = domain->triclinic;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixAddForceNonConstant>(0,nlocal),*this,foriginal_kk);
    copymode = 0;
  }

  atomKK->modified(execution_space, F_MASK);

  foriginal[0] = foriginal_kk.values[0];
  foriginal[1] = foriginal_kk.values[1];
  foriginal[2] = foriginal_kk.values[2];
  foriginal[3] = foriginal_kk.values[3];

}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixAddForceKokkos<DeviceType>::operator()(TagFixAddForceConstant, const int &i, double_4& foriginal_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    auto unwrap = DomainKokkos::unmap(domainKK_prd,domainKK_h,domainKK_triclinic,x_i,image(i));
    foriginal_kk.values[0] -= xvalue*unwrap[0] + yvalue*unwrap[1] + zvalue*unwrap[2];
    foriginal_kk.values[1] += f(i,0);
    foriginal_kk.values[2] += f(i,1);
    foriginal_kk.values[3] += f(i,2);
    f(i,0) += xvalue;
    f(i,1) += yvalue;
    f(i,2) += zvalue;
    /*
      if (evflag){
          v[0] = xvalue * unwrap[0];
          v[1] = yvalue * unwrap[1];
          v[2] = zvalue * unwrap[2];
          v[3] = xvalue * unwrap[1];
          v[4] = xvalue * unwrap[2];
          v[5] = yvalue * unwrap[2];
          v_tally(i,v);
      }
   */
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixAddForceKokkos<DeviceType>::operator()(TagFixAddForceNonConstant, const int &i, double_4& foriginal_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    
    auto unwrap = DomainKokkos::unmap(domainKK_prd,domainKK_h,domainKK_triclinic,x_i,image(i));
    //if (xstyle == ATOM) xvalue = d_sforce(i,0);
    //if (ystyle == ATOM) yvalue = d_sforce(i,1);
    //if (zstyle == ATOM) zvalue = d_sforce(i,2);

    if (estyle == ATOM) {
       foriginal_kk.values[0] += d_sforce(i,3);
    } else {
      if (xstyle) foriginal_kk.values[0] -= xvalue*unwrap[0];
      if (ystyle) foriginal_kk.values[0] -= yvalue*unwrap[1];
      if (zstyle) foriginal_kk.values[0] -= zvalue*unwrap[2];
    }
    foriginal_kk.values[1] += f(i,0);
    foriginal_kk.values[2] += f(i,1);
    foriginal_kk.values[3] += f(i,2);

    if (xstyle == ATOM) f(i,0) += d_sforce(i,0);
    else {
       if (xstyle) f(i,0) += xvalue;
    }
    if (ystyle == ATOM) f(i,1) += d_sforce(i,1);
    else {
       if (ystyle) f(i,1) += yvalue;
    }
    if (zstyle == ATOM) f(i,2) += d_sforce(i,2);
    else {
       if (zstyle) f(i,2) += zvalue;
    }

    //if (xstyle) f(i,0) += xvalue;
    //if (ystyle) f(i,1) += yvalue;
    //if (zstyle) f(i,2) += zvalue;
/*
    if (evflag) {
       v[0] = xstyle ? xvalue*unwrap[0] : 0.0;
       v[1] = ystyle ? yvalue*unwrap[1] : 0.0;
       v[2] = zstyle ? zvalue*unwrap[2] : 0.0;
       v[3] = xstyle ? xvalue*unwrap[1] : 0.0;
       v[4] = xstyle ? xvalue*unwrap[2] : 0.0;
       v[5] = ystyle ? yvalue*unwrap[2] : 0.0;
       v_tally(i, v);
    }
*/
  }
}

namespace LAMMPS_NS {
template class FixAddForceKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixAddForceKokkos<LMPHostType>;
#endif
}

