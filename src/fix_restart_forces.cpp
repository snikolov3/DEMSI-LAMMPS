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

#include <stdlib.h>
#include <string.h>
#include "fix_restart_forces.h"
#include "atom.h"
#include "group.h"
#include "modify.h"
#include "domain.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "update.h"
#include "output.h"
#include "error.h"
#include "fix_store.h"

using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixRestartForces::FixRestartForces(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  id_fix(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal fix restart/forces command");

  // if triclinic is set, in the run setup calls to x2lamda and lamda2x will ruin bfb even with this fix
  // Note, the addition of other fixes can easily also ruin bfb restarts for instance ones with RNG or
  // fix/deform determines ho it strains the box based on the initial size 
  if(domain->triclinic) 
      error->warning(FLERR,"Triclinic box prevents bit-for-bit restarts");

}

/* ---------------------------------------------------------------------- */

FixRestartForces::~FixRestartForces()
{
  if (modify->nfix) modify->delete_fix(id_fix);    
}

/* ---------------------------------------------------------------------- */

void FixRestartForces::post_constructor()
{
  //Alternatively, one could use a local array but atoms could be sorted in setup
  //One would also need to register grow and restart callbacks
  int n = strlen(id) + strlen("_restart_forces") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_restart_forces");

  char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "1";

  if (atom->torque_flag) {
    torque_flag = 1;
    newarg[5] = (char *) "6";
  } else {
    newarg[5] = (char *) "3";
  }

  modify->add_fix(6,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;  
}


/* ---------------------------------------------------------------------- */

int FixRestartForces::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;  
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRestartForces::setup(int vflag)
{
  // Check if values in fix store were loaded from  restart file
  if(fix->restart_reset){ 
    double **f = atom->f;
    double **torque = atom->torque;
    double **values = fix->astore;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int i,m;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (m = 0; m < 3; m++) {
          f[i][m] = values[i][m];
          if(torque_flag) torque[i][m] = values[i][m+3];
        }
      }    
    }  
  }
  
}

/* ---------------------------------------------------------------------- */

void FixRestartForces::end_of_step()
{
  if (not output->restart_flag) return;
   
  if(update->ntimestep == output->next_restart) {
    double **f = atom->f;
    double **torque = atom->torque;
    double **values = fix->astore;     
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    for(int i = 0; i < nlocal; i ++){
      if (mask[i] & groupbit) {
        for(int m = 0; m < 3; m++){
          values[i][m] = f[i][m];
          if(torque_flag) values[i][m+3] = torque[i][m];
        }
      }
    } 
  }
}


