/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(restart/forces,FixRestartForces)

#else

#ifndef LMP_FIX_RESTART_FORCES_H
#define LMP_FIX_RESTART_FORCES_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRestartForces : public Fix {
 public:
  FixRestartForces(class LAMMPS *, int, char **);
  ~FixRestartForces();
  void post_constructor();  
  int setmask();
  void setup(int);
  void end_of_step();

 private:  
  char *id_fix;
  class FixStore *fix;
  int torque_flag;
};

}

#endif
#endif
