/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(ptmmN,FixPTMMN);
// clang-format on
#else

#ifndef LMP_FIX_PTMMN_H
#define LMP_FIX_PTMMN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPTMMN : public Fix {
 public:
  FixPTMMN(class LAMMPS *, int, char **);
//  virtual void initial_integrate(int);
//  virtual void end_of_step();
  int setmask();
  virtual ~FixPTMMN();
  void init();
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
//  void initial_integrate(int) override;
//  void end_of_step() override;
//  void setup(int) override;
  void post_force(int);
//  void pre_force(int);
  void grow_arrays(int);
  double IronBCCMag(double);
  virtual double memory_usage();
//  void compute_vector(int);
//  double compute_scalar();
 private:
  char *id_ptm, *id_stress, *fast, *id_voro;
  class Compute *c_ptm, *c_stress, *c_voro;
  class PairSpin **spin_pair;
  class Pair *pair;
  int npairs, npairspin;
  double **output;
  int nmax;
  int ave_steps;
  double rneigh;
  double fcc_trans;
  double hcp_trans;
  int voro_flag;
  class RanMars *random;
  int seed;
  double **ptm;
  double **voro;
  double **stress;
};
}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix setforce does not exist

Self-explanatory.

E: Variable name for fix setforce does not exist

Self-explanatory.

E: Variable for fix setforce is invalid style

Only equal-style variables can be used.

E: Cannot use non-zero forces in an energy minimization

Fix setforce cannot be used in this manner.  Use fix addforce
instead.

*/
