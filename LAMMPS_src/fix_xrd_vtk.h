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

#ifdef FIX_CLASS

FixStyle(xrd/vtk,FixXRDvtk)

#else

#ifndef LMP_FIX_XRD_VTK_H
#define LMP_FIX_XRD_VTK_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixXRDvtk : public Fix {
 public:
  FixXRDvtk(class LAMMPS *, int, char **);
  ~FixXRDvtk();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_vector(int);
  void reset_timestep(bigint);


 private:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  bigint nvalid;
  int which;
  char *ids;
  FILE *fp;
  int nrows;
  int ncols;

  int ave,nwindow,nsum,startstep;
  int overwrite;
  long filepos;

  int norm,iwindow,window_limit;
  double *vector;
  double *vector_total;
  double **vector_list;

  void invoke_scalar(bigint);
  void invoke_vector(bigint);
  void options(int, char **);

  bigint nextvalid();
  
  class ComputeXRD *compute_xrd;

  double  Zone[3];           // Zone axis to view XRD
  double  R_Ewald;           // Radius of Ewald sphere (distance units)
  double  lambda;            // Radiation wavelength (distance units)
  double  dK[3];             // spacing of reciprocal points in each dimension
  double  Kmax;              // Maximum reciprocal distance to explore    
  int     Knmax[3];          // maximum integer value for K points in each dimension
  int     Knmin[3];          // minimum integer value for K points in each dimension

//int     KnSlice[6];        // min 0-2 max 2-5 hkl index using zone
  double  Max2Theta;         
  double  Min2Theta;
  double  c[3];              // Parameters controlling resolution of reciprocal space explored
  double  prd_inv[3];        // Inverse spacing of unit cell

  int     argindex; 
  char    *filename;         // user-specified file
  int     nOutput; 
  int     Dim[3];
  bool    manual;            // Turn on manual reciprocal map

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for fix ave/time does not exist

Self-explanatory.

E: Fix ID for fix ave/time does not exist

Self-explanatory.

E: Invalid fix ave/time off column

Self-explantory.

E: Fix ave/time compute does not calculate a scalar

Self-explantory.

E: Fix ave/time compute does not calculate a vector

Self-explantory.

E: Fix ave/time compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ave/time compute does not calculate an array

Self-explanatory.

E: Fix ave/time compute array is accessed out-of-range

An index for the array is out of bounds.

E: Fix ave/time fix does not calculate a scalar

Self-explanatory.

E: Fix ave/time fix does not calculate a vector

Self-explanatory.

E: Fix ave/time fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix for fix ave/time not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/time
is requesting a value on a non-allowed timestep.

E: Fix ave/time fix does not calculate an array

Self-explanatory.

E: Fix ave/time fix array is accessed out-of-range

An index for the array is out of bounds.

E: Variable name for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time variable is not equal-style variable

Self-explanatory.

E: Fix ave/time cannot use variable with vector mode

Variables produce scalar values.

E: Fix ave/time columns are inconsistent lengths

Self-explanatory.

E: Fix ave/time cannot set output array intensive/extensive from these inputs

One of more of the vector inputs has individual elements which are
flagged as intensive or extensive.  Such an input cannot be flagged as
all intensive/extensive when turned into an array by fix ave/time.

E: Cannot open fix ave/time file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Fix ave/time missed timestep

You cannot reset the timestep to a value beyond where the fix
expects to next perform averaging.

*/
