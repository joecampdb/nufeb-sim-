/* -*- c++ -*- ----------------------------------------------------------
   NUFEB - Individual-based Modelling for Microbial Communities
   Custom fix: Anammox growth with Hg2+ non-competitive inhibition

   Extends fix_growth_anammox to include:
     mu = mu_max * [NH4]/(K_NH4+[NH4]) * [NO2]/(K_NO2+[NO2])
                 * K_O2/(K_O2+[O2])
                 * K_I/(K_I+[Hg2+])          <-- inhibition term
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nufeb/growth/anammox_hg,FixGrowthAnammoxHg)

#else

#ifndef LMP_FIX_GROWTH_ANAMMOX_HG_H
#define LMP_FIX_GROWTH_ANAMMOX_HG_H

#include "fix_growth.h"

namespace LAMMPS_NS {

class FixGrowthAnammoxHg: public FixGrowth {
 public:
  FixGrowthAnammoxHg(class LAMMPS *, int, char **);
  virtual ~FixGrowthAnammoxHg() {}

  virtual void update_atoms();
  virtual void update_cells();

 protected:
  // standard anammox substrates
  int inh4;
  int io2;
  int ino2;
  int ino3;
  int in2;

  // mercury substrate index
  int ihg;

  // half-saturation (affinity) constants
  double nh4_affinity;
  double o2_affinity;
  double no2_affinity;

  // Hg2+ inhibition constant (kg/m3)
  // growth = 0 when [Hg2+] >> K_I
  // growth = mu_max when [Hg2+] << K_I
  double hg_inhibition;

  // kinetic parameters
  double growth;
  double yield;
  double maintain;
  double decay;
};

}

#endif
#endif
