/* -*- c++ -*- ----------------------------------------------------------
   NUFEB - Individual-based Modelling for Microbial Communities
   Custom fix: Pseudomonas pseudoalcaligenes CECT5344

   Cyanotrophic bacterium that:
     1. Degrades SCN- as sole nitrogen source (Monod kinetics)
     2. Requires O2 (obligate aerobe)
     3. Produces NH3 from SCN- hydrolysis (feeds anammox/AOB)
     4. Accumulates Hg2+ intracellularly (mercury sink)
     5. Has high Hg2+ tolerance (K_I >> anammox)

   Kinetics:
     mu = mu_max * [SCN]/(K_SCN+[SCN]) * [O2]/(K_O2+[O2])
                 * K_I/(K_I+[Hg2+])

   Based on literature:
     Ibanez et al. 2023, Microbiology Spectrum (proteomics)
     Luque-Almagro et al. 2010, J Hazard Mater (kinetics)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nufeb/growth/cect5344,FixGrowthCECT5344)

#else

#ifndef LMP_FIX_GROWTH_CECT5344_H
#define LMP_FIX_GROWTH_CECT5344_H

#include "fix_growth.h"

namespace LAMMPS_NS {

class FixGrowthCECT5344: public FixGrowth {
 public:
  FixGrowthCECT5344(class LAMMPS *, int, char **);
  virtual ~FixGrowthCECT5344() {}

  virtual void update_atoms();
  virtual void update_cells();

 protected:
  // substrate indices on the grid
  int iscn;           // thiocyanate (electron donor / N source)
  int io2;            // dissolved oxygen (electron acceptor)
  int inh4;           // ammonium (product of SCN hydrolysis)
  int ihg;            // mercury(II) ions (inhibitor + bioaccumulated)

  // half-saturation constants (kg/m3)
  double scn_affinity;
  double o2_affinity;

  // Hg2+ inhibition constant (kg/m3)
  double hg_inhibition;

  // Hg2+ intracellular accumulation rate constant (m3 kg-1 s-1)
  // removal flux: reac[hg] -= hg_uptake * [Hg2+] * dens_cect
  double hg_uptake;

  // NH3 yield from SCN- (kg NH4 produced per kg SCN consumed)
  // stoichiometric: SCN-(58) + 2H2O -> NH3(17) + CO2 + HS-
  // mass ratio = 17/58 = 0.293
  double nh4_yield;

  // kinetic parameters
  double growth;      // max specific growth rate (s-1)
  double yield;       // biomass yield on SCN- (kg biomass / kg SCN)
  double maintain;    // maintenance rate (s-1)
  double decay;       // decay/death rate (s-1)
};

}

#endif
#endif
