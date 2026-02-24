/* ----------------------------------------------------------------------
   NUFEB - Individual-based Modelling for Microbial Communities
   Custom fix: Pseudomonas pseudoalcaligenes CECT5344

   Usage in input script:
     fix ID group nufeb/growth/cect5344 scn K_scn o2 K_o2 nh4 hg K_I &
       growth MU yield Y maintain M decay D hg_uptake KU nh4_yield YN

   Example:
     fix growth_cect CECT nufeb/growth/cect5344 scn 5e-3 o2 2e-4 nh4 hg 1e-2 &
       growth 5.5e-5 yield 0.5 maintain 3e-6 decay 5e-7 &
       hg_uptake 5e-4 nh4_yield 0.293

   Reaction network:
     SCN- + 2H2O --> NH3 + CO2 + HS-     (thiocyanate hydrolase)
     Hg2+(aq) --> Hg2+(intracellular)     (bioaccumulation, mercury sink)

   Growth kinetics:
     mu = mu_max * [SCN]/(K_SCN+[SCN]) * [O2]/(K_O2+[O2])
                 * K_I/(K_I+[Hg2+])

   Parameter basis:
     Hg MIC with cyanide:    200 uM  --> K_I ~ 1e-2 kg/m3
     Hg MIC without cyanide: 10  uM  --> K_I ~ 1e-3 kg/m3
     CN degradation rate:    1.9-2.31 mg CN/L/OD/h at pH 9
     (Ibanez et al. 2023; Luque-Almagro et al. 2010)
------------------------------------------------------------------------- */

#include "fix_growth_cect5344.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "group.h"
#include "grid_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixGrowthCECT5344::FixGrowthCECT5344(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  // minimum: fix ID group style scn K o2 K nh4 hg K_I = 12 positional
  if (narg < 12)
    error->all(FLERR, "Illegal fix nufeb/growth/cect5344 command");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/cect5344 requires grid_style nufeb/chemostat");

  iscn = -1;
  io2 = -1;
  inh4 = -1;
  ihg = -1;

  scn_affinity = 0.0;
  o2_affinity = 0.0;
  hg_inhibition = 0.0;
  hg_uptake = 0.0;
  nh4_yield = 0.293;   // stoichiometric default: 17/58

  growth = 0.0;
  yield = 1.0;
  maintain = 0.0;
  decay = 0.0;

  // parse positional substrate arguments: scn K_scn o2 K_o2 nh4 hg K_I
  iscn = grid->find(arg[3]);
  if (iscn < 0)
    error->all(FLERR, "Can't find substrate name: scn");
  scn_affinity = utils::numeric(FLERR,arg[4],true,lmp);
  if (scn_affinity <= 0)
    error->all(FLERR, "SCN affinity must be greater than zero");

  io2 = grid->find(arg[5]);
  if (io2 < 0)
    error->all(FLERR, "Can't find substrate name: o2");
  o2_affinity = utils::numeric(FLERR,arg[6],true,lmp);
  if (o2_affinity <= 0)
    error->all(FLERR, "O2 affinity must be greater than zero");

  inh4 = grid->find(arg[7]);
  if (inh4 < 0)
    error->all(FLERR, "Can't find substrate name: nh4");

  ihg = grid->find(arg[8]);
  if (ihg < 0)
    error->all(FLERR, "Can't find substrate name: hg");
  hg_inhibition = utils::numeric(FLERR,arg[9],true,lmp);
  if (hg_inhibition <= 0)
    error->all(FLERR, "Hg inhibition constant must be greater than zero");

  // parse keyword arguments
  int iarg = 10;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "maintain") == 0) {
      maintain = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "hg_uptake") == 0) {
      hg_uptake = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "nh4_yield") == 0) {
      nh4_yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/cect5344 command");
    }
  }
}

/* ----------------------------------------------------------------------
   update_cells: compute reaction rates on the Eulerian grid
   called during the chemistry step of NUFEB run_style
------------------------------------------------------------------------- */

void FixGrowthCECT5344::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      // Monod growth on SCN- with O2 and Hg2+ inhibition
      double mu = growth
        * conc[iscn][i] / (scn_affinity + conc[iscn][i])
        * conc[io2][i] / (o2_affinity + conc[io2][i])
        * hg_inhibition / (hg_inhibition + conc[ihg][i]);

      double biomass_rate = mu * dens[igroup][i];

      // SCN- consumption (substrate for growth)
      reac[iscn][i] -= (1.0 / yield) * biomass_rate;

      // NH4+ production from SCN- hydrolysis
      // stoichiometry: each kg SCN- consumed yields nh4_yield kg NH4+
      reac[inh4][i] += nh4_yield * (1.0 / yield) * biomass_rate;

      // Hg2+ intracellular accumulation (mercury sink)
      // first-order uptake proportional to [Hg2+] and CECT5344 biomass density
      // removes free Hg2+ from solution, protecting neighboring cells
      if (hg_uptake > 0.0 && conc[ihg][i] > 0.0) {
        reac[ihg][i] -= hg_uptake * conc[ihg][i] * dens[igroup][i];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   update_atoms: compute per-cell growth rates for the Lagrangian particles
   called during the biology step of NUFEB run_style
------------------------------------------------------------------------- */

void FixGrowthCECT5344::update_atoms()
{
  double **conc = grid->conc;

  for (int i = 0; i < grid->ncells; i++) {
    double mu = growth
      * conc[iscn][i] / (scn_affinity + conc[iscn][i])
      * conc[io2][i] / (o2_affinity + conc[io2][i])
      * hg_inhibition / (hg_inhibition + conc[ihg][i]);

    double maint = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);

    grid->growth[igroup][i][0] = mu - maint - decay;
  }

  update_atoms_coccus();
}
