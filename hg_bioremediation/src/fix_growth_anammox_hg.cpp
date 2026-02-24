/* ----------------------------------------------------------------------
   NUFEB - Individual-based Modelling for Microbial Communities
   Custom fix: Anammox growth with Hg2+ non-competitive inhibition

   Usage in input script:
     fix ID group nufeb/growth/anammox_hg nh4 K_nh4 o2 K_o2 no2 K_no2 no3 n2 hg K_I &
       growth MU yield Y maintain M decay D

   Example:
     fix growth_ana ANA nufeb/growth/anammox_hg nh4 7e-5 o2 1e-5 no2 5e-5 no3 n2 hg 1e-5 &
       growth 9.26e-7 yield 0.159 maintain 3.5e-8 decay 3e-8

   The Hg2+ inhibition uses non-competitive inhibition:
     mu_eff = mu_max * M_nh4 * M_no2 * I_o2 * I_hg
   where:
     M_nh4 = [NH4]/(K_NH4 + [NH4])          Monod for ammonium
     M_no2 = [NO2]/(K_NO2 + [NO2])          Monod for nitrite
     I_o2  = K_O2/(K_O2 + [O2])             O2 inhibition (anammox is anaerobic)
     I_hg  = K_I/(K_I + [Hg2+])             Hg2+ non-competitive inhibition
------------------------------------------------------------------------- */

#include "fix_growth_anammox_hg.h"

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

FixGrowthAnammoxHg::FixGrowthAnammoxHg(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  // minimum args: fix ID group style nh4 K o2 K no2 K no3 n2 hg K_I = 14
  if (narg < 14)
    error->all(FLERR, "Illegal fix nufeb/growth/anammox_hg command");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/anammox_hg requires grid_style nufeb/chemostat");

  inh4 = -1;
  io2 = -1;
  ino2 = -1;
  ino3 = -1;
  in2 = -1;
  ihg = -1;

  nh4_affinity = 0.0;
  o2_affinity = 0.0;
  no2_affinity = 0.0;
  hg_inhibition = 0.0;

  growth = 0.0;
  yield = 1.0;
  maintain = 0.0;
  decay = 0.0;

  // parse positional substrate arguments: nh4 K_nh4 o2 K_o2 no2 K_no2 no3 n2 hg K_I
  inh4 = grid->find(arg[3]);
  if (inh4 < 0)
    error->all(FLERR, "Can't find substrate name: nh4");
  nh4_affinity = utils::numeric(FLERR,arg[4],true,lmp);

  io2 = grid->find(arg[5]);
  if (io2 < 0)
    error->all(FLERR, "Can't find substrate name: o2");
  o2_affinity = utils::numeric(FLERR,arg[6],true,lmp);

  ino2 = grid->find(arg[7]);
  if (ino2 < 0)
    error->all(FLERR, "Can't find substrate name: no2");
  no2_affinity = utils::numeric(FLERR,arg[8],true,lmp);

  ino3 = grid->find(arg[9]);
  if (ino3 < 0)
    error->all(FLERR, "Can't find substrate name: no3");

  in2 = grid->find(arg[10]);
  if (in2 < 0)
    error->all(FLERR, "Can't find substrate name: n2");

  ihg = grid->find(arg[11]);
  if (ihg < 0)
    error->all(FLERR, "Can't find substrate name: hg");
  hg_inhibition = utils::numeric(FLERR,arg[12],true,lmp);
  if (hg_inhibition <= 0)
    error->all(FLERR, "Hg inhibition constant must be greater than zero");

  // parse keyword arguments
  int iarg = 13;
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
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/anammox_hg command");
    }
  }
}

/* ----------------------------------------------------------------------
   update_cells: compute reaction rates on the Eulerian grid
   called during the chemistry step of NUFEB run_style
------------------------------------------------------------------------- */

void FixGrowthAnammoxHg::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      // Monod kinetics with O2 inhibition and Hg2+ inhibition
      double mu = growth
        * conc[inh4][i] / (nh4_affinity + conc[inh4][i])
        * conc[ino2][i] / (no2_affinity + conc[ino2][i])
        * o2_affinity / (o2_affinity + conc[io2][i])
        * hg_inhibition / (hg_inhibition + conc[ihg][i]);   // <-- Hg inhibition

      // substrate consumption/production stoichiometry (same as standard anammox)
      reac[inh4][i] -= 1.0 / yield * mu * dens[igroup][i];
      reac[ino2][i] -= (1.0 / yield + 1.0 / 1.14) * mu * dens[igroup][i];
      reac[ino3][i] += (1.0 / 1.14) * mu * dens[igroup][i];
      reac[in2][i]  += (2.0 / yield) * mu * dens[igroup][i];

      // Hg2+ is not consumed by anammox -- it inhibits but is not metabolized
      // (Tier 2 could add biosorption here: reac[ihg][i] -= k_sorp * dens * conc[ihg])
    }
  }
}

/* ----------------------------------------------------------------------
   update_atoms: compute per-cell growth rates for the Lagrangian particles
   called during the biology step of NUFEB run_style
------------------------------------------------------------------------- */

void FixGrowthAnammoxHg::update_atoms()
{
  double **conc = grid->conc;

  for (int i = 0; i < grid->ncells; i++) {
    double mu = growth
      * conc[inh4][i] / (nh4_affinity + conc[inh4][i])
      * conc[ino2][i] / (no2_affinity + conc[ino2][i])
      * o2_affinity / (o2_affinity + conc[io2][i])
      * hg_inhibition / (hg_inhibition + conc[ihg][i]);   // <-- Hg inhibition

    double maint = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);

    grid->growth[igroup][i][0] = mu - maint - decay;
  }

  update_atoms_coccus();
}
