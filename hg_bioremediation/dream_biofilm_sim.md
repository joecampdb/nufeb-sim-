# The Dream Biofilm Simulation

## A Technical Design for Maximizing Both Beauty and Scientific Impact in NUFEB-2

---

## 1. Motivation: The Optimization Problem

Biofilm simulations face a fundamental tension. The parameters that produce the
most visually striking morphologies (sparse seeding, strong diffusion limitation,
low shear, heterogeneous species) also tend to produce the most scientifically
informative results -- but only up to a point. Beyond that point, computational
cost explodes, numerical stability degrades, and the simulation becomes an
expensive screensaver rather than a research tool.

This document explores the design space for a "dream" NUFEB-2 simulation that
simultaneously maximizes:

- **Morphological complexity** (finger structures, mushroom caps, channels,
  stratified layers, EPS scaffolding)
- **Visual beauty** (species contrast, structural drama, translucent EPS matrix)
- **Scientific utility** (testable hypotheses about mercury bioremediation,
  consortia design, and biofilm ecology)
- **Computational tractability** (completable on a workstation in < 48 hours)

The key insight from the literature is that these objectives are not orthogonal.
The most scientifically interesting biofilms *are* the most morphologically
complex ones, because structure determines function. The challenge is finding the
parameter regime where complexity is maximized without computational blowup.

---

## 2. Lessons from the Mercury Reef

The ongoing Mercury Reef simulation (400x400x800um, 7 atom types, 4 seed
clusters, 365-day run) has revealed several critical dynamics:

| Observation | Implication |
|-------------|-------------|
| DEAD cells reach 90% of total atoms by day 70 | Dead cell accumulation is the primary computational bottleneck |
| HET dominates the living population (~100k steady state) | Heterotroph carrying capacity is reached early; further growth is zero-sum |
| NOB went extinct | The nitrogen cycle is incomplete without NOB; community is unstable |
| ANA collapsed to 12 cells | Hg2+ inhibition is too strong at 50 ug/L; CECT cannot rescue fast enough |
| CECT expanded 15x between day 70 and day 124 | CECT growth lag allows Hg damage before remediation kicks in |
| EPS peaked at 161k then declined to ~40k | EPS decay (2e-6) exceeds secretion at steady state; scaffolding is dissolving |
| 4-lobe morphology with CECT eruption zones | Asymmetric seeding produces scientifically meaningful spatial heterogeneity |
| Steps slow from 1.5s to 25s as atoms pass 1.5M | Neighbor list and pair calculations scale with dead cell count |

### Critical Design Failures to Correct

1. **No dead cell management.** The simulation spends 80%+ of compute on inert
   particles. A periodic dead cell culling strategy (or the `nufeb/merge_eps`
   fix repurposed for dead cells) would reduce atom count by 5-10x.

2. **EPS net loss at steady state.** With `epsyield 0.18` and `eps decay 2e-6`,
   the scaffolding dissolves faster than it forms once HET growth slows. The
   dream simulation needs either higher epsyield, lower eps decay, or both.

3. **ANA rescue too slow.** CECT needs ~60 days to build sufficient population
   to meaningfully draw down Hg2+. By then, ANA is functionally extinct. Either
   start with more CECT, lower initial Hg2+, or increase hg_uptake rate.

4. **NOB extinction.** Without NOB, the nitrification cascade breaks. NOB needs
   a protected niche (e.g., co-located with AOB in a dedicated seed cluster with
   favorable NO2- supply).

---

## 3. Technical Design: The Dream Simulation

### 3.1 Domain and Grid

| Parameter | Mercury Reef | Dream Setup | Rationale |
|-----------|-------------|-------------|-----------|
| Domain | 400x400x800um | 600x600x1000um | Wider for more colony separation; taller for vertical stratification |
| Grid spacing | 12.5um (32x32x64) | 20um (30x30x50) | Coarser grid: 45k cells vs 65k cells; faster diffusion solver |
| Boundary (z) | ff (fixed-fixed) | ff | Non-periodic top for bulk liquid interface |
| Boundary (xy) | pp (periodic) | pp | Periodic lateral boundaries avoid edge effects |

The coarser grid (20um) sacrifices some diffusion resolution but keeps the solver
tractable. At 30x30x50 = 45,000 grid cells, each diffusion iteration is ~30%
faster than the current 65,536-cell grid. The wider domain (600um) provides room
for 8-12 seed clusters with sufficient separation for independent colony
morphology before merger.

```lammps
region simu_domain block 0.0 6e-4 0.0 6e-4 0.0 1e-3 units box
create_box 8 simu_domain

grid_style nufeb/chemostat 8 sub nh4 o2 no2 no3 n2 hg scn 2e-5
```

### 3.2 Species Palette (8 Atom Types)

| Type | Species | Role | Color | mu_max (s-1) | Key Trait |
|------|---------|------|-------|-------------|-----------|
| 1 | HET | Heterotroph, EPS producer | Tan | 6.94e-5 | Structural backbone |
| 2 | AOB | NH4+ -> NO2- | Dark Olive Green | 2.37e-5 | Nitrification step 1 |
| 3 | NOB | NO2- -> NO3- | Sienna | 1.68e-5 | Nitrification step 2 |
| 4 | ANA | Anammox (Hg-inhibited) | Gold | 9.26e-7 | N2 production, Hg "prey" |
| 5 | CECT | SCN- degrader, Hg sink | Teal | 2.50e-5 | Mercury "protector" |
| 6 | CYANO | Cyanobacteria (phototroph) | Chartreuse | ~1.5e-5 | O2 producer at surface, visual crown |
| 7 | EPS | Extracellular matrix | Wheat | -- | Structural scaffolding |
| 8 | DEAD | Lysed cells | Dark Gray | -- | Necrotic core |

**New addition: CYANO (Cyanobacteria).** NUFEB has a built-in
`nufeb/growth/cyano` fix for phototrophs that grow on light + CO2 and produce O2
and sucrose. Adding cyanobacteria to the surface creates:

- A visually distinct **chartreuse crown** on top of the biofilm (phototrophs
  must be at the surface for light access)
- An **endogenous O2 source** that intensifies the oxic/anoxic gradient within
  the biofilm, driving stronger vertical stratification
- A **sucrose cross-feeding** opportunity: CYANO exports sucrose, which could
  support HET growth in a syntrophic loop (implementable via `nufeb/growth/ecoli`
  or a custom Monod fix)
- A **color contrast anchor**: chartreuse at the top, gold ANA and teal CECT in
  the deep interior, darkgray DEAD at the core

The phototroph layer also has direct research significance: cyanobacterial mats
are the dominant mercury-cycling biofilms in freshwater periphyton systems, where
the oxic/anoxic boundary determines the balance between Hg(II) reduction
(detoxification) and Hg methylation (toxification).

### 3.3 Seeding Strategy: 12 Scattered Micro-Patches

Instead of 4 corner clusters, use 12 small seed regions scattered across the
substratum with different species compositions. This creates:

- Multiple independent colony morphologies before merger
- Competition boundaries where colonies collide (the most structurally complex
  zones in real biofilms)
- Spatially heterogeneous mercury remediation (some zones protected by CECT,
  others exposed)

```lammps
lattice sc 1e-5 origin 0 0 0

# 12 seed patches (4x3 grid with jitter, each ~40um)
region s01 block  4  8  4  8 0 1    # HET + CYANO
region s02 block 14 18  4  8 0 1    # AOB + NOB + HET
region s03 block 24 28  4  8 0 1    # CECT + HET
region s04 block 34 38  4  8 0 1    # ANA + CECT (rescue pair)
region s05 block  4  8 14 18 0 1    # HET + AOB
region s06 block 14 18 14 18 0 1    # CECT + ANA + HET (core rescue)
region s07 block 24 28 14 18 0 1    # HET + NOB + CYANO
region s08 block 34 38 14 18 0 1    # AOB + NOB
region s09 block  4  8 24 28 0 1    # HET + CECT + CYANO
region s10 block 14 18 24 28 0 1    # ANA + HET + AOB
region s11 block 24 28 24 28 0 1    # CECT + HET + NOB
region s12 block 34 38 24 28 0 1    # HET + CYANO + CECT
```

Each patch gets 10-20 cells, with composition chosen to create functional
diversity:

- **Rescue pairs** (s04, s06): ANA + CECT co-seeded so CECT can locally
  protect ANA from Hg2+ before the global concentration drops
- **Nitrifier hubs** (s02, s08): AOB + NOB co-seeded to establish the
  nitrification cascade before either species goes extinct
- **Phototroph anchors** (s01, s07, s09, s12): CYANO seeded at 4 positions
  to ensure surface colonization from multiple origins
- **Remediation nodes** (s03, s11): CECT-heavy patches as mercury sinks

### 3.4 EPS System (Strengthened)

The Mercury Reef showed EPS declining at steady state. The dream simulation
adjusts parameters to maintain a persistent, visually prominent EPS matrix:

| Parameter | Mercury Reef | Dream | Rationale |
|-----------|-------------|-------|-----------|
| epsyield | 0.18 | 0.25 | More EPS per unit HET growth |
| epsdens | 30 | 30 | Standard |
| eps_secretion ratio | 1.3 | 1.2 | Earlier secretion = more frequent EPS bead ejection |
| eps decay | 2e-6 | 8e-7 | Slower dissolution preserves scaffolding |
| adhesion ke | 5e2 | 1e3 | Stronger EPS-mediated cohesion for structural integrity |
| wall adhesion kn | 1e3 | 2e3 | Firmer substratum attachment |

Additionally, **EPS merging** (`nufeb/merge_eps`) should be added to prevent EPS
bead count explosion. Small EPS particles below 1.2um merge with their nearest
neighbor, maintaining visual fidelity while halving EPS atom count:

```lammps
fix merge_eps EPS nufeb/merge_eps 1.2e-6 9876 epsdens 30
```

### 3.5 Shear: Carving Channels and Streamers

The single most impactful addition for visual complexity is `fix nufeb/shear`.
Without shear, biofilms grow into smooth mounds. With shear, they develop:

- **Water channels** carved between microcolonies
- **Streamers** extending downstream from colony tips
- **Asymmetric morphology** with upstream-facing flat profiles and downstream
  towers/mushrooms
- **Detachment events** that expose internal structure

NUFEB's shear fix applies a linear Couette flow profile:

```lammps
# Gentle shear: 0.05 s-1 in +x direction, water viscosity
# Layer height = 500um (half-domain, so deep biofilm sees less shear)
fix shear1 all nufeb/shear 0.05 1e-3 +x layer 5e-4
```

The shear rate of 0.05 s-1 is deliberately low. Literature shows that high shear
(> 0.5 s-1) flattens biofilms into featureless mats, while the intermediate
regime (0.01 - 0.1 s-1) maximizes morphological complexity by selectively
removing weakly attached cells while allowing towers and mushrooms to persist.

**Research hypothesis 1:** *Shear-induced channel formation enhances mercury
bioremediation efficiency by increasing the effective surface area between the
biofilm and the bulk liquid, allowing deeper Hg2+ penetration to interior CECT
colonies.*

### 3.6 Dead Cell Management

The Mercury Reef's computational bottleneck is its 1.35 million dead cells. The
dream simulation implements a periodic dead cell culling strategy:

**Option A: Buried dead cell deletion.** Every N biological timesteps, delete
DEAD atoms that are more than 3 cell diameters below the biofilm surface and
completely surrounded by other dead cells. These particles contribute nothing to
structure or visual appearance.

This requires a custom compute or a periodic `delete_atoms` command. In NUFEB,
the simplest approach is a region-based deletion of dead cells below a z-threshold
that tracks the biofilm base:

```lammps
# Delete deeply buried dead cells every 240 steps (10 days)
# Keep only the top 100um of dead cells for visual necrotic core
variable zsurf equal 1e-4   # adjust dynamically or use fixed threshold
region dead_cull block INF INF INF INF 0.0 ${zsurf} units box
# (executed via a run-loop or post-processing script)
```

**Option B: `nufeb/merge_eps` adapted for DEAD.** Merge neighboring dead cells
into larger inert aggregates, reducing particle count while preserving visual
bulk. This is less invasive than deletion.

**Target:** Keep total atom count below 500,000 at all times. This ensures
steps remain under 10 seconds even at day 365.

### 3.7 Substrate Chemistry

| Substrate | Initial (kg/m3) | BC (z-top) | Role |
|-----------|-----------------|------------|------|
| sub | 1.2e-3 | Dirichlet | Organic carbon for HET |
| nh4 | 1e-3 | Dirichlet | Ammonium for AOB/ANA |
| o2 | 0.5e-3 | Dirichlet | Dissolved oxygen |
| no2 | 1e-8 | Dirichlet | Nitrite (AOB product, NOB/ANA substrate) |
| no3 | 1e-8 | Dirichlet | Nitrate (NOB product) |
| n2 | 1e-8 | Dirichlet | Dinitrogen (ANA product) |
| hg | 2.5e-5 | Dirichlet | Hg2+ -- 25 ug/L (halved from Mercury Reef) |
| scn | 1.5e-3 | Dirichlet | Thiocyanate for CECT |

**Key change: Hg2+ reduced to 25 ug/L** (from 50 ug/L). The Mercury Reef showed
that 50 ug/L suppresses ANA to functional extinction before CECT can respond.
At 25 ug/L with the same K_I = 10 ug/L, ANA inhibition is still strong (71%)
but survivable, creating a dynamic tension where CECT activity gradually shifts
the balance.

**Research hypothesis 2:** *There exists a critical Hg2+ concentration threshold
below which anammox bacteria can persist long enough for CECT-mediated
remediation to rescue the population, and above which ANA extinction is
irreversible regardless of CECT density.*

*Note: If CYANO is included, additional substrates (light, co2, sucrose, gco2)
would be needed. This adds complexity to the grid solver. An alternative is to
approximate CYANO as a simple O2-producing autotroph using `nufeb/growth/monod`
on a generic "light" substrate with fixed Dirichlet BC at the top surface.*

### 3.8 Diffusion Coefficients

| Substrate | D_bulk (m2/s) | Biofilm ratio | D_biofilm (m2/s) |
|-----------|--------------|---------------|------------------|
| sub | 1.16e-9 | 0.75 | 8.7e-10 |
| nh4 | 1.97e-9 | 0.75 | 1.48e-9 |
| o2 | 2.30e-9 | 0.75 | 1.73e-9 |
| no2 | 1.85e-9 | 0.75 | 1.39e-9 |
| no3 | 1.85e-9 | 0.75 | 1.39e-9 |
| n2 | 2.30e-9 | 0.75 | 1.73e-9 |
| hg | 9.10e-10 | 0.40 | 3.64e-10 |
| scn | 1.85e-9 | 0.60 | 1.11e-9 |

**Key change: Hg diffusion ratio lowered to 0.40** (from 0.50). Mercury ions
bind strongly to EPS and cell surfaces, further retarding diffusion within the
biofilm matrix. This creates steeper Hg2+ gradients, which:

- Makes CECT's local Hg depletion zones more pronounced and visually mappable
  in VTK concentration fields
- Creates stronger spatial protection for ANA cells near CECT colonies
- Is more physically realistic (Hg2+ has strong affinity for thiol groups in EPS)

**Research hypothesis 3:** *Reduced Hg2+ diffusivity within the EPS matrix
creates localized "mercury shadow" zones around CECT colonies where ANA can
survive, leading to emergent spatial co-localization of these two species -- a
form of niche construction mediated by EPS.*

---

## 4. Hypotheses at the Beauty-Utility Frontier

The dream simulation is designed to test hypotheses that are simultaneously
visually dramatic and scientifically novel. Each hypothesis predicts an emergent
spatial pattern that would be visible in both LAMMPS renders and Blender
visualizations.

### H1: The Mercury Shadow Effect

**Statement:** CECT colonies create local Hg2+ depletion zones ("mercury
shadows") within which ANA can survive, leading to emergent spatial
co-localization of the Hg-sensitive and Hg-accumulating species.

**Visual signature:** In VTK concentration fields, Hg2+ depletion zones appear
as dark halos around teal CECT clusters. Gold ANA cells should preferentially
survive within these halos, creating teal-gold spatial pairs against a tan HET
background.

**Scientific significance:** This would demonstrate a form of biologically
mediated niche construction -- one species altering the local chemical
environment to create habitat for another. This has direct implications for
designing bioremediation consortia: spatial arrangement matters as much as species
composition.

**Testable in simulation:** Compare ANA survival rates as a function of distance
from nearest CECT colony centroid. If the mercury shadow hypothesis holds, ANA
survival should decay exponentially with distance from CECT.

### H2: Shear-Enhanced Remediation

**Statement:** Gentle shear flow carves water channels through the biofilm that
increase the effective contact surface between bulk liquid and interior CECT
colonies, enhancing overall Hg2+ uptake rate compared to a quiescent (no-shear)
biofilm of equal biomass.

**Visual signature:** Channels visible as voids in the biofilm structure, with
CECT colonies lining channel walls (analogous to vascular endothelium). The
biofilm transitions from a smooth mound to a ridged, channelized architecture.

**Scientific significance:** This connects biofilm structural mechanics to
bioremediation performance -- a largely unexplored link in the mercury literature.
If confirmed, it suggests that bioreactor design should optimize for intermediate
shear rather than maximizing biofilm density.

**Testable in simulation:** Run two identical simulations with and without
`fix nufeb/shear`. Compare total Hg2+ removal at day 180. Measure effective
biofilm surface area in both cases using VTK isosurface analysis.

### H3: The Phototroph Crown

**Statement:** Surface-colonizing cyanobacteria create an O2-producing "crown"
that intensifies the oxic/anoxic gradient within the biofilm, pushing the anoxic
zone (where ANA operates) deeper into the biofilm interior and creating a thicker
aerobic barrier that retards Hg2+ diffusion to ANA.

**Visual signature:** A distinct chartreuse cap on the biofilm surface, with
oxygen concentration fields showing a bright oxic zone in the top 50-100um and
deep anoxia below. ANA (gold) cells stratified into a thin band at the
oxic/anoxic interface.

**Scientific significance:** Cyanobacterial mats are the dominant periphyton
biofilm in freshwater systems where mercury cycling occurs. The interplay between
photosynthetic O2 production and mercury methylation/demethylation is a central
open question in environmental mercury science. Simulating this interface
directly addresses the "periphyton paradox" identified in recent reviews.

**Testable in simulation:** Compare O2 penetration depth and ANA vertical
position with and without CYANO. Measure whether CYANO's O2 production shifts
the ANA population deeper or eliminates it entirely (O2 inhibition of anammox).

### H4: Nutrient-Complexity Sweet Spot

**Statement:** There exists an intermediate nutrient supply regime that maximizes
biofilm morphological complexity (quantified by fractal dimension, roughness
coefficient, and surface-to-volume ratio). Too much nutrient produces flat,
dense mats; too little produces sparse, fragmented clusters.

**Visual signature:** At the optimal nutrient level, the biofilm exhibits
simultaneous mushroom caps, finger protrusions, internal voids, and surface
rugosity -- the full repertoire of biofilm architecture.

**Scientific significance:** This has been demonstrated computationally for
single-species biofilms (Royal Society Interface 2023), but never for
multi-species mercury-cycling consortia. Demonstrating that the
nutrient-complexity relationship holds for functional consortia would bridge
the gap between theoretical biofilm ecology and practical bioreactor design.

**Testable in simulation:** Run 3 simulations at sub = {0.4e-3, 1.2e-3, 3.6e-3}
(low/medium/high carbon). Quantify morphological complexity from VTK outputs
using surface roughness, height variance, and fractal dimension metrics.

### H5: Emergent Division of Labor Through Spatial Segregation

**Statement:** In a multi-species biofilm with metabolic cross-feeding (AOB
produces NO2- for NOB and ANA; CECT produces NH4+ for AOB; HET produces EPS
for all), the species self-organize into spatially distinct functional layers
without any explicit positioning -- an emergent division of labor driven solely
by substrate gradients and growth kinetics.

**Visual signature:** Vertical stratification visible in cross-sections: CYANO
(chartreuse) at the surface, HET (tan) in the upper aerobic zone, AOB (olive)
and NOB (sienna) in the microaerobic transition, ANA (gold) and CECT (teal) in
the deep anoxic zone, DEAD (gray) at the core.

**Scientific significance:** Self-organized spatial structure in multi-species
biofilms is theoretically predicted but rarely demonstrated in silico with more
than 3 species. A 6-species (+ EPS + DEAD) simulation showing clear vertical
zonation would be a significant contribution to IbM biofilm modeling.

**Testable in simulation:** Compute species-specific mean z-position as a
function of time. If spatial segregation is emergent, species z-positions should
separate monotonically over the first 60-90 days and then stabilize.

---

## 5. Computational Optimization Strategy

### 5.1 Dead Cell Budget

The single most important optimization. Target: **< 500k total atoms at all
times.**

| Strategy | Implementation | Atom Reduction |
|----------|---------------|---------------|
| Buried dead cell deletion | Periodic `delete_atoms` in z < biofilm_base | ~60-80% of DEAD |
| Dead cell merging | `nufeb/merge_eps`-style fix for DEAD | ~50% of DEAD |
| Higher death diameter | Increase from 2um to 2.5um | Fewer cells die (they divide instead) |
| Lower decay rates | Reduce maintenance costs | Fewer cells starve to death |

### 5.2 Grid Optimization

| Strategy | Implementation | Speedup |
|----------|---------------|---------|
| Coarser grid (20um) | `grid_style ... 2e-5` | ~30% per diffusion step |
| Lower diffmax | `diffmax 500` (from 1000) | Caps worst-case diffusion iterations |
| Higher difftol | `difftol 1e-5` (from 1e-6) | Fewer iterations to converge |

### 5.3 Output Optimization

| Output | Mercury Reef | Dream | Rationale |
|--------|-------------|-------|-----------|
| PNG every | 40 steps | 60 steps | Fewer renders, still smooth video |
| VTK every | 80 steps | 120 steps | Sufficient for Blender animation |
| Thermo every | 1 step | 1 step | Keep for monitoring |

### 5.4 Estimated Run Time

With dead cell management keeping atoms under 500k:

- Steps: 6480 (270 days)
- Estimated time/step: 5-10 seconds (with <500k atoms and coarser grid)
- **Total: 9-18 hours** (vs. 40-55+ hours for Mercury Reef at 365 days)

---

## 6. Blender Rendering: Where Beauty Meets Data

The dream simulation's VTK output is specifically designed for the Blender Cycles
pipeline described previously. Key rendering considerations:

### 6.1 Material Assignments

| Type | Blender Material | Key Properties |
|------|-----------------|----------------|
| HET (Tan) | Principled BSDF | SSS 0.1, clearcoat 0.15 (wet membrane) |
| AOB (Olive) | Principled BSDF | SSS 0.15, green-biased scattering |
| NOB (Sienna) | Principled BSDF | SSS 0.1, warm reddish scattering |
| ANA (Gold) | Principled BSDF | SSS 0.2, high clearcoat (the "prey" glows) |
| CECT (Teal) | Principled BSDF | SSS 0.15, cool scattering (the "protector") |
| CYANO (Chartreuse) | Principled BSDF | SSS 0.2, slight emission 0.01 (photosynthetic glow) |
| EPS (Wheat) | Principled BSDF | **Full SSS + transmission 0.4, IOR 1.33** (translucent gel) |
| DEAD (Dark Gray) | Principled BSDF | Matte roughness 0.85, no SSS (lifeless) |

### 6.2 Concentration Field Overlays

The `.vti` grid data enables volumetric rendering of substrate fields:

- **Hg2+ concentration:** Render as a red-to-transparent volume overlay.
  Mercury shadows around CECT colonies become visible as clear zones in a red
  haze -- a direct visualization of Hypothesis H1.
- **O2 concentration:** Blue-to-transparent gradient showing the oxic/anoxic
  boundary. The CYANO crown's O2 production is visible as a bright blue cap.
- **Reaction rate fields:** Growth rate isosurfaces show where active metabolism
  is occurring vs. dormant zones.

### 6.3 Camera and Lighting

- **HDRI:** Laboratory or underwater HDRI from Poly Haven for neutral,
  scientific lighting with subtle caustics through EPS
- **Depth of field:** f/4.0 focused on the biofilm-liquid interface, blurring
  the substratum and upper bulk liquid
- **Camera path:** Slow 360-degree orbit for the animation, with a dolly zoom
  into a channel cross-section at the midpoint
- **Cross-section reveal:** A clipping plane sweeps through the biofilm at
  day 270 to expose internal stratification

---

## 7. Validation Pathway

The dream simulation's predictions can be validated against published
experimental data:

| Hypothesis | Validation Data Source |
|------------|----------------------|
| H1 (Mercury shadow) | Fluorescence microscopy of Hg-reporter biofilms; Hg microelectrode profiling |
| H2 (Shear-enhanced remediation) | Packed-bed bioreactor Hg removal efficiency at varying flow rates (AEM 2002, spatially oscillating Hg activity) |
| H3 (Phototroph crown) | Periphyton O2 microsensor profiling; confocal imaging of cyanobacterial mats |
| H4 (Nutrient-complexity) | COMSTAT/BiofilmQ analysis of confocal z-stacks at varying C/N ratios |
| H5 (Emergent stratification) | FISH (fluorescence in situ hybridization) cross-sections of multi-species biofilms |

### Morphological Metrics (Extractable from VTK)

Standard biofilm morphology quantification metrics applicable to our outputs:

- **Biovolume:** total volume occupied by living cells
- **Mean/max thickness:** biofilm height profile statistics
- **Roughness coefficient (Ra):** standard deviation of height; higher Ra =
  more complex surface topology
- **Surface area to biovolume ratio:** higher = more channelized, more
  gas/liquid exchange surface
- **Fractal dimension:** 2.0 = flat film, ~2.6 = maximally complex, <2.0 =
  sparse clusters
- **Species-specific z-profiles:** vertical distribution of each species type
- **Nearest-neighbor species co-occurrence:** spatial correlation between
  ANA and CECT (H1 test)

---

## 8. Summary: The Pareto Front

The dream simulation sits at the intersection of four optimization axes:

```
                    Beauty
                      |
                      |     * Dream Simulation
                      |    /
                      |   /
                      |  /
                      | /
 Tractability --------+----------- Complexity
                      |
                      |
                      |
                      |
                    Utility
```

**Beauty** is maximized by: 8 species with distinct colors, translucent EPS,
shear-carved channels, phototroph crown, Blender Cycles rendering with SSS and
volumetric Hg overlays.

**Complexity** is maximized by: 12 scattered seed patches, gentle shear flow,
intermediate nutrient supply, strong diffusion limitation for Hg2+, multiple
metabolic cross-feeding loops.

**Utility** is maximized by: 5 testable hypotheses (mercury shadow, shear
enhancement, phototroph crown, nutrient-complexity sweet spot, emergent
stratification), direct connection to open questions in mercury bioremediation
science, validation pathway against experimental techniques.

**Tractability** is maintained by: dead cell management (<500k atoms), coarser
grid (20um), EPS merging, 270-day run (sufficient for mature morphology), output
every 60-120 steps.

The key insight is that these objectives are largely **aligned, not opposed**.
The most scientifically interesting biofilm (one with emergent spatial structure
driven by metabolic cross-feeding and environmental gradients) is also the most
visually striking one. The only true trade-off is computational cost, and that
is addressed by engineering the simulation to be efficient (dead cell management,
grid optimization) rather than by sacrificing biological realism.

---

## References

- Li, B. et al. (2019). NUFEB: A massively parallel simulator for individual-based
  modelling of microbial communities. *PLOS Computational Biology.*
- Buttner, L. et al. (2024). Is it selfish to be filamentous in biofilms?
  *PLOS Computational Biology.*
- Drescher, K. et al. (2022). Extracellular DNA in biofilm streamers. *PNAS.*
- Picioreanu, C. et al. (2004). Particle-based multidimensional multispecies
  biofilm model. *Applied and Environmental Microbiology.*
- Royal Society Interface (2023). Nutrient-driven biofilm stratification and
  fractal dimension optimization.
- AEM (2002). Spatially oscillating mercury-reducing activity in packed-bed
  bioreactors.
- Frontiers in Microbiology (2023). Methylmercury formation in *Geobacter*
  biofilms.
- Frontiers in Environmental Science (2025). Periphyton-mercury interactions:
  the periphyton paradox.
- Nature Communications (2024). Engineering microbiomes for bioremediation.
- ISME Communications (2025). Ecological design of synthetic microbial communities.
- PMC (2022). Gradients and consequences of heterogeneity in biofilms.
- PMC (2022). CFD-DEM coupled EPS deformation modeling.
- PMC (2024). Global priority questions in biofilm research.
