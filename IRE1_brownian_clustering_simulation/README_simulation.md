# IRE1α Brownian Clustering Simulation — Figure 1

Code accompanying:

> **A topology-dependent IRE1α signaling axis generates a cryptic XBP1 isoform present in human cancers**
> Rodrigues Reyes JR, Vergel De Dios C, Datta D, Batjargal T, Pecchia S, Papa FR, Wilson MZ.

---

## Overview

This script implements a coarse-grained Brownian dynamics simulation that models how the subcellular location of IRE1α — ER membrane vs. cytoplasm — determines its oligomerization kinetics and cluster morphology. The key insight is that ER membrane confinement produces anomalous (subdiffusive) motion, which limits cluster coalescence and maintains many small, persistent foci. Cytoplasmic IRE1α undergoes normal 3D Brownian diffusion, leading to rapid coalescence into large hyperclusters within hours.

The simulation reproduces the phenotypic difference between ER-opto-IRE1 and Cyto-opto-IRE1 observed in live-cell imaging experiments and generates the panels shown in **Figure 1** of the paper.

---

## Model Description

Two conditions are simulated in parallel, starting from **N = 200** IRE1 monomers:

### ER-tethered (2D membrane diffusion)
- Particles diffuse on a **12 × 12 μm²** surface representing the ER membrane (z = 0).
- Diffusion coefficient follows the **Saffman–Delbrück** model for membrane proteins:
  `D_mem(n) = D₀ / n`, where `D₀ = 0.005 μm²/s`.
- The ER membrane is partitioned into **corrals** of side length `Lc = 1.2 μm`, representing compartmentalization by ER tubular structures.
- Particles crossing a corral boundary are reflected with probability `1 − p_hop` (`p_hop = 3 × 10⁻⁴`), producing **subdiffusive motion (α < 1)**.
- Outcome: many small, persistent clusters — matching the ER-opto-IRE1 phenotype.

### Cytoplasmic (3D free Brownian motion)
- Particles diffuse freely in a **12 × 12 × 8 μm³** volume.
- Diffusion coefficient follows the **Stokes–Einstein** relation:
  `D_cyto(n) = D₀ / n^(1/3)`, where `D₀ = 0.15 μm²/s` (30× higher than membrane).
- No corral boundaries; only reflecting domain walls.
- Outcome: normal diffusion (α ≈ 1), rapid coalescence into large hyperclusters — matching the Cyto-opto-IRE1 phenotype.

### Cluster merging
Clusters merge irreversibly when their separation falls below the sum of their effective radii:
`r_eff(n) = r₀ · n^(1/3)`, with `r₀ = 0.05 μm` (50 nm interaction range, representing CRY2Olig photobody contact).
Upon merging, the surviving cluster is placed at the mass-weighted centroid.

### MSD analysis
To extract the anomalous diffusion exponent α, short tracking simulations (180 s, 5 s intervals) are run at 1 h and 4 h timepoints — mimicking the experimental imaging protocol. α is obtained by log–log linear regression of MSD vs. lag time τ.

---

## Simulation Parameters

| Parameter | Symbol | ER-tethered | Cytoplasmic |
|---|---|---|---|
| Number of monomers | N | 200 | 200 |
| Domain (x × y) | Lx × Ly | 12 × 12 μm² | 12 × 12 μm² |
| Domain height | Lz | n/a (z = 0) | 8 μm |
| Monomer diffusion coeff. | D₀ | 0.005 μm²/s | 0.15 μm²/s |
| Time step | Δt | 1 s | 1 s |
| Total simulation time | T | 36,000 s (10 h) | 36,000 s (10 h) |
| Base merge radius | r₀ | 0.05 μm | 0.05 μm |
| Corral side length | Lc | 1.2 μm | n/a |
| Hopping probability | p_hop | 3 × 10⁻⁴ | n/a |
| Snapshot interval | — | 120 s | 120 s |

---

## Outputs

Running the script produces:

| File | Description |
|---|---|
| `figure1_simulation.png/.pdf` | Main 8-panel figure (snapshots, MSD curves, α boxplot, cluster count and size over time) |
| `figure_S1_cluster_size_dist.png/.pdf` | Cluster size distributions at 6 timepoints |
| `figure_S2_alpha_vs_size.png/.pdf` | Anomalous diffusion exponent α vs. cluster size |
| `figure_S3_msd_fits.png/.pdf` | Individual MSD curves with power-law fits |
| `figure_S4_diffusion_scaling.png/.pdf` | Diffusion coefficient scaling and cluster growth kinetics |
| `figure_S5_spatial.png/.pdf` | Nearest-neighbor distances and rank-size plots |
| `supplementary_video.mp4` (or `.gif`) | Side-by-side animation of both conditions over 10 hours |

---

## Requirements

```
python >= 3.9
numpy
scipy
matplotlib
```

Install dependencies:
```bash
pip install numpy scipy matplotlib
```

For video export, [FFmpeg](https://ffmpeg.org/) must be installed and accessible on your `PATH`. If FFmpeg is unavailable, the script falls back to a GIF using Pillow.

---

## Usage

1. Set the output directory at the top of `__main__`:
   ```python
   OUT_DIR = "/path/to/your/output/directory"
   ```

2. Run the script:
   ```bash
   python ire1_figure1_simulation.py
   ```

Runtime is approximately **5–15 minutes** on a modern laptop (the ER simulation is the bottleneck due to per-particle corral boundary checks at each of the 36,000 time steps).

---

## Code Structure

| Function | Description |
|---|---|
| `run_clustering()` | Core simulation engine — Brownian dynamics with corral boundaries and irreversible merging |
| `run_msd_tracking()` | Short tracking simulation for MSD and α extraction |
| `make_figure()` | Generates the main 8-panel Figure 1 |
| `make_supp_csd()` | Supplementary Figure S1: cluster size distributions |
| `make_supp_alpha_scatter()` | Supplementary Figure S2: α vs. cluster size scatter |
| `make_supp_msd_fits()` | Supplementary Figure S3: individual MSD fits |
| `make_supp_diffusion_scaling()` | Supplementary Figure S4: diffusion scaling and mass kinetics |
| `make_supp_spatial()` | Supplementary Figure S5: spatial organization |
| `make_video()` | Supplementary video: 3D animation of both conditions |
