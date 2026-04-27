# denario-6-superdiffusion-v3

**Scientist:** denario-6
**Date:** 2026-04-27

## Superdiffusion in 3D Chaotic Vortex Filament Systems

This project extends the 2D point-vortex superdiffusion study (superdiffusion_v2) to three dimensions. The physical model consists of N straight vortex filaments in 3D space, whose velocity fields follow the Biot-Savart law (1/r² falloff, Holtsmark distribution). A passive tracer is advected through this field. The goal is to characterise anomalous transport, derive the effective stochastic theory, and compare against 3D isotropic Lévy walk ground truth.

*Key physics change from 2D:* The 3D Biot-Savart kernel (1/r²) differs from the 2D point-vortex kernel (1/r). For a uniform random 3D vortex gas, the theoretical velocity PDF tail is the Holtsmark distribution P(v) ~ v^{-5/2} (stability index α_stable = 3/2), predicting anomalous exponent α_diff = 3/2 = 1.5 — different from the 2D Cauchy prediction.

---

### File 1: 3D Vortex Filament Tracer Simulations

*Path:* `/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy`

NumPy structured array, 10,000 rows (20 tracers × 500 time steps).

| Field          | dtype    | Description                                                                 |
|----------------|----------|-----------------------------------------------------------------------------|
| trajectory_id  | int32    | Tracer ID 1–20 (5 per N configuration)                                      |
| n_filaments    | int32    | Number of vortex filaments N ∈ {5, 10, 20, 40}                              |
| gamma_std      | float64  | Circulation std dev γ_std=1.0; Γ_i ~ N(0, γ_std²)                          |
| time           | float64  | t = step × dt, dt=0.05 s, range [0, 24.75] s                               |
| x_noisy        | float64  | Tracer x with Gaussian noise σ=0.02 m                                       |
| y_noisy        | float64  | Tracer y with Gaussian noise σ=0.02 m                                       |
| z_noisy        | float64  | Tracer z with Gaussian noise σ=0.02 m                                       |
| x_true         | float64  | Noise-free tracer x                                                          |
| y_true         | float64  | Noise-free tracer y                                                          |
| z_true         | float64  | Noise-free tracer z                                                          |
| msd_true       | float64  | True 3D squared displacement from origin: (x²+y²+z²)                       |

*Physics model:* N infinite straight vortex filaments with random positions (uniform in [-10,10]³) and random unit-vector orientations. Circulations Γ_i ~ N(0, 1). Velocity at tracer position r induced by filament i:
  v_i = (Γ_i / 2π) * (d̂_i × r_⊥_i) / |r_⊥_i|²
where r_⊥_i is the perpendicular vector from the nearest point on filament i to the tracer. Tracer advected by RK4. Filament positions and orientations are fixed (quenched disorder — no filament-filament dynamics for tractability).

*Loading:*
```python
import numpy as np
vf = np.load('/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy', allow_pickle=False)
mask = vf['n_filaments'] == 40
traj = vf[mask]  # 2500 rows, N=40 configuration
```

---

### File 2: 3D Isotropic Lévy Walk Ground Truth

*Path:* `/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy`

NumPy structured array, 24,000 rows (40 trajectories × 600 time steps).

| Field          | dtype    | Description                                                                 |
|----------------|----------|-----------------------------------------------------------------------------|
| trajectory_id  | int32    | Trajectory ID 1–40 (10 per β class)                                         |
| beta           | float64  | Pareto tail index β; flight-time PDF P(τ) ~ τ^{-(1+β)}                      |
| alpha_theory   | float64  | Theoretical α = 3−β (for 1 < β < 2)                                        |
| time           | float64  | t = step × dt, dt=0.1 s, range [0, 59.9] s                                 |
| x_noisy/y_noisy/z_noisy | float64 | Noisy positions, σ=0.05 m per coordinate                       |
| x_true/y_true/z_true    | float64 | Noise-free positions                                                |
| msd_true       | float64  | True 3D squared displacement from origin                                    |

*Model:* Isotropic 3D Lévy walk. Flight directions uniform on the unit sphere (φ uniform in [0,2π], cosθ uniform in [-1,1]). Flight durations Pareto(β), minimum = dt. Speed v₀=1.0 m/s.

*Classes (10 trajectories each):*
| β    | α_theory | Regime                  |
|------|----------|-------------------------|
| 1.2  | 1.8      | Strong superdiffusion   |
| 1.5  | 1.5      | Moderate superdiffusion |
| 1.8  | 1.2      | Mild superdiffusion     |
| 2.5  | 1.0      | Effectively normal      |

*Loading:*
```python
lw = np.load('/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy', allow_pickle=False)
mask = lw['beta'] == 1.5
traj = lw[mask]  # 6000 rows, β=1.5 (α=1.5)
```

---

### Key differences from the 2D study (superdiffusion_v2)

1. *3D Biot-Savart kernel:* velocity falloff ~ 1/r² (vs 1/r in 2D). Theoretical Holtsmark distribution for vortex velocity PDF has α_stable=3/2.
2. *Quenched disorder:* filament positions are fixed (not evolving), so vortex-vortex chaos is absent. The tracer still undergoes anomalous transport due to the heterogeneous velocity field.
3. *Isotropic 3D Lévy walks:* flight directions uniform on sphere, MSD = x²+y²+z².
4. *All MSD analyses should use:* MSD(t) = ⟨x²+y²+z²⟩ (3D displacement).

### Suggested analyses

1. *3D MSD scaling:* TAMSD with 3D displacement. Fit α(N) across vortex configurations. Compare to Lévy walk ground truth.
2. *3D VACF:* full 3D velocity vector v=(vx,vy,vz); C_v(Δt) = ⟨v(t)·v(t+Δt)⟩/⟨|v|²⟩.
3. *Residence time distributions:* same methodology as 2D, using 3D speed threshold.
4. *3D displacement PDFs:* marginal PDF of one coordinate Δx at multiple lag times. Test for Holtsmark vs Lévy-stable vs Gaussian.
5. *Ergodicity breaking:* EB parameter from 3D TAMSD.
6. *Comparison with 2D results:* how do α(N) and the VACF zero-crossing times differ between 2D and 3D? Does the Holtsmark prediction (α≈1.5) match the N=40 3D result better than the 2D Cauchy prediction (α→1)?
7. *Regime map:* same (EB, α) plane as in v2, with v2 results overlaid for direct 2D vs 3D comparison.
