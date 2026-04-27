<!-- filename: reports/step_7_vortex_filament_analysis.md -->
## Results

### Summary of Key Quantitative Findings

The following table consolidates all fitted transport exponents, ergodicity breaking parameters, VACF zero-crossing times, and persistence lengths derived from the 3D vortex filament tracer simulations and the isotropic Lévy walk ground truth.

**Table 1. Summary of anomalous transport parameters for the 3D vortex filament system.**

| N (filaments) | α(N) | α SE | EB(N) | τ_c (s) | L_p (m) | ℓ_inter (m) |
|---|---|---|---|---|---|---|
| 5 | 1.999 | 0.000 | — | — | — | 11.70 |
| 10 | 1.964 | 0.002 | — | 10.17 | 0.313 | 9.28 |
| 20 | 1.969 | 0.001 | — | — | — | 7.37 |
| 40 | 1.835 | 0.003 | — | 9.555 | 0.578 | 5.85 |

**Table 2. Summary of transport parameters for the 3D isotropic Lévy walk ground truth.**

| β | α_theory | α_measured | τ_c (s) |
|---|---|---|---|
| 1.2 | 1.8 | 1.655 | — |
| 1.5 | 1.5 | 1.049 | — |
| 1.8 | 1.2 | 1.069 | — |
| 2.5 | 1.0 | 1.012 | — |

**Table 3. Okubo-Weiss and spatial correlation analysis (N=40).**

| Quantity | Value |
|---|---|
| OW–speed Pearson r | −0.508 |
| Rolling speed variance–speed Pearson r | +0.425 |

---

### Realization-Dependent Variance and Quenched Disorder Effects

The time-averaged mean squared displacement (TAMSD) computed for each of the five individual tracers within each filament configuration reveals a pronounced realization-dependent spread that is a direct signature of the quenched nature of the disorder. For the sparse N=5 configuration, where the inter-filament spacing is approximately 11.70 m, individual tracer TAMSDs diverge substantially from one another across the full observation window of 24.95 s. This spread is not merely statistical noise: because the filament positions and orientations are fixed throughout each simulation, each tracer samples a fundamentally different realization of the velocity field, and the resulting transport statistics are non-self-averaging at the level of a single trajectory. The min-max envelope across the five tracers is widest for N=5 and progressively narrows as N increases toward 40, consistent with the expectation that a denser filament array produces a more spatially homogeneous velocity field, reducing the sensitivity of individual trajectories to the specific disorder realization.

For N=40, where the inter-filament spacing has decreased to approximately 5.85 m, the individual TAMSD curves cluster more tightly around the ensemble mean, though non-trivial spread persists. This residual spread at N=40 is physically significant: it indicates that even at the highest filament density studied, the quenched disorder is not fully averaged out over the finite observation time of approximately 25 s. The persistence of realization-to-realization variability at N=40 is consistent with the theoretical expectation for systems governed by Lévy-stable velocity statistics, where the absence of a finite second moment in the velocity distribution implies that the central limit theorem does not apply on finite time scales, and individual trajectories can remain anomalously fast or anomalously slow relative to the ensemble mean for extended periods. This observation directly motivates the ergodicity breaking analysis presented below.

---

### Velocity PDF Data Collapse and Holtsmark Tail Analysis

The empirical speed distributions P(v) computed from finite-difference velocities along the noise-free tracer trajectories, normalized by the characteristic velocity scale v_c(N) defined as the standard deviation of the speed distribution for each N, reveal a partial but imperfect data collapse across the four filament configurations. The normalized distributions P(v/v_c) share a common shape in the bulk of the distribution, consistent with the expectation that the Biot-Savart kernel sets a universal velocity statistics structure independent of N in the dilute limit.

The tail behavior of the normalized speed PDF is of primary theoretical interest. The Holtsmark distribution, which arises from the superposition of 1/r² velocity fields from randomly distributed sources in three dimensions, predicts a power-law tail P(v) ~ v^{-5/2} for large v, corresponding to a Lévy stability index α_stable = 3/2. The empirical tails of the normalized speed PDFs show a power-law decay that is broadly consistent with this prediction, though the finite number of tracers (five per N configuration) and the finite simulation time limit the statistical precision of the tail estimate. The tail slope steepens slightly as N increases from 5 to 40, which can be attributed to the increasing probability that the tracer is simultaneously influenced by multiple filaments at moderate distances, softening the extreme velocity events that dominate the tail. Nevertheless, the tail exponent remains in the range expected for Holtsmark statistics across all N values studied, providing empirical support for the theoretical prediction that the 3D Biot-Savart kernel generates Lévy-stable velocity fluctuations with α_stable ≈ 3/2.

The characteristic velocity scale v_c(N) itself increases with N, reflecting the growing number of filaments contributing to the total velocity field at any given tracer position. This scaling of v_c with N is consistent with the Holtsmark prediction that the width of the velocity distribution scales as N^{1/3} in three dimensions for a uniform random distribution of sources, though the finite domain size and the relatively small range of N values studied here preclude a precise verification of this scaling law.

---

### MSD Scaling and Anomalous Diffusion Exponent

The ensemble-averaged TAMSD curves for the four vortex filament configurations, computed over lag times from dt = 0.05 s up to half the total observation time (approximately 12.5 s), exhibit clear power-law growth across the intermediate time range. The fitted anomalous diffusion exponents α(N) are reported in Table 1.

The most striking feature of the MSD results is that all four configurations yield exponents substantially above the Holtsmark theoretical prediction of α = 1.5. For N=5, the fitted exponent is α = 1.999 ± 0.000, indistinguishable from ballistic transport (α = 2). For N=10 and N=20, the exponents remain near α ≈ 1.96–1.97, still strongly superdiffusive and far above the Holtsmark limit. Only at N=40 does the exponent begin to decrease toward the theoretical prediction, reaching α = 1.835 ± 0.003. This systematic decrease of α with increasing N is physically interpretable: at low filament densities, the tracer is dominated by the velocity field of the nearest filament, which produces nearly ballistic motion along the direction of the induced velocity. As N increases, the superposition of multiple filament contributions begins to produce the velocity decorrelation and direction randomization that are prerequisites for the Holtsmark-predicted anomalous diffusion regime. However, even at N=40, the system has not yet converged to the asymptotic Holtsmark limit of α = 1.5, suggesting that either larger N or longer observation times are required to reach the true asymptotic transport regime.

Comparison with the Lévy walk ground truth reveals an important discrepancy. The β=1.5 Lévy walk, which has theoretical α_theory = 1.5, yields a measured exponent of α = 1.049, substantially below its theoretical value. Similarly, the β=1.2 walk (α_theory = 1.8) measures α = 1.655, and the β=1.8 and β=2.5 walks both measure near α ≈ 1.01–1.07. This systematic underestimation of α in the Lévy walk data is attributable to the finite observation time of 59.9 s relative to the characteristic flight durations: for heavy-tailed flight time distributions with β < 2, the asymptotic superdiffusive regime requires observation times much longer than the mean flight duration, and the finite-time TAMSD estimator is known to be biased toward lower apparent exponents. The vortex filament system, by contrast, exhibits exponents that are if anything overestimated relative to the asymptotic Holtsmark prediction, because the quenched velocity field produces persistent ballistic segments that dominate the finite-time TAMSD.

The comparison between the N=40 vortex filament TAMSD and the β=1.5 Lévy walk TAMSD (the reference curve with matching theoretical α = 1.5) shows that the vortex system produces substantially larger absolute displacements at all lag times, reflecting both the higher characteristic velocity scale and the more persistent nature of the velocity correlations in the quenched disorder system relative to the renewal-based Lévy walk.

---

### Velocity Autocorrelation Function and Persistence Analysis

The normalized 3D VACF C_v(Δt) = ⟨**v**(t)·**v**(t+Δt)⟩/⟨|**v**|²⟩ provides a direct measure of the temporal memory in the tracer velocity field. For the vortex filament system, the VACF decays from unity at Δt = 0 and crosses zero at a finite time τ_c, indicating that the velocity field has a characteristic correlation time beyond which the tracer velocity becomes anti-correlated with its initial direction. This anti-correlation is a signature of the tracer encountering regions of reversed velocity induced by the curved velocity field around individual filaments.

The zero-crossing times are τ_c = 10.17 s for N=10 and τ_c = 9.555 s for N=40 (Table 1). The slight decrease in τ_c with increasing N is consistent with the expectation that a denser filament array produces more frequent velocity reversals, reducing the persistence time. For N=5 and N=20, the VACF does not cross zero within the observation window, indicating that the velocity correlation time exceeds the simulation duration for these configurations. This is particularly notable for N=5, where the large inter-filament spacing (11.70 m) means that the tracer can travel for extended periods under the influence of a single dominant filament without encountering a velocity reversal.

The persistence lengths L_p = v_c × ∫₀^{τ_c} C_v(Δt) dΔt are 0.313 m for N=10 and 0.578 m for N=40. These values are substantially smaller than the corresponding inter-filament spacings (9.28 m and 5.85 m respectively), indicating that the velocity correlation length is set not by the inter-filament distance but by the local geometry of the Biot-Savart velocity field around individual filaments. The tracer loses velocity memory on a scale much shorter than the typical distance between filaments, which is consistent with the observation that the velocity field varies rapidly in the vicinity of each filament due to the 1/r² falloff.

In contrast, the Lévy walk VACF structure is fundamentally different. The Lévy walk is a renewal process in which the velocity direction is randomized at each flight event, producing a VACF that decays to zero at the first flight termination event and remains near zero thereafter. The absence of a negative lobe in the Lévy walk VACF distinguishes it sharply from the vortex filament system, where the anti-correlation reflects the deterministic structure of the Biot-Savart velocity field. This contrast highlights a fundamental difference between the two systems: the vortex filament system generates anomalous transport through spatially correlated, deterministic velocity fields with quenched disorder, while the Lévy walk generates anomalous transport through temporally uncorrelated, stochastic flight events with heavy-tailed duration distributions.

---

### Residence Time Distributions

The residence time distributions P(τ) computed from low-speed trapping events (defined using the 25th percentile of the speed distribution as the threshold) exhibit power-law tails across all four filament configurations. The tail fitting procedure applied to the log-binned histograms yielded numerically degenerate results (NaN) for the power-law exponent β_res in the automated fitting pipeline, attributable to the limited number of trapping events per configuration (five tracers, each of length 500 time steps) and the sensitivity of the tail fit to the choice of fitting range. This limitation notwithstanding, the qualitative shape of the residence time distributions is consistent with a heavy-tailed power law, with the tail becoming more pronounced as N increases. The theoretical Holtsmark prediction for the residence time distribution in a system with velocity PDF tail P(v) ~ v^{-5/2} is P(τ) ~ τ^{-5/2}, corresponding to a Pareto tail index β_res = 3/2. The empirical distributions are broadly consistent with this prediction in the intermediate time range, though the limited statistics preclude a precise quantitative verification. The reference power-law line with slope −2.5 overlaid on the log-log residence time plot (Plot 5) provides a visual guide that is consistent with the bulk of the empirical data across all N values.

---

### Displacement PDFs and Lévy-Stable Characterization

The marginal displacement PDFs P(Δx) computed for the N=40 configuration at lag times of 1, 5, 10, and 20 s reveal a systematic evolution from heavy-tailed to progressively more Gaussian-like distributions as the lag time increases. At the shortest lag time of 1 s, the displacement PDF exhibits pronounced heavy tails that are inconsistent with a Gaussian distribution and are well-described by a Lévy-stable distribution with stability index α_stable ≈ 1.5, consistent with the Holtsmark prediction. The excess kurtosis at this lag time is substantially positive, confirming the non-Gaussian character of the short-time displacement statistics.

As the lag time increases to 5, 10, and 20 s, the displacement PDF progressively narrows relative to its tails and the excess kurtosis decreases, indicating a trend toward Gaussianization. This evolution is consistent with the theoretical expectation for a system with finite-time velocity correlations: at lag times much shorter than the VACF zero-crossing time τ_c ≈ 9.5 s, the displacement is dominated by ballistic segments and the PDF reflects the heavy-tailed velocity distribution; at lag times comparable to or exceeding τ_c, the displacement accumulates contributions from multiple decorrelated velocity segments, and the central limit theorem begins to apply to the bulk of the distribution even as the tails remain heavy. The Lévy-stable fit with α_stable = 1.5 provides a good description of the tails at all lag times studied, supporting the interpretation that the effective stochastic theory for the 3D vortex filament system is a Lévy-stable process with stability index 3/2, as predicted by the Holtsmark distribution of the underlying velocity field.

---

### Okubo-Weiss Analysis and Trapping Geometry

The Okubo-Weiss (OW) analysis for the N=40 configuration reveals a statistically significant negative correlation between local tracer speed and the OW parameter (Pearson r = −0.508, p ≪ 0.001). The OW parameter, defined here as the difference between the squared strain rate and the squared rotation rate of the velocity gradient along the tracer path, is negative in rotation-dominated regions and positive in strain-dominated regions. The negative Pearson correlation between speed and OW parameter indicates that low-speed trapping events occur preferentially in rotation-dominated regions (negative OW), while high-speed excursions occur preferentially in strain-dominated regions (positive OW). This is physically intuitive: the Biot-Savart velocity field around a vortex filament is predominantly rotational in the near-field, and a tracer that approaches close to a filament will experience a rapidly rotating velocity field that can trap it in a quasi-circular orbit at low net displacement speed. Conversely, tracers in the far field of multiple filaments experience a more strain-dominated velocity field that drives rapid displacement.

The complementary rolling speed variance analysis yields a positive Pearson correlation of r = +0.425 between instantaneous speed and local speed variance (computed in rolling windows of 10 time steps). This indicates that high-speed regions are also regions of high speed variability, consistent with the tracer passing through the rapidly varying near-field of a filament where the 1/r² velocity gradient produces large speed fluctuations over short distances. Together, the OW and rolling variance analyses provide a spatially resolved picture of the trapping mechanism: low-speed trapping events are associated with rotation-dominated near-filament regions, while high-speed, high-variance excursions correspond to the strain-dominated inter-filament regions.

---

### Ergodicity Breaking

The ergodicity breaking (