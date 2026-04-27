

Iteration 0:
**Project Summary: 3D Quenched Vortex Filament Superdiffusion**

**1. Objective & Model**
Characterize anomalous transport of passive tracers in a 3D quenched vortex filament field (Biot-Savart $1/r^2$ kernel). Compare against 3D isotropic Lévy walk ground truth to evaluate convergence toward the theoretical Holtsmark limit ($\alpha_{stable}=1.5$).

**2. Key Findings**
*   **Transport Regime:** All configurations ($N \in \{5, 10, 20, 40\}$) exhibit strong superdiffusion ($\alpha > 1.5$). Exponents decrease with $N$ ($\alpha \approx 2.0$ at $N=5$ to $\alpha \approx 1.84$ at $N=40$), indicating a slow approach to the Holtsmark-predicted $\alpha=1.5$.
*   **Quenched Disorder:** Unlike Lévy walks, the quenched (fixed) filament positions prevent VACF zero-crossings, leading to persistent velocity correlations and near-ballistic transport.
*   **Velocity PDF:** Hill estimation of the speed tail yields $\alpha_{stable} \in [2.7, 24.8]$, significantly higher than the theoretical 1.5. The system remains in a pre-asymptotic regime where finite-size effects and specific filament realizations dominate.
*   **Trapping:** Residence time distributions show heavy tails at $N=40$ ($\beta_{trap} \approx 0.81$), suggesting divergent mean trapping times, consistent with the observed superdiffusive scaling.

**3. Limitations & Uncertainties**
*   **Statistical Power:** Extremely low tracer count ($N_{tracers}=5$ per configuration) results in wide confidence intervals and high sensitivity to specific filament realizations.
*   **Finite-Time Effects:** The simulation duration is insufficient to reach the asymptotic diffusive limit, causing measured $\alpha$ to remain biased toward ballistic values.
*   **Hill Estimator:** Optimal $k$ values suggest the estimator is sampling the bulk rather than the true asymptotic tail.

**4. Future Directions**
*   **Increase Sample Size:** Increase tracer count to $O(100)$ per configuration to narrow confidence intervals and resolve the convergence trend of $\alpha(N)$.
*   **Domain Scaling:** Evaluate sensitivity to simulation box size to distinguish between finite-domain truncation and true asymptotic behavior.
*   **Annealed Dynamics:** Introduce filament-filament dynamics to test if evolving disorder recovers the expected Holtsmark/Lévy-stable transport limits.
*   **Refined Tail Analysis:** Implement more robust tail-fitting methods (e.g., Maximum Likelihood Estimation for power laws) to better estimate $\alpha_{stable}$ and $\beta_{trap}$.
        

Iteration 1:
**Methodological Evolution**
- **Analysis Scope:** This iteration introduces a 3D quenched vortex filament model, transitioning from the 2D point-vortex framework (Iteration 0) to a 3D Biot-Savart kernel ($1/r^2$ falloff).
- **Metric Refinement:** The MSD calculation was updated to 3D displacement ($\sqrt{x^2+y^2+z^2}$).
- **New Diagnostics:** Added Okubo-Weiss (OW) criterion analysis to correlate tracer trapping with local velocity gradient topology (rotation vs. strain) and computed the Ergodicity Breaking (EB) parameter to quantify non-self-averaging behavior.
- **Ground Truth Comparison:** Introduced a 3D isotropic Lévy walk as a baseline to contrast renewal-based stochastic transport with the deterministic, spatially correlated transport of the quenched filament system.

**Performance Delta**
- **Anomalous Exponent ($\alpha$):** The 3D system exhibits strong superdiffusion ($\alpha \approx 1.8–2.0$), significantly higher than the theoretical Holtsmark prediction ($\alpha=1.5$). This represents a regression in model accuracy relative to the theoretical limit compared to 2D results, which aligned more closely with Cauchy-based predictions.
- **Robustness:** The system shows high realization-dependent variance, particularly at low $N$. While increasing $N$ to 40 improved convergence toward the Holtsmark limit, the system remains non-ergodic on the observed time scales.
- **Interpretability:** The OW analysis successfully identified that low-speed trapping occurs in rotation-dominated regions, providing a clearer fluid-dynamical mechanism for anomalous transport than the purely statistical 2D approach.

**Synthesis**
- **Causal Attribution:** The observed overestimation of $\alpha$ (relative to $\alpha=1.5$) is attributed to the quenched nature of the filaments. Unlike the Lévy walk, which randomizes velocity directions at each renewal event, the Biot-Savart field induces persistent ballistic segments as tracers traverse the near-field of individual filaments.
- **Validity and Limits:** The 3D system is not yet in the asymptotic regime. The discrepancy between the measured $\alpha$ and the Holtsmark prediction suggests that the "trapping" interpretation of 2D systems must be replaced by a "geometric manifold" interpretation in 3D, where transport is dictated by the topology of the velocity gradient tensor.
- **Next Steps:** The failure of the finite-time TAMSD to reach the theoretical $\alpha=1.5$ limit indicates that future iterations must either increase the simulation duration or employ ensemble-averaging over multiple filament configurations to mitigate the bias introduced by quenched disorder.
        