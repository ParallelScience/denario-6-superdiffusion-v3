

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
        