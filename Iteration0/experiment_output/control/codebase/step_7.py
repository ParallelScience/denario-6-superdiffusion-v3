# filename: codebase/step_7.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import stats
from scipy.stats import levy_stable, kstest, norm
import warnings
import os
from step_1 import extract_tracer_trajectories_vortex, extract_tracer_trajectories_levy, compute_velocities_central_diff
from step_5 import compute_tamsd_single, get_per_tracer_tamsd, compute_eb_at_lag, fit_powerlaw_window
VF_PATH = '/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy'
LW_PATH = '/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy'
DATA_DIR = 'data/'
N_BOOTSTRAP = 500
RNG_SEED = 42
N_VALS = [5, 10, 20, 40]
BETA_VALS = [1.2, 1.5, 1.8, 2.5]
LAG_STEPS = [1, 5, 10, 20]
def bootstrap_alpha_eb_full(per_tracer_tamsd, lag_times, n_bootstrap=500, rng=None, win_start_frac=0.1, win_end_frac=0.5):
    if rng is None:
        rng = np.random.default_rng(RNG_SEED)
    n_tracers = len(per_tracer_tamsd)
    lag_idx_fixed = max(0, int(len(lag_times) * 0.25) - 1)
    tamsd_arr = np.array(per_tracer_tamsd)
    alpha_boot = np.empty(n_bootstrap)
    eb_boot = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n_tracers, size=n_tracers)
        sample = tamsd_arr[idx]
        ensemble = np.mean(sample, axis=0)
        alpha_b, _ = fit_powerlaw_window(lag_times, ensemble, win_start_frac, win_end_frac)
        alpha_boot[b] = alpha_b
        vals = sample[:, lag_idx_fixed]
        mean_sq = np.mean(vals ** 2)
        sq_mean = np.mean(vals) ** 2
        eb_boot[b] = (mean_sq - sq_mean) / sq_mean if sq_mean > 0 else np.nan
    valid_alpha = alpha_boot[np.isfinite(alpha_boot)]
    valid_eb = eb_boot[np.isfinite(eb_boot)]
    alpha_mean = float(np.mean(valid_alpha)) if len(valid_alpha) > 0 else np.nan
    alpha_ci = (float(np.percentile(valid_alpha, 2.5)), float(np.percentile(valid_alpha, 97.5))) if len(valid_alpha) > 0 else (np.nan, np.nan)
    alpha_width = alpha_ci[1] - alpha_ci[0] if np.isfinite(alpha_ci[0]) else np.nan
    eb_mean = float(np.mean(valid_eb)) if len(valid_eb) > 0 else np.nan
    eb_ci = (float(np.percentile(valid_eb, 2.5)), float(np.percentile(valid_eb, 97.5))) if len(valid_eb) > 0 else (np.nan, np.nan)
    eb_width = eb_ci[1] - eb_ci[0] if np.isfinite(eb_ci[0]) else np.nan
    return alpha_mean, alpha_ci, alpha_width, eb_mean, eb_ci, eb_width, alpha_boot, eb_boot
def fit_holtsmark_params(speeds):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            params = levy_stable.fit(speeds, f0=1.5, f1=0.0)
            scale = float(params[3])
            loc = float(params[2])
        except Exception:
            scale = np.nan
            loc = np.nan
    return scale, loc
def ks_test_holtsmark(speeds, scale, loc):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            ks_stat, p_value = kstest(speeds, lambda x: levy_stable.cdf(x, 1.5, 0.0, loc=loc, scale=scale))
        except Exception:
            ks_stat, p_value = np.nan, np.nan
    return float(ks_stat), float(p_value)
def ks_test_gaussian(dx_arr):
    mu = float(np.mean(dx_arr))
    sigma = float(np.std(dx_arr, ddof=1))
    if sigma == 0.0:
        return np.nan, np.nan, mu, sigma
    ks_stat, p_value = kstest(dx_arr, lambda x: norm.cdf(x, loc=mu, scale=sigma))
    return float(ks_stat), float(p_value), mu, sigma
if __name__ == '__main__':
    rng = np.random.default_rng(RNG_SEED)
    vf_raw = np.load(VF_PATH, allow_pickle=False)
    lw_raw = np.load(LW_PATH, allow_pickle=False)
    vf_tracers = extract_tracer_trajectories_vortex(vf_raw)
    lw_tracers = extract_tracer_trajectories_levy(lw_raw)
    print('=' * 70)
    print('STEP 7: STATISTICAL VALIDATION — BOOTSTRAPPING AND KS TESTS')
    print('=' * 70)
    for N in N_VALS:
        keys_N = [(n, tid) for (n, tid) in vf_tracers.keys() if n == N]
        per_tracer_tamsd, lag_times, dt = get_per_tracer_tamsd(vf_tracers, keys_N, max_lag_fraction=0.5)
        alpha_mean, alpha_ci, alpha_width, eb_mean, eb_ci, eb_width, _, _ = bootstrap_alpha_eb_full(per_tracer_tamsd, lag_times, n_bootstrap=N_BOOTSTRAP, rng=rng)
        print('N=' + str(N) + ' alpha=' + str(round(alpha_mean, 4)) + ' CI=' + str(alpha_ci) + ' EB=' + str(round(eb_mean, 4)))
    for beta in BETA_VALS:
        keys_b = [(b, tid) for (b, tid) in lw_tracers.keys() if b == beta]
        per_tracer_tamsd, lag_times, dt = get_per_tracer_tamsd(lw_tracers, keys_b, max_lag_fraction=0.5)
        alpha_mean, alpha_ci, alpha_width, eb_mean, eb_ci, eb_width, _, _ = bootstrap_alpha_eb_full(per_tracer_tamsd, lag_times, n_bootstrap=N_BOOTSTRAP, rng=rng)
        print('Beta=' + str(beta) + ' alpha=' + str(round(alpha_mean, 4)) + ' CI=' + str(alpha_ci) + ' EB=' + str(round(eb_mean, 4)))