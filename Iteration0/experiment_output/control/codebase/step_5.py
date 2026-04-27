# filename: codebase/step_5.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import stats
import os
from step_1 import (extract_tracer_trajectories_vortex, extract_tracer_trajectories_levy, load_and_verify_vortex, load_and_verify_levy)
VF_PATH = '/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy'
LW_PATH = '/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy'
DATA_DIR = 'data/'
N_BOOTSTRAP = 500
RNG_SEED = 42
def compute_tamsd_single(x, y, z, max_lag_fraction=0.5):
    T = len(x)
    max_lag = int(T * max_lag_fraction)
    tamsd = np.empty(max_lag)
    for lag in range(1, max_lag + 1):
        dx = x[lag:] - x[:-lag]
        dy = y[lag:] - y[:-lag]
        dz = z[lag:] - z[:-lag]
        tamsd[lag - 1] = np.mean(dx**2 + dy**2 + dz**2)
    return tamsd, np.arange(1, max_lag + 1)
def get_per_tracer_tamsd(tracer_dict, keys_list, max_lag_fraction=0.5):
    per_tracer_tamsd = []
    dt = None
    lag_times_ref = None
    for key in keys_list:
        traj = tracer_dict[key]
        x, y, z = traj['x_true'], traj['y_true'], traj['z_true']
        time = traj['time']
        if dt is None:
            dt = float(time[1] - time[0])
        tamsd, lag_idx = compute_tamsd_single(x, y, z, max_lag_fraction=max_lag_fraction)
        per_tracer_tamsd.append(tamsd)
        if lag_times_ref is None:
            lag_times_ref = lag_idx * dt
    min_len = min(len(t) for t in per_tracer_tamsd)
    per_tracer_tamsd = [t[:min_len] for t in per_tracer_tamsd]
    lag_times = lag_times_ref[:min_len]
    return per_tracer_tamsd, lag_times, dt
def compute_eb_at_lag(per_tracer_tamsd, lag_idx):
    vals = np.array([t[lag_idx] for t in per_tracer_tamsd])
    mean_sq = np.mean(vals**2)
    sq_mean = np.mean(vals)**2
    if sq_mean == 0.0:
        return np.nan
    return float((mean_sq - sq_mean) / sq_mean)
def compute_eb_vs_lag(per_tracer_tamsd, lag_times):
    n_lags = len(lag_times)
    eb_curve = np.empty(n_lags)
    for i in range(n_lags):
        eb_curve[i] = compute_eb_at_lag(per_tracer_tamsd, i)
    return eb_curve
def fit_powerlaw_window(lag_times, ensemble_tamsd, win_start_frac=0.1, win_end_frac=0.5):
    T_max = lag_times[-1]
    mask = (lag_times >= win_start_frac * T_max) & (lag_times <= win_end_frac * T_max)
    lt = lag_times[mask]
    ms = ensemble_tamsd[mask]
    valid = (lt > 0) & (ms > 0) & np.isfinite(ms)
    lt = lt[valid]
    ms = ms[valid]
    if len(lt) < 2:
        return np.nan, np.nan
    slope, _, _, _, se = stats.linregress(np.log(lt), np.log(ms))
    return float(slope), float(se)
def bootstrap_alpha_eb(per_tracer_tamsd, lag_times, n_bootstrap=500, rng=None, win_start_frac=0.1, win_end_frac=0.5):
    if rng is None:
        rng = np.random.default_rng(RNG_SEED)
    n_tracers = len(per_tracer_tamsd)
    lag_idx_fixed = int(len(lag_times) * 0.25)
    alpha_boot = np.empty(n_bootstrap)
    eb_boot = np.empty(n_bootstrap)
    tamsd_arr = np.array(per_tracer_tamsd)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n_tracers, size=n_tracers)
        sample = tamsd_arr[idx]
        ensemble = np.mean(sample, axis=0)
        alpha_b, _ = fit_powerlaw_window(lag_times, ensemble, win_start_frac, win_end_frac)
        alpha_boot[b] = alpha_b
        vals = sample[:, lag_idx_fixed]
        mean_sq = np.mean(vals**2)
        sq_mean = np.mean(vals)**2
        eb_boot[b] = (mean_sq - sq_mean) / sq_mean if sq_mean > 0 else np.nan
    valid_alpha = alpha_boot[np.isfinite(alpha_boot)]
    valid_eb = eb_boot[np.isfinite(eb_boot)]
    alpha_ci = (float(np.percentile(valid_alpha, 2.5)), float(np.percentile(valid_alpha, 97.5))) if len(valid_alpha) > 0 else (np.nan, np.nan)
    eb_ci = (float(np.percentile(valid_eb, 2.5)), float(np.percentile(valid_eb, 97.5))) if len(valid_eb) > 0 else (np.nan, np.nan)
    return alpha_ci, eb_ci, alpha_boot, eb_boot
if __name__ == '__main__':
    rng = np.random.default_rng(RNG_SEED)
    vf_data = load_and_verify_vortex(VF_PATH)
    lw_data = load_and_verify_levy(LW_PATH)
    vf_tracers = extract_tracer_trajectories_vortex(vf_data)
    lw_tracers = extract_tracer_trajectories_levy(lw_data)
    N_vals = [5, 10, 20, 40]
    beta_vals = [1.2, 1.5, 1.8, 2.5]
    for N in N_vals:
        keys_N = [(n, tid) for (n, tid) in vf_tracers.keys() if n == N]
        per_tracer_tamsd, lag_times, dt = get_per_tracer_tamsd(vf_tracers, keys_N, max_lag_fraction=0.5)
        alpha_ci, eb_ci, _, _ = bootstrap_alpha_eb(per_tracer_tamsd, lag_times, n_bootstrap=N_BOOTSTRAP, rng=rng)
        print('N=' + str(N) + ': alpha_CI=' + str(alpha_ci) + ', EB_CI=' + str(eb_ci))
    for beta in beta_vals:
        keys_b = [(b, tid) for (b, tid) in lw_tracers.keys() if b == beta]
        per_tracer_tamsd, lag_times, dt = get_per_tracer_tamsd(lw_tracers, keys_b, max_lag_fraction=0.5)
        alpha_ci, eb_ci, _, _ = bootstrap_alpha_eb(per_tracer_tamsd, lag_times, n_bootstrap=N_BOOTSTRAP, rng=rng)
        print('beta=' + str(beta) + ': alpha_CI=' + str(alpha_ci) + ', EB_CI=' + str(eb_ci))