# filename: codebase/step_2.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import stats
import os
from step_1 import (
    extract_tracer_trajectories_vortex,
    extract_tracer_trajectories_levy,
    load_and_verify_vortex,
    load_and_verify_levy,
)

VF_PATH = '/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy'
LW_PATH = '/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy'
DATA_DIR = 'data/'

def compute_tamsd_vectorized(x, y, z, max_lag_fraction=0.5):
    T = len(x)
    max_lag = int(T * max_lag_fraction)
    tamsd = np.empty(max_lag)
    for lag in range(1, max_lag + 1):
        dx = x[lag:] - x[:-lag]
        dy = y[lag:] - y[:-lag]
        dz = z[lag:] - z[:-lag]
        tamsd[lag - 1] = np.mean(dx**2 + dy**2 + dz**2)
    return tamsd, np.arange(1, max_lag + 1)

def compute_ensemble_tamsd(tracer_dict, keys_list, max_lag_fraction=0.5):
    per_tracer_tamsd = []
    lag_indices_ref = None
    dt = None
    for key in keys_list:
        traj = tracer_dict[key]
        x, y, z = traj['x_true'], traj['y_true'], traj['z_true']
        time = traj['time']
        if dt is None:
            dt = float(time[1] - time[0])
        tamsd, lag_idx = compute_tamsd_vectorized(x, y, z, max_lag_fraction=max_lag_fraction)
        per_tracer_tamsd.append(tamsd)
        if lag_indices_ref is None:
            lag_indices_ref = lag_idx
    min_len = min(len(t) for t in per_tracer_tamsd)
    per_tracer_arr = np.array([t[:min_len] for t in per_tracer_tamsd])
    ensemble_tamsd = np.mean(per_tracer_arr, axis=0)
    return per_tracer_tamsd, ensemble_tamsd, lag_indices_ref[:min_len], dt

def compute_local_slope(lag_times, tamsd, window=7):
    log_t = np.log(lag_times)
    log_msd = np.log(tamsd)
    n = len(log_t)
    half = window // 2
    slopes = np.full(n, np.nan)
    for i in range(n):
        i0 = max(0, i - half)
        i1 = min(n, i + half + 1)
        if i1 - i0 < 2:
            continue
        slope, _, _, _, _ = stats.linregress(log_t[i0:i1], log_msd[i0:i1])
        slopes[i] = slope
    return lag_times, slopes

def find_stable_scaling_window(lag_times, slopes, tolerance=0.1, min_points=5):
    T_max = lag_times[-1]
    valid = (lag_times >= 0.05 * T_max) & (lag_times <= 0.6 * T_max) & np.isfinite(slopes)
    valid_idx = np.where(valid)[0]
    if len(valid_idx) < min_points:
        mask = valid.copy()
        med = float(np.nanmedian(slopes[valid])) if valid.any() else 1.0
        if valid_idx.size > 0:
            return lag_times[valid_idx[0]], lag_times[valid_idx[-1]], med, mask
        else:
            return lag_times[0], lag_times[-1], med, mask
    best_start = valid_idx[0]
    best_end = valid_idx[0]
    best_len = 0
    n_valid = len(valid_idx)
    for ii in range(n_valid):
        i = valid_idx[ii]
        run_slopes = [slopes[i]]
        for jj in range(ii + 1, n_valid):
            j = valid_idx[jj]
            if j != valid_idx[jj - 1] + 1:
                break
            run_slopes.append(slopes[j])
            seg = np.array(run_slopes)
            med = np.median(seg)
            if np.all(np.abs(seg - med) <= tolerance):
                if len(seg) > best_len:
                    best_len = len(seg)
                    best_start = i
                    best_end = j
    mask = np.zeros(len(lag_times), dtype=bool)
    mask[best_start:best_end + 1] = True
    med_slope = float(np.median(slopes[best_start:best_end + 1]))
    return float(lag_times[best_start]), float(lag_times[best_end]), med_slope, mask

def fit_powerlaw_in_window(lag_times, tamsd, window_mask):
    lt = lag_times[window_mask]
    ms = tamsd[window_mask]
    valid = (lt > 0) & (ms > 0) & np.isfinite(ms)
    lt = lt[valid]
    ms = ms[valid]
    if len(lt) < 2:
        return np.nan, np.nan, np.nan, np.nan
    slope, intercept, r_value, _, stderr = stats.linregress(np.log(lt), np.log(ms))
    return float(slope), float(stderr), float(intercept), float(r_value**2)

def process_dataset(tracer_dict, group_key, group_vals, max_lag_fraction=0.5):
    results = {}
    for gval in group_vals:
        keys_g = [k for k in tracer_dict.keys() if k[0] == gval]
        per_tracer_tamsd, ensemble_tamsd, lag_idx, dt = compute_ensemble_tamsd(
            tracer_dict, keys_g, max_lag_fraction=max_lag_fraction
        )
        lag_times = lag_idx * dt
        _, slopes = compute_local_slope(lag_times, ensemble_tamsd, window=7)
        win_start, win_end, med_slope, win_mask = find_stable_scaling_window(
            lag_times, slopes, tolerance=0.1, min_points=5
        )
        alpha, alpha_se, log_A, r2 = fit_powerlaw_in_window(lag_times, ensemble_tamsd, win_mask)
        results[gval] = {
            'per_tracer_tamsd': per_tracer_tamsd,
            'ensemble_tamsd': ensemble_tamsd,
            'lag_times': lag_times,
            'lag_indices': lag_idx,
            'slopes': slopes,
            'win_start': win_start,
            'win_end': win_end,
            'win_mask': win_mask,
            'alpha': alpha,
            'alpha_se': alpha_se,
            'log_A': log_A,
            'r2': r2,
            'dt': dt,
            'n_tracers': len(keys_g),
        }
    return results

if __name__ == '__main__':
    vf_data = load_and_verify_vortex(VF_PATH)
    lw_data = load_and_verify_levy(LW_PATH)
    vf_tracers = extract_tracer_trajectories_vortex(vf_data)
    lw_tracers = extract_tracer_trajectories_levy(lw_data)
    N_vals = sorted(set(k[0] for k in vf_tracers.keys()))
    beta_vals = sorted(set(k[0] for k in lw_tracers.keys()))
    vf_results = process_dataset(vf_tracers, 0, N_vals, max_lag_fraction=0.5)
    lw_results = process_dataset(lw_tracers, 0, beta_vals, max_lag_fraction=0.5)
    for N in N_vals:
        r = vf_results[N]
        print('N = ' + str(N) + ' filaments: alpha = ' + str(round(r['alpha'], 4)))
    for b in beta_vals:
        r = lw_results[b]
        print('beta = ' + str(b) + ': alpha = ' + str(round(r['alpha'], 4)))