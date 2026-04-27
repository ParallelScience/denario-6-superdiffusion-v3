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

def compute_local_slope(lag_times, tamsd, window=5):
    log_t, log_msd = np.log(lag_times), np.log(tamsd)
    n = len(log_t)
    half = window // 2
    slopes = np.zeros(n)
    for i in range(n):
        i0, i1 = max(0, i - half), min(n, i + half + 1)
        if i1 - i0 < 2:
            slopes[i] = np.nan
            continue
        slope, _, _, _, _ = stats.linregress(log_t[i0:i1], log_msd[i0:i1])
        slopes[i] = slope
    return lag_times, slopes

def find_stable_scaling_window(lag_times, slopes, tolerance=0.1, min_points=5):
    T_max = lag_times[-1]
    valid = (lag_times >= 0.05 * T_max) & (lag_times <= 0.6 * T_max) & np.isfinite(slopes)
    valid_idx = np.where(valid)[0]
    if len(valid_idx) < min_points:
        return lag_times[0], lag_times[-1], np.nanmedian(slopes[valid]), valid
    best_start, best_end, best_len = valid_idx[0], valid_idx[0], 0
    for i in valid_idx:
        for j in valid_idx:
            if j < i: continue
            seg = slopes[i:j+1]
            if len(seg) < min_points: continue
            med = np.median(seg)
            if np.all(np.abs(seg - med) <= tolerance):
                if (j - i + 1) > best_len:
                    best_len, best_start, best_end = j - i + 1, i, j
    mask = np.zeros(len(lag_times), dtype=bool)
    mask[best_start:best_end+1] = True
    return lag_times[best_start], lag_times[best_end], np.median(slopes[best_start:best_end+1]), mask

def fit_powerlaw_in_window(lag_times, tamsd, window_mask):
    lt, ms = lag_times[window_mask], tamsd[window_mask]
    valid = (lt > 0) & (ms > 0)
    slope, intercept, r_value, _, stderr = stats.linregress(np.log(lt[valid]), np.log(ms[valid]))
    return slope, stderr, intercept, r_value**2

if __name__ == '__main__':
    vf_data = load_and_verify_vortex()
    lw_data = load_and_verify_levy()
    print("Processing complete.")