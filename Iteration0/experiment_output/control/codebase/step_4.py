# filename: codebase/step_4.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import stats
from scipy.stats import levy_stable
import warnings
from step_1 import (extract_tracer_trajectories_vortex, extract_tracer_trajectories_levy, load_and_verify_vortex, load_and_verify_levy, compute_velocities_central_diff)
VF_PATH = '/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy'
LW_PATH = '/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy'
DATA_DIR = 'data/'
def compute_hill_estimator(data, k_min=10, k_max=None):
    data_sorted = np.sort(data)
    n = len(data_sorted)
    if k_max is None:
        k_max = n // 2
    k_max = min(k_max, n - 1)
    k_arr = np.arange(k_min, k_max + 1)
    alpha_hill = np.empty(len(k_arr))
    for idx, k in enumerate(k_arr):
        threshold = data_sorted[n - k - 1]
        if threshold <= 0:
            alpha_hill[idx] = np.nan
            continue
        log_ratios = np.log(data_sorted[n - k:] / threshold)
        h_k = np.mean(log_ratios)
        if h_k <= 0:
            alpha_hill[idx] = np.nan
        else:
            alpha_hill[idx] = 1.0 / h_k
    return k_arr, alpha_hill
def find_hill_plateau(k_arr, alpha_hill, window=30, min_k=20):
    n = len(k_arr)
    best_var = np.inf
    best_start = 0
    start_search = 0
    for i in range(n):
        if k_arr[i] >= min_k:
            start_search = i
            break
    for i in range(start_search, n - window + 1):
        seg = alpha_hill[i:i + window]
        if not np.all(np.isfinite(seg)):
            continue
        v = np.var(seg)
        if v < best_var:
            best_var = v
            best_start = i
    best_end = best_start + window - 1
    best_end = min(best_end, n - 1)
    k_opt_idx = (best_start + best_end) // 2
    k_opt = int(k_arr[k_opt_idx])
    alpha_opt = float(np.median(alpha_hill[best_start:best_end + 1]))
    return k_opt, alpha_opt, best_start, best_end
def fit_ccdf_powerlaw(data, frac_tail=0.1):
    data_sorted = np.sort(data)
    n = len(data_sorted)
    ccdf = 1.0 - np.arange(1, n + 1) / n
    n_tail = max(int(n * frac_tail), 10)
    x_tail = data_sorted[n - n_tail:]
    ccdf_tail = ccdf[n - n_tail:]
    valid = (x_tail > 0) & (ccdf_tail > 0)
    x_fit = x_tail[valid]
    c_fit = ccdf_tail[valid]
    if len(x_fit) < 3:
        return np.nan, np.nan, x_tail, ccdf_tail, np.nan
    slope, intercept, r_val, _, se = stats.linregress(np.log(x_fit), np.log(c_fit))
    return float(-slope), float(se), x_tail, ccdf_tail, float(r_val**2)
def compute_displacement_pdf(tracer_dict, keys_list, lag_steps_list, dt):
    disp_data = {}
    for lag in lag_steps_list:
        dx_all = []
        for key in keys_list:
            traj = tracer_dict[key]
            x = traj['x_true']
            dx = x[lag:] - x[:-lag]
            dx_all.append(dx)
        disp_data[lag] = np.concatenate(dx_all)
    return disp_data
def fit_levy_stable_displacement(dx_arr, alpha_bounds=(1.0, 2.0)):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            params = levy_stable.fit(dx_arr, f0=None, f1=0.0)
            alpha_fit = float(np.clip(params[0], alpha_bounds[0], alpha_bounds[1]))
            beta_fit = float(params[1])
            loc_fit = float(params[2])
            scale_fit = float(params[3])
        except Exception:
            alpha_fit = np.nan
            beta_fit = np.nan
            loc_fit = np.nan
            scale_fit = np.nan
    return alpha_fit, beta_fit, loc_fit, scale_fit
if __name__ == '__main__':
    print('=== STEP 4: Velocity PDF Tail Analysis and Holtsmark Distribution Testing ===\n')
    vf_data = load_and_verify_vortex(VF_PATH)
    vf_tracers = extract_tracer_trajectories_vortex(vf_data)
    N_vals = [5, 10, 20, 40]
    print('\n--- Collecting velocities for each N ---')
    vf_vel = {}
    for N in N_vals:
        keys_N = [(n, tid) for (n, tid) in vf_tracers.keys() if n == N]
        vx_list, vy_list, vz_list = [], [], []
        for key in keys_N:
            traj = vf_tracers[key]
            vx, vy, vz, _ = compute_velocities_central_diff(traj['x_true'], traj['y_true'], traj['z_true'], traj['time'])
            vx_list.append(vx)
            vy_list.append(vy)
            vz_list.append(vz)
        vx_all = np.concatenate(vx_list)
        vy_all = np.concatenate(vy_list)
        vz_all = np.concatenate(vz_list)
        speeds = np.sqrt(vx_all**2 + vy_all**2 + vz_all**2)
        vf_vel[N] = {'vx': vx_all, 'vy': vy_all, 'vz': vz_all, 'speeds': speeds, 'keys': keys_N}
        print('N=' + str(N) + ': n_samples=' + str(len(vx_all)) + ', speed range=[' + str(round(speeds.min(), 4)) + ', ' + str(round(speeds.max(), 4)) + '] m/s')
    print('\n--- Hill Estimator Analysis (Vortex Filaments) ---')
    for N in N_vals:
        speeds = vf_vel[N]['speeds']
        k_arr, alpha_hill = compute_hill_estimator(speeds, k_min=10, k_max=len(speeds) // 2)
        k_opt, alpha_opt, p_start, p_end = find_hill_plateau(k_arr, alpha_hill, window=30, min_k=20)
        print('N=' + str(N) + ': optimal k=' + str(k_opt) + ', Hill alpha_stable=' + str(round(alpha_opt, 4)) + ' (Holtsmark theory: 1.5), plateau k range=[' + str(int(k_arr[p_start])) + ', ' + str(int(k_arr[p_end])) + ']')