# filename: codebase/step_4.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy.signal import savgol_filter
from scipy import stats
from step_1 import (
    extract_tracer_trajectories_vortex,
    extract_tracer_trajectories_levy,
    load_and_verify_vortex,
    load_and_verify_levy,
    compute_velocities_central_diff,
)
VF_PATH = '/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy'
LW_PATH = '/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy'
DATA_DIR = 'data/'
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
def compute_all_tamsd(tracer_dict, keys_list, max_lag_fraction=0.5):
    per_tracer = []
    dt = None
    for key in keys_list:
        traj = tracer_dict[key]
        x, y, z, time = traj['x_true'], traj['y_true'], traj['z_true'], traj['time']
        if dt is None:
            dt = float(time[1] - time[0])
        tamsd, _ = compute_tamsd_single(x, y, z, max_lag_fraction)
        per_tracer.append(tamsd)
    min_len = min(len(t) for t in per_tracer)
    per_tracer_arr = np.array([t[:min_len] for t in per_tracer])
    ensemble_mean = np.mean(per_tracer_arr, axis=0)
    lag_times = np.arange(1, min_len + 1) * dt
    return per_tracer_arr, ensemble_mean, lag_times, dt
def compute_local_slope_savgol(lag_times, tamsd, window_length=15, polyorder=2):
    valid = (lag_times > 0) & (tamsd > 0) & np.isfinite(tamsd)
    log_t = np.where(valid, np.log(lag_times), 0.0)
    log_msd = np.where(valid, np.log(tamsd), 0.0)
    n = len(log_t)
    wl = min(window_length, n)
    if wl % 2 == 0:
        wl -= 1
    if wl < 3:
        return np.full(n, np.nan)
    smooth_log_msd = savgol_filter(log_msd, window_length=wl, polyorder=min(polyorder, wl - 1))
    dlog_t = np.gradient(log_t)
    dlog_t_safe = np.where(np.abs(dlog_t) > 1e-15, dlog_t, 1e-15)
    slopes = np.gradient(smooth_log_msd) / dlog_t_safe
    slopes[~valid] = np.nan
    return slopes
def find_stable_window(lag_times, slopes, t_min=1.0, last_frac=0.2, tolerance=0.15, min_points=5):
    T_max = lag_times[-1]
    t_cutoff = T_max * (1.0 - last_frac)
    candidate = (lag_times >= t_min) & (lag_times <= t_cutoff) & np.isfinite(slopes)
    candidate_idx = np.where(candidate)[0]
    if len(candidate_idx) < min_points:
        mask = candidate.copy()
        med = float(np.nanmedian(slopes[candidate])) if candidate.any() else 1.0
        t_s = float(lag_times[candidate_idx[0]]) if len(candidate_idx) > 0 else t_min
        t_e = float(lag_times[candidate_idx[-1]]) if len(candidate_idx) > 0 else T_max
        return mask, med, t_s, t_e
    best_start = candidate_idx[0]
    best_end = candidate_idx[0]
    best_len = 1
    n = len(candidate_idx)
    for ii in range(n):
        i = candidate_idx[ii]
        run = [i]
        for jj in range(ii + 1, n):
            j = candidate_idx[jj]
            if j != candidate_idx[jj - 1] + 1:
                break
            run.append(j)
            seg = slopes[np.array(run)]
            med = np.median(seg)
            if np.all(np.abs(seg - med) <= tolerance):
                if len(run) > best_len:
                    best_len = len(run)
                    best_start = run[0]
                    best_end = run[-1]
    win_mask = np.zeros(len(lag_times), dtype=bool)
    win_mask[best_start:best_end + 1] = True
    alpha_plateau = float(np.median(slopes[best_start:best_end + 1]))
    return win_mask, alpha_plateau, float(lag_times[best_start]), float(lag_times[best_end])
def fit_powerlaw(lag_times, tamsd, win_mask):
    lt = lag_times[win_mask]
    ms = tamsd[win_mask]
    valid = (lt > 0) & (ms > 0) & np.isfinite(ms)
    lt, ms = lt[valid], ms[valid]
    if len(lt) < 2:
        return np.nan, np.nan, np.nan, np.nan
    slope, intercept, r_val, _, stderr = stats.linregress(np.log(lt), np.log(ms))
    return float(slope), float(stderr), float(intercept), float(r_val**2)
def compute_vacf_vectorized(tracer_dict, keys_list, dt, max_lag_fraction=0.5):
    all_vx, all_vy, all_vz = [], [], []
    for key in keys_list:
        traj = tracer_dict[key]
        vx, vy, vz, _ = compute_velocities_central_diff(traj['x_true'], traj['y_true'], traj['z_true'], traj['time'])
        all_vx.append(vx)
        all_vy.append(vy)
        all_vz.append(vz)
    speeds_all = np.concatenate([np.sqrt(vx**2 + vy**2 + vz**2) for vx, vy, vz in zip(all_vx, all_vy, all_vz)])
    v_c = float(np.std(speeds_all))
    T_vel = len(all_vx[0])
    max_lag = int(T_vel * max_lag_fraction)
    v2_all = np.concatenate([vx**2 + vy**2 + vz**2 for vx, vy, vz in zip(all_vx, all_vy, all_vz)])
    v2_mean = float(np.mean(v2_all))
    vacf_sum = np.zeros(max_lag)
    vacf_cnt = np.zeros(max_lag, dtype=np.int64)
    for vx, vy, vz in zip(all_vx, all_vy, all_vz):
        T = len(vx)
        ml = min(max_lag, T - 1)
        for lag in range(1, ml + 1):
            n_pairs = T - lag
            vacf_sum[lag - 1] += (np.dot(vx[:n_pairs], vx[lag:lag + n_pairs]) + np.dot(vy[:n_pairs], vy[lag:lag + n_pairs]) + np.dot(vz[:n_pairs], vz[lag:lag + n_pairs]))
            vacf_cnt[lag - 1] += n_pairs
    safe_cnt = np.where(vacf_cnt > 0, vacf_cnt, 1)
    vacf_norm = (vacf_sum / safe_cnt) / (v2_mean if v2_mean > 0 else 1.0)
    lag_times_vacf = np.arange(1, max_lag + 1) * dt
    return vacf_norm, lag_times_vacf, v_c
def find_zero_crossing(vacf, lag_times):
    for i in range(len(vacf) - 1):
        if np.isfinite(vacf[i]) and np.isfinite(vacf[i + 1]):
            if vacf[i] >= 0 and vacf[i + 1] < 0:
                frac = vacf[i] / (vacf[i] - vacf[i + 1])
                return float(lag_times[i] + frac * (lag_times[i + 1] - lag_times[i]))
    return np.nan
def compute_persistence_length(vacf, lag_times, v_c, tau_c):
    if not np.isfinite(tau_c):
        return np.nan
    mask = lag_times <= tau_c
    if mask.sum() < 2:
        return np.nan
    lt_cut = np.concatenate([[0.0], lag_times[mask]])
    vacf_cut = np.concatenate([[1.0], vacf[mask]])
    return float(v_c * np.trapezoid(vacf_cut, lt_cut))
def compute_eb_curve(per_tracer_arr, ensemble_mean):
    var_tamsd = np.var(per_tracer_arr, axis=0, ddof=1)
    mean_sq = ensemble_mean**2
    safe_mean_sq = np.where(mean_sq > 0, mean_sq, np.nan)
    eb = var_tamsd / safe_mean_sq
    return eb
if __name__ == '__main__':
    vf_data = load_and_verify_vortex(VF_PATH)
    lw_data = load_and_verify_levy(LW_PATH)
    vortex_tracers = extract_tracer_trajectories_vortex(vf_data)
    levy_tracers = extract_tracer_trajectories_levy(lw_data)
    N_vals = sorted(set(k[0] for k in vortex_tracers.keys()))
    vortex_results = {}
    for N in N_vals:
        keys_N = [k for k in vortex_tracers.keys() if k[0] == N]
        per_tracer_arr, ensemble_mean, lag_times, dt = compute_all_tamsd(vortex_tracers, keys_N)
        slopes = compute_local_slope_savgol(lag_times, ensemble_mean)
        win_mask, alpha_plateau, t_start, t_end = find_stable_window(lag_times, slopes)
        alpha, alpha_se, log_A, r2 = fit_powerlaw(lag_times, ensemble_mean, win_mask)
        vacf, vacf_lags, v_c = compute_vacf_vectorized(vortex_tracers, keys_N, dt)
        tau_c = find_zero_crossing(vacf, vacf_lags)
        L_p = compute_persistence_length(vacf, vacf_lags, v_c, tau_c)
        inter_filament = 20.0 * (N ** (-1.0 / 3.0))
        vortex_results[str(N)] = {'alpha': alpha, 'alpha_se': alpha_se, 'tau_c': tau_c, 'L_p': L_p, 'inter_filament': inter_filament}
        print('N=' + str(N) + ': alpha=' + str(round(alpha, 3)) + ' +/- ' + str(round(alpha_se, 3)) + ', tau_c=' + str(round(tau_c, 3)) + 's, L_p=' + str(round(L_p, 3)) + 'm, inter_filament=' + str(round(inter_filament, 3)) + 'm')
    np.savez(os.path.join(DATA_DIR, 'vortex_analysis_results.npz'), **vortex_results)
    print('Analysis complete.')