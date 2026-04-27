# filename: codebase/step_3.py
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
    compute_velocities_central_diff,
)

VF_PATH = '/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy'
LW_PATH = '/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy'
DATA_DIR = 'data/'

def compute_vacf_fft(vx, vy, vz):
    T = len(vx)
    n_fft = 2 * T
    def autocorr_component(v):
        V = np.fft.rfft(v, n=n_fft)
        S = V * np.conj(V)
        ac = np.fft.irfft(S)[:T]
        lags = np.arange(T)
        counts = T - lags
        return ac / counts
    return autocorr_component(vx) + autocorr_component(vy) + autocorr_component(vz)

def compute_normalized_vacf_per_tracer(vx, vy, vz):
    vacf_raw = compute_vacf_fft(vx, vy, vz)
    c0 = vacf_raw[0]
    if c0 == 0.0:
        return np.zeros_like(vacf_raw)
    return vacf_raw / c0

def ensemble_vacf(per_tracer_velocities, max_lag_fraction=0.5):
    all_vacf = []
    dt = None
    for traj in per_tracer_velocities:
        vx, vy, vz, t = traj['vx'], traj['vy'], traj['vz'], traj['t']
        if dt is None:
            dt = float(t[1] - t[0]) if len(t) > 1 else 0.05
        vacf_norm = compute_normalized_vacf_per_tracer(vx, vy, vz)
        all_vacf.append(vacf_norm)
    min_len = min(len(v) for v in all_vacf)
    max_lag = int(min_len * max_lag_fraction)
    arr = np.array([v[:max_lag] for v in all_vacf])
    vacf_mean = np.mean(arr, axis=0)
    vacf_std = np.std(arr, axis=0, ddof=1) if arr.shape[0] > 1 else np.zeros(max_lag)
    lag_times = np.arange(max_lag) * dt
    return lag_times, vacf_mean, vacf_std, dt

def find_zero_crossing(lag_times, vacf):
    for i in range(1, len(vacf)):
        if vacf[i - 1] > 0 and vacf[i] <= 0:
            frac = vacf[i - 1] / (vacf[i - 1] - vacf[i])
            return float(lag_times[i - 1] + frac * (lag_times[i] - lag_times[i - 1]))
    return None

def compute_residence_times(vx, vy, vz, dt, threshold):
    speed = np.sqrt(vx**2 + vy**2 + vz**2)
    below = speed < threshold
    durations = []
    in_event = False
    count = 0
    for b in below:
        if b:
            in_event = True
            count += 1
        else:
            if in_event and count > 0:
                durations.append(count * dt)
            in_event = False
            count = 0
    if in_event and count > 0:
        durations.append(count * dt)
    return np.array(durations)

def pareto_mle_beta(durations, tau_min=None):
    if tau_min is None:
        tau_min = float(np.min(durations))
    mask = durations >= tau_min
    d = durations[mask]
    n = len(d)
    if n < 2:
        return np.nan, np.nan, n, tau_min
    log_ratios = np.log(d / tau_min)
    beta_hat = float(n / np.sum(log_ratios))
    beta_se = beta_hat / np.sqrt(n)
    return beta_hat, beta_se, n, tau_min

def process_vacf_and_residence_vortex(tracer_dict, N_vals):
    results = {}
    for N in N_vals:
        keys_N = [(n, tid) for (n, tid) in tracer_dict.keys() if n == N]
        per_tracer_velocities = []
        all_speeds = []
        dt_val = None
        for key in keys_N:
            traj = tracer_dict[key]
            vx, vy, vz, t_int = compute_velocities_central_diff(traj['x_true'], traj['y_true'], traj['z_true'], traj['time'])
            if dt_val is None:
                dt_val = float(traj['time'][1] - traj['time'][0])
            per_tracer_velocities.append({'vx': vx, 'vy': vy, 'vz': vz, 't': t_int})
            speed = np.sqrt(vx**2 + vy**2 + vz**2)
            all_speeds.append(speed)
        all_speeds_concat = np.concatenate(all_speeds)
        threshold = float(np.percentile(all_speeds_concat, 25))
        lag_times, vacf_mean, vacf_std, dt = ensemble_vacf(per_tracer_velocities, max_lag_fraction=0.5)
        tau_c = find_zero_crossing(lag_times, vacf_mean)
        all_durations = []
        for ptv in per_tracer_velocities:
            durs = compute_residence_times(ptv['vx'], ptv['vy'], ptv['vz'], dt_val, threshold)
            all_durations.append(durs)
        all_durations_concat = np.concatenate(all_durations) if all_durations else np.array([])
        if len(all_durations_concat) >= 2:
            beta_hat, beta_se, n_used, tau_min_used = pareto_mle_beta(all_durations_concat)
        else:
            beta_hat, beta_se, n_used, tau_min_used = np.nan, np.nan, 0, np.nan
        results[N] = {'lag_times': lag_times, 'vacf_mean': vacf_mean, 'vacf_std': vacf_std, 'tau_c': tau_c, 'threshold': threshold, 'all_durations': all_durations_concat, 'beta_trap': beta_hat, 'beta_trap_se': beta_se, 'n_residence': n_used, 'tau_min': tau_min_used, 'dt': dt_val, 'n_tracers': len(keys_N)}
    return results

def process_vacf_and_residence_levy(tracer_dict, beta_vals):
    results = {}
    for beta in beta_vals:
        keys_b = [(b, tid) for (b, tid) in tracer_dict.keys() if b == beta]
        per_tracer_velocities = []
        all_speeds = []
        dt_val = None
        for key in keys_b:
            traj = tracer_dict[key]
            vx, vy, vz, t_int = compute_velocities_central_diff(traj['x_true'], traj['y_true'], traj['z_true'], traj['time'])
            if dt_val is None:
                dt_val = float(traj['time'][1] - traj['time'][0])
            per_tracer_velocities.append({'vx': vx, 'vy': vy, 'vz': vz, 't': t_int})
            speed = np.sqrt(vx**2 + vy**2 + vz**2)
            all_speeds.append(speed)
        all_speeds_concat = np.concatenate(all_speeds)
        threshold = float(np.percentile(all_speeds_concat, 25))
        lag_times, vacf_mean, vacf_std, dt = ensemble_vacf(per_tracer_velocities, max_lag_fraction=0.5)
        tau_c = find_zero_crossing(lag_times, vacf_mean)
        all_durations = []
        for ptv in per_tracer_velocities:
            durs = compute_residence_times(ptv['vx'], ptv['vy'], ptv['vz'], dt_val, threshold)
            all_durations.append(durs)
        all_durations_concat = np.concatenate(all_durations) if all_durations else np.array([])
        if len(all_durations_concat) >= 2:
            beta_hat, beta_se, n_used, tau_min_used = pareto_mle_beta(all_durations_concat)
        else:
            beta_hat, beta_se, n_used, tau_min_used = np.nan, np.nan, 0, np.nan
        results[beta] = {'lag_times': lag_times, 'vacf_mean': vacf_mean, 'vacf_std': vacf_std, 'tau_c': tau_c, 'threshold': threshold, 'all_durations': all_durations_concat, 'beta_trap': beta_hat, 'beta_trap_se': beta_se, 'n_residence': n_used, 'tau_min': tau_min_used, 'dt': dt_val}
    return results

if __name__ == '__main__':
    vf_data = load_and_verify_vortex(VF_PATH)
    vf_tracers = extract_tracer_trajectories_vortex(vf_data)
    vortex_results = process_vacf_and_residence_vortex(vf_tracers, [5, 10, 20, 40])
    for N, res in vortex_results.items():
        print('Vortex N=' + str(N) + ': tau_c=' + str(res['tau_c']) + ', beta_trap=' + str(res['beta_trap']))
    lw_data = load_and_verify_levy(LW_PATH)
    lw_tracers = extract_tracer_trajectories_levy(lw_data)
    levy_results = process_vacf_and_residence_levy(lw_tracers, [1.2, 1.5, 1.8, 2.5])
    for b, res in levy_results.items():
        print('Lévy beta=' + str(b) + ': tau_c=' + str(res['tau_c']) + ', beta_trap=' + str(res['beta_trap']))
    np.save(os.path.join(DATA_DIR, 'vortex_vacf_residence.npy'), vortex_results)
    np.save(os.path.join(DATA_DIR, 'levy_vacf_residence.npy'), levy_results)