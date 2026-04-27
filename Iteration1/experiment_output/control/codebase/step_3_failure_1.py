# filename: codebase/step_3.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import stats
from scipy.ndimage import uniform_filter1d
import os
from step_1 import extract_tracer_trajectories_vortex, load_and_verify_vortex, compute_velocities_central_diff
VF_PATH = '/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy'
DATA_DIR = 'data/'
def compute_speed_per_tracer(tracer_dict, N, dt=0.05):
    keys_N = [(n, tid) for (n, tid) in tracer_dict.keys() if n == N]
    results = []
    for key in keys_N:
        traj = tracer_dict[key]
        vx, vy, vz, t_int = compute_velocities_central_diff(traj['x_true'], traj['y_true'], traj['z_true'], traj['time'])
        speed = np.sqrt(vx**2 + vy**2 + vz**2)
        results.append({'speed': speed, 'vx': vx, 'vy': vy, 'vz': vz, 't_interior': t_int, 'trajectory_id': key[1], 'x_true': traj['x_true'][1:-1], 'y_true': traj['y_true'][1:-1], 'z_true': traj['z_true'][1:-1]})
    return results
def compute_residence_times(speed, threshold):
    below = speed < threshold
    durations = []
    count = 0
    for val in below:
        if val:
            count += 1
        else:
            if count > 0:
                durations.append(count)
                count = 0
    if count > 0:
        durations.append(count)
    return np.array(durations, dtype=float)
def compute_residence_time_distributions(tracer_list, dt, n_bins=30):
    all_speeds = np.concatenate([t['speed'] for t in tracer_list])
    median_speed = float(np.median(all_speeds))
    threshold = median_speed / 2.0
    all_durations = []
    for t in tracer_list:
        dur = compute_residence_times(t['speed'], threshold)
        all_durations.append(dur * dt)
    all_dur_arr = np.concatenate(all_durations) if all_durations else np.array([])
    if len(all_dur_arr) < 2:
        return {'bin_centers': np.array([]), 'hist': np.array([]), 'threshold': threshold, 'all_durations': all_dur_arr, 'median_speed': median_speed, 'n_events': 0}
    min_dur = all_dur_arr.min()
    max_dur = all_dur_arr.max()
    if min_dur <= 0:
        min_dur = dt
    bins = np.logspace(np.log10(min_dur), np.log10(max_dur), n_bins + 1)
    hist, edges = np.histogram(all_dur_arr, bins=bins, density=True)
    bin_centers = np.sqrt(edges[:-1] * edges[1:])
    return {'bin_centers': bin_centers, 'hist': hist, 'threshold': threshold, 'all_durations': all_dur_arr, 'median_speed': median_speed, 'n_events': len(all_dur_arr)}
def compute_displacement_pdfs_N40(tracer_dict, lag_times_s, dt=0.05, n_bins=80):
    N = 40
    keys_N = [(n, tid) for (n, tid) in tracer_dict.keys() if n == N]
    results = {}
    for lag_s in lag_times_s:
        lag_steps = int(round(lag_s / dt))
        delta_x_all = []
        for key in keys_N:
            traj = tracer_dict[key]
            x = traj['x_true']
            if len(x) > lag_steps:
                dx = x[lag_steps:] - x[:-lag_steps]
                delta_x_all.append(dx)
        if not delta_x_all:
            results[lag_s] = {'bin_centers': np.array([]), 'pdf': np.array([]), 'delta_x': np.array([])}
            continue
        dx_arr = np.concatenate(delta_x_all)
        std_dx = float(np.std(dx_arr))
        kurt = float(stats.kurtosis(dx_arr, fisher=True))
        abs_max = np.percentile(np.abs(dx_arr), 99.5)
        bins = np.linspace(-abs_max, abs_max, n_bins + 1)
        hist, edges = np.histogram(dx_arr, bins=bins, density=True)
        bin_centers = 0.5 * (edges[:-1] + edges[1:])
        results[lag_s] = {'bin_centers': bin_centers, 'pdf': hist, 'delta_x': dx_arr, 'std_dx': std_dx, 'kurtosis': kurt, 'n_samples': len(dx_arr)}
    return results
def compute_okubo_weiss_proxy(tracer_list, dt):
    ow_all = []
    speed_all = []
    strain_all = []
    rot_all = []
    for t in tracer_list:
        vx, vy, vz, speed = t['vx'], t['vy'], t['vz'], t['speed']
        n = len(vx)
        if n < 3: continue
        dvx, dvy, dvz = np.gradient(vx, dt), np.gradient(vy, dt), np.gradient(vz, dt)
        strain_sq = dvx**2 + dvy**2 + dvz**2
        cross_x, cross_y, cross_z = vy * dvz - vz * dvy, vz * dvx - vx * dvz, vx * dvy - vy * dvx
        cross_mag_sq = cross_x**2 + cross_y**2 + cross_z**2
        speed_sq = speed**2
        safe_speed_sq = np.where(speed_sq > 1e-30, speed_sq, 1e-30)
        rot_sq = cross_mag_sq / safe_speed_sq**2
        ow = strain_sq - rot_sq
        ow_all.append(ow)
        speed_all.append(speed)
        strain_all.append(np.sqrt(strain_sq))
        rot_all.append(np.sqrt(rot_sq))
    ow_arr, speed_arr, strain_arr, rot_arr = np.concatenate(ow_all), np.concatenate(speed_all), np.concatenate(strain_all), np.concatenate(rot_all)
    finite_mask = np.isfinite(ow_arr) & np.isfinite(speed_arr)
    r, p = stats.pearsonr(speed_arr[finite_mask], ow_arr[finite_mask])
    return {'ow_proxy': ow_arr[finite_mask], 'speed_interior': speed_arr[finite_mask], 'pearson_r': float(r), 'pearson_p': float(p), 'strain_rate': strain_arr[finite_mask], 'rotation_rate': rot_arr[finite_mask]}
def compute_rolling_speed_variance_correlation(tracer_list, window=10):
    all_speed, all_rolling_var = [], []
    for t in tracer_list:
        speed = t['speed']
        n = len(speed)
        if n < window + 1: continue
        speed_sq = speed**2
        rolling_mean_sq = uniform_filter1d(speed_sq, size=window, mode='reflect')
        rolling_mean = uniform_filter1d(speed, size=window, mode='reflect')
        rolling_var = np.maximum(rolling_mean_sq - rolling_mean**2, 0.0)
        half = window // 2
        valid = slice(half, n - half)
        all_speed.append(speed[valid])
        all_rolling_var.append(rolling_var[valid])
    if not all_speed: return {'pearson_r': np.nan, 'pearson_p': np.nan}
    speed_arr, var_arr = np.concatenate(all_speed), np.concatenate(all_rolling_var)
    finite_mask = np.isfinite(speed_arr) & np.isfinite(var_arr)
    r, p = stats.pearsonr(speed_arr[finite_mask], var_arr[finite_mask])
    return {'pearson_r': float(r), 'pearson_p': float(p), 'all_speed': speed_arr[finite_mask], 'all_rolling_var': var_arr[finite_mask]}
if __name__ == '__main__':
    vf = load_and_verify_vortex(VF_PATH)
    tracer_dict = extract_tracer_trajectories_vortex(vf)
    for N in [5, 10, 20, 40]:
        tracers = compute_speed_per_tracer(tracer_dict, N)
        res_dist = compute_residence_time_distributions(tracers, 0.05)
        np.savez(os.path.join(DATA_DIR, 'residence_times_N' + str(N) + '.npz'), **res_dist)
        print('Saved residence times for N=' + str(N))
    pdfs = compute_displacement_pdfs_N40(tracer_dict, [1.0, 5.0, 10.0, 20.0])
    np.savez(os.path.join(DATA_DIR, 'displacement_pdfs_N40.npz'), **pdfs)
    print('Saved displacement PDFs for N=40')
    tracers_40 = compute_speed_per_tracer(tracer_dict, 40)
    ow_results = compute_okubo_weiss_proxy(tracers_40, 0.05)
    np.savez(os.path.join(DATA_DIR, 'ow_results_N40.npz'), **ow_results)
    print('Saved OW results for N=40. Pearson r: ' + str(ow_results['pearson_r']))
    corr_results = compute_rolling_speed_variance_correlation(tracers_40)
    np.savez(os.path.join(DATA_DIR, 'spatial_corr_N40.npz'), **corr_results)
    print('Saved spatial correlation for N=40. Pearson r: ' + str(corr_results['pearson_r']))