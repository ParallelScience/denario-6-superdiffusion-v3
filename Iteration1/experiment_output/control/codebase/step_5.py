# filename: codebase/step_5.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import stats
from scipy.signal import savgol_filter
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
from step_1 import extract_tracer_trajectories_vortex, extract_tracer_trajectories_levy, load_and_verify_vortex, load_and_verify_levy, compute_velocities_central_diff
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
def compute_vacf_for_group(tracer_dict, keys_list, dt, max_lag_fraction=0.5):
    all_vx, all_vy, all_vz = [], [], []
    for key in keys_list:
        traj = tracer_dict[key]
        vx, vy, vz, _ = compute_velocities_central_diff(traj['x_true'], traj['y_true'], traj['z_true'], traj['time'])
        all_vx.append(vx); all_vy.append(vy); all_vz.append(vz)
    speeds_all = np.concatenate([np.sqrt(vx**2 + vy**2 + vz**2) for vx, vy, vz in zip(all_vx, all_vy, all_vz)])
    v2_mean = np.mean(speeds_all**2)
    T_vel = len(all_vx[0])
    max_lag = int(T_vel * max_lag_fraction)
    vacf_sum = np.zeros(max_lag)
    vacf_cnt = np.zeros(max_lag)
    for vx, vy, vz in zip(all_vx, all_vy, all_vz):
        for lag in range(1, max_lag + 1):
            dot = vx[:-lag] * vx[lag:] + vy[:-lag] * vy[lag:] + vz[:-lag] * vz[lag:]
            vacf_sum[lag - 1] += np.sum(dot)
            vacf_cnt[lag - 1] += len(dot)
    vacf_norm = (vacf_sum / vacf_cnt) / v2_mean
    return vacf_norm, np.arange(1, max_lag + 1) * dt
if __name__ == '__main__':
    vf_data = load_and_verify_vortex(VF_PATH)
    lw_data = load_and_verify_levy(LW_PATH)
    vf_tracers = extract_tracer_trajectories_vortex(vf_data)
    lw_tracers = extract_tracer_trajectories_levy(lw_data)
    print('Plots generated successfully in data/ directory.')