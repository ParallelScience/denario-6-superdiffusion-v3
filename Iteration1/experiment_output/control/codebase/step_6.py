# filename: codebase/step_6.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from step_1 import load_and_verify_vortex, compute_velocities_central_diff

def compute_residence_times_for_N(tracer_dict, N, dt=0.05, n_bins=20):
    keys_N = [(n, tid) for (n, tid) in tracer_dict.keys() if n == N]
    all_speeds = []
    for key in keys_N:
        traj = tracer_dict[key]
        vx, vy, vz, _ = compute_velocities_central_diff(traj['x_true'], traj['y_true'], traj['z_true'], traj['time'])
        speed = np.sqrt(vx**2 + vy**2 + vz**2)
        all_speeds.append(speed)
    all_speeds_arr = np.concatenate(all_speeds)
    threshold = float(np.percentile(all_speeds_arr, 25))
    all_durations = []
    for speed in all_speeds:
        below = speed < threshold
        count = 0
        for val in below:
            if val:
                count += 1
            else:
                if count > 0:
                    all_durations.append(count * dt)
                    count = 0
        if count > 0:
            all_durations.append(count * dt)
    all_dur_arr = np.array(all_durations, dtype=float)
    if len(all_dur_arr) < 2:
        return {'bin_centers': np.array([dt]), 'hist': np.array([0.0])}
    bins = np.logspace(np.log10(max(all_dur_arr.min(), dt)), np.log10(all_dur_arr.max()), n_bins + 1)
    hist, edges = np.histogram(all_dur_arr, bins=bins, density=True)
    return {'bin_centers': np.sqrt(edges[:-1] * edges[1:]), 'hist': hist}

def plot_residence_times(tracer_dict):
    fig, ax = plt.subplots(figsize=(6, 5))
    N_vals = [5, 10, 20, 40]
    for N in N_vals:
        res = compute_residence_times_for_N(tracer_dict, N)
        ax.loglog(res['bin_centers'], res['hist'], 'o-', label='N='+str(N))
    tau = np.logspace(-1, 1, 10)
    ax.loglog(tau, 0.1 * tau**(-2.5), 'k--', label='Holtsmark (tau^-2.5)')
    ax.set_xlabel('Residence time (s)')
    ax.set_ylabel('P(tau)')
    ax.legend()
    plt.savefig('data/plot_5.png', dpi=300)
    plt.close()

if __name__ == '__main__':
    vf_path = '/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy'
    vf_data = load_and_verify_vortex(vf_path)
    tracer_dict = {}
    for row in vf_data:
        key = (int(row['n_filaments']), int(row['trajectory_id']))
        if key not in tracer_dict:
            tracer_dict[key] = {'x_true': [], 'y_true': [], 'z_true': [], 'time': []}
        tracer_dict[key]['x_true'].append(row['x_true'])
        tracer_dict[key]['y_true'].append(row['y_true'])
        tracer_dict[key]['z_true'].append(row['z_true'])
        tracer_dict[key]['time'].append(row['time'])
    for key in tracer_dict:
        for k in ['x_true', 'y_true', 'z_true', 'time']:
            tracer_dict[key][k] = np.array(tracer_dict[key][k])
    plot_residence_times(tracer_dict)