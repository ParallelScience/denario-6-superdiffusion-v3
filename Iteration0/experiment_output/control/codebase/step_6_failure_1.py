# filename: codebase/step_6.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import matplotlib.pyplot as plt
from step_1 import compute_velocities_central_diff
from step_2 import process_dataset
from step_3 import ensemble_vacf, find_zero_crossing

DATA_DIR = 'data/'
N_VALS = [5, 10, 20, 40]
BETA_VALS = [1.2, 1.5, 1.8, 2.5]
N_COLORS = {5: '#1f77b4', 10: '#ff7f0e', 20: '#2ca02c', 40: '#d62728'}
BETA_COLORS = {1.2: '#9467bd', 1.5: '#8c564b', 1.8: '#e377c2', 2.5: '#7f7f7f'}

def plot_msd_scaling(vf_tracers, lw_tracers):
    vf_results = process_dataset(vf_tracers, 0, N_VALS, max_lag_fraction=0.5)
    lw_results = process_dataset(lw_tracers, 0, BETA_VALS, max_lag_fraction=0.5)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for N in N_VALS:
        r = vf_results[N]
        axes[0, 0].loglog(r['lag_times'], r['ensemble_tamsd'], color=N_COLORS[N], label='N=' + str(N))
        axes[1, 0].semilogx(r['lag_times'], r['slopes'], color=N_COLORS[N])
    for beta in BETA_VALS:
        r = lw_results[beta]
        axes[0, 1].loglog(r['lag_times'], r['ensemble_tamsd'], color=BETA_COLORS[beta], label='b=' + str(beta))
        axes[1, 1].semilogx(r['lag_times'], r['slopes'], color=BETA_COLORS[beta])
    plt.savefig(os.path.join(DATA_DIR, 'msd_scaling.png'))

if __name__ == '__main__':
    vf = np.load('/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy', allow_pickle=False)
    lw = np.load('/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy', allow_pickle=False)
    plot_msd_scaling(vf, lw)