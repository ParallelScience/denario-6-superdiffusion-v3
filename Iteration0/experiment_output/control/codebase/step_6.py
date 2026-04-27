# filename: codebase/step_6.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
from step_1 import extract_tracer_trajectories_vortex, extract_tracer_trajectories_levy
from step_2 import process_dataset

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
    axes[0, 0].set_xlabel('Lag time (s)')
    axes[0, 0].set_ylabel('TAMSD (m^2)')
    axes[0, 0].set_title('Vortex Filaments: TAMSD')
    axes[0, 0].legend()
    axes[1, 0].set_xlabel('Lag time (s)')
    axes[1, 0].set_ylabel('Local slope')
    axes[1, 0].set_title('Vortex Filaments: Local Slope')
    for beta in BETA_VALS:
        r = lw_results[beta]
        axes[0, 1].loglog(r['lag_times'], r['ensemble_tamsd'], color=BETA_COLORS[beta], label='β=' + str(beta))
        axes[1, 1].semilogx(r['lag_times'], r['slopes'], color=BETA_COLORS[beta])
    axes[0, 1].set_xlabel('Lag time (s)')
    axes[0, 1].set_ylabel('TAMSD (m^2)')
    axes[0, 1].set_title('Lévy Walk: TAMSD')
    axes[0, 1].legend()
    axes[1, 1].set_xlabel('Lag time (s)')
    axes[1, 1].set_ylabel('Local slope')
    axes[1, 1].set_title('Lévy Walk: Local Slope')
    plt.tight_layout()
    plt.savefig(os.path.join(DATA_DIR, 'msd_scaling.png'))

if __name__ == '__main__':
    vf_raw = np.load('/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy', allow_pickle=False)
    lw_raw = np.load('/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy', allow_pickle=False)
    vf_tracers = extract_tracer_trajectories_vortex(vf_raw)
    lw_tracers = extract_tracer_trajectories_levy(lw_raw)
    plot_msd_scaling(vf_tracers, lw_tracers)