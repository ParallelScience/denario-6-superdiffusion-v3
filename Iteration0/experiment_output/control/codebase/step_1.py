# filename: codebase/step_1.py
import sys
import os
sys.path.insert(0, os.path.abspath("codebase"))
sys.path.insert(0, "/home/node/data/compsep_data/")
import numpy as np
import os

def load_and_verify_vortex(path):
    vf = np.load(path, allow_pickle=False)
    print('=== VORTEX FILAMENT DATASET ===')
    print('Total rows:', vf.shape[0])
    print('Fields and dtypes:')
    for name in vf.dtype.names:
        print('  ', name, ':', vf.dtype[name])
    print('Unique trajectory_id:', np.unique(vf['trajectory_id']))
    print('Unique n_filaments:', np.unique(vf['n_filaments']))
    print('Unique gamma_std:', np.unique(vf['gamma_std']))
    print('Time range: [', vf['time'].min(), ',', vf['time'].max(), '] s')
    dt_vals = np.diff(np.unique(vf['time']))
    print('dt (from unique times):', dt_vals[0], 's')
    return vf

def load_and_verify_levy(path):
    lw = np.load(path, allow_pickle=False)
    print('\n=== LEVY WALK 3D DATASET ===')
    print('Total rows:', lw.shape[0])
    print('Fields and dtypes:')
    for name in lw.dtype.names:
        print('  ', name, ':', lw.dtype[name])
    print('Unique trajectory_id:', np.unique(lw['trajectory_id']))
    print('Unique beta:', np.unique(lw['beta']))
    print('Unique alpha_theory:', np.unique(lw['alpha_theory']))
    print('Time range: [', lw['time'].min(), ',', lw['time'].max(), '] s')
    dt_vals = np.diff(np.unique(lw['time']))
    print('dt (from unique times):', dt_vals[0], 's')
    return lw

def extract_tracer_trajectories_vortex(vf):
    tracer_dict = {}
    n_filaments_vals = np.unique(vf['n_filaments'])
    for N in n_filaments_vals:
        mask_N = vf['n_filaments'] == N
        sub = vf[mask_N]
        traj_ids = np.unique(sub['trajectory_id'])
        for tid in traj_ids:
            mask_t = sub['trajectory_id'] == tid
            rows = sub[mask_t]
            sort_idx = np.argsort(rows['time'])
            rows = rows[sort_idx]
            tracer_dict[(int(N), int(tid))] = {
                'time': rows['time'].copy(),
                'x_true': rows['x_true'].copy(),
                'y_true': rows['y_true'].copy(),
                'z_true': rows['z_true'].copy(),
                'x_noisy': rows['x_noisy'].copy(),
                'y_noisy': rows['y_noisy'].copy(),
                'z_noisy': rows['z_noisy'].copy(),
                'msd_true': rows['msd_true'].copy(),
            }
    return tracer_dict

def extract_tracer_trajectories_levy(lw):
    tracer_dict = {}
    beta_vals = np.unique(lw['beta'])
    for beta in beta_vals:
        mask_b = lw['beta'] == beta
        sub = lw[mask_b]
        traj_ids = np.unique(sub['trajectory_id'])
        for tid in traj_ids:
            mask_t = sub['trajectory_id'] == tid
            rows = sub[mask_t]
            sort_idx = np.argsort(rows['time'])
            rows = rows[sort_idx]
            tracer_dict[(float(beta), int(tid))] = {
                'time': rows['time'].copy(),
                'x_true': rows['x_true'].copy(),
                'y_true': rows['y_true'].copy(),
                'z_true': rows['z_true'].copy(),
                'x_noisy': rows['x_noisy'].copy(),
                'y_noisy': rows['y_noisy'].copy(),
                'z_noisy': rows['z_noisy'].copy(),
                'msd_true': rows['msd_true'].copy(),
                'alpha_theory': float(rows['alpha_theory'][0]),
            }
    return tracer_dict

def compute_velocities_central_diff(x, y, z, time):
    dt = time[1:] - time[:-1]
    dt_fwd = dt[1:]
    dt_bwd = dt[:-1]
    dt_sum = dt_fwd + dt_bwd
    vx = (x[2:] - x[:-2]) / dt_sum
    vy = (y[2:] - y[:-2]) / dt_sum
    vz = (z[2:] - z[:-2]) / dt_sum
    t_interior = time[1:-1]
    return vx, vy, vz, t_interior

def compute_velocity_stats_vortex(tracer_dict):
    N_vals = sorted(set(k[0] for k in tracer_dict.keys()))
    stats = {}
    for N in N_vals:
        keys_N = [(n, tid) for (n, tid) in tracer_dict.keys() if n == N]
        vx_list, vy_list, vz_list = [], [], []
        per_tracer = []
        time_range = None
        dt_val = None
        for key in keys_N:
            traj = tracer_dict[key]
            vx, vy, vz, t_int = compute_velocities_central_diff(traj['x_true'], traj['y_true'], traj['z_true'], traj['time'])
            vx_list.append(vx)
            vy_list.append(vy)
            vz_list.append(vz)
            per_tracer.append({'vx': vx, 'vy': vy, 'vz': vz, 't': t_int, 'trajectory_id': key[1]})
            if time_range is None:
                time_range = (traj['time'].min(), traj['time'].max())
                dt_val = traj['time'][1] - traj['time'][0]
        vx_all = np.concatenate(vx_list)
        vy_all = np.concatenate(vy_list)
        vz_all = np.concatenate(vz_list)
        speeds = np.sqrt(vx_all**2 + vy_all**2 + vz_all**2)
        sigma_v2 = np.var(vx_all) + np.var(vy_all) + np.var(vz_all)
        mean_speed = np.mean(speeds)
        stats[N] = {'sigma_v2': sigma_v2, 'mean_speed': mean_speed, 'speeds': speeds, 'vx_all': vx_all, 'vy_all': vy_all, 'vz_all': vz_all, 'n_tracers': len(keys_N), 'time_range': time_range, 'dt': dt_val, 'per_tracer_velocities': per_tracer}
    return stats

def compute_velocity_stats_levy(tracer_dict):
    beta_vals = sorted(set(k[0] for k in tracer_dict.keys()))
    stats = {}
    for beta in beta_vals:
        keys_b = [(b, tid) for (b, tid) in tracer_dict.keys() if b == beta]
        vx_list, vy_list, vz_list = [], [], []
        per_tracer = []
        time_range = None
        dt_val = None
        alpha_theory = None
        for key in keys_b:
            traj = tracer_dict[key]
            vx, vy, vz, t_int = compute_velocities_central_diff(traj['x_true'], traj['y_true'], traj['z_true'], traj['time'])
            vx_list.append(vx)
            vy_list.append(vy)
            vz_list.append(vz)
            per_tracer.append({'vx': vx, 'vy': vy, 'vz': vz, 't': t_int, 'trajectory_id': key[1]})
            if time_range is None:
                time_range = (traj['time'].min(), traj['time'].max())
                dt_val = traj['time'][1] - traj['time'][0]
                alpha_theory = traj['alpha_theory']
        vx_all = np.concatenate(vx_list)
        vy_all = np.concatenate(vy_list)
        vz_all = np.concatenate(vz_list)
        speeds = np.sqrt(vx_all**2 + vy_all**2 + vz_all**2)
        sigma_v2 = np.var(vx_all) + np.var(vy_all) + np.var(vz_all)
        mean_speed = np.mean(speeds)
        stats[beta] = {'sigma_v2': sigma_v2, 'mean_speed': mean_speed, 'speeds': speeds, 'vx_all': vx_all, 'vy_all': vy_all, 'vz_all': vz_all, 'n_tracers': len(keys_b), 'time_range': time_range, 'dt': dt_val, 'alpha_theory': alpha_theory, 'per_tracer_velocities': per_tracer}
    return stats

if __name__ == '__main__':
    vf_data = load_and_verify_vortex('/home/node/work/projects/superdiffusion_v3/vortex_filament_tracers.npy')
    lw_data = load_and_verify_levy('/home/node/work/projects/superdiffusion_v3/levy_walk_3d_trajectories.npy')
    vf_tracers = extract_tracer_trajectories_vortex(vf_data)
    lw_tracers = extract_tracer_trajectories_levy(lw_data)
    vf_stats = compute_velocity_stats_vortex(vf_tracers)
    lw_stats = compute_velocity_stats_levy(lw_tracers)
    np.savez('data/vortex_stats.npz', stats=vf_stats)
    np.savez('data/levy_stats.npz', stats=lw_stats)
    print('Saved velocity statistics to data/vortex_stats.npz and data/levy_stats.npz')