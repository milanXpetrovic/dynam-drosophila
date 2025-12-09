# %%
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from natsort import natsort_keygen, natsorted

import settings
from utils import fileio

# ego(nodesAttr$negative_inf_weighted) + #! distances between all flies matrix
# ego(nodesAttr$negative_inf_weighted) + #! number of uniques fly met between interaction
# ego(nodesAttr$negative_inf_weighted) + #! difference in distance traveled between interactions

TIME_WINDOW_SIZE_SEC = 36
NUM_FLIES = 12
MAX_TIME = 28800
TIME_WINDOW = 24 * TIME_WINDOW_SIZE_SEC
FLY_IDS = range(1, NUM_FLIES + 1)

SOC_SPACE_DISTANCE = {
    'CS_10D': 7.75,
    'Cs_5DIZ': 9.25,
    'CsCh': 8.25
}

INTERACTION_SPACE_DISTANCE = {
    'CS_10D': 2.0,
    'Cs_5DIZ': 2.5,
    'CsCh': 2.25
}


def number_of_flies_in_soc_space(df):
    """
    Computes the number of neighboring flies within social space per time step and saves the result in long format.
    """

    all_flies = natsorted(pd.unique([x for col in df.columns for x in col.split()]))
    fly_index = {fly: i for i, fly in enumerate(all_flies)}

    df = (df < SOC_SPACE_DISTANCE[treatment]).astype(int)
    df_long = df.melt(var_name='pair', value_name='value')
    df_long[['fly_a', 'fly_b']] = df_long['pair'].str.split(' ', expand=True)

    index_pairs = [tuple(col.split()) for col in df.columns]
    row_col_indices = [(fly_index[a], fly_index[b]) for a, b in index_pairs]

    matrices = np.zeros((df.shape[0], len(all_flies), len(all_flies)))

    for col_idx, (i, j) in enumerate(row_col_indices):
        matrices[:, i, j] = df.iloc[:, col_idx].values

    clean_flies = [fly.replace('.csv', '') for fly in all_flies]
    clean_flies = natsorted(clean_flies)

    df_neighbors = pd.DataFrame(matrices.sum(axis=2).astype(int), columns=clean_flies)
    df_neighbors['time'] = df_neighbors.index

    df_long = df_neighbors.melt(id_vars='time', var_name='actor', value_name='value')
    natsort_key = natsort_keygen()
    df_long = df_long.sort_values(by=['time', 'actor'], key=lambda col: col.map(
        natsort_key) if col.name == 'actor' else col)
    df_long = df_long.reset_index(drop=True)

    save_path = output_dir / group_name / "number_of_flies_in_soc_space.csv"
    df_long.to_csv(save_path)


def get_unique_flies(df_distances, df_edgelist, threshold):
    df_distances = (df_distances <= threshold).astype(int)

    df_long = df_distances.melt(var_name='pair', value_name='value')
    df_long[['fly_a', 'fly_b']] = df_long['pair'].str.split(' ', expand=True)

    all_flies = natsorted(pd.unique([x for col in df_distances.columns for x in col.split()]))
    fly_index = {fly: i for i, fly in enumerate(all_flies)}
    index_pairs = [tuple(col.split()) for col in df_distances.columns]
    row_col_indices = [(fly_index[a], fly_index[b]) for a, b in index_pairs]
    all_flies = [fly.replace('.csv', '') for fly in all_flies]

    matrices = np.zeros((df_distances.shape[0], len(all_flies), len(all_flies)))
    for col_idx, (i, j) in enumerate(row_col_indices):
        matrices[:, i, j] = df_distances.iloc[:, col_idx].values

    grouped = df_edgelist.groupby('sender')

    res = []
    for sender, group in grouped:
        sender_number = int(re.search(r'\d+', sender).group())
        times = group['time'].sort_values().to_numpy()

        res.append({"time": times[0], "actor": sender, "value": 0})
        for i in range(len(times) - 1):
            t1, t2 = times[i], times[i + 1]

            sliced_matrices = matrices[t1:t2]  # ! we take time - 1, i.e. count unique before new interaction
            summed_matrix = sliced_matrices.sum(axis=0)
            summed_matrix = (summed_matrix > 0).astype(int)

            value = sum(summed_matrix[sender_number-1])
            res.append({"time": t2, "actor": sender, "value": value})

    df_results = pd.DataFrame(res)
    df_results = df_results.sort_values('time').reset_index(drop=True)

    return df_results


def distance_traveled_between_interactions(df_edgelist, group_trajectories):
    grouped = df_edgelist.groupby('sender')
    res = []

    for sender, group in grouped:
        trajectory = pd.read_csv(group_trajectories[f'{sender}.csv'], usecols=['pos x', 'pos y'])
        times = group['time'].sort_values().to_numpy()

        res.append({"time": times[0], "actor": sender, "value": 0})

        pos_x = trajectory['pos x'].to_numpy()
        pos_y = trajectory['pos y'].to_numpy()
        dx = np.diff(pos_x)
        dy = np.diff(pos_y)
        step_distances = np.sqrt(dx**2 + dy**2)
        step_distances = np.insert(step_distances, 0, 0)
        cumulative_distance = np.cumsum(step_distances)

        for i in range(len(times) - 1):
            t1, t2 = times[i], times[i + 1]

            dist = np.round(cumulative_distance[t2] - cumulative_distance[t1], 3)
            res.append({"time": t2, "actor": sender, "value": dist})

    df_results = pd.DataFrame(res)
    df_results = df_results.sort_values('time').reset_index(drop=True)

    return df_results


for treatment in settings.TREATMENTS:
    print(treatment)
    input_distances_dir = settings.DISTANCE_MATRIX_DIR / treatment
    input_edgelists_dir = settings.EDGELISTS_DIR / treatment
    input_traj_dir = settings.TRAJECTORIES_DIR / treatment

    output_dir = settings.COVARIANCES_DIR / treatment
    output_dir.mkdir(parents=True, exist_ok=True)

    groups_distances = fileio.load_files_from_folder(input_distances_dir)
    groups_edgelists = fileio.load_files_from_folder(input_edgelists_dir)
    all_trajectories = fileio.load_multiple_folders(input_traj_dir)

    for (group_name_traj, group_path_traj), (_, group_path_edges) in zip(all_trajectories.items(), groups_edgelists.items()):
        group_name = Path(group_name_traj).stem
        print(group_name)

        df_edgelist = pd.read_csv(group_path_edges, index_col=0)
        group_trajectories = fileio.load_files_from_folder(all_trajectories[group_name])

        df_results = distance_traveled_between_interactions(df_edgelist, group_trajectories)
        save_path = output_dir / group_name / "distance_traveled_between_interactions.csv"
        df_results.to_csv(save_path)

    for (group_name_dist, group_path_dist), (group_name_edges, group_path_edges) in zip(groups_distances.items(), groups_edgelists.items()):
        group_name = Path(group_name_dist).stem
        print(group_name)

        df_distances = pd.read_csv(group_path_dist, index_col=0)
        df_edgelist = pd.read_csv(group_path_edges, index_col=0)

        number_of_flies_in_soc_space(df_distances)

        threshold = SOC_SPACE_DISTANCE[treatment]
        df_results = get_unique_flies(df_distances, df_edgelist, threshold)
        save_path = output_dir / group_name / "unique_partners_met_social_space.csv"
        df_results.to_csv(save_path)

        threshold = INTERACTION_SPACE_DISTANCE[treatment]
        df_results = get_unique_flies(df_distances, df_edgelist, threshold)
        save_path = output_dir / group_name / "unique_partners_met_interaction_space.csv"
        df_results.to_csv(save_path)
