# %%
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

import settings
from utils import fileio

# --- Constants ---
TIME_WINDOW_SIZE_SEC = 36 # 3, 12
NUM_FLIES = 12
MAX_TIME = 28800
TIME_WINDOW = 24 * TIME_WINDOW_SIZE_SEC
FLY_IDS = range(1, NUM_FLIES + 1)


def generate_column_names(fly_ids):
    """Generates various column names based on fly IDs."""
    cols = {}
    cols['in_deg'] = [f'fly{i}_in_degree' for i in fly_ids]
    cols['out_deg'] = [f'fly{i}_out_degree' for i in fly_ids]
    cols['in_weighted'] = [f'fly{i}_in_weighted' for i in fly_ids]
    cols['out_weighted'] = [f'fly{i}_out_weighted' for i in fly_ids]
    cols['influence'] = [f'fly{i}_influence' for i in fly_ids]
    cols['inf_weighted'] = [f'fly{i}_inf_weighted' for i in fly_ids]
    cols['all_degrees'] = cols['in_deg'] + cols['out_deg'] + cols['in_weighted'] + cols['out_weighted']
    return cols


def calculate_degrees(df, time_window, max_time, all_degree_cols):
    """
    Calculates degrees for each interaction.
    This version fixes the bug where values were being incremented (+=) instead of assigned (=).
    """
    res = pd.DataFrame(0, index=range(max_time), columns=all_degree_cols)
    times = df['time'].to_numpy()

    for idx, row in df.iterrows():
        interaction_time = row['time']
        sender = row['sender']
        receiver = row['receiver']

        start_time = interaction_time - time_window
        start_idx = np.searchsorted(times, start_time, side='left')

        df_time_window = df.iloc[start_idx: idx + 1]

        if df_time_window.empty:
            continue

        out_degree_w = df_time_window['sender'].value_counts()
        in_degree_w = df_time_window['receiver'].value_counts()
        res.loc[interaction_time, f'{receiver}_in_weighted'] = in_degree_w.get(receiver, 0)
        res.loc[interaction_time, f'{sender}_out_weighted'] = out_degree_w.get(sender, 0)

        unique_interactions = df_time_window.drop_duplicates(subset=['sender', 'receiver'])
        out_degree = unique_interactions['sender'].value_counts()
        in_degree = unique_interactions['receiver'].value_counts()
        res.loc[interaction_time, f'{receiver}_in_degree'] = in_degree.get(receiver, 0)
        res.loc[interaction_time, f'{sender}_out_degree'] = out_degree.get(sender, 0)

    last_non_zero_indices = res.apply(lambda col: col[col > 0].last_valid_index()).fillna(0).astype(int)
    res = res.replace(0, np.nan).ffill().fillna(0).astype(int)

    for col in res.columns:
        last_index = last_non_zero_indices[col]
        cutoff_index = last_index + time_window
        if cutoff_index < max_time:
            res.loc[cutoff_index + 1:, col] = 0

    return res


def calculate_influence(df, degrees_df, time_window, influence_cols, inf_weighted_cols):
    """Calculates influence scores based on pre-computed degrees."""
    influence_df = pd.DataFrame(0, index=degrees_df.index, columns=influence_cols + inf_weighted_cols)

    for time_val, sender, receiver in zip(df['time'], df['sender'], df['receiver']):
        out_deg_sender = degrees_df.at[time_val, f'{sender}_out_degree']
        out_deg_receiver = degrees_df.at[time_val, f'{receiver}_out_degree']
        out_deg_sender_w = degrees_df.at[time_val, f'{sender}_out_weighted']
        out_deg_receiver_w = degrees_df.at[time_val, f'{receiver}_out_weighted']

        influence_df.at[time_val, f'{receiver}_influence'] = out_deg_receiver - out_deg_sender
        influence_df.at[time_val, f'{receiver}_inf_weighted'] = out_deg_receiver_w - out_deg_sender_w

    influence_tw = influence_df.replace(0, np.nan).rolling(window=time_window, min_periods=1).mean()
    return influence_tw.fillna(0)


def calculate_burstiness(degrees_df, time_window, out_degree_cols, out_weighted_cols):
    """
    Calculates burstiness metrics.
    This version fixes the ValueError by using a dictionary comprehension to build the DataFrame.
    """
    rolling_sum_df = degrees_df.rolling(window=time_window, min_periods=1).sum()
    res_burstines_unw = rolling_sum_df.filter(like='_out_degree')
    max_rolling = res_burstines_unw.rolling(window=time_window, min_periods=1).max().max(axis=1)

    burstiness_max_values = pd.DataFrame({col: max_rolling for col in out_degree_cols},
                                         index=res_burstines_unw.index)

    res_burstines_weighted = rolling_sum_df.filter(like='_out_weighted')
    max_weight_rolling = res_burstines_weighted.rolling(window=time_window, min_periods=1).max().max(axis=1)

    burstiness_weighted_max_values = pd.DataFrame({col: max_weight_rolling for col in out_weighted_cols},
                                                  index=res_burstines_weighted.index)

    return burstiness_max_values, burstiness_weighted_max_values


def melt_and_save(df, value_suffix, output_dir, file_prefix=""):
    """Melts a DataFrame to long format and saves it to a CSV."""
    if df.empty:
        return

    df_reset = df.reset_index().rename(columns={'index': 'time'})
    df_long = df_reset.melt(id_vars=['time'], var_name='actor', value_name='value')
    df_long['actor'] = df_long['actor'].str.replace(value_suffix, '', regex=False)
    df_long = df_long.sort_values(by=['time', 'actor']).reset_index(drop=True)

    save_path = output_dir / f'{file_prefix}{value_suffix[1:]}.csv'
    df_long.to_csv(save_path, index=False)
    print(f"Saved {save_path}")


def process_group(group_name, group_path, output_dir_base):
    """Main processing pipeline for a single group file."""
    print(f"\n--- Processing group: {group_name} ---")
    start_time = time.time()

    group_output_dir = output_dir_base / group_name
    group_output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(group_path, index_col=0)
    df = df.sort_values('time').reset_index(drop=True)

    col_names = generate_column_names(FLY_IDS)
    degrees_df = calculate_degrees(df, TIME_WINDOW, MAX_TIME, col_names['all_degrees'])
    influence_df = calculate_influence(df, degrees_df, TIME_WINDOW, col_names['influence'], col_names['inf_weighted'])
    burstiness_df, burstiness_w_df = calculate_burstiness(
        degrees_df, TIME_WINDOW, col_names['out_deg'], col_names['out_weighted'])

    print("Saving all results...")
    for cov_value in ['_out_degree', '_out_weighted', '_in_degree', '_in_weighted']:
        melt_and_save(degrees_df.filter(like=cov_value), cov_value, group_output_dir)

    melt_and_save(burstiness_df, '_out_degree', group_output_dir, file_prefix="burstines_")
    melt_and_save(burstiness_w_df, '_out_weighted', group_output_dir, file_prefix="burstines_")

    for cov_value in ['_influence', '_inf_weighted']:
        df_to_save = influence_df.filter(like=cov_value)
        positive_df = df_to_save.where(df_to_save > 0, 0)
        melt_and_save(positive_df, cov_value, group_output_dir, file_prefix="positive")
        negative_df = df_to_save.where(df_to_save < 0, 0)
        melt_and_save(negative_df, cov_value, group_output_dir, file_prefix="negative")

    end_time = time.time()
    print(f"--- Finished {group_name} in {end_time - start_time:.2f} seconds ---")


def main():
    """Main script execution."""
    for treatment in settings.TREATMENTS:
        input_dir = settings.EDGELISTS_DIR / treatment
        output_dir = settings.COVARIANCES_DIR / treatment
        output_dir.mkdir(parents=True, exist_ok=True)

        groups = fileio.load_files_from_folder(input_dir)

        for group_name_csv, group_path in groups.items():
            group_name = Path(group_name_csv).stem
            process_group(group_name, group_path, output_dir)


if __name__ == '__main__':
    main()
