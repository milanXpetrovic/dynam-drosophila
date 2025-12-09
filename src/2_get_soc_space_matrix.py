# %%
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import settings
from utils import fileio

SOC_SPACE_DISTANCE = {
    'CS_10D': 7.75,
    'Cs_5DIZ': 9.25,
    'CsCh': 8.25
}


def make_soc_matrix(df_row, treatment_soc_dist):
    n = 12
    dist_matrix = np.zeros((n, n))
    columns = df.columns

    for col_name, value in zip(columns, df_row):
        matches = re.findall(r'fly(\d+).csv', col_name)
        if len(matches) == 2:
            fly1 = int(matches[0]) - 1
            fly2 = int(matches[1]) - 1

            dist_matrix[fly1, fly2] = value
            dist_matrix[fly2, fly1] = value

    soc_space_matrix = (dist_matrix < treatment_soc_dist).astype(int)
    np.fill_diagonal(soc_space_matrix, 0)

    return soc_space_matrix


for treatment in settings.TREATMENTS:
    treatment_soc_dist = SOC_SPACE_DISTANCE[treatment]
    input_dir = settings.DISTANCE_MATRIX_DIR / treatment
    output_dir = settings.SOC_SPACE_MATRIX_DIR / treatment
    output_dir.mkdir(parents=True, exist_ok=True)
    groups = fileio.load_files_from_folder(input_dir)

    for group_name_csv, group_path in groups.items():
        group_name_csv = 'CTRL10_14_02_2024_11_04_A1.csv'
        print(group_name_csv)
        group_name = Path(group_name_csv).stem
        # group_save_path = output_dir / group_name
        # group_save_path.mkdir(parents=True, exist_ok=True)

        edgelists = pd.read_csv(group_path.replace('distance_matrix', 'edgelists'), index_col=0)
        df = pd.read_csv(group_path, index_col=0)
        soc_space_matrix = (df < treatment_soc_dist).astype(int)

        df_row = df.iloc[0].values
        m_first = make_soc_matrix(df_row, treatment_soc_dist)

        flies = [f"fly{i+1}" for i in range(m_first.shape[0])]
        m_first_df = pd.DataFrame(m_first, index=flies, columns=flies)

        # m_first_path = group_save_path / f"m_at_t0.csv"
        # m_first_df.to_csv(m_first_path)

        updates = []
        matrices = []
        for i in range(1, len(df)):
            df_row = df.iloc[i].values
            curr_matrix = make_soc_matrix(df_row, treatment_soc_dist)
            matrices.append(curr_matrix)

            diffs = np.where(curr_matrix != m_first)
            prev_matrix = curr_matrix

            # if len(updates):
            #     sys.exit()

            for fly_a, fly_b in zip(*diffs):
                if fly_a < fly_b:
                    updates.append({
                        'time': i,
                        'i': fly_a,
                        'j': fly_b,
                        'update': curr_matrix[fly_a, fly_b]
                    })
                    # updates.append({
                    #     'time': i,
                    #     'i': fly_b,
                    #     'j': fly_a,
                    #     'update': curr_matrix[fly_a, fly_b]
                    # })

        df_updateds = pd.DataFrame(updates)
        sys.exit()
# %%
# CTRL10_14_02_2024_11_04_A1
# 109,fly2,fly1,2004,2035,1

for row in df_updateds.iterrows():
    idx, data = row
    t, i, j, update = data
    if t < 109:
        m_first[i][j] = update

    if t > 109:
        sys.exit()

# %%

    # updated_m = m_first

    # %%
    # for row in edgelists.iterrows():
    # idx, data = row
    # sender, receiver, t, _, _ = data
    # sender_idx = int(sender.replace('fly', '')) - 1
    # receiver_idx = int(receiver.replace('fly', '')) - 1

    # if not matrices[t][receiver_idx][sender_idx] and matrices[t][sender_idx][receiver_idx]:
    #     print(group_name)
    #     print(data)
    #     sys.exit()

    # %%

    # sys.exit()
    # %%

    # diffs = np.where(curr_matrix != m_first)
    # prev_matrix = curr_matrix

    # for fly_a, fly_b in zip(*diffs):
    #     if fly_a < fly_b:
    #         updates.append({
    #             'time': i,
    #             'i': fly_a,
    #             'j': fly_b,
    #             'update': curr_matrix[fly_a, fly_b]
    #         })

    # df_updateds = pd.DataFrame(updates)
    # df_save_path = group_save_path / 'updates.csv'
    # df_updateds.to_csv(df_save_path, index=False)
