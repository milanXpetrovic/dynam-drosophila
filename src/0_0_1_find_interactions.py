# %%
import os
import sys

import numpy as np
import pandas as pd
import toml

import settings
from utils import fileio

VIDEO_FPS = 24

for treatment in settings.TREATMENTS:
    SCRIPT_OUTPUT = settings.EDGELISTS_DIR / treatment
    SCRIPT_OUTPUT.mkdir(parents=True, exist_ok=True)

    input_distances_dir = settings.DISTANCE_MATRIX_DIR / treatment
    input_angles_dir = settings.ANGLE_MATRIX_DIR / treatment

    groups_distances = fileio.load_files_from_folder(input_distances_dir)
    groups_angles = fileio.load_files_from_folder(input_angles_dir)

    TREATMENT_CONFIG = settings.CONFIGS_DIR / "interaction_criteria" / f"{treatment}.toml"
    with TREATMENT_CONFIG.open() as f:
        treatment_config = toml.load(f)

    ANGLE = treatment_config["ANGLE"]
    DISTANCE = treatment_config["DISTANCE"]
    TIME = treatment_config["TIME"]

    for angles_tuple, distances_tuple in zip(groups_angles.items(), groups_distances.items()):
        angles_name, angles_path = angles_tuple
        print(angles_name)
        distances_name, distances_path = distances_tuple

        if angles_name != distances_name:
            sys.exit()

        df_angles = pd.read_csv(angles_path, index_col=0)
        df_distances = pd.read_csv(distances_path, index_col=0)

        edgelist = pd.DataFrame(
            columns=[
                "sender",
                "receiver",
                "time",
                "end_time",
            ])

        for angles_col, distances_col in zip(df_angles.columns, df_distances.columns):
            if angles_col != distances_col:
                sys.exit()

            df = pd.concat([df_angles[angles_col], df_distances[distances_col]], axis=1)
            df.columns = ["angle", "distance"]

            distance_mask = df["distance"] <= DISTANCE
            angle_mask = (df["angle"] >= ANGLE[0]) & (df["angle"] <= ANGLE[1])
            df = df[distance_mask & angle_mask]

            timecut_frames = int(TIME * VIDEO_FPS)
            clear_list_of_df = [d for _, d in df.groupby(df.index - np.arange(len(df))) if len(d) >= timecut_frames]

            node_1, node_2 = angles_col.split(" ")
            node_1, node_2 = node_1.replace(".csv", ""), node_2.replace(".csv", "")

            for interaction in clear_list_of_df:
                data = {
                    "sender": node_1,
                    "receiver": node_2,
                    "time": int(interaction.index[0]),
                    "end_time": int(interaction.index[-1]),
                }

                row = pd.DataFrame.from_dict(data, orient="index").T
                edgelist = pd.concat([edgelist, row], ignore_index=True)

        edgelist = edgelist.sort_values("time").reset_index(drop=True)
        edgelist.to_csv(os.path.join(SCRIPT_OUTPUT, angles_name))
