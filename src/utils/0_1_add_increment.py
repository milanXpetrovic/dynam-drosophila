# %%
import sys

import pandas as pd

import settings
from utils import fileio

for TREATMENT in settings.TREATMENTS:

    treatment_dir = settings.EDGELISTS_DIR / TREATMENT
    treatment_dir.mkdir(parents=True, exist_ok=True)

    group = fileio.load_files_from_folder(treatment_dir)

    for group_name, group_path in group.items():
        print(group_name)
        df = pd.read_csv(group_path, index_col=0)
        df['increment'] = 1

        df.to_csv(group_path)
