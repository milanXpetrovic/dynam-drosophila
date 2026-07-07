# %%
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

        while df['time'].duplicated().any():
            duplicated_mask = df['time'].duplicated(keep='first')
            df.loc[duplicated_mask, 'time'] += 1

        save_name = treatment_dir / group_name
        df.to_csv(save_name)
