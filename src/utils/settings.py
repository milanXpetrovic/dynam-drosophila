from pathlib import Path

TREATMENTS = ['CS_10D', 'Cs_5DIZ', 'CsCh']

ROOT_DIR = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT_DIR / 'data'
CONFIGS_DIR = ROOT_DIR / 'configs'
EDGELISTS_DIR = DATA_DIR / 'edgelists'
COVARIANCES_DIR = DATA_DIR / 'covariances'
DISTANCE_MATRIX_DIR = DATA_DIR / 'distance_matrix'
ANGLE_MATRIX_DIR = DATA_DIR / 'angle_matrix'
SOC_SPACE_MATRIX_DIR = DATA_DIR / 'soc_space_matrix'
TRAJECTORIES_DIR = DATA_DIR / 'trajectories'
