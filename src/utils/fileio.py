import os
import re
import sys


def natural_sort(l):
    def convert(text):
        return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key):
        return [convert(c) for c in re.split("([0-9]+)", key)]

    return sorted(l, key=alphanum_key)


def load_multiple_folders(path):
    """
    Loads subfolders from a specified directory into a dictionary.

    Parameters:
    path (str): Directory path to load subfolders from.

    Returns:
    dict: Subfolder names as keys and their full paths as values.

    Raises:
    SystemExit: If the path is invalid or the directory is empty.
    """

    if not os.path.exists(path) or not os.path.isdir(path):
        sys.exit(f"{path} is invalid path!")

    subfolders = [f.name for f in os.scandir(path) if f.is_dir()]

    if not subfolders:
        sys.exit("Directory is empty!")

    files = {}

    for folder in subfolders:
        files[folder] = os.path.join(path, folder)

    return files


def load_files_from_folder(path, file_format=".csv", n_sort=False):
    """
    Loads files of a specific format from a folder into a dictionary.

    Parameters:
    path (str): Directory path to load files from.
    file_format (str, optional): File format to filter by (default is ".csv").
    n_sort (bool, optional): Apply natural sorting if True (default is False).

    Returns:
    dict: Filenames as keys and their full paths as values.

    Raises:
    SystemExit: If the directory is empty.
    """

    if not os.listdir(path):
        sys.exit("Directory is empty")

    files_dict = {}

    for r, d, f in os.walk(path):
        if n_sort:
            f = natural_sort(f)
        for file in f:
            if file_format in file:
                files_dict.update({file: os.path.join(r, file)})

    return files_dict
