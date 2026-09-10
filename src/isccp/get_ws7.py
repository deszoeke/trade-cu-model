import os
import xarray as xr
import numpy as np

data_dir = "./data/global/"
def xtr7(file_path):
    with xr.open_dataset(file_path) as ds:
        # 1. Extract out the 42-element histogram vector for Weather State 7 (index 6)
        # The variable 'ws' contains the cloud fractions (in %) for each cluster centroid
        ws7_flat_histogram = ds['ws'].isel(n_centroids=6).values

        # 2. Reshape the flat 42-bin vector back into the physical 2D histogram matrix
        # ISCCP joint histograms are stored as (6 Optical Thickness ranges, 7 Pressure ranges)
        # and transposed to (7 Pressure ranges, 6 Optical Thickness ranges), per H1.py
        ws7_2d_matrix = ws7_flat_histogram.reshape(6, 7).T

    print("Extracted 2D Cloud Fraction Matrix for Weather State 7:")
    print(ws7_2d_matrix)
    return ws7_2d_matrix

for YEAR in range(1983, 1984):
    file = os.path.join(data_dir, f"{YEAR}.nc")
    xtr7(file)
