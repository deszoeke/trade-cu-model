import os
import glob
import numpy as np
import xarray as xr

# Paths
DATA_DIR = "./data/global"  # Adjust to your download folder path
OUTPUT_FILE = "isccp_weather_states_monthly_climatology_counts_1983_2017.nc"

# Map out monthly variable names expected across files
MONTHS = [
    "january", "february", "march", "april", "may", "june",
    "july", "august", "september", "october", "november", "december"
]

# Grid parameters from metadata
NLON, NLAT = 360, 180
NCENTROIDS = 11  # Weather States 1 through 11; WS7 is shallow cumulus
NMONTHS = 12     # 12 Calendar Months

def main():
    # Find all annual netcdf files
    file_paths = sorted(glob.glob(os.path.join(DATA_DIR, "*.nc")))
    if not file_paths:
        print(f"No NetCDF files found in {DATA_DIR}. Please check your path.")
        return

    # 1. Initialize arrays to accumulate counts across separate months
    # Shape: (12 months, 11 centroids, 360 longitude, 180 latitude)
    total_timesteps_accum = np.zeros((NMONTHS, NCENTROIDS, NLON, NLAT), dtype=np.int64)
    # Shape: (12 months, 360 longitude, 180 latitude) - tracking valid steps per grid box per month
    total_valid_timesteps = np.zeros((NMONTHS, NLON, NLAT), dtype=np.int64)

    # Coordinates placeholders to build output later
    lon_coords = None
    lat_coords = None

    print(f"Found {len(file_paths)} files. Starting processing loop...")

    for path in file_paths:
        print(f"Processing file: {os.path.basename(path)}")
        
        # Open individual file using xarray
        with xr.open_dataset(path) as ds:
            # Capture spatial coordinates from the first valid file
            if lon_coords is None and 'longitude' in ds.coords:
                lon_coords = ds.longitude.values
                lat_coords = ds.latitude.values

            # Loop through each month index (0 to 11) and variable name
            for m_idx, month_var in enumerate(MONTHS):
                if month_var in ds.variables:
                    # Load 3D matrix into memory (Time, Longitude, Latitude)
                    month_data = ds[month_var].values
                    
                    # Track valid, non-missing timesteps for the denominator for this month
                    valid_mask = (month_data >= 1) & (month_data <= NCENTROIDS)
                    total_valid_timesteps[m_idx, :, :] += np.sum(valid_mask, axis=0)
                    
                    # Compute counts for each centroid (1 to 11) for this month
                    for c_idx in range(1, NCENTROIDS + 1):
                        centroid_counts_this_month = np.sum(month_data == c_idx, axis=0)
                        total_timesteps_accum[m_idx, c_idx - 1, :, :] += centroid_counts_this_month

    # Fallback coordinates if the NetCDF metadata lacks explicit lon/lat coordinate arrays
    if lon_coords is None:
        lon_coords = np.arange(0.5, 360.5, 1.0)  # Standard 1-deg centers
        lat_coords = np.arange(-89.5, 90.5, 1.0)

    # 2. Package accumulation into a clean xarray Dataset
    output_ds = xr.Dataset(
        data_vars={
            "sum_of_timesteps": (["month", "centroid", "longitude", "latitude"], total_timesteps_accum),
            "total_valid_timesteps": (["month", "longitude", "latitude"], total_valid_timesteps)
        },
        coords={
            "month": np.arange(1, NMONTHS + 1),       # 1 = January, 12 = December
            "centroid": np.arange(1, NCENTROIDS + 1), # 1 to 11
            "longitude": lon_coords,
            "latitude": lat_coords
        },
        attrs={
            "description": "Total sum of timesteps per centroid and valid timesteps clustered by calendar month (1983-2017).",
            "history": "Aggregated from raw 3-hourly observations. Separated into 12 month buckets across all target years."
        }
    )

    # 3. Save compilation to disk
    print(f"Saving monthly climatology dataset to {OUTPUT_FILE}...")
    output_ds.to_netcdf(OUTPUT_FILE)
    print("Done!")

if __name__ == "__main__":
    main()

