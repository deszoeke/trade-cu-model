import xarray as xr

#import numpy as np
#months_standard = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
#months_leap     = np.array([31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

# Load accumulated file
ds = xr.open_dataset("./data/isccp_weather_states_monthly_climatology_counts_1983_2017.nc")

# Calculate Relative Frequency of Occurrence (RFO) 
# Resulting shape will be (month: 12, centroid: 11, longitude: 360, latitude: 180)
rfo_climatology = ds["sum_of_timesteps"] / ds["total_valid_timesteps"]

# Quick save of the finished RFO profiles
rfo_climatology.to_netcdf("isccp_weather_states_monthly_RFO_climatology.nc")

