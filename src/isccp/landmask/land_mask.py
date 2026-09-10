import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import regionmask
import xarray as xr

# 1. Create a 1x1 degree global coordinate grid
# Latitudes from -89.5 to 89.5 (centers of the 1-degree bins)
# Longitudes from -179.5 to 179.5
lats = np.arange(-89.5, 90, 1.0)
lons = np.arange(-179.5, 180, 1.0)

# 2. Fetch the land mask from regionmask (Natural Earth land dataset)
# scale='110m' is perfect for a 1-degree grid. Use '50m' for higher precision.
land_mask_dataset = regionmask.defined_regions.natural_earth_v5_0_0.land_110

# 3. Create the 2D boolean mask over our grid coordinates
# Returns: 0 for land, NaN for ocean (or vice versa depending on your needs)
mask = land_mask_dataset.mask(lons, lats)

# Convert NaN to 0 (Ocean) and 0 to 1 (Land) for a traditional binary mask
binary_mask = xr.where(mask.isnull(), 0, 1)

# 4. Save to a standard NetCDF file
binary_mask.to_netcdf("land_mask_1x1.nc")
print("Successfully generated and saved 'land_mask_1x1.nc'!")

# 5. Visual Check: Plot the 1x1 grid mask
fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot the mask data
binary_mask.plot(
    ax=ax,
    transform=ccrs.PlateCarree(),
    cmap="Blues",
    cbar_kwargs={"label": "0 = Ocean, 1 = Land", "orientation": "horizontal"},
)

# Add coastlines to check alignment
ax.coastlines(color="red", linewidth=1)
ax.set_title("1x1 Degree Global Land Mask")

plt.show()

