import os
import xarray as xr
import numpy as np

# calculate SW and LW kernels by multiplying the pc-tau 
# histogram by the zelinka pc-tau histogram kernel for each latitude and month.

# kernelfile = "/home/deszoeks/projects/ATOMIC/cloud-radiative-kernels/data/obs_cloud_kernels4.nc"
kernelfile = "../../../cloud-radiative-kernels/data/obs_cloud_kernels4.nc"
kernel = xr.open_dataset(kernelfile)
# kernel has dims (time=12, tau=7, plev=7, lat=72[, albcs=3 for SW only]);
# variables are 'SWkernel' and 'LWkernel' (not 'sw_kernel'/'lw_kernel')


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

file = os.path.join(data_dir, "1983.nc") # same every year
ws7_tau_pc_hist = xtr7(file) # (7 pc bins, 6 tau bins)

# Sanity check (per the original comment): the ws7 pc axis turns out to run
# TOA -> surface, i.e. REVERSED relative to the kernel's plev axis (which runs
# surface, 1000-800 hPa, -> TOA). Confirmed here: fraction is largest in the
# LAST two rows, which is physically correct for shallow cumulus (low-topped).
print("ws7 pc-row sums before reversal (expect largest LAST):", ws7_tau_pc_hist.sum(axis=1))
ws7_pc_tau = ws7_tau_pc_hist[::-1, :]  # now index 0 = near-surface (1000-800 hPa), matches kernel plev

# The kernel splits out a sub-visible tau bin (0.01-0.3) that ws7's 6-bin tau
# axis doesn't have, so ws7 has one fewer tau bin than the kernel (6 vs 7).
# Per plan: reuse ws7's thinnest tau bin (0.3-1.3) for both of the kernel's
# first two tau bins (0.01-0.3 and 0.3-1.3).
ws7_pc_tau7 = np.concatenate([ws7_pc_tau[:, :1], ws7_pc_tau], axis=1)  # (7 pc, 7 tau)

# Kernel arrays are (tau, plev[, albcs]); transpose ws7 to (tau, pc) to match.
ws7_hist_fraction = ws7_pc_tau7.T  # (7 tau, 7 pc), a 0-1 fraction (confirmed:
# ds['ws'] sums to ~0.05-0.99 per centroid, not ~5-99, so it is NOT already in %)

# The kernel is in W/m^2 per PERCENT cloud fraction (ds['SWkernel'].attrs['units']
# == 'W/m^2/%'), so the histogram must be in percentage points, not a 0-1 fraction.
ws7_hist = ws7_hist_fraction * 100.0  # (7 tau, 7 pc), percentage points

# for each calendar month and latitude,
# multiply the pc-tau histogram by the kernel for that month and latitude
months = np.arange(1, 13)
lats = kernel.lat.values  # kernel's native 2.5 deg grid (-88.75..88.75)
sw_by_month_lat = np.full((len(months), len(lats)), np.nan)
lw_by_month_lat = np.full((len(months), len(lats)), np.nan)

# interpolate the kernel to clear sky albedo of 0.086
alb = 0.086
wgts = np.array([alb, 1.0-alb, 0])/0.5 # alb <= 0.5

for mi, month in enumerate(months):
    for li, lat in enumerate(lats):
        kernel_sw = kernel['SWkernel'].sel(time=month, lat=lat).values  # (tau, plev, albcs)
        kernel_lw = kernel['LWkernel'].sel(time=month, lat=lat).values  # (tau, plev)

        # Surface albedo (albcs=0, 0.5, 1.0): use albcs=0 (dark/clear ocean).
        # kernel_sw_alb0 = kernel_sw[:, :, 0]  # (tau, plev)
        # inner product of weights and kernel along the albcs axis (axis=-1) interpolates
        kernel_sw_alb = np.tensordot(wgts, kernel_sw, axes=([-1], [-1]))  # (tau, plev)
        # kernel_sw_alb = interp_along_axis([0.086], [0, 0.5, 1.0], kernel_sw, axis=-1)[..., 0]

        # multiply the pc-tau histogram by the kernel
        sw_kernel_weighted = ws7_hist * kernel_sw_alb
        lw_kernel_weighted = ws7_hist * kernel_lw

        # sum over the pc-tau bins to get the total SW and LW kernels for this month and latitude
        swcre = sw_kernel_weighted.sum()
        lwcre = lw_kernel_weighted.sum()

        sw_by_month_lat[mi, li] = swcre
        lw_by_month_lat[mi, li] = lwcre
        print(f"Month: {month}, Latitude: {lat}, SW Kernel: {swcre}, LW Kernel: {lwcre}")

# This is the WS7 CRE per unit occurrence -- the histogram is
# constant, only the kernel varies by month/lat. Next step (separate/later):
# multiply and sum these by the actual WS7 RFO(lat, lon, month) -- e.g. from
# rfo_ws7_marine_tropics.py's counts file -- to get the real radiative effect
# at each (lat, lon, month).
out = xr.Dataset(
    {
        "SW_kernel": (("month", "lat"), sw_by_month_lat),
        "LW_kernel": (("month", "lat"), lw_by_month_lat),
    },
    coords={"month": months, "lat": lats},
)
out["SW_kernel"].attrs["long_name"] = "SW cloud radiative kernel weighted by the WS7 pc-tau histogram (albcs=0)"
out["LW_kernel"].attrs["long_name"] = "LW cloud radiative kernel weighted by the WS7 pc-tau histogram"
out["SW_kernel"].attrs["units"] = "W m-2 per unit WS7 occurrence (histogram already in cloud-fraction percent)"
out["LW_kernel"].attrs["units"] = "W m-2 per unit WS7 occurrence (histogram already in cloud-fraction percent)"

out_file = "ws7_cre_by_month_lat.nc"
out.to_netcdf(out_file)
print(f"Saved SW/LW WS7 kernel (month, lat) to {out_file}")
