#!/usr/bin/env python3
"""
Compute the cos(lat)-weighted mean clear-sky shortwave (TOA) albedo over
ocean grid points between 40S and 40N, from CERES EBAF-TOA climatology.

Albedo is computed as a FLUX RATIO (not a mean-of-ratios):
    alpha = sum(w * F_clr_up) / sum(w * F_down)
where w = cos(lat) * ocean_mask, summed over all qualifying lat/lon points.
This is the physically meaningful area-weighted albedo (energy-conserving),
as opposed to averaging the per-gridpoint ratios, which would overweight
low-insolation points.

Assumes:
    ./ceres-ebaf/CERES_EBAF-TOA_Edition4.2_200003-202306.nc
    ./landmask/land_mask_1x1.nc

Adjust CERES_FILE / MASK_FILE below if paths differ.
"""

import numpy as np
import xarray as xr

CERES_FILE = "./ceres-ebaf/CERES_EBAF-TOA_Edition4.2_200003-202306.nc"
MASK_FILE = "./landmask/land_mask_1x1.nc"

LAT_MIN, LAT_MAX = -40.0, 40.0

# Which mask value means "ocean"? We auto-detect below by checking which
# choice gives a global ocean fraction close to the known ~0.71. Override
# here if you already know the convention (0 or 1) for your file.
OCEAN_MASK_VALUE = None  # None = auto-detect


def load_data():
    ceres = xr.open_dataset(CERES_FILE)
    mask_ds = xr.open_dataset(MASK_FILE)
    return ceres, mask_ds


def align_mask_to_ceres(ceres, mask_ds):
    """Return the land/ocean mask on the CERES lat/lon grid."""
    mask = mask_ds["mask"]

    same_grid = (
        mask.sizes.get("lat") == ceres.sizes.get("lat")
        and mask.sizes.get("lon") == ceres.sizes.get("lon")
        and np.allclose(mask["lat"].values, ceres["lat"].values, atol=1e-3)
        and np.allclose(
            np.mod(mask["lon"].values, 360),
            np.mod(ceres["lon"].values, 360),
            atol=1e-3,
        )
    )

    if same_grid:
        mask = mask.assign_coords(lat=ceres["lat"].values, lon=ceres["lon"].values)
        return mask

    # Otherwise, reindex/interpolate the mask onto the CERES grid (nearest
    # neighbor, since this is a categorical 0/1 field).
    mask = mask.reindex(
        lat=ceres["lat"].values, lon=ceres["lon"].values, method="nearest"
    )
    return mask


def pick_ocean_value(mask, lat):
    """Auto-detect whether ocean == 0 or ocean == 1 using the known global
    ocean fraction (~0.71), area-weighted by cos(lat)."""
    if OCEAN_MASK_VALUE is not None:
        return OCEAN_MASK_VALUE

    w = np.cos(np.deg2rad(lat))
    w2d = xr.DataArray(w, coords={"lat": lat}, dims=["lat"]).broadcast_like(mask)

    frac_is_0 = float((w2d.where(mask == 0)).sum() / w2d.sum())
    frac_is_1 = float((w2d.where(mask == 1)).sum() / w2d.sum())

    target = 0.71
    choice = 0 if abs(frac_is_0 - target) < abs(frac_is_1 - target) else 1
    print(
        f"[mask autodetect] cos(lat)-weighted fraction where mask==0: {frac_is_0:.3f}\n"
        f"[mask autodetect] cos(lat)-weighted fraction where mask==1: {frac_is_1:.3f}\n"
        f"[mask autodetect] -> treating mask == {choice} as OCEAN "
        f"(expected global ocean fraction ~0.71)"
    )
    return choice


def clean(da):
    """Replace CERES fill values / out-of-range values with NaN."""
    fv = da.encoding.get("_FillValue", da.attrs.get("_FillValue", None))
    da = da.where(da != -999.0)
    vmin = da.attrs.get("valid_min", None)
    vmax = da.attrs.get("valid_max", None)
    if vmin is not None:
        da = da.where(da >= vmin)
    if vmax is not None:
        da = da.where(da <= vmax)
    return da


def main():
    ceres, mask_ds = load_data()

    sw_clr = clean(ceres["toa_sw_clr_c_clim"])  # (ctime, lat, lon)
    sw_down = clean(ceres["solar_clim"])         # (ctime, lat, lon)

    lat = ceres["lat"].values

    mask = align_mask_to_ceres(ceres, mask_ds)
    ocean_val = pick_ocean_value(mask, lat)
    ocean_mask = (mask == ocean_val)

    # Restrict to 40S-40N
    lat_band = (ceres["lat"] >= LAT_MIN) & (ceres["lat"] <= LAT_MAX)

    sw_clr = sw_clr.sel(lat=lat_band)
    sw_down = sw_down.sel(lat=lat_band)
    ocean_mask_band = ocean_mask.sel(lat=lat_band)

    lat_band_vals = sw_clr["lat"].values
    coslat = xr.DataArray(
        np.cos(np.deg2rad(lat_band_vals)), coords={"lat": lat_band_vals}, dims=["lat"]
    )

    # Combined weight: cos(lat) x ocean mask (broadcasts over ctime, lon)
    weight = coslat * ocean_mask_band.astype(float)

    # --- Per-month (climatological month) flux-weighted albedo ---
    num = (sw_clr * weight).sum(dim=["lat", "lon"], skipna=True)
    den = (sw_down * weight).sum(dim=["lat", "lon"], skipna=True)
    monthly_albedo = num / den

    print("\nMonthly climatological ocean clear-sky albedo, 40S-40N:")
    for i, val in enumerate(monthly_albedo.values, start=1):
        print(f"  Month {i:2d}: alpha_clr = {val:.4f}")

    # --- Annual mean: ratio of annual-mean fluxes (energy-consistent),
    #     not a mean of the monthly ratios ---
    annual_num = num.mean(dim="ctime")
    annual_den = den.mean(dim="ctime")
    annual_albedo_fluxratio = float(annual_num / annual_den)

    # For comparison: simple mean of the monthly ratios
    annual_albedo_meanofratios = float(monthly_albedo.mean(dim="ctime"))

    print(f"\nAnnual mean clear-sky ocean albedo (40S-40N):")
    print(f"  Ratio of annual-mean fluxes : {annual_albedo_fluxratio:.4f}")
    print(f"  Mean of monthly ratios      : {annual_albedo_meanofratios:.4f}")

    # Sanity-check: fraction of area in the 40S-40N band that is ocean
    coslat_full = coslat.broadcast_like(ocean_mask_band).astype(float)
    ocean_weight_sum = float((coslat_full * ocean_mask_band.astype(float)).sum())
    total_weight_sum = float(coslat_full.sum())
    band_ocean_frac = ocean_weight_sum / total_weight_sum
    print(
        f"\n[sanity check] cos(lat)-weighted ocean fraction within 40S-40N: "
        f"{band_ocean_frac:.3f} (typically ~0.75-0.80 for this latitude band)"
    )

    # ------------------------------------------------------------------
    # MEAN OF THE QUOTIENT (pointwise albedo at each lat/lon/month, then
    # averaged), using the full monthly time series rather than the
    # 12-month climatology. This differs from the flux-ratio (quotient of
    # the mean) computed above whenever clear-sky flux and insolation
    # covary across space/time (e.g. AOD or residual cloud contamination
    # correlated with insolation) - exactly the effect being tested here.
    # ------------------------------------------------------------------
    sw_clr_mon = clean(ceres["toa_sw_clr_c_mon"]).sel(lat=lat_band)   # (time, lat, lon)
    sw_down_mon = clean(ceres["solar_mon"]).sel(lat=lat_band)         # (time, lat, lon)

    pointwise_albedo = (sw_clr_mon / sw_down_mon.where(sw_down_mon > 0)).where(
        ocean_mask_band
    )

    weight_sum = float(weight.sum())  # same at every time step (mask is static)

    # cos(lat)*ocean-weighted spatial mean at each individual month
    albedo_by_time = (pointwise_albedo * weight).sum(dim=["lat", "lon"], skipna=True) / weight_sum

    mean_of_quotient = float(albedo_by_time.mean(dim="time", skipna=True))
    std_of_quotient = float(albedo_by_time.std(dim="time", skipna=True))

    print(f"\nMean of the quotient (pointwise ratio, full monthly record, 40S-40N ocean):")
    print(f"  mean(alpha_clr) over all months = {mean_of_quotient:.4f}  "
          f"(std across months = {std_of_quotient:.4f}, n_months = {albedo_by_time.sizes['time']})")
    print(f"\nFor comparison, quotient of the mean (climatology, computed above):")
    print(f"  alpha_clr = mean(F_clr) / mean(F_down) = {annual_albedo_fluxratio:.4f}")
    print(f"  difference (mean-of-quotient minus quotient-of-mean) = "
          f"{mean_of_quotient - annual_albedo_fluxratio:+.4f}")


if __name__ == "__main__":
    main()

