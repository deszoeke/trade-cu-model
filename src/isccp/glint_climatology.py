#!/usr/bin/env python3
"""
tropical_ocean_glint_climatology.py

Data pipeline / climatology layer: computes solar geometry and area+flux-
weights the Cox-Munk physics from coxmunk_physics.py over 40S-40N ocean.

Wind-speed dependence in the Cox-Munk glint+whitecap albedo is weak over
the range of trade-wind-relevant conditions (see coxmunk_physics.py
__main__ output: <0.01 albedo change across U=0-25 m/s except at sun
angles >60 deg, which contribute little to insolation-weighted means).
Given that, this script uses a single representative wind speed
(FIXED_WIND_SPEED, default 8 m/s - typical trade-wind value) rather than
loading a full wind climatology. All averaging is done over solar
geometry only, flux-weighted by cos(theta0) (the actual incident solar
radiation), area-weighted by cos(lat), and restricted to ocean.

If you want to reintroduce a real wind climatology later (e.g. to check
whether the "weak dependence" assumption breaks down somewhere, such as
strong dust-driven wind events or Southern Ocean fringe of the band),
swap FIXED_WIND_SPEED for a loaded U(lat,lon,month) field - the rest of
the pipeline (theta0 computation, flux-weighting, spatial aggregation)
is unchanged either way.

INPUT FILE (edit path as needed)
----------------------------------
    ./landmask/land_mask_1x1.nc   - land/ocean mask (0/1, convention
                                     auto-detected against the known ~0.71
                                     global ocean fraction). Also supplies
                                     the lat/lon grid, since no wind file
                                     is loaded.
"""

import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

import coxmunk_physics as cm

# ----------------------------------------------------------------------
# Configuration - EDIT THESE
# ----------------------------------------------------------------------
MASK_FILE = "./landmask/land_mask_1x1.nc"

LAT_MIN, LAT_MAX = -40.0, 40.0
N_HOUR_ANGLE_NODES = 8   # Gauss-Legendre nodes for diurnal flux-weighting
OCEAN_MASK_VALUE = None  # None = auto-detect
FIXED_WIND_SPEED = 8.0   # m/s - representative trade-wind value; Cox-Munk
                          # glint+whitecap albedo is only weakly sensitive
                          # to this (see coxmunk_physics.py __main__ output)


# ----------------------------------------------------------------------
# Solar geometry (pure functions, no I/O - could also live in a separate
# "solar_geometry.py" if it grows; kept here since it's climatology-glue,
# not surface-reflectance physics)
# ----------------------------------------------------------------------
MID_MONTH_DOY = np.array([17, 47, 75, 105, 136, 166, 197, 228, 258, 289, 319, 345])


def solar_declination(doy):
    return np.deg2rad(23.44) * np.sin(2 * np.pi * (284 + doy) / 365.0)


def hour_angle_nodes_weights(lat_rad, decl, n_nodes=N_HOUR_ANGLE_NODES):
    """Gauss-Legendre nodes/weights over daylight hour angle [-h0, h0]."""
    cos_h0 = np.clip(-np.tan(lat_rad) * np.tan(decl), -1.0, 1.0)
    h0 = np.arccos(cos_h0)

    x, w = np.polynomial.legendre.leggauss(n_nodes)
    h = 0.5 * (x + 1) * (2 * h0) - h0
    wh = w * h0
    return h, wh


# ----------------------------------------------------------------------
# File loading (auto-detect grid conventions)
# ----------------------------------------------------------------------
def load_ocean_mask_and_grid():
    """Load the mask file and return (ocean_mask, lat, lon) restricted to
    LAT_MIN..LAT_MAX. The mask file also supplies the working lat/lon grid
    since no wind file is loaded."""
    mask_ds = xr.open_dataset(MASK_FILE)
    mask = mask_ds["mask"]

    lat_name = "lat" if "lat" in mask.coords else "latitude"
    lon_name = "lon" if "lon" in mask.coords else "longitude"
    mask = mask.rename({lat_name: "lat", lon_name: "lon"})

    band = (mask["lat"].values >= LAT_MIN) & (mask["lat"].values <= LAT_MAX)
    mask = mask.isel(lat=np.where(band)[0])

    lat = mask["lat"].values
    lon = mask["lon"].values

    w = np.cos(np.deg2rad(lat))
    w2d = xr.DataArray(w, coords={"lat": lat}, dims=["lat"]).broadcast_like(mask)
    frac0 = float((w2d.where(mask == 0)).sum() / w2d.sum())
    frac1 = float((w2d.where(mask == 1)).sum() / w2d.sum())
    choice = OCEAN_MASK_VALUE
    if choice is None:
        # NOTE: this 0.71 target is the *global* ocean fraction; restricted
        # to 40S-40N the true fraction is somewhat higher (~0.75-0.80), so
        # this autodetect is only meant to pick the right mask convention
        # (0 vs 1), not to validate the band-restricted fraction itself.
        choice = 0 if abs(frac0 - 0.71) < abs(frac1 - 0.71) else 1
    print(f"[mask] treating mask == {choice} as OCEAN "
          f"(cos(lat)-weighted fractions within {LAT_MIN}-{LAT_MAX}: "
          f"mask==0 -> {frac0:.3f}, mask==1 -> {frac1:.3f})")

    return (mask == choice), lat, lon


# ----------------------------------------------------------------------
# Main pipeline
# ----------------------------------------------------------------------
def main():
    print("Building Cox-Munk glint lookup table (see coxmunk_physics.py)...")
    theta0_rad, U_grid, glint_table = cm.build_glint_lookup()
    if not cm.self_test(theta0_rad, U_grid, glint_table):
        raise SystemExit("coxmunk_physics self-test failed - aborting.")

    interp = RegularGridInterpolator(
        (theta0_rad, U_grid), glint_table, bounds_error=False, fill_value=None
    )
    glint_lookup_fn = lambda th0, U: interp(np.stack([th0.ravel(), U.ravel()], axis=-1)).reshape(th0.shape)

    ocean_mask, lat_band, lon = load_ocean_mask_and_grid()

    coslat = xr.DataArray(np.cos(np.deg2rad(lat_band)), coords={"lat": lat_band}, dims=["lat"])
    spatial_weight = (coslat * ocean_mask.astype(float)).values  # (lat, lon)

    n_lat, n_lon = len(lat_band), len(lon)
    lat_rad_col = np.deg2rad(lat_band)[:, None]

    U_month = np.full((n_lat, n_lon), FIXED_WIND_SPEED)  # constant, per module docstring

    monthly_flat = np.zeros(12)
    monthly_total = np.zeros(12)
    monthly_weight = np.zeros(12)

    for m in range(12):
        decl = solar_declination(MID_MONTH_DOY[m])
        h_nodes, h_weights = hour_angle_nodes_weights(lat_rad_col, decl)  # (n_lat, N_nodes)

        cos_theta0 = (
            np.sin(lat_rad_col) * np.sin(decl)
            + np.cos(lat_rad_col) * np.cos(decl) * np.cos(h_nodes)
        )
        cos_theta0 = np.clip(cos_theta0, 0.0, 1.0)  # daylight only
        theta0 = np.arccos(cos_theta0)  # (n_lat, N_nodes)

        num_total = np.zeros((n_lat, n_lon))
        num_flat = np.zeros((n_lat, n_lon))
        den = np.zeros((n_lat, n_lon))

        for k in range(h_nodes.shape[1]):
            th0_k = theta0[:, k]
            w_k = h_weights[:, k]
            cos_k = cos_theta0[:, k]

            th0_grid = np.repeat(th0_k[:, None], n_lon, axis=1)
            alpha_total_k = cm.total_albedo(th0_grid, U_month, glint_lookup_fn)
            alpha_flat_k = cm.fresnel_unpolarized(th0_grid)  # U=0, no whitecap baseline

            flux_w = (cos_k * w_k)[:, None]  # this IS "average incident solar
                                               # radiation cos(theta0)" - the
                                               # weight used throughout
            num_total += alpha_total_k * flux_w
            num_flat += alpha_flat_k * flux_w
            den += np.broadcast_to(flux_w, (n_lat, n_lon))

        albedo_total_month = num_total / den
        albedo_flat_month = num_flat / den
        insolation_proxy = den

        w_full = spatial_weight * insolation_proxy
        monthly_total[m] = np.nansum(albedo_total_month * w_full) / np.nansum(w_full)
        monthly_flat[m] = np.nansum(albedo_flat_month * w_full) / np.nansum(w_full)
        monthly_weight[m] = np.nansum(w_full)

    print(f"\n(fixed wind speed = {FIXED_WIND_SPEED} m/s throughout)")
    print("Monthly climatological ocean albedo, 40S-40N "
          "(flux-weighted diurnally by cos(theta0), area-weighted spatially):")
    print(f"{'Month':>6} {'Flat (U=0)':>12} {'Glint+foam':>12} {'Excess':>10}")
    for m in range(12):
        print(f"{m+1:6d} {monthly_flat[m]:12.4f} {monthly_total[m]:12.4f} "
              f"{monthly_total[m]-monthly_flat[m]:10.4f}")

    w = monthly_weight
    annual_flat = np.sum(monthly_flat * w) / np.sum(w)
    annual_total = np.sum(monthly_total * w) / np.sum(w)

    print(f"\nAnnual (insolation-weighted) mean, 40S-40N ocean:")
    print(f"  Flat-ocean baseline (U=0, no whitecap):        {annual_flat:.4f}")
    print(f"  Cox-Munk glint + whitecap (U={FIXED_WIND_SPEED:.0f} m/s):        {annual_total:.4f}")
    print(f"  Wind-driven excess:                             {annual_total - annual_flat:+.4f}")
    print(
        "\nCompare this excess to the EBAF-vs-Zheng gap discussed earlier: "
        "if the excess here is small relative to that gap, it supports "
        "undetected cloud contamination (not glint sampling) as the dominant "
        "explanation; if comparable, glint-driven sampling bias in satellite "
        "clear-sky screening deserves more weight."
    )


if __name__ == "__main__":
    main()

