"""
Relative frequency of occurrence (RFO) of ISCCP Weather State 7 (shallow
cumulus) over marine (ocean) pixels equatorward of 40 deg latitude, expressed
as a fraction of the WHOLE Earth's area: pixels poleward of 40 deg (and land
pixels) contribute 0 to the numerator, but the full global area is still the
denominator.

Uses the pre-aggregated per-gridpoint counts in
isccp_weather_states_monthly_climatology_counts_1983_2017.nc
(built by sum_centroid_counts.py) and the 1x1 deg land/ocean mask in
landmask/land_mask_1x1.nc.

Outputs:
  - the total-record RFO (all months and years pooled)
  - the RFO for each calendar month
  - a proportional-area ("egg", Mollweide) map of the record (annual) average
    climatology, wrapped at 30 deg E longitude
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

COUNTS_FILE = "./data/isccp_weather_states_monthly_climatology_counts_1983_2017.nc"
LAND_MASK_FILE = "./landmask/land_mask_1x1.nc"
CRE_FILE = "./ws7_cre_by_month_lat.nc"  # from ws7_kernel.py: SW/LW kernel x WS7 histogram, by (month, lat)
WS7_CENTROID = 7          # shallow cumulus weather state
LAT_LIMIT = 40.0          # deg, equatorward of this latitude counts toward RFO
WRAP_LON = 30.0           # deg E, map seam / wrap longitude for the egg plot
MAP_FILE = "ws7_rfo_annual_climatology_egg.png"
SWCRE_MAP_FILE = "ws7_swcre_annual_climatology_egg.png"

MONTH_NAMES = [
    "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December",
]


def load_ocean_mask(lon, lat):
    """Boolean ocean mask (True = ocean) on the counts file's (longitude, latitude) grid."""
    land = xr.open_dataset(LAND_MASK_FILE)["mask"]
    # land_mask lon is -179.5..179.5; counts lon is 0.5..359.5 -- put both on 0-360 and sort
    land = land.assign_coords(lon=((land.lon + 360) % 360)).sortby("lon")
    land = land.rename({"lon": "longitude", "lat": "latitude"})
    land = land.reindex(longitude=lon, latitude=lat, method="nearest")
    return land == 0  # 0 = ocean, 1 = land


def load_cre_kernel(lat):
    """
    SW/LW WS7 cloud radiative kernel (ws7_cre_by_month_lat.nc), interpolated
    from its native 2.5 deg latitude grid onto `lat` (the counts file's 1 deg
    grid). Returned DataArrays are (month, latitude); they broadcast over
    longitude automatically since they carry no longitude dimension.
    """
    cre = xr.open_dataset(CRE_FILE).rename({"lat": "latitude"})
    cre = cre.interp(latitude=lat, kwargs={"fill_value": "extrapolate"})
    return cre["SW_kernel"], cre["LW_kernel"]


def global_area_weighted_rfo(pixel_rfo, valid, area_weight_global, dims=("longitude", "latitude")):
    """
    Area-weighted RFO using the WHOLE global area (all lat/lon) as the
    denominator. Pixels where `valid` is False (land, or poleward of
    LAT_LIMIT) contribute 0 to the numerator, as if their RFO were 0.
    """
    numerator = (pixel_rfo.where(valid, 0.0) * area_weight_global).sum(dim=list(dims))
    denominator = area_weight_global.sum(dim=list(dims))
    return numerator / denominator


def plot_egg_map(field, lon, lat, out_file, vmin, vmax, title, cmap="Blues_r"):
    """
    Proportional-area (Mollweide, "egg") map of a record (annual) average
    field over ALL ocean pixels (full globe, not restricted to |lat| <=
    LAT_LIMIT), wrapped at WRAP_LON so the map seam falls at that meridian
    instead of the antimeridian. Land pixels (NaN in `field`) are filled
    black.
    """
    # Mollweide wraps at central_longitude +/- 180, so put the seam at
    # WRAP_LON by centering the map on WRAP_LON - 180.
    central_longitude = WRAP_LON - 180.0
    proj = ccrs.Mollweide(central_longitude=central_longitude)

    fig = plt.figure(figsize=(11, 6))
    ax = plt.axes(projection=proj)
    ax.set_global()

    cmap = plt.get_cmap(cmap, lut=20).copy()

    data = field.transpose("latitude", "longitude")
    mesh = ax.pcolormesh(
        lon, lat, data,
        transform=ccrs.PlateCarree(),
        cmap=cmap, vmin=vmin, vmax=vmax, shading="auto",
        zorder=0,
    )

    # Draw land ON TOP of the data (not via the colormap's "bad" color, and
    # not underneath): reprojecting the pcolormesh mesh from PlateCarree to
    # Mollweide smears colors across masked (NaN, land) cells instead of
    # leaving them transparent, so an opaque feature drawn underneath would
    # still get painted over. Covering land afterward hides that smearing.
    ax.add_feature(cfeature.LAND, facecolor="black", edgecolor="k", linewidth=0.5, zorder=1)

    # Reference lines at the +/- LAT_LIMIT latitudes used for the RFO/CRE stats.
    line_lons = np.linspace(-180, 180, 361)
    for lat_val in (-LAT_LIMIT, LAT_LIMIT):
        ax.plot(
            line_lons, np.full_like(line_lons, lat_val),
            transform=ccrs.PlateCarree(),
            color="k", linestyle="--", linewidth=1, zorder=2,
        )

    cb = plt.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.06, shrink=0.7)
    ax.set_title(title)

    fig.savefig(out_file, dpi=150, bbox_inches="tight")
    print(f"Saved annual climatology map to {out_file}")
    return fig, ax


def main():
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    ds = xr.open_dataset(COUNTS_FILE)
    counts7 = ds["sum_of_timesteps"].sel(centroid=WS7_CENTROID)  # (month, longitude, latitude)
    total_valid = ds["total_valid_timesteps"]                    # (month, longitude, latitude)

    ocean_mask = load_ocean_mask(ds.longitude, ds.latitude)      # (longitude, latitude)
    tropics_mask = np.abs(ds.latitude) <= LAT_LIMIT              # (latitude,)
    marine_tropics_mask = ocean_mask & tropics_mask              # (longitude, latitude)

    # Area weight: grid-cell area on a regular lat/lon grid is proportional to cos(lat).
    # This weight is computed over the WHOLE globe (used as the RFO denominator),
    # independent of the ocean/tropics mask.
    area_weight_global = np.cos(np.deg2rad(ds.latitude))
    area_weight_global = area_weight_global.broadcast_like(counts7.isel(month=0))

    # --- Per-calendar-month RFO ---
    valid_month = marine_tropics_mask & (total_valid > 0)
    pixel_rfo_month = counts7 / total_valid.where(total_valid > 0)
    monthly_rfo = global_area_weighted_rfo(pixel_rfo_month, valid_month, area_weight_global)

    # --- Total-record (annual) RFO: pool numerator and denominator across all months first ---
    counts7_total = counts7.sum(dim="month")
    total_valid_total = total_valid.sum(dim="month")
    valid_total = marine_tropics_mask & (total_valid_total > 0)
    pixel_rfo_total = counts7_total / total_valid_total.where(total_valid_total > 0)
    overall_rfo = global_area_weighted_rfo(pixel_rfo_total, valid_total, area_weight_global)

    print("RFO of WS7 (shallow cumulus), ocean & |lat| <= 40 deg, global-Earth-area denominator")
    print("=" * 85)
    for m, name in enumerate(MONTH_NAMES, start=1):
        val = float(monthly_rfo.sel(month=m))
        print(f"  {name:<10s}: {val:.4f}  ({val * 100:.2f}%)")
    print("-" * 85)
    print(f"  {'Total record':<10s}: {float(overall_rfo):.4f}  ({float(overall_rfo) * 100:.2f}%)")

    # --- Planetary-mean CRE from WS7 clouds, same domain/weighting as the RFO above ---
    # counts7 * kernel(month, lat) is the WS7-occurrence-weighted CRE contribution
    # per timestep; pooling this (instead of counts7 alone) across months before
    # dividing by total_valid_total reuses the RFO's own pooling weights, so months
    # with more valid observations are weighted accordingly, exactly as for RFO.
    sw_kernel, lw_kernel = load_cre_kernel(ds.latitude)
    counts7_sw = counts7 * sw_kernel      # (month, longitude, latitude)
    counts7_lw = counts7 * lw_kernel

    pixel_swcre_month = counts7_sw / total_valid.where(total_valid > 0)
    pixel_lwcre_month = counts7_lw / total_valid.where(total_valid > 0)
    monthly_swcre = global_area_weighted_rfo(pixel_swcre_month, valid_month, area_weight_global)
    monthly_lwcre = global_area_weighted_rfo(pixel_lwcre_month, valid_month, area_weight_global)

    counts7_sw_total = counts7_sw.sum(dim="month")
    counts7_lw_total = counts7_lw.sum(dim="month")
    pixel_swcre_total = counts7_sw_total / total_valid_total.where(total_valid_total > 0)
    pixel_lwcre_total = counts7_lw_total / total_valid_total.where(total_valid_total > 0)
    overall_swcre = global_area_weighted_rfo(pixel_swcre_total, valid_total, area_weight_global)
    overall_lwcre = global_area_weighted_rfo(pixel_lwcre_total, valid_total, area_weight_global)

    print()
    print("Planetary-mean CRE from WS7 clouds, ocean & |lat| <= 40 deg, global-Earth-area denominator")
    print("=" * 85)
    print(f"  {'Month':<10s}  {'SW (W/m2)':>12s}  {'LW (W/m2)':>12s}  {'Net (W/m2)':>12s}")
    for m, name in enumerate(MONTH_NAMES, start=1):
        sw = float(monthly_swcre.sel(month=m))
        lw = float(monthly_lwcre.sel(month=m))
        print(f"  {name:<10s}  {sw:12.4f}  {lw:12.4f}  {sw + lw:12.4f}")
    print("-" * 85)
    sw_tot, lw_tot = float(overall_swcre), float(overall_lwcre)
    print(f"  {'Total record':<10s}  {sw_tot:12.4f}  {lw_tot:12.4f}  {sw_tot + lw_tot:12.4f}")

    # Record (annual) average climatology maps: show all ocean pixels globally
    # (not restricted to |lat| <= LAT_LIMIT), with land set to NaN so it plots
    # black. The printed RFO/CRE statistics above remain equatorward of 40 deg.
    valid_ocean = ocean_mask & (total_valid_total > 0)

    annual_rfo_field = pixel_rfo_total.where(valid_ocean)
    plot_egg_map(
        annual_rfo_field, ds.longitude, ds.latitude, out_file=MAP_FILE,
        vmin=0, vmax=1, cmap="Blues_r",
        title="ISCCP-H shallow cumulus regime (WS7) frequency of occurrence",
    )

    annual_swcre_field = pixel_swcre_total.where(valid_ocean)
    plot_egg_map(
        annual_swcre_field, ds.longitude, ds.latitude, out_file=SWCRE_MAP_FILE,
        vmin=-40, vmax=0, cmap="Blues_r",
        title="Record average SW cloud radiative effect from WS7 (shallow cumulus), W m$^{-2}$",
    )

    return monthly_rfo, overall_rfo, monthly_swcre, monthly_lwcre, overall_swcre, overall_lwcre


if __name__ == "__main__":
    main()
