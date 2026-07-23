using Pkg; Pkg.activate(".")
using Revise
using NCDatasets
using Statistics
# using Interpolations
using PythonPlot
using Printf

# Apply your global Arial font rules cleanly across the global scope
rc("font", family="sans-serif")
rc("font"; Symbol("sans-serif") => "Helvetica")

datadir = "../../../ATOMIC_GOES/data/"

# ==============================================================================
# 1. RADIATIVE TRANSFER LOOKUPS
# ==============================================================================

# not used. Every cloudy pixel counts as 100 % cloudy in the ISCCP binning
# no subpixel mixed cloud fraction
# function compute_subpixel_cloud_fraction(alpha_obs, sza_deg, tau, cth_km, alpha_surface)
#     if isnan(tau) || isnan(alpha_obs) || isnan(sza_deg) || isnan(cth_km) || isnan(alpha_surface)
#         return NaN
#     end
#     if tau <= 0.0 || alpha_obs <= alpha_surface
#         return 0.0
#     end
#     mu_0 = cos(deg2rad(sza_deg))
#     if mu_0 <= 0.0 return 0.0 end

#     g = cth_km > 6.0 ? 0.75 : 0.85
#     # gamma_1 = 0.75 * (1.0 - (g * mu_0))
#     # Coakley-Chylek 2-stream approximation for cloud albedo:
#     u = 0.75 * (1.0 - g)
#     alpha_cloud = alpha_cloud = (u * tau) / (u * tau + mu_0)
#     # ISCCP method for multiple scattering
#     b = 6.8 * mu_0          # approximate ISCCP coefficient
#     tp = sqrt(sqrt(tau^3)) # tau^0.75
#     alpha_cloud = tp / (tp + b) # power is 0.75 for liquid clouds?
#     # ISCCP from Zelinka
#     # tp = τ^0.895 # why is it different from 0.75 above?
#     # αc = tp∕(tp + 6.82),
#     if alpha_cloud <= alpha_surface return 0.0 end
#     return clamp((alpha_obs - alpha_surface) / (alpha_cloud - alpha_surface), 0.0, 1.0)
# end

# radiative functions for cloud optical thickness and albedo

"""
asymmetry parameter for the scattering phase function of cloud droplets, 
as a function of effective radius (microns)
"""
asym_parameter(r_e) = 0.823 + 0.0035 * r_e
asym_parameter_kernel = asym_parameter(10.0) # 0.858 for 10 micron effective radius
tau_fac(r_e) = (1-asym_parameter(r_e)) / (1-asym_parameter_kernel)
calc_tau_scaled(r_e, tau) = tau * tau_fac(r_e)
"albedo asymptotically matching radiative calculations of kernel"
function albedo_kernel(tau_scaled, sza=90.0, r_e=10.0)
    omgts = (1-asym_parameter(r_e)) * tau_scaled
    omgts / (2*cosd(sza) + omgts)
end


function get_bin_index(value, edges)
    if isnan(value) || value < edges[1] || value > edges[end] return 0 end
    return searchsortedlast(edges, value)
end

# ==============================================================================
# 2. FUNCTIONS TO ACCUMULATE PIXELS IN-PLACE 
# ==============================================================================

"""
    update_albedo_profile!(acc_cloudy, acc_clear_profile, acc_clear_scalar, global_denominator, nc_file, ...)
Processes a single granule file. Allocates 100% clear space to a single scalar accumulator,
while splitting subpixel cloud structures into the 2D joint matrix and 1D pressure profile.
All trackers are modified strictly in place.
"""
function update_albedo_profile!(
    albedo_profile_accumulator::Vector{Float64},
    reflec_profile_accumulator::Vector{Float64},
    cloudy_profile_count::Vector{Int64},
    cloudy_all_count::Vector{Int64},
    clear_all_count::Vector{Int64},
    total_all_count::Vector{Int64}, # <- updates
    nc_file::String,
    lat_bounds::Tuple{Float64, Float64},
    lon_bounds::Tuple{Float64, Float64},
    vza_thr = 60.0, 
    sza_thr = 60.0 )

    if !isfile(nc_file)
        println("Warning: file not found: ", nc_file)
        return
    end

    lat_min, lat_max = lat_bounds
    lon_min, lon_max = lon_bounds

    ds = NCDatasets.Dataset(nc_file, "r")
    
    lats            = coalesce.(ds["latitude"][:,:], NaN32)
    lons            = coalesce.(ds["longitude"][:,:], NaN32)
    reflectance_vis = coalesce.(ds["reflectance_vis"][:,:], NaN32)
    tau             = coalesce.(ds["cloud_visible_optical_depth"][:,:], NaN32)
    particle_size   = coalesce.(ds["cloud_particle_size"][:,:], NaN32)
    pixel_sza       = coalesce.(ds["pixel_sza"][:,:], NaN32)
    pixel_vza       = coalesce.(ds["pixel_vza"][:,:], NaN32)
    # pc              = coalesce.(ds["cloud_top_pressure"][:,:], NaN32)
    cloud_height    = coalesce.(ds["cloud_top_height"][:,:], NaN32)

    close(ds)

    # radiative calculations for each pixel
    tau_scaled = calc_tau_scaled.(particle_size, tau)
    albedo = albedo_kernel.(tau_scaled, pixel_sza, particle_size)

    # ISCCP-consistent binary classification: cloudy if retrieval produced a valid tau
    cloud_mask = .!isnan.(tau) .& (tau .> 0.0)

    # Universal geographic / solar zenith condition mask
    spatial_mask = (lats .>= lat_min) .& (lats .<= lat_max) .& 
                   (lons .>= lon_min) .& (lons .<= lon_max) .&
                   (pixel_vza .<= vza_thr) .& (pixel_sza .<= sza_thr) .& 
                   (.!isnan.(reflectance_vis))

    # valid subsets
    v_cloud = cloud_mask[spatial_mask]
    v_height = cloud_height[spatial_mask]
    v_tau   = tau_scaled[spatial_mask]
    v_albedo = albedo[spatial_mask]
    v_reflectance = reflectance_vis[spatial_mask]

    file_valid_footprints = length(v_cloud)
    if file_valid_footprints == 0 return end

    # Accumulate footprint contents down the track
    for i in 1:file_valid_footprints
        total_all_count[1] += 1
        if !v_cloud[i] # clear
            clear_all_count[1] += 1
        else           # cloud
            cloudy_all_count[1] += 1
            h_bin = get_bin_index(1e3*v_height[i], height_bins)
            if 1 <= h_bin <= nh
                cloudy_profile_count[h_bin] += 1
                albedo_profile_accumulator[h_bin] += v_albedo[i]
                reflec_profile_accumulator[h_bin] += v_reflectance[i]
            end
        end
    end

    println("$nc_file ($file_valid_footprints pixels)")
end

"""
count cloudy and clear pixels as a function of height.
count pixels, reflectance weighted pixels, and albedo weighted pixels.
"""
function compile_albedo_profile(lat_bounds, lon_bounds, data_file_list)
    height_bins = 500.0:10.0:4000.0
    nh = length(height_bins)

    # --- PREALLOCATE ACCUMULATORS ---
    albedo_profile_accumulator = zeros(Float64, nh)
    reflec_profile_accumulator = zeros(Float64, nh)
    cloudy_profile_count = zeros(Int64, nh)
    cloudy_all_count = zeros(Int64, 1)
    clear_all_count  = zeros(Int64, 1)
    total_all_count  = zeros(Int64, 1) # denominator of all scanned pixels

    # Run state ingestion loop sequence
    for file in data_file_list
        update_albedo_profile!(
            albedo_profile_accumulator,
            reflec_profile_accumulator,
            cloudy_profile_count,
            cloudy_all_count,
            clear_all_count,
            total_all_count,
            joinpath(datadir, "GOES/all", file), 
            lat_bounds, 
            lon_bounds )
    end

    return (albedo_profile_accumulator, 
        reflec_profile_accumulator, 
        cloudy_profile_count, 
        cloudy_all_count, 
        clear_all_count, 
        total_all_count)
end

# ==============================================================================
# 3. RUN PARAMETERS
# ==============================================================================

# EUREC4A region
lat_bounds = (12.5, 16.0)
lon_bounds = (-60.0, -49.0)
# files matching data
daylight_file(s) = 1200 <= parse(Int,match(r"(\d{4})\.PX\.02K\.NC$", s).captures[1]) <= 1920
data_file_list = filter(daylight_file, readdir(joinpath(datadir, "GOES/all"))) # [1:3]

# Compile albedo and reflectance-weighted clouds
(   albedo_profile_accumulator,
    reflec_profile_accumulator,
    cloudy_profile_count,
    cloudy_all_count,
    clear_all_count,
    total_all_count ) = compile_albedo_profile(
    lat_bounds, lon_bounds, data_file_list[:] )

# ==============================================================================
# 4. NORMALIZE FRACTIONS
# ==============================================================================
albedo_profile = albedo_profile_accumulator ./ total_all_count[1]
reflec_profile = reflec_profile_accumulator ./ total_all_count[1]
cloud_profile  = cloudy_profile_count       ./ total_all_count[1] # fraction of full pixels
cloud_total_frac = cloudy_all_count[1] / total_all_count[1]
clear_total_frac = clear_all_count[1]  / total_all_count[1]

begin
    println("Total pixel count: ", total_all_count[1])
    println("Regional pure clear pixel fraction: ", round(100*clear_total_frac, digits=2), " %")
    println("Regional cloudy pixel fraction: ", round(100*cloud_total_frac, digits=2), " %")
    println("Regional mean albedo: ", round(sum(albedo_profile), digits=3))
    println("Regional mean reflectance: ", round(sum(reflec_profile), digits=3))
end

#=
# Save the profiles in a netcdf file.
# first dump and edit the obs kernel cdl file, then
# ncgen -o shcu_isccp_cloud_pct.nc shcu_isccp_cloud_pct.cdl
# copy the variables we want
# ncks -A -C -v plev_bnds,tau_bnds,plev,tau obs_cloud_kernels4.nc shcu_isccp_cloud_pct.nc
NCDatasets.Dataset(joinpath(datadir, "shcu_cloud_albedo_refl_profile.nc"), "a") do ald
    ald["albedo_profile"][:] .= albedo_profile
    ald["reflec_profile"][:] .= reflec_profile
    ald["cloud_profile"][:]  .= cloud_profile
    ald["cloud_total_frac"][1] = cloud_total_frac
    ald["clear_total_frac"][1] = clear_total_frac
    ald["total_all_count"][1] = total_all_count[1]
end
=#

# plot the profiles
clf()
plot(albedo_profile, height_bins/1e3, label="albedo-fraction")
plot(reflec_profile, height_bins/1e3, label="reflectance-fraction")
plot(cloud_profile, height_bins/1e3, label="cloud fraction")
legend(frameon=false)
xlim([0, 0.003])
plot(1.4*albedo_profile, height_bins/1e3, label="albedo-fraction")
