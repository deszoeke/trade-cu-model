using Pkg; Pkg.activate(".")
using Revise
using NCDatasets
using Dates
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

function yyyydoy_HHMM2dt(s::AbstractString)
    # Slice out the numeric components based on string position
    yr   = parse(Int, s[1:4])       # "2020"
    doy  = parse(Int, s[5:7])       # "365"
    hour = parse(Int, s[9:10])      # "12"
    min  = parse(Int, s[11:12])     # "30"
    
    # Add Day of Year to the first day of that year (Subtract 1 since January 1st is Day 1)
    base_date = Date(yr) + Day(doy - 1)
    
    # Combine the Date and Time into a final DateTime object
    return DateTime(base_date, Time(hour, min))
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
    cloud_profile_count::Vector{Int64},
    cloud_all_count::Vector{Int64},
    clear_all_count::Vector{Int64},
    total_all_count::Vector{Int64}, # <- updates
    nc_file::String,
    lat_bounds::Tuple{Float64, Float64},
    lon_bounds::Tuple{Float64, Float64},
    vza_thr=75.0, 
    sza_thr=75.0,
    height_bins=0.0:10.0:4000.0 )

    nh = length(height_bins)

    if !isfile(nc_file)
        println("Warning: file not found: ", nc_file)
        return
    end

    lat_min, lat_max = lat_bounds
    lon_min, lon_max = lon_bounds

    ds = NCDatasets.Dataset(nc_file, "r");
    lats            = coalesce.(ds["latitude"][:,:], NaN32)
    lons            = coalesce.(ds["longitude"][:,:], NaN32)
    reflectance_vis = coalesce.(ds["reflectance_vis"][:,:], NaN32)
    tau             = coalesce.(ds["cloud_visible_optical_depth"][:,:], NaN32)
    # sw_albedo       = coalesce.(ds["broadband_shortwave_albedo"][:,:], NaN32)
    particle_size   = coalesce.(ds["cloud_particle_size"][:,:], NaN32)
    pixel_sza       = coalesce.(ds["pixel_sza"][:,:], NaN32)
    pixel_vza       = coalesce.(ds["pixel_vza"][:,:], NaN32)
    cloud_height    = coalesce.(ds["cloud_top_height"][:,:], NaN32)
    close(ds)

    # radiative calculations for each pixel
    tau_scaled = calc_tau_scaled.(particle_size, tau)
    albedo = albedo_kernel.(tau_scaled, pixel_sza, particle_size) # cloud albedo

    # ISCCP-consistent binary classification: cloudy if retrieval produced a valid tau
    cloud_mask = .!isnan.(tau) .& (tau .> 0.0)

    # geographic / solar zenith condition mask; include clear and cloudy pixels
    spatial_mask = ( (lats .>= lat_min) .& (lats .<= lat_max) .& 
                     (lons .>= lon_min) .& (lons .<= lon_max) .&
                     (pixel_vza .<= vza_thr) .& (pixel_sza .<= sza_thr) .& 
                     (.!isnan.(reflectance_vis)) )
    # valid subsets
    v_cloud = cloud_mask[spatial_mask]
    v_height = cloud_height[spatial_mask]
    v_tau   = tau_scaled[spatial_mask]
    v_albedo = albedo[spatial_mask]
    v_reflectance = reflectance_vis[spatial_mask]

    file_valid_footprints = length(v_cloud)
    if file_valid_footprints == 0 return end

    # Accumulate footprint contents
    for i in 1:file_valid_footprints
        total_all_count[1] += 1
        if !v_cloud[i] # clear
            clear_all_count[1] += 1
        else           # cloud
            cloud_all_count[1] += 1
            h_bin = get_bin_index(1e3*v_height[i], height_bins)
            if 1 <= h_bin <= nh # only accumulates below 4 km
                cloud_profile_count[h_bin] += 1
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
function compile_albedo_profile(lat_bounds, lon_bounds, data_file_list, fileout)
    height_bins = 0.0:10.0:4000.0
    nh = length(height_bins) # 401

    sza_thr = 75.0
    vza_thr = 75.0

    # --- PREALLOCATE ACCUMULATORS ---
    albedo_profile_accumulator = zeros(Float64, nh)
    reflec_profile_accumulator = zeros(Float64, nh)
    cloud_profile_count = zeros(Int64, nh)
    cloud_all_count = zeros(Int64, 1)
    clear_all_count  = zeros(Int64, 1)
    total_all_count  = zeros(Int64, 1) # denominator of all scanned pixels

    # with open output file...
    NCDatasets.Dataset(fileout, "a") do dsout
        dsout["cloud_top_height"][:] = height_bins
        # ...loop over input files
        for (fi, file) in enumerate(data_file_list)
            # zero accumulators each time
            albedo_profile_accumulator .= 0
            reflec_profile_accumulator .= 0
            cloud_profile_count .= 0
            cloud_all_count .= 0
            clear_all_count .= 0
            total_all_count .= 0 # denominator of all scanned pixels

            # update counts for one image file
            update_albedo_profile!(
                albedo_profile_accumulator,
                reflec_profile_accumulator,
                cloud_profile_count,
                cloud_all_count,
                clear_all_count,
                total_all_count,
                joinpath(datadir, "GOES/all", file), 
                lat_bounds, 
                lon_bounds,
                sza_thr, vza_thr, 
                height_bins )

            # println("clear_all_count: ", clear_all_count[1], " cloud_all_count: ", cloud_all_count[1], " total_all_count: ", total_all_count[1])

            # calculate fractions
            cloud_low_count = sum(cloud_profile_count)
            clear_and_low_count = clear_all_count[1] + cloud_low_count # not obscured
            albedo_profile = albedo_profile_accumulator ./ clear_and_low_count
            reflec_profile = reflec_profile_accumulator ./ clear_and_low_count
            cloud_profile  = cloud_profile_count        ./ clear_and_low_count
            cloud_total_frac = cloud_all_count[1] / total_all_count[1]
            clear_total_frac = clear_all_count[1] / total_all_count[1]

            # save one record each time to file
            dsout["time"][fi] = yyyydoy_HHMM2dt(data_file_list[fi][17:28])
            dsout["albedo_profile"][:, fi] .= albedo_profile
            dsout["reflec_profile"][:, fi] .= reflec_profile
            dsout["cloud_profile"][ :, fi] .=  cloud_profile
            dsout["cloud_total_frac"][fi] = cloud_total_frac
            dsout["clear_total_frac"][fi] = clear_total_frac
            dsout["total_all_count"][fi] = total_all_count[1]
        end
    end

    # return only the last record
    return (albedo_profile_accumulator, 
        reflec_profile_accumulator, 
        cloud_profile_count, 
        cloud_all_count, 
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
sza_thr = 75.0
vza_thr = 75.0
data_file_list = filter(daylight_file, readdir(joinpath(datadir, "GOES/all"))) # [1:3]
height_bins=0.0:10.0:4000.0

# Compile albedo and reflectance-weighted cloud counts
fileout =joinpath(datadir, "shcu_cloud_albedo_refl_profile_ts.nc")
(   albedo_profile_accumulator,
    reflec_profile_accumulator,
    cloud_profile_count,
    cloud_all_count,
    clear_all_count,
    total_all_count ) = compile_albedo_profile(
        lat_bounds, lon_bounds, data_file_list[:], fileout )

# ==============================================================================
# 4. NORMALIZE FRACTIONS below 4 km for last file
# ==============================================================================
cloud_low_count = sum(cloud_profile_count)
clear_and_low_count = clear_all_count[1] + cloud_low_count # not obscured
albedo_profile = albedo_profile_accumulator ./ clear_and_low_count
reflec_profile = reflec_profile_accumulator ./ clear_and_low_count
cloud_profile  = cloud_profile_count        ./ clear_and_low_count
cloud_total_frac = cloud_all_count[1] / total_all_count[1]
clear_total_frac = clear_all_count[1] / total_all_count[1]

begin
    println("Total pixel count: ", total_all_count[1])
    println("Regional pure clear pixel fraction: ", round(100*clear_total_frac, digits=2), " %")
    println("Regional cloudy pixel fraction: ", round(100*cloud_total_frac, digits=2), " %")
    println("Regional mean albedo: ", round(sum(albedo_profile), digits=3))
    println("Regional mean reflectance: ", round(sum(reflec_profile), digits=3))
    println("albedo/reflectance fraction ratio: ", round(sum(albedo_profile)/sum(reflec_profile), digits=3))
end

# record accumulations not stored or saved here !

height_bins = 0.0:10.0:4000.0
nh = length(height_bins)

# plot the first few profiles
clf()
for fi = 1:5
    Dataset(fileout, "r") do dsin
        # Get the index you want to read (e.g., the first record)

        # --- Read a single record matching your variables ---
        # (Note: time reads back directly as a Julia DateTime object)
        dt_parsed        = dsin["time"][fi]
        albedo_profile   = dsin["albedo_profile"][:, fi]
        reflec_profile   = dsin["reflec_profile"][:, fi]
        cloud_profile    = dsin["cloud_profile"][:, fi]
        cloud_total_frac = dsin["cloud_total_frac"][fi]
        clear_total_frac = dsin["clear_total_frac"][fi]
        total_all_count  = dsin["total_all_count"][fi]

        subplot(1,2,1)
        plot(100*albedo_profile, height_bins/1e3, linewidth=0.3, color="tab:blue", label="albedo-fraction")
        plot(100*reflec_profile, height_bins/1e3, linewidth=0.3, color="tab:orange", label="reflectance-fraction")
        plot(100*cloud_profile, height_bins/1e3, linewidth=0.3, color="k", label="cloud fraction")
        # legend(frameon=false)
        xlim([-0.0002, 0.6]); ylim([0, 4])
        xlabel("cloud amount per 10 m height bin (%)")
        ylabel("height (km)")
        subplot(1,2,2)
        plot(cumsum(reverse(albedo_profile)), reverse(height_bins)/1e3, linewidth=0.3, color="tab:blue", label="albedo-fraction")
        plot(cumsum(reverse(reflec_profile)), reverse(height_bins)/1e3, linewidth=0.3, color="tab:orange", label="reflectance-fraction")
        plot(cumsum(reverse(cloud_profile)), reverse(height_bins)/1e3, linewidth=0.3, color="k", label="0.1x cloud fraction")
        xlabel("cumulative cloud amount")
        xlim([-0.0005, 0.7]); ylim([0, 4])
    end
end

# [ savefig(joinpath(datadir, "shcu_cloud_albedo_refl_profile.$f")) for f in ["png", "pdf", "svg", "eps" ] ]

# read the whole file with the same variable names
(  time,
   cloud_top_height,
   albedo_profile,
   reflec_profile,
   cloud_profile,
   cloud_total_frac,
   clear_total_frac,
   total_all_count ) = NCDatasets.Dataset(fileout, "r") do dsin
        dsin["time"][:],
        dsin["cloud_top_height"][:],
        dsin["albedo_profile"][:,:],
        dsin["reflec_profile"][:,:],
        dsin["cloud_profile"][:,:],
        dsin["cloud_total_frac"][:],
        dsin["clear_total_frac"][:],
        dsin["total_all_count"][:]
end

# full record means, weighted by valid pixel counts
# for the appropriate denominators
clear_and_low_fraction = clear_total_frac .+ sum(cloud_profile, dims=1)[:]
lowrecmean(x) = x*clear_and_low_fraction / sum(clear_and_low_fraction)
albedo_profile_mean = lowrecmean(albedo_profile) # record mean profiles
reflec_profile_mean = lowrecmean(reflec_profile)
cloud_profile_mean  = lowrecmean(cloud_profile)

allrecmean(x) = x'*total_all_count / sum(total_all_count)
cloud_total_frac_mean = allrecmean(cloud_total_frac) # record mean fractions
clear_total_frac_mean = allrecmean(clear_total_frac)

# variance
lowrecvar(x) = lowrecmean(x.^2) .- lowrecmean(x).^2
albedo_profile_std = sqrt.(lowrecvar(albedo_profile))
reflec_profile_std = sqrt.(lowrecvar(reflec_profile))
cloud_profile_std  = sqrt.(lowrecvar(cloud_profile))
norm = 1/sqrt(length(time))

let albedo_profile = lowrecmean(albedo_profile),
    reflec_profile = lowrecmean(reflec_profile),
    cloud_profile  = lowrecmean(cloud_profile),
    albedo_profile_err = norm*sqrt.(lowrecvar(albedo_profile)),
    reflec_profile_err = norm*sqrt.(lowrecvar(reflec_profile)),
    cloud_profile_err  = norm*sqrt.(lowrecvar(cloud_profile))

    clf()
    subplot(1,2,1)
    plot(100*albedo_profile_mean, height_bins/1e3, linewidth=1, color="tab:blue", label="albedo-weighted")
    plot(100*reflec_profile_mean, height_bins/1e3, linewidth=1, color="tab:orange", label="reflectance-weighted")
    plot(0.5*100*cloud_profile_mean, height_bins/1e3, linewidth=1, color="k", label="0.5x cloud fraction")
    plot(100*(albedo_profile_mean .- albedo_profile_err), height_bins/1e3, linewidth=0.4, color="tab:blue")
    plot(100*(reflec_profile_mean .- reflec_profile_err), height_bins/1e3, linewidth=0.4, color="tab:orange")
    plot(0.5*100*(cloud_profile_mean .- cloud_profile_err), height_bins/1e3, linewidth=0.4, color="k")
    plot(100*(albedo_profile_mean .+ albedo_profile_err), height_bins/1e3, linewidth=0.4, color="tab:blue")
    plot(100*(reflec_profile_mean .+ reflec_profile_err), height_bins/1e3, linewidth=0.4, color="tab:orange")
    plot(0.5*100*(cloud_profile_mean .+ cloud_profile_err), height_bins/1e3, linewidth=0.4, color="k")
    legend(frameon=false)
    xlim([-0.0002, 0.16]); ylim([0, 4])
    xlabel("cloud amount per 10 m height bin (%)")
    ylabel("height (km)")
    subplot(1,2,2)
    plot(cumsum(reverse(albedo_profile_mean)), reverse(height_bins)/1e3, linewidth=1, color="tab:blue", label="albedo-weighted")
    plot(cumsum(reverse(reflec_profile_mean)), reverse(height_bins)/1e3, linewidth=1, color="tab:orange", label="reflectance-weighted")
    plot(0.5*cumsum(reverse(cloud_profile_mean)), reverse(height_bins)/1e3, linewidth=1, color="k", label="0.5x cloud fraction")
    plot(cumsum(reverse(albedo_profile_mean .- albedo_profile_err)), reverse(height_bins)/1e3, linewidth=0.4, color="tab:blue")
    plot(cumsum(reverse(reflec_profile_mean .- reflec_profile_err)), reverse(height_bins)/1e3, linewidth=0.4, color="tab:orange")
    plot(0.5*cumsum(reverse(cloud_profile_mean .- cloud_profile_err)), reverse(height_bins)/1e3, linewidth=0.4, color="k")
    plot(cumsum(reverse(albedo_profile_mean .+ albedo_profile_err)), reverse(height_bins)/1e3, linewidth=0.4, color="tab:blue")
    plot(cumsum(reverse(reflec_profile_mean .+ reflec_profile_err)), reverse(height_bins)/1e3, linewidth=0.4, color="tab:orange")
    plot(0.5*cumsum(reverse(cloud_profile_mean .+ cloud_profile_err)), reverse(height_bins)/1e3, linewidth=0.4, color="k")
    xlabel("cumulative cloud amount")
    xlim([-0.0005, 0.25]); ylim([0, 4])
end