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
function albedo_kernel(r_e, tau_scaled)
    omgts = (1-asym_parameter(r_e)) * tau_scaled
    omgts / (2 + omgts)
end


function get_bin_index(value, edges)
    if isnan(value) || value < edges[1] || value > edges[end] return 0 end
    return searchsortedlast(edges, value)
end

# ==============================================================================
# 2. FUNCTIONS TO ACCUMULATE PIXELS IN-PLACE 
# ==============================================================================

"""
    update_isccp_accumulation!(acc_cloudy, acc_clear_profile, acc_clear_scalar, global_denominator, nc_file, ...)
Processes a single granule file. Allocates 100% clear space to a single scalar accumulator,
while splitting subpixel cloud structures into the 2D joint matrix and 1D pressure profile.
All trackers are modified strictly in place.
"""
function update_isccp_accumulation!(
    acc_sum_albedo::Vector{Float64},
    acc_cloudy::Matrix{Float64}, 
    acc_clear_profile::Vector{Float64}, 
    acc_clear_scalar::Vector{Float64}, # Single scalar inside a 1-element container
    global_denominator::Vector{Int64}, # Single element container for global normalization base
    nc_file::String,
    lat_bounds::Tuple{Float64, Float64},
    lon_bounds::Tuple{Float64, Float64},
    tau_edges::Vector{T},
    pc_edges::Vector{T};
    vza_thr = 60.0, 
    sza_thr = 60.0 ) where T <: Real

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
    pc              = coalesce.(ds["cloud_top_pressure"][:,:], NaN32)

    close(ds)

    tau_scaled = calc_tau_scaled.(particle_size, tau)
    albedo = albedo_kernel.(particle_size, tau_scaled)

    # ISCCP-consistent binary classification: cloudy if retrieval produced a valid tau
    cloud_mask = .!isnan.(tau) .& (tau .> 0.0)

    # Universal geographic / solar zenith condition mask
    spatial_mask = (lats .>= lat_min) .& (lats .<= lat_max) .& 
                   (lons .>= lon_min) .& (lons .<= lon_max) .&
                   (pixel_vza .<= vza_thr) .& (pixel_sza .<= sza_thr) .& 
                   (.!isnan.(reflectance_vis))

    # valid subsets
    v_cloud = cloud_mask[spatial_mask]
    v_pc    = pc[spatial_mask]
    v_tau   = tau_scaled[spatial_mask]
    v_albedo = albedo[spatial_mask]

    file_valid_footprints = length(v_cloud)
    if file_valid_footprints == 0 return end

    # Accumulate footprint contents down the track
    for i in 1:file_valid_footprints

        if !v_cloud[i] # clear
            acc_clear_scalar[1] += 1.0
        else           # cloud retrieval
            acc_sum_albedo[1] += v_albedo[i]
            p_bin = get_bin_index(v_pc[i], pc_edges)
            if p_bin > 0
                tau_bin = get_bin_index(v_tau[i], tau_edges)
                if tau_bin > 0
                    acc_cloudy[p_bin, tau_bin] += 1.0
                end
            end
        end
    end

    # Increment master footprint tracker in place
    global_denominator[1] += file_valid_footprints
    println("$nc_file ($file_valid_footprints pixels)")
end

function compile_isccp_histogram(lat_bounds, lon_bounds, data_file_list)
    #tau_edges = [0.0, 0.3, 1.3, 3.6, 9.4, 23.0, 60.0, 1000.0]
    #pc_edges  = [10.0, 180.0, 310.0, 440.0, 560.0, 680.0, 800.0, 1000.0]
    # Zelinka kernel bin edges:
    tau_edges = [0.01, 0.3, 1.3, 3.6, 9.4, 23, 60, 380]
    pc_edges = reverse(float([1000, 800, 680, 560, 440, 310, 180, 50])) 
    # ! pc_edges are reversed relative to Zelinka kernel !

    # --- PREALLOCATE ACCUMULATORS ---
    albedo_mean_accumulator = zeros(Float64, 1)
    isccp_cloudy_accumulator = zeros(Float64, 7, 7)
    isccp_clear_profile_acc  = zeros(Float64, 7)
    global_clear_pixel_counter = zeros(Float64, 1) # Single explicit scalar counter for completely clear space
    global_total_footprints    = zeros(Int64, 1)   # Global denominator tracking all scanned pixels

    # Run state ingestion loop sequence
    for file in data_file_list
        update_isccp_accumulation!(
            albedo_mean_accumulator,
            isccp_cloudy_accumulator, 
            isccp_clear_profile_acc,
            global_clear_pixel_counter,
            global_total_footprints, 
            joinpath(datadir, "GOES/all", file), 
            lat_bounds, 
            lon_bounds, 
            tau_edges, 
            pc_edges )
    end

    return (albedo_mean_accumulator, isccp_cloudy_accumulator, isccp_clear_profile_acc, global_clear_pixel_counter, global_total_footprints)
end

# ZELINKA FEEDBACK KERNEL FUNCTIONS
"""
K_sw_cloudy, K_lw_cloudy = load_zelinka_kernel(kernel_file)
plev_bnds = [1000.0  800.0  680.0  560.0  440.0  310.0  180.0  50.0]
tau_bnds = [0.01  0.3  1.3  3.6  9.4  23.0  60.0  380.0]
"""
function load_zelinka_kernel( kernel_file=joinpath(datadir, "obs_cloud_kernels4.nc") )
    NCDatasets.Dataset(kernel_file, "r") do kd
    K_sw_cloudy = kd["SWkernel"][:,:,:,:,:] # -> (albcs, lat, plev, tau, time)
    K_lw_cloudy = kd["LWkernel"][:,:,:,:]
    return (K_sw_cloudy, K_lw_cloudy)
    end
end
# lowest 2 cloud levels are indexed plev[1:2], height incresing, plev is decreasing

# ==============================================================================
# 3. ENVIRONMENT STATE PREALLOCATION
# ==============================================================================

# EUREC4A region
lat_bounds = (12.5, 16.0)
lon_bounds = (-60.0, -49.0)
# files matching data
daylight_file(s) = 1200 <= parse(Int,match(r"(\d{4})\.PX\.02K\.NC$", s).captures[1]) <= 1920
data_file_list = filter(daylight_file, readdir(joinpath(datadir, "GOES/all"))) # [1:3]

# Set up the histogram netcdf file.
# first dump and edit the obs kernel cdl file, then
# ncgen -o shcu_isccp_cloud_pct.nc shcu_isccp_cloud_pct.cdl
# copy the variables we want
# ncks -A -C -v plev_bnds,tau_bnds,plev,tau obs_cloud_kernels4.nc shcu_isccp_cloud_pct.nc
filename = joinpath(datadir, "shcu_isccp_cloud_pct4.nc")
if false
    # Compile the ISCCP histograms across the domain track
    (   albedo_mean_accumulator,
        isccp_cloudy_accumulator, 
        isccp_clear_profile_acc, 
        global_clear_pixel_counter, 
        global_total_footprints ) = compile_isccp_histogram(
            lat_bounds, lon_bounds, data_file_list )

    # ==============================================================================
    # 4. NORMALIZATION FOR RAD KERNELS
    # ==============================================================================
    denom = global_total_footprints[1]
    albedo_mean = albedo_mean_accumulator[1] / denom

    # histogram percentage matrices
    isccp_cloudy_histogram_pct = (isccp_cloudy_accumulator ./ denom) .* 100.0
    isccp_clear_profile_pct    = (isccp_clear_profile_acc  ./ denom) .* 100.0 # all zeros
    total_domain_clear_sky_pct = (global_clear_pixel_counter / denom) * 100.0

    println("Total pixel count: ", denom)
    println("Regional pure clear pixel fraction: ", round(total_domain_clear_sky_pct[1], digits=2), " %")
    println("Regional mean albedo: ", round(albedo_mean, digits=3))

    # write the cloud histogram data
    # plev is reversed relative to kernel data!
    NCDatasets.Dataset(filename, "a") do ihd
        ihd["albedo"][1]           = albedo_mean
        ihd["cloud_hist"][:,:]    .= isccp_cloudy_histogram_pct
        ihd["clear_prof"][:]      .= isccp_clear_profile_pct
        ihd["clear_pixel_frac"][1] = total_domain_clear_sky_pct
    end
else
    # load data from the netcdf file
    (albedo_mean, 
    isccp_cloudy_histogram_pct, 
     isccp_clear_profile_pct, 
     total_domain_clear_sky_pct) = NCDatasets.Dataset(filename, "r") do ihd
        ihd["albedo"][1],
        ihd["cloud_hist"][:,:],
        ihd["clear_prof"][:],
        ihd["clear_pixel_frac"][1] 
    end
end

# now use kernels to compute CRE or feedbacks across the cloudy and clear histogram structures
K_sw_cloudy, K_lw_cloudy = load_zelinka_kernel();

# average kernels in Jan-Feb, lat_range, lon_range
# SW(time, tau, plev, lat, albcs) 
# netcdf structure reverses to (albcs, lat, plev, tau, time) in julia
ll = 42:43  # searchsortedfirst(lat,12.0):searchsortedlast(lat,16.5) # lat range for EUREC4A
aa = 1   # searchsortedlast(albcs,0.05) # clear sky albedo approx 0
tt = 1:2
# K_sw_sc_cloudy = mean(K_sw_cloudy[tt,:,:,ll,aa], dims=(1,4)) |> dropdims
# K_lw_sc_cloudy = mean(K_lw_cloudy[tt,:,:,ll   ], dims=(1,4)) |> dropdims
K_sw_sc_cloudy = mean(K_sw_cloudy[aa, ll, :,:, tt], dims=(1,4))[1,:,:,1] # plev x tau
K_lw_sc_cloudy = mean(K_lw_cloudy[    ll, :,:, tt], dims=(1,4))[1,:,:,1]
# scale kernel bins by observed cloud fractions
ich = coalesce.(reverse(isccp_cloudy_histogram_pct, dims=1), NaN)
sw_cre_hist = ich .* K_sw_sc_cloudy
lw_cre_hist = ich .* K_lw_sc_cloudy
# would sum to get total CRE

# low cloud CRE (lowest 2 height levels, greatest pressure levels indexed 1:2)
sw_low_cre, lw_low_cre = sum(sw_cre_hist[1:2, :]), sum(lw_cre_hist[1:2, :]) # W/m^2

# compute Delta CRE W/m^2/K for +5 % decrease in lowest 2 cloud levels 
dR_sw_ShCu, dR_lw_ShCu = -0.05 .* (sw_low_cre, lw_low_cre)

# plots
tau_edges = [0.01, 0.3, 1.3, 3.6, 9.4, 23, 60, 380]
pc_edges = float([1000.0, 800, 680, 560, 440, 310, 180, 50]) # correct order

function format_prec(val; sigdigits=2)
    r = round(val; sigdigits=sigdigits)
    # If it ends in .0 (like 26.0 or 380.0), convert to Int to drop the decimal
    return isinteger(r) ? string(Int(r)) : string(r)
end

"""
    plot_isccp_matrix(data, tau_edges, pc_edges)

Plots an ISCCP 2D histogram where each cell is rendered with equal visual size,
regardless of the underlying physical values of the bins.

- `data`: A 2D Matrix of dimensions (length(pc_edges)-1, length(tau_edges)-1)
- `tau_edges`: Array of optical thickness boundaries
- `pc_edges`: Array of cloud top pressure boundaries (sorted high-to-low or low-to-high)
"""
function plot_isccp_matrix(ax, data::Matrix, tau_edges=tau_edges, pc_edges=pc_edges; 
    cmap=ColorMap("Blues"), kwargs...)

    fig = ax.figure
    # Get dimensions based on the edge vectors
    num_tau_bins = length(tau_edges) - 1
    num_press_bins = length(pc_edges) - 1
    
    # Validate data matrix dimensions
    if size(data) != (num_press_bins, num_tau_bins)
        error("Data dimensions $(size(data)) must match matrix bins ($num_press_bins, $num_tau_bins)")
    end

    # Force equal visual spacing by creating an identical, sequential grid index
    x_grid = 0:num_tau_bins
    y_grid = 0:num_press_bins

    # Plot using the mock sequential grid to force uniform pixel/square sizing
    pcm = ax.pcolormesh(
        x_grid, 
        y_grid, 
        data, 
        shading="flat", 
        edgecolors="white", 
        linewidth=1,
        cmap=cmap,
        kwargs...
    )

    cm   = pcm.cmap
    norm = pcm.norm

    # Add text annotations centered inside each bin
    # Determine maximum value to dynamically scale text color contrast
    max_val = maximum(data)
    for i in 1:num_press_bins
        for j in 1:num_tau_bins
            val = data[i, j]

            r, g, b, _ = PythonPlot.PythonCall.pyconvert(
                Tuple{Float64, Float64, Float64, Float64}, cm(norm(val)) )
            lum = 0.299*r + 0.587*g + 0.114*b
            # Dynamically set text color based on matrix shading depth for legibility
            #text_color = val > (max_val * 0.5) ? "white" : "black"
            text_color = (lum > 0.5) ? "black" : "white"

            # Find the center coordinates for the current cell
            x_center = (x_grid[j] + x_grid[j+1]) / 2
            y_center = (y_grid[i] + y_grid[i+1]) / 2
            
            # Formats value to 2 decimal places (adjust as needed)
            rv = round(val, digits=1)
            label_text = rv==0.0 ? "" : @sprintf("%.2g", rv)

            ax.text(
                x_center, 
                y_center, 
                label_text,
                ha="center", 
                va="center", 
                color=text_color,
                fontsize=9,
                # weight="bold"
            )
        end
    end

    # Label the ticks exactly at the bin edges
    ax.set_xticks(x_grid)
    ax.set_xticklabels(string.(format_prec.(tau_edges)))
    ax.set_yticks(y_grid)
    ax.set_yticklabels(string.(round.(Int, pc_edges)))

    # Format labels
    ax.set_xlabel("cloud optical thickness (τ)")
    ax.set_ylabel("cloud top pressure (hPa)")
    # ax.set_title("ISCCP joint cloud histogram")

    # Add contextual color bar
    cbar = fig.colorbar(pcm, ax=ax)
    # cbar.set_label("Cloud Fraction / Frequency")

    fig.tight_layout()
    return pcm, cbar
end
function plot_isccp_matrix(data::Matrix, tau_edges=tau_edges, pc_edges=pc_edges;
    cmap=ColorMap("Blues"), kwargs...)
    ax = gca()
    return plot_isccp_matrix(ax, data, tau_edges, pc_edges, cmap=cmap, kwargs...)
end

# initialize figure
fig, axs = subplots(2, 1, figsize=(5, 6))
axs[0].invert_yaxis() # Invert y-axis to have higher pressures at the bottom
axs[1].invert_yaxis() # Invert y-axis to have higher pressures at the
# axs[2].invert_yaxis() # Invert y-axis to have higher pressures at the
pclf, cb1 = plot_isccp_matrix(axs[0], reverse(ich, dims=1), tau_edges, reverse(pc_edges),
    cmap=ColorMap("Blues").resampled(10))
pcre, cb2 = plot_isccp_matrix(axs[1], reverse(sw_cre_hist, dims=1), tau_edges, reverse(pc_edges),
    cmap=ColorMap("Blues_r").resampled(10))
# pcre.set_cmap("Blues_r")
axs[0].set_title("GOES cloud fraction (%)")
axs[1].set_title("SW CRE (W m⁻²)")
axs[0].set_xlabel(nothing)
# axs[1].set_ylabel(nothing)
tight_layout()
display(fig)

for fmt in ["png", "svg", "pdf"]
    fig.savefig("goes_isccp_cloud_histogram_and_cre.$fmt")
end

# print out the CREs
println("SW, LW low cloud radiative effects: $(sw_low_cre), $(lw_low_cre)")
println("SW, LW low cloud radiative effect change: $(dR_sw_ShCu), $(dR_lw_ShCu)")

# get CERES FBCT histogram for comparison
# /Users/deszoeks/Data/ATOMIC/satellite/ISCCP/CER_FluxByCldTyp-Month_Terra-Aqua-MODIS_Edition4A-20260521_024852"
# compile_fbct_hist.jl
fbct_histogram = [
 7.70293   10.2463      3.19888      0.521623     0.0567247   0.00132773
 1.50099    5.29244     5.71525      1.81545      0.304576    0.0157222
 0.642293   0.861229    1.28077      0.967155     0.422698    0.042823
 0.304027   0.252883    0.557532     0.511711     0.314445    0.042039
 2.05153    1.4502      1.62918      1.10968      0.59592     0.0858304
 3.54012    1.54909     0.831532     0.686792     0.704872    0.291805
 0.159219   0.00309147  0.000219069  0.000448272  0.00240926  0.0101675 ]
fbct_od = [ 0.02, 1.27, 3.55, 9.38, 22.63, 60.36, 378.65 ]
fbct_press = [ 1000.0, 800.0, 680.0, 560.0, 440.0, 310.0, 180.0, 10.0 ]
fbct_sw_cre_kernel = [
 10.2211   31.0894  59.0229   91.5348  118.163  139.476
 14.8665   36.8831  64.7321   95.9388  123.27   144.842
 12.5879   39.0152  70.0865  102.226   128.906  155.732
  7.91026  50.8679  81.5521  110.729   134.987  161.091
 15.5884   57.0098  91.7088  123.456   147.007  167.925
 13.1715   53.8972  92.0407  130.11    155.798  174.39
  4.31043  38.7438  80.809   125.422   154.941  173.838 ]
fbct_sw_cre_histogram = -0.01*fbct_histogram.*fbct_sw_cre_kernel

# add the FBCT histogram to the plot
fig, axs = subplots(2, 2, figsize=(8, 6))
axs[0,0].invert_yaxis() # Invert y-axis to have higher pressures at the bottom
axs[0,1].invert_yaxis() # Invert y-axis to have higher pressures at the
axs[1,0].invert_yaxis() # Invert y-axis to have higher pressures at the
axs[1,1].invert_yaxis() # Invert y-axis to have higher pressures at the
pclf, cb1 = plot_isccp_matrix(axs[0,0], reverse(ich, dims=1), tau_edges, reverse(pc_edges),
    cmap=ColorMap("Blues").resampled(10))
pcre, cb2 = plot_isccp_matrix(axs[1,0], reverse(sw_cre_hist, dims=1), tau_edges, reverse(pc_edges),
    cmap=ColorMap("Blues_r").resampled(10))
pfbc, cb3 = plot_isccp_matrix(axs[0,1], reverse(fbct_histogram, dims=1), fbct_od, reverse(fbct_press),
    cmap=ColorMap("Blues").resampled(10))
pfsw, cb4 = plot_isccp_matrix(axs[1,1], reverse(fbct_sw_cre_histogram, dims=1), fbct_od, reverse(fbct_press),
    cmap=ColorMap("Blues_r").resampled(10))

axs[0,0].set_title("GOES cloud fraction (%)")
axs[1,0].set_title("GOES SW CRE (W m⁻²)")
axs[1,1].set_title("CERES FBCT SW CRE (W m⁻²)")
axs[0,1].set_title("CERES FBCT cloud fraction (%)")
axs[0,1].set_ylabel(nothing)
axs[0,1].set_xlabel(nothing)
axs[0,0].set_xlabel(nothing)
[ ax.set_aspect(0.8) for ax in axs.flatten() ]  # Ensure equal aspect ratio for all subplots

# panel labels
axs[0,0].text(6.4, 0.6, "a", size=12)
axs[1,0].text(6.4, 0.6, "c", size=12)
axs[0,1].text(5.4, 0.6, "b", size=12)
axs[1,1].text(5.4, 0.6, "d", size=12)

tight_layout()
display(fig)

[ fig.savefig("goes_and_fbct_cloud_and_cre_histograms.$fmt") for fmt in ["png", "svg", "pdf"] ]