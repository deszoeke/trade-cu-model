using Pkg; Pkg.activate(".")
using NCDatasets
using Statistics
using PythonPlot
using Printf

datadir = "../../../ATOMIC_GOES/data/"

# ==============================================================================
# 1. RADIATIVE TRANSFER LOOKUPS
# ==============================================================================

function compute_subpixel_cloud_fraction(alpha_obs, sza_deg, tau, cth_km, alpha_surface)
    if isnan(tau) || isnan(alpha_obs) || isnan(sza_deg) || isnan(cth_km) || isnan(alpha_surface)
        return NaN
    end
    if tau <= 0.0 || alpha_obs <= alpha_surface
        return 0.0
    end
    mu_0 = cos(deg2rad(sza_deg))
    if mu_0 <= 0.0 return 0.0 end

    g = cth_km > 6.0 ? 0.75 : 0.85
    gamma_1 = 0.75 * (1.0 - (g * mu_0))
    alpha_cloud = (gamma_1 * tau) / (1.0 + (gamma_1 * tau))
    
    if alpha_cloud <= alpha_surface return 0.0 end
    return clamp((alpha_obs - alpha_surface) / (alpha_cloud - alpha_surface), 0.0, 1.0)
end

function get_bin_index(value, edges)
    if isnan(value) || value < edges[1] || value > edges[end] return 0 end
    return searchsortedlast(edges, value)
end

# ==============================================================================
# 2. FUNCTIOMS TO ACCUMULATE PIXELS IN-PLACE 
# ==============================================================================

"""
    update_isccp_accumulation!(acc_cloudy, acc_clear_profile, acc_clear_scalar, global_denominator, nc_file, ...)
    
Processes a single granule file. Allocates 100% clear space to a single scalar accumulator,
while splitting subpixel cloud structures into the 2D joint matrix and 1D pressure profile.
All trackers are modified strictly in place.
"""
function update_isccp_accumulation!(
    acc_cloudy::Matrix{Float64}, 
    acc_clear_profile::Vector{Float64}, 
    acc_clear_scalar::Vector{Float64}, # Single scalar inside a 1-element container
    global_denominator::Vector{Int64}, # Single element container for global normalization base
    nc_file::String,
    lat_bounds::Tuple{Float64, Float64},
    lon_bounds::Tuple{Float64, Float64},
    tau_edges::Vector{T},
    pc_edges::Vector{T} ) where T <: Real

    if !isfile(nc_file)
        println("Warning: file not found: ", nc_file)
        return
    end

    lat_min, lat_max = lat_bounds
    lon_min, lon_max = lon_bounds

    ds = NCDatasets.Dataset(nc_file, "r")
    
    lats            = coalesce.(ds["latitude"][:,:], NaN)
    lons            = coalesce.(ds["longitude"][:,:], NaN)
    reflectance_vis = coalesce.(ds["reflectance_vis"][:,:], NaN)
    alpha_obs       = coalesce.(ds["broadband_shortwave_albedo"][:,:], NaN)
    tau             = coalesce.(ds["cloud_visible_optical_depth"][:,:], NaN)
    cth_km          = coalesce.(ds["cloud_top_height"][:,:], NaN)
    pixel_sza       = coalesce.(ds["pixel_sza"][:,:], NaN)
    alpha_surface   = coalesce.(ds["clearsky_vis_reflectance"][:,:], NaN) # different grid
    pixel_vza       = coalesce.(ds["pixel_vza"][:,:], NaN)
    pc              = coalesce.(ds["cloud_top_pressure"][:,:], NaN)

    close(ds) 

    # Compute cloud fraction
    cloud_fraction = compute_subpixel_cloud_fraction.(alpha_obs, pixel_sza, tau, cth_km, 0.0)

    # Universal geographic / solar zenith condition mask
    spatial_mask = (lats .>= lat_min) .& (lats .<= lat_max) .& 
                   (lons .>= lon_min) .& (lons .<= lon_max) .&
                   (pixel_vza .<= 60.0) .& (pixel_sza .<= 60.0) .& 
                   (.!isnan.(reflectance_vis))

    v_cf   = cloud_fraction[spatial_mask]
    v_pc   = pc[spatial_mask]
    v_tau  = tau[spatial_mask]
    
    file_valid_footprints = length(v_cf)
    if file_valid_footprints == 0 return end

    # Accumulate footprint contents down the track
    for i in 1:file_valid_footprints
        f = v_cf[i]
        
        # --- BRANCH 1: Completely Clear Pixel (Skipped by retrieval / has no cloud properties) ---
        if isnan(f) || f <= 0.0
            acc_clear_scalar[1] += 1.0  # Accumulate cleanly into the single standalone scalar tracker
            
        # --- BRANCH 2: Footprint contains retrieved cloud layer columns ---
        else
            clear_fraction  = 1.0 - f
            cloudy_fraction = f
            
            p_bin = get_bin_index(v_pc[i], pc_edges)
            if p_bin > 0
                # Record the clear sky background remnant of this cloud footprint at its cloud top pressure
                acc_clear_profile[p_bin] += clear_fraction
                
                # Record the cloudy fraction into the 2D joint histogram matrix
                tau_bin = get_bin_index(v_tau[i], tau_edges)
                if tau_bin > 0 
                    acc_cloudy[p_bin, tau_bin] += cloudy_fraction 
                end
            end
        end
    end

    # Increment master footprint tracker in place
    global_denominator[1] += file_valid_footprints
    println("Processed file: $nc_file (Added $file_valid_footprints footprints to global pool).")
end

function compile_isccp_histogram(lats_range, lons_range, data_file_list)
    #tau_edges = [0.0, 0.3, 1.3, 3.6, 9.4, 23.0, 60.0, 1000.0]
    #pc_edges  = [10.0, 180.0, 310.0, 440.0, 560.0, 680.0, 800.0, 1000.0]
    # Zelinke kernel bin edges:
    tau_edges = [0.01, 0.3, 1.3, 3.6, 9.4, 23, 60, 380]
    pc_edges = reverse(float([1000, 800, 680, 560, 440, 310, 180, 50]))

    # --- PREALLOCATE ACCUMULATORS ---
    isccp_cloudy_accumulator = zeros(Float64, 7, 7)
    isccp_clear_profile_acc  = zeros(Float64, 7)
    global_clear_pixel_counter = zeros(Float64, 1) # Single explicit scalar counter for completely clear space
    global_total_footprints    = zeros(Int64, 1)   # Global denominator tracking all scanned pixels

    # Run state ingestion loop sequence
    for file in data_file_list
        update_isccp_accumulation!(
            isccp_cloudy_accumulator, 
            isccp_clear_profile_acc,
            global_clear_pixel_counter,
            global_total_footprints, 
            joinpath(datadir, "GOES/all", file), 
            lats_range, 
            lons_range, 
            tau_edges, 
            pc_edges )

    end

    return (isccp_cloudy_accumulator, isccp_clear_profile_acc, global_clear_pixel_counter, global_total_footprints)
end

# ZELINKA FEEDBACK KERNEL FUNCTIONS
function load_zelinka_kernel( kernel_file=joinpath(datadir, "obs_cloud_kernels4.nc") )
    NCDatasets.Dataset(kernel_file, "r") do kd
    K_sw_cloudy = kd["SWkernel"][:,:,:,:,:] # -> (albcs, lat, plev, tau, time)
    K_lw_cloudy = kd["LWkernel"][:,:,:,:]
    return (K_sw_cloudy, K_lw_cloudy)
    end
end

# not used:
function compute_zelinka_feedback(
    dR_cloudy::Matrix{Float64}, 
    dR_clear_prof::Vector{Float64}, 
    dR_clear_scalar::Float64,
    K_sw_cloudy, 
    K_lw_cloudy )
    
    # 1. Component sums from the subpixel profiles
    sw_cloudy = sum(dR_cloudy .* K_sw_cloudy)
    lw_cloudy = sum(dR_cloudy .* K_lw_cloudy)
    
    return (SW=sw_cloudy, LW=sw_cloudy, NET=total_sw + total_lw)
end

# ==============================================================================
# 3. ENVIRONMENT STATE PREALLOCATION
# ==============================================================================

# EUREC4A region
lats_range = (12.5, 16.0)
lons_range = (-60.0, -49.0)
# Array playlist matching your data structure
data_file_list = filter(endswith(".NC"), readdir(joinpath(datadir, "GOES/all")))[37:39]
# Compile the ISCCP histograms across the domain track
(   isccp_cloudy_accumulator, 
    isccp_clear_profile_acc, 
    global_clear_pixel_counter, 
    global_total_footprints ) = compile_isccp_histogram(
        lats_range, lons_range, data_file_list )

# ==============================================================================
# 4. NORMALIZATION FOR RAD KERNELS
# ==============================================================================
denom = global_total_footprints[1]

# histogram percentage matrices
isccp_cloudy_histogram_pct = (isccp_cloudy_accumulator ./ denom) .* 100.0
isccp_clear_profile_pct    = (isccp_clear_profile_acc  ./ denom) .* 100.0
total_domain_clear_sky_pct = (global_clear_pixel_counter / denom) * 100.0

println("Total pixel count: ", denom)
println("Regional pure clear pixel fraction: ", round(total_domain_clear_sky_pct[1], digits=2), " %")

# now use kernels to compute CRE or feedbacks across the cloudy and clear histogram structures
K_sw_cloudy, K_lw_cloudy = load_zelinka_kernel()

# average kernels in Jan-Feb, lat_range, lon_range
# SW(time, tau, plev, lat, albcs) 
# netcdf structure reverses to (albcs, lat, plev, tau, time) in julia
ll = 41 # searchsortedlast(lat,12.0)
aa = 1 # searchsortedlast(albcs,0.05) # clear sky albedo approx 0
tt = 1:2
# K_sw_sc_cloudy = mean(K_sw_cloudy[tt,:,:,ll,aa], dims=(1,4)) |> dropdims
# K_lw_sc_cloudy = mean(K_lw_cloudy[tt,:,:,ll   ], dims=(1,4)) |> dropdims
K_sw_sc_cloudy = mean(K_sw_cloudy[aa, ll, :,:, tt], dims=3)
K_lw_sc_cloudy = mean(K_lw_cloudy[    ll, :,:, tt], dims=3)
# scale kernel bins by observed cloud fractions
sw_cre_hist = isccp_cloudy_histogram_pct .* K_sw_sc_cloudy
lw_cre_hist = isccp_cloudy_histogram_pct .* K_lw_sc_cloudy
# would sum to get total CRE

# compute Delta CRE W/m^2/K for +5 % decrease in lowest 2 cloud levels 
dR_sw_ShCu = -0.05 * sw_cre_hist[(end-1):end, :] # apply to lowest 2 pressure levels
dR_lw_ShCu = -0.05 * lw_cre_hist[(end-1):end, :]

dSWCRE_ShCu = sum(dR_sw_ShCu) # local W/m^2/K sum over lowest 2 pressure levels and all tau bins

# plots
tau_edges = [0.01, 0.3, 1.3, 3.6, 9.4, 23, 60, 380]
pc_edges = reverse(float([1000, 800, 680, 560, 440, 310, 180, 50]))

"""
    plot_isccp_matrix(data, tau_edges, pc_edges)

Plots an ISCCP 2D histogram where each cell is rendered with equal visual size,
regardless of the underlying physical values of the bins.

- `data`: A 2D Matrix of dimensions (length(pc_edges)-1, length(tau_edges)-1)
- `tau_edges`: Array of optical thickness boundaries
- `pc_edges`: Array of cloud top pressure boundaries (sorted high-to-low or low-to-high)
"""

function plot_isccp_matrix(data::Matrix, tau_edges=tau_edges, pc_edges=pc_edges;
    cmap=ColorMap("Blues"), kwargs...)
    ax = gca()
    return plot_isccp_matrix(ax, data, tau_edges, pc_edges, cmap=cmap, kwargs...)
end
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
            label_text = @sprintf("%.2f", val) 

            ax.text(
                x_center, 
                y_center, 
                label_text,
                ha="center", 
                va="center", 
                color=text_color,
                fontsize=9,
                weight="bold"
            )
        end
    end

    # Label the ticks exactly at the bin edges
    ax.set_xticks(x_grid)
    ax.set_xticklabels(string.(tau_edges))
    ax.set_yticks(y_grid)
    ax.set_yticklabels(string.(pc_edges))

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

# initialize figure
fig, axs = subplots(1, 2, figsize=(8, 3))
axs[0].invert_yaxis() # Invert y-axis to have higher pressures at the bottom
axs[1].invert_yaxis() # Invert y-axis to have higher pressures at the
pclf, cb1 = plot_isccp_matrix(axs[0], isccp_cloudy_histogram_pct, tau_edges, pc_edges)
pcre, cb2 = plot_isccp_matrix(axs[1], sw_cre_hist[:,:,1], tau_edges, pc_edges,
    cmap=ColorMap("Blues_r"))
# pcre.set_cmap("Blues_r")
axs[0].set_title("GOES cloud fraction (%)")
axs[1].set_title("SW CRE (W m⁻²)")
axs[1].set_ylabel(nothing)
display(fig)