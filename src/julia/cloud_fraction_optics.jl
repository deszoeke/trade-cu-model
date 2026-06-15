using NCDatasets

"""
    compute_subpixel_cloud_fraction(alpha_obs, sza_deg, tau, cth_km, alpha_surface)
"""
function compute_subpixel_cloud_fraction(alpha_obs, sza_deg, tau, cth_km, alpha_surface)
    if isnan(tau) || isnan(alpha_obs) || isnan(sza_deg) || isnan(cth_km) || isnan(alpha_surface)
        return NaN
    end
    if tau <= 0.0 || alpha_obs <= alpha_surface
        return 0.0
    end
    mu_0 = cos(deg2rad(sza_deg))
    if mu_0 <= 0.0
        return 0.0
    end

    g = cth_km > 6.0 ? 0.75 : 0.85
    gamma_1 = 0.75 * (1.0 - (g * mu_0))
    alpha_cloud = (gamma_1 * tau) / (1.0 + (gamma_1 * tau))
    
    if alpha_cloud <= alpha_surface
        return 0.0
    end
        
    f = (alpha_obs - alpha_surface) / (alpha_cloud - alpha_surface)
    return clamp(f, 0.0, 1.0)
end

"""
    get_bin_index(value, edges)
"""
function get_bin_index(value, edges)
    if isnan(value) || value < edges || value > edges[end]
        return 0
    end
    return searchsortedlast(edges, value)
end

"""
    estimate_clear_sky_pressure(skin_temp_k)
"""
function estimate_clear_sky_pressure(skin_temp_k)
    if isnan(skin_temp_k)
        return 950.0 
    end
    dt = skin_temp_k - 288.15
    estimated_p = 1013.25 * (1.0 - (0.0065 * dt / 288.15))^(-5.257)
    return clamp(estimated_p, 15.0, 1010.0)
end

# --- Configuration & Constraints ---
nc_file = "goes16_nasa_retreival.nc"

LAT_MIN, LAT_MAX = 20.0, 45.0
LON_MIN, LON_MAX = -90.0, -60.0

# Zelinka 2012 / ISCCP Bin Edges
tau_edges = [0.0, 0.3, 1.3, 3.6, 9.4, 23.0, 60.0, 1000.0]
pc_edges  = [10.0, 180.0, 310.0, 440.0, 560.0, 680.0, 800.0, 1000.0]

# Preallocate accumulation grids
isccp_cloudy_histogram = zeros(Float64, 7, 7)
isccp_clear_histogram  = zeros(Float64, 7) 

NCDatasets.Dataset(nc_file, "r") do ds
    println("Reading curvilinear satellite coordinate arrays...")
    lats = coalesce.(ds["latitude"][:,:], NaN)
    lons = coalesce.(ds["longitude"][:,:], NaN)
    
    println("Loading radiometric and cloud variables...")
    reflectance_vis = coalesce.(ds["reflectance_vis"][:,:], NaN)
    alpha_obs       = coalesce.(ds["broadband_shortwave_albedo"][:,:], NaN)
    tau             = coalesce.(ds["cloud_visible_optical_depth"][:,:], NaN)
    cth_km          = coalesce.(ds["cloud_top_height"][:,:], NaN)
    pixel_sza       = coalesce.(ds["pixel_sza"][:,:], NaN)
    alpha_surface   = coalesce.(ds["clearsky_vis_reflectance"][:,:], NaN)
    pixel_vza       = coalesce.(ds["pixel_vza"][:,:], NaN)
    pc              = coalesce.(ds["cloud_top_pressure"][:,:], NaN)
    skin_temp       = coalesce.(ds["pixel_skin_temperature"][:,:], NaN)

    println("Calculating subpixel cloud fractions...")
    cloud_fraction = compute_subpixel_cloud_fraction.(alpha_obs, pixel_sza, tau, cth_km, alpha_surface)

    println("Applying strict geometry constraints to find all physical footprints...")
    # The absolute master filter. It acts as our complete bounding domain, 
    # capturing all successful instrument scans within our target bounds.
    spatial_mask = (lats .>= LAT_MIN) .& (lats .<= LAT_MAX) .& 
                   (lons .>= LON_MIN) .& (lons .<= LON_MAX) .&
                   (pixel_vza .<= 60.0) .& 
                   (pixel_sza .<= 60.0) .&       # Filter solar zenith angle > 60° as missing
                   (.!isnan.(reflectance_vis))   # Ensure data is not a dropped/corrupted scan

    # Unroll coordinates down to 1D vectors for lightning-fast loop ingestion
    v_cloud_fraction = cloud_fraction[spatial_mask]
    v_pc             = pc[spatial_mask]
    v_tau            = tau[spatial_mask]
    v_skin_temp      = skin_temp[spatial_mask]
    
    # Corrected Global Denominator: Captures EVERY valid operational pixel in the grid box
    total_valid_pixels = length(v_cloud_fraction)
    println("Domain Total Footprints (Global Denominator): ", total_valid_pixels)

    if total_valid_pixels == 0
        return
    end

    println("Accumulating fractional mass across the domain track...")
    for i in 1:total_valid_pixels
        f = v_cloud_fraction[i]
        
        # --- CASE 1: Footprint is completely clear (Skipped by cloud retrieval algorithm) ---
        if isnan(f) || f <= 0.0
            clear_fraction = 1.0
            
            # Map clear sky to a surface pressure using its skin temperature
            p_val = estimate_clear_sky_pressure(v_skin_temp[i])
            p_bin = get_bin_index(p_val, pc_edges)
            
            if p_bin > 0
                isccp_clear_histogram[p_bin] += clear_fraction
            end
            
        # --- CASE 2: Footprint contains subpixel cloud structures ---
        else
            clear_fraction  = 1.0 - f
            cloudy_fraction = f
            
            p_bin   = get_bin_index(v_pc[i], pc_edges)
            tau_bin = get_bin_index(v_tau[i], tau_edges)
            
            if p_bin > 0
                # Record the clear portion at this pixel's cloud top height
                isccp_clear_histogram[p_bin] += clear_fraction
                
                # Record the cloudy portion in the 2D joint histogram
                if tau_bin > 0
                    isccp_cloudy_histogram[p_bin, tau_bin] += cloudy_fraction
                end
            end
        end
    end
    
    # Standard Regional Normalization: Divide every element by the universal global denominator
    cloudy_pct = (isccp_cloudy_histogram ./ total_valid_pixels) .* 100.0
    clear_pct  = (isccp_clear_histogram  ./ total_valid_pixels) .* 100.0

    println("\n--- Final 7x7 Cloudy Fraction Histogram (%) ---")
    display(cloudy_pct)
    
    println("\n--- Final 1x7 Clear Fraction Vector by Pressure Bin (%) ---")
    for b in 1:7
        println("Bin $b ($(pc_edges[b]) to $(pc_edges[b+1]) hPa): ", round(clear_pct[b], digits=2), "%")
    end
    
    println("\nVerification: Total Domain Area Sum (Cloudy + Clear) = ", 
            round(sum(cloudy_pct) + sum(clear_pct), digits=1), "%")
end
