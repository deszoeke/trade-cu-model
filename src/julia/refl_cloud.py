import numpy as np

def compute_subpixel_cloud_fraction(alpha_obs, sza_deg, tau, cth_km, alpha_surface=0.05):
    """
    Computes subpixel cloud fraction using a two-stream radiative transfer inversion.
    
    Parameters:
    - alpha_obs: Observed pixel albedo from GOES
    - sza_deg: Solar Zenith Angle in degrees
    - tau: Cloud Optical Depth reported for the pixel
    - cth_km: Cloud Top Height in km (used to determine cloud phase)
    - alpha_surface: Clear-sky background albedo (default 0.05 for ocean)
    """
    # 1. Handle clear or edge cases
    if tau <= 0 or alpha_obs <= alpha_surface:
        return 0.0
    
    # 2. Geometry: Cosine of Solar Zenith Angle
    mu_0 = np.cos(np.radians(sza_deg))
    if mu_0 <= 0:  # Nighttime/Terminator safety check
        return 0.0

    # 3. Microphysics: Determine asymmetry parameter based on Cloud Top Height
    # High clouds (>6 km) are treated as ice; low clouds as liquid water
    g = 0.75 if cth_km > 6.0 else 0.85
    
    # 4. Two-Stream Parameterization
    gamma_1 = 0.75 * (1.0 - (g * mu_0))
    
    # 5. Calculate theoretical 100% overcast cloud albedo
    alpha_cloud = (gamma_1 * tau) / (1.0 + (gamma_1 * tau))
    
    # 6. Apply Linear Mixing Model to extract cloud fraction
    if alpha_cloud <= alpha_surface:
        return 0.0
        
    f = (alpha_obs - alpha_surface) / (alpha_cloud - alpha_surface)
    
    # Bound the results physically between 0 and 1
    return np.clip(f, 0.0, 1.0)

# Example usage:
# Pixel observed albedo = 0.35, SZA = 30°, Optical Depth = 4.0, Cloud Top Height = 3.5 km
f_proxy = compute_subpixel_cloud_fraction(alpha_obs=0.35, sza_deg=30.0, tau=4.0, cth_km=3.5)
print(f"Estimated subpixel cloud fraction proxy: {f_proxy:.2f}")
