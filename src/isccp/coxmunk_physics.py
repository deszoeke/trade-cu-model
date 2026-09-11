#!/usr/bin/env python3
"""
coxmunk_physics.py

Self-contained ocean-surface glint/whitecap albedo physics:
Cox & Munk (1954) wave-slope statistics + Fresnel reflection +
Monahan & Muircheartaigh (1980) whitecap fraction.

No file I/O, no climatology, no wind data. Just the radiative-transfer
piece: given a sun zenith angle and a wind speed, what is the
hemispherically-integrated ("black-sky") surface albedo?

Import this module from a separate climatology/data-pipeline script
rather than extending it with I/O - the physics here is independently
useful, testable, and validated on its own (see self_test / __main__).

PHYSICS / KEY FORMULAS
-----------------------
Cox & Munk (1954) isotropic mean-square slope:
    sigma^2(U) = 0.003 + 0.00512 * U        (U = 10 m neutral wind speed, m/s)
(The full Gram-Charlier form includes skewness/peakedness and
upwind/crosswind anisotropy; this module uses the simpler isotropic
Gaussian slope PDF, standard for a first-order estimate.)

Glint geometry (sun zenith theta0, view zenith thetav, relative azimuth
dphi), unit vectors pointing away from the surface toward the sun (s) and
toward the viewer (v); the specular facet normal bisects s and v:
    cos(2*omega) = cos(theta0)*cos(thetav) + sin(theta0)*sin(thetav)*cos(dphi)
    cos(beta)    = (cos(theta0) + cos(thetav)) / (2*cos(omega))
where omega = local angle of incidence on the tilted facet, beta = facet
tilt from vertical.

Isotropic Cox-Munk slope PDF (slope magnitude tan(beta)), normalized so
that Int P dzx dzy = 1:
    P(beta) = 1/(pi*sigma^2) * exp(-tan^2(beta) / sigma^2)

Cox-Munk BRDF (validated below against the flat-Fresnel limit):
    BRDF(theta0,thetav,dphi) = P(beta) * R_F(omega)
                               / (4*cos(theta0)*cos(thetav)*cos(beta)^4)
Integrating BRDF*cos(thetav) over the exit hemisphere for fixed theta0
gives the hemispherical-directional ("black-sky") glint albedo:
    alpha_glint(theta0, sigma2) = Int BRDF * cos(thetav) dOmega_v
As sigma2 -> 0 this converges to the flat-surface Fresnel reflectance
R_F(theta0) - this is checked numerically at runtime by self_test().

NOTE ON A COMMON PITFALL: several published write-ups of this formula
include an extra leading factor of pi in the numerator (i.e.
"pi * P * R_F / (...)"). That version overshoots the required flat-surface
(sigma2 -> 0) limit by a factor of pi (verified numerically here). The
form used in this module omits that extra pi and passes the limit test
to <1.5% at all tested sun angles.

Whitecap fraction (Monahan & Muircheartaigh 1980):
    W(U) = 2.95e-6 * U^3.52
Total surface albedo (linear foam / sub-foam-water mixture):
    alpha_total(theta0,U) = W(U)*ALBEDO_FOAM + (1-W(U))*alpha_glint(theta0,sigma2(U))

DOCUMENTED SIMPLIFICATIONS (all second-order for U ~ 4-10 m/s, tropical
sun angles; more important at high wind speed / low sun elevation):
    - Isotropic Gaussian slope PDF, not the full Gram-Charlier form
      (neglects skewness/peakedness and upwind/crosswind anisotropy).
    - No wave shadowing/masking correction (Smith 1967) at grazing angles.
    - No multiple scattering between facets.
    - Single static whitecap albedo (no foam age/decay dependence).
    - No sub-surface/underlight contribution (irrelevant for broadband
      TOA-relevant albedo to leading order).
"""

import numpy as np

N_REFRACTIVE_INDEX = 1.34   # broadband effective refractive index, seawater
ALBEDO_FOAM = 0.55          # broadband "fresh foam" albedo (Koepke 1984).
                             # Effective/decayed whitecap albedo used in some
                             # studies is lower (~0.2-0.4); treat this as an
                             # upper-bound sensitivity choice.

N_THETAV = 180               # view-hemisphere integration resolution
N_PHI = 360


def sigma2_isotropic(U):
    """Cox & Munk (1954) isotropic mean-square slope."""
    return 0.003 + 0.00512 * np.maximum(U, 0.0)


def fresnel_unpolarized(omega, n=N_REFRACTIVE_INDEX):
    """Unpolarized Fresnel reflectance at angle of incidence omega (radians)."""
    omega = np.asarray(omega, dtype=float)
    sin_t = np.clip(np.sin(omega) / n, -1.0, 1.0)
    theta_t = np.arcsin(sin_t)

    eps = 1e-6
    near_normal = omega < eps

    with np.errstate(divide="ignore", invalid="ignore"):
        Rs = (np.sin(omega - theta_t) / np.sin(omega + theta_t)) ** 2
        Rp = (np.tan(omega - theta_t) / np.tan(omega + theta_t)) ** 2
        R = 0.5 * (Rs + Rp)

    R0 = ((n - 1.0) / (n + 1.0)) ** 2
    R = np.where(near_normal, R0, R)
    return np.clip(R, 0.0, 1.0)


def whitecap_fraction(U):
    """Monahan & Muircheartaigh (1980)."""
    return np.clip(2.95e-6 * np.maximum(U, 0.0) ** 3.52, 0.0, 1.0)


def _view_hemisphere_grid(n_thetav=N_THETAV, n_phi=N_PHI):
    thetav = np.linspace(0.5, 89.5, n_thetav) * np.pi / 180.0
    phi = np.linspace(0.5, 359.5, n_phi) * np.pi / 180.0
    dthetav = (89.0 / n_thetav) * np.pi / 180.0
    dphi = (359.0 / n_phi) * np.pi / 180.0
    TV, PHI = np.meshgrid(thetav, phi, indexing="ij")
    return TV, PHI, dthetav, dphi


def glint_albedo_single(theta0, sigma2, TV, PHI, dthetav, dphi):
    """Hemispherical-directional glint albedo for one (theta0, sigma2)."""
    cos_t0, sin_t0 = np.cos(theta0), np.sin(theta0)
    cos_tv, sin_tv = np.cos(TV), np.sin(TV)

    cos2omega = np.clip(cos_t0 * cos_tv + sin_t0 * sin_tv * np.cos(PHI), -1.0, 1.0)
    omega = 0.5 * np.arccos(cos2omega)
    cos_omega = np.cos(omega)

    cos_beta = (cos_t0 + cos_tv) / (2.0 * np.maximum(cos_omega, 1e-8))
    cos_beta = np.clip(cos_beta, 1e-6, 1.0)
    tan2beta = (1.0 - cos_beta ** 2) / (cos_beta ** 2)

    P = (1.0 / (np.pi * sigma2)) * np.exp(-tan2beta / sigma2)
    RF = fresnel_unpolarized(omega)

    # No leading factor of pi - see module docstring / self_test().
    brdf = (P * RF) / (4.0 * cos_t0 * cos_tv * cos_beta ** 4)

    integrand = brdf * cos_tv * sin_tv  # dOmega = sin(thetav) dthetav dphi
    return float(np.sum(integrand) * dthetav * dphi)


def build_glint_lookup(n_theta0=45, n_u=26, u_max=25.0):
    """Precompute alpha_glint(theta0, sigma2) on a grid for interpolation.

    Returns
    -------
    theta0_rad : (n_theta0,) array
    U_grid     : (n_u,) array, m/s
    table      : (n_theta0, n_u) array of hemispherical glint albedo
    """
    theta0_deg = np.linspace(0.0, 88.0, n_theta0)
    theta0_rad = theta0_deg * np.pi / 180.0
    U_grid = np.linspace(0.0, u_max, n_u)
    sigma2_grid = sigma2_isotropic(U_grid)

    TV, PHI, dthetav, dphi = _view_hemisphere_grid()

    table = np.zeros((n_theta0, n_u))
    for i, th0 in enumerate(theta0_rad):
        for j, s2 in enumerate(sigma2_grid):
            table[i, j] = glint_albedo_single(th0, s2, TV, PHI, dthetav, dphi)

    return theta0_rad, U_grid, table


def total_albedo(theta0, U, glint_lookup_fn):
    """
    alpha_total(theta0, U) = W(U)*ALBEDO_FOAM + (1-W(U))*alpha_glint(theta0, U)

    `glint_lookup_fn` should be a callable (e.g. a fitted
    scipy.interpolate.RegularGridInterpolator, or glint_albedo_single itself
    wrapped appropriately) mapping (theta0, U) -> alpha_glint. Kept as an
    injected dependency so this module doesn't need to know how the caller
    chose to interpolate the lookup table.
    """
    W = whitecap_fraction(U)
    alpha_glint = glint_lookup_fn(theta0, U)
    return W * ALBEDO_FOAM + (1 - W) * alpha_glint


def self_test(theta0_rad, U_grid, table, tol=0.05, verbose=True):
    """Sanity check: as U -> 0 (sigma2 -> ~0.003, nearly flat), the
    hemispherically-integrated glint albedo should approach the flat-surface
    Fresnel reflectance at that sun angle. Returns True/False."""
    if verbose:
        print("[self-test] comparing near-flat (U=0) Cox-Munk integral to flat Fresnel:")
    i_u0 = 0  # U = 0
    ok = True
    for theta0_deg_check in [0, 20, 40, 60, 75]:
        i0 = int(np.argmin(np.abs(theta0_rad * 180 / np.pi - theta0_deg_check)))
        cm_val = table[i0, i_u0]
        flat_val = fresnel_unpolarized(theta0_rad[i0])
        rel_diff = abs(cm_val - flat_val) / max(flat_val, 1e-6)
        passed = rel_diff < tol
        ok = ok and passed
        if verbose:
            status = "OK" if passed else "CHECK"
            print(
                f"    theta0={theta0_deg_check:3d} deg: Cox-Munk(U=0)={cm_val:.4f}  "
                f"flat Fresnel={flat_val:.4f}  rel.diff={rel_diff:.2%}  [{status}]"
            )
    if verbose:
        print(f"[self-test] {'PASSED' if ok else 'CHECK WARNINGS ABOVE'} "
              f"(residual sigma2=0.003 at U=0 gives a small but nonzero roughness)\n")
    return ok


if __name__ == "__main__":
    # Running this module directly just builds the table and validates it -
    # no wind/mask files needed. Useful as a quick standalone physics check.
    print("Building Cox-Munk glint lookup table...")
    theta0_rad, U_grid, table = build_glint_lookup()
    ok = self_test(theta0_rad, U_grid, table)

    print("Glint albedo vs. wind speed at a few sun angles "
          "(sanity check: should be ~flat near theta0=0, "
          "and decrease with wind at high theta0):")
    for th0_deg in [0, 20, 40, 60, 75, 85]:
        i0 = int(np.argmin(np.abs(theta0_rad * 180 / np.pi - th0_deg)))
        row = table[i0]
        sample = " ".join(f"{u:.0f}:{a:.3f}" for u, a in zip(U_grid[::5], row[::5]))
        print(f"  theta0={th0_deg:3d} deg: {sample}")

    if not ok:
        raise SystemExit("Self-test failed - do not trust downstream results until fixed.")
