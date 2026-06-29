# Context Plan - 2026-06-28

## Purpose
Pause handoff for supervised code changes tomorrow.
This captures what was learned, what was reverted, and what is still unresolved.

## Current State (as of pause)
- File state is mixed after iterative edits and user undo actions.
- `src/julia/TradeCuModel.jl` currently has:
  - `compute_normalized_fluxes` masking above cloud top in normalized flux outputs.
  - `cloudflux_allsky(..., allskyeddyflux, ...)` currently partitions by `Feddy_i = allskyeddyflux / a_i`.
- `src/julia/TradeCuExperiments.jl` currently uses:
  - `F2z, G_ls = calc_Ftot(...)`
  - `cloudflux_allsky(..., allskyeddyflux=G_ls, ...)`
  - `M = F_cld ./ (qcld - qm)` and `w = copy(M)` (in-cloud convention).
- User reverted at least one previous assistant change in `TradeCuExperiments.jl`.

## Agreed Physical Workflow (from user)
1. Run cloud model ensemble first to get cloud tops `h_i` (or `ztop_i`) from solved `q_{c,i}`.
2. Project cloud fractions onto solved `h_i`.
3. Compute category fluxes and decomposition in `cloudflux_allsky`:
   - total category moisture flux `F_i`
   - updraft flux `F_cld,i`
   - precipitation flux `F_p,i`
4. Sometimes perform a first cloud-ensemble run, then interpolate resulting `h_i` to the 10 m grid, and use that as input for a refined control run.

## Key Inconsistencies Identified
1. Flux partition target mismatch
- Current code path in `integrate_experiment!` passes `G_ls` into `cloudflux_allsky` and divides by `a_i` there.
- But partitioning should use category flux by cloud-top category (derived from `calcF2` / `F2z`) and then map/evaluate with solved `h_i`.

2. Mass-flux convention mismatch
- User convention: stored matrix `M` should mean all-sky mass flux by category and height:
  - `M_i = a_i * w_i`
  - `M_tot = sum_i(a_i * w_i)`
- Current code stores `M = F_cld / Delta q` and `w = M`, which corresponds to in-cloud mass flux (not all-sky `M_i`).

3. Symbol naming ambiguity
- In code, `F_cld` is treated as in-cloud updraft moisture flux.
- User notation can read `Fcld_i = M_i Delta q_i`; if `M_i` is all-sky, then this would imply an all-sky flux symbol.
- Need explicit naming convention for in-cloud vs all-sky variables.

## What Was Confirmed Previously
- Two real bug mechanisms had been found earlier in the iteration:
  1. Extra `1/a` amplification path in velocity calculation (historically produced order `10^6` m/s spikes).
  2. Missing/insufficient masking above cloud top in normalized flux path.
- Current high values are lower than before but still physically too large, indicating remaining arithmetic/partition consistency issues.

## Supervised Change Plan for Tomorrow
1. Lock conventions first (required before code edits)
- Decide and document:
  - `w_i`: in-cloud plume velocity.
  - `M_i`: all-sky category mass flux (`a_i * w_i`) or in-cloud mass flux.
  - `F_cld,i` and `G_cld,i` naming split (in-cloud vs all-sky).

2. Rewire flux partition in one place
- In `integrate_experiment!` + `cloudflux_allsky`:
  - Keep cloud-top-first workflow.
  - Ensure `calc_Ftot` provides large-scale all-sky reference and/or category partition driver.
  - Use category-consistent `F_i` for decomposition after solved `h_i` are available.

3. Implement mass-flux outputs per chosen convention
- If user convention is kept:
  - `w = F_cld / Delta q` (in-cloud)
  - `M = a_i * w` (all-sky stored matrix)
  - Keep all-sky moisture outputs as `G_cld = a_i * F_cld`, `G_pcp = a_i * F_pcp`.

4. Add explicit diagnostics in code path
- Diagnostic checks during development:
  - `max|G - sum_i(a_i * F_i)|`
  - `max|G_cld + G_pcp - reconstructed total|`
  - per-category constancy of `F_i` from cloud base to `h_i`.

5. Validate with control and one DIM-style run
- Control: baseline `cth_*`.
- DIM/mesoscale: altered `cth_*`.
- Verify conventions and magnitudes before broader runs.

## Minimal Test Snippets to Re-run Tomorrow
- Recompute control and print:
  - max finite `w`
  - location `(z, i)` of max
  - local `F_cld`, `Delta q`, `a_i`
- Reconstruct all-sky moisture flux and compare against large-scale target profile.

## Notes
- User requested supervised edits only after pause.
- Do not auto-apply additional arithmetic changes without confirming conventions first.
