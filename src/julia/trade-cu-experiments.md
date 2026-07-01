# Intent for Trade Cumulus Experiments

# Cloud categories for different sink rates
 - There are two separate ensembles of cloud categories $i$, one corresponding to `control` (and its experiments) with 600 prior `tot_sink` categories, and one to `control-sink` (and its experiments) with ~331 `sinkz` categories.
 - Both control simulations (C: `control`, `control-sink`) are designed to have cloud top height area fraction distributions $a_i$ that match the observed GOES cloud top height (`cth_bin, rfv_acc, rfv_nrm`).
 - The `control` simulation is prescribed 600 `tot_sink` and simulates `ztop` for each category $i$.
 - The `control-sink` simulation interpolates, from the `tot_sink` and `ztop` of the `control` simulation, the (331) sink rates `sinkz` at the GOES cloud top height bins `cth_bin`, so the cloud fractions `rfv_nrm` can be used directly for cloud fractions in the corresponding cloud category with `ztop=cth_bin`.
 - Since the cloud top height changes in every experiment, it might be just as good to conservatively interpolate the cloud area from the `control`'s `ztop` and `tot_sink` categories (rather than using `sinkz` in the `control-sink`, etc. experiments), but we maintain both ensembles of cloud categories.
 
 # Clouds and moisture fluxes
 - The model computes cloud total water and cloud top height. (Cloud top height was used to get the area fraction for each category above.)
 - The model computes moisture fluxes by the clouds (updrafts and precipitation), by partitioning a prescribed vertical profile of (all-sky mean) fluxes $G(z)$ to the cloud categories $i$ according to their resulting cloud top height.
 - The prescribed flux profile $G$ can be perturbed in experiments. $G$ and the simulated cloud top height affect the flux partitioned to each cloud category $G_i$ in each experiment.

# Experiments
 - Experiments in batch E1 of perturbation experiments are given the same sink rates as the control experiments, and due to other perturbations, new clouds and fluxes are computed for each sink rate category $i$.
 - Experiments in a second batch E2 keep the same cloud categories $i$ as their respective C and E1 experiments, but scale the sink rate in each category.
 - Experiments simulate perturbations to clouds and fluxes due to perturbations prescribed in E1 and E2.
 - Mass fluxes for each category $M_i = G_i / (\Delta q)_i$ in each experiment.

 # Cloud fraction predictions
 - Cloud fraction $a_i$ in the controls C is set to match observations.
 - We use the mass flux $M_i = a_i w_i$ to predict the change to the area fraction $a_i$ in the (non-control) experiments, assuming $w_i$ is invariant from control to experiment, so $d\ln(a_i) = d\ln(M_i)$.



