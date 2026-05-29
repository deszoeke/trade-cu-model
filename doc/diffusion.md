# Diffusion and mass flux

Low cloud response to climate warming in GCMs is often described in terms of the effect of the
turbulence or the cumulus clouds mixing the moisture,
If turbulent and cloud mixing dry the cloud layer, then they reduce the cloud fraction.
Models where low cloud fraction is dominated by clouds near cloud base height
are particularly sensitive to changes in humidity at cloud base height.
Tendencies from the turbulence and cumulus vertical mixing processes 
change the moisture profile in a few days.
On seasonal and climate time scales they establish a climatological moisture profile.
Sherwood et al. (2014), for example, measures cumulus mixing as reduction of vertical gradient of humidity.

In an equilibrium climate, the moisture profile adjusts so that the fluxes of moisture balance.

We parameterize the cumulus moisture flux by a mass flux model. This cumulus moisture flux can be written as a diffusive downgradient flux mixing the vertical moisture gradient:
$$
F_{Cu} = M \Delta q = -\mu (\partial q / \partial z),
$$
where $\Delta q = q_{cld} - q$, and $\mu$ is the diffusivity of moisture by the cumulus flux.

The diffusivity can be estimated directly from $\mu = -F_{Cu}/(\partial q / \partial z)$.

The entraining plume mass flux model is used to estimate the structure of $q_{cld}$:
$$
\partial q_{cld} / \partial z = -(\epsilon + \alpha) (q_{cld}-q),
$$
where lateral entrainment $\epsilon$ and precipitation autoconversion $\alpha$ rates are inverse length scales.
A third inverse length scale $\chi(z)$ is devised such that
$$
\partial(\Delta q)/\partial z = -(\epsilon + \alpha + \chi),
$$
where $\chi$ relates the cloud moisture difference to the vertical gradient,
$$
\chi = -(\partial q / \partial z) / \Delta q > 0 \\
\chi = \partial \ln(\Delta q)/ \partial z + \epsilon + \alpha
$$
The diffusivity can also be written $\mu = M / \chi$.

The diffusivity can be estimated directly from $\mu = -F_{Cu}/(\partial q / \partial z)$.
Ignoring precipitation, this can be estimated from large scale considerations alone, because the all-sky 
average cumulus moisture flux balances the large scale flux $0 = F_{Cu} + F_{LS} - F_0$
where the constant of integration $F_0$ is equal to the surface evaporation. Moisture from the surface evaporation is transmitted (minus a subcloud moistening tendency) by the subcloud layer turbulence to the cloud base.
$$
\mu = F_{LS}/(\partial q / \partial z).
$$
Differentiating between one climate and another,
$$
d\ln \mu + d\ln(\partial q / \partial z) = d\ln F_{LS}.
$$