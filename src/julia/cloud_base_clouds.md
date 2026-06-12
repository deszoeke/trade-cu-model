## Cloud-base clouds

Observations, high-resolution models, single-column models, and climate models represent different cloud fraction at cloud base. Though cloud fraction is maximum at cloud base by construction assuming that all clouds extend from their tops to cloud base, the /cloud top/ height distribution in GOES observations is uniform with no maximum near cloud base. In many representations, however, the maximum cloud fraction is about 5 % at cloud base (Bretherton et al. XXXX, Blossey et al. XXXX, Vogel et al. 2022). These cloud-base clouds contribute strongly to the cloud radiative effect, and knowing the existing fraction of cloud base clouds is a challenging and necessary condition for understanding how they might change.

These lowest of the low clouds are degenerate to our cumulus mass flux model. Rather than representing moisture flux in coherent cloudy updraft plumes, cloud-base clouds feed on moisture from turbulent eddies in the subcloud layer. We present an alternative subcloud layer turbulence model for diagnosing the cloud-base clouds in this process.

Cloud-base cloud fraction is diagnosed from the fraction of subcloud layer updrafts that reach saturation (Sommeria and Deardorff, 1977; Mellor 1977). The mean saturation specific humidity at cloud base $\overline{q}_s(T_{CB})$ is taken to be constant, and the turbulent subcloud specific humidity $q = \overline{q} + q'$ where the fluctuations $q'$ are normally distributed with width 
$\sigma_q = \overline{q'^2}^{1/2}$. The cloud fraction 
$$
C_{CB} = (1 - \mathrm{erf}(s) ) / 2
$$
is the integral of the humidity distribution greater than saturation, which depends on ratio of the mean saturation depression to the fluctuation width $s = (\overline{q}_{s}(T_{CB}) - \overline{q})/\sigma_q$.

To evaluate the cloud-base fraction, we need to define the width of the fluctuations due to turbulent eddies. We conceptualize this as thermals that lift moisture from moist lower subcloud layer to the cloud base. The width is thus simply written in terms of the surface - cloud base specific humidity difference $\sigma_q = q_0 - q_{CB}$. Using this definition for the /Ron Brown/ soundings, $\sigma_q = $ 1.5 g kg$^{-1}$ and the cloud-base cloud fraction is 4 % (absolute), which agrees broadly with EUREC4A cloud base cloud fraction observations.

TABLE. Standard deviation of specific humidity estimated from subcloud humidity scales.

| $\sigma_q$ (g kg$^{-1}$) | source |
| :--- | :--- |
| 0.17 | flux $\overline{w'q'}/\sqrt{w_*^2 + u_*^2}$ |
| 0.19 | flux $\overline{w'q'}/w_*$ |
| 0.46 | flux $\overline{w'q'}/u_*$ |
| 0.44 | standard deviation below $z_{CB}$ within each sounding |
| 0.99 | standard deviation among soundings |
| 1.5  | vertical $\overline{q_0}-\overline{q_{CB}}$
| 0.53 | turbulent mixed layer scalar variance $\sigma_{q\mathrm{turb}}^2 = (q_0-q_{CB}) \overline{w'q'}/w^*$
[USE 3x THIS ONE: It has the flux and the gradient within the parameterization]


The mean width of the specific humidity distribution $\sigma_{q\mathrm{turb}} = $\overline{q'^2}^{1/2}$ is parameterized from a production-dissipation balance of specific humidity variance,
$$
0 = -\overline{w'q'}(\partial \overline{q}/\partial z) - \chi_q.
$$
where 
$
\chi_q = \kappa \overline{q'(\partial^2q'/\partial x_i^2)}
$
is the humidity variance dissipation
for the subcloud layer from the surface to cloud base $z_{CB}$ 
with uniform humidity gradient
$
\partial \overline{q}/\partial z.
$
The dissipation (dimensions $[ \chi_q ] = [ q'^2 ]/\mathrm{time}$) is parameterized
$$
\chi_q = K w^* \overline{q'^2} / z_{CB}.
$$
The convetive velocity scale
$w^* = (z_{CB} B)^{1/3}$
is larger than the friction velocity $u^*$,
where
$B = g \overline{w'T_v'}_0 / T_{v0}$
is the surface buoyancy flux.
The variance resulting from the production-dissipation balance,
$$
\sigma_{q\mathrm{turb}}^2 = (q_0-q_{CB}) \overline{w'q'}/ (K w^*)
$$
depends on the surface to cloud base humidity difference 
$q_0-q_{CB} = -z_{CB} (\partial \overline{q} / \partial z)$
and the moisture flux 
$\overline{w'q'}$. 
Choosing the empirical constant $K=0.1$, 
we get 
$\sigma_{q\mathrm{turb}} = $ 1.7 g kg $^{-1}$
and cloud-base cloud fraction 
$$
C_{CB} = \{ 1 - \mathrm{erf}(
  [q_s(z_{CB})-q(z_{CB})] / 
  [\sigma_{q\mathrm{turb}}]) \} / 2
  = 5.7 \%
$$
matching the observed EUREC4A cloud-base cloud fraction 
(5.4 $\pm$ 3.1 \% from BASTALIAS, Bony et al., 2022; Vogel et al., 2022).

### Scaling for subcloud humidity under warming

In a warmer climate, the saturation vapor pressure increases +7 \% 
but precpitation increases only +2 \%. 
The equal cloud-base moisture flux and surface evaporation 
that balance the decrease of precipitation 
also increase +2 \%
(resulting in no moisture convergence and steady humidity).

### Surface flux

If wind and the heat coefficient stay the same