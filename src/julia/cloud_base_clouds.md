## Cloud-base clouds

Observations, high-resolution models, single-column models, and climate models represent different cloud fraction at cloud base. Though cloud fraction is maximum at cloud base by construction assuming that all clouds extend from their tops to cloud base, the /cloud top/ height distribution in GOES observations is uniform with no maximum near cloud base. In many representations, however, the maximum cloud fraction is about 5 % at cloud base (Bretherton et al. XXXX, Blossey et al. XXXX, Vogel et al. 2022). These cloud-base clouds contribute strongly to the cloud radiative effect. Knowing the existing fraction of cloud base clouds is a challenging and necessary condition for understanding how they might change.

These lowest of the low clouds are degenerate to our cumulus mass flux model, because they are not represented by coherent cloudy updraft plumes. Rather, cloud-base clouds feed on moisture from turbulent eddies in the subcloud layer. We present subcloud layer turbulence model for diagnosing the cloud-base clouds.

Cloud-base cloud fraction is diagnosed from the fraction of subcloud layer updrafts that reach saturation (Sommeria and Deardorff, 1977; Mellor 1977). The mean saturation specific humidity at cloud base $\overline{q}_s(T_{CB})$ is taken to be constant, and the turbulent subcloud specific humidity $q = \overline{q} + q'$ where the fluctuations $q'$ are normally distributed with width 
$\sigma_q = \overline{q'^2}^{1/2}$. The cloud fraction 
$$
C_{CB} = (1 - \mathrm{erf}(s) ) / 2
$$
is the integral of the humidity distribution greater than saturation, which depends on ratio of the mean saturation depression to the fluctuation width $s = (\overline{q}_{s}(T_{CB}) - \overline{q})/\sigma_q$.

To evaluate the cloud-base cloud fraction, we need to define the width of the fluctuations due to turbulent eddies. We conceptualize this process as thermals that lift moisture from moist lower subcloud layer to the cloud base.
<!--
TABLE. Standard deviation of specific humidity estimated from subcloud humidity scales.

| $\sigma_q$ (g kg$^{-1}$) | source |
| :--- | :--- |
| 0.17 | flux $\overline{w'q'}/\sqrt{w_*^2 + u_*^2}$ |
| 0.19 | flux $\overline{w'q'}/w_*$ |
| 0.46 | flux $q^* = \overline{w'q'}/u_*$ |
| 0.44 | standard deviation below $z_{CB}$ within each sounding |
| 0.99 | standard deviation among soundings |
| 1.5  | vertical $\overline{q_0}-\overline{q_{CB}}$
| 0.53 | turbulent mixed layer scalar variance $\sigma_{q\mathrm{turb}}^2 = (q_0-q_{CB}) \overline{w'q'}/w^*$
[USE 3x THIS ONE: It has the flux and the gradient within the parameterization]
!-->

The mean width of the specific humidity distribution 
$
\sigma_{q\mathrm{turb}} = \overline{q'^2}^{1/2}
$ 
is found from production-dissipation balance of specific humidity variance,
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
\chi_q = K u^* \overline{q'^2} / z_{CB}.
$$
The friction velocity is
$u^* = (\tau/\rho)^{1/2}$.

<!--
$w^* = (z_{CB} B)^{1/3}$
is larger than the friction velocity $u^*$,
where
$B = g \overline{w'T_v'}_0 / T_{v0}$
is the surface buoyancy flux.
!-->
The humidity variance resulting from the production-dissipation balance,
$$
\sigma_{q\mathrm{turb}}^2 = (q_0-q_{CB}) \overline{w'q'}/ (K u^*)
$$
depends on the surface to cloud base humidity difference 
$q_0-q_{CB} = -z_{CB} (\partial \overline{q} / \partial z)$
and the moisture flux 
$\overline{w'q'}$. 
<!--Choosing the empirical constant $K=0.3$, 
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
For constant humidities, the cloud fraction decreses exponentially with on the inverse humidity variance $\sigma_{q\mathrm{turb}}^{-2}.
!-->

### Relative humidity
The cloud-base cloud fraction is expressed as the fraction of parcels whose relative humidity reaches or exceeds 1 at the cloud base. We write the specific humidity standard deviation in terms of the friction velocity $u^*$,
$$
\sigma_q = (k_1\overline{w'q'}(q_0-q_{CB})/u^*)^{1/2}
 = (k_1 u^* q^* (q_0-q_{CB})/u^*)^{1/2}
 = (k_1     q^* (q_0-q_{CB}))^{1/2}.
$$
[the factor $k_1$ absorbs differences between the convective velocity scale and the friction velocity]
We write the standard deviation of relative humidity RH at cloud base,
$$
\sigma_{RH} = (q^* (q_0 - q_{CB}) / q_{sCB}^2)^{1/2}.
$$
The cloud fraction is
$$
C_{CB} = (1-\mathrm{erf}(s)) / 2
$$
with the argument
$
s = (1-\overline{RH}_{CB})/\sigma_{RH},
$
which measures the mean saturation deficit at cloud base, relative to the turbulent fluctuations of RH.

### Scaling subcloud humidity with warming

If wind speed and the vapor transfer coefficient stay the same
and $q_s$ increases by +7 \% /K,
the sensitivity of the surface evaporation flux 
$
E_0 = \rho L C_E U q_s (1-\mathrm{RH})
$
is limited to +2 \% /K by quenching the surface relative saturation deficit
$1-\mathrm{RH}$ by -5 \% /K.
Taking the mean surface relative humidity at which vapor evaporates 
in the tropical circulation to be RH = 0.73, the relative saturation 
deficit 1-RH = 0.27 decreases by -5 \% /K to 1-RH = 0.26, 
an absolute increase of RH of about +0.01 /K.

In a warmer climate, the saturation vapor pressure increases +7 \% /K 
but precpitation increases only +2 \% /K. 
The surface evaporation (and cloud-base moisture flux)
that balances the decrease of precipitation 
also increase +2 \% /K
(resulting in no moisture convergence and steady humidity).
Factoring the surface flux as
$$
\overline{w'q'} = u^* q^*
$$
we expect the stress and friction velocity $u^*$ to remain constant, so the humidity scale $q^*$ increases +2 \% /K, and the absolute humidities scale with +0.07 /K of the saturation humidity.
$$
d\ln\sigma_{RH} = \frac{1}{2}[d\ln q^* + d\ln(q_0-q_{CB})] - d\ln q_{sCB} \\
= \frac{1}{2}[0.02 + 0.07] - 0.07 = -0.025.
$$
We scale the argument of the error function,
$$
d\ln s = d\ln(1-\overline{RH}_{CB}) - d\ln\sigma_{RH} \\
  = -0.05 - (-0.025) = -0.025
$$
EUREC4A observed cloud-base cloud fraction is 
5.4 $\pm$ 3.1 \% 
(from BASTALIAS, Bony et al., 2022; Vogel et al., 2022)
this is near
$$
C_{CB} = (1-\mathrm{erf}(s_0))/2 = 0.055,
$$
which is obtained for $s_0 = 1.13$.
Scaling $s$ by 0.975 /K gives cloud-base cloud fraction
$$
C_{CB}(0.975s_0) = 0.06,
$$
a +9 \% /K increase in cloud-base clouds. 
This 0.5 \% (absolute) /K increase of the cloud-base clouds closely offsets the -0.5 \% /K decrease of the shallow cumulus clouds.

<!--
### IGNORE...
Cloud base occurs near the lifting condensation level (LCL) of the surface air, where $q(0) \approx q_s(CB)$ (TABLE HUMIES). The subcloud humidity gradient decreases with height, $q(CB) < q(0)$, but cloud base relative humidity increases with height $RH(0) < RH(CB) < 1$. Saturation specific humidity increases by +7 \%/K. The cloud base humidity is observed in the present climate. Projecting for constant cloud base relative humidity and increased width $\sigma_q$ due to the flux, the cloud fraction increases from 5.7 \% to 6.9 \% (+19 \%). Projecting RH increases by 1-RH decreasing uniformly by -5 \%, results in cloud-base cloud fraction increasing to 7.3 \% (+28 \%). Constant diffusivity estimated from the present gradient and flux
$
k = -(\partial\overline{q}/\partial z) / \overline{w'q'}
$
projects cloud-base clouds to increase to 8.8 \% /K (+54 \%).

Cloud-base clouds increase in all the experiments.  Increasing cloud base cloud fraction in a warmer climate would contribute a negative effect to the low cloud feedback. The results are implicitly proportional to the cloud-base cloud fraction in the present climate, where we have tuned the parameters of the subcloud turbulence to the EUREC4A cloud-base cloud fraction. They are furthermore sensitve to the parameterization of the mean humidity at cloud base.

TABLE HUMIES. Observed subcloud humidities (g kg $^{-1}$) TOO MUCH PRECISION

| name | value |
| :--- | :--- |
| $q(CB)$   | 13.24  |
| $q_s(CB)$ | 15.12  |
| RH(CB)    | 0.8756 |
| $q(0)$    | 14.74  |
| $q_s(0)$  | 20.22  |
| RH(0)     | 0.7289 |

!-->