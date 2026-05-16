# On the factor for mean vertical velocity within dropsonde circles

Net (cloud and subsidence) mass flux $M_{net}$ goes to zero for a large domain, especially approaching infinite size.
Dropsonde circles have finite radius $r$=110 km. 
The net mass flux for a finite area is
$$
M_{net} = M_{up} - M_{dn},
$$
with $M_{up}>=0$ in the cloud and $M_{dn} >= 0$ in the compensating subsiding area.
If the finite area contains much cloud then $M_{net}>0$. $M_{net}$ is measured over the whole area (fraction unity) so $M_{net} = W$ is just the finite area-mean velocity.

Let us parameterize $M_{net} = \gamma M_{up}$. Then
$$
\gamma M_{up} = M_{up} - M_{dn} \\
\gamma = 1 - M_{dn}/M_{up} \\
M_{dn} = (1-\gamma) M_{up}.
$$
Defining
$\alpha = W a_{up}/M_{up}$, we have
$ \alpha = a_{up} \gamma$.

We found for $r$=110 km, and negligibly weak mass flux $M_{dn} = w_{dn}(1-a_{up})$ suitable for balancing over very large area, that 
$ W_\mathrm{JOANNE} = \alpha W_\mathrm{GOES} = \alpha M_{up}/a_{up}$ with $\alpha \approx$ 0.1.


Limiting behaviors:

$\gamma = 1; M_{dn}=0$

$\gamma = 0; M_{dn}=M_{up}$

When $a_{up}=1$ then $\alpha=\gamma=1$, $M_{dn}=0$, and $w_{dn}$ is undefined.

A critical point is at $\alpha = a_{up}$. This is not impossible since both are near 0.1. In this case $\gamma = 1$, $M_{dn} = 0$, $W = M_{net} = M_{up}$.

I think that $\gamma$ and $\alpha$ should vary with the size $r$ of the averaging domain. The smaller the averaging area, the stronger the compensating subsidence. This effect also reduces the subsidence away from the clouds. 

If the linear scaling factor is a good model, then it's as if all the clouds are fluxing the moisture up with less total mass flux. This might represent some small inner-cloud process. For example, if turbulent updrafts and downdrafts transport moisture against the mean gradient, then these carry some of the moisture flux, reducing the need for bulk mass flux in the plume. It would be an interesting story if only ~1/10 of the plume mass flux were sufficient to make fluxes to balance the moisture budget of the free troposphere.