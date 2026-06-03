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

## Water diffusion equation for the shallow cumulus region
If we write the flux as a diffusion then the Eulerian moisture equation, which balances in a mean climate,
$$
0 = \frac{\partial q}{\partial t} = 
    \frac{\partial}{\partial z} \left ( \mu \frac{\partial q}{\partial z} \right )
    -W \frac{\partial q}{\partial z}
    -V_H \cdot \nabla q,
$$
can be written as a linear operator of the vertical moisture gradient
(using $d(\mu dq) = d\mu dq + \mu d^2q = (d\mu d + \mu d^2)q $)
$$
0 =
\left ( \frac{\partial \mu}{\partial z  }\frac{\partial}{\partial z} + 
  \mu \frac{\partial^2  }{\partial z^2} 
  - W   \frac{\partial}{\partial z}
  - V_H \cdot \nabla   \right ) q = \mathcal{L_0} q
$$
Note the advection source $-V_H \cdot \nabla q < 0$ dries the cumulus layer. Writing this as a downstream moistening (of the total derivative) on the left hand side.
$$
V_H \cdot \nabla q =
\left ( \frac{\partial \mu}{\partial z}
+ (\mu\frac{\partial}{\partial z} - W) \right ) 
\frac{\partial q}{\partial z} =
\mathcal{L} \frac{\partial q}{\partial z}
$$
The small moistening downstream towards the deep convective region is balanced by the linear operator
$
\mathcal{L} = ( \frac{\partial \mu}{\partial z} + (\mu\frac{\partial}{\partial z} - W)  )
$
on the vertical moisture gradient.
Under climate change, the downstream moistening, 
the linear operator, and the vertical moisture gradient 
$\partial q/\partial z$ all change.

Take $x$ to be the streamwise horizontal direction of the mean flow $U = V_H$ so that the advection is $U\frac{\partial q}{\partial x}$. For this mean circulation 
$$
\frac{\partial U}{\partial x} = -\frac{\partial W}{\partial z}.
$$

Differentiating from one climate to another, we have
$$
0 = d(\mathcal{L}_0 q) = d\mathcal{L}_0 q + \mathcal{L}_0 dq.
$$
where
$$
\mathcal{L}_0 = \frac{\partial \mu}{\partial z  }\frac{\partial}{\partial z} + 
  \mu \frac{\partial^2  }{\partial z^2} 
  - W   \frac{\partial}{\partial z}
  - V_H \cdot \nabla
$$
Differentiating the linear operator,
$$
d\mathcal{L}_0 = d(\frac{\partial \mu}{\partial z  }\frac{\partial}{\partial z}) + 
 d(\mu \frac{\partial^2  }{\partial z^2})
  - d(W   \frac{\partial}{\partial z})
  - d(V_H \cdot \nabla )
= \frac{\partial d\mu}{\partial z  }\frac{\partial}{\partial z}
  +      \frac{\partial \mu}{\partial z  }\frac{\partial}{\partial z}d 
  + d\mu \frac{\partial^2  }{\partial z^2} 
  +  \mu \frac{\partial^2  }{\partial z^2}
  - dW   \frac{\partial}{\partial z}
  -  W   \frac{\partial}{\partial z}
  - dV_H \cdot \nabla
  -  V_H \cdot \nabla
$$
and grouping by operator parameter changes $d\mu$, $dW$, and $dV_H$
apart from humidity changes $dq$,
$$
d\mathcal{L}_0 = 
         \frac{\partial d\mu}{\partial z  }\frac{\partial}{\partial z}
  + d\mu \frac{\partial^2  }{\partial z^2} 
  - dW   \frac{\partial}{\partial z}
  - dV_H \cdot \nabla
  +  \left ( \frac{\partial \mu}{\partial z  }\frac{\partial}{\partial z}
      +  \mu \frac{\partial^2  }{\partial z^2}
      -    W \frac{\partial}{\partial z}
      -  V_H \cdot \nabla \right ).
$$
So
$$
d(\mathcal{L}_0 q) =
  \left (\frac{\partial d\mu}{\partial z  }\frac{\partial}{\partial z}
  + d\mu \frac{\partial^2  }{\partial z^2} 
  - dW   \frac{\partial}{\partial z}
  - dV_H \cdot \nabla \right ) q
  + \left ( \frac{\partial \mu}{\partial z  }\frac{\partial}{\partial z}
      +  \mu \frac{\partial^2  }{\partial z^2}
      -  W \frac{\partial}{\partial z}
      -  V_H \cdot \nabla \right ) dq.
$$

