# cloud optical properties

g_asym = 0.85 # for liquid phase, 0.75 for ice phase

refl_cloud(tau; g_asym=g_asym) = (1-g_asym)*tau / (2 + (1-g_asym)*tau)
mu0(angle) = cosd(angle)


function albedo_cloud(tau, mu0; g_asym=g_asym)
   gamma1 = 0.75 * (1 - g_asym * mu0)
   return gamma1 * tau / (1 + gamma1 * tau)
end