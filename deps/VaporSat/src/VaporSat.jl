# module Vaporsat for qsat for basic cloud calculations
module VaporSat

using ForwardDiff

export qsat, dqsdT, dlnqsdT
export Twet_autodiff

# constants
Cp = 1005.7  # from my Davies-Jones function, was 1005.
Cpv = 1870.0 # J/kg/K
Cw  = 4190.0
L0 = 2.501e6 # J/kg

C = 273.15 # K
Rd = 287.04
Rv = 461.5
RdoRv=Rd/Rv

grav = 9.8 # m/s2

"latent heat of water vapor"
LvK(TempK) = L0 + (Cpv-Cw) * (TempK-273.0)

# functions
"""
es(T,p) = is saturation vapor pressure based on Wexler's formula,
with enhancement factor for moist air rather than water vapor.
The enhancement factor requires a pressure.
T [degrees C], p [Pa] (note the reversed input order), es [Pa];
keyword TK [Kelvin] overwrites positional T[C].
es(T,p; TK=tk[Kelvin], P=pr[hPa])
From A. L. Buck 1981: JAM, 20, 1527-1532.
SPdeS 7 July 2004
"""
#function es(T,p=1e5; P=p*1e-2)
#    esat = 1e2 * 6.1121*(1.0007 + 3.46e-8*P)*exp((17.502*T)/(240.97 + T)) # convert es to Pa
#end
"T [C]; supplying TK [Kelvin] by keyword ignores positional T"
function es(T,p=1e5; TK=T+C, P=p*1e-2)
    T = TK - C
    esat = 1e2 * 6.1121*(1.0007 + 3.46e-8*P)*exp((17.502*T)/(240.97 + T)) # converts es to Pa
end

"""
qsat(p,T) is saturation specific humidity based on Wexler's formula for es
with enhancement factor (see es.m).
p [Pa], T [degrees C], qsat [kg/kg]
From A. L. Buck 1981: JAM, 20, 1527-1532.
SPdeS 7 July 2004
"""
qsat(p,T) = RdoRv*es(T,p) / (p + (RdoRv-1)*es(T,p))

"dqsdT(p,T[C]) derivative of qsat with respect to T at p,T by autodiff of Bolton's qsat"
dqsdT(p,T) = ForwardDiff.derivative(t -> qsat(p,t), T)
dlnqsdT(p,T) = ForwardDiff.derivative(t -> log(qsat(p,t)), T)

# wet bulb temperature methods
# for approximating the evap process

"General single Newton iteration to update x toward f(x) = fhat for a univariate function f"
updatex(f, x, fhat) = x + (fhat-f(x)) / ForwardDiff.derivative(f, x)

"""
Twet_autodiff(T[K], q[kg/kg], p[Pa]; niter=2) wet bulb temperature using Newton's method
for target specific humidity q[kg/kg]. Uses automatic differntiation.
"""
function Twet_autodiff(T, q, p; niter=2)
    b(t) = f(t) = (t - T) + LvK((T+t)/2)/Cp * (qs(p,t-C) - q)
    t=T
    for i in 1:niter
        t = updatex(f, t, 0)
    end
    t
end
# 2 iterations converges to ~0.001 K

# call as...
# q = rh*qsat(pa, Ta)
# Twet_autodiff(Ta, rh*qsat(pa, Ta-C), pa) 

#=
# dew point temperature
function dewpoint_from_q_buck1981(q, P=101325.0)
    # Convert specific humidity to vapor pressure
    e = (q * P) / (0.622 + 0.378 * q)
    
    # Invert Buck (1981): es = 611.21 * exp((18.678 - T/234.5) * (T / (257.14 + T)))    
    b = log(e / 611.21)
    Td = (257.14 * b) / (18.678 - b)
end
=#

# integrate pressure from the hydrostatic equation
dlnp(Tv, dz) = -gravity * dz / (Rd * Tv)
press_hydrostatic(Tv, dz, p0) = p0 * exp(cumsum(dlnp(Tv, dz)))

# reference pressure [Pa] profile
p_q(Tm) = interpolate((Tm,), press_hydrostatic.((1+(1/RdoRv-1)*qm)*Tm, dz, p0))

"""Tcond_autodiff(T[K], q[kg/kg], p[Pa]; niter=2) 
Temperature at which condensation likely occurred 
for target specific humidity q[kg/kg] and reference pressure p[Pa].
Uses automatic differentiation.
"""
function Tcond_autodiff(q; niter=5, tol=1e-2)
    f(t) = p_q(t)/esat(t-C) - (RdoRv*(1/q - 1) + 1)
    t = 280.0 # first guess
    for i in 1:niter
        tmp = updatex(f, t, 0)
        if abs(tmp - t) < tol
            t = tmp
            break
        end
        t = tmp
    end
    t
end

end # module VaporSat
