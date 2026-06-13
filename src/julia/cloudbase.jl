# load modules
using Revise
using Pkg
# Pkg.activate(joinpath(homedir(), "Projects/ATOMIC/ATOMIC_GOES/julia/ATOMIC_GOES"))
cd("/Users/deszoeks/Projects/ATOMIC/trade-cu-model/src/julia")
Pkg.activate(".")

using Printf, Dates, NCDatasets, JLD2
using Statistics, SpecialFunctions, Interpolations, LinearAlgebra

# using PyPlot -> not thread safe in IJulia with Julia 1.12.4+
using PythonPlot
using VaporSat # dev ../../deps/VaporSat
includet("TradeCuModel.jl")
using .TradeCuModel

using NaNStatistics
using Missings

# parameters and initialization
KelvinCelsius=273.15 # K
zi   = 4.0e3    # m
ztop = 4.0e3
zcb  = 700      # m
divg = 1.5e-6   # 1/s
x = 0.53 # parameter precipitation efficiency

# ensemble of sink rates
ns = 600 # number of sink rates
tot_sink = range(6.3523e-4, 5.7e-3, length=ns) # min tuned for x=0.53 to get the highest possible cloud top
# tot_sink = (1 .+tanh.(range(-8*pi, 0, length=ns))) .* (5e-3 - 1e-4) .+ 1e-4
# tot_sink = range(6.1716e-4, 5.8e-3, length=ns) # min tuned for x=0.53 to get the highest possible cloud top

# get sounding data
z, tam, thm, qm, pm = get_mean_soundings() # missing in, NaN out
#m  K       kg/kg Pa 
qs  = qsat.(pm, tam.-KelvinCelsius) # kg/kg
tvm  = virtual_temp.(tam, qm) # virtual temperature, K
thvm = virtual_temp.(thm, qm) # virtual potential temperature, K
dz = z[2]-z[1]
icb = findfirst(z .>= zcb) # cloud base index

# standard deviation among soundings
ds = get_sounding_dataset()
qdata = ds[:q][:,:]
qstd = tstd(qdata)
qstdsubcloud = Vector{eltype(qdata)}(undef, size(qdata,2))
for i in 1:size(qdata,2)
    qstdsubcloud[i] = std(skipmissing(qdata[1:icb,i]))
end
# qstdsubcloud = 0.44 g/kg -- much smaller than the mean subcloud q difference of 1.5 g/kg
# subcloud vertical difference of Dq 1.5 g/kg, but standard deviation is 0.44 g/kg.
# while std among soundings is 0.99 g/kg.
# get An-Yi's aggregated GOES cloud fraction vs height data
rfv_nrm, rfv_acc, cth_bin = NCDataset("../../data/goes16_binned_low4km_20200115_20200219.nc") do dsa
    rfv_nrm = mean(dsa[:rfv_nrm][:,:], dims=2)
    rfv_acc = mean(dsa[:rfv_acc][:,:], dims=2)
    cth_bin = dsa[:cth_bin][:]
    return rfv_nrm, rfv_acc, cth_bin
end
#size(rfv_nrm), size(cth_bin) # (351,1)
# line up model grid and GOES cloud top height bins
offset = findfirst(x->x≈cth_bin[1]*1e3, z) - 1 # 50; cth_bin starts at z=500

# get PSL data for turbulence parameters
psl = NCDataset("../../data/flux/EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc")
ustar = mean(skipmissing(psl[:ustar][:])) # 0.3 m/s
# wstar = mean(skipmissing(psl[:wstar][:]))
tair = mean(skipmissing(psl[:tair][:])) # 300 K
qair = mean(skipmissing(psl[:qair][:]))
hs = mean(skipmissing(-psl[:hs_bulk][:])) # 300 K
hl = mean(skipmissing(-psl[:hl_bulk][:]))

Rd = 287.0 # J/kg/K
invrho(t=tair+273.15,q=qair*1e-3 ,p=101800) = Rd*virtual_temp(t, q) / p
"compute sfc buoyancy flux in m^2 s^-3 from PSL met data"
function calc_bflx(t=tair+273.15, q=qair*1e-3, hs=hs, hl=hl)
    Cp = 1005
    g = 9.8
    inv_rho = invrho(t, q)
    wt = inv_rho * hs / Cp
    wq = inv_rho * hl / LvK(t)
    fac = 1 + 0.608*q
    fac2 = 0.608 / fac
    g * ( wt/t + fac2*wq ) # buoyancy flux, m^2 s^-3
end
# WARNING PSL hb_bulk is B * (ρCpT)/g [W/m^2]
# psl[:hb_bulk][:]/3.5e4

# surface buoyancy flux
bflx0 = calc_bflx()
wstar = (700.0 * bflx0)^(1/3) # 0.72 m/s
Uscl = sqrt(ustar^2 + wstar^2) # 0.77 m/s

# cloud base water flux from vapor and precipitation flux
# rhb_prate = mean(skipmissing(psl["prate"][:])) / 3600
rhb_prate = 8.88e-6 # kg s-1  # L*rhb_prate = 27.5 W/m^2
Ecb = 153 # W/m^2; 180.6 - 27.5; just the turbulent vapor flux at cloud base
rhoL = mean(filter(isfinite, @. pm / (Rd*tvm) * LvK(tam))) # 2.41e6 J/m^3
wq = Ecb/rhoL # not used when Dq is used directly !!
Dq = qm[isfinite.(qm)][1] - qm[icb] # 0.0015 kg/kg
# standard deviation of induced by flux 
sigmaq_flux = wq ./ [ustar, wstar, Uscl] # 0.46, 0.19, 0.17 g/kg

# cloud base condensation fraction
Ccb = 0.5 * (1 - erf((qs[icb]-qm[icb])/Dq))

function cloudbase_cloudfrac(; qm=qm, qs=qs, sigmaq=qm[isfinite.(qm)][1] - qm[icb], icb=icb)
    Ccb = 0.5 * (1 - erf((qs[icb]-qm[icb])/sigmaq))
end

cloudbase_cloudfrac(icb=icb  ) # 0.038 at 700 m
cloudbase_cloudfrac(icb=icb+5) # 0.065 at 750 m

[ cloudbase_cloudfrac(icb=i) for i in icb.+(-10:20) ] # 0.11 at 700 m
minimum( (qs .- qm) |> filter(isfinite) ) # 0.0015 kg/kg at 700 m

function safe_argmin(x)
    idx = findall(isfinite, x)
    idx[argmin(x[idx])]
end
safe_argmin((qs .- qm)[1:500])

# look at humidity profiles in SCL
Dq = qm[isfinite.(qm)][1] - qm[icb] # 1.5 g/kg mean diff across subcloud
sigmaq = nanmean(qstd[1:icb]) # 0.99 g/kg for subcloud
plot(1e3*qstd[1:200], 1e-3*z[1:200])
plot(1e3*qm[1:200], 1e-3*z[1:200])
plot(1e3*qs[1:200], 1e-3*z[1:200])

display(gcf())

# From production-disspation balance of the q variance
# sigmaq^2 = wq * Dq / wstar
# strictly this is a proportionality, not an equality, but it gives the right order of magnitude for sigmaq in the SCL
fac = 0.1:0.01:0.5
Ccb_fac = zeros(length(fac))
for (i,f) in enumerate(fac)
    sigmaq_MLturb = sqrt(wq * Dq / (f * wstar) ) # .53 g/kg, close to the observed 0.44 g/kg
    Ccb_fac[i] = 0.5 * (1 - erf((qs[icb]-qm[icb])/(sigmaq_MLturb))) # 5.7 %
end
plot(fac, Ccb_fac)

# overload the method
"""
cloudbase_cloudfrac(; Dq, satdef_cb, wq, wstar=wstar)
"""
function cloudbase_cloudfrac(; Dq, satdef_cb, wq, wstar=wstar, fac=0.1)
    sigmaq = sqrt(wq * Dq / (fac * wstar) ) # g/kg
    Ccb = 0.5 * (1 - erf((satdef_cb)/sigmaq)) # cloud base cloud fraction
end

q0 = filter(isfinite, qm)[1]
qs0 = filter(isfinite, qs)[1]
qcb = qm[icb]
qs_cb = qs[icb]
RH_0 = q0 / qs0 # 0.73
RH_cb = qcb / qs_cb # 0.8756

# CONTROL, present climate
Ccb_ctl = cloudbase_cloudfrac(Dq = q0-qcb, 
                    satdef_cb = qs_cb-qcb,
                    wq = Ecb/rhoL) # 0.057
# increase Ecb +2 % to balance precip increase of +2 %.
Ccb_fluxonly = cloudbase_cloudfrac(Dq = q0-qcb,
                    satdef_cb = qs_cb-qcb,
                    wq = 1.02 * Ecb/rhoL) # 0.05903, +3.2%
# 0.05903 : just increasing Ecb +2 % increases cloud base cloud fraction to 5.9 %,
# a relative change of +3.1 %.

# adjust the q, qs profile for climate change, too.
# they are solved in terms of the original qs for RH0=0.8, RH1=0.81

cloudbase_cloudfrac(Dq = 1.07 * (q0-qcb),
                    satdef_cb = 1.07 * (qs_cb-qcb),
                    wq = 1.02 * Ecb/rhoL) # 0.0530

# try surface saturation at LCL condition
cloudbase_cloudfrac(Dq = 1.07 * (q0-qcb), # supposing uniform RH1
                    satdef_cb = 0.95*1.07 * (qs_cb-qcb),
                    wq = 1.02 * Ecb/rhoL) # 0.0623

# Project the new q_CB.
# Several different alternative assumptions can be made about how it changes with climate change:
# In all of them, (1-RH) quenches by *0.95 at the surface to reduce the evaporative vapor flux (to *1.02).
# 0. RH_cb is constant at cloud base so satdef_cb = (1-RH_cb) * qs_cb
Ccb_constRHcb = let qs_cb=1.07*qs_cb, qcb=RH_cb*qs_cb, qs0=1.07*qs0, q0=(1-0.95(1-RH_0))*qs0
    cloudbase_cloudfrac(Dq = q0-qcb,
                    satdef_cb = (1-RH_cb)*qs_cb,
                    wq = 1.02 * Ecb/rhoL) # 0.06854, +19.8 %
end
# 1. RH increases uniformly by 1 % (absolute) on [0 zcb]. 
# DO NOT USE: the ratio 0.95*(1-RH) is better!
#   q_cb increases by 1.07*(RH0+1%).
#   qs_cb - q_cb = [1-(RH0+1%)] * qs_cb = [1-(RH0+1%)] * 1.07*qs_cb0
Ccb_RHp1pct = let qs_cb=1.07*qs_cb, qcb=(RH_cb+0.01)*qs_cb, q0=1.07*(RH_0+0.01)*qs0
    cloudbase_cloudfrac(Dq = q0-qcb,
                    satdef_cb = qs_cb-qcb,
                    wq = 1.02 * Ecb/rhoL) # 0.07186, +25.5 %
end
# 2. (1-RH) decreases by the ratio *0.95 at all z, so q_CB increases by 0.95*1.07.
#   qs_cb - q_cb = 0.95 * (1-RH0) * qs_cb = 0.95 * 1.07 * (1-RH0) * qs_cb0
Ccb_RHp95 = let qs_cb=1.07*qs_cb, qcb=(1-0.95*(1-RH_cb))*qs_cb, qs0=1.07*qs0, q0=(1-0.95*(1-RH_0))*qs0
    cloudbase_cloudfrac(Dq = q0-qcb,
                    satdef_cb = 0.95*(1-RH_cb)*qs_cb,
                    wq = 1.02 * Ecb/rhoL) # 0.07334, +28.1 %
end

# 3. Diffusivity k is constant and uniform
#  zcb will appear in \sigma_q instead of the vertical q difference now.
#  wq = -k dq/dz = k (q_sfc - q_cb) / zcb
#  dq/dz = - 1/k * wq
#  q0 - q_cb = zcb/k * wq
#  sigma_q^2 = wq * (q_sfc - q_cb) / (0.1 * wstar) = wq^2 * zcb/(0.1*k*wstar)
Ccb_k = let k = Ecb/rhoL * zcb / (-(qcb-q0)), wq = 1.02 * wq, q0 = (1-0.95*(1-RH_0)) * 1.07*qs0
    q_cb = q0 - wq * zcb/k 
    satdef_cb = 1.07*qs_cb - q_cb
    sigma_q = wq * sqrt(zcb/(0.1*k*wstar))
    ccf = 0.5 * (1 - erf(satdef_cb/sigma_q)) # 0.08808, +54 %
end

# experiment with cloud base RH
fac = 0.9:0.01:1.1
Ccb_varRHcb = zeros(length(fac))
let qs_cb=1.07*qs_cb, qcb=RH_cb*qs_cb, qs0=1.07*qs0, q0=(1-0.95(1-RH_0))*qs0
    for (i, qcb) in enumerate(fac .* RH_cb*qs_cb)
    Ccb_varRHcb[i] = cloudbase_cloudfrac(Dq = q0-qcb,
                    satdef_cb = qs_cb-qcb,
                    wq = 1.02 * Ecb/rhoL)
    end
end
plot(100*fac.*RH_cb, 100*Ccb_varRHcb)
ylabel("cloud base cloud fraction (%)")
xlabel("cloud base RH (%)")

# climate changes:
# E_0  +2 %/K
# q_s  +7 %/K
# 1-RH -5 %/K
# RH0 = 0.8, RH1 = 0.81