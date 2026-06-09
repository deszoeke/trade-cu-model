# load modules
using Revise
using Pkg
# Pkg.activate(joinpath(homedir(), "Projects/ATOMIC/ATOMIC_GOES/julia/ATOMIC_GOES"))
cd("/Users/deszoeks/Projects/ATOMIC/trade-cu-model/src/julia")
Pkg.activate(".")

using Printf, Dates, NCDatasets, JLD2
using Statistics, Interpolations, LinearAlgebra

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