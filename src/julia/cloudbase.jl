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
qstd = tstd(ds[:q][:,:])

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

rhoL = mean(filter(isfinite, calc_rhoL.(tvm, pm)[z.<=ztop])) # 2.41e6 J/m^3

# cloud base water flux from vapor and precipitation flux
# rhb_prate = mean(skipmissing(psl["prate"][:])) / 3600
rhb_prate = 8.88e-6 # kg s-1
Ecb = 180 # W/m^2; E0 - 35; just the turbulent cloud vapor flux
wq = Ecb/rhoL
Dq = qm[isfinite.(qm)][1] - qm[icb] # 0.0015 kg/kg
# cloud base condensation fraction
Ccb = 0.5 * (1 - erf((qs[icb]-qm[icb])/Dq))

function cloudbase_cloudfrac(; qm=qm, qs=qs, wq=Ecb/rhoL, icb=icb)
    Dq = qm[isfinite.(qm)][1] - qm[icb]
    Ccb = 0.5 * (1 - erf((qs[icb]-qm[icb])/Dq))
    return Ccb
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