using Revise
using Pkg; Pkg.activate(".")

using Dates
using NCDatasets
using PyPlot
using Interpolations

using PyCall
# allow for plotting with missing values
function PyCall.PyObject(a::Array{Union{T,Missing},N}) where {T,N}
    numpy_ma = PyCall.pyimport("numpy").ma
    pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
end

"""
get the indices of y that correspond for each variable x
so y[index[i]] <= x[i] < y[index[i+1]]
"""
coarseind(x, y) = [searchsortedlast(y, xi) for xi in x]
#= "with type stability and preallocation"
function coarseind(x, y)
   indices = Vector{Int}(undef, length(x))
   for i in eachindex(x, indices)
       indices[i] = searchsortedlast(y, x[i])
   end
   return indices
end =#

#url = "https://www.ncei.noaa.gov/thredds-ocean/catalog/psl/atomic/p3/Picarro/catalog.html" # thredds
#url = "https://www.ncei.noaa.gov/data/oceans/oar/psl/atomic-2020/p3/Picarro/"
# ds = NCDataset(url)
p3dir = "../../data/p3"
p3files = readdir(joinpath(p3dir, "Picarro")) |> filter(startswith("EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-Isotope-Ratios-1Hz"))
p3h2ofiles = readdir(joinpath(p3dir, "Picarro")) |> filter(startswith("EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-1Hz"))
p3metfiles = readdir(joinpath(p3dir, "meteorology")) |> filter(startswith("EUREC4A_ATOMIC_P3_Flight-Level"))

#dsi = NCDataset( joinpath(p3dir, p3files[1]) )

"filter out takeoff and landing"
function filt_flight( dsi )
    # exclude take-off and landing to avoid sampling island air
    takeoff = findfirst(dsi[:alt][:] .>= 2500 )
    landing = findlast( dsi[:alt][:] .>= 2500 )
    t_h = (dsi[:time][:].-dsi[:time][takeoff])./Hour(1) # hours since first ascent
    ii = findall(    eachindex(t_h) .>= takeoff 
                 .&& eachindex(t_h) .<= landing 
                 .&& dsi[:lon][:] .>= -58.5     )
    return ii, t_h
end

"read flight isotope file, filter takeoff and landing"
function get_iso_ds( f )
    dsi = NCDataset( f )
    ii, t_h = filt_flight( dsi )
    return dsi, ii, t_h
end

m2n(x) = replace(x, missing => NaN)

function get_h2o_var_at_iso_times( dsw, varname, dsi, ii,
        epoch_ms_i = Dates.datetime2epochms.(dsi[:time][:]),
        epoch_ms_w = Dates.datetime2epochms.(dsw[:time][:]) )
    itp = linear_interpolation(epoch_ms_w, m2n(dsw[varname][:]))
    itp( epoch_ms_i[ii] )
end

function plot_flight_alt( dsi, ii )
    subplot(2,1,1)
    plot(dsi[:lon][ii], dsi[:alt][ii]/1e3)
    subplot(2,1,2)
    plot(t_h[ii], dsi[:alt][ii]/1e3)
end

function plot_iso_prof( dsi, ii )
    plot(dsi[:dD][ii], dsi[:alt][ii]/1e3, linewidth=0.2)
    ylabel("height (km)")
    xlabel(L"\delta"*"D "*L"(10^{-3})")
end

spechum(mixratio) = mixratio / (1+mixratio)
const RDvsmow = 155.76e-6     # unitless, deuterium
const R18Ovsmow = 2005.2e-6   # unitless, oxygen-18
delt(R, Rvsmow=RDvsmow) = 1e3*(R/Rvsmow - 1)     # permil
Rati(d, Rvsmow=RDvsmow) = Rvsmow * (1 + d*1e-3); # d in permil
alphae_d(T) = exp( 1158.8e-12 .*T.^3 - 1620.1e-9 .*T.^2 + 794.84e-6 .*T - 161.04e-3 + 2.9992e6./T.^3 )

"mixing line"
#mix(a, x0, x1) = a.*x1 + (1-a).*x0
function mix(x0,y0, x1,y1, n::Integer=20)
    a = (0:n)./n
    x = a.*x1 .+ (1.0 .-a).*x0
    y = a.*y1 .+ (1.0 .-a).*y0
    return x,y
end

"Rayleigh model isotope ratio R"
function Rayleigh(T=305.0:-5.0:220.0, R0=1.5e-4, q0=2e-2,
        RR=similar(T), q=similar(T) )
    Cpv = 1870.0 # J/kg/K
    Cw  = 4190.0
    Cp  = 1005.
    Rd  = 287.0
    C   = 273.15
    L0  = 2.501e6
    Rv  = 461.0
    "Lv(T) [J/kg] Latent heat of vaporization of water as a function of temperature [K]"
    Lv(T) = L0 + (Cpv-Cw) * (T-C)

    RR[1] = R0
    q[ 1] = q0
    for i in eachindex(RR[1:end-1])
        dlnT = log( T[i+1] / T[i] )
        dlnq = ( Lv(T[i])/(Rv*T[i]) + Cp/Rd ) * dlnT
        dlnR = ( alphae_d(T[i])-1 ) * dlnq
        RR[i+1] = RR[i] * exp(dlnR)
        q[i+1]  = q[i]  * exp(dlnq)
    end
    return RR, q # Rayleigh isotope ratio, specific humidity
end
RR, qR = Rayleigh()

dstring = [ match(r"\d{4}\d{2}\d{2}", f).match for f in joinpath.(p3dir, "Picarro", p3files) ]

# plot profiles for all flights
clf()
for (i,f) in enumerate( joinpath.(p3dir, "Picarro", p3files) )
    dsi, ii, t_h = get_iso_ds( f )
    plot_iso_prof( dsi, ii )
end
xlim([-500, -65])
ylim([0, 4])

# loglog plot q*dD vs q for all flights
clf()
for (i,f) in enumerate( joinpath.(p3dir, "Picarro", p3files) )
    dsi, ii, t_h = get_iso_ds( f )
    h2ofile = filter(contains(dstring[i]), p3h2ofiles)[1] # matching water vapor file
    dsw = NCDataset( joinpath(p3dir, "Picarro", h2ofile) )
    q  = spechum.( get_h2o_var_at_iso_times( dsw, :mmr, dsi, ii ).*1e-3 )
    rh =           get_h2o_var_at_iso_times( dsw, :rh_iso, dsi, ii )
    loglog(q, -dsi[:dD][ii].*q, linestyle="none", marker=".", markersize=0.2) # weird
end
plot([1e-5, 2e-2], (5e-3.*[1e-5, 2e-2]./1e-5).^0.8 , color="k", linewidth=0.3)
plot([1e-5, 2e-2], (7e-4.*[1e-5, 2e-2]./1e-5).^0.8 , color="k", linewidth=0.3)
plot([1e-5, 2e-2], (5e-3.*[1e-5, 2e-2]./1e-5).^0.6 , color="k", linewidth=0.3)
plot(6e-3.+[0,0], [2e-1, 4e-1], color="k", linewidth=0.3)
text(7e-3, 2e-1, "trade cumulus")
xlabel("specific humidity q")
ylabel("q*dD "*L"10^{-3}")

ylim([0.1, 6.5])
xlim([1e-3, 2e-2])

# loglog plot qi vs q for all flights
clf()
for (i,f) in enumerate( joinpath.(p3dir, "Picarro", p3files) )
    dsi, ii, t_h = get_iso_ds( f )
    h2ofile = filter(contains(dstring[i]), p3h2ofiles)[1] # matching water vapor file
    dsw = NCDataset( joinpath(p3dir, "Picarro", h2ofile) )
    q = spechum.( get_h2o_var_at_iso_times( dsw, :mmr, dsi, ii ).*1e-3 )
    rh =          get_h2o_var_at_iso_times( dsw, :rh_iso, dsi, ii )
    icld = rh .> 90.0
    icldp = rh .> 97.0
    qi = Rati.(dsi[:dD][ii]) .* q
    semilogx(q, log.(Rati.(dsi[:dD][ii])./RDvsmow), linestyle="none", marker=".", markersize=0.2, color="r")
    semilogx(q[icld], log.(Rati.(dsi[:dD][ii[icld]])./RDvsmow), linestyle="none", marker=".", markersize=0.2, color="b")
    semilogx(q[icldp], log.(Rati.(dsi[:dD][ii[icldp]])./RDvsmow), linestyle="none", marker=".", markersize=0.4, color="c")
end
semilogx(qR, log.(1.08*RR/RDvsmow), linewidth=0.4, color="k") # ~ +80 permil
semilogx(qR, log.(0.97*RR/RDvsmow), linewidth=0.4, color="k") # ~ -30 permil
qm,Rm = mix(qR[1],0.97*RR[1], qR[6], 1.08*RR[6])
semilogx(qm, log.(Rm/RDvsmow), linewidth=0.4, color="k", linestyle="--")# 
qm,Rm = mix(qR[5],0.97*RR[5], qR[8], 1.08*RR[8])
semilogx(qm, log.(Rm/RDvsmow), linewidth=0.4, color="k", linestyle="--")# 
xlim([9e-4, 2e-2])
ylim([-0.25, 0])
title("deuterium")
ylabel(L"\ln(R/R_s)\approx\delta")
xlabel("specific humidity (kg kg\$^{-1}\$)")
tight_layout()

fmts = ["png", "pdf", "eps", "svg"]
[ savefig("P3_RD_vs_q.$(f)") for f in fmts ]

# Normalize each profile by a maximum R0, rather than RDvsmow
# to collapse profiles.
clf()
for (i,f) in enumerate( joinpath.(p3dir, "Picarro", p3files) )
    dsi, ii, t_h = get_iso_ds( f )
    h2ofile = filter(contains(dstring[i]), p3h2ofiles)[1] # matching water vapor file
    dsw = NCDataset( joinpath(p3dir, "Picarro", h2ofile) )
    q = spechum.( get_h2o_var_at_iso_times( dsw, :mmr, dsi, ii ).*1e-3 )
    rh =          get_h2o_var_at_iso_times( dsw, :rh_iso, dsi, ii )
    icld = rh .> 90.0
    icldp = rh .> 97.0
    qi = Rati.(dsi[:dD][ii]) .* q
    R0 = Rati(maximum(skipmissing(dsi[:dD][ii])))
    print("$(R0), ")
    semilogx(q, log.(Rati.(dsi[:dD][ii])./R0), linestyle="none", marker=".", markersize=0.2, color="r")
    semilogx(q[icld], log.(Rati.(dsi[:dD][ii[icld]])./R0), linestyle="none", marker=".", markersize=0.2, color="b")
    semilogx(q[icldp], log.(Rati.(dsi[:dD][ii[icldp]])./R0), linestyle="none", marker=".", markersize=0.4, color="c")
end
semilogx(qR, log.(1.12*RR/RR[1]), linewidth=0.4, color="k") # ~ +80 permil
semilogx(qR, log.(RR/RR[1]), linewidth=0.4, color="k") # ~ -30 permil
qm,Rm = mix(qR[1],RR[1], qR[6], 1.12*RR[6])
semilogx(qm, log.(Rm/RR[1]), linewidth=0.4, color="k", linestyle="--")#
qm,Rm = mix(qR[5],RR[5], qR[8], 1.12*RR[8])
semilogx(qm, log.(Rm/RR[1]), linewidth=0.4, color="k", linestyle="--")#
xlim([1e-3, 2e-2])
ylim([-0.13, 0])
title("deuterium")
ylabel(L"\ln(R/R_s)\approx\delta")
xlabel("specific humidity (kg kg\$^{-1}\$)")
tight_layout()

[ savefig("P3_RDoR0_vs_q.$(f)") for f in fmts ]
