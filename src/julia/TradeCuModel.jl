module TradeCuModel

# load modules
using Revise # does using Revise first allow for revising code in ATOMIC_GOES???
using Pkg
# Pkg.activate(joinpath(homedir(), "Projects/ATOMIC/ATOMIC_GOES/julia/ATOMIC_GOES"))
cd("/Users/deszoeks/Projects/ATOMIC/trade-cu-model/src/julia")
Pkg.activate(".")

using Printf
using Dates
using NCDatasets
using JLD2

using Statistics
using Interpolations
using LinearAlgebra

# using PyPlot -> not thread safe in IJulia with Julia 1.12.4+
using PythonPlot

using VaporSat # dev ../../deps/VaporSat

# exported functions
export similarmissing, n2m, nisf2m
export mid, m2z, anom, runningmean, recursef, upsample1, upsample, ξ
export ddz, q_total, cloudflux!
export precipflux_down!, precipflux_down, precipflux_down_sfc
export calcF2
export cloudflux_1x
export updraft_w_dq
export dadsinkrate, find_contour!
export tmean, tstd
export get_sounding_dataset, get_mean_soundings, get_goes_cloud_data
export virtual_temp, calc_rhoL, Lv, LvK
# export largescale_drying

# utility functions; extend matplotlib

if @isdefined(PythonPlot)
    using PythonCall
    "Recast julia Arrays with missing as python masked arrays."
    function PythonCall.Py(a::Array{Union{T,Missing},N}) where {T,N}
        np = PythonCall.pyimport("numpy")
        np.ma.masked_invalid(np.array(replace(a, missing => NaN)))
    end

    "close plots with close(\"all\")"
    close(s::String) = PythonCall.pyimport("matplotlib.pyplot").close(s)
elseif @isdefined(PyPlot)
    using PyCall
    using PyCall: PyObject
    # allow for plotting with missing values
    function PyCall.PyObject(a::Array{Union{T,Missing},N}) where {T,N}
        numpy_ma = PyCall.pyimport("numpy").ma
        pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
    end
end

"convert arrays with missings to masked python float64 arrays"
# msk(y) = ma.masked_invalid(np.array(replace(y, missing => NaN), dtype="float64"))


"plotf(f, x, args...)  function plot wrapper"
plotf(f, x, args...) = plot(f.(x), x, args...)

pd = permutedims

"midpoint averager"
mid(x) = 0.5*(x[1:end-1].+x[2:end]) # doesn't work for Int16

"missing to zero"
# m2z(x) = ismissing(x) ? 0 : x
m2z(x) = replace(x, missing=>0)

# statistics fcns
"anomaly"
anom(r, kwargs...) = r .- mean(r, kwargs...)

"runningmean(x, n=7, f=isfinite) n-point running mean points satisfying f(x)==true"
function runningmean(x, n=7, f=isfinite)
    notf2zero(x) = f(x) ? x : 0
    cs = cumsum(notf2zero.(x))
    ct = cumsum(f.(x))
    (cs[n+1:end] .- cs[1:end-n]) ./ (ct[n+1:end] .- ct[1:end-n])
end

recurse(f, x, n) = n>0 ? f(recurse(f, x, n-1)) : x
"abstractly composes f n times; recursef(f,3) = f∘f∘f = x->f(f(f(x))), without evaluating"
recursef(f,n) = reduce(∘, ntuple(_ -> f, n)) 

"declare arrays with elements of Union{Missing,T} unions"
similarmissing(s::Tuple, T::Type=Float64) = Array{Union{Missing, T}}(missing, s)
similarmissing(A::Array) = similarmissing(size(A), eltype(A))

"moving_average(A, m) m-point moving average of A"
function moving_average(A::AbstractArray, m::Int, cond=x->true)
    out = similar(A)
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = m÷2 * oneunit(Ifirst)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I-I1):min(Ilast, I+I1)
            if cond(A[J]) # add condition, e.g. isfinite
                s += A[J]
                n += 1
            end
        end
        out[I] = s/n
    end
    return out
end

# means and helper functions

n2m(x) = isnan(x) ? missing : x
skipnan(x) = skipmissing(n2m.(x))

"replace nonmissing nonfinite by missing"
nisf2m(x) = !ismissing(x) && isfinite(x) ? x : missing
# usem e.g., mean( skipmissing(nisf2m.(x)) )

"upsample by factor n, between 2 points"
upsample1(x0, x1; n=10, a=(0:n)/n) = x0*(1-a) + x1*a

"upsample X by factor n, between points of vector x."
function upsample(x; n=10, a=(0:n)/n)
    X = zeros((length(x)-1)*n + 1)
    for i in eachindex( x[1:end-1] )
        X[n*(i-1)+1 .+ (0:n)] .= x[i]*(1 .- a) + x[i+1]*a
    end
    X
end

"\\xi(x) wrapper for upsampling linearly in log space"
ξ(x) = exp.(upsample(log.(x)))

"vertical derivative"
ddz(q,z) = diff(q) ./ diff(z)

# functions for thermodynamics

# using .VaporSat -> dev ../../deps/VaporSat
# VaporSat exports dqsdT, qsat
# here we overload some unexported methods in VaporSat

# thermodynamic constants
Rd=287.04 # J/K/kg Bolton
Cp=1005.7 # J/K/kg Bolton
Rv=461.5 # Bolton
# specific heats of vapor and liquid
Cpv=1870 # Bolton
Cw=4190  # Bolton
gravity=9.8 #
KelvinCelsius=273.15 # K

# Lv(T) = 2.501e6 + (Cpv-Cw) * T # Celsius
"Lv(T) [J/kg] Latent heat of vaporization of water as a function of temperature [K]."
Lv(T) = 2.501e6 + (Cpv-Cw) * (T-KelvinCelsius)
LvK(T) = Lv(T-273.15)

"""
qv(p/ev) = Rd / (Rv * (p/ev + (Cp/Rv-1)))
specific humidity [kg/kg]
"""
qv(poev) = Rd / (Rv * (poev + (Cp/Rv-1.)))
qv(p,ev) = qv(p/ev)

"vapor pressure in units of p, qv unitless"
ev(p, qv) = p*qv / (Rd/Rv + qv)

"""
esat(T,p) [hPa] computes saturation vapor pressure based on Wexler's formula,
with enhancement factor for moist air rather than water vapor.
The enhancement factor requires a pressure.
T [degrees C], p [hPa] (note the reversed input order), es [hPa]
From A. L. Buck 1981: JAM, 20, 1527-1532.
SPdeS 7 July 2004
"""
esat(T, p) = 6.1121*(1.0007 + 3.46e-8*p)*exp((17.502*T)/(240.97 + T))
#              hPa

# use method from VaporSat qsat(p[Pa], T[C])
# note order of arguments and units
# "qsat(T[C], p[hPa])"
# qsat(T, p) = qv(p/esat(T, p))

"specific humidity q [unitless] from mixing ratio w [unitless]"
q_from_w(w) = w/(1+w)

"""
Tlcl(T[K], ev[Pa])
Temperature at the LCL [K]. From Bolton, 1980, MWR, 108, 1046-1053.
"""
Tlcl(T, ev) = 2840. / (3.5*log(T) - log( 0.01*ev ) - 4.805) + 55.

# adiabatic lapse rate
lapse_ad = gravity/Cp

"moist adiabatic lapse rate p[Pa], T[K]"
lapse_moist(T,p) = lapse_ad * (1 + LvK(T)*qsat(p,T-KelvinCelsius)/(Rd*T)) / (1 + Lv(T)^2*Rd/Rv*qsat(p,T-KelvinCelsius)/(Rd*Cp*T*T))

"parcel process adiabatic change in saturation vapor pressure, p[Pa], T[K]"
dqsatdz_moistad(p,T) = dqsdT(p,T-KelvinCelsius) * -lapse_moist(T,p) # negative
# Tolerate small inconsistencies between dqsdT (Bolton) 
# and lapse_moist (Rogers and Yau) methods.

"density temperature for liquid water loading (exact, or take default approximation)"
liqldfac(ql; rhoa_o_rhoL=0.001) = 1 - (1 - rhoa_o_rhoL) * ql
Tliqld(Tv, ql; p=1e5, rhoa=p/(Rd*Tv)) = Tv * liqldfac(ql; rhoa_o_rhoL=rhoa/1000.0)

calc_rho( T, p) = p/(Rd*T)
calc_rhoL(T, p) = p/(Rd*T) * LvK(T) # should use Tv in p/(Rd*T) but T in LvK
# Tmean = mean(skipnan(tmean(:ta)[ds[:alt].<=ztop])) # 286.6K = mean Cu layer temperature
# pmean = mean(skipnan(tmean(:p)[ds[:alt].<=ztop]))
# rhoL = calc_rhoL(Tmean, pmean)

# functions compute clouds, moisture and mass flux
# from q, qs sounding, cloud fraction, 
# Tmean, pmean

"virtual temperature, K"
virtual_temp(t,q) = t * (1 + 0.608*q)

# functions for getting sounding data
function get_sounding_dataset(; datapath = joinpath(homedir(),"Data/ATOMIC/radiosonde/level2"))
    # readdir(datapath)
    ncfile = filter(
        x->startswith(x,"EUREC4A_RonBrown") && endswith(x,"Vaisala-RS_L2_v3.0.0.nc"), 
        readdir(datapath))
    ds = NCDataset( joinpath(datapath, ncfile[1]) )
end

"average all sounding times. preserve levels"
# function tmean(var, ds=get_sounding_dataset())
#     [mean(skipmissing(ds[var][i,:])) for i in axes(ds[var][:,:],1)]
# end

tmean(var::Symbol, ds=get_sounding_dataset()) = tmean(ds[var][:,:])
tmean(var::String, ds=get_sounding_dataset()) = tmean(ds[var][:,:])
tstd(var::Symbol, ds=get_sounding_dataset()) = tstd(ds[var][:,:])
tstd(var::String, ds=get_sounding_dataset()) = tstd(ds[var][:,:])

function tmean(data::AbstractArray)
    # Preallocate the result matrix with the correct element type excluding Missing
    T = nonmissingtype(eltype(data)) 
    rows = size(data, 1)
    result = Vector{T}(undef, rows)
    
    # Loop over rows efficiently using @view
    @inbounds for i in 1:rows
        result[i] = mean(skipmissing(@view data[i, :]))
    end
    return result
end

function tstd(data::AbstractArray)
    # Preallocate the result matrix with the correct element type excluding Missing
    T = nonmissingtype(eltype(data)) 
    rows = size(data, 1)
    result = Vector{T}(undef, rows)
    
    # Loop over rows efficiently using @view
    @inbounds for i in 1:rows
        result[i] = std(skipmissing(@view data[i, :]))
    end
    return result
end

# ds[:alt][:] # 3100, 10-m vert resolution
# ds[:ta][:,:] # 3100×277 Matrix{Union{Missing, Float32}}

# get mean sounding data
function get_mean_soundings(ds=get_sounding_dataset())
    z = ds[:alt][:]
    tam = tmean(:ta, ds)
    thm = tmean(:theta, ds)
    qm  = tmean(:q, ds)
    pm  = tmean(:p, ds)
    return z, tam, thm, qm, pm
end

# get An-Yi's aggregated GOES cloud fraction vs height data
function get_goes_cloud_data()
    NCDataset("../../data/goes16_binned_low4km_20200115_20200219.nc") do dsa
        rfv_nrm = mean(dsa[:rfv_nrm][:,:], dims=2)[:]
        rfv_acc = mean(dsa[:rfv_acc][:,:], dims=2)[:]
        cth_bin = dsa[:cth_bin][:]
        return rfv_nrm, rfv_acc, cth_bin
    end
end

# functions for the cloud model

# cloud model - updraft total water and cloud water functions
# Updraft q_total, ql do not depend on F or Fcld.
# Below, mass flux M depends on q_total and Fcld.

"qtilde=(1-x)*qm .+ x*qs intermediate environmental moisture function"
function q_tilde(qm, qs, x)
    (1-x)*qm .+ x*qs
end

"updraft q_total function"
function q_total(aedz, x, qs, qm; i0, qt0, stopattop=false)
    qtilde=(1-x)*qm .+ x*qs
    qt = copy(qtilde)
    qt[1:i0] .= qt0 # cloud base initial condition
    # in cloud
    for i in i0:length(qt)-1
        if qt[i] > qs[i] # in cloud
            a = aedz
            Q = 0.5*(qtilde[i]+qtilde[i+1])
            qt[i+1] = ((1-0.5*a)*qt[i] + a*Q) / (1+0.5*a)
        elseif !stopattop   # above cloud
            a = (1-x)*aedz  # no precipitation sink
            Q = 0.5*(qm[i]+qm[i+1])
            qt[i+1] = ((1-0.5*a)*qt[i] + a*Q) / (1+0.5*a)
        elseif stopattop  # above cloud, stop at cloud top
            qt[i+1] = -999
        end
    end
    return qt
end

# cloud flux Fcld functions
# ql needs to be calculate by the cloud functions already
"""
Iterates up cloud flux modifying Fcld with precipitation efficiency x and total flux F.
Must initialize cloud base Fcld[icb] first!
"""
function cloudflux!(Fcld, x, ae, F, ql, qt, qm, icb=icb, itop=min(lastindex(F),findfirst(.!ismissing.(F) .&& F.<0)+20))
    Δq = qt .- qm
    c = @. -x*ae * ql / Δq * dz/2 # Q*dz stepping quantity # /2 added 2025 Aug 16
    # step Fcld through the trade Cu layer
    for i in icb:itop
        Fcld[i+1] = (F[i+1]-F[i] + Fcld[i]*(1+c[i])) / (1-c[i+1])                
    end
end

"integrate precipitation down to cloud base by stepping"
function precipflux_down_sfc(x, ae, F, ql, qt, qm; istart=maximum(findall(ql.>0))+1, icb=icb)
    c = @. x*ae * ql/(2*(qt-qm)) * dz
    #P[istart] = 0.0
    P = 0.0 # at cloud top istart
    for i in (istart-1):-1:icb
        # update P by 1 step
        # P[i] = ( (1+c[i+1])*P[i+1] + c[i]*F[i] + c[i+1]*F[i+1] ) / ( 1-c[i] )
        P = ( (1+c[i+1])*P + c[i]*F[i] + c[i+1]*F[i+1] ) / ( 1-c[i] )
    end
    return P # at cloud base
end

"precipitation profile integrated downward from cloud top"
function precipflux_down!(P, x, ae, F, ql, qt, qm; istart=maximum(findall(ql.>0))+1, icb=icb)
    c = @. x*ae * ql/(2*(qt-qm)) * dz
    # P = zeros(size(F)) # at cloud top istart
    P[istart] = 0.0
    for i in (istart-1):-1:icb
        # update P by 1 step
        P[i] = ( (1+c[i+1])*P[i+1] + c[i]*F[i] + c[i+1]*F[i+1] ) / ( 1-c[i] )
    end
    P[1:icb-1] .= NaN
end

"precipitation profile integrated downward from cloud top"
function precipflux_down(x, ae, F, ql, qt, qm; istart=maximum(findall(ql.>0))+1, icb=icb, dz)
    c = @. x*ae * ql/(2*(qt-qm)) * dz
    P = zeros(size(F)) # at cloud top istart
    for i in (istart-1):-1:icb
        # update P by 1 step
        P[i] = ( (1+c[i+1])*P[i+1] + c[i]*F[i] + c[i+1]*F[i+1] ) / ( 1-c[i] )
    end
    P[1:icb-1] .= NaN # below cloud base
    return P # precip flux vector
end


# large-scale flux
# calculate F2

# ΔG = -diff(G) # z grid
# prototype
# F2 = -diff(G)[offset.+eachindex(rfv_nrm)] ./ rfv_nrm[:] # dimension category_i
# rfv_acc[end] == 0.0

# line up model grid and cloud top height bins
# offset = findfirst(x->x≈cth_bin[1]*1e3, z) - 1 # 50; cth_bin starts at z=500

calcF2_(G, acc; sk=1) = diff(G[1:sk:end]) ./ diff(acc[1:sk:end]) # dimension category_i, inputs asligned!

"""
calcF2(G, rfv_acc)
Vertically uniform cloud eddy flux 
The all-sky mean flux profile ``G(z)`` has a flux difference 
at each height ``h_i`` of ``(ΔG)_i=-[G(h_{i+1})-G(h_i)]``. 
This flux difference is assumed to be carried by the clouds with 
area fraction ``a_i`` and height 
``h_i``, so the mean flux in cloud category ``i`` is vertically uniform, 
``F_i=a_i^{-1} (ΔG)_i``.

This is the flux for each cloud category ``i``. It is uniform from cloud base to ``h_i`` and goes to zero above ``h_i``.
"""
function calcF2(G, rfv_acc, offset; sk)
    align2i(G) = G[offset.+eachindex(rfv_acc)] # dimension category_
    return calcF2_(align2i(G), rfv_acc; sk=sk) # inputs not aligned
end


# iterate cloud model for F2z_i

# "terminate the cloud at ql<=0"
# findcloudtop(ql, z; zcb=2000) = findfirst(ql .<= 0 .&& z .> zcb )

"terminate the cloud either at ql=0 or at the elevated minimum ql in the trade inversion"
function findcloudtop(ql, z; zcb, ist=findfirst(z .>= 2000), ien=findfirst(z .>= 5000))
    ql = coalesce.(ql, NaN)
    if all(isnan, ql)
        return 0
    else
        itop = findfirst(ql .<= 0 .&& z .> zcb )
        if isnothing(itop) # choose the minimum in 2000<=z<=5000
            itop = argmin(ql[ist:ien]) + ist - 1
        end
    end
    itop
end

# detect deep convective clouds
"max ql of low clouds, below 3600 m"
low_argmax(ql,z) = argmax(skipmissing(replace(ql, NaN=>missing) .* (z.<3600)))

"Find the minimum in cloud above the trade inversion. Shallow Cu cloud top should be below this."
function high_argmin(ql,z) 
    lax = low_argmax(ql,z)
    lst = findlast( z .< 3600 )
    argmin(skipmissing(replace(ql[lax:lst], NaN=>missing))) + lax - 1
end

"return true for deep clouds having no ql <= 0 between max ql of low cloud and 4 km"
deep_cloud(ql,z) = minimum(ql[low_argmax(ql,z):findfirst(z.>=4000)]) > 0.0

# compute a measure of cloud fraction da/dsink = da/dh * dh/dsink

"return the first unique indices of a"
uniqueind(a) = unique(i -> a[i], eachindex(a))

function sortunique(ztop)
    ui_ = uniqueind(ztop[:])
    oi = sortperm(ztop[ui_])
    ui = ui_[oi]
end

"""
dadsinkrate(ztop, tot_sink, rfv_nrm, dh=10.0)
Compute cloud fraction density per unit sink rate 
da/dsinkrate = da/dh * dh/dsinkrate.
Interpolates satellite coordinate to model sinkrate coordiante.
"""
function dadsinkrate(ztop, tot_sink, cth_bin, rfv_nrm, dh=10.0)
    # condition vertical coordinate
    ui = sortunique(ztop)
    h = ztop[ui] # sorted ascending
    hm_ = (h[1:end-1] .+ h[2:end]) / 2
    ii = hm_ .<= 4000 # hardwired !!
    hm = hm_[ii]

    dh_dsink = diff(h) ./ diff(tot_sink[ui])

    # satellite coordinate --> da/dh
    da_dh = rfv_nrm[:] / dh
    # interpolate satellite a(h) to unique h_sink coordinate
    da_dh_to_sink = linear_interpolation((cth_bin*1e3,), da_dh).(hm) 
    da_dsink = da_dh_to_sink .* -dh_dsink[ii] # truncates just to useful sink_rate bins
    # println("sum(ii) = $(sum(ii))")
    da_dsink, findall(ii) # return the indices too
end

# usage later:
# da_dsink = dadsinkrate(ztop, tot_sink, cth_bin, rfv_nrm)

# refine cloud top for x=0.53

# nx = 500
# pcp_effic = range(0, 1, length=nx) #0:0.005:0.9         # precipitation efficiency
# ns = 600 # number of sink rates to test
# tot_sink = (1 .+tanh.(range(-8*pi, 0, length=ns))) .* (5e-3 - 1e-4) .+ 1e-4
# tot_sink = range(6.1716e-4, 5.8e-3, length=ns) # min tuned for x=0.53 to get the highest possible cloud top
# tot_sink = range(6.3523e-4, 5.7e-3, length=ns) # min tuned for x=0.53 to get the highest possible cloud top

"Compute clouds and fluxes for a single precipitation efficiency x and a range of sink rates tot_sink."
function cloudflux_1x(tot_sink=tot_sink; x=x, 
    z, nz=length(z), dz=z[2]-z[1],
    qm=qm, qs=qs, F2z=F2z, icb=icb, qcb=qcb)

    ns = length(tot_sink)
    qtc  = Array{ Union{Missing, Float64},2}(missing, nz, ns)
    Fcld = Array{ Union{Missing, Float64},2}(missing, nz, ns)
    Fp   = Array{ Union{Missing, Float64},2}(missing, nz, ns)
    # Pcb  = Matrix{Union{Missing, Float64}}(missing, ns)
    ztop = Vector{Union{Missing, Float64}}(missing, ns)
    # zmin = Matrix{Union{Missing, Float64}}(missing, ns)
    ff   = Vector{Union{Missing, Float64}}(missing, nz)

    # ix = findfirst( pcp_effic .>= x)
    qtilde = q_tilde(qm, qs, x) # depends on precipitation efficiency
    for ia in eachindex(tot_sink)
        ae = tot_sink[ia]
        # compute clouds
        qt = q_total(dz*ae, x, qs, qm; i0=icb, qt0=qcb)
        qd = qt.-qs
        ql = max.(0, qd)
        itop = findcloudtop(ql,z; zcb=z[icb])
        if !isnothing(itop) 
            ztop[ia] = z[itop] # ztop can be up to 20 km
            if !deep_cloud(ql,z)
                qtc[:,ia] .= qt
                fill!(ff, F2z[itop]) # use uniform flux for the experiment's cloud top height
                Fp[:,ia] .= -precipflux_down( x, ae, ff, ql, qt, qm, istart=itop, icb=icb, dz=dz )
                # Pcb[ia] = -Fp[icb,ia]
                # cloud updraft flux
                Fcld[:,ia] .= -Fp[:,ia] .+ F2z[itop]
                # cloud        = -precip       eddy
            end
        end
    end
    ztop, Fcld, Fp, qtc
end


"compute w for a single x and range of sink rates"
function updraft_w_dq(Fcld, qtc, qm, z, ztop)
    dq = qtc .- qm
    w = convert(Matrix{Union{Missing,Float64}}, Fcld ./ dq)
    for i in axes(w,2) # blank out w above clouds
        ii = z .> ztop[i]
        w[ii,i] .= missing
    end
    return w, dq
end

# example usage:
# w53, dq53 = updraft_w_dq(Fcld53, qtc53, qm, z, zt53)

# compute da/dsink for x=0.53
# da_dsink = dadsinkrate(zt53, tot_sink, cth_bin, rfv_nrm)

# Find the total sink rate at which the cloud top reaches each z in the grid.
"""
    find_contour!(X_out, x, y, q, c)
where q=c for each y in the grid. X_out has dimensions of y.
"""
function find_contour!(X_out, x, q, c)
    Nx, Ny = size(q)
    peaks = [argmax([isnan(v) || ismissing(v) ? -Inf : v for v in q[i, :]]) for i in 1:Nx]

    for j in 1:Ny
        idx = searchsortedlast(reverse(@view q[:, j]), c)
        (idx == 0 || idx == Nx) && continue
        i = Nx - idx # index in verso
        if j > peaks[i] && !isnan(q[i,j]) && !isnan(q[i+1,j])
            X_out[j] = x[i] + (x[i+1] - x[i]) * (c - q[i,j]) / (q[i+1,j] - q[i,j])
        end
    end
end
function find_contour(x,y, q, c=0) # returning wrapper
    X_out = Vector{Float64}(missing, length(y))
    find_contour!(X_out, x, q, c)
    return X_out
end
# interpolate the total sink rate that gives cloud top at each z.
# iz = findall(zcb .< z .<= ztop)
# sinkz = Vector{Float64}(missing, length(z))
# find_contour!(sinkz[iz], tot_sink, z[iz], pd((qtc.-qs)[iz,:]), 0)
# plot(sinkz, z, "k")

"Plot cloud profiles for a single precipitation efficiency x and a range of sink rates tot_sink."
function plot_qv_qc_w(tot_sink, z, da_dsink, ztop, qtc, Fcld, w, qm, qs; x=0.53, ntp=500)
    if @isdefined(PythonPlot) PyPlot=PythonPlot end

    # vertical coordinate
    ui = sortunique(ztop)
    h = ztop[ui] # sorted ascending
    hm_ = (h[1:end-1] .+ h[2:end]) / 2
    ii = hm_ .<= 4000
    hm = hm_[ii]

    dq = qtc .- qm
    
    figure(figsize=(9, 11))

    subplot(5,2,1)
    contourf(tot_sink*1e3, z[1:ntp]/1e3, qtc[1:ntp,:]*1e3, cmap=ColorMap("RdYlBu_r"))
    cb = colorbar()
    plot(tot_sink*1e3, ztop[:]/1e3, "k", label="cloud top height")
    # plot(tot_sink[ui][1:end-1][ii]*1e3, 0.02*da_dsink, "w", linewidth=0.5, label="da/d(sink)")
    ylim([0, 4])
    title("total cloud specific humidity (g kg\$^{-1}\$)")
    text(0.05, 3.4, "a")

    subplot(5,2,3)
    contourf(tot_sink*1e3, z[1:ntp]/1e3, dq[1:ntp,:]*1e3, levels=0:0.5:10, vmin=0, vmax=10, cmap=ColorMap("RdYlBu_r"))
    cb = colorbar()
    contour(tot_sink*1e3, z[1:ntp]/1e3, max.(0,qtc[1:ntp,:].-qs[1:ntp])*1e3, levels=0:0.2:3, vmin=0, vmax=3, colors="k", linewidths=0.5)
    plot(tot_sink*1e3, ztop[:]/1e3, "w")
    ylim([0, 4])
    title("updraft\$-\$environment humidity,\ncloud water (g kg\$^{-1}\$)")
    # xlabel("sink rate (km\$^{-1}\$)")
    ylabel("z (km)")
    text(0.05, 3.4, "b")

    subplot(5,2,5)
    levs = 0:0.05:0.5
    norm = PyPlot.matplotlib.colors.BoundaryNorm(levs, 256)
    pcolormesh(tot_sink*1e3, z[1:ntp]/1e3, w[1:ntp,:], norm=norm, cmap=ColorMap("RdYlBu_r"))
    colorbar()
    plot(tot_sink*1e3, ztop[:]/1e3, "k", label="cloud top height")
    ylim([0, 4])
    xlabel("sink rate "*L"\alpha + \epsilon"*"(km\$^{-1}\$)")
    title("cloud-mean vertical velocity (m s\$^{-1}\$)")
    text(0.05, 3.4, "c")

    tight_layout()
    return gcf()
end

end # module 
