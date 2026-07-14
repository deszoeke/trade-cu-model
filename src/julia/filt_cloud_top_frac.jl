includet("TradeCuModel.jl")

using PythonPlot
using Printf
using Statistics

if @isdefined(PythonPlot)
    using PythonCall
    "Convert arrays with missing to numpy masked arrays without overriding global conversion."
    function as_masked_array(a::AbstractArray{Union{T,Missing},N}) where {T<:Real,N}
        np = PythonCall.pyimport("numpy")
        np.ma.array(coalesce.(a, NaN), mask=ismissing.(a))
    end

    "close plots with close(\"all\")"
    close(x) = PythonCall.pyimport("matplotlib.pyplot").close(x)
elseif @isdefined(PyPlot)
    using PyCall
    using PyCall: PyObject
    # allow for plotting with missing values
    function PyCall.PyObject(a::Array{Union{T,Missing},N}) where {T,N}
        numpy_ma = PyCall.pyimport("numpy").ma
        pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
    end
end

moving_average = TradeCuModel.moving_average
recurse        = TradeCuModel.recurse


# get An-Yi's aggregated GOES cloud fraction vs height data
# function get_goes_cloud_data()
#     NCDataset("../../data/satellite/GOES-16/goes16_binned_low4km_20200115_20200219.nc") do dsa
#         rfv_nrm = tmean(dsa[:rfv_nrm][:,:]) # skipmissing: don't let fill values on individual days blank low-height bins
#         rfv_acc = tmean(dsa[:rfv_acc][:,:])
#         cth_bin = 1e3 * dsa[:cth_bin][:] # km -> m
#         return rfv_nrm, rfv_acc, cth_bin
#     end
# end
rfv_nrm, rfv_acc, cth_bin = get_goes_cloud_data()

flt(x) = recurse(x->moving_average(x, m, good), x, n) # nx moving average filter

"bilinear interpolation of y(x) between (x1,y1) and (x2,y2)"
bilinear(x1,x2, y1,y2, x) = ( x2 == x1 ? y1 : y1 + (y2-y1) * (x-x1) / (x2-x1) )

"returns a function that interpolates y(x) from vectors X, Y with X descending"
function interp_extrap_descending( X::AbstractVector{<:Real}, Y::AbstractVector{<:Real} )
    function itp(x)
        # strictly NaN true out-of-bounds inputs
        # (x > X[1] || x < X[end]) && return NaN
        # search for descending vector order
        j = clamp(searchsortedfirst(X, x, rev=true), 2, length(X)) # clamped to data intervals [2, N]
        bilinear(X[j-1],X[j], Y[j-1],Y[j], x)
    end
    return itp
end

"returns a function that interpolates y(x) from vectors X, Y with X ascending"
function interp_extrap_ascending( X::AbstractVector{<:Real}, Y::AbstractVector{<:Real} )
    function itp(x)
        # strictly NaN true out-of-bounds inputs
        # (x < X[1] || x > X[end]) && return NaN
        # search for ascending vector order
        j = clamp(searchsortedfirst(X, x), 2, length(X)) # clamped to data intervals [2, N]
        bilinear(X[j-1],X[j], Y[j-1],Y[j], x)
    end
    return itp
end

# filter rfv_nrm and rfv_acc to reduce noise in cloud fraction
"""
filter rfv_nrm with moving average filter.
extrapolate at bottom, and set top cloud fraction to 0 (avoids undershooting 0 cloud fraction)
"""
function filt_rfv(rfv_nrm; n=3, m=5, zerotop=false)
    mp =  m÷2 * n
    idx = clamp.(1-mp:lastindex(rfv_nrm)+mp, 1,lastindex(rfv_nrm)) # pad ends
    flt(x) = recurse(x->moving_average(x, m, good), x, n) # nx moving average filter
    # extrapolate at bottom
    r = interp_extrap_ascending( eachindex(rfv_nrm), rfv_nrm ).(1-mp:lastindex(rfv_nrm)+mp) # [mp+1:lastindex(rfv_nrm)+mp]
    # set last value to last data value = 0
    r[lastindex(rfv_nrm)+mp:end] .= rfv_nrm[end] 
    fltidx = mp+1:lastindex(rfv_nrm)+mp
    y = flt(r)[fltidx]
    zerotop && (y[end] = 0.0) # set top cloud fraction to 0
    return y
end

n = 3 # number of times to apply moving average filter  
m = 5 # moving average window size
mp =  m÷2 * n
idx = clamp.(1-mp:lastindex(rfv_nrm)+mp, 1,lastindex(rfv_nrm)) # pad ends

clf()
plot(rfv_nrm, cth_bin/1e3, label="rfv_nrm")
# plot(flt(rfv_nrm[idx])[mp+1:end-mp], cth_bin/1e3, label="filtered")
plot(filt_rfv(rfv_nrm), cth_bin/1e3, label="filtered & extrapolated")
legend()

clf()
plot(rfv_acc, cth_bin/1e3, label="rfv_acc")
# plot(flt(rfv_nrm[idx])[mp+1:end-mp], cth_bin/1e3, label="filtered")
plot(filt_rfv(rfv_acc), cth_bin/1e3, label="filtered & extrapolated")
legend()

rfv_nrm_f = filt_rfv(rfv_nrm)
rfv_acc_f = filt_rfv(rfv_acc; zerotop=true) # set top cloud fraction to 0

