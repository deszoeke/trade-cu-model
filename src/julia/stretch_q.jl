using Revise
using Debugger

using PythonPlot, PythonCall
using Printf
using Statistics

# facilitate and customize plots
if @isdefined(PythonPlot) && @isdefined(PythonCall)
    "Convert arrays with missing to numpy masked arrays without overriding global conversion."
    function as_masked_array(a::AbstractArray{Union{T,Missing},N}) where {T<:Real,N}
        np = PythonCall.pyimport("numpy")
        np.ma.array(coalesce.(a, NaN), mask=ismissing.(a))
    end

    "close plots with close(\"all\")"
    close(x) = PythonCall.pyimport("matplotlib.pyplot").close(x)
elseif @isdefined(PyPlot) && @isdefined(PyCall)
    # using PyCall: PyObject
    # allow for plotting with missing values
    function PyCall.PyObject(a::Array{Union{T,Missing},N}) where {T,N}
        numpy_ma = PyCall.pyimport("numpy").ma
        pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
    end
end

includet("TradeCuExperiments.jl")
using .TradeCuExperiments
import .TradeCuModel.interpolate_ascending
import .TradeCuModel.interpolate_descending
import .TradeCuModel.interp_extrap_ascending
import .TradeCuModel.interp_extrap_descending
import .TradeCuModel.moving_average
import .TradeCuModel.recurse

pd = permutedims
good(x) = !ismissing(x) && isfinite(x)
f0(x) = good(x) ? x : 0

"""
filter x with moving average filter.
"""
function filt_q(x; n=5, m=7, )
    mp =  m÷2 * n
    good(x) = !ismissing(x) && isfinite(x)
    flt(x) = recurse(x->moving_average(x, m, good), x, n) # nx moving average filter
    y = flt(x)
    return y
end

# set up plot defaults
PythonPlot.matplotlib.rcParams["font.family"] = "sans-serif"
PythonPlot.matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "OpenSans"]
# make fonts bigger by mutating rcParams
font_settings = Dict(
    "font.size" => 14,       # Base size
    "axes.titlesize" => 18,  # Subplot titles
    "axes.labelsize" => 16,  # X/Y labels
    "xtick.labelsize" => 14, # X-axis numbers
    "ytick.labelsize" => 14, # Y-axis numbers
    "legend.fontsize" => 14  # Legend text
)
matplotlib.rcParams.update(font_settings)

( qm, qs, zcb, qcb, E_cb, x, divg, sfc_adv,
    tot_sink, cth_bin, rfv_acc, rfv_nrm, 
    rhoL, ns, nz ) = setup_experiments(ctx=ctx);
ctx = init_context()
z = ctx.z;

clf()
plot(qm*1e3,z/1e3)
ylim([0, 6])

x = log.(filt_q(qm)./maximum(f0,qm))
plot(x,z/1e3)
discrx = 0:-0.05:-3.5
discrz = interpolate_descending(x,z).(discrx)

plot(discrx, discrz)

stretch(z; a=300) = z + a*sin(2*pi*z/6e3)
stretch(z; a=300) = a*sqrt(sin(pi*z/6e3))
# stretch by multiplying
stretch(z; a=0.1) = a*max( 0, sin(pi*(z-700)/(4e3-700)) )

clf()
plot(discrx, discrz)
# plot(discrx, discrz .+ stretch.(discrz))
# plot(discrx, discrz .- stretch.(discrz))
plot(discrx, @. discrz * (1 + 0.1*sin(pi*discrz/6e3)))
plot(discrx, @. discrz * (1 - 0.1*sin(pi*discrz/6e3)))

# stretch by multiplying
clf()
plot(exp.(discrx),    discrz )
plot(exp.(discrx), @. discrz * (1 + 0.1*max(0, sin(pi*(discrz-700)/(4e3-700)))))
plot(exp.(discrx), @. discrz * (1 - 0.1*max(0, sin(pi*(discrz-700)/(4e3-700)))))

# look at derivatives
clf()
plot(diff(exp.(discrx))./diff(discrz),    discrz[2:end] )
plot( diff(exp.(discrx))./diff(@. discrz * (1 + stretch(discrz))), 
                              (@. discrz * (1 + stretch(discrz)))[2:end] )
plot( diff(exp.(discrx))./diff(@. discrz * (1 - stretch(discrz))), 
                              (@. discrz * (1 - stretch(discrz)))[2:end] )