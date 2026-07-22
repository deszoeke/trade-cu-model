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

qmax = maximum(f0,qm)
x = log.(filt_q(qm)./qmax)
plot(x,z/1e3)

# discretize x = ln humidity for stretching
discrx = 0:-0.05:-3.5
discrz = interpolate_descending(x,z).(discrx)
discrz[1] = z[findfirst(isfinite,qm)]

plot(discrx, discrz)

stretch(z; a=300) = z + a*sin(2*pi*z/6e3)
stretch(z; a=300) = a*sqrt(sin(pi*z/6e3))
# stretch below 4 km by multiplying
stretch(z; zcb=700.0, zct=4000.0) = max( 0, sin(pi*(z-zcb)/(zct-zcb)) )

# clf()
# plot(discrx, discrz)
# plot(discrx, @. discrz * (1 + 0.1*sin(pi*discrz/6e3)))
# plot(discrx, @. discrz * (1 - 0.1*sin(pi*discrz/6e3)))

# stretch by multiplying
clf()
plot(exp.(discrx),    discrz )
plot(exp.(discrx), @. discrz * (1 + 0.1*stretch(discrz)))
plot(exp.(discrx), @. discrz * (1 - 0.1*stretch(discrz)))

# look at derivatives
clf()
plot( diff(exp.(discrx))./diff(discrz),    discrz[2:end] )
plot( diff(exp.(discrx))./diff(@. discrz * (1 + 0.1*stretch(discrz))), 
                              (@. discrz * (1 + 0.1*stretch(discrz)))[2:end] )
plot( diff(exp.(discrx))./diff(@. discrz * (1 - 0.1*stretch(discrz))), 
                              (@. discrz * (1 - 0.1*stretch(discrz)))[2:end] )

# reinterpolate discrx from stretched back to the original z grid
xint = interpolate_ascending(discrz[2:end], discrx[2:end]).(z)
xint = interpolate_ascending((@. discrz * (1.0 + 0.1*stretch(discrz)))[2:end], discrx[2:end]).(z) # z stretched up
xint = interpolate_ascending((@. discrz * (1.0 - 0.1*stretch(discrz)))[2:end], discrx[2:end]).(z) # z stretched down

# clf()
plot(1e3*qmax*exp.(interpolate_ascending(discrz, discrx).(z)), z/1e3, "k", label="control")
plot(1e3*qmax*exp.(interpolate_ascending((@. discrz * (1.0 + 0.1*stretch(discrz))), discrx).(z)), z/1e3, label="up")
plot(1e3*qmax*exp.(interpolate_ascending((@. discrz * (1.0 - 0.1*stretch(discrz))), discrx).(z)), z/1e3, label="down")
plot(1e3*qmax*exp.(interpolate_ascending(discrz, discrx).(z)), z/1e3, "k")
ylim([0, 4]); xlim([0, 15]); 
xlabel("q (g/kg)"); ylabel("z (km)")
title("stretched humidity profiles")
[ savefig("stretch_q.$(f)") for f in ["png", "pdf", "eps", "svg"] ]

function stretch_humidity(qm, z; a=0.1)
    qmax = maximum(f0,qm)
    x = log.(filt_q(qm)./qmax)
    # discretize x = ln humidity for stretching
    discrx = 0:-0.05:-3.5
    discrz = interpolate_descending(x,z).(discrx)
    discrz[1] = z[findfirst(isfinite,qm)]
    # reinterpolate discrx from stretched back to the original z grid
    xint = interpolate_ascending((@. discrz * (1.0 + a*stretch(discrz))), discrx).(z) # z stretched up
    return qmax*exp.(xint)
end

clf()
plot((@. stretch(discrz)), discrz/1e3, color="tab:blue", label="up")
# plot((@. 1.0 - 0.7*stretch(discrz)), discrz/1e3, color="tab:orange", label="down")
plot(1e3*stretch_humidity(qm, z; a= 0.0), z/1e3, "k", label="control")
plot(1e3*stretch_humidity(qm, z; a=+0.1), z/1e3, color="tab:blue", label="up")
plot(1e3*stretch_humidity(qm, z; a=-0.1), z/1e3, color="tab:orange", label="down")
plot(1e3*stretch_humidity(qm, z; a= 0.0), z/1e3, "k")
ylim([0, 4]); xlim([0, 15]); 
xlabel("q (g/kg)"); ylabel("z (km)")
title("stretched humidity profiles")
[ savefig("stretch_q.$(f)") for f in ["png", "pdf", "eps", "svg"] ]

# A general way to make matching humidity and eddy flux profiles is to stretch/squash
# qm, ztop, G, G_i, F_i. Areas a_i stay the same. 
# Squashing increases and stretching decreasesthe sink rate α+ϵ.