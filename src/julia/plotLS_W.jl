using Revise
using Printf
using Statistics

using PythonPlot
using PythonCall
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

includet("TradeCuExperiments.jl")
using .TradeCuExperiments
interpolate_ascending = TradeCuModel.interpolate_ascending

# get data (runs model too)
ctx, ExpDict, controlsink, sinkm5, sinkp5 = test_control_sink();

# load standard parameters
( qm, qs, zcb, qcb, E_cb, x, divg, sfc_adv,
    tot_sink, cth_bin, rfv_acc, rfv_nrm, 
    rhoL, ns, nz ) = setup_experiments(ctx=ctx);
divg = ExpDict["control"].input.divg
sfc_adv = ExpDict["control"].input.sfc_adv
W = TradeCuModel.subsidence.(ctx.z; divg=divg, zi=ctx.zi)

# large-scale subsidence vertical velocity (m/s)
subsidence(z_; divg=divg, zi=zi) = -min(z_, zi) * divg
"interpolate the derivative dq/dz at z_"
ddz(q=qm, z=ctx.z) = interpolate_ascending( z[1:end-1].+0.5*diff(z), diff(q) ./ diff(z) ) # returns a function of z_
dqdz(z) = ddz(qm, ctx.z)(z)
"""
large-scale drying profile S_ls(z) = -w*dq/dz - 1.7e-8*(zdivg-z)/zdivg
kg/kg s^-1
"""
function largescale_drying(z_; q=qm, z=ctx.z, divg=divg, sfc_adv=sfc_adv, zdivg=4e3)
    # zdivg = 4e3 # m; top height
    dqdz(z_) = ddz(q, z)(z_)
    wdqdz = subsidence(z_; divg=divg, zi=zdivg) * dqdz(z_)
    -wdqdz - sfc_adv * max(0, (zdivg - z_) ./ zdivg) # advection increases to sfc_adv at surface
end

# plot the subsidence W
clf()
plot([0, 0], [0, 5], "k-")
plot( 1e3*W, ctx.z/1e3, label="W" )
plot( 1e3*86400*largescale_drying.(ctx.z), ctx.z/1e3, label="S_ls" )
ylim([0, 5])
ylabel("height (km)")
xlabel("large-scale subsidence W (mm/s)\nlarge-scale drying (g/kg d\$^{-1}\$)")
title("control")
tight_layout()

[ savefig("subs_drying.$f") for f in ["png", "pdf", "svg"] ]
