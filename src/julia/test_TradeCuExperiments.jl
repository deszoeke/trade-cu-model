# using REPL
# REPL.activate(TradeCuExperiments) # steps in to module scope
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

PythonPlot.matplotlib.rcParams["font.family"] = "sans-serif"
PythonPlot.matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "OpenSans"]

include("TradeCuExperiments.jl")
using .TradeCuExperiments


# experiment with sink rate fixed by the control
# ctx, ExpDict, controlsink, sinkm5, sinkp5 = TradeCuExperiments.test_control_sink() # not needed in REPL
ctx, ExpDict, controlsink, sinkm5, sinkp5 = test_control_sink()

# load standard parameters
( qm, qs, zcb, qcb, E_cb, x, divg, 
    tot_sink, cth_bin, rfv_acc, rfv_nrm, 
    rhoL, E_cb, qcb, ns, nz ) = setup_experiments(ctx=ctx)
sinkz = controlsink.input.tot_sink
controlsink.output.acld
ztop = ctx.z/1e3
acc = cumsum(controlsink.output.acld)

# mesoscale DIM experiments
DIMsink = define_experiment(; 
    name="DIMsink",
    description="DIM sink rate, cloud tops align with z grid",
    qm= (1 .- 0.95*(1 .- qm./qs)).*qs*1.07,
    qs=1.07*qs,
    zcb=zcb,
    qcb=1.07*qcb,
    E_cb=1.02*E_cb,
    x=x,
    divg=0.95*divg,
    tot_sink=sinkz,
    cth_bin=ztop,
    rfv_acc=acc,
    rfv_nrm=controlsink.output.acld );
DIMsinkm5 = define_experiment(; 
    name="DIMsink-5%", description="DIM sink rate -5%",
    qm= (1 .- 0.95*(1 .- qm./qs)).*qs*1.07, qs=1.07*qs, zcb=zcb,
    qcb=1.07*qcb, E_cb=1.02*E_cb, x=x, divg=0.95*divg,
    tot_sink=sinkz*0.95,
    cth_bin=ztop, rfv_acc=acc, rfv_nrm=controlsink.output.acld ); # how to specify the cloud fraction?
DIMsinkp5 = define_experiment(; 
    name="DIMsink+5%", description="DIM sink rate +5%",
    qm= (1 .- 0.95*(1 .- qm./qs)).*qs*1.07, qs=1.07*qs, zcb=zcb,
    qcb=1.07*qcb, E_cb=1.02*E_cb, x=x, divg=0.95*divg,
    tot_sink=sinkz*1.05,
    cth_bin=ztop, rfv_acc=acc, rfv_nrm=controlsink.output.acld );
# integrate
for exp in [DIMsink, DIMsinkm5, DIMsinkp5]
    println(exp.name)
    integrate_experiment!(exp, ctx=ctx)
    push!(ExpDict, exp.name => exp)
end


"plot cloud vertical velocity for an experiment"
function plot_exp_w(e, ctx=ctx)
    pcolormesh(ctx.z[51:350]/1e3,ctx.z[51:350]/1e3, e.output.w[51:350,51:350], 
        cmap=ColorMap("BuPu")) #, vmin=0.0, vmax=5)
    colorbar()
    ylim([0.5, 3.5]); xlim([0.5, 3.5])
    xlabel("cloud top height (km)")
    ylabel("z coordinate (km)")
    title("cloud vertical velocity (m/s)")
    gca().set_aspect("equal")
    tight_layout()
end

"plot cloud liquid for an experiment"
function plot_exp_qc(e, ctx=ctx)
    ql = calc_ql(e)
    pcolormesh(ctx.z[51:350]/1e3,ctx.z[51:350]/1e3, ql[51:350,51:350]*1e3, 
        cmap=get_cmap("BuPu",10), vmin=0.0, vmax=2.0)
    colorbar()
    ylim([0.5, 3.5]); xlim([0.5, 3.5])
    xlabel("cloud top height (km)")
    ylabel("z coordinate (km)")
    title("cloud liquid (g/kg)")
    gca().set_aspect("equal")
    tight_layout()
end

# plot cloud w
clf() 
subplot(2,3,1)
plot_exp_w(controlsink, ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("control")
gca().set_aspect("equal")
subplot(2,3,2)
plot_exp_w(sinkm5, ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("sink rate -5%")
gca().set_aspect("equal")
subplot(2,3,3)
plot_exp_w(sinkp5, ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("sink rate +5%")
gca().set_aspect("equal")
suptitle("cloud vertical velocity (m/s)")
tight_layout()

for f in ["png", "pdf", "svg"]
    savefig("experiment_cloud_vertical_velocity.$f")
end

# plot cloud qc liquid
# clf() 
subplot(2,3,4)
plot_exp_qc(controlsink)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("control sink rate")
gca().set_aspect("equal")
subplot(2,3,5)
plot_exp_qc(sinkm5)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("sink rate -5%")
gca().set_aspect("equal")
subplot(2,3,6)
plot_exp_qc(sinkp5)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("sink rate +5%")
gca().set_aspect("equal")
suptitle("cloud vertical velocity (m/s), cloud liquid (g/kg)")
tight_layout()
for f in ["png", "pdf", "svg"]
    savefig("experiment_cloud_liquid.$f")
end

"""
power law on log-log plot
y = y0 * (x/x0)^p
"""
plotl(x, y0, p) = loglog(x, y0 * (x./x[1]).^p) 
plotlz(z, x0, p) = loglog(x0 * (z./z[1]).^p, z) 

clf()
loglog(sinkz, a_i, linestyle="none", marker=".")
plotl([6.55e-4, 6e-4], 1.5e-4, 10)
xlim([6e-4, 1e-3]); ylim([1e-4, 2e-4])
# so at low tot_sink (hi clouds), a = 1.5e-4 * (tot_sink/6.55e-4)^10

clf()
loglog( a_i, ctx.z/1e3, marker=".")
ylim([2.5,4.1])
plotlz([2.5, 4.1], 3.1e-4, -5.5)
# above 3 km, extrapolate that cloud fraction decreases as 
# a = (3e-4/10 m) (z/2.5 km)^-5.5

# graph shifting the sink rate by ±5%
clf()
loglog(     sinkz, a_i)
loglog(0.95*sinkz, a_i, ":")
loglog(1.05*sinkz, a_i, ":")
# Cloud fraction is so small at the low-sink end that no extrapolation is recommended.

# plot sink rate parameter ensembles for mesoscale experiments
clf()
subplot(2,2,1)
plot(1e3*controlsink.input.tot_sink, ctx.z/1e3, marker=".", linestyle="none", label="control")
plot(1e3*sinkm5.input.tot_sink, ctx.z/1e3, marker=".", linestyle="none", label="total sink -5%")
plot(1e3*sinkp5.input.tot_sink, ctx.z/1e3, marker=".", linestyle="none", label="total sink +5%")
plot(1e4*a_i[1:401], ctx.z[1:401]/1e3, "k", label=L"a_i")
legend( frameon=false )
ylabel("cloud top height (km)")
xlabel("total sink rate (km\$^{-1}\$)")
xlim([0, 6]); ylim([0, 3.5])
title("cloud top height vs. sink rate")

subplot(2,2,2)
plot(1e3*controlsink.input.tot_sink, 1e3*a_i, label="control")
plot(1e3*sinkm5.input.tot_sink, 1e3*a_i, label="total sink -5%")
plot(1e3*sinkp5.input.tot_sink, 1e3*a_i, label="total sink +5%")
title("cloud fraction vs. sink rate")
xlabel("total sink rate (km\$^{-1}\$)")
ylabel("cloud fraction (10\$^{-3}\$)")
xlim([0, 6])
tight_layout()
suptitle("sink rates for mesoscale experiments")
