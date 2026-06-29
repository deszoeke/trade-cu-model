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

includet("TradeCuExperiments.jl")
using .TradeCuExperiments

interp_cloudtop_height = TradeCuModel.interp_cloudtop_height # testing this

# experiment with sink rate fixed by the control
# ctx, ExpDict, controlsink, sinkm5, sinkp5 = TradeCuExperiments.test_control_sink() # not needed in REPL
ctx, ExpDict, controlsink, sinkm5, sinkp5 = test_control_sink();

# load standard parameters
( qm, qs, zcb, qcb, E_cb, x, divg, 
    tot_sink, cth_bin, rfv_acc, rfv_nrm, 
    rhoL, E_cb, qcb, ns, nz ) = setup_experiments(ctx=ctx)
# Get matched total sink rate and cloud top ensemble from control experiment
# aligned with the z grid for use in experiments.
# controlsink.output.acld # cloud fractions interpolated to ztop = z grid
function good_sinks(controlsink; zb=700.0, zt=4000.0)
    ii = findall(x-> zb<=x<=zt, ctx.z)
    (ctx.z[ii], controlsink.input.tot_sink[ii], 
     cumsum(controlsink.output.acld)[ii], controlsink.output.acld[ii])
end
ztop, sinkz, acc, acld = good_sinks(controlsink)

# DIM experiments with different sink rates, to simulate mesoscale organization changes
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
    tot_sink=reverse(sinkz), # sinkz increases
    cth_bin= ztop,
    rfv_acc= acc,
    rfv_nrm= acld );
DIMsinkm5 = define_experiment(; 
    name="DIMsink-5%", description="DIM sink rate -5%",
    qm= (1 .- 0.95*(1 .- qm./qs)).*qs*1.07, qs=1.07*qs, zcb=zcb,
    qcb=1.07*qcb, E_cb=1.02*E_cb, x=x, divg=0.95*divg,
    tot_sink=reverse(sinkz*0.95),
    cth_bin= ztop,
    rfv_acc= acc, 
    rfv_nrm= acld ); # how to specify the cloud fraction?
DIMsinkp5 = define_experiment(;
    name="DIMsink+5%", description="DIM sink rate +5%",
    qm= (1 .- 0.95*(1 .- qm./qs)).*qs*1.07, qs=1.07*qs, zcb=zcb,
    qcb=1.07*qcb, E_cb=1.02*E_cb, x=x, divg=0.95*divg,
    tot_sink=reverse(sinkz*1.05),
    cth_bin= ztop,
    rfv_acc= acc, 
    rfv_nrm= acld );
# integrate
for exp in [DIMsink, DIMsinkm5, DIMsinkp5]
    println(exp.name)
    integrate_experiment!(exp, ctx=ctx)
    push!(ExpDict, exp.name => exp)
end

"plot cloud vertical velocity for an experiment"
function plot_exp(e, var, ctx=ctx)
    ztop = interp_cloudtop_height(ctx.z, e.output.qc .- e.input.qs)
    iz = findall(500.0 .<= ctx.z .<= 3500.0) # 10 m bins, z=500 m
    ik = findall(x-> !ismissing(x) && isfinite(x), ztop) # tot_sink -> ztop
    q = getfield(e.output, var)
    pcolormesh(ztop[ik]/1e3, ctx.z[iz]/1e3, q[iz,ik], 
        cmap=ColorMap("BuPu")) #, vmin=0.0, vmax=5)
    colorbar()
    ylim([0.5, 3.5]); xlim([0.5, 3.5])
    xlabel("cloud top height (km)")
    ylabel("z coordinate (km)")
    title("$(e.name) $(string(var))")
    gca().set_aspect("equal")
end

"plot cloud vertical velocity for an experiment"
function plot_exp_w(e, ctx=ctx; kwargs...)
    ztop = interp_cloudtop_height(ctx.z, e.output.qc .- e.input.qs)
    iz = findall(500.0 .<= ctx.z .<= 3500.0) # 10 m bins, z=500 m
    ik = findall(x-> !ismissing(x) && isfinite(x), ztop) # tot_sink -> ztop
    log_norm = PythonPlot.matplotlib.colors.LogNorm(vmin=minimum(0.001), vmax=maximum(5))
    pcolormesh(ztop[ik]/1e3, ctx.z[iz]/1e3, e.output.w[iz,ik], 
        cmap=ColorMap("BuPu")) #, norm=log_norm, kwargs...)
    colorbar()
    ylim([0.5, 3.5]); xlim([0.5, 3.5])
    xlabel("cloud top height (km)")
    ylabel("z coordinate (km)")
    title("cloud vertical velocity (m/s)")
    gca().set_aspect("equal")
    tight_layout()
end

"plot cloud liquid for an experiment"
function plot_exp_qc(e, ctx=ctx; kwargs...)
    ztop = interp_cloudtop_height(ctx.z, e.output.qc .- e.input.qs) # recomputing here works
    ql = calc_ql(e)
    iz = findall(500.0 .<= ctx.z .<= 3500.0) # 10 m bins, z=500 m
    ik = findall(x-> !ismissing(x) && isfinite(x), ztop) # tot_sink -> ztop

    pcolormesh(ztop[ik]/1e3, ctx.z[iz]/1e3, ql[iz,ik]*1e3, 
        cmap=get_cmap("BuPu",10), vmin=0.0, vmax=2.0, kwargs...)
    colorbar()
    ylim([0.5, 3.5]); xlim([0.5, 3.5])
    xlabel("cloud top height (km)")
    ylabel("z coordinate (km)")
    title("cloud liquid (g/kg)")
    gca().set_aspect("equal")
    # tight_layout()
end

# plot cloud w
clf(); fig = figure()
fig.set_size_inches([10, 5])
subplot(2,3,1)
plot_exp_w(ExpDict["control"], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("control")
gca().set_aspect("equal")
subplot(2,3,2)
plot_exp_w(ExpDict["subsidence-5%"], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("subsidence -5%")
gca().set_aspect("equal")
subplot(2,3,3)
plot_exp_w(ExpDict["q&qs+7%"], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("q and qs +7%")
gca().set_aspect("equal")
suptitle("cloud vertical velocity (m/s)")
tight_layout()
# for f in ["png", "pdf", "svg"]
#     savefig("experiment_cloud_vertical_velocity.$f")
# end

# plot cloud qc liquid
# clf() 
subplot(2,3,4)
plot_exp_qc(ExpDict["control"], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("control")
gca().set_aspect("equal")
subplot(2,3,5)
plot_exp_qc(ExpDict["subsidence-5%"], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("(same clouds) subsidence -5%")
gca().set_aspect("equal")
subplot(2,3,6)
plot_exp_qc(ExpDict["q&qs+7%"], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title("q and qs +7%")
gca().set_aspect("equal")
suptitle("cloud vertical velocity (m/s), cloud liquid (g/kg)")
tight_layout()
# for f in ["png", "pdf", "svg"]
#     savefig("experiment_cloud_liquid.$f")
# end


# plot the experiments where the cloud top height change
# propagates to the flux partition.
# plot cloud w
expmts = ["control-sink", "DIMsink", "DIMsink-5%"]

clf(); fig = figure(); fig.set_size_inches([10, 5])
subplot(2,3,1)
plot_exp_w(ExpDict[expmts[1]], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title(expmts[1])
gca().set_aspect("equal")
subplot(2,3,2)
plot_exp_w(ExpDict[expmts[2]], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title(expmts[2])
gca().set_aspect("equal")
subplot(2,3,3)
plot_exp_w(ExpDict[expmts[3]], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title(expmts[3])
gca().set_aspect("equal")
suptitle("cloud vertical velocity (m/s)")
tight_layout()
# for f in ["png", "pdf", "svg"]
#     savefig("experiment_cloud_vertical_velocity.$f")
# end

# plot cloud qc liquid
# clf() 
subplot(2,3,4)
plot_exp_qc(ExpDict[expmts[1]], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title(expmts[1])
gca().set_aspect("equal")
subplot(2,3,5)
plot_exp_qc(ExpDict[expmts[2]], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title(expmts[2])
gca().set_aspect("equal")
subplot(2,3,6)
plot_exp_qc(ExpDict[expmts[3]], ctx)
plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
title(expmts[3])
gca().set_aspect("equal")
suptitle("cloud vertical velocity (m/s), cloud liquid (g/kg)")
tight_layout()
# for f in ["png", "pdf", "svg"]
#     savefig("experiment_cloud_liquid.$f")
# end




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
