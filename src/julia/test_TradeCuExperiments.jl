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

interp_cloudtop_height = TradeCuModel.interp_cloudtop_height # needed for plotting output

# experiment with sink rate fixed by the control
ctx, ExpDict, controlsink, sinkm5, sinkp5 = test_control_sink();

# load standard parameters
( qm, qs, zcb, qcb, E_cb, x, divg, sfc_adv,
    tot_sink, cth_bin, rfv_acc, rfv_nrm, 
    rhoL, ns, nz ) = setup_experiments(ctx=ctx);
# Get matched total sink rate and cloud top ensemble from control experiment
# aligned with the z grid for use in experiments.
# controlsink.output.acld; cloud fractions interpolated to ztop = z grid
function good_sinks(controlsink; zb=700.0, zt=4000.0)
    ii     = findall(x-> zb<=x<=zt, ctx.z)       # z-grid indices for 700–4000 m
    ii_cth = findall(x-> zb<=x<=zt, ctx.cth_bin)  # GOES CTH indices for 700–4000 m
    acld  = controlsink.output.acld               # model cloud fractions (331 elements)
    sinkz = controlsink.input.tot_sink            # 331 elements
    acc   = ctx.rfv_acc[ii_cth]                   # GOES survival function (decreasing)
    (ctx.z[ii], sinkz, acc, acld)
end
ztop, sinkz, acc, acld = good_sinks(controlsink);

# deprecate control-sink and DIMsink experiments.

# DIM experiments with different sink rates, to simulate mesoscale organization changes
DIMsink = define_experiment(controlsink; 
    name="DIMsink", description="DIM sink rate, cloud tops align with z grid",
    qm=1.07*qm, qs=1.07*qs, qcb=1.07*qcb,
    E_cb=1.02*E_cb,
    divg=0.95*divg, sfc_adv=0.95*sfc_adv,
    tot_sink=sinkz,               # sinkz decreases with index
    cth_bin= ztop, rfv_nrm= acld, # probably ignored
    control=false, a_i_control=controlsink.output.acld, M_i_control=controlsink.output.M );
DIMsinkm5 = define_experiment(DIMsink; name="DIMsink-5%", description="DIM sink rate -5%",
    tot_sink=sinkz*0.95 );
DIMsinkp5 = define_experiment(DIMsink; name="DIMsink+5%", description="DIM sink rate +5%",
    tot_sink=sinkz*1.05 );

# integrate
for exp in [DIMsink, DIMsinkm5, DIMsinkp5]
    println(exp.name)
    integrate_experiment!(exp, ctx=ctx)
    push!(ExpDict, exp.name => exp)
end

# recompute the cloud fraction for the DIM experiments, assuming
# w_i is invariant from the control experiment, stretching w profile from old ztop to new ztop.
# a_i = new_area( ExpDict["DIMsink"], ExpDict["control-sink"], ctx )
# possibly only cloud BASE mass flux is a good indicator of cloud fraction.

# quick count of valid fields
good(x) = !ismissing(x) && isfinite(x)
function inventory(e)
    for f in fieldnames(typeof(e))
        v = getfield(e, f)
        if eltype(v) <: Union{Missing, T} where T <: Number
            println("$(string(f)) has $(sum(good,v))")
        else
            println("$(string(f)) is nonnumeric.")
        end
    end
end
# inventory(ExpDict["control-sink"].input)
# inventory(ExpDict["control-sink"].output)

# analyze total cloud fraction from experiments
totcld(e) =  sum(filter(isfinite,skipmissing(e.output.acld)))
# totcld(ExpDict["control-sink"]), totcld(ExpDict["DIMsink"]), totcld(ExpDict["DIMsink-5%"]), totcld(ExpDict["DIMsink+5%"])
dlna(e,c) = log(totcld(e)/totcld(c))
function dlna_samelevs(e,c) 
    ii = good.(e.output.acld) .& good.(c.output.acld)
    log(sum(e.output.acld[ii])/sum(c.output.acld[ii]))
end
# limiting to same sink rate has no effect on E1 experiments
function dlna_limit_index(e,c) 
    ik = findall(x-> !ismissing(x) && isfinite(x), c.output.ztop)
    totcld(e) =  sum(filter(isfinite,skipmissing(e.output.acld[ik])))
    log(totcld(e)/totcld(c))
end
function dlna_limit_ztop(e,c) # do not go above the highest control cloud top height
    ik = findall(x-> !ismissing(x) && isfinite(x), c.output.ztop)
    mz = maximum(c.output.ztop[ik])
    ij = findall(x-> !ismissing(x) && isfinite(x) && x<=mz, e.output.ztop)
    # println("max ztop for control = $(mz/1e3) km, max ztop for experiment = $(maximum(e.output.ztop[ij])/1e3) km")
    log(sum(e.output.acld[ij])/sum(c.output.acld[ik]))
end

"""
dlna_itp_ztop(e,c)   interpolate at max control.output.ztop
compares the same valid cloud levels as for the control experiment
"""
function dlna_itp_ztop(e,c)
    ik = findall(x-> !ismissing(x) && isfinite(x), c.output.ztop)
    mz = maximum(c.output.ztop[ik])
    ca = sum(c.output.acld[ik])
    ij = findall(x-> !ismissing(x) && isfinite(x), e.output.ztop)
    ea = TradeCuModel.interpolate_ascending(
        coalesce.(e.output.ztop[ij],NaN), 
        cumsum(coalesce.(e.output.acld[ij],NaN)) )(mz)
    log(ea/ca)
end

# print readable table for E1 experiments
open("cloud_frac_table.txt", "w") do io
    redirect_stdout(io) do
        println("cloud fraction, % change from control")
        println(@sprintf("%-15s | %10s", "experiment", "dlna")) # interpolating to ztop
        println("-"^30)
        for exp in ["subsidence-5%", "Ecb+2%", "q&qs+7%", "lclRH+0.003", "lclRH+0.006"]
            println(@sprintf("%-15s | %10.2f", exp, 
                100*dlna_itp_ztop(   ExpDict[exp], ExpDict["control"]) ) )
        end
    end
end
# Inverse dependence of mass flux on Δq, M = G/(Δq) 
# increases Δq and decreases M = a*w and thus cloud fraction.

# mesoscale organization experiments with different sink rates
# cloud tops lie on the same curve when plotted vs sink rate
plot(ExpDict["control"].input.tot_sink*1e3, ExpDict["control"].output.ztop/1e3)
plot(ExpDict["sink+5%"].input.tot_sink*1e3, ExpDict["sink+5%"].output.ztop/1e3)
plot(ExpDict["sink-5%"].input.tot_sink*1e3, ExpDict["sink-5%"].output.ztop/1e3)
clf() # plot vs. category subscripts - shifts same height clouds to different sink rates
plot( ExpDict["control"].output.ztop/1e3)
plot( ExpDict["sink+5%"].output.ztop/1e3)
plot( ExpDict["sink-5%"].output.ztop/1e3)

Experiment = TradeCuModel.Experiment
function plot_exp_vs_control(c::Experiment, e1::Experiment, e2::Experiment, var::Symbol; f=identity, kwargs...)
    cvar = f(getfield(c.output, var))
    e1var = f(getfield(e1.output, var))
    e2var = f(getfield(e2.output, var))

    clf() # plot experiment vs control cloud top heights
    plot( cvar, e1var, label=e1.name; kwargs...)
    plot( cvar, e2var, label=e2.name; kwargs...)
    xlabel("control"); ylabel("experiment")
    legend(frameon=false)
    return gcf()
end

plot_exp_vs_control(ExpDict["control"], ExpDict["sink-5%"], ExpDict["sink+5%"], 
    :ztop; f=(x->x*1e-3))
plot([0.7, 3.5], [0.7, 3.5], "k-", linewidth=0.5, label="1:1")
ylim([0.7, 3.2]); xlim([0.7, 3.2])

plot_exp_vs_control(ExpDict["control"], ExpDict["sink-5%"], ExpDict["sink+5%"], 
    :M; f=(M->M[71,:]) ) # ack

plot_exp_vs_control(ExpDict["control"], ExpDict["sink-5%"], ExpDict["sink+5%"], 
    :acld; linestyle="none", marker=".", markersize=1 )

# print readable table for E2 sink rate experiments
begin
    println("cloud fraction, % change from control")
    println(@sprintf("%-15s | %10s", "experiment", "dlna")) # interpolating to ztop
    println("-"^30)
    for exp in ["DIMsink", "DIMsink-5%", "DIMsink+5%"]
        println(@sprintf("%-15s | %10.2f", exp, 
            100*dlna(   ExpDict[exp], ExpDict["control-sink"]) ) )
    end
end
# interpolating to the original ztop gives very different results
# Clouds with weaker sink rates rise above the original ztop.
# Stronger sink rates give missing clouds at original ztop.

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


"plot effect of experiments on cloud fraction profiles"
function plot_cld_tot_profs(expmts)
    ax = fig.add_subplot(1,1, 1)
    for (i, exp) in enumerate(expmts)
        ztop = ExpDict[exp].output.ztop
        ik = findall(x-> !ismissing(x) && isfinite(x), ztop)
        ac = cumsum(ExpDict[exp].output.acld[ik])
        ax.plot(ac*1e2, ztop[ik]/1e3, 
            linewidth=0.5, marker=".", markersize=2,
            label="$(exp), a=$(round(maximum(ac*100), digits=1))%")
    end
    ax.set_xlabel("cumulative cloud fraction (%)")
    ax.set_ylabel("cloud top height (km)")
    ax.set_xlim([0, 10]); ax.set_ylim([0.5, 3.5])
    ax.legend(frameon=false)
    return ax
end

expmts = ["control", "subsidence-5%", "q&qs+7%", "lclRH+0.00875", "lclRH+0.0175"]
fig = gcf(); fig.clf()
ax = plot_cld_tot_profs(expmts)

expmts = ["control-sink", "DIMsink", "DIMsink-5%", "DIMsink+5%"]
fig = gcf(); fig.clf()
ax = plot_cld_tot_profs(expmts)

# plotting experiment structure unpackers and calculators
unpack(v,e) = foldl(getfield, v; init=e)
multop(f, v, T...) = f.(broadcast(x->unpack(v,x), T)...)
# computes dlnM output between experiments e-c
f(x,y) = log(x/y)
dlnM = multop(f, [:output, :M], ExpDict["DIMsink"], ExpDict["control-sink"]);

# plotting functions for experiment output
Experiment = TradeCuModel.Experiment # for dispatch

"pcolormesh q::Matrix for experiment e"
function plot_exp_var(e::Experiment, q::Matrix, ctx=ctx; ncolor=10, ax=nothing, kwargs...)
    ax = isnothing(ax) ? gca() : ax
    ztop = e.output.ztop
    iz = findall(500.0 .<= ctx.z .<= 3500.0) # 10 m bins, z=500 m
    ik = findall(x-> !ismissing(x) && isfinite(x), ztop)
    pcm = ax.pcolormesh(ztop[ik]/1e3, ctx.z[iz]/1e3, q[iz,ik]; cmap=get_cmap("BuPu", ncolor), kwargs...)
    ax.get_figure().colorbar(pcm, ax=ax)
    ax.plot(ztop[ik]/1e3, e.output.acld[ik], "w-", linewidth=0.5, label="cloud_i fraction")
    ax.plot([0.5, 3.5], [0.5, 3.5], "k-", linewidth=0.5)
    ax.set_ylim([0.5, 3.5]); ax.set_xlim([0.5, 3.5]); ax.set_yticks([1.0, 2.0, 3.0],["1", "2", "3"])
    ax.set_aspect("equal")
end
"pcolormesh an experiment output variable var"
plot_exp_var(e::Experiment, var::Symbol, ctx=ctx; ax=nothing, kwargs...) = plot_exp_var(e::Experiment, getfield(e.output, var), ctx; ax=ax, kwargs...)
plot_exp_var(f, e::Experiment, var::Symbol, ctx=ctx; ax=nothing, kwargs...) = plot_exp_var(e::Experiment, f.(getfield(e.output, var)), ctx; ax=ax, kwargs...)

# plot Δq = qc-qm for surface (1-RH)-5% experiment
dq_cld(e) = e.output.qc .- e.input.qm
# expmts = ["control", "subsidence-5%", "q&qs+7%", "sfc(1-RH)-5%"]
exp = "sfc(1-RH)-5%"
e = ExpDict[exp]
c = ExpDict["control"]
fig = gcf(); # fig.set_size_inches([5, 5]); 
fig.clf()
ax = fig.subplots(1, 1)
qx = (100*log.(max.(0, dq_cld(e)./dq_cld(c)))); # % change in Δq relative to control
qx[ctx.z .< ctx.zcb, :] .= NaN;
for (i, zt) in enumerate(e.output.ztop); qx[ctx.z.>=zt.+10.0, i] .= NaN; end
qx[qx.<1.0] .= NaN; qx[qx.>4.25] .= NaN;
plot_exp_var(e, qx, ctx; ax=ax, cmap=get_cmap("BuPu", 13), vmin=1.0, vmax=4.25)
ax.plot([0.7, 3.5], [0.7, 3.5], "k-", linewidth=0.5)
ax.plot([0.7, 3.5], [0.7, 0.7], "k-", linewidth=0.5)
ax.set_title("\$d\$ln\$(q_c-q)\$ [%] for $(exp) - control")
ax.set_ylabel("z coordinate (km)")
ax.set_xlabel("cloud top height (km)")
fig.tight_layout()
[ fig.savefig("experiment_dq.$f") for f in ["png", "pdf", "svg"] ];

# plot cloud w
expmts = ["control-sink", "DIMsink", "DIMsink-5%"]
fig = gcf(); fig.set_size_inches([10, 5]); fig.clf()
axs = fig.subplots(2, 3)
for (i, exp) in enumerate(expmts)
    ax = axs[0, i-1]
    plot_exp_var(ExpDict[exp], :w, ctx; vmax=1.0, ax=ax)
    ax.set_title(exp)
    (i-1)%3 == 0 && ax.set_ylabel("in-cloud vertical velocity (m/s)\n\nz coordinate (km)", size=12)
end
for (i, exp) in enumerate(expmts)
    ax = axs[1, i-1]
    plot_exp_var(ExpDict[exp], 1e3*calc_ql(ExpDict[exp]), ctx; vmin=0.0, vmax=2.0, ax=ax)
    # ax.set_title(exp)
    (i-1)%3 == 0 && ax.set_ylabel("cloud liquid water (g/kg)\n\nz coordinate (km)", size=12)
    ax.set_xlabel("cloud top height (km)")
end
fig.tight_layout()
[ fig.savefig("experiment_cloud_vel_liquid.$f") for f in ["png", "pdf", "svg"] ]


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
