# using REPL
# REPL.activate(TradeCuExperiments) # steps in to module

# functions added to the module
# function define_control_sink_experiment(; ctx::ModelContext, sinkz=sinkz)
#     ( qm, qs, zcb, qcb, E_cb, x, divg, 
#         tot_sink, cth_bin, rfv_acc, rfv_nrm, 
#         rhoL, E_cb, qcb, ns, nz ) = TradeCuExperiments.setup_experiments(ctx=ctx)

#     return Experiment(
#             "control-sink", "Control for setting compact set of sink rates",
#             ModelInput(qm, qs, zcb, qcb, E_cb, x, divg, sinkz, cth_bin, rfv_acc, rfv_nrm),
#             TradeCuExperiments.allocate_output(nz,ns) )
# end

# function cloud_i_area( ctx )
#     z=ctx.z
#     nz = length(ctx.z)
#     cth_bin=ctx.cth_bin
#     rfv_acc=ctx.rfv_acc
#     # align z grid with fractions
#     offset = findfirst(x->x≈cth_bin[1]*1e3, z) - 1 # offset for 10 m bins
#     a_i = zeros(nz)
#     # cloud with qc=0 at h_i diverges below between h_(i-1) and h_i
#     a_i[offset+1 .+ eachindex(rfv_acc[1:end-1])] .= -diff(rfv_acc)
#     return a_i
# end

using Revise
using Debugger

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

include("TradeCuExperiments.jl")
using .TradeCuExperiments

# run all the experiments and fill the output structures
ctx = init_context()
ExpDict = define_experiments(ctx=ctx)
for exp in values(ExpDict) # run all the defined experiments
    println(exp.name)
    integrate_experiment!(exp, ctx=ctx)
end
# text output
println("Experiment names: $(keys(ExpDict))")
println("ctx: $(fieldnames(typeof(ctx)))")
control = ExpDict["control"]
println("control.input: $(fieldnames(typeof(control.input)))")
println("control.output: $(fieldnames(typeof(control.output)))")

# experiment hierarchy
keyorder = [
  "control", 
  "subsidence-5%",
  "qs+7%",
  "q&qs+7%",
  "Ecb+2%",
  "(1-RH)-5%",
  "DIM" ]

#=
# plot figures of results
close("all")
figure()
colrs = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
for (i,k) in enumerate(keyorder)
    exp = ExpDict[k]
    plot( 1e3*exp.input.qm, 1e-3*ctx.z, linewidth=0.5, label=k, color=colrs[i] )
    plot( 1e3*exp.input.qs, 1e-3*ctx.z, linewidth=0.5, color=colrs[i] )
end
legend(frameon=false)
ylim([0, 4])
xlabel("specific humidity (g/kg)")
ylabel("z (km)")
display(gcf())

# contour the cloud profiles as a function of sink rate
PythonPlot.matplotlib.pyplot.close("all")
figure()
all_ql = vcat([vec(1e3 .* coalesce.(max.(0, ExpDict[k].output.qc .- ExpDict[k].input.qs), NaN)) for k in keyorder]...)
finite_vals = filter(isfinite, all_ql)
vmax = isempty(finite_vals) ? 1.0 : maximum(finite_vals)
norm = PythonPlot.matplotlib.colors.Normalize(vmin=0.0, vmax=vmax)
for (i,k) in enumerate(keyorder[[1, 2, 3, 4, 5, 6]])
    exp = ExpDict[k]
    ql = max.(0, exp.output.qc .- exp.input.qs)
    ql_plot = 1e3 .* coalesce.(ql, NaN)
    contour(exp.input.tot_sink, 1e-3*ctx.z, ql_plot,
        colors=colrs[i], norm=norm, vmin=0.0, vmax=vmax)
    plot([-1, -2], [-1, -2], color=colrs[i], label=k)
end
# colorbar()
ylim([0, 4])
xlim([minimum(control.input.tot_sink), maximum(control.input.tot_sink)])
xlabel("sink rate (km\$^{-1}\$)")
ylabel("z (km)")
legend(frameon=false)
display(gcf())
=#

# The clouds don't depend on the fluxes at all.
# Another way to experiment is to keep the control distribution of
# sink rates, and recompute the clouds

# So we can integrate the clouds separately,
# find the sink rates that give clouds of each height,
# assign fractions to these from the obseved 
# cloud top height and 
# partition the large-scale moisture flux to the cloud flux.
# Then finally calculate the cloud fluxes, mass fluxes, and velocities.

# get sink rate as a function of cloud top height
sinkz = get_sinkrate( ExpDict["control"]; ctx=ctx )
a_i = cloud_i_area( ctx )
# a_i, sinkz aligned with cloud top heights z
# plot(sinkz, ctx.z, "k"); display(gcf())

# Now sinkz and ctx.z are 1 to 1, and 75:302 are valid low clouds.
# Map these to the cloud top height distribution 
# to get the cloud fraction for each of the interpolated sink rates in sinkz.

# test the compact sink rates
# insert sinkz in as the total sink rate
controlsink = define_control_sink_experiment(ctx=ctx, sinkz=sinkz)
# this fails:
integrate_experiment!(controlsink, ctx=ctx) # doesn't work with sinkz

# function integrate_experiment_sinkz!(exp::Experiment; ctx::ModelContext)
#     println("integrate_experiment! $(exp.name)")
#     z = ctx.z
#     dz = ctx.dz

#     # compute cloud+precipitation moisture flux that balances large scale drying
#     # ! independent of clouds !
#     F2z, G_ls = calc_Ftot( ctx=ctx, E_cb=exp.input.E_cb, divg=exp.input.divg,
#         qm=exp.input.qm, 
#         cth_bin=exp.input.cth_bin, cth_acc=exp.input.cth_acc, 
#         dz=dz, icb=findfirst(z .>= exp.input.zcb) )
#     # Flux is distributed to the clouds by their cloud top height distribution.

#     # run the cloud model for many sink rates
#     # ! qcld is independent of fluxes !
#     zt, F_cld, F_pcp, qcld = cloudflux_1x(
#         exp.input.tot_sink; x=exp.input.x, 
#         z=z, nz=length(z), 
#         dz=z[2]-z[1],
#         qm=exp.input.qm, qs=exp.input.qs, 
#         F2z=F2z, icb=findfirst(z .>= exp.input.zcb), 
#         qcb=exp.input.qcb )

#     # get sink rate as a function of cloud top height
#     # postprocess to get w, a, M, ...
#     w, _ = updraft_w_dq(F_cld, qcld, exp.input.qm, z, zt)
#     # print("size(w): $(size(w))") # nz, ns = (3100, 600)

#     # Interpolate satellite coordinate to model sinkrate coordiante.
#     # cloud fraction density per unit sink rate
#     # da/dsinkrate = da/dh * dh/dsinkrate.

#     da_dsink, da_ind = dadsinkrate(zt, exp.input.tot_sink, exp.input.cth_bin, exp.input.cth_nrm)
#     acld = ctx.dsink * da_dsink # cloud area fraction in sink rate bin

#     println("size(w)=$(size(w))")
#     exp.output.w .= w
#     exp.output.acld[da_ind] .= acld
#     exp.output.M[:,da_ind] .= w[:,da_ind] .* acld' # mass flux per sink rate
#     exp.output.qc .= qcld
#     exp.output.F_cld .= F_cld
#     exp.output.F_pcp .= F_pcp
#     exp.output.G_cld[:,da_ind] .= F_cld[:,da_ind] .* acld'
#     exp.output.G_pcp[:,da_ind] .= F_pcp[:,da_ind] .* acld'

#     return nothing
# end
# let
#     integrate_experiment_sinkz!(control_sink; ctx=ctx)
# end

# area fractions, not testing TradeCuExperiments
# maximum(ctx.rfv_acc) - ctx.rfv_acc[findlast(ctx.cth_bin .< 1.2)] # 1.4 %
# maximum(ctx.rfv_acc) - ctx.rfv_acc[findlast(ctx.cth_bin .< 1.0)] # 1.9 %
