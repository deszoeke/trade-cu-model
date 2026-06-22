# environment and code loading header
using Revise
cd(joinpath(homedir(), "Projects/ATOMIC/trade-cu-model/src/julia"))
using Pkg; Pkg.activate(".")
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

"Trade Cumulus Model Experiments Module"
module TradeCuExperiments

using ..TradeCuModel
using Statistics
using VaporSat # dev ../../deps/VaporSat

export ModelInput, ModelOutput, ModelContext, Experiment # types for experiment data
export ExpDict # dictionary contains defined experiments
export init_context, define_experiments
export integrate_experiment!

# Define the Inputs Container
struct ModelInput
    qm::Vector
    qs::Vector
    zcb::Number
    qcb::Number
    E_cb::Number
    x::Number
    divg::Number
    tot_sink::Vector
    cth_bin::Vector
    cth_acc::Vector # accumulated cloud fraction below cth_bin
    cth_nrm::Vector # cloud fraction within cth_bin
end

# Define the Outputs Container
struct ModelOutput
    M::Matrix{Union{U, Missing}} where U<:Number
    w::Matrix{Union{U, Missing}} where U<:Number
    acld::Vector{Union{U, Missing}} where U<:Number
    qc::Matrix{Union{U, Missing}} where U<:Number
    F_cld::Matrix{Union{U, Missing}} where U<:Number
    F_pcp::Matrix{Union{U, Missing}} where U<:Number
    G_cld::Matrix{Union{U, Missing}} where U<:Number
    G_pcp::Matrix{Union{U, Missing}} where U<:Number
end

# Define Experiments of input-output pairs
struct Experiment
    name::String
    description::String
    input::ModelInput
    output::ModelOutput
end

"Shared, cached data loaded once for all experiments."
struct ModelContext
    z::Vector{Float64}
    qm::Vector{Float64}
    qs::Vector{Float64}
    cth_bin::Vector{Float64}
    rfv_acc::Vector{Float64}
    rfv_nrm::Vector{Float64}
    dz::Float64
    dsink::Float64
    zcb::Float64
    zi::Float64
    ztop::Float64
    rhoL::Float64
end

"Load sounding and GOES data once and return initialized context."
function init_context()
    KelvinCelsius = 273.15 # K
    zi = 4.0e3             # m
    ztop = 4.0e3           # m
    zcb = 700.0            # m

    tot_sink = range(6.3523e-4, 5.7e-3, length=600)
    dsink = tot_sink[2] - tot_sink[1]

    z, tam, _, qm, pm = get_mean_soundings()
    qs = qsat.(pm, tam .- KelvinCelsius)
    tvm = virtual_temp.(tam, qm)
    dz = z[2] - z[1]

    rfv_nrm, rfv_acc, cth_bin = get_goes_cloud_data()
    rhoL = mean(filter(isfinite, calc_rhoL.(tvm, pm)[z .<= ztop]))

    return ModelContext(
        collect(z),
        collect(qm),
        collect(qs),
        vec(collect(cth_bin)),
        vec(collect(rfv_acc)),
        vec(collect(rfv_nrm)),
        dz,
        dsink,
        zcb,
        zi,
        ztop,
        rhoL,
    )
end


"initialize parameters, grids, data for experiments"
function setup_experiments(; ctx::ModelContext)
    #   parameters and initialization
    zi   = ctx.zi    # m  # inversion for tapering subsidence
    ztop = ctx.ztop  # m  # top of cloud model integration domain
    zcb  = ctx.zcb   # m
    divg = 1.5e-6   # 1/s
    x = 0.53 # parameter precipitation efficiency

    # ensemble of sink rates
    tot_sink = range(6.3523e-4, 5.7e-3, length=600) # min tuned for x=0.53 to get the highest possible cloud top
    dsink = ctx.dsink
    # tot_sink = (1 .+tanh.(range(-8*pi, 0, length=600))) .* (5e-3 - 1e-4) .+ 1e-4
    # tot_sink = range(6.1716e-4, 5.8e-3, length=600) # min tuned for x=0.53 to get the highest possible cloud top

    z = ctx.z
    qm = ctx.qm
    qs = ctx.qs
    dz = ctx.dz
    # icb = findfirst(z .>= zcb) # cloud base index

    cth_bin = ctx.cth_bin
    rfv_acc = ctx.rfv_acc
    rfv_nrm = ctx.rfv_nrm
    #size(rfv_nrm), size(cth_bin) # (351,1)

    rhoL = ctx.rhoL # 2.41e6 J/m^3
    E_cb = 180.0 # W/m^2; just the cloud vapor flux; E0 - 35?
    qcb = qs[findfirst(z .>= zcb)] # kg/kg; cloud base specific humidity
    ns = length(tot_sink)
    nz = length(z)
    return ( qm, qs, zcb, qcb, E_cb, x, divg, 
        tot_sink, cth_bin, rfv_acc, rfv_nrm, 
        rhoL, E_cb, qcb, ns, nz )
end

"initialize experiments with input parameters and empty output structures"
function define_experiments(; ctx::ModelContext)
    ( qm, qs, zcb, qcb, E_cb, x, divg, 
        tot_sink, cth_bin, rfv_acc, rfv_nrm, 
        rhoL, E_cb, qcb, ns, nz ) = setup_experiments(ctx=ctx)

    allocate_output() = ModelOutput( 
        Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
        Vector{Union{Missing, Float64}}(missing,    ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
        Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
        Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns) )
        
    # initialize experiment input and output structures
    control = Experiment(
        "control",
        "Control",
        ModelInput(qm, qs, zcb, qcb, E_cb, x, divg, tot_sink, cth_bin, rfv_acc, rfv_nrm),
        allocate_output()
    )

    subsminus5pct = Experiment(
        "subsidence-5%",
        "LS subsidence - 5%",
        ModelInput(qm, qs, zcb, qcb, E_cb, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
        allocate_output()
    )

    qsplus7pct = Experiment(
        "qs+7%",
        "qs + 7%, LS subsidence - 5%, LOW RH!",
        ModelInput(qm, qs*1.07, zcb, qcb*1.07, E_cb, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
        allocate_output()
    )

    qplus7pct = Experiment(
        "q&qs+7%",
        "q and qs + 7%, subsidence - 5%, RH=control",
        ModelInput(qm*1.07, qs*1.07, zcb, qcb*1.07, E_cb, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
        allocate_output()
    )

    ecbplus2pct = Experiment(
        "Ecb+2%",
        "E_cb + 2%, q,qs + 7%, subsidence - 5%",
        ModelInput(qm*1.07, qs*1.07, zcb, qcb*1.07, E_cb*1.02, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
        allocate_output()
    )

    qm_ = @. (1 - 0.95*(1-qm/qs)) * qs * 1.07
    cRHminus5pct = Experiment(
        "(1-RH)-5%",
        "subcloud (1-RH) - 5%, E_cb + 2%, q,qs + 7%, subsidence - 5%",
        ModelInput(qm_, qs*1.07, zcb, qcb*1.07, E_cb*1.02, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
        allocate_output()
    )

    # experiment dictionary for looping, and defining short names
    ExpDict = Dict(
        "control" => control,
        "subsidence-5%" => subsminus5pct,
        "qs+7%" => qsplus7pct,
        "q&qs+7%" => qplus7pct,
        "Ecb+2%" => ecbplus2pct,
        "(1-RH)-5%" => cRHminus5pct
    )
end

"large scale moisture flux profile distributed to clouds top height bins"
function calc_Ftot(; ctx::ModelContext,
                     E_cb, # W/m^2; just the cloud vapor flux; E0 - 35?
                     rhb_prate=8.88e-6, # kg/s
                     divg,
                     qm,
                     cth_bin,
                     cth_acc,
                     dz=ctx.dz,
                     icb=findfirst(ctx.z .>= ctx.zcb) )

    z = ctx.z
    zi = ctx.zi
    ztop = ctx.ztop
    rhoL = ctx.rhoL

    # align cth data; cth_bin starts at z=500
    offset = findfirst(x->x≈cth_bin[1]*1e3, z) - 1 # 50

    # cloud base water flux from vapor and precipitation flux
    # rhb_prate = mean(skipmissing(psl["prate"][:])) / 3600
    # all sky _total_ flux at cloud base
    calc_cloud_mean(s,z) = mean(skipmissing(tmean(s)[z.<=ztop])) # 286.6K = mean Cu layer temperature
    G_cb = E_cb/rhoL - rhb_prate

    # specify the LS drying source profile due to
    # subsidence and advection
    "large scale subsidence vertical velocity (m/s)"
    subsidence(z; divg=divg, zi=zi) = -min(z,zi) * divg # <0
    # ^ negate on rhs
    # large scale drying moisture source
    function largescale_drying(q,z)
        wdqdz = subsidence.(z[1:end-1]) .* ddz(q,z)
        S_ls = -wdqdz .- 1.7e-8*max.(0, (4e3.-z[1:end-1])/4e3) 
    end

    # define LS source
    S_ls = largescale_drying(qm,z)

    # Total all-sky flux G with cloud base BC
    G = Array{Union{Float64, Missing}}(missing, size(qm))
    G[icb:end] = G_cb .+ vcat(0, cumsum( S_ls[icb:end] ) * dz)

    F2 = calcF2(G, cth_acc, offset; sk=1) # divvies flux into cloud top height bins

    # coarsen bin size to 100 m
    # F2c = calcF2(G, cth_acc, offset; sk=10)
    # mean not weighted by area:
    # F2m = mean(skipmissing(nisf2m.(F2)))
    # mean weighted by area:
    # F2m = sum( skipmissing(nisf2m.(F2) .* rfv_nrm[1:end-1]) ) / sum( rfv_nrm )

    # align F2 to cloud model z grid
    F2z = Array{Union{Missing, Float64}}(missing, size(qm))
    F2z[offset.+eachindex(F2)] .= F2[:] 
    # needs to be masked with 0 above cloud top
    return F2z, G
end

"integrate an experiment based on inputs, modify output in place"
function integrate_experiment!(exp::Experiment; ctx::ModelContext)
    z = ctx.z
    dz = ctx.dz

    # compute cloud+precipitation moisture flux that balances large scale drying
    F2z, G_ls = calc_Ftot( ctx=ctx, E_cb=exp.input.E_cb, divg=exp.input.divg,
        qm=exp.input.qm, 
        cth_bin=exp.input.cth_bin, cth_acc=exp.input.cth_acc, 
        dz=dz, icb=findfirst(z .>= exp.input.zcb) )
    # Flux is distributed to the clouds by their cloud top height distribution.

    # run the cloud model for many sink rates
    zt, F_cld, F_pcp, qcld = cloudflux_1x(
        exp.input.tot_sink; x=exp.input.x, 
        z=z, nz=length(z), 
        dz=z[2]-z[1],
        qm=exp.input.qm, qs=exp.input.qs, 
        F2z=F2z, icb=findfirst(z .>= exp.input.zcb), 
        qcb=exp.input.qcb )

    # postprocess to get w, a, M, ...
    w, _ = updraft_w_dq(F_cld, qcld, exp.input.qm, z, zt)
    # print("size(w): $(size(w))") # nz, ns = (3100, 600)

    # Interpolate satellite coordinate to model sinkrate coordiante.
    # cloud fraction density per unit sink rate
    # da/dsinkrate = da/dh * dh/dsinkrate.
    da_dsink, da_ind = dadsinkrate(zt, exp.input.tot_sink, exp.input.cth_bin, exp.input.cth_nrm)
    acld = ctx.dsink * da_dsink # cloud area fraction in sink rate bin

    exp.output.w .= w
    exp.output.acld[da_ind] .= acld
    exp.output.M[:,da_ind] .= w[:,da_ind] .* acld' # mass flux per sink rate
    exp.output.qc .= qcld
    exp.output.F_cld .= F_cld
    exp.output.F_pcp .= F_pcp
    exp.output.G_cld[:,da_ind] .= F_cld[:,da_ind] .* acld'
    exp.output.G_pcp[:,da_ind] .= F_pcp[:,da_ind] .* acld'

    return nothing
end

end # module TradeCuExperiments

using .TradeCuExperiments

# run all the experiments and fill the output structures
ctx = init_context()
ExpDict = define_experiments(ctx=ctx)
for exp in values(ExpDict)
    println(exp.name)
    integrate_experiment!(exp, ctx=ctx)
end

println("Experiment names: $(keys(ExpDict))")
println("ctx: $(fieldnames(typeof(ctx)))")
control = ExpDict["control"]
println("control.input: $(fieldnames(typeof(control.input)))")
println("control.output: $(fieldnames(typeof(control.output)))")

keyorder = [
  "control", 
  "subsidence-5%",
  "qs+7%",
  "q&qs+7%",
  "Ecb+2%",
  "(1-RH)-5%",
  ]

PythonPlot.matplotlib.pyplot.close("all")
figure()
colrs = ["C0", "C1", "C2", "C3", "C4", "C5"]
for (i,k) in enumerate(keyorder)
    exp = ExpDict[k]
    plot( 1e3*exp.input.qm, 1e-3*ctx.z, label=k, color=colrs[i] )
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

# anther way to experiment is to keep the control distribution of
# sink rates, and recompute the clouds

# The clouds don't depend on the fluxes at all.
# So we can integrate the clouds separately,
# find the sink rates that give clouds of each height,
# assign fractions to these from the obseved 
# cloud top height and 
# partition the large-scale moisture flux to the cloud flux.
# Then finally calculate the cloud fluxes, mass fluxes, and velocities.

# get sink rate as a function of cloud top height
exp = ExpDict["control"]
iz = findall(ctx.zcb .<= ctx.z .<= ctx.ztop)
sinkz = NaN .+ zeros(length(ctx.z)) # fill a vector with NaN
TradeCuModel.find_contour!(sinkz[iz], exp.input.tot_sink,
    permutedims((exp.output.qc.-exp.input.qs)[iz,:]), 0.0)
plot(sinkz, ctx.z, "k"); display(gcf())

q = permutedims(exp.output.qc.-exp.input.qs)
# sinkz = similar(ctx.z, Union{eltype(ctx.z), Missing})
sinkz = NaN .+ zeros(length(ctx.z)) # fill a vector with NaN
iz = findall(ctx.zcb .<= ctx.z .<= ctx.ztop)
TradeCuModel.find_contour!(@view(sinkz[iz]), exp.input.tot_sink, q[:,iz], 0.0)
plot(sinkz, ctx.z, "k"); display(gcf())
# Now sinkz and ctx.z are 1 to 1, and 75:302 are valid low clouds.
# Map these to the cloud top height distribution 
# to get the cloud fraction for each of the interpolated sink rates in sinkz.
