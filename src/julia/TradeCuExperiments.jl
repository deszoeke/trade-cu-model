# environment and code loading header
# using Revise
# joinpath(homedir(), "projects/ATOMIC/trade-cu-model/src/julia") |> p-> ispath(p) && cd(p)
# # joinpath(homedir(), "Projects/ATOMIC/trade-cu-model/src/julia") |> p-> ispath(p) && cd(p)
# using Pkg; Pkg.activate(".") # may be superfluous if opening julia --project in local folder

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

export ExpDict # dictionary contains defined experiments
export init_context, define_experiments, define_experiment
export interp_sinkrate #, get_sinkrate
export integrate_experiment!
# temporarily export more to experiment in outside environment
export setup_experiments
export cloud_i_area
export test_control_sink
export calc_ql
export new_area


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
    # reverse so that tot_sink decreases and ztop increases
    tot_sink = reverse(range(6.3523e-4, 5.7e-3, length=600)) # min tuned for x=0.53 to get the highest possible cloud top
    dsink = ctx.dsink
    # tot_sink = (1 .+tanh.(range(-8*pi, 0, length=600))) .* (5e-3 - 1e-4) .+ 1e-4
    # tot_sink = range(6.1716e-4, 5.8e-3, length=600) # min tuned for x=0.53 to get the highest possible cloud top

    z = ctx.z
    qm = ctx.qm
    qs = ctx.qs
    dz = ctx.dz
    # icb = findfirst(z .>= zcb) # cloud base index

    cth_bin = ctx.cth_bin # preliminary context gets clobbered in some experiments
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

allocate_output(nz,ns) = ModelOutput( 
    Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
    Vector{Union{Missing, Float64}}(missing,    ns), Vector{Union{Missing, Float64}}(missing,    ns), 
    Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
    Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
    Matrix{Union{Missing, Float64}}(missing, nz,ns) )

function define_experiment(; name, description, 
    qm, qs, zcb, qcb, E_cb, x, divg, tot_sink, cth_bin, rfv_acc, rfv_nrm, control=true, a_i_control=nothing, M_i_control=nothing)
    Experiment( name, description,
        ModelInput(qm, qs, zcb, qcb, E_cb, x, divg, tot_sink, cth_bin, rfv_acc, rfv_nrm, control, a_i_control, M_i_control),
        allocate_output(length(qm),length(tot_sink)) )
end

"initialize experiments with input parameters and empty output structures"
function define_experiments(; ctx::ModelContext)
    ( qm, qs, zcb, qcb, E_cb, x, divg, 
        tot_sink, cth_bin, rfv_acc, rfv_nrm, 
        rhoL, E_cb, qcb, ns, nz ) = setup_experiments(ctx=ctx)
        
    # initialize experiment input and output structures
    control = define_experiment(; name="control", description="Control",
        qm=qm, qs=qs, zcb=zcb, qcb=qcb, E_cb=E_cb, x=x, divg=divg, tot_sink=tot_sink, 
        cth_bin=cth_bin, rfv_acc=rfv_acc, rfv_nrm=rfv_nrm, control=true, a_i_control=nothing, M_i_control=nothing )
    # integrate control to pass control a_i, M_i to experiments
    integrate_experiment!(control, ctx=ctx)
    # control a_i, M_i scale updraft w for experiments
    a_i_control = control.output.acld # cloud area fraction for each cloud top height bin
    M_i_control = control.output.M # mass flux for each cloud top height bin

    subsminus5pct = define_experiment(
        name="subsidence-5%", description="LS subsidence-5%",
        qm=qm, qs=qs, zcb=zcb, qcb=qcb, E_cb=E_cb, x=x, divg=divg*0.95, tot_sink=tot_sink, 
        cth_bin=cth_bin, rfv_acc=rfv_acc, rfv_nrm=rfv_nrm, control=false, a_i_control=a_i_control, M_i_control=M_i_control )

    qsplus7pct = define_experiment(
        name="qs+7%", description="qs+7%, LS subsidence-5%, LOW RH!",
        qm=qm, qs=qs*1.07, zcb=zcb, qcb=qcb*1.07, E_cb=E_cb, x=x, divg=divg*0.95, tot_sink=tot_sink, 
        cth_bin=cth_bin, rfv_acc=rfv_acc, rfv_nrm=rfv_nrm, control=false, a_i_control=a_i_control, M_i_control=M_i_control )

    qplus7pct = define_experiment(
        name="q&qs+7%", description="q and qs +7%, subsidence-5%, RH=control",
        qm=qm*1.07, qs=qs*1.07, zcb=zcb, qcb=qcb*1.07, E_cb=E_cb, x=x, divg=divg*0.95, tot_sink=tot_sink, 
        cth_bin=cth_bin, rfv_acc=rfv_acc, rfv_nrm=rfv_nrm, control=false, a_i_control=a_i_control, M_i_control=M_i_control )

    ecbplus2pct = define_experiment(
        name="Ecb+2%", description="E_cb + 2%, q&qs+7%, subsidence-5%",
        qm=qm*1.07, qs=qs*1.07, zcb=zcb, qcb=qcb*1.07, E_cb=E_cb*1.02, x=x, divg=divg*0.95, tot_sink=tot_sink, 
        cth_bin=cth_bin, rfv_acc=rfv_acc, rfv_nrm=rfv_nrm, control=false, a_i_control=a_i_control, M_i_control=M_i_control )

    qm_new = @. (1 - 0.95*(1-qm/qs)) * qs * 1.07
    cRHminus5pct = define_experiment(
        name="(1-RH)-5%", description="subcloud (1-RH)-5%, E_cb+2%, q&qs+7%, subsidence-5%",
        qm=qm_new, qs=qs*1.07, zcb=zcb, qcb=qcb*1.07, E_cb=E_cb*1.02, x=x, divg=divg*0.95, tot_sink=tot_sink, 
        cth_bin=cth_bin, rfv_acc=rfv_acc, rfv_nrm=rfv_nrm, control=false, a_i_control=a_i_control, M_i_control=M_i_control )

    # "DIM" is exactly as "(1-RH)-5%" above
    DIM = define_experiment(
        name="DIM", description="Descent Inhibited Moisture Flux; a_i-5%, (1-RH)-5%, E_cb+2%, q&qs+7%, subsidence-5%",
        qm=qm_new, qs=qs*1.07, zcb=zcb, qcb=qcb*1.07, E_cb=E_cb*1.02, x=x, divg=divg*0.95, tot_sink=tot_sink, 
        cth_bin=cth_bin, rfv_acc=rfv_acc, rfv_nrm=rfv_nrm, control=false, a_i_control=a_i_control, M_i_control=M_i_control )

    # experiment dictionary for looping, and defining short names
    ExpDict = Dict(
        "control" => control,
        "subsidence-5%" => subsminus5pct,
        "qs+7%" => qsplus7pct,
        "q&qs+7%" => qplus7pct,
        "Ecb+2%" => ecbplus2pct,
        "(1-RH)-5%" => cRHminus5pct,
        "DIM" => DIM
    )
end

# should just use define_experiment to define a sink experiment
function define_control_sink_experiment(; ctx::ModelContext, sinkz=sinkz, 
    name="control-sink", description="control sink rate experiment",
    control=true, a_i_control=nothing, M_i_control=nothing)
    ( qm, qs, zcb, qcb, E_cb, x, divg, 
        tot_sink, cth_bin, rfv_acc, rfv_nrm, 
        rhoL, E_cb, qcb, ns, nz ) = setup_experiments(ctx=ctx)

    return Experiment( name, description,
            ModelInput(qm, qs, zcb, qcb, E_cb, x, divg, sinkz, cth_bin, rfv_acc, rfv_nrm,
            control, a_i_control, M_i_control),
            allocate_output(nz, length(sinkz)) )
end

# function define_DIM_sink_experiment(; ctx::ModelContext, sinkz=sinkz)
#     ( qm, qs, zcb, qcb, E_cb, x, divg, 
#         tot_sink, cth_bin, rfv_acc, rfv_nrm, 
#         rhoL, E_cb, qcb, ns, nz ) = setup_experiments(ctx=ctx)

#     qm_new = @. (1 - 0.95*(1-qm/qs)) * qs * 1.07
#     return Experiment(
#             "DIM-z", "DIM with z≈ztop",
#             ModelInput(qm_new, qs*1.07, zcb, qcb*1.07, E_cb*1.02, x, divg*0.95, sinkz, cth_bin, 0.95*rfv_acc, 0.95*rfv_nrm),
#             allocate_output(nz,nz) ) # note nz,nz
# end

"gets GOES area for nearest cloud top height"
function cloud_i_area( ctx::ModelContext )
    z=ctx.z
    nz = length(ctx.z)
    cth_bin=ctx.cth_bin
    rfv_acc=ctx.rfv_acc
    # align z grid with fractions
    offset = findfirst(x-> abs(x-cth_bin[1]) < 1.0, z) - 1 # offset for 10 m bins
    a_i = zeros(nz)
    # cloud with qc=0 at h_i diverges below between h_(i-1) and h_i
    a_i[offset+1 .+ eachindex(rfv_acc[1:end-1])] .= -diff(rfv_acc)
    return a_i
end

"gets GOES area for nearest cloud top height"
function cloud_i_area( exp::Experiment; ctx::ModelContext )
    z=ctx.z
    nz = length(ctx.z)
    cth_bin=exp.input.cth_bin # experiment cloud top height supersedes context
    cth_acc=exp.input.cth_acc
    # align z grid with fractions
    offset = findfirst(x-> abs(x-cth_bin[1]) < 1.0, z) - 1 # offset for 10 m bins
    a_i = zeros(nz)
    # cloud with qc=0 at h_i diverges below between h_(i-1) and h_i
    a_i[offset+1 .+ eachindex(cth_acc[1:end-1])] .= -diff(cth_acc)
    return a_i
end

# Use calc_Ftot(...) directly for flux diagnostics.

# The in-cloud moisture flux F_i (due to cloud updraft + precipitation) for cloud i with area fraction a_i is constant from cloud base to cloud top height h_i. The all-sky eddy flux is G = sum{F_i a_i}. The total cloud area is a = sum{a_i}. The updraft flux for cloud i is Fcld_i = M_i (\Delta q)_i. The mass flux is M_i = w_i a_i. The total moisture flux F_i = Fcld_i + Fp_i is the sum of the updraft (positive, upward) and precipitation (negative, downward) fluxes. We have the all sky flux from large-scale balance; partition it to the total eddy flux F_i according to a_i of each cloud top height h_i, then decompose to M_i and w_i. 
# refine the sink rate and cloud top height interpolation

"binear interpolation of y(x) between (x1,y1) and (x2,y2)"
bilinear(x1,x2, y1,y2, x) = ( x2 == x1 ? y1 : y1 + (y2-y1) * (x-x1) / (x2-x1) )

# functions for computing new cloud area a_i for a new experiment 
# assuming the w_i are constant.
"normz(ztop) normalized z coordinate between zcb and ztop"
normz(z, ztop, zcb) = (z - zcb) / (ztop - zcb)

function stretch_to_new_ztop( oldQ::Vector, oldztop::Number, newztop::Number, ctx::ModelContext )
    z = ctx.z
    zcb = ctx.zcb
    newQ = fill(NaN, size(z))
    good(x) = !ismissing(x) && isfinite(x)
    if !ismissing(oldztop) && isfinite(oldztop) && !ismissing(newztop) && isfinite(newztop)
        Q = coalesce.(oldQ, NaN)
        ii = good.(oldQ)
        jj = zcb .<= z .<= newztop # interpolate at z between cloud base and NEW cloud top height
        newQ[jj] = interpolate_ascending( normz.(z, oldztop, zcb)[ii], Q[ii] ).(normz.(z[jj], newztop, zcb))
    end
    return newQ
end

"""
new_area(M::Vector, ztop::Union{Missing,Float64}, w_old, ztop_old, ctx)
    new experiment cloud area, assuming stretching invariant w_old
"""
function new_area(M::Vector, ztop::Union{Missing,Float64}, w_old, ztop_old, ctx)
    # new cloud area fraction for new ztop, given mass flux and old velocity
    M ./ stretch_to_new_ztop( w_old, ztop_old, ztop, ctx )
end
"""
new_area(expt::Experiment, ctl::Experiment)
    compute the cloud area for a new experiment expt mass flux and cloud top height
    given the old experiment ctl w and ztop.
"""
function new_area(expt::Experiment, ctl::Experiment, ctx::ModelContext)
    z = ctx.z
    zcb = ctx.zcb
    ztop    = expt.output.ztop
    ztop_old = ctl.output.ztop
    a_i_new = fill(NaN, size(expt.output.M))
    for i in eachindex(ztop)
        if !ismissing(ztop[i]) && isfinite(ztop[i]) && !ismissing(ztop_old[i]) && isfinite(ztop_old[i])
            @assert ztop[i] > zcb
            tmp = new_area( expt.output.M[:,i], ztop[i], ctl.output.w[:,i], ztop_old[i], ctx )
            a_i_new[:,i] = coalesce.(tmp, NaN)
        end
    end
    return a_i_new
end

# other potentially helpful functions for computing cloud properties

calc_ql(exp) = @. max(0, exp.output.qc - exp.input.qs)

"average cloud properties s[zt] in cloud top height bins zbin[inbin-1] < zt <= zbin[inbin]"
function meanbin_scalar!(a::Vector, c::Vector, s::Vector, zt::Vector, zbin::Vector)
    eachindex(zt)[ findfirst( zt.>=700.0 ):end ]
    for i in eachindex(zt)
        # find the bin index for cloud top height zbin[inbin-1] < zt <= zbin[inbin]
        inbin = searchsortedfirst( zbin, zt[i] )
        if !isempty(inbin) && inbin <= length(a)
            a[inbin] += s[i]
            c[inbin] += 1
        end
    end
    return a ./ c # return the mean
end
function meanbin_scalar(s,zt, zbin)
    a = zeros(length(zbin))
    c = zeros(length(zbin))
    meanbin_scalar!(a,c, s,zt, zbin)
    return a
end

"average cloud properties s[:,zt] in cloud top height bins zbin[inbin-1] < zt <= zbin[inbin]"
function meanbin_vec!(a::Matrix, c::Matrix, s::Matrix, zt::Vector, zbin::Vector)
    eachindex(zt)[ findfirst( zt.>=700.0 ):end ]
    for i in eachindex(zt)
        # find the bin index for cloud top height zbin[inbin-1] < zt <= zbin[inbin]
        inbin = searchsortedfirst( zbin, zt[i] )
        if !isempty(inbin) && inbin <= length(a)
            a[:,inbin] .+= s[:,i]
            c[:,inbin] .+= 1
        end
    end
    return a ./ c # return the mean
end
function meanbin_vec(s,zt, zbin)
    a = zeros(size(s,1), length(zbin))
    c = zeros(size(s,1), length(zbin))
    meanbin_vec!(a,c, s,zt, zbin)
    return a
end

# EXPERIMENT TESTS
"""
inject this exeriment function within TradeCuExperiments module
to test the sink rate setting experiment.
"""
function test_control_sink()
    # run all the experiments and fill the output structures
    ctx = init_context()
    ExpDict = define_experiments(ctx=ctx)
    for exp in values(ExpDict) # run all the defined experiments
        println(exp.name)
        integrate_experiment!(exp, ctx=ctx)
    end

    # experiment hierarchy
    keyorder = [
    "control", 
    "subsidence-5%",
    "qs+7%",
    "q&qs+7%",
    "Ecb+2%",
    "(1-RH)-5%",
    "DIM" ]

    # get sink rate as a function of cloud top height from control
    sinkz_full = interp_sinkrate( ExpDict["control"]; ctx=ctx ) # nz-length, extrapolates!
    ii_cld = findall(700 .≤ ctx.z .≤ 4000)  # 331 elements matching GOES CTH range
    sinkz = sinkz_full[ii_cld]              # sink rates for cloud tops at 700–4000 m

    # a_i = cloud_i_area( ctx )
    controlsink = define_control_sink_experiment(ctx=ctx, sinkz=sinkz,
        name="control-sink", description="control sink rate experiment",
        control=true, a_i_control=nothing, M_i_control=nothing)
    a_i_control = controlsink.output.acld # sized for 331 sink rate bins
    M_i_control = controlsink.output.M
    
    sinkm5 = define_control_sink_experiment(ctx=ctx, sinkz=0.95*sinkz,
        name="control-sink-5%", description="control sink rate experiment -5%",
        control=false, a_i_control=a_i_control, M_i_control=M_i_control)
    sinkp5 = define_control_sink_experiment(ctx=ctx, sinkz=1.05*sinkz,
        name="control-sink+5%", description="control sink rate experiment +5%",
        control=false, a_i_control=a_i_control, M_i_control=M_i_control)
    integrate_experiment!(controlsink, ctx=ctx)
    integrate_experiment!(sinkm5,      ctx=ctx)
    integrate_experiment!(sinkp5,      ctx=ctx)
    # append the controlsink experiments to the dictionary
    for exp in [controlsink, sinkm5, sinkp5]
        push!(ExpDict, exp.name => exp)
    end


    # The clouds don't depend on the fluxes at all.
    # Another way to experiment is to keep the control distribution of
    # sink rates, and recompute the clouds

    # So we can integrate the clouds separately,
    # find the sink rates that give clouds of each height,
    # assign fractions to these from the obseved 
    # cloud top height and 
    # partition the large-scale moisture flux to the cloud flux.
    # Then finally calculate the cloud fluxes, mass fluxes, and velocities.
    return ctx, ExpDict, controlsink, sinkm5, sinkp5
end

end # module TradeCuExperiments
