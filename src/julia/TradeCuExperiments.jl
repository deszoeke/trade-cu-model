# environment and code loading header
using Revise 
using Pkg; Pkg.activate(".")

includet("TradeCuModel.jl")

using PythonPlot
using Statistics

if @isdefined(PythonPlot)
    using PythonCall
    "Recast julia Arrays with missing as python masked arrays."
    function PythonCall.Py(a::Array{Union{T,Missing},N}) where {T,N}
        np = PythonCall.pyimport("numpy")
        np.ma.masked_invalid(np.array(replace(a, missing => NaN)))
    end

    "close plots with close(\"all\")"
    close(s::String) = PythonCall.pyimport("matplotlib.pyplot").close(s)
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
using PythonPlot
using Statistics
using VaporSat # dev ../../deps/VaporSat

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

# initialize model inputs

# parameters and initialization
KelvinCelsius=273.15 # K
zi   = 4.0e3    # m  # inversion for tapering subsidence
ztop = 4.0e3    # m  # top of cloud model integration domain
zcb  = 700.0    # m
divg = 1.5e-6   # 1/s
x = 0.53 # parameter precipitation efficiency

# ensemble of sink rates
tot_sink = range(6.3523e-4, 5.7e-3, length=600) # min tuned for x=0.53 to get the highest possible cloud top
dsink = tot_sink[2] - tot_sink[1]
# tot_sink = (1 .+tanh.(range(-8*pi, 0, length=600))) .* (5e-3 - 1e-4) .+ 1e-4
# tot_sink = range(6.1716e-4, 5.8e-3, length=600) # min tuned for x=0.53 to get the highest possible cloud top

z, tam, thm, qm, pm = get_mean_soundings()
#m  K       kg/kg Pa 
qs  = qsat.(pm, tam.-KelvinCelsius) # kg/kg
tvm  = virtual_temp.(tam, qm) # virtual temperature, K
# thvm = virtual_temp.(thm, qm) # virtual potential temperature, K
dz = z[2]-z[1]
# icb = findfirst(z .>= zcb) # cloud base index

rfv_nrm, rfv_acc, cth_bin = get_goes_cloud_data() # reflectivity data
#size(rfv_nrm), size(cth_bin) # (351,1)

rhoL = mean(filter(isfinite, calc_rhoL.(tvm, pm)[z.<=ztop])) # 2.41e6 J/m^3
E_cb = 180.0 # W/m^2; just the cloud vapor flux; E0 - 35?
qcb = qs[findfirst(z .>= zcb)] # kg/kg; cloud base specific humidity
ns = length(tot_sink)
nz = length(z)

# initialize experiment input and output structures
control = Experiment(
    "control",
    "Control",
    ModelInput(qm, qs, zcb, qcb, E_cb, x, divg, tot_sink, cth_bin, rfv_acc, rfv_nrm),
    ModelOutput( Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Vector{Union{Missing, Float64}}(missing,    ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns) )
)

subsminus5pct = Experiment(
    "subsidence-5%",
    "LS subsidence - 5%",
    ModelInput(qm, qs, zcb, qcb, E_cb, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
    ModelOutput( Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Vector{Union{Missing, Float64}}(missing,    ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns) )
)

qsplus7pct = Experiment(
    "qs+7%",
    "qs + 7%, LS subsidence - 5%, LOW RH!",
    ModelInput(qm, qs*1.07, zcb, qcb*1.07, E_cb, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
    ModelOutput( Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Vector{Union{Missing, Float64}}(missing,    ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns) )
)

qplus7pct = Experiment(
    "q&qs+7%",
    "q and qs + 7%, subsidence - 5%, RH=control",
    ModelInput(qm*1.07, qs*1.07, zcb, qcb*1.07, E_cb, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
    ModelOutput( Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Vector{Union{Missing, Float64}}(missing,    ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns) )
)

ecbplus2pct = Experiment(
    "Ecb+2%",
    "E_cb + 2%, q,qs + 7%, subsidence - 5%",
    ModelInput(qm*1.07, qs*1.07, zcb, qcb*1.07, E_cb*1.02, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
    ModelOutput( Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Vector{Union{Missing, Float64}}(missing,    ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns) )
)

qm_ = @. (1 - 0.95*(1-qm/qs)) * qs * 1.07
cRHminus5pct = Experiment(
    "(1-RH)-5%",
    "subcloud (1-RH) - 5%, E_cb + 2%, q,qs + 7%, subsidence - 5%",
    ModelInput(qm_, qs*1.07, zcb, qcb*1.07, E_cb*1.02, x, divg*0.95, tot_sink, cth_bin, rfv_acc, rfv_nrm),
    ModelOutput( Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Vector{Union{Missing, Float64}}(missing,    ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns),
                 Matrix{Union{Missing, Float64}}(missing, nz,ns), Matrix{Union{Missing, Float64}}(missing, nz,ns) )
)

# experiment dictionary for looping, and defining short names
ExpDict = Dict(
    "control" => control,
    "subsidence-5%" => subsminus5pct,
    "LSqs+7%" => LSqsplus7pct,
    "q&qs+7%" => qplus7pct,
    "Ecb+2%" => ecbplus2pct,
    "(1-RH)-5%" => cRHminus5pct
)

"large scale moisture flux profile distributed to clouds top height bins"
function calc_Ftot(; E_cb=E_cb, # W/m^2; just the cloud vapor flux; E0 - 35?
                     rhb_prate=8.88e-6, # kg/s
                     divg=divg,
                     qm=qm,
                     cth_bin=cth_bin,
                     cth_acc=rfv_acc,
                     dz=dz,
                     icb=findfirst(z .>= zcb) )

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
function integrate_experiment!(exp::Experiment)
    # compute cloud+precipitation moisture flux that balances large scale drying
    F2z, G_ls = calc_Ftot( E_cb=exp.input.E_cb, divg=exp.input.divg, 
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
    print("size(da_dsink): $(size(da_dsink))") # <600 
    acld = dsink * da_dsink # cloud area fraction in sink rate bin
    print("size(acld): $(size(acld))")

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

# run all the experiments and fill the output structures
for exp in values(ExpDict)
    integrate_experiment!(exp)
end