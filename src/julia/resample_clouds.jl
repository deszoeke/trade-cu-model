using Revise
# cd("/Users/deszoeks/projects/ATOMIC/trade-cu-model/src/julia")
# using Pkg; Pkg.activate(".")
using NCDatasets
using Dates
using Printf
using Statistics
using LinearAlgebra
using PythonPlot
using PythonCall

includet("TradeCuExperiments.jl")
using .TradeCuExperiments

# utility functions; extend matplotlib

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

"plotf(f, x, args...)  function plot wrapper"
plotf(f, x, args...) = plot(f.(x), x, args...)

pd = permutedims

# set up plot defaults
PythonPlot.matplotlib.rcParams["font.family"] = "sans-serif"
PythonPlot.matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "OpenSans"]
# make fonts bigger by mutating rcParams
fsz = 12
font_settings = Dict(
    "font.size" => fsz,         # Base size
    "axes.titlesize" =>fsz+2,   # Subplot titles
    "axes.labelsize" => fsz+2,  # X/Y labels
    "xtick.labelsize" => fsz,   # X-axis numbers
    "ytick.labelsize" => fsz,   # Y-axis numbers
    "legend.fontsize" => fsz    # Legend text
)
matplotlib.rcParams.update(font_settings)

# get Simon's GOES albedo-weighted cloud fraction
function get_goes_cloud_ts()
    NCDataset("../../data/satellite/GOES-16/shcu_cloud_albedo_refl_profile_ts.nc") do dsa
        time = dsa[:time][:]
        rfv_nrm = dsa[:albedo_c_profile][:,:] # skipmissing: don't let fill values on individual days blank low-height bins
        rfv_acc = reverse(cumsum(reverse(rfv_nrm, dims=1), dims=1), dims=1)
        cth_bin = dsa[:cloud_top_height][:] # m
        return time, cth_bin,rfv_nrm, rfv_acc
    end
end
time, cth_bin, rfv_nrm, rfv_acc = get_goes_cloud_ts()

# first indices in each day (or after a gap of 6 hours)
iday = [0; findall(diff(time) .> Hour(6))] .+ 1

# temporal division
# mean first 10 days
ii = 1:iday[10+1]-1
ii = findall(time .<= time[1] + Day(10))
# last 10 days
jj = findall(time .>= time[end] - Day(10))
clf()
plot( mean(rfv_nrm[:,jj], dims=2), cth_bin)
plot( mean(rfv_nrm[:,ii], dims=2), cth_bin)

# EOFs of temporal cloud fraction variability
anom(x) = x .- mean(x, dims=2)

function cf_modes(rfv_nrm)
    A = anom(rfv_nrm)
    n = size(A, 2)
    # total variance profile directly from the matrix
    # var(A, dims=2) uses the standard (n - 1) normalization automatically
    spatial_variance = var(A, dims=2)
    total_matrix_variance = sum(spatial_variance)

    F = svd(A)
    U = F.U  # Spatial patterns (length = 401)
    S = F.S  # Singular values
    V = F.V  # Temporal patterns (length = 835)
    Vt = F.Vt # transpose for plotting
    U[:,1] .*= -1; V[1,:] .*= -1; Vt[:,1] .*= -1 # flip sign of first mode to match the mean profile

    mode_variances = S.^2 ./ (n - 1)
    U_scaled = U * Diagonal(S) # spatial modes
    mode_stds = S ./ sqrt(n - 1)
    U_physical = U * Diagonal(mode_stds) # spatial modes in physical units, fraction/10 m bin (std of cloud fraction)

    return U_physical, U_scaled, mode_stds, Vt, spatial_variance, total_matrix_variance
end
# log transform for normalizing relative cloud changes
cf_tx(c) = log(c + 1e-6)
cf_rx(x) = exp(x) - 1e-6
zz = 500.0 .<= cth_bin .<= 3100.0 # 261
z = cth_bin[zz]
c = rfv_nrm[zz,:]
x = cf_tx.(c)
cx0 = mean(x, dims=2)
# cf0 = cf_rx.(cx0)
cf0 = mean(rfv_nrm[zz,:], dims=2) # mean cloud fraction profile in the 500-3100 m range
U_physical, U_scaled, mode_stds, Vt, spatial_variance, total_matrix_variance = cf_modes(c)
n = size(rfv_nrm, 2)
# modes scaled to be anomalies of the mean cloud fraction profile
# scale mode amplitude by 1/sqrt(n-1)
U_mean = U_physical ./ sqrt( n - 1 )

# resample with mean +- 1 std of the first 2 modes: 9 experiments
"perturbation of the mean cloud fraction profile ∑ᵢ( coefsᵢ * U[:,i] )"
cf_prime(coefs, U=U_mean) = U[:, 1:length(coefs)] * coefs

clf()
plot([0,0], [0, 4], "k-", linewidth=0.5)
plot(cf0, z/1e3, "k--", linewidth=2)
for i in -1:1, j in -1:1
    plot((10*cf_prime([i, j]) .+ cf0)[:], z/1e3, label=@sprintf("%+d%+d", i, j))
    i>=0 && j>=0 && plot((10*cf_prime([i, j]))[:], z/1e3)
end
legend(frameon=false)
# 3 anomalous resamplings have negative cloud fraction: 0,+1; +1,-1; +1,+1
# nevertheless, the linear modes of cf show more variety in cloud top height.
# the log(cf) anomalies vary mostly in amplitude, in contrast.
# excluding these is fine
# 12
# --  0-  +-X
# -0  00  +0
# -+  0+X ++X
spatial_modes_variance = (U_scaled .^ 2) ./ (n - 1)
total_variance_profile = sum(spatial_modes_variance, dims=2)
cumulative_mode_variance = cumsum(spatial_modes_variance, dims=2)

# plot the vertical modes and their cumulative variance explained
clf()
subplot(1,2,1)
plot( [0,0], [0, 4], "k-", linewidth=0.5)
plot( cf0, z/1e3, linewidth=1, color="k", label="mean", linestyle="--")
plot( U_physical[:,1:3] .* [-1, -1, -1]', z/1e3, linewidth=2, label=["mode 1" "mode 2" "mode 3"])
plot( U_physical[:,4:6], z/1e3, linewidth=0.5, color="gray")
xlabel("cloud fraction modes")
ylabel("cloud top height (km)")
ylim([0, 3.5])
xlim([-2e-4, 6e-4])
legend(frameon=false)

subplot(1,2,2)
# plot( [0,0], [0, 4], "k-", linewidth=0.5)
plot( sqrt.(spatial_variance), z/1e3, "k--", linewidth=1, label="total std")
plot( sqrt.(cumulative_mode_variance[:,1:3]), z/1e3, linewidth=1.4)
xlabel("cloud fraction\ncumulative standard deviation")
ylim([0, 3.5])
xlim([0, 3.5e-4])
tight_layout()
# for f in ["svg", "png", "pdf", "eps"]
#     savefig("cf_modes." * f)
# end

# amplitude converges to the total std at each height
# first 3 modes explain 0.8-0.9 of the total variance

# mode 1 "more" is trade inversion cloud fraction (slightly obscures) low clouds.
# mode 2 "higher" trade inversion clouds (2-3 km), adjusts the height of the trade cu cloud top (inversion height); moves cloud higher.
# mode 3 "low" has more 0.5-1 km clouds, and vertically limited trade inversion clouds



# ( qm, qs, zcb, qcb, E_cb, x, divg, sfc_adv,
#     tot_sink, cth_bin, rfv_acc, rfv_nrm, 
#     rho, rhoL, ns, nz ) = setup_experiments(ctx=ctx)
ModelContext = TradeCuModel.ModelContext
ctx = init_context()

"initialize cloud fraction resampling experiments"
function define_cf_experiments(; ctx::ModelContext, cf0=cf0, U=U_mean, stds=10.0)

    ( qm, qs, zcb, qcb, E_cb, x, divg, sfc_adv,
      tot_sink, cth_bin, rfv_acc, rfv_nrm, 
      rho, rhoL, ns, nz ) = setup_experiments(ctx=ctx)
      
    icb = findfirst(ctx.z .>= zcb) # cloud base index
    
    # initialize control experiment input and output structures
    control = define_experiment(; name="control", description="Control",
        qm=qm, qs=qs, zcb=zcb, qcb=qcb, E_cb=E_cb, x=x, divg=divg, sfc_adv=sfc_adv, tot_sink=tot_sink, 
        cth_bin=cth_bin, rfv_acc=rfv_acc, rfv_nrm=rfv_nrm, 
        control=true, a_i_control=nothing, M_i_control=nothing )

    # integrate control
    integrate_experiment!(control, ctx=ctx)
    # control a_i and M_i will scale a for experiments
    a_i_control = control.output.acld # cloud area fraction for each cloud top height bin
    M_i_control = control.output.M # mass flux for each cloud top height bin

    # cf experiments: control experiments resample cloud fraction profiles,
    # ±1 std of the first 2 modes

    CfExpDict=Dict{String, TradeCuModel.Experiment}()
    for i in -1:1, j in -1:1
        # new cloud top height distribution resamples cloud fractions responsible
        # for the moisture flux.
        expname = @sprintf("cf %2d %2d", i, j)
        cf = zeros(size(cth_bin)) # initialize cloud fraction profile
        zz = 500.0 .<= cth_bin .<= 3100.0 # 261
        cf[zz] = TradeCuModel.filt_rfv( (stds*cf_prime([i, j]) .+ cf0)[:] )
        cumcf = reverse( cumsum(reverse( cf )) )

        control = define_experiment( control; name=expname, 
            description="cloud fraction resampled, cf0 + $(i)*mode1 + $(j)*mode2",
            qm=qm, qs=qs, zcb=zcb, qcb=qcb, E_cb=E_cb, x=x, divg=divg, 
            sfc_adv=sfc_adv, tot_sink=tot_sink, cth_bin=cth_bin, 
            rfv_acc=cumcf, rfv_nrm=cf,
            control=true, a_i_control=nothing, M_i_control=nothing )
        integrate_experiment!(control, ctx=ctx)
        # control a_i and M_i will scale a for experiments
        a_i_control = control.output.acld # cloud area fraction for each cloud top height bin
        M_i_control = control.output.M # mass flux for each cloud top height bin
        
        dimexp = define_experiment(control; name="DIM-" * expname, 
            description="DIM subsidence-5% Ecb+2% q,qs+7% " * @sprintf("cf%+1dM1%+1dM2",i,j),
            divg=divg*0.95, sfc_adv=0.95*sfc_adv,
            E_cb=E_cb*1.02,
            qs=qs*1.07, qcb=qcb*1.07, qm=qm*1.07,
            control=false, a_i_control=a_i_control, M_i_control=M_i_control )
        integrate_experiment!(dimexp, ctx=ctx)

        push!(CfExpDict, expname => control)
        push!(CfExpDict, "dim-" * expname => dimexp)
    end

    return CfExpDict
end

# initialize cf experiments
# run control simulations with resampled cf.
# use the control cloud fraction and mass flux for DIM experiments.
CfExpDict = define_cf_experiments(; ctx=ctx, cf0=cf0, U=U_physical)
# already done by define_cf_experiments:
# for expmt in values(CfExpDict) # run all the defined experiments
#     integrate_experiment!(expmt, ctx=ctx)
# end

# plot
# 12
# --  0-  +-X
# -0  00  +0
# -+  0+X ++X
clf()
fig = gcf()
fig.subplots_adjust(left=0.10, right=0.82, bottom=0.12, top=0.90, wspace=0.18, hspace=0.26)
axs = fig.subplots(3, 3)
log_norm = matplotlib.colors.LogNorm(vmin=1e-1, vmax=1e1)
levels = 10 .^ range(log10(1e-1), log10(1e1), length=13)
right_mappable = nothing
for i in -1:1, j in -1:1
    ax = axs[i+1, j+1]
    req = @sprintf("cf %2d %2d", i, j)
    data = CfExpDict[req].output.w
    cm = ax.contourf(CfExpDict[req].input.tot_sink*1e3, ctx.z/1e3,
        data, levels=levels, cmap=get_cmap("RdYlBu_r", 13), norm=log_norm, extend="both")
    if j == 1
        right_mappable = cm
    end
    ax.plot(2e3*CfExpDict[req].input.cth_nrm, CfExpDict[req].input.cth_bin/1e3, color="k", label="cloud fraction")
    ax.set_ylim([0, 3.5])
    ax.set_xlim([0, 5.6])
    j == -1 && ax.set_ylabel("height (km)", fontsize=12)
    i ==  1 && ax.set_xlabel("sink rate (1/km)", fontsize=12)
    ax.text(2, 2.5, @sprintf("%+d,%+d σ", 10*i, 10*j), fontsize=12)
    ax.tick_params(labelsize=12)
end
if right_mappable !== nothing
    cax = fig.add_axes([0.84, 0.23, 0.018, 0.42])
    fig.colorbar(right_mappable, cax=cax, extend="both")
end
fig.suptitle("cloud fraction mean ±10σ of modes 1 and 2", fontsize=14)
# for f in ["svg", "png", "pdf", "eps"]
#     savefig("resample_cf_w." * f)
# end

# plot cloud fraction ratio of DIM to control experiments
good(x) = !ismissing(x) && isfinite(x)
req_exp(i,j) = @sprintf("cf %2d %2d", i, j)
dim_exp_pair(D, req) = (D[req], D["dim-" * req])
function cf_plot_acld_ratio(c, e)
    plot(c.input.tot_sink*1e3, e.output.acld ./ c.output.acld .-1)
end

clf()
fig = gcf()
for i in -1:1, j in -1:1
    c,e = dim_exp_pair(CfExpDict, req_exp(i,j))
    cf_plot_acld_ratio(c, e)
end
# if a sink rate generates a cloud top height in the 
# observed cloud top height distribution, then its cloud fraction
# is always reduced by the same ratio in the DIM experiment, 
# regardless of the cloud fraction profile.
#
# Different cloud fraction are reduced by different fractions

c,e = dim_exp_pair(CfExpDict, req_exp(0,0))
icb = searchsortedlast(ctx.z, 700)
ii = @. good(e.output.M[icb,:]) & good(c.output.M[icb,:])
sum(e.output.M[icb,ii]) / sum(c.output.M[icb,ii]) # !almost! the same as acld
sum(e.output.acld[ii]) / sum(c.output.acld[ii])

for i in -1:1, j in -1:1
    c,e = dim_exp_pair(CfExpDict, req_exp(i,j))
    l = plot(c.output.acld, c.output.ztop/1e3)
    plot(e.output.acld, e.output.ztop/1e3, "--", color=l[0].get_color())
end
xlabel("cloud fraction per sink rate bin")
ylabel("cloud top height for sink rate (km)")
title("cloud fraction change for each moisture sink rate")

# look at diurnal cycle?
date = floor.(time, Day(1))
hourofday = Dates.value.(time.-date) / 86_400_000 * 24