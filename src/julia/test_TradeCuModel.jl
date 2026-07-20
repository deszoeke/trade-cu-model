using Revise
using Debugger
using Printf
using Statistics

using PythonPlot, PythonCall
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

include("TradeCuExperiments.jl") # nuts and bolts to get into experiments
using .TradeCuExperiments
calc_G_allsky = TradeCuModel.calc_G_allsky
import .TradeCuModel.interpolate_ascending

# set data
ctx = init_context()
( qm, qs, zcb, qcb, E_cb, x, divg, sfc_adv,
        tot_sink, cth_bin, rfv_acc, rfv_nrm, rho, rhoL, ns, nz ) = setup_experiments(ctx=ctx)
z = ctx.z
ztop = collect(710.0:50:6000.0) # cloud-top height grid
G_i, G_tot = calc_G_allsky(ztop; z=z, E_cb=E_cb, rhb_prate=0, divg=divg, sfc_adv=1.5e-8,
    qm=qm, rho=rho, rhoL=rhoL, zcb=zcb)
clf()
subplot(1,2,1)
plot(G_i/10, ztop/1e3)
xlabel("(kg/kg /s)\nall-sky eddy source ")
ylabel("height (km)")
subplot(1,2,2)
plot(G_tot, ztop/1e3)
xlabel("(kg/kg m/s)\nall-sky eddy flux G_tot")
xlim([0, 7e-5])
tight_layout()

rhb_prate=8.88e-6, # kg/m^2/s
sfc_adv=1.7e-8,
zi=4000.0, zcb=700.0

# all-sky total flux at cloud base
G_cb = E_cb/rhoL + rhb_prate/rho # 5e-5 kg/kg m/s about 5x the magnitude of the precipitation
# G_cb .*= 1.2

# large-scale subsidence vertical velocity (m/s)
subsidence(z_; divg=divg, zi=zi) = -min(z_, zi) * divg
"interpolate the derivative dq/dz at z_"
ddz(z_, q=qm,z=z) = interpolate_ascending( z[1:end-1].+0.5*diff(z), diff(q) ./ diff(z) )(z_)
"""
large-scale drying profile S_ls(z) = -w*dq/dz - 1.7e-8*(zdivg-z)/zdivg
kg/kg s^-1
"""
function largescale_drying(q, z_; divg=divg, sfc_adv=sfc_adv)
    zdivg = 4e3 # m; top height hardwired
    wdqdz = subsidence(z_; divg=divg, zi=zdivg) .* ddz(z_, q)
    -wdqdz - sfc_adv * max(0, (zdivg - z_) ./ zdivg) # advection increases to sfc_adv at surface
end
"integral of large scale drying between z0 and z1 (kg/kg m/s)"
S_ls_dz(z0,z1) = ( largescale_drying(qm, z0, divg=divg, sfc_adv=sfc_adv) + 
                   largescale_drying(qm, z1, divg=divg, sfc_adv=sfc_adv) ) * 0.5 * (z1-z0)

    nztop = length(ztop)
    nztop == 0 && return Float64[], Float64[]

    zt = coalesce.(ztop, NaN)
    finite = isfinite.(zt)
    valid = finite .& (zt .>= zcb)
    any(finite .& (zt .< zcb)) && @warn "Some zt < zcb."

    if !any(valid)
        return fill(NaN, nztop), fill(NaN, nztop)
    end
    # Total all-sky flux G with cloud-base boundary condition
    nztop = length(ztop)
    G_i   = Array{Union{Float64, Missing}}(missing, nztop)
    G_tot = Array{Union{Float64, Missing}}(missing, nztop)

    ii = findall(valid)
    iis = ii[sortperm(zt[ii])] # indices of ztop

    # G_i is eddy flux converging between z_i - z_i-1
    # G_i(zcb) = 0.0 # at cloud base
    dG = similarmissing((length(iis),), Float64)
    dG[1] = S_ls_dz(zcb, zt[iis[1]]) # first cloud
    dG[2:end] = S_ls_dz.(zt[iis[1:end-1]], zt[iis[2:end]]) # subsequent clouds
    # G_tot is total all-sky eddy flux
    # G_tot(zcb) = G_cb
    G_i[iis] = -dG
    G_tot[iis] = G_cb .+ cumsum(dG) # vertically integrated all-sky flux
