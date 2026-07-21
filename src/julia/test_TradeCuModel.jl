using Revise
using Debugger
using Printf
using Statistics
using NCDatasets

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

good(x) = !ismissing(x) && isfinite(x)
psl = NCDataset("../../data/flux/EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc");
rhb_prate = mean(skipmissing(psl[:prate][:])/3600*rhoL) # -> 22 W/m^2
rhb_evap  = mean(skipmissing(psl[:evap ][:])/3600*rhoL) # -> 182 W/m^2

# set data
ctx = init_context()
( qm, qs, zcb, qcb, E_cb, x, divg, sfc_adv,
        tot_sink, cth_bin, rfv_acc, rfv_nrm, rho, rhoL, ns, nz ) = setup_experiments(ctx=ctx)
z = ctx.z
ztop = collect(710.0:50:8000.0) # cloud-top height grid
G_i, G_tot = calc_G_allsky(ztop; z=z, E_cb=rhb_evap, rhb_prate=rhb_prate+5,
    divg=divg, sfc_adv=1.7e-8,
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
