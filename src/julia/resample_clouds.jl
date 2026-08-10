using Revise
# cd("/Users/deszoeks/projects/ATOMIC/trade-cu-model/src/julia")
# using Pkg; Pkg.activate(".")
using NCDatasets
using Statistics
using PythonPlot
using PythonCall


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
function get_goes_cloud_data()
    NCDataset("../../data/satellite/GOES-16/shcu_cloud_albedo_refl_profile_mean.nc") do dsa
        rfv_nrm = dsa[:albedo_profile][:] # skipmissing: don't let fill values on individual days blank low-height bins
        rfv_acc = reverse(cumsum(reverse(rfv_nrm)))
        cth_bin = dsa[:cloud_top_height][:] # m
        return rfv_nrm, rfv_acc, cth_bin
    end
end