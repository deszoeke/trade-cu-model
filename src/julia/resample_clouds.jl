using Revise
# cd("/Users/deszoeks/projects/ATOMIC/trade-cu-model/src/julia")
# using Pkg; Pkg.activate(".")
using NCDatasets
using Dates
using Statistics
using LinearAlgebra
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
A = anom(rfv_nrm)
n = size(A, 2)
# total variance profile directly from the matrix
# var(A, dims=2) uses the standard (n - 1) normalization automatically
spatial_variance = var(A, dims=2)
total_matrix_variance = sum(spatial_variance)

F = svd(A)
U = F.U  # Spatial patterns (length = 401)
S = F.S  # Singular values
V = F.V  # Temporal patterns (length = 730)
Vt = F.Vt # transpose for plotting
U[:,1] .*= -1; V[1,:] .*= -1; Vt[:,1] .*= -1 # flip sign of first mode to match the mean profile

mode_variances = S.^2 ./ (n - 1)
U_scaled = U * Diagonal(S) # spatial modes
mode_stds = S ./ sqrt(n - 1)
U_physical = U * Diagonal(mode_stds)
spatial_modes_variance = (U_scaled .^ 2) ./ (n - 1)
total_variance_profile = sum(spatial_modes_variance, dims=2)
cumulative_mode_variance = cumsum(spatial_modes_variance, dims=2)

clf()
subplot(1,2,1)
plot([0,0], [0, 4], "k-", linewidth=0.5)
plot(mean(rfv_nrm, dims=2), cth_bin/1e3, linewidth=1, color="k", label="10mean", linestyle="--")
plot( U_physical[:,1:3], cth_bin/1e3, linewidth=2, label=["1" "2" "3"])
plot( U_physical[:,4:6], cth_bin/1e3, linewidth=0.7)
xlabel("cloud fraction modes")
ylabel("cloud top height (km)")
ylim([0, 4])
# legend(frameon=false)

subplot(1,2,2)
plot([0,0], [0, 4], "k-", linewidth=0.5)
plot(sqrt.(spatial_variance), cth_bin/1e3, "k--", linewidth=1, label="total std")
plot(sqrt.(cumulative_mode_variance[:,1:3]), cth_bin/1e3, linewidth=1.4)
xlabel("cloud fraction\ncumulative amplitude")
ylim([0, 4])
xlim([0, 3.5e-4])
tight_layout()
# amplitude converges to the total std at each height
# first 3 modes explain 0.8-0.9 of the total variance

# mode 1 is trade inversion cloud fraction (slightly obscures) low clouds.
# mode 2 adjusts the height of the trade cu cloud top (inversion height); moves cloud higher.
# mode 3 is slightly higher 0.5-1 km clouds, and slightly lower trade inversion clouds

# resample with mean +- 1 std of the first 2 modes