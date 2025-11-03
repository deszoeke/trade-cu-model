using Revise

using NCDatasets
using PyPlot
using Interpolations

using PyCall
# allow for plotting with missing values
function PyCall.PyObject(a::Array{Union{T,Missing},N}) where {T,N}
    numpy_ma = PyCall.pyimport("numpy").ma
    pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
end


"""
get the indices of y that correspond for each variable x
so y[index[i]] <= x[i] < y[index[i+1]]
"""
coarseind(x, y) = [searchsortedlast(y, xi) for xi in x]
#= "with type stability and preallocation"
function coarseind(x, y)
   indices = Vector{Int}(undef, length(x))
   for i in eachindex(x, indices)
       indices[i] = searchsortedlast(y, x[i])
   end
   return indices
end =#

#url = "https://www.ncei.noaa.gov/thredds-ocean/catalog/psl/atomic/p3/Picarro/catalog.html" # thredds
#url = "https://www.ncei.noaa.gov/data/oceans/oar/psl/atomic-2020/p3/Picarro/"
# ds = NCDataset(url)
p3dir = "../../data/p3"
p3files = readdir(joinpath(p3dir, "Picarro")) |> filter(startswith("EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-Isotope-Ratios-1Hz"))
p3h2ofiles = readdir(joinpath(p3dir, "Picarro")) |> filter(startswith("EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-1Hz"))
p3metfiles = readdir(joinpath(p3dir, "meteorology")) |> filter(startswith("EUREC4A_ATOMIC_P3_Flight-Level"))

#dsi = NCDataset( joinpath(p3dir, p3files[1]) )

"filter out takeoff and landing"
function filt_flight( dsi )
    # exclude take-off and landing to avoid sampling island air
    takeoff = findfirst(dsi[:alt][:] .>= 2500 )
    landing = findlast( dsi[:alt][:] .>= 2500 )
    t_h = (dsi[:time][:].-dsi[:time][takeoff])./Hour(1) # hours since first ascent
    ii = findall(    eachindex(t_h) .>= takeoff 
                 .&& eachindex(t_h) .<= landing 
                 .&& dsi[:lon][:] .>= -58.5     )
    return ii, t_h
end

"read flight isotope file, filter takeoff and landing"
function get_iso_ds( f )
    dsi = NCDataset( f )
    ii, t_h = filt_flight( dsi )
    return dsi, ii, t_h
end

function plot_flight_alt( dsi, ii )
    subplot(2,1,1)
    plot(dsi[:lon][ii], dsi[:alt][ii]/1e3)
    subplot(2,1,2)
    plot(t_h[ii], dsi[:alt][ii]/1e3)
end

function plot_iso_prof( dsi, ii )
    plot(dsi[:dD][ii], dsi[:alt][ii]/1e3, linewidth=0.2)
    ylabel("height (km)")
    xlabel(L"\delta"*"D "*L"(10^{-3})")
end

# plot profiles for all flights
clf()
fig, ax = subplots(2,1,1)
for f in joinpath.(p3dir, "Picarro", p3files[1])
    dsi, ii, t_h = get_iso_ds( f )
    dstring = match(r"\d{4}\d{2}\d{2}", f).match
    h2ofile = filter(contains(dstring), p3h2ofiles)[1] # matching water vapor file
    dsw = NCDataset( h2ofile )
    plot_iso_prof( dsi, ii )
end
xlim([-500, -65])
ylim([0, 4])

