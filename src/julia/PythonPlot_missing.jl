using PythonPlot
using PythonCall

"Recast julia Arrays with missing as python masked arrays."
function PythonCall.Py(a::Array{Union{T,Missing},N}) where {T,N}
    np = PythonCall.pyimport("numpy")
    np.ma.masked_invalid(np.array(replace(a, missing => NaN)))
end

# test
plot([1,2,missing,4,5], marker="o")
