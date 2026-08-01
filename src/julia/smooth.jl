using Revise
using Statistics


"moving_average(A, m) m-point moving average of A"
function moving_average(A::AbstractArray, m::Union{Int,Vector{Int},Tuple{Int}}, cond=x->true)
    out = similar(A)
    R = CartesianIndices(A)
    Ifirst, Ilast = first(R), last(R)
    I1 = CartesianIndex(Tuple(m.÷2) .* Tuple(oneunit(Ifirst))) # index 1/2 window
    for I in R
        n, s = 0, zero(eltype(out))
        for J in max(Ifirst, I-I1):min(Ilast, I+I1)
            if cond(A[J]) # add condition, e.g. isfinite
                s += A[J]
                n += 1
            end
        end
        out[I] = s/n
    end
    return out
end

recurse(f, x, n) = n>0 ? f(recurse(f, x, n-1)) : x


"bilinear interpolation of y(x) between (x1,y1) and (x2,y2)"
bilinear(x1,x2, y1,y2, x) = ( x2 == x1 ? y1 : y1 + (y2-y1) * (x-x1) / (x2-x1) )

"returns a function that interpolates y(x) from vectors X, Y with X ascending"
function interp_extrap_ascending( X::AbstractVector{<:Real}, Y::AbstractVector{<:Real} )
    function itp(x)
        # strictly NaN true out-of-bounds inputs
        # (x < X[1] || x > X[end]) && return NaN
        # search for ascending vector order
        j = clamp(searchsortedfirst(X, x), 2, length(X)) # clamped to data intervals [2, N]
        bilinear(X[j-1],X[j], Y[j-1],Y[j], x)
    end
    return itp
end

"""
filter rfv_nrm with moving average filter.
extrapolate at bottom, and set top cloud fraction to 0 (avoids undershooting 0 cloud fraction)
"""
function filt_rfv(rfv_nrm; n=3, m=5, zerotop=false)
    mp =  m÷2 * n
    idx = clamp.(1-mp:lastindex(rfv_nrm)+mp, 1,lastindex(rfv_nrm)) # pad ends
    good(x) = !ismissing(x) && isfinite(x)
    flt(x) = recurse(x->moving_average(x, m, good), x, n) # nx moving average filter
    # extrapolate at bottom
    r = interp_extrap_ascending( eachindex(rfv_nrm), rfv_nrm ).(1-mp:lastindex(rfv_nrm)+mp) # [mp+1:lastindex(rfv_nrm)+mp]
    # set last value to last data value = 0
    r[lastindex(rfv_nrm)+mp:end] .= rfv_nrm[end]
    fltidx = mp+1:lastindex(rfv_nrm)+mp
    y = flt(r)[fltidx]
    zerotop && (y[end] = 0.0) # set top cloud fraction to 0
    return y
end

# test
x = (x -> x + mod(x,3)).(1:400)
filt_rfv(x)
