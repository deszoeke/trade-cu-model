import numpy as np


def moving_average(A, m, cond=None):
    """m-point moving average of A"""
    if cond is None:
        cond = lambda x: True
    A = np.asarray(A, dtype=float)
    if np.isscalar(m):
        m = (m,)
    half = tuple(mi // 2 for mi in m)
    out = np.empty_like(A, dtype=float)
    it = np.nditer(A, flags=["multi_index"])
    while not it.finished:
        idx = it.multi_index
        slices = tuple(
            slice(max(0, idx[d] - half[d]), min(A.shape[d], idx[d] + half[d] + 1))
            for d in range(A.ndim)
        )
        patch = A[slices]
        mask = np.vectorize(cond)(patch)
        vals = patch[mask]
        out[idx] = vals.mean() if vals.size > 0 else np.nan
        it.iternext()
    return out


def recurse(f, x, n):
    return f(recurse(f, x, n - 1)) if n > 0 else x


def bilinear(x1, x2, y1, y2, x):
    """bilinear interpolation of y(x) between (x1,y1) and (x2,y2)"""
    if x2 == x1:
        return y1
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)


def interp_extrap_ascending(X, Y):
    """returns a function that interpolates y(x) from vectors X, Y with X ascending"""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)

    def itp(x):
        j = np.searchsorted(X, x, side="left")
        j = np.clip(j, 1, len(X) - 1)
        return bilinear(X[j - 1], X[j], Y[j - 1], Y[j], x)

    return np.vectorize(itp)


def filt_rfv(rfv_nrm, n=3, m=5, zerotop=False):
    """
    filter rfv_nrm with moving average filter.
    extrapolate at bottom, and set top cloud fraction to 0 (avoids undershooting 0 cloud fraction)
    """
    rfv_nrm = np.asarray(rfv_nrm, dtype=float)
    N = len(rfv_nrm)
    mp = (m // 2) * n

    def good(x):
        return not np.isnan(x) and np.isfinite(x)

    def flt(x):
        return recurse(lambda a: moving_average(a, m, good), x, n)

    # extrapolate at bottom using 1-based index matching Julia's eachindex
    indices = np.arange(1, N + 1, dtype=float)
    pad_indices = np.arange(1 - mp, N + mp + 1, dtype=float)
    itp = interp_extrap_ascending(indices, rfv_nrm)
    r = itp(pad_indices).astype(float)

    # set last values to last data value
    r[N + mp - 1 :] = rfv_nrm[-1]

    flt_slice = slice(mp, N + mp)
    y = flt(r)[flt_slice]

    if zerotop:
        y[-1] = 0.0
    return y


if __name__ == "__main__":
    x = np.array([i + (i % 3) for i in range(1, 401)], dtype=float)
    result = filt_rfv(x)
    print(result)
