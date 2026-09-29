mata:
// B-spline basis and derivative (Cox-de Boor), matching R splines2 intercept=FALSE.
// public: bspline_basis(x, degree, interior_knots) -> n x K, K = nk + degree
//         bspline_deriv(x, degree, interior_knots) -> n x K

real colvector _extknots(real scalar p, real colvector interior, real scalar lo, real scalar hi)
{
    return(J(p+1, 1, lo) \ interior \ J(p+1, 1, hi))
}

// full basis (1 x Kfull) at scalar x; Kfull = degree + 1 + nk
real rowvector _bsp_full(real scalar x, real scalar p, real colvector t)
{
    real scalar L, Kfull, i, k, denom1, denom2
    real rowvector B, Bnew
    L = rows(t)
    Kfull = L - p - 1
    B = J(1, Kfull + p, 0)
    for (i = 1; i <= Kfull + p - 1; i++) {
        if (t[i] <= x & x < t[i+1]) B[i] = 1
    }
    i = Kfull + p                                  // right boundary: include x == t[L]
    if (t[i] <= x & x <= t[i+1]) B[i] = 1
    for (k = 1; k <= p; k++) {
        Bnew = J(1, Kfull + p - k, 0)
        for (i = 1; i <= Kfull + p - k; i++) {
            denom1 = t[i+k] - t[i]
            denom2 = t[i+k+1] - t[i+1]
            if (denom1 > 0) Bnew[i] = (x - t[i]) / denom1 * B[i]
            if (denom2 > 0) Bnew[i] = Bnew[i] + (t[i+k+1] - x) / denom2 * B[i+1]
        }
        B = Bnew
    }
    return(B)
}

// full derivative (1 x Kfull)
real rowvector _bsp_deriv_full(real scalar x, real scalar p, real colvector t)
{
    real scalar L, Kfull, i, denom1, denom2
    real rowvector Blow, D
    L = rows(t)
    Kfull = L - p - 1
    if (p == 0) return(J(1, Kfull, 0))
    Blow = _bsp_full(x, p-1, t)                    // length Kfull + 1
    D = J(1, Kfull, 0)
    for (i = 1; i <= Kfull; i++) {
        denom1 = t[i+p] - t[i]
        denom2 = t[i+p+1] - t[i+1]
        if (denom1 > 0) D[i] = D[i] + p / denom1 * Blow[i]
        if (denom2 > 0) D[i] = D[i] - p / denom2 * Blow[i+1]
    }
    return(D)
}

real matrix bspline_basis(real colvector x, real scalar p, real colvector interior)
{
    real scalar lo, hi, n, i, K, Kfull
    real colvector t
    real matrix out
    real rowvector row
    n = rows(x)
    lo = min(x); hi = max(x)
    K = rows(interior) + p
    Kfull = p + 1 + rows(interior)
    t = _extknots(p, interior, lo, hi)
    out = J(n, K, 0)
    for (i = 1; i <= n; i++) {
        row = _bsp_full(x[i], p, t)
        out[i, .] = row[| 2 \ Kfull |]             // drop first basis fn (intercept)
    }
    return(out)
}

real matrix bspline_deriv(real colvector x, real scalar p, real colvector interior)
{
    real scalar lo, hi, n, i, K, Kfull
    real colvector t
    real matrix out
    real rowvector row
    n = rows(x)
    lo = min(x); hi = max(x)
    K = rows(interior) + p
    Kfull = p + 1 + rows(interior)
    t = _extknots(p, interior, lo, hi)
    out = J(n, K, 0)
    for (i = 1; i <= n; i++) {
        row = _bsp_deriv_full(x[i], p, t)
        out[i, .] = row[| 2 \ Kfull |]
    }
    return(out)
}
end
