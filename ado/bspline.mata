mata:
// B-spline basis/derivative (Cox-de Boor) + continuous-DiD fit, matching R splines2/contdid.

real colvector _extknots(real scalar p, real colvector interior, real scalar lo, real scalar hi)
{
    return(J(p+1, 1, lo) \ interior \ J(p+1, 1, hi))
}

real rowvector _bsp_full(real scalar x, real scalar p, real colvector t)
{
    real scalar L, Kfull, i, k, denom1, denom2
    real rowvector B, Bnew
    L = rows(t)
    Kfull = L - p - 1
    if (x == t[L]) {
        B = J(1, Kfull, 0)
        B[Kfull] = 1
        return(B)
    }
    B = J(1, Kfull + p, 0)
    for (i = 1; i <= Kfull + p; i++) {
        if (t[i] <= x & x < t[i+1]) B[i] = 1
    }
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

real rowvector _bsp_deriv_full(real scalar x, real scalar p, real colvector t)
{
    real scalar L, Kfull, i, denom1, denom2, val
    real rowvector Blow, D
    L = rows(t)
    Kfull = L - p - 1
    if (p == 0) return(J(1, Kfull, 0))
    if (x == t[L]) {
        D = J(1, Kfull, 0)
        val = p / (t[L-1] - t[L-p-1])
        D[Kfull-1] = -val
        D[Kfull]   =  val
        return(D)
    }
    Blow = _bsp_full(x, p-1, t)
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
        out[i, .] = row[| 2 \ Kfull |]
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

real colvector _quantile_knots(real colvector dt, real scalar nk)
{
    real colvector s, k
    real scalar n, i
    s = sort(dt, 1)
    n = rows(s)
    k = J(nk, 1, 0)
    for (i = 1; i <= nk; i++) {
        k[i] = s[round(n * i / (nk + 1))]
    }
    return(k)
}

void _contdid_fit(real colvector d, real colvector dy, real scalar degree,
                  real colvector knots, real colvector deval,
                  real colvector att, real colvector acrt)
{
    real matrix Bt, Xt, Be, Xe, Bd
    real colvector dt, dyt, beta, m0
    dt = select(d, d :> 0)
    dyt = select(dy, d :> 0)
    m0 = mean(select(dy, d :== 0))
    Bt = bspline_basis(dt, degree, knots)
    Xt = J(rows(Bt), 1, 1), Bt
    beta = qrsolve(Xt, dyt)
    Be = bspline_basis(deval, degree, knots)
    Xe = J(rows(Be), 1, 1), Be
    att = Xe * beta :- m0
    Bd = bspline_deriv(deval, degree, knots)
    acrt = Bd * beta[2..rows(beta)]
}

void _contdid_run(string scalar dvname, string scalar dyvname, real scalar degree,
                  string scalar knotmatname, real scalar nk, real scalar npoints,
                  real scalar dmin, real scalar dmax,
                  string scalar attname, string scalar acrtname)
{
    real colvector d, dy, dt, knots, deval, att, acrt
    real matrix km
    real scalar i
    d  = st_data(., dvname)
    dy = st_data(., dyvname)
    if (nk > 0) {
        dt = select(d, d :> 0)
        knots = _quantile_knots(dt, nk)
    }
    else if (knotmatname != "") {
        km = st_matrix(knotmatname)
        knots = km'
    }
    else {
        knots = J(0, 1, 0)
    }
    deval = J(npoints, 1, 0)
    for (i = 1; i <= npoints; i++) {
        deval[i] = dmin + (dmax - dmin) * (i - 1) / (npoints - 1)
    }
    _contdid_fit(d, dy, degree, knots, deval, att, acrt)
    st_matrix(attname, att')
    st_matrix(acrtname, acrt')
}
end
