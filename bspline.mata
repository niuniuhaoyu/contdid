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

// uniform confidence band via multiplier bootstrap (iid units).
void _contdid_cband(string scalar dvname, string scalar dyvname, real scalar degree,
                    string scalar knotmatname, real scalar nk, real scalar npoints,
                    real scalar dmin, real scalar dmax, real scalar B, real scalar level)
{
    real matrix Bt, Xt, Be, Bd, Xe, bread, IFb, psi_a, psi_c, cbatt, cbacrt
    real colvector d, dy, dt, dyt, knots, deval, beta, e, m0, xi, G
    real colvector sigma_a, sigma_c, supb_a, supb_c, att, acrt, se_a, se_c
    real scalar n, nt, n0, p0, i, j, b, K, np, alpha, qa, qc
    real colvector tridx

    d  = st_data(., dvname)
    dy = st_data(., dyvname)
    if (nk > 0) {
        dt = select(d, d :> 0)
        knots = _quantile_knots(dt, nk)
    }
    else if (knotmatname != "") {
        knots = st_matrix(knotmatname)'
    }
    else {
        knots = J(0, 1, 0)
    }
    deval = J(npoints, 1, 0)
    for (i = 1; i <= npoints; i++) {
        deval[i] = dmin + (dmax - dmin) * (i - 1) / (npoints - 1)
    }

    n = rows(d)
    dt = select(d, d :> 0)
    dyt = select(dy, d :> 0)
    nt = rows(dt)
    n0 = n - nt
    p0 = n0 / n
    m0 = mean(select(dy, d :== 0))

    Bt = bspline_basis(dt, degree, knots)
    Xt = J(nt, 1, 1), Bt
    K = cols(Bt)
    beta = qrsolve(Xt, dyt)
    e = dyt - Xt * beta
    bread = invsym(Xt'Xt / nt)

    Be = bspline_basis(deval, degree, knots)
    Xe = J(rows(Be), 1, 1), Be
    Bd = bspline_deriv(deval, degree, knots)
    np = rows(deval)

    // point estimates
    att = Xe * beta :- m0
    acrt = Bd * beta[2..rows(beta)]

    IFb = J(nt, K+1, 0)
    for (i = 1; i <= nt; i++) {
        IFb[i, .] = (e[i] * (bread * Xt[i,.]'))'
    }
    tridx = selectindex(d :> 0)

    psi_a = J(n, np, 0)
    psi_c = J(n, np, 0)
    for (j = 1; j <= np; j++) {
        for (i = 1; i <= nt; i++) {
            psi_a[tridx[i], j] = Xe[j, .] * IFb[i, .]'
            psi_c[tridx[i], j] = Bd[j, .] * IFb[i, 2..cols(IFb)]'
        }
        for (i = 1; i <= n; i++) {
            if (d[i] == 0) psi_a[i, j] = -(dy[i] - m0) / p0
        }
    }

    sigma_a = J(1, np, 0)
    sigma_c = J(1, np, 0)
    for (j = 1; j <= np; j++) {
        sigma_a[1, j] = sqrt(mean(psi_a[., j] :^ 2))
        sigma_c[1, j] = sqrt(mean(psi_c[., j] :^ 2))
    }

    supb_a = J(B, 1, 0)
    supb_c = J(B, 1, 0)
    for (b = 1; b <= B; b++) {
        xi = (runiform(n, 1) :> 0.5) :* 2 :- 1   // Rademacher +/-1 multipliers
        G = (xi' * psi_a) / sqrt(n)
        supb_a[b] = max(abs(G) :/ sigma_a)
        G = (xi' * psi_c) / sqrt(n)
        supb_c[b] = max(abs(G) :/ sigma_c)
    }

    alpha = 1 - level / 100
    supb_a = sort(supb_a, 1)
    supb_c = sort(supb_c, 1)
    qa = ceil((1 - alpha) * B)
    st_numscalar("crit_att", supb_a[qa])
    st_numscalar("crit_acrt", supb_c[qa])

    se_a = sigma_a' / sqrt(n)
    se_c = sigma_c' / sqrt(n)

    // uniform band matrices: [d, est, cb_lb, cb_ub]
    cbatt = J(np, 4, 0)
    cbacrt = J(np, 4, 0)
    for (j = 1; j <= np; j++) {
        cbatt[j, .]  = (deval[j], att[j], att[j] - supb_a[qa] * se_a[j], att[j] + supb_a[qa] * se_a[j])
        cbacrt[j, .] = (deval[j], acrt[j], acrt[j] - supb_c[qa] * se_c[j], acrt[j] + supb_c[qa] * se_c[j])
    }
    st_matrix("cb_att", cbatt)
    st_matrix("cb_acrt", cbacrt)
}

// staggered adoption: group-time dose-response + aggregation.
// d (n), g (n, 0=never), Y (n x T wide, cols in time order), tvals (sorted time values),
// deval (eval doses). returns att, acrt aggregated over post-treatment (g,t).
void _contdid_stag(real colvector d, real colvector g, real matrix Y, real colvector tvals,
                   real scalar degree, real colvector knots, real colvector deval,
                   real colvector att, real colvector acrt)
{
    real colvector cohorts, dy, dact, gtatt, gtacrt
    real scalar n, T, np, c, ti, ti2, base, nc, wsum
    real colvector subidx

    n = rows(d)
    T = rows(tvals)
    np = rows(deval)
    cohorts = uniqrows(select(g, g :> 0))

    att = J(np, 1, 0)
    acrt = J(np, 1, 0)
    wsum = 0

    for (c = 1; c <= rows(cohorts); c++) {
        ti = selectindex(tvals :== cohorts[c])
        if (ti == 0) continue
        base = ti - 1
        if (base < 1) continue
        for (ti2 = ti; ti2 <= T; ti2++) {
            subidx = selectindex((g :== cohorts[c]) :| (g :> tvals[ti2]) :| (g :== 0))
            dy = Y[subidx, ti2] - Y[subidx, base]
            dact = d[subidx] :* (g[subidx] :== cohorts[c])
            _contdid_fit(dact, dy, degree, knots, deval, gtatt, gtacrt)
            nc = sum(g :== cohorts[c])
            att = att + nc * gtatt
            acrt = acrt + nc * gtacrt
            wsum = wsum + nc
        }
    }
    att = att / wsum
    acrt = acrt / wsum
}

void _contdid_stag_run(string scalar dvname, string scalar gvname, string scalar yvars,
                       string scalar tvalsname, real scalar degree,
                       string scalar knotmatname, real scalar nk, real scalar npoints,
                       real scalar dmin, real scalar dmax,
                       string scalar attname, string scalar acrtname)
{
    real colvector d, g, dt, knots, deval, att, acrt, tvals
    real matrix Y, km
    real scalar i
    d = st_data(., dvname)
    g = st_data(., gvname)
    Y = st_data(., yvars)
    tvals = st_matrix(tvalsname)'
    if (nk > 0) {
        dt = select(d, d :> 0)
        knots = _quantile_knots(dt, nk)
    }
    else if (knotmatname != "") {
        knots = st_matrix(knotmatname)'
    }
    else {
        knots = J(0, 1, 0)
    }
    deval = J(npoints, 1, 0)
    for (i = 1; i <= npoints; i++) {
        deval[i] = dmin + (dmax - dmin) * (i - 1) / (npoints - 1)
    }
    _contdid_stag(d, g, Y, tvals, degree, knots, deval, att, acrt)
    st_matrix(attname, att')
    st_matrix(acrtname, acrt')
}

// ---------------------------------------------------------------------------
// v3a: covariates (Conditional Strong Parallel Trends; CGBS SI.3, Prop S3)
//   ATT_x(d) = E[dY | X=x, D=d] - E[dY | X=x, D=0]
//   ATT(d)   = E_X[ ATT_x(d) | D>0 ]  (aggregate over treated X)
// Model: dY = b0 + B(D)'b + X'g + sum_j (B(D)*X_j)'d_j + e
// Requires explicit boundary knots [0, dmax] so that B(0) is well-defined.
// ---------------------------------------------------------------------------

real matrix bspline_basis_bk(real colvector x, real scalar p, real colvector interior,
                             real scalar lo, real scalar hi)
{
    real scalar n, i, K, Kfull
    real colvector t
    real matrix out
    real rowvector row
    n = rows(x)
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

real matrix bspline_deriv_bk(real colvector x, real scalar p, real colvector interior,
                             real scalar lo, real scalar hi)
{
    real scalar n, i, K, Kfull
    real colvector t
    real matrix out
    real rowvector row
    n = rows(x)
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

void _contdid_fit_cov(real colvector d, real colvector dy, real matrix Xm,
                      real scalar degree, real colvector knots,
                      real colvector deval, real colvector att, real colvector acrt)
{
    real scalar n, J, K, dmax, np, i, j, off
    real matrix B, W, beta, Bn, Bd
    real rowvector B0, g, xbar

    n = rows(d); J = cols(Xm); dmax = max(d)
    B = bspline_basis_bk(d, degree, knots, 0, dmax)
    K = cols(B)
    W = (J(n,1,1), B, Xm)
    for (j = 1; j <= J; j++) W = W, (B :* Xm[.,j])
    beta = qrsolve(W, dy)

    // treatment-group covariate means
    xbar = J(1, J, 0)
    for (j = 1; j <= J; j++) xbar[1,j] = mean(select(Xm[.,j], d :> 0))

    // g[k] = beta_B[k] + sum_j delta_{j,k} * xbar_j
    g = J(1, K, 0)
    for (i = 1; i <= K; i++) g[1,i] = beta[1+i]
    for (j = 1; j <= J; j++) {
        off = 2 + K + J + (j - 1) * K
        for (i = 1; i <= K; i++) g[1,i] = g[1,i] + beta[off + i - 1] * xbar[1,j]
    }

    np = rows(deval)
    Bn = bspline_basis_bk(deval, degree, knots, 0, dmax)
    B0 = bspline_basis_bk(J(1,1,0), degree, knots, 0, dmax)
    Bd = bspline_deriv_bk(deval, degree, knots, 0, dmax)
    att = J(np,1,0); acrt = J(np,1,0)
    for (i = 1; i <= np; i++) {
        att[i]  = (Bn[i,.] :- B0) * g'
        acrt[i] = Bd[i,.] * g'
    }
}

void _contdid_run_cov(string scalar dvname, string scalar dyvname, string scalar xnames,
                      real scalar degree, string scalar knotmatname, real scalar nk,
                      real scalar npoints, real scalar dmin, real scalar dmax,
                      string scalar attname, string scalar acrtname)
{
    real colvector d, dy, dt, knots, deval, att, acrt
    real matrix Xm, km
    real scalar i
    d  = st_data(., dvname)
    dy = st_data(., dyvname)
    Xm = st_data(., xnames)
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
    _contdid_fit_cov(d, dy, Xm, degree, knots, deval, att, acrt)
    st_matrix(attname, att')
    st_matrix(acrtname, acrt')
}

// v3a cband: influence-function multiplier bootstrap with covariates.
// IF includes IF(theta) (OLS among all units) and IF(xbar) (treated covariate means).
void _contdid_cband_cov(string scalar dvname, string scalar dyvname, string scalar xnames,
                        real scalar degree, string scalar knotmatname, real scalar nk,
                        real scalar npoints, real scalar dmin, real scalar dmax,
                        real scalar Bdraws, real scalar level)
{
    real scalar n, J, K, dmaxv, np, i, j, jj, k, off, b, alpha, qa, qc, va, vc, sa, sc, p0
    real colvector d, dy, dt, knots, deval, e, supb_a, supb_c, sigma_a, sigma_c, xi, G, se_a, se_c, att, acrt
    real rowvector  B0, cA, cC, g, xbar
    real matrix Xm, B, W, theta, bread, IFtheta, Bn, Bd, km, IFxbar, psi_a, psi_c, cbatt, cbacrt
    real colvector rows_a, rows_c

    d  = st_data(., dvname)
    dy = st_data(., dyvname)
    Xm = st_data(., xnames)
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

    n = rows(d); J = cols(Xm); dmaxv = max(d)
    B = bspline_basis_bk(d, degree, knots, 0, dmaxv)
    K = cols(B)
    W = (J(n,1,1), B, Xm)
    for (jj = 1; jj <= J; jj++) W = W, (B :* Xm[.,jj])
    theta = qrsolve(W, dy)
    e = dy - W * theta
    bread = invsym(W'W / n)
    IFtheta = J(n, cols(W), 0)
    for (i = 1; i <= n; i++) IFtheta[i,.] = (e[i] * (bread * W[i,.]'))'

    xbar = J(1, J, 0)
    for (jj = 1; jj <= J; jj++) xbar[1,jj] = mean(select(Xm[.,jj], d :> 0))
    p0 = mean(d :> 0)
    IFxbar = J(n, J, 0)
    for (jj = 1; jj <= J; jj++) {
        for (i = 1; i <= n; i++) {
            if (d[i] > 0) IFxbar[i,jj] = (Xm[i,jj] - xbar[1,jj]) / p0
        }
    }

    deval = J(npoints, 1, 0)
    for (i = 1; i <= npoints; i++) deval[i] = dmin + (dmax - dmin) * (i - 1) / (npoints - 1)
    np = npoints
    Bn = bspline_basis_bk(deval, degree, knots, 0, dmaxv)
    B0 = bspline_basis_bk(J(1,1,0), degree, knots, 0, dmaxv)
    Bd = bspline_deriv_bk(deval, degree, knots, 0, dmaxv)

    psi_a = J(n, np, 0)
    psi_c = J(n, np, 0)
    for (j = 1; j <= np; j++) {
        cA = Bn[j,.] :- B0
        cC = Bd[j,.]
        for (i = 1; i <= n; i++) {
            va = 0; vc = 0
            for (k = 1; k <= K; k++) {
                va = va + cA[k] * IFtheta[i, 1+k]
                vc = vc + cC[k] * IFtheta[i, 1+k]
            }
            for (jj = 1; jj <= J; jj++) {
                off = 2 + K + J + (jj - 1) * K
                for (k = 1; k <= K; k++) {
                    va = va + xbar[1,jj] * cA[k] * IFtheta[i, off + k - 1]
                    vc = vc + xbar[1,jj] * cC[k] * IFtheta[i, off + k - 1]
                }
            }
            for (jj = 1; jj <= J; jj++) {
                off = 2 + K + J + (jj - 1) * K
                sa = 0; sc = 0
                for (k = 1; k <= K; k++) {
                    sa = sa + theta[off + k - 1] * cA[k]
                    sc = sc + theta[off + k - 1] * cC[k]
                }
                va = va + sa * IFxbar[i,jj]
                vc = vc + sc * IFxbar[i,jj]
            }
            psi_a[i,j] = va
            psi_c[i,j] = vc
        }
    }

    sigma_a = J(1, np, 0)
    sigma_c = J(1, np, 0)
    for (j = 1; j <= np; j++) {
        sigma_a[1,j] = sqrt(mean(psi_a[.,j] :^ 2))
        sigma_c[1,j] = sqrt(mean(psi_c[.,j] :^ 2))
    }

    supb_a = J(Bdraws, 1, 0)
    supb_c = J(Bdraws, 1, 0)
    for (b = 1; b <= Bdraws; b++) {
        xi = (runiform(n, 1) :> 0.5) :* 2 :- 1
        G = (xi' * psi_a) / sqrt(n)
        supb_a[b] = max(abs(G) :/ sigma_a)
        G = (xi' * psi_c) / sqrt(n)
        supb_c[b] = max(abs(G) :/ sigma_c)
    }
    alpha = 1 - level / 100
    supb_a = sort(supb_a, 1)
    supb_c = sort(supb_c, 1)
    qa = ceil((1 - alpha) * Bdraws)
    qc = qa
    st_numscalar("crit_att", supb_a[qa])
    st_numscalar("crit_acrt", supb_c[qc])

    se_a = sigma_a' / sqrt(n)
    se_c = sigma_c' / sqrt(n)

    att = J(np, 1, 0)
    acrt = J(np, 1, 0)
    for (j = 1; j <= np; j++) {
        g = J(1, K, 0)
        for (k = 1; k <= K; k++) g[1,k] = theta[1+k]
        for (jj = 1; jj <= J; jj++) {
            off = 2 + K + J + (jj - 1) * K
            for (k = 1; k <= K; k++) g[1,k] = g[1,k] + theta[off + k - 1] * xbar[1,jj]
        }
        att[j]  = (Bn[j,.] :- B0) * g'
        acrt[j] = Bd[j,.] * g'
    }

    cbatt = J(np, 4, 0)
    cbacrt = J(np, 4, 0)
    for (j = 1; j <= np; j++) {
        cbatt[j,.]  = (deval[j], att[j], att[j] - supb_a[qa] * se_a[j], att[j] + supb_a[qa] * se_a[j])
        cbacrt[j,.] = (deval[j], acrt[j], acrt[j] - supb_c[qc] * se_c[j], acrt[j] + supb_c[qc] * se_c[j])
    }
    st_matrix("cb_att", cbatt)
    st_matrix("cb_acrt", cbacrt)
}
end
