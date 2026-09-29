*! _test_bspline.do — verify Mata B-spline matches R splines2 golden reference
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
run "ado/bspline.mata"

mata:
x = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)'
knots = (0.4, 0.7)'
B = bspline_basis(x, 3, knots)
D = bspline_deriv(x, 3, knots)
st_matrix("mB", B)
st_matrix("mD", D)
end

* write to csv
clear
svmat mB
export delimited using "examples/reference/mata_basis.csv", replace
clear
svmat mD
export delimited using "examples/reference/mata_deriv.csv", replace
di "done"
