# dump_splines2.R — generate golden-reference B-spline basis + derivative
# for cross-checking the Mata implementation.
library(splines2)

# fixed evaluation grid and a fixed set of interior knots
x <- seq(0.1, 0.9, length.out = 9)
degree <- 3
knots <- c(0.4, 0.7)          # 2 interior knots (quantiles of U(0,1) approx)

B <- bSpline(x, degree = degree, knots = knots, intercept = FALSE)
D <- dbs(x, derivs = 1, degree = degree, knots = knots, intercept = FALSE)

cat("K (ncol, intercept=FALSE):", ncol(B), "\n")
cat("expect K = length(knots)+degree =", length(knots)+degree, "\n")

write.csv(cbind(x, B), "splines2_basis.csv", row.names = FALSE)
write.csv(cbind(x, D), "splines2_deriv.csv",  row.names = FALSE)

cat("basis (first 3 rows):\n")
print(round(B[1:3, ], 6))
cat("deriv (first 3 rows):\n")
print(round(D[1:3, ], 6))
