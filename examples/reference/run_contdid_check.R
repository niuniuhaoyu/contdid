# run_contdid_check.R — independent base-R reimplementation of the contdid
# linear-in-dose estimator, to cross-check the Stata .ado numerically.
# Matches R contdid's default (num_knots=0, degree=1):
#   ATT(d) = (a + b*d) - mean(dy | D=0),  where a,b from OLS of dy on d among D>0.

df <- read.csv("contdid_sim.csv", header = TRUE)

pre  <- df[df$t == 0, c("id", "y")]
post <- df[df$t == 1, c("id", "y", "d")]
names(pre)[2]  <- "y0"
names(post)[2] <- "y1"
m <- merge(pre, post, by = "id")
m$dy <- m$y1 - m$y0

tr  <- m[m$d > 0, ]
fit <- lm(dy ~ d, data = tr)
a   <- coef(fit)[["(Intercept)"]]
b   <- coef(fit)[["d"]]
m0  <- mean(m$dy[m$d == 0])
dmin <- min(m$d[m$d > 0])
dmax <- max(m$d[m$d > 0])

cat(sprintf("a=%.8f b=%.8f m0=%.8f dmin=%.8f dmax=%.8f\n", a, b, m0, dmin, dmax))
for (k in 1:5) {
  dk  <- dmin + (dmax - dmin) * (k - 1) / 4
  att <- a + b * dk - m0
  cat(sprintf("k=%d d=%.4f ATT=%.6f\n", k, dk, att))
}
