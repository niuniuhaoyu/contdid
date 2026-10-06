# v3a_prototype.R — 原型验证协变量估计量（Prop S3 设计）
# 真实 ATT(d) = (0.5 + 0.2 * mean(x|D>0)) * d
suppressMessages(library(splines2))
df <- read.csv("D:/OpenCode/contdid/examples/reference/contdid_cov_sim.csv", header = TRUE)
pre  <- df[df$t == 0, c("id","y")]; post <- df[df$t == 1, c("id","y","d","x")]
names(pre)[2] <- "y0"; names(post)[2] <- "y1"
m <- merge(pre, post, by = "id"); m$dy <- m$y1 - m$y0
d <- m$d; x <- m$x; dy <- m$dy
xtr <- mean(x[d > 0])
dmax <- max(d)
dvals <- seq(min(d[d > 0]), max(d[d > 0]), length.out = 9)
true_att <- (0.5 + 0.2 * xtr) * dvals

# --- (1) 无条件估计（v1 口径）---
tr <- d > 0
fit_u <- lm(dy ~ d, data = data.frame(dy = dy[tr], d = d[tr]))
att_u <- predict(fit_u, newdata = data.frame(d = dvals)) - mean(dy[d == 0])

# --- (2) 含 X×D 交互的联合 sieve（推荐设计）---
B  <- bSpline(d,    degree = 1, knots = numeric(0), Boundary.knots = c(0, dmax))
Bn <- bSpline(dvals,degree = 1, knots = numeric(0), Boundary.knots = c(0, dmax))
B0 <- bSpline(0,    degree = 1, knots = numeric(0), Boundary.knots = c(0, dmax))
X  <- cbind(1, B, x, B * x)
beta <- qr.solve(X, dy)
# 系数: (Intercept), B, x, B:x  -> ATT(d) = (bB + bBx * xbar)*(B(d)-B(0))
bB  <- beta[2]; bBx <- beta[4]
att_c <- as.numeric((bB + bBx * xtr) * (Bn - as.numeric(B0)))

cat(sprintf("xbar(treated) = %.4f ; dmax = %.4f\n", xtr, dmax))
cat(sprintf("true slope = %.6f\n\n", 0.5 + 0.2 * xtr))
out <- data.frame(d = round(dvals,4), true = true_att,
                  uncond = att_u, cond = att_c,
                  bias_uncond = att_u - true_att, bias_cond = att_c - true_att)
print(round(out, 6))
cat("\nmax|bias| uncond =", round(max(abs(out$bias_uncond)), 6), "\n")
cat("max|bias| cond   =", round(max(abs(out$bias_cond)), 6), "\n")
