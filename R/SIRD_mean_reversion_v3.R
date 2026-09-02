############################################################
# SAME PLOT:
# deterministic + stochastic mean-reversion together
############################################################

rm(list = ls())

#############################
# 1) PARAMETERS
#############################
beta0 <- 0.15
theta <- 0.40
xi    <- 0.035
T_end <- 30
dt    <- 0.01

beta_bar_vec <- c(0.22, 0.35, 0.50)

# time grid
t <- seq(0, T_end, by = dt)

#############################
# 2) FUNCTIONS
#############################

# deterministic path
beta_path <- function(t, beta_bar, beta0, theta) {
  beta_bar + (beta0 - beta_bar) * exp(-theta * t)
}

# stochastic path
simulate_beta_path <- function(beta0, beta_bar, theta, xi, T_end, dt) {
  n <- floor(T_end / dt)
  t <- seq(0, T_end, by = dt)
  beta <- numeric(length(t))
  beta[1] <- beta0
  
  for (i in 1:n) {
    z <- rnorm(1)
    beta[i + 1] <- beta[i] + theta * (beta_bar - beta[i]) * dt + xi * sqrt(dt) * z
    beta[i + 1] <- max(beta[i + 1], 0)
  }
  
  beta
}

#############################
# 3) BUILD PATHS
#############################

# deterministic paths
beta_det_mat <- sapply(beta_bar_vec, function(bb) {
  beta_path(t, bb, beta0, theta)
})
beta_det_mat <- as.matrix(beta_det_mat)

# stochastic paths
set.seed(123)
beta_sto_mat <- sapply(beta_bar_vec, function(bb) {
  simulate_beta_path(beta0, bb, theta, xi, T_end, dt)
})
beta_sto_mat <- as.matrix(beta_sto_mat)

#############################
# 4) PLOT
#############################

# choose line colors automatically
cols <- cols <- c("steelblue", "firebrick", "darkgreen")

par(mar = c(5, 5, 4, 2))

# empty frame first
plot(
  t, beta_det_mat[, 1],
  type = "n",
  xlab = "Time",
  ylab = expression(beta(t)),
  main = expression(paste("Deterministic and Stochastic Mean-Reversion of ", beta(t))),
  cex.main = 1.0,
  ylim = range(c(beta_det_mat, beta_sto_mat))
)

# stochastic paths first (thin)
for (j in seq_along(beta_bar_vec)) {
  lines(
    t, beta_sto_mat[, j],
    lwd = 1.2,
    lty = 1,
    col = cols[j]
  )
}

# deterministic paths on top (thick)
for (j in seq_along(beta_bar_vec)) {
  lines(
    t, beta_det_mat[, j],
    lwd = 3,
    lty = 2,
    col = cols[j]
  )
}

# horizontal lines for long-run means
for (j in seq_along(beta_bar_vec)) {
  abline(h = beta_bar_vec[j], lty = 3, col = cols[j])
}

# initial value line and point
abline(h = beta0, lty = 2, col = "black")
points(0, beta0, pch = 16, cex = 1.2)

#############################
# 5) LEGEND
#############################
legend(
  "topleft",
  legend = c(
    paste0("stochastic, beta_bar = ", beta_bar_vec),
    paste0("deterministic, beta_bar = ", beta_bar_vec),
    paste0("beta0 = ", beta0)
  ),
  col = c(cols, cols, "black"),
  lty = c(rep(1, length(beta_bar_vec)),
          rep(3, length(beta_bar_vec)),
          2),
  lwd = c(rep(1.2, length(beta_bar_vec)),
          rep(3, length(beta_bar_vec)),
          1.5),
  pch = c(rep(NA, 2 * length(beta_bar_vec)), 16),
  bty = "n",
  cex = 0.85
)