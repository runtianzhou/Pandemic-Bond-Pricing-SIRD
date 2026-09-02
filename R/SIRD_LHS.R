############################################################
# LATIN HYPERCUBE SENSITIVITY FOR SIRD BOND PRICE
# Latin hypercube sampling (LHS)
############################################################

library(lhs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

set.seed(123)

############################################################
# 1) USER SETTINGS
############################################################

n_lhs <- 500
M_sim <- 1500

# Fixed bond / trigger settings
PH_fixed      <- 0.1468
death_trigger <- 2500
maI_trigger   <- 250
N_pop         <- 8.3e9
I0            <- 100
R0_init       <- 0
D0            <- 0
beta0_fixed   <- 0.15

# Parameter ranges for LHS
param_ranges <- list(
  beta_bar = c(0.20, 0.40),
  gamma    = c(0.01, 0.10),
  xi       = c(0.10, 1.20),
  alpha    = c(0.01, 0.10),
  theta    = c(0.10, 1.30)
)

############################################################
# 1A) CONSOLE TABLE OF PARAMETER RANGES
############################################################

range_table <- data.frame(
  Parameter = c("beta_bar", "gamma", "xi", "alpha", "theta"),
  Lower     = c(
    param_ranges$beta_bar[1],
    param_ranges$gamma[1],
    param_ranges$xi[1],
    param_ranges$alpha[1],
    param_ranges$theta[1]
  ),
  Upper     = c(
    param_ranges$beta_bar[2],
    param_ranges$gamma[2],
    param_ranges$xi[2],
    param_ranges$alpha[2],
    param_ranges$theta[2]
  ),
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("LHS PARAMETER RANGES\n")
cat("========================================\n")
print(format(range_table, digits = 4), row.names = FALSE)

############################################################
# 2) MAP LHS(0,1) TO PARAMETER RANGES
############################################################

map_to_range <- function(u, lower, upper) {
  lower + u * (upper - lower)
}

lhs_u <- randomLHS(n_lhs, length(param_ranges))
colnames(lhs_u) <- names(param_ranges)

lhs_params <- as.data.frame(lhs_u)

for (nm in names(param_ranges)) {
  lhs_params[[nm]] <- map_to_range(
    lhs_params[[nm]],
    param_ranges[[nm]][1],
    param_ranges[[nm]][2]
  )
}

lhs_params$beta0 <- beta0_fixed

############################################################
# 3) REAL SIRD-OU MONTE CARLO PRICING CODE
############################################################

beta_bar_fn <- function(t, beta, beta0, theta) {
  beta + (beta0 - beta) * exp(-theta * t)
}

sigma_t <- function(t, theta, xi) {
  (xi / sqrt(2 * theta)) * sqrt(1 - exp(-2 * theta * t))
}

roll_mean_7_matrix <- function(X) {
  n <- nrow(X)
  M <- ncol(X)
  out <- matrix(NA_real_, nrow = n, ncol = M)
  if (n < 7) return(out)
  
  for (i in 7:n) {
    out[i, ] <- colMeans(X[(i - 6):i, , drop = FALSE])
  }
  out
}

# Placeholder bond valuation
bond_value_full_lhs <- function() {
  225
}

bond_value_triggered_lhs <- function(tau_day, T_end = 1104) {
  max(0, 225 * (tau_day / T_end))
}

simulate_sird_paths_lhs <- function(
    T_end = 1104,
    M = 1000,
    N_pop = 8.3e9,
    I0_count = 100,
    R0_count = 0,
    D0_count = 0,
    beta,
    beta0,
    theta,
    xi,
    alpha,
    gamma,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  days <- 0:T_end
  n    <- length(days)
  
  I0 <- I0_count / N_pop
  R0 <- R0_count / N_pop
  D0 <- D0_count / N_pop
  S0 <- 1 - I0 - R0 - D0
  
  if (S0 <= 0) stop("Initial S0 must be positive.")
  
  S <- matrix(NA_real_, nrow = n, ncol = M)
  I <- matrix(NA_real_, nrow = n, ncol = M)
  R <- matrix(NA_real_, nrow = n, ncol = M)
  D <- matrix(NA_real_, nrow = n, ncol = M)
  
  S[1, ] <- S0
  I[1, ] <- I0
  R[1, ] <- R0
  D[1, ] <- D0
  
  for (k in 1:T_end) {
    tk  <- days[k]
    btk <- beta_bar_fn(tk, beta, beta0, theta)
    sig <- sigma_t(tk, theta, xi)
    
    Zk  <- rnorm(M)
    dBk <- Zk
    
    S_prev <- S[k, ]
    I_prev <- I[k, ]
    
    eta_I <- sig * S_prev
    drift_logI <- btk * S_prev - (alpha + gamma) - 0.5 * eta_I^2
    I_new <- I_prev * exp(drift_logI + eta_I * dBk)
    
    eta_S <- -sig * I_prev
    drift_logS <- -btk * I_prev - 0.5 * eta_S^2
    S_new <- S_prev * exp(drift_logS + eta_S * dBk)
    
    S_new <- pmax(S_new, 0)
    I_new <- pmax(I_new, 0)
    
    total_SI <- S_new + I_new
    bad_SI <- total_SI > 1
    if (any(bad_SI)) {
      S_new[bad_SI] <- S_new[bad_SI] / total_SI[bad_SI]
      I_new[bad_SI] <- I_new[bad_SI] / total_SI[bad_SI]
    }
    
    remaining_RD <- 1 - S_new - I_new
    R_new <- (alpha / (alpha + gamma)) * remaining_RD
    D_new <- (gamma / (alpha + gamma)) * remaining_RD
    
    R_new <- pmax(R_new, 0)
    D_new <- pmax(D_new, 0)
    
    S[k + 1, ] <- S_new
    I[k + 1, ] <- I_new
    R[k + 1, ] <- R_new
    D[k + 1, ] <- D_new
  }
  
  list(days = days, S = S, I = I, R = R, D = D)
}

build_trigger_series_lhs <- function(sim, N_pop) {
  I_count <- sim$I * N_pop
  R_count <- sim$R * N_pop
  D_total <- sim$D * N_pop
  
  C_total <- I_count + R_count + D_total
  
  n <- nrow(C_total)
  M <- ncol(C_total)
  
  daily_new_inf <- matrix(NA_real_, nrow = n, ncol = M)
  daily_new_inf[1, ] <- NA_real_
  
  if (n >= 2) {
    daily_new_inf[2:n, ] <- pmax(
      0,
      C_total[2:n, , drop = FALSE] - C_total[1:(n - 1), , drop = FALSE]
    )
  }
  
  maI <- roll_mean_7_matrix(daily_new_inf)
  
  list(
    I_count       = I_count,
    R_count       = R_count,
    D_total       = D_total,
    C_total       = C_total,
    daily_new_inf = daily_new_inf,
    maI           = maI
  )
}

calc_trigger_tau_lhs <- function(days, D_total, maI,
                                 death_trigger = 2500,
                                 maI_trigger = 250) {
  M <- ncol(D_total)
  tau <- rep(NA_integer_, M)
  
  for (j in 1:M) {
    cond <- (!is.na(D_total[, j])) &
      (!is.na(maI[, j])) &
      (D_total[, j] >= death_trigger) &
      (maI[, j]     >= maI_trigger)
    
    hit <- which(cond)
    if (length(hit) > 0) {
      tau[j] <- days[hit[1]]
    }
  }
  
  tau
}

price_bond_sird <- function(beta_bar, gamma, xi, alpha, theta,
                            beta0 = 0.15,
                            PH = 0.1468,
                            death_trigger = 2500,
                            maI_trigger = 250,
                            N_pop = 8.3e9,
                            I0 = 100, R0_init = 0, D0 = 0,
                            M = 1000) {
  
  sim <- simulate_sird_paths_lhs(
    T_end     = 1104,
    M         = M,
    N_pop     = N_pop,
    I0_count  = I0,
    R0_count  = R0_init,
    D0_count  = D0,
    beta      = beta_bar,
    beta0     = beta0,
    theta     = theta,
    xi        = xi,
    alpha     = alpha,
    gamma     = gamma
  )
  
  series <- build_trigger_series_lhs(sim, N_pop = N_pop)
  
  tau <- calc_trigger_tau_lhs(
    days          = sim$days,
    D_total       = series$D_total,
    maI           = series$maI,
    death_trigger = death_trigger,
    maI_trigger   = maI_trigger
  )
  
  triggered <- which(is.finite(tau))
  p_cond <- length(triggered) / M
  PC <- PH * p_cond
  
  full_price <- bond_value_full_lhs()
  
  if (length(triggered) > 0) {
    partial_vals <- vapply(
      tau[triggered],
      function(x) bond_value_triggered_lhs(x, T_end = 1104),
      FUN.VALUE = numeric(1)
    )
    trigger_price <- mean(partial_vals)
  } else {
    trigger_price <- full_price
  }
  
  bond_price <- trigger_price * PC + full_price * (1 - PC)
  
  list(
    PC = PC,
    bond_price = bond_price
  )
}

############################################################
# 4) RUN LHS EXPERIMENT
############################################################

results_list <- vector("list", nrow(lhs_params))

for (i in seq_len(nrow(lhs_params))) {
  pars <- lhs_params[i, ]
  
  out <- price_bond_sird(
    beta_bar = pars$beta_bar,
    gamma    = pars$gamma,
    xi       = pars$xi,
    alpha    = pars$alpha,
    theta    = pars$theta,
    beta0    = pars$beta0,
    PH       = PH_fixed,
    death_trigger = death_trigger,
    maI_trigger   = maI_trigger,
    N_pop    = N_pop,
    I0       = I0,
    R0_init  = R0_init,
    D0       = D0,
    M        = M_sim
  )
  
  results_list[[i]] <- data.frame(
    beta_bar   = pars$beta_bar,
    gamma      = pars$gamma,
    xi         = pars$xi,
    alpha      = pars$alpha,
    theta      = pars$theta,
    beta0      = pars$beta0,
    PC         = out$PC,
    bond_price = out$bond_price
  )
}

results <- bind_rows(results_list)

############################################################
# 5) GLOBAL SENSITIVITY VIA STANDARDIZED REGRESSION
############################################################

fit <- lm(
  scale(bond_price) ~ scale(beta_bar) + scale(gamma) +
    scale(xi) + scale(alpha) + scale(theta),
  data = results
)

coef_df <- data.frame(
  parameter = names(coef(fit))[-1],
  effect    = as.numeric(coef(fit)[-1])
)

coef_df$parameter <- gsub("scale\\(|\\)", "", coef_df$parameter)
coef_df <- coef_df %>% mutate(abs_effect = abs(effect))

############################################################
# 6) GRAPH 1: BAR PLOT
############################################################

p1 <- ggplot(coef_df, aes(x = reorder(parameter, abs_effect), y = effect, fill = effect)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#B22222",
    mid = "#F2F2F2",
    high = "#2166AC",
    midpoint = 0,
    name = "Effect on\nbond price"
  ) +
  scale_x_discrete(
    labels = c(
      beta_bar = expression(bar(beta)),
      gamma    = expression(gamma),
      xi       = expression(xi),
      alpha    = expression(alpha),
      theta    = expression(theta)
    )
  ) +
  labs(
    title = "Global sensitivity of bond price under LHS",
    x = "Parameter",
    y = "Effect on bond price"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.title = element_text(face = "bold")
  )

print(p1)

############################################################
# 7) COMBINED FACET SCATTER PLOT
############################################################

plot_data <- results %>%
  select(bond_price, PC, beta_bar, gamma, xi, alpha, theta) %>%
  pivot_longer(
    cols = c(beta_bar, gamma, xi, alpha, theta),
    names_to = "parameter",
    values_to = "value"
  )

plot_data$parameter <- factor(
  plot_data$parameter,
  levels = c("beta_bar", "gamma", "xi", "alpha", "theta")
)

facet_labs <- as_labeller(c(
  beta_bar = "bar(beta)",
  gamma    = "gamma",
  xi       = "xi",
  alpha    = "alpha",
  theta    = "theta"
), label_parsed)

p_all <- ggplot(plot_data, aes(x = value, y = bond_price, color = bond_price)) +
  geom_point(alpha = 0.60, size = 1.8) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1) +
  facet_wrap(~ parameter, scales = "free_x", ncol = 2, labeller = facet_labs) +
  scale_color_gradient(
    low = "skyblue",
    high = "red",
    name = "Bond price"
  ) +
  labs(
    title = "LHS results: bond price versus SIRD parameters",
    x = "Parameter value",
    y = "Bond price"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 12)
  )

print(p_all)

############################################################
# 8) SAVE RESULTS
############################################################

write.csv(results, "lhs_sird_bond_price_results.csv", row.names = FALSE)
write.csv(range_table, "lhs_parameter_ranges.csv", row.names = FALSE)

ggsave("lhs_sensitivity_barplot.png", p1, width = 8, height = 5, dpi = 300)
ggsave("lhs_all_parameter_scatter.png", p_all, width = 10, height = 8, dpi = 300)