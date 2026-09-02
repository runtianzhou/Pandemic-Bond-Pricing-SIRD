###############################################################################
# NEW TABLE 9 ONLY
# Average / minimum / maximum pandemic bond price
# under the SIRD parameter grid
###############################################################################

#############################
# 1) Build parameter grid
#############################
set.seed(123)

beta_bar_grid <- sort(unique(scenario_params$beta_bar))   # 0.22, 0.35
gamma_grid    <- sort(unique(scenario_params$gamma))      # 0.01, 0.09, 0.10
xi_grid       <- sort(unique(scenario_params$xi))         # 0.12, 0.25, 1.20
alpha_grid    <- sort(unique(scenario_params$alpha))      # 0.01, 0.10
beta0_grid    <- sort(unique(scenario_params$beta0))      # 0.15
theta_grid    <- sort(unique(scenario_params$theta))      # 0.13, 0.40, 1.30

d1 <- expand.grid(
  beta_bar = beta_bar_grid,
  gamma    = gamma_grid,
  xi       = xi_grid,
  alpha    = alpha_grid,
  beta0    = beta0_grid,
  theta    = theta_grid,
  stringsAsFactors = FALSE
)

len_d1 <- nrow(d1)

#############################
# 2) Run pricing on the full grid
#############################
Bond_Price <- numeric(len_d1)
P_C        <- numeric(len_d1)
trigger_n  <- numeric(len_d1)

for (j in 1:len_d1) {
  if (j %% 25 == 0 || j == len_d1) {
    cat("Running grid case", j, "out of", len_d1, "\n")
  }
  
  out <- run_sird_pricing(
    beta           = d1$beta_bar[j],
    gamma          = d1$gamma[j],
    xi             = d1$xi[j],
    alpha          = d1$alpha[j],
    beta0          = d1$beta0[j],
    theta          = d1$theta[j],
    numsim         = 100,
    n              = 1104,
    N_pop          = 8.30e9,
    I0_count       = 100,
    R0_count       = 0,
    D0_count       = 0,
    theta_D_total  = 2500,
    theta_I_ma7    = 250,
    pandemic_prob  = 0.1468,
    seed           = 123
  )
  
  Bond_Price[j] <- out$Bond_Price
  P_C[j]        <- out$P_C
  trigger_n[j]  <- out$trigger
}

#############################
# 3) Store full results
#############################
df_expand <- data.frame(
  beta_bar  = d1$beta_bar,
  gamma     = d1$gamma,
  xi        = d1$xi,
  alpha     = d1$alpha,
  beta0     = d1$beta0,
  theta     = d1$theta,
  Price     = Bond_Price,
  Trig_Prob = P_C,
  Trigger_n = trigger_n,
  stringsAsFactors = FALSE
)

#############################
# 4) Summary statistics
#############################
avg_price_grid <- mean(df_expand$Price, na.rm = TRUE)
avg_prob_grid  <- mean(df_expand$Trig_Prob, na.rm = TRUE)

min_idx <- which.min(df_expand$Price)
max_idx <- which.max(df_expand$Price)

min_row <- df_expand[min_idx, , drop = FALSE]
max_row <- df_expand[max_idx, , drop = FALSE]

cat("\n====================================\n")
cat("AVERAGE BOND PRICE UNDER GRID\n")
cat("====================================\n")
cat(sprintf("USD %.4f million\n", avg_price_grid))

cat("\n====================================\n")
cat("MINIMUM BOND PRICE PARAMETER SET\n")
cat("====================================\n")
print(min_row)

cat("\n====================================\n")
cat("MAXIMUM BOND PRICE PARAMETER SET\n")
cat("====================================\n")
print(max_row)

#############################
# 5) Build main Table 9
#############################
table9_df <- data.frame(
  Item = c("beta_bar",
           "gamma",
           "xi",
           "alpha",
           "beta0",
           "theta",
           "Death trigger",
           "maI trigger",
           "# of parameter combinations",
           "# Simulations",
           "Average bond price under this grid",
           "Minimum bond price under this grid",
           "Parameters for minimum bond price",
           "Maximum bond price under this grid",
           "Parameters for maximum bond price",
           "Average probability for trigger activation"),
  Value = c(
    paste(beta_bar_grid, collapse = ", "),
    paste(sprintf("%.2f", gamma_grid), collapse = ", "),
    paste(sprintf("%.2f", xi_grid), collapse = ", "),
    paste(sprintf("%.2f", alpha_grid), collapse = ", "),
    paste(sprintf("%.2f", beta0_grid), collapse = ", "),
    paste(sprintf("%.2f", theta_grid), collapse = ", "),
    "2500",
    "250",
    format(nrow(d1), big.mark = ","),
    "each combination 100 times",
    paste0("USD ", sprintf("%.4f", avg_price_grid), " million"),
    paste0("USD ", sprintf("%.4f", min_row$Price), " million"),
    paste0(
      "beta_bar=", sprintf("%.2f", min_row$beta_bar),
      ", gamma=", sprintf("%.2f", min_row$gamma),
      ", xi=", sprintf("%.2f", min_row$xi),
      ", alpha=", sprintf("%.2f", min_row$alpha),
      ", beta0=", sprintf("%.2f", min_row$beta0),
      ", theta=", sprintf("%.2f", min_row$theta)
    ),
    paste0("USD ", sprintf("%.4f", max_row$Price), " million"),
    paste0(
      "beta_bar=", sprintf("%.2f", max_row$beta_bar),
      ", gamma=", sprintf("%.2f", max_row$gamma),
      ", xi=", sprintf("%.2f", max_row$xi),
      ", alpha=", sprintf("%.2f", max_row$alpha),
      ", beta0=", sprintf("%.2f", max_row$beta0),
      ", theta=", sprintf("%.2f", max_row$theta)
    ),
    sprintf("%.4f", avg_prob_grid)
  ),
  stringsAsFactors = FALSE
)

#############################
# 6) Build compact min-max table
#############################
table9_minmax_df <- data.frame(
  Case = c("Minimum bond price", "Maximum bond price"),
  BondPrice = c(
    sprintf("%.4f", min_row$Price),
    sprintf("%.4f", max_row$Price)
  ),
  beta_bar = c(
    sprintf("%.2f", min_row$beta_bar),
    sprintf("%.2f", max_row$beta_bar)
  ),
  gamma = c(
    sprintf("%.2f", min_row$gamma),
    sprintf("%.2f", max_row$gamma)
  ),
  xi = c(
    sprintf("%.2f", min_row$xi),
    sprintf("%.2f", max_row$xi)
  ),
  alpha = c(
    sprintf("%.2f", min_row$alpha),
    sprintf("%.2f", max_row$alpha)
  ),
  beta0 = c(
    sprintf("%.2f", min_row$beta0),
    sprintf("%.2f", max_row$beta0)
  ),
  theta = c(
    sprintf("%.2f", min_row$theta),
    sprintf("%.2f", max_row$theta)
  ),
  TrigProb = c(
    sprintf("%.4f", min_row$Trig_Prob),
    sprintf("%.4f", max_row$Trig_Prob)
  ),
  stringsAsFactors = FALSE
)

#############################
# 7) Save raw results
#############################
saveRDS(df_expand, file = "df_expand_SIRD.rds")
saveRDS(d1,        file = "d1_SIRD.rds")

write.csv(df_expand,         "df_expand_SIRD_full.csv", row.names = FALSE)
write.csv(table9_df,         "Table9_SIRD_data.csv", row.names = FALSE)
write.csv(table9_minmax_df,  "Table9_SIRD_minmax_data.csv", row.names = FALSE)

#############################
# 8) Save paper-style Table 9
#############################
make_paper_table(
  df           = table9_df,
  file_name    = "Table9_SIRD_paper_style.png",
  title_text   = "Table. Pandemic bond price summary under various SIRD parameter combinations.",
  footer_text  = NULL,
  width_px     = 2600,
  height_px    = 1400,
  table_y      = 0.52,
  table_width  = 0.92,
  table_height = 0.78
)

#############################
# 9) Save compact min-max table
#############################
make_paper_table(
  df           = table9_minmax_df,
  file_name    = "Table9_SIRD_minmax_paper_style.png",
  title_text   = "Table. Minimum and maximum pandemic bond prices under the SIRD parameter grid.",
  footer_text  = NULL,
  width_px     = 2600,
  height_px    = 900,
  table_y      = 0.52,
  table_width  = 0.95,
  table_height = 0.55
)

#############################
# 10) Final message
#############################
message("Done. Files created:")
message(" - df_expand_SIRD.rds")
message(" - d1_SIRD.rds")
message(" - df_expand_SIRD_full.csv")
message(" - Table9_SIRD_data.csv")
message(" - Table9_SIRD_minmax_data.csv")
message(" - Table9_SIRD_paper_style.png")
message(" - Table9_SIRD_minmax_paper_style.png")