###############################################################################
# Numerical Example 2
# SIRD replacement version for NEW Table 8
#
# New Table 8:
#   Keep all parameters fixed at Scenario I values
#   Only vary beta_bar across:
#     Scenario I   = 0.22
#     Scenario II  = 0.35
#     Scenario III = 0.50
#
# Trigger used here:
#   1) D_total(t) >= 2500
#   2) maI(t)     >= 250
#
# where
#   C_total(t)            = I(t) + R(t) + D(t)
#   daily_new_infections  = max(0, C_total(t) - C_total(t-1))
#   maI(t)                = 7-day moving average of daily new infections
###############################################################################

# =========================
# 0) Packages
# =========================
# install.packages(c("ggpubr","gridExtra","COVID19","demodelr","matrixStats",
#                    "ggplot2","forecast","tseries","RQuantLib","quantmod",
#                    "dplyr","tidyverse","zoo","Sim.DiffProc","rstudioapi"))

library(ggpubr)
library(gridExtra)
library(grid)
library(gtable)
library(COVID19)
library(demodelr)
library(matrixStats)
library(ggplot2)
library(forecast)
library(tseries)
library(RQuantLib)
library(quantmod)
library(dplyr)
library(tidyverse)
library(zoo)
library(Sim.DiffProc)
library(rstudioapi)

# =========================
# 1) Working directory
# =========================
script_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
if (!nzchar(script_path)) stop("Please save this script in RStudio first, then run it again.")
setwd(dirname(script_path))
message("Working directory: ", getwd())

# =========================
# 2) Download and clean COVID data
# =========================
USdata <- covid19(
  country = "US",
  level   = 1,
  start   = "2020-03-01",
  end     = "2022-12-31"
)

USdata[is.na(USdata)] <- 0

covid_data <- USdata[, c("date", "confirmed", "deaths")]

y1 <- USdata$confirmed
y2 <- USdata$deaths

n_obs <- length(y1)
days_obs <- 0:(n_obs - 1)

# =========================
# 3) Daily growth-rate plots
# =========================

# Infection growth rate
daily_infected <- diff(y1, lag = 1)
daily_growth_I <- daily_infected / pmax(head(y1, -1), 1)

dat1 <- data.frame(
  Date        = 1:length(daily_growth_I),
  Growth.rate = daily_growth_I
)

p1 <- ggplot(dat1, aes(x = Date, y = Growth.rate)) +
  geom_line(color = "blue", linewidth = 0.75) +
  labs(
    title   = "",
    x       = "Day",
    y       = "Daily Growth rate: Infections",
    caption = ""
  ) +
  theme_bw(base_family = "Times")

ggsave(filename = "growth_rate_I.png", plot = p1, width = 8, height = 4.5, dpi = 300)

# Death growth rate
daily_death <- diff(y2, lag = 1)
daily_growth_D <- daily_death / pmax(head(y2, -1), 1)

dat2 <- data.frame(
  Date        = 1:length(daily_growth_D),
  Growth.rate = daily_growth_D
)

p2 <- ggplot(dat2, aes(x = Date, y = Growth.rate)) +
  geom_line(color = "blue", linewidth = 0.75) +
  labs(
    title   = "",
    x       = "Day",
    y       = "Daily Growth rate: Death",
    caption = ""
  ) +
  theme_bw(base_family = "Times")

ggsave(filename = "growth_rate_D.png", plot = p2, width = 8, height = 4.5, dpi = 300)

p3 <- ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave(filename = "growth_rate_combine.png", plot = p3, width = 12, height = 4.5, dpi = 300)

# =========================
# 4) NEW Table 8 scenario parameters
#    Keep all Scenario I parameters fixed
#    Only vary beta_bar
# =========================
scenario_params <- data.frame(
  Scenario = c("I", "II", "III"),
  beta_bar = c(0.22, 0.35, 0.50),   # only this changes
  gamma    = c(0.10, 0.10, 0.10),
  xi       = c(0.25, 0.25, 0.25),
  alpha    = c(0.10, 0.10, 0.10),
  beta0    = c(0.15, 0.15, 0.15),
  theta    = c(0.40, 0.40, 0.40),
  stringsAsFactors = FALSE
)

# =========================
# 5) Bond helper functions
# =========================
coupon_data <- readRDS(file = "coupon_example2.rds")

pv <- function(i, cf, t) {
  sum(cf * (i)^t)
}

# full bond value if no trigger
full_bond_value <- function() {
  cf     <- as.numeric(coupon_data$CF)
  y_rate <- as.numeric(coupon_data$Disc.rate)
  t      <- coupon_data$Time / 360
  pv(y_rate, cf, t)
}

# partial value if trigger occurs at trigger_time
partial_bond_value <- function(trigger_time) {
  coupondf <- subset(coupon_data, coupon_data$Time <= trigger_time)
  
  if (nrow(coupondf) == 0) {
    return(0)
  }
  
  cf     <- as.numeric(coupondf$CF)
  y_rate <- as.numeric(coupondf$Disc.rate)
  t      <- coupondf$Time / 360
  
  pv(y_rate, cf, t)
}

# =========================
# 6) SIRD helper functions
# =========================

beta_bar_fn <- function(t, beta, beta0, theta) {
  beta + (beta0 - beta) * exp(-theta * t)
}

sigma_t <- function(t, theta, xi) {
  (xi / sqrt(2 * theta)) * sqrt(1 - exp(-2 * theta * t))
}

roll_mean_7 <- function(x) {
  n <- length(x)
  out <- rep(NA_real_, n)
  if (n < 7) return(out)
  for (i in 7:n) {
    out[i] <- mean(x[(i - 6):i], na.rm = FALSE)
  }
  out
}

simulate_sird_paths <- function(
    T_end,
    M,
    N_pop,
    I0_count,
    R0_count,
    D0_count,
    beta,
    beta0,
    theta,
    xi,
    alpha,
    gamma,
    seed = 123
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
    dBk <- Zk   # dt = 1
    
    S_prev <- S[k, ]
    I_prev <- I[k, ]
    R_prev <- R[k, ]
    D_prev <- D[k, ]
    
    eta_I <- sig * S_prev
    drift_logI <- btk * S_prev - (alpha + gamma) - 0.5 * eta_I^2
    I_new <- I_prev * exp(drift_logI + eta_I * dBk)
    
    eta_S <- -sig * I_prev
    drift_logS <- -btk * I_prev - 0.5 * eta_S^2
    S_new <- S_prev * exp(drift_logS + eta_S * dBk)
    
    R_new <- R_prev + alpha * I_prev
    D_new <- D_prev + gamma * I_prev
    
    S_new <- pmax(S_new, 0)
    I_new <- pmax(I_new, 0)
    R_new <- pmax(R_new, 0)
    D_new <- pmax(D_new, 0)
    
    total_new <- S_new + I_new + R_new + D_new
    bad <- total_new > 1
    if (any(bad)) {
      S_new[bad] <- S_new[bad] / total_new[bad]
      I_new[bad] <- I_new[bad] / total_new[bad]
      R_new[bad] <- R_new[bad] / total_new[bad]
      D_new[bad] <- D_new[bad] / total_new[bad]
    }
    
    S[k + 1, ] <- S_new
    I[k + 1, ] <- I_new
    R[k + 1, ] <- R_new
    D[k + 1, ] <- D_new
  }
  
  list(
    days = days,
    S = S,
    I = I,
    R = R,
    D = D
  )
}

build_trigger_series_sird <- function(sim, N_pop) {
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
  
  maI <- matrix(NA_real_, nrow = n, ncol = M)
  for (j in 1:M) {
    maI[, j] <- roll_mean_7(daily_new_inf[, j])
  }
  
  list(
    I_count       = I_count,
    R_count       = R_count,
    D_total       = D_total,
    C_total       = C_total,
    daily_new_inf = daily_new_inf,
    maI           = maI
  )
}

calc_trigger_tau_sird <- function(days, D_total, maI, theta_D_total = 2500, theta_I_ma7 = 250) {
  M <- ncol(D_total)
  tau <- rep(NA_integer_, M)
  
  for (j in 1:M) {
    cond <- (!is.na(D_total[, j])) &
      (!is.na(maI[, j])) &
      (D_total[, j] >= theta_D_total) &
      (maI[, j]     >= theta_I_ma7)
    
    hit <- which(cond)
    if (length(hit) > 0) {
      tau[j] <- days[hit[1]]
    }
  }
  
  tau
}

# =========================
# 7) Common pricing runner for one SIRD parameter set
# =========================
run_sird_pricing <- function(
    beta,
    gamma,
    xi,
    alpha = 0.10,
    beta0 = 0.15,
    theta = 0.40,
    numsim = 5000,
    n = 1104,
    N_pop = 8.30e9,
    I0_count = 100,
    R0_count = 0,
    D0_count = 0,
    theta_D_total = 2500,
    theta_I_ma7 = 250,
    pandemic_prob = 0.1468,
    seed = 123
) {
  sim <- simulate_sird_paths(
    T_end     = n,
    M         = numsim,
    N_pop     = N_pop,
    I0_count  = I0_count,
    R0_count  = R0_count,
    D0_count  = D0_count,
    beta      = beta,
    beta0     = beta0,
    theta     = theta,
    xi        = xi,
    alpha     = alpha,
    gamma     = gamma,
    seed      = seed
  )
  
  series <- build_trigger_series_sird(sim, N_pop = N_pop)
  
  tau <- calc_trigger_tau_sird(
    days          = sim$days,
    D_total       = series$D_total,
    maI           = series$maI,
    theta_D_total = theta_D_total,
    theta_I_ma7   = theta_I_ma7
  )
  
  triggered <- which(is.finite(tau))
  trigger   <- length(triggered)
  
  BondValue_partial <- c()
  
  if (length(triggered) > 0) {
    for (ii in triggered) {
      trigger_time <- tau[ii]
      presentV <- partial_bond_value(trigger_time)
      BondValue_partial <- c(BondValue_partial, presentV)
    }
  }
  
  price1 <- mean(BondValue_partial, na.rm = TRUE)
  price2 <- full_bond_value()
  
  if (!is.finite(price1)) {
    price1 <- price2
  }
  
  prob1 <- pandemic_prob
  prob2 <- trigger / numsim
  prob3 <- prob1 * prob2
  prob4 <- 1 - prob3
  
  price3 <- (price1 * prob3) + (price2 * prob4)
  
  list(
    Bond_Price = price3,
    P_C        = prob3,
    trigger    = trigger,
    tau        = tau
  )
}

# =========================
# 8) Scenario pricing (NEW Table 8)
# =========================
set.seed(123)

Total_infected <- max(y1, na.rm = TRUE)
Total_death    <- max(y2, na.rm = TRUE)

Bond_Price <- c()
P_C        <- c()

for (j in 1:nrow(scenario_params)) {
  print(paste0("Simulating SIRD scenario ", scenario_params$Scenario[j]))
  
  out <- run_sird_pricing(
    beta           = scenario_params$beta_bar[j],
    gamma          = scenario_params$gamma[j],
    xi             = scenario_params$xi[j],
    alpha          = scenario_params$alpha[j],
    beta0          = scenario_params$beta0[j],
    theta          = scenario_params$theta[j],
    numsim         = 5000,
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
  
  Bond_Price <- c(Bond_Price, out$Bond_Price)
  P_C        <- c(P_C, out$P_C)
}

print(round(Bond_Price, 4))
print(round(P_C, 4))
print(mean(Bond_Price, na.rm = TRUE))

Bond_Price_3 <- Bond_Price
P_C_3        <- P_C

beta_bar_3vec <- scenario_params$beta_bar
gamma_3vec    <- scenario_params$gamma
xi_3vec       <- scenario_params$xi
alpha_3vec    <- scenario_params$alpha
beta0_3vec    <- scenario_params$beta0
theta_3vec    <- scenario_params$theta

###############################################################################
# Save simulation results
###############################################################################
saveRDS(Bond_Price_3,   file = "Bond_Price_3_SIRD_newbeta.rds")
saveRDS(P_C_3,          file = "P_C_3_SIRD_newbeta.rds")

saveRDS(beta_bar_3vec,  file = "beta_bar_3vec_SIRD_newbeta.rds")
saveRDS(gamma_3vec,     file = "gamma_3vec_SIRD_newbeta.rds")
saveRDS(xi_3vec,        file = "xi_3vec_SIRD_newbeta.rds")
saveRDS(alpha_3vec,     file = "alpha_3vec_SIRD_newbeta.rds")
saveRDS(beta0_3vec,     file = "beta0_3vec_SIRD_newbeta.rds")
saveRDS(theta_3vec,     file = "theta_3vec_SIRD_newbeta.rds")

saveRDS(Total_infected, file = "Total_infected_SIRD_newbeta.rds")
saveRDS(Total_death,    file = "Total_death_SIRD_newbeta.rds")

message("SIRD simulation results saved successfully.")

# =========================
# 9) Build NEW Table 8 data
# =========================
table8_df <- data.frame(
  Scenario = c("# Simulations",
               "N_pop",
               "I(0)",
               "R(0)",
               "D(0)",
               "beta_bar",
               "gamma",
               "xi",
               "alpha",
               "beta0",
               "theta",
               "Death trigger",
               "maI trigger",
               "Bond Price",
               "P(C)",
               "P(H)"),
  I = c("5000",
        format(8.30e9, scientific = TRUE),
        "100",
        "0",
        "0",
        sprintf("%.4f", beta_bar_3vec[1]),
        sprintf("%.4f", gamma_3vec[1]),
        sprintf("%.4f", xi_3vec[1]),
        sprintf("%.4f", alpha_3vec[1]),
        sprintf("%.4f", beta0_3vec[1]),
        sprintf("%.4f", theta_3vec[1]),
        "2500",
        "250",
        sprintf("%.4f", Bond_Price_3[1]),
        sprintf("%.4f", P_C_3[1]),
        "0.1468"),
  II = c("5000",
         format(8.30e9, scientific = TRUE),
         "100",
         "0",
         "0",
         sprintf("%.4f", beta_bar_3vec[2]),
         sprintf("%.4f", gamma_3vec[2]),
         sprintf("%.4f", xi_3vec[2]),
         sprintf("%.4f", alpha_3vec[2]),
         sprintf("%.4f", beta0_3vec[2]),
         sprintf("%.4f", theta_3vec[2]),
         "2500",
         "250",
         sprintf("%.4f", Bond_Price_3[2]),
         sprintf("%.4f", P_C_3[2]),
         "0.1468"),
  III = c("5000",
          format(8.30e9, scientific = TRUE),
          "100",
          "0",
          "0",
          sprintf("%.4f", beta_bar_3vec[3]),
          sprintf("%.4f", gamma_3vec[3]),
          sprintf("%.4f", xi_3vec[3]),
          sprintf("%.4f", alpha_3vec[3]),
          sprintf("%.4f", beta0_3vec[3]),
          sprintf("%.4f", theta_3vec[3]),
          "2500",
          "250",
          sprintf("%.4f", Bond_Price_3[3]),
          sprintf("%.4f", P_C_3[3]),
          "0.1468"),
  stringsAsFactors = FALSE
)

avg_price_table8 <- mean(Bond_Price_3, na.rm = TRUE)

print(table8_df)
print(avg_price_table8)

# =========================
# 10) Helper function for paper-style tables
# =========================
make_paper_table <- function(df, file_name, title_text, footer_text = NULL,
                             width_px = 1800, height_px = 1100, res_dpi = 200,
                             table_x = 0.5, table_y = 0.52,
                             table_width = 0.90, table_height = 0.70,
                             base_fontsize = 16, header_fontsize = 17,
                             title_fontsize = 20, footer_fontsize = 18) {
  
  tg <- tableGrob(
    df,
    rows = NULL,
    theme = ttheme_minimal(
      core = list(
        fg_params = list(fontsize = base_fontsize, fontfamily = "Times"),
        bg_params = list(fill = "white", col = NA)
      ),
      colhead = list(
        fg_params = list(fontsize = header_fontsize, fontface = "bold", fontfamily = "Times"),
        bg_params = list(fill = "white", col = NA)
      )
    )
  )
  
  for (i in seq_len(nrow(df) + 1)) {
    tg <- gtable_add_grob(
      tg,
      grobs = segmentsGrob(
        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
        y0 = unit(0, "npc"), y1 = unit(0, "npc"),
        gp = gpar(lwd = 1.1)
      ),
      t = i, l = 1, r = ncol(tg)
    )
  }
  
  tg <- gtable_add_grob(
    tg,
    grobs = segmentsGrob(
      x0 = unit(0, "npc"), x1 = unit(1, "npc"),
      y0 = unit(1, "npc"), y1 = unit(1, "npc"),
      gp = gpar(lwd = 1.6)
    ),
    t = 1, l = 1, r = ncol(tg)
  )
  
  png(filename = file_name, width = width_px, height = height_px, res = res_dpi)
  grid.newpage()
  
  grid.text(
    title_text,
    x = 0.5, y = 0.95,
    gp = gpar(fontsize = title_fontsize, fontfamily = "Times", fontface = "bold")
  )
  
  grid.draw(
    editGrob(
      tg,
      vp = viewport(
        x = table_x, y = table_y,
        width = table_width, height = table_height
      )
    )
  )
  
  if (!is.null(footer_text)) {
    grid.text(
      footer_text,
      x = 0.5, y = 0.06,
      gp = gpar(fontsize = footer_fontsize, fontfamily = "Times")
    )
  }
  
  dev.off()
}

# =========================
# 11) Save NEW Table 8 like the paper
# =========================
make_paper_table(
  df           = table8_df,
  file_name    = "Table8_SIRD_newbeta_paper_style.png",
  title_text   = "Sensitivity of the fair value of the bond price under SIRD model with Scenario I parameters fixed and varying beta_bar.",
  footer_text  = paste0("Average price: USD ", sprintf("%.3f", avg_price_table8), " million"),
  width_px     = 3000,
  height_px    = 1300,
  table_y      = 0.52,
  table_width  = 0.92,
  table_height = 0.76
)

# =========================
# 12) Save raw table
# =========================
write.csv(table8_df, "Table8_SIRD_newbeta_data.csv", row.names = FALSE)

message("Done. Files created:")
message(" - growth_rate_I.png")
message(" - growth_rate_D.png")
message(" - growth_rate_combine.png")
message(" - Bond_Price_3_SIRD_newbeta.rds")
message(" - P_C_3_SIRD_newbeta.rds")
message(" - beta_bar_3vec_SIRD_newbeta.rds")
message(" - gamma_3vec_SIRD_newbeta.rds")
message(" - xi_3vec_SIRD_newbeta.rds")
message(" - alpha_3vec_SIRD_newbeta.rds")
message(" - beta0_3vec_SIRD_newbeta.rds")
message(" - theta_3vec_SIRD_newbeta.rds")
message(" - Total_infected_SIRD_newbeta.rds")
message(" - Total_death_SIRD_newbeta.rds")
message(" - Table8_SIRD_newbeta_paper_style.png")
message(" - Table8_SIRD_newbeta_data.csv")