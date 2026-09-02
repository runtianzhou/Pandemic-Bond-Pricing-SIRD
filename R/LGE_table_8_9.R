###############################################################################
# Numerical Example 2
# Paper-replication version for Figure 8, Table 8, and Table 9
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

# graphics.off()

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

n    <- length(y1)
days <- 0:(n - 1)
x1   <- days
x2   <- days

# =========================
# 3) Daily growth-rate plots
# =========================

# Infection growth rate
daily_infected <- diff(y1, lag = 1)
daily_growth_I <- daily_infected / head(y1, -1)

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

# IMPORTANT: paper-style implementation
# This is intentionally kept as in the original code path to reproduce the paper
daily_growth_D <- daily_death / head(y1, -1)

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

# Combined plot
p3 <- ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave(filename = "growth_rate_combine.png", plot = p3, width = 12, height = 4.5, dpi = 300)

# =========================
# 4) Compute candidate growth / volatility parameters
# =========================

# Growth parameters
gI_1 <- round(max(daily_growth_I, na.rm = TRUE), 4)
gD_1 <- round(max(daily_growth_D, na.rm = TRUE), 4)

gI_2 <- round(mean(daily_growth_I[daily_growth_I > 0], na.rm = TRUE), 4)
gD_2 <- round(mean(daily_growth_D[daily_growth_D > 0], na.rm = TRUE), 4)

gI_3 <- round(mean(daily_growth_I[1:50], na.rm = TRUE), 4)
gD_3 <- round(mean(daily_growth_D[1:50], na.rm = TRUE), 4)

# Volatility parameters
sI_1 <- round(sd(daily_growth_I, na.rm = TRUE), 4)
sD_1 <- round(sd(daily_growth_D, na.rm = TRUE), 4)

sI_2 <- round(sd(daily_growth_I[daily_growth_I > 0], na.rm = TRUE), 4)
sD_2 <- round(sd(daily_growth_D[daily_growth_D > 0], na.rm = TRUE), 4)

sI_3 <- round(sd(daily_growth_I[1:50], na.rm = TRUE), 4)
sD_3 <- round(sd(daily_growth_D[1:50], na.rm = TRUE), 4)

# =========================
# 5) Scenario pricing (Table 8)
# =========================
set.seed(123)

y1 <- USdata$confirmed
y2 <- USdata$deaths

Total_infected <- max(y1, na.rm = TRUE)
Total_death    <- max(y2, na.rm = TRUE)

g_I <- c(gI_1, gI_2, gI_3)
g_D <- c(gD_1, gD_2, gD_3)
s_I <- c(sI_1, sI_2, sI_3)
s_D <- c(sD_1, sD_2, sD_3)

Bond_Price <- c()
P_C        <- c()

# Read coupon data once
coupon_data <- readRDS(file = "coupon_example2.rds")

pv <- function(i, cf, t) {
  sum(cf * (i)^t)
}

for (j in 1:3) {
  print(paste0("Simulating scenario ", j))
  
  y1_max  <- Total_infected
  y2_max  <- Total_death
  y1_min  <- 1
  y2_min  <- 1
  numsim  <- 5000
  n       <- 1104
  theta_I <- 5000
  theta_D <- 2500
  
  s0_I <- s_I[j]
  r0_I <- g_I[j]
  K0_I <- y1_max
  x0_I <- y1_min
  
  s0_D <- s_D[j]
  r0_D <- g_I[j]    # IMPORTANT: paper-style choice, same growth level
  K0_D <- y2_max
  x0_D <- y2_min
  
  dt  <- 1
  m   <- n
  M   <- numsim
  rho <- 0.5
  
  fun_st_logistic <- function(dt, m, rho, M) {
    fx <- expression(r0_I * x * (1 - x / K0_I), r0_D * y * (1 - y / K0_D))
    gx <- expression(s0_I * x * (1 - x / K0_I), s0_D * y * (1 - y / K0_D))
    x0 <- c(x0_I, x0_D)
    t0 <- 0
    Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
    
    set.seed(123)
    res <- snssde2d(
      N         = m,
      M         = M,
      x0        = x0,
      t0        = t0,
      Dt        = dt,
      drift     = fx,
      diffusion = gx,
      corr      = Sigma
    )
    return(res)
  }
  
  res <- fun_st_logistic(dt, m, rho, M)
  
  simulated_I <- matrix(unlist(res$X), ncol = numsim, byrow = FALSE)
  simulated_D <- matrix(unlist(res$Y), ncol = numsim, byrow = FALSE)
  
  Increment_I <- colDiffs(simulated_I, lag = 1)
  Increment_I <- pmax(Increment_I, 0)
  nn <- nrow(Increment_I)
  
  index <- 1:nn
  
  sim_avg7_I <- matrix(ncol = numsim, nrow = length(7:nn))
  sim_GR_I   <- matrix(ncol = numsim, nrow = length(7:nn))
  
  for (i in 1:numsim) {
    obj1 <- zoo(Increment_I[, i], index)
    
    seven_avg_I <- rollmean(obj1, 7, align = "right", fill = 0)
    seven_sd_I  <- rollapply(data = obj1, width = 7, FUN = sd, align = "right", fill = 0)
    
    y11 <- as.data.frame(seven_avg_I[7:nn])[[1]]
    y33 <- as.data.frame(seven_sd_I[7:nn])[[1]]
    
    sim_avg7_I[, i] <- y11
    sim_GR_I[, i]   <- (y11 - (1.533 * y33))
  }
  
  Increment_D <- colDiffs(simulated_D, lag = 1)
  Increment_D <- pmax(Increment_D, 0)
  sim_avg7_D  <- matrix(ncol = numsim, nrow = length(7:nn))
  
  for (i in 1:numsim) {
    obj1 <- zoo(Increment_D[, i], index)
    seven_avg_D <- rollmean(obj1, 7, align = "right", fill = 0)
    y22 <- as.data.frame(seven_avg_D[7:nn])[[1]]
    sim_avg7_D[, i] <- y22
  }
  
  trigger      <- 0
  common_index <- c()
  common_time  <- c()
  
  for (i in 1:numsim) {
    set_I  <- which(sim_avg7_I[, i] > theta_I)
    set_D  <- which(sim_avg7_D[, i] > theta_D)
    set_GR <- which(sim_GR_I[, i] > 0)
    
    common <- Reduce(intersect, list(set_I, set_D, set_GR))
    
    if (length(common) > 0) {
      trigger      <- trigger + 1
      common_index <- append(common_index, i)
      common_time  <- append(common_time, common[1])
    }
  }
  
  BondValue_partial <- c()
  
  for (i in 1:length(common_index)) {
    trigger_time <- common_time[i] + 6
    
    coupondf <- subset(coupon_data, coupon_data$Time <= trigger_time)
    
    cf     <- as.numeric(coupondf$CF)
    y_rate <- as.numeric(coupondf$Disc.rate)
    t      <- coupondf$Time / 360
    
    presentV <- pv(y_rate, cf, t)
    BondValue_partial <- append(BondValue_partial, presentV)
  }
  
  price1 <- mean(BondValue_partial, na.rm = TRUE)
  price2 <- 225
  
  prob1 <- 0.1468
  prob2 <- trigger / numsim
  prob3 <- prob1 * prob2
  prob4 <- 1 - prob3
  
  price3 <- (price1 * prob3) + (price2 * prob4)
  
  Bond_Price <- append(Bond_Price, price3)
  P_C        <- append(P_C, prob3)
}

print(round(Bond_Price, 4))
print(round(P_C, 4))
print(mean(Bond_Price, na.rm = TRUE))

Bond_Price_3 <- Bond_Price
P_C_3        <- P_C
g_I_3vec     <- g_I
g_D_3vec     <- g_D
s_I_3vec     <- s_I
s_D_3vec     <- s_D

# =========================
# 6) Expanded grid pricing (Table 9)
# =========================
set.seed(123)

y1 <- USdata$confirmed
y2 <- USdata$deaths

Total_infected <- max(y1, na.rm = TRUE)
Total_death    <- max(y2, na.rm = TRUE)

g_I <- seq(0.1, 0.4, 0.1)
g_D <- seq(0.1, 0.4, 0.1)
s_I <- seq(0.02, 0.2, 0.02)
s_D <- seq(0.02, 0.2, 0.02)
rho_grid <- seq(0.1, 0.9, 0.1)

d1 <- expand.grid(g_I = g_I, g_D = g_D, s_I = s_I, s_D = s_D, rho = rho_grid)
len_d1 <- nrow(d1)

Bond_Price <- c()
P_C        <- c()

# Read once
coupon_data <- readRDS(file = "coupon_example2.rds")

for (j in 1:len_d1) {
  if (j %% 100 == 0) {
    print(paste0("Simulating scenario ", j, " out of ", len_d1))
  }
  
  y1_max  <- Total_infected
  y2_max  <- Total_death
  y1_min  <- 1
  y2_min  <- 1
  numsim  <- 100
  n       <- 1104
  theta_I <- 5000
  theta_D <- 2500
  
  s0_I <- d1$s_I[j]
  r0_I <- d1$g_I[j]
  K0_I <- y1_max
  x0_I <- y1_min
  
  s0_D <- d1$s_D[j]
  r0_D <- d1$g_D[j]
  K0_D <- y2_max
  x0_D <- y2_min
  
  dt  <- 1
  m   <- n
  M   <- numsim
  rho <- d1$rho[j]
  
  fun_st_logistic <- function(dt, m, rho, M) {
    fx <- expression(r0_I * x * (1 - x / K0_I), r0_D * y * (1 - y / K0_D))
    gx <- expression(s0_I * x * (1 - x / K0_I), s0_D * y * (1 - y / K0_D))
    x0 <- c(x0_I, x0_D)
    t0 <- 0
    Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
    
    set.seed(123)
    res <- snssde2d(
      N         = m,
      M         = M,
      x0        = x0,
      t0        = t0,
      Dt        = dt,
      drift     = fx,
      diffusion = gx,
      corr      = Sigma
    )
    return(res)
  }
  
  res <- fun_st_logistic(dt, m, rho, M)
  
  simulated_I <- matrix(unlist(res$X), ncol = numsim, byrow = FALSE)
  simulated_D <- matrix(unlist(res$Y), ncol = numsim, byrow = FALSE)
  
  Increment_I <- colDiffs(simulated_I, lag = 1)
  Increment_I <- pmax(Increment_I, 0)
  nn <- nrow(Increment_I)
  
  index <- 1:nn
  
  sim_avg7_I <- matrix(ncol = numsim, nrow = length(7:nn))
  sim_GR_I   <- matrix(ncol = numsim, nrow = length(7:nn))
  
  for (i in 1:numsim) {
    obj1 <- zoo(Increment_I[, i], index)
    
    seven_avg_I <- rollmean(obj1, 7, align = "right", fill = 0)
    seven_sd_I  <- rollapply(data = obj1, width = 7, FUN = sd, align = "right", fill = 0)
    
    y11 <- as.data.frame(seven_avg_I[7:nn])[[1]]
    y33 <- as.data.frame(seven_sd_I[7:nn])[[1]]
    
    sim_avg7_I[, i] <- y11
    sim_GR_I[, i]   <- (y11 - (1.533 * y33))
  }
  
  Increment_D <- colDiffs(simulated_D, lag = 1)
  Increment_D <- pmax(Increment_D, 0)
  sim_avg7_D  <- matrix(ncol = numsim, nrow = length(7:nn))
  
  for (i in 1:numsim) {
    obj1 <- zoo(Increment_D[, i], index)
    seven_avg_D <- rollmean(obj1, 7, align = "right", fill = 0)
    y22 <- as.data.frame(seven_avg_D[7:nn])[[1]]
    sim_avg7_D[, i] <- y22
  }
  
  trigger      <- 0
  common_index <- c()
  common_time  <- c()
  
  for (i in 1:numsim) {
    set_I  <- which(sim_avg7_I[, i] > theta_I)
    set_D  <- which(sim_avg7_D[, i] > theta_D)
    set_GR <- which(sim_GR_I[, i] > 0)
    
    common <- Reduce(intersect, list(set_I, set_D, set_GR))
    
    if (length(common) > 0) {
      trigger      <- trigger + 1
      common_index <- append(common_index, i)
      common_time  <- append(common_time, common[1])
    }
  }
  
  BondValue_partial <- c()
  
  for (i in 1:length(common_index)) {
    trigger_time <- common_time[i] + 6
    
    coupondf <- subset(coupon_data, coupon_data$Time <= trigger_time)
    
    cf     <- as.numeric(coupondf$CF)
    y_rate <- as.numeric(coupondf$Disc.rate)
    t      <- coupondf$Time / 360
    
    presentV <- pv(y_rate, cf, t)
    BondValue_partial <- append(BondValue_partial, presentV)
  }
  
  price1 <- mean(BondValue_partial, na.rm = TRUE)
  price2 <- 225
  
  prob1 <- 0.1468
  prob2 <- trigger / numsim
  prob3 <- prob1 * prob2
  prob4 <- 1 - prob3
  
  price3 <- (price1 * prob3) + (price2 * prob4)
  
  Bond_Price <- append(Bond_Price, price3)
  P_C        <- append(P_C, prob3)
}

df_expand <- data.frame(
  Price     = Bond_Price,
  Trig.Prob = P_C
)

print(round(Bond_Price, 4))
print(round(P_C, 4))
print(mean(Bond_Price, na.rm = TRUE))

###############################################################################
# Save simulation results
###############################################################################
saveRDS(Bond_Price_3,   file = "Bond_Price_3.rds")
saveRDS(P_C_3,          file = "P_C_3.rds")
saveRDS(g_I_3vec,       file = "g_I_3vec.rds")
saveRDS(g_D_3vec,       file = "g_D_3vec.rds")
saveRDS(s_I_3vec,       file = "s_I_3vec.rds")
saveRDS(s_D_3vec,       file = "s_D_3vec.rds")

saveRDS(df_expand,      file = "df_expand.rds")
saveRDS(d1,             file = "d1.rds")

saveRDS(Total_infected, file = "Total_infected.rds")
saveRDS(Total_death,    file = "Total_death.rds")

message("Simulation results saved successfully.")

# =========================
# 7) Build Table 8 data
# =========================
table8_df <- data.frame(
  Scenario = c("# Simulations",
               "K_I",
               "K_D",
               "I(0)",
               "D(0)",
               "g_I, g_D",
               "sigma_I",
               "sigma_D",
               "rho",
               "Bond Price",
               "P(C)",
               "P(H)"),
  I = c("5000",
        paste0(round(Total_infected / 1e6, 2), " M"),
        paste0(round(Total_death / 1e6, 1), " M"),
        "1",
        "1",
        sprintf("%.4f", g_I_3vec[1]),
        sprintf("%.4f", s_I_3vec[1]),
        sprintf("%.4f", s_D_3vec[1]),
        "0.5",
        sprintf("%.4f", Bond_Price_3[1]),
        sprintf("%.4f", P_C_3[1]),
        "0.1468"),
  II = c("5000",
         paste0(round(Total_infected / 1e6, 2), " M"),
         paste0(round(Total_death / 1e6, 1), " M"),
         "1",
         "1",
         sprintf("%.4f", g_I_3vec[2]),
         sprintf("%.4f", s_I_3vec[2]),
         sprintf("%.4f", s_D_3vec[2]),
         "0.5",
         sprintf("%.4f", Bond_Price_3[2]),
         sprintf("%.4f", P_C_3[2]),
         "0.1468"),
  III = c("5000",
          paste0(round(Total_infected / 1e6, 2), " M"),
          paste0(round(Total_death / 1e6, 1), " M"),
          "1",
          "1",
          sprintf("%.4f", g_I_3vec[3]),
          sprintf("%.4f", s_I_3vec[3]),
          sprintf("%.4f", s_D_3vec[3]),
          "0.5",
          sprintf("%.4f", Bond_Price_3[3]),
          sprintf("%.4f", P_C_3[3]),
          "0.1468"),
  stringsAsFactors = FALSE
)

avg_price_table8 <- mean(Bond_Price_3, na.rm = TRUE)

# =========================
# 8) Build Table 9 data
# =========================
table9_df <- data.frame(
  Item = c("g_I",
           "g_D",
           "s_I",
           "s_D",
           "rho",
           "# of parameter combinations",
           "# Simulations",
           "Average bond price under this grid",
           "Average probability for trigger activation"),
  Value = c("0.1, 0.2, 0.3, 0.4",
            "0.1, 0.2, 0.3, 0.4",
            "0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20",
            "0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20",
            "0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9",
            format(nrow(d1), big.mark = ","),
            "each combination 100 times",
            paste0("USD ", sprintf("%.4f", mean(df_expand$Price, na.rm = TRUE)), " million"),
            sprintf("%.4f", mean(df_expand$Trig.Prob, na.rm = TRUE))),
  stringsAsFactors = FALSE
)

# =========================
# 9) Helper function for paper-style tables
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
# 10) Save Table 8 like the paper
# =========================
make_paper_table(
  df           = table8_df,
  file_name    = "Table8_paper_style.png",
  title_text   = "Table 8. Numeric example 2: The fair value of the bond price under proposed model.",
  footer_text  = paste0("Average price: USD ", sprintf("%.3f", avg_price_table8), " million"),
  width_px     = 2600,
  height_px    = 1200,
  table_y      = 0.52,
  table_width  = 0.92,
  table_height = 0.72
)

# =========================
# 11) Save Table 9 like the paper
# =========================
make_paper_table(
  df           = table9_df,
  file_name    = "Table9_paper_style.png",
  title_text   = "Table 9. Numerical example 2: average pandemic bond price under various parameter combinations.",
  footer_text  = NULL,
  width_px     = 2600,
  height_px    = 1000,
  table_y      = 0.52,
  table_width  = 0.92,
  table_height = 0.68
)

# =========================
# 12) Save raw tables
# =========================
write.csv(table8_df, "Table8_data.csv", row.names = FALSE)
write.csv(table9_df, "Table9_data.csv", row.names = FALSE)

message("Done. Files created:")
message(" - growth_rate_I.png")
message(" - growth_rate_D.png")
message(" - growth_rate_combine.png")
message(" - Bond_Price_3.rds")
message(" - P_C_3.rds")
message(" - g_I_3vec.rds")
message(" - g_D_3vec.rds")
message(" - s_I_3vec.rds")
message(" - s_D_3vec.rds")
message(" - df_expand.rds")
message(" - d1.rds")
message(" - Total_infected.rds")
message(" - Total_death.rds")
message(" - Table8_paper_style.png")
message(" - Table9_paper_style.png")
message(" - Table8_data.csv")
message(" - Table9_data.csv")