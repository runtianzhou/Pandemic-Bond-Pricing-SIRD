############################################################
# FULL CORRECT R CODING (GR-only trigger, World Bank 84-day GR(t))
# + Mark LAST time GR(t) > 0 with a vertical line + day label
# + DO NOT put graphs together: one plot per path (sequential)
############################################################

#############################
# 0) ONE-PLACE CONFIG
#############################
CFG <- list(
  N     = 1104,
  dt    = 1,
  M     = 1,      # set to 50 if you want; plots will be sequential
  seed  = 123,
  
  # If you already have I_mat (N x M) from your Algorithm 1 simulation,
  # set simulate=FALSE and define I_mat yourself.
  simulate = TRUE,
  
  # Demo simulator params (only if simulate=TRUE)
  I0     = 1,
  KI     = 833e6,
  gI     = 0.320,
  sigmaI = 0.10,
  
  # Plot window
  x_min = 50,
  x_max = 150,
  
  # y-limits for GR plot (NULL = auto per path)
  y_GR = NULL
)

#############################
# 1) World Bank GR(t) (TCA-based, 5×14-day schedule)
#############################
calc_worldbank_GR_TCA <- function(TCA, undefined_value = NA_real_) {
  TCA <- as.numeric(TCA)
  n <- length(TCA)
  
  GR  <- rep(undefined_value, n)
  mu  <- rep(NA_real_, n)
  st  <- rep(NA_real_, n)
  se  <- rep(NA_real_, n)
  
  NCRC <- matrix(NA_real_, nrow = n, ncol = 5)
  colnames(NCRC) <- paste0("NCRC", 1:5)
  
  safe_log_ratio <- function(a, b) {
    if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(NA_real_)
    log(a / b)
  }
  
  for (t in seq_len(n)) {
    if (t <= 84) next  # GR undefined before 84-day lookback exists
    
    d0 <- TCA[t]      - TCA[t - 14]
    d1 <- TCA[t - 14] - TCA[t - 28]
    d2 <- TCA[t - 28] - TCA[t - 42]
    d3 <- TCA[t - 42] - TCA[t - 56]
    d4 <- TCA[t - 56] - TCA[t - 70]
    d5 <- TCA[t - 70] - TCA[t - 84]
    
    ncrc1 <- safe_log_ratio(d0, d1)
    ncrc2 <- safe_log_ratio(d1, d2)
    ncrc3 <- safe_log_ratio(d2, d3)
    ncrc4 <- safe_log_ratio(d3, d4)
    ncrc5 <- safe_log_ratio(d4, d5)
    
    NCRC[t, ] <- c(ncrc1, ncrc2, ncrc3, ncrc4, ncrc5)
    if (anyNA(NCRC[t, ])) next
    
    mu_t <- mean(NCRC[t, ])
    st_t <- sqrt(sum((NCRC[t, ] - mu_t)^2) / 4)  # sample SD, df=4
    se_t <- st_t / sqrt(5)
    
    mu[t] <- mu_t
    st[t] <- st_t
    se[t] <- se_t
    GR[t] <- mu_t - 1.533 * se_t
  }
  
  list(GR = GR, mu = mu, st = st, se = se, NCRC = NCRC)
}

#############################
# 2) Trigger time using ONLY GR(t) > 0 (t>84)
#############################
trigger_time_GR_only <- function(GR_vec) {
  idx <- which(seq_along(GR_vec) > 84 & is.finite(GR_vec) & GR_vec > 0)
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

#############################
# 3) Find + draw LAST time GR(t) > 0 (vertical line + label)
#############################
mark_last_GR_positive <- function(x, GR,
                                  threshold = 0,
                                  line_col = "darkgreen",
                                  line_lty = 3,
                                  line_lwd = 2,
                                  label_cex = 0.85,
                                  label_pos = 1,
                                  label_y = NULL) {
  stopifnot(length(x) == length(GR))
  
  ok <- is.finite(x) & is.finite(GR)
  x_ok <- x[ok]
  g_ok <- GR[ok]
  if (length(x_ok) < 1) return(invisible(NA_real_))
  
  idx_pos <- which(g_ok > threshold)
  if (length(idx_pos) == 0) return(invisible(NA_real_))
  
  i_last <- idx_pos[length(idx_pos)]
  x_last <- x_ok[i_last]
  
  abline(v = x_last, col = line_col, lty = line_lty, lwd = line_lwd)
  
  if (is.null(label_y)) {
    usr <- par("usr")
    label_y <- usr[4] - 0.05 * (usr[4] - usr[3])
  }
  
  text(x = x_last, y = label_y,
       labels = paste0("last GR>0: day ", round(x_last, 2)),
       col = line_col, cex = label_cex, pos = label_pos)
  
  invisible(x_last)
}

#############################
# 4) Plot ONE path only (no combined graphs)
#############################
plot_GR_single <- function(days, GR_vec, M, tau_first = NA_integer_,
                           xlim = NULL, ylim = NULL,
                           main = paste0("World Bank GR(t), M=", M, " with GR>0 intervals")) {
  ok <- is.finite(GR_vec)
  if (!any(ok)) stop("GR has no finite values to plot (need N>=85 and positive 14-day increments).")
  
  if (is.null(xlim)) xlim <- range(days)
  if (is.null(ylim)) ylim <- range(GR_vec[ok])
  
  plot(days, GR_vec, type = "l", lwd = 2,
       xlab = "Day", ylab = "GR(t)",
       main = main, xlim = xlim, ylim = ylim)
  
  abline(h = 0, lty = 2, col = "red")   # GR threshold
  abline(v = 85, lty = 3, col = "blue") # first possible day where t>84
  
  # First trigger time (optional)
  if (is.finite(tau_first)) {
    abline(v = tau_first, lty = 3, col = "purple")
    text(tau_first, ylim[2], labels = paste0("tau(first)>0=", tau_first),
         pos = 4, cex = 0.8, col = "purple")
  }
  
  # LAST GR>0 line + label
  mark_last_GR_positive(days, GR_vec, threshold = 0, line_col = "darkgreen")
}

#############################
# 5) OPTIONAL demo simulator for cumulative I(t) (if needed)
#############################
simulate_logistic_I_EM <- function(CFG) {
  set.seed(CFG$seed)
  N <- CFG$N
  M <- CFG$M
  dt <- CFG$dt
  
  I <- matrix(NA_real_, nrow = N, ncol = M)
  I[1, ] <- CFG$I0
  
  for (m in seq_len(M)) {
    dW <- sqrt(dt) * rnorm(N - 1)
    for (k in 2:N) {
      Ik <- I[k - 1, m]
      a  <- CFG$gI * Ik * (1 - Ik / CFG$KI)
      b  <- CFG$sigmaI * Ik * (1 - Ik / CFG$KI)
      I_next <- Ik + a * dt + b * dW[k - 1]
      I[k, m] <- max(I_next, 0)
    }
  }
  
  list(t = seq_len(N), I = I)
}

#############################
# 6) RUN (GR-only setting)
#############################
graphics.off()  # clear any broken graphics state

# --- Input cumulative series matrix (TCA) ---
if (CFG$simulate) {
  sim <- simulate_logistic_I_EM(CFG)
  days  <- sim$t
  I_mat <- sim$I               # use as TCA for GR(t)
} else {
  stop("Set CFG$simulate=FALSE and provide I_mat (N x M) and days.")
}

# --- Compute GR and trigger times per path ---
N <- nrow(I_mat)
M <- ncol(I_mat)
GR_mat <- matrix(NA_real_, N, M)
tau_first_pos <- rep(NA_integer_, M)
tau_last_pos  <- rep(NA_real_, M)

for (m in seq_len(M)) {
  wb <- calc_worldbank_GR_TCA(TCA = I_mat[, m], undefined_value = NA_real_)
  GR_mat[, m] <- wb$GR
  tau_first_pos[m] <- trigger_time_GR_only(GR_mat[, m])
  
  # compute last GR>0 day (no plotting here)
  ok <- is.finite(GR_mat[, m])
  idx <- which(ok & GR_mat[, m] > 0 & days > 84)
  tau_last_pos[m] <- if (length(idx) == 0) NA_real_ else days[idx[length(idx)]]
}

print(data.frame(path = 1:M, tau_first_pos = tau_first_pos, tau_last_pos = tau_last_pos))

# --- Plot ONE graph at a time (no combined plots) ---
xlim_use <- c(CFG$x_min, CFG$x_max)

for (m in seq_len(M)) {
  plot_GR_single(days, GR_mat[, m], tau_first = tau_first_pos[m],
                 xlim = xlim_use, ylim = CFG$y_GR,
                 main = paste0("World Bank GR(t), M=", M, " with GR>0 intervals"))
  if (m < M) {
    cat("\nShowing Path", m, "of", M, " — press Enter for next plot...\n")
    readline()
  }
}