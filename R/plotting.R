# ------------------------------------------------------------
# Plotting utilities for BP-GAS-FZ results
# Expected fields in `result` (depending on plot):
# - result$in_sample$f_ts, result$in_sample$rho_ts (optional)
# - result$in_sample$v_ts, result$in_sample$e_t
# - result$Yp_train, result$Yp_oos
# - result$oos_results$v_oos, result$oos_results$e_oos, result$oos_results$rho_oos (optional)
# - result$eval_results$breaches (optional)
# - result$Y1_train, result$Y2_train, result$Y1_oos, result$Y2_oos
# ------------------------------------------------------------

# Helper: safe range padding
.pad_range <- function(x, pad = 0.05) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(r + c(-1, 1))
  r + c(-pad, pad) * diff(r)
}

# Plot 1: Latent factors (f1,f2,f12) + dynamic correlation rho
plot_latent_factors <- function(result) {

  if (is.null(result$in_sample) || is.null(result$in_sample$f_ts)) {
    stop("result$in_sample$f_ts not found. Run the model with return_full=TRUE in the loss.")
  }

  f_ts <- result$in_sample$f_ts

  # Prefer stored rho_ts; otherwise compute from f12
  rho_ts <- result$in_sample$rho_ts
  if (is.null(rho_ts) && ncol(f_ts) >= 3) {
    rho_ts <- tanh(f_ts[, 3])
  }

  ticker1 <- result$ticker1
  ticker2 <- result$ticker2

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  plot(f_ts[, 1], type = "l",
       main = paste("Latent Factor f1 (", ticker1, "vs.", ticker2, ")"),
       xlab = "Time Index", ylab = "f1")

  plot(f_ts[, 2], type = "l",
       main = paste("Latent Factor f2 (", ticker1, "vs.", ticker2, ")"),
       xlab = "Time Index", ylab = "f2")

  plot(f_ts[, 3], type = "l",
       main = paste("Correlation Driver f12 (", ticker1, "vs.", ticker2, ")"),
       xlab = "Time Index", ylab = "f12")

  if (!is.null(rho_ts)) {
    plot(rho_ts, type = "l",
         main = paste("Dynamic Rho = tanh(f12) (", ticker1, "vs.", ticker2, ")"),
         xlab = "Time Index", ylab = "Rho")
  } else {
    plot.new()
    title("Dynamic Rho unavailable")
  }
}

# Plot 2: OOS VaR/ES vs portfolio returns + breach markers
plot_oos_var_es <- function(result) {

  Yp_oos <- result$Yp_oos
  v_oos <- result$oos_results$v_oos
  e_oos <- result$oos_results$e_oos

  ticker1 <- result$ticker1
  ticker2 <- result$ticker2

  breach_indices <- which(Yp_oos < v_oos)

  ylim <- .pad_range(c(Yp_oos, v_oos, e_oos), pad = 0.05)

  plot(Yp_oos, type = "l",
       main = paste("OOS Returns vs VaR/ES (", ticker1, "vs.", ticker2, ")"),
       xlab = "Time Index", ylab = "Value", ylim = ylim)

  lines(v_oos, lwd = 1.5)
  lines(e_oos, lwd = 1.5)

  if (length(breach_indices) > 0) {
    points(breach_indices, Yp_oos[breach_indices], pch = 4)
  }

  legend("topright",
         legend = c("Portfolio Returns", "VaR", "ES", "Breaches"),
         lty = c(1, 1, 1, NA),
         pch = c(NA, NA, NA, 4),
         bty = "n")
}

# Plot 3: In-sample VaR/ES vs portfolio returns
plot_is_var_es <- function(result) {

  if (is.null(result$in_sample) || is.null(result$in_sample$v_ts) || is.null(result$in_sample$e_t)) {
    stop("Need result$in_sample$v_ts and result$in_sample$e_t. Run loss with return_full=TRUE.")
  }

  Yp_train <- result$Yp_train
  v_ts <- result$in_sample$v_ts
  e_ts <- result$in_sample$e_t

  ticker1 <- result$ticker1
  ticker2 <- result$ticker2

  ylim <- .pad_range(c(Yp_train, v_ts, e_ts), pad = 0.05)

  plot(Yp_train, type = "l",
       main = paste("In-sample Returns vs VaR/ES (", ticker1, "vs.", ticker2, ")"),
       xlab = "Time Index", ylab = "Value", ylim = ylim)

  lines(v_ts, lwd = 1.5)
  lines(e_ts, lwd = 1.5)

  legend("topright",
         legend = c("Portfolio Returns", "VaR", "ES"),
         lty = 1, bty = "n")
}

# Plot 4: OOS breach timeline
plot_oos_breaches <- function(result) {

  Yp_oos <- result$Yp_oos
  ticker1 <- result$ticker1
  ticker2 <- result$ticker2

  # Prefer stored breach indices if available
  breach_indices <- NULL
  if (!is.null(result$eval_results) && !is.null(result$eval_results$breaches)) {
    breach_indices <- result$eval_results$breaches
  } else if (!is.null(result$oos_results$v_oos)) {
    breach_indices <- which(Yp_oos < result$oos_results$v_oos)
  } else {
    breach_indices <- integer(0)
  }

  ylim <- .pad_range(Yp_oos, pad = 0.05)

  plot(Yp_oos, type = "l",
       main = paste("OOS Breach Timeline (", ticker1, "vs.", ticker2, ")"),
       xlab = "Time Index", ylab = "Portfolio Returns", ylim = ylim)

  if (length(breach_indices) > 0) {
    points(breach_indices, Yp_oos[breach_indices], pch = 4)
  }

  legend("topright",
         legend = c("Portfolio Returns", "Breaches"),
         lty = c(1, NA), pch = c(NA, 4),
         bty = "n")
}

# Plot 5: Dynamic rho vs rolling correlation (train + OOS)
plot_dynamic_rho <- function(result, window = 30) {

  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("Package 'zoo' is required for rolling correlation. Please install it.")
  }

  Y1_train <- result$Y1_train
  Y2_train <- result$Y2_train

  Y1_oos <- result$Y1_oos
  Y2_oos <- result$Y2_oos

  # in-sample rho
  rho_ts <- NULL
  if (!is.null(result$in_sample$rho_ts)) {
    rho_ts <- result$in_sample$rho_ts
  } else if (!is.null(result$in_sample$f_ts) && ncol(result$in_sample$f_ts) >= 3) {
    rho_ts <- tanh(result$in_sample$f_ts[, 3])
  } else {
    stop("Need in-sample rho_ts or f_ts to compute rho.")
  }

  # oos rho (only Block-B forecast returns this; others may not)
  rho_oos <- result$oos_results$rho_oos
  if (is.null(rho_oos)) {
    rho_oos <- rep(NA_real_, length(Y1_oos))
  }

  rolling_cor <- zoo::rollapplyr(
    cbind(Y1_train, Y2_train),
    width = window,
    FUN = function(x) cor(x[, 1], x[, 2]),
    by.column = FALSE,
    fill = NA_real_
  )

  rolling_cor_oos <- zoo::rollapplyr(
    cbind(Y1_oos, Y2_oos),
    width = window,
    FUN = function(x) cor(x[, 1], x[, 2]),
    by.column = FALSE,
    fill = NA_real_
  )

  rho_combined <- c(rho_ts, rho_oos)
  rolling_cor_combined <- c(as.numeric(rolling_cor), as.numeric(rolling_cor_oos))

  ticker1 <- result$ticker1
  ticker2 <- result$ticker2

  ylim <- .pad_range(c(rho_combined, rolling_cor_combined), pad = 0.1)

  plot(rho_combined, type = "l",
       main = paste("Dynamic Rho vs Rolling Correlation (", ticker1, "vs.", ticker2, ")"),
       xlab = "Time Index (In-sample + OOS)", ylab = "Correlation",
       ylim = ylim)

  lines(rolling_cor_combined, lwd = 1.5)

  abline(v = length(rho_ts), lty = 2)

  legend("topright",
         legend = c("Dynamic Rho", paste0(window, "-day Rolling Correlation"), "OOS Start"),
         lty = c(1, 1, 2),
         bty = "n")
}
