# ------------------------------------------------------------
# Block-B BP-GAS-FZ runner
# - B is block-diagonal: 2x2 block + separate 3rd diagonal element
# - A is full 3x3
# Estimation: nloptr global search + BFGS refinement
# Dependencies (must be sourced/available):
#   - load_asset_data()   [from R/data_load.R]
#   - sampleve()          [from R/utils.R]
#   - BP_GAS_FZ_Loss_BlockB_Final()   [from R/model_blockB_loss.R]
#   - BP_GAS_FZ_BlockB_Forecast()     [from R/model_blockB_forecast.R]
# ------------------------------------------------------------

print_summary_blockB <- function(result, ticker1, ticker2) {

  cat("----------------------------------------------------\n")
  cat(ticker1, "vs.", ticker2, "\n")
  cat("----------------------------------------------------\n")

  cat("Estimated parameters:\n")
  cat(sprintf(" a = %.4f, b = %.4f (c = %.4f)\n",
              result$parameters$a, result$parameters$b, result$parameters$c))

  cat(" B matrix:\n")
  print(round(result$parameters$B, 4))

  cat(" A matrix:\n")
  print(round(result$parameters$A, 4))

  cat("Final FZ loss:", result$fz_loss, "\n")

  cat("\nKupiec Test Results:\n")
  cat(sprintf(" LR = %.4f\n p-value = %.4f\n breaches = %d\n expected_breaches = %.1f\n",
              result$eval_results$kupiec$LR,
              result$eval_results$kupiec$p_value,
              result$eval_results$kupiec$breaches,
              result$eval_results$kupiec$expected_breaches))

  cat("\nConditional Coverage Test:\n")
  if (!is.null(result$eval_results$cc_test) &&
      !inherits(result$eval_results$cc_test, "try-error")) {

    cc <- result$eval_results$cc_test

    cat(sprintf(
      paste0(
        " expected.exceed = %d\n actual.exceed = %d\n",
        " uc.H0 = %s\n uc.LRstat = %.4f\n uc.critical = %.4f\n uc.LRp = %.4f\n uc.Decision = %s\n",
        " cc.H0 = %s\n cc.LRstat = %.4f\n cc.critical = %.4f\n cc.LRp = %.4f\n cc.Decision = %s\n"
      ),
      cc$expected.exceed, cc$actual.exceed,
      cc$uc.H0, cc$uc.LRstat, cc$uc.critical, cc$uc.LRp, cc$uc.Decision,
      cc$cc.H0, cc$cc.LRstat, cc$cc.critical, cc$cc.LRp, cc$cc.Decision
    ))

  } else {
    cat("Conditional coverage test failed.\n")
  }

  cat("\nBerkowitz Test p-value:", result$eval_results$berkowitz$p_value, "\n")
  if (is.na(result$eval_results$berkowitz$p_value)) {
    cat("Berkowitz test could not be computed due to insufficient/degenerate breaches.\n")
  } else if (result$eval_results$berkowitz$p_value > 0.05) {
    cat("Berkowitz test: No evidence of VaR-ES misspecification (p > 0.05).\n")
  } else {
    cat("Berkowitz test: Evidence of VaR-ES misspecification (p <= 0.05).\n")
  }

  cat("Number of OOS breaches:", length(result$eval_results$breaches), "\n")

  if (length(result$eval_results$breaches) > 0) {
    Z1_vals <- result$Yp_oos[result$eval_results$breaches] /
      result$oos_results$e_oos[result$eval_results$breaches]
    Z1_stat <- mean(Z1_vals, na.rm = TRUE)
    cat("Z1 Statistic:", Z1_stat, "Expected ~1 for calibration\n")
  } else {
    cat("No breaches occurred; cannot compute Z1 statistic.\n")
  }

  cat("Volatility (standard deviation):\n")
  cat(ticker1, ":", sd(result$Y1_train), ticker2, ":", sd(result$Y2_train), "\n")
  cat("Correlation:", cor(result$Y1_train, result$Y2_train), "\n")

  losses <- result$oos_results$loss_oos
  mean_loss <- mean(losses, na.rm = TRUE)
  sd_loss <- sd(losses, na.rm = TRUE)
  n_loss <- sum(!is.na(losses))

  cat("\n--- Out-of-Sample Loss ---\n")
  cat(sprintf("Average OOS Loss: %.6f\n", mean_loss))

  if (n_loss > 1 && is.finite(sd_loss) && sd_loss > 0) {
    error <- qt(0.975, df = n_loss - 1) * sd_loss / sqrt(n_loss)
    cat(sprintf("95%% CI for Loss: [%.6f, %.6f]\n", mean_loss - error, mean_loss + error))
  } else {
    cat("95% CI for Loss: not available (too few non-NA observations).\n")
  }
}


run_evaluation_tests <- function(Yp_oos, v_oos, e_oos, alpha, ticker1, ticker2) {

  VaR_Kupiec_Test <- function(Yp, VaR, alpha) {
    breaches <- sum(Yp < VaR, na.rm = TRUE)
    N <- length(Yp)
    p <- breaches / N

    # Guard: p=0 or p=1 makes log((1-p)^(...)) valid, but be cautious numerically
    LR <- -2 * (
      log((1 - alpha)^(N - breaches) * alpha^breaches) -
        log((1 - p)^(N - breaches) * p^breaches)
    )

    p_value <- 1 - pchisq(LR, df = 1)
    list(LR = LR, p_value = p_value, breaches = breaches, expected_breaches = alpha * N)
  }

  kupiec_result <- VaR_Kupiec_Test(Yp = Yp_oos, VaR = v_oos, alpha = alpha)

  cc_test <- suppressWarnings(
    try(VaRTest(alpha = alpha, actual = Yp_oos, VaR = v_oos, conf.level = 0.95), silent = TRUE)
  )

  berkowitz_test <- function(Yp, v, e) {
    breaches <- Yp < v
    xb <- Yp[breaches]
    eb <- e[breaches]

    # Need at least 2 breaches and non-degenerate sd to do anything meaningful
    if (length(xb) < 2) return(list(p_value = NA))

    s <- sd(xb, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) return(list(p_value = NA))

    z <- (xb - eb) / s
    if (any(!is.finite(z))) return(list(p_value = NA))

    lm_test <- lm(z ~ 1)
    p_value <- summary(lm_test)$coefficients[1, 4]
    list(p_value = p_value)
  }

  berkowitz_result <- berkowitz_test(Yp_oos, v_oos, e_oos)
  breaches <- which(Yp_oos < v_oos)

  list(
    kupiec = kupiec_result,
    cc_test = cc_test,
    berkowitz = berkowitz_result,
    breaches = breaches
  )
}


run_blockB_model <- function(ticker1, ticker2,
                             start_date, end_date,
                             w1, w2,
                             alpha,
                             maxeval = 8000) {

  data <- load_asset_data(ticker1, ticker2, start_date, end_date)
  Y1 <- data$Y1
  Y2 <- data$Y2

  T_total <- length(Y1)
  T_train <- floor(0.75 * T_total)

  Y1_train <- Y1[1:T_train]
  Y2_train <- Y2[1:T_train]
  Y1_oos <- Y1[(T_train + 1):T_total]
  Y2_oos <- Y2[(T_train + 1):T_total]

  cat("Train:", T_train, "OOS:", T_total - T_train, "\n")

  # GARCH initialization for b11/b22 (optional but useful)
  garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0)),
    distribution.model = "std"
  )

  garch_fit1 <- ugarchfit(spec = garch_spec, data = Y1_train)
  garch_fit2 <- ugarchfit(spec = garch_spec, data = Y2_train)

  beta1_init <- as.numeric(coef(garch_fit1)["beta1"])
  beta2_init <- as.numeric(coef(garch_fit2)["beta1"])

  # clamp to the model's admissible region for diagonals (>=0.7)
  beta1_init <- max(0.7, min(0.999, beta1_init))
  beta2_init <- max(0.7, min(0.999, beta2_init))

  # Starting parameters
  Yp_train <- w1 * Y1_train + w2 * Y2_train
  vebar <- sampleve(Yp_train, alpha)

  theta_start <- c(
    log(-vebar[1, 2]),                 # theta_b
    qnorm(vebar[1, 1] / vebar[1, 2]),  # theta_c
    qnorm((beta1_init + 1) / 2),       # b11
    0.3,                               # b12 (tanh-mapped)
    0.3,                               # b21 (tanh-mapped)
    qnorm((beta2_init + 1) / 2),       # b22
    qnorm((0.75 + 1) / 2),             # b33
    rep(log(0.02), 9)                  # A entries
  )

  lower <- c(-5, -2, 1.0364, -3, -3, 1.0364, 1.0364, rep(-6, 9))
  upper <- c(-1.5, 0.524, 2, 3, 3, 2, 2, rep(-2, 9))
  theta_start <- pmin(pmax(theta_start, lower), upper)

  eval_f_wrapper <- function(x, Y1, Y2, w1, w2, alpha, tau) {
    BP_GAS_FZ_Loss_BlockB_Final(
      x,
      Y1 = Y1, Y2 = Y2,
      w1 = w1, w2 = w2,
      alpha = alpha,
      ticker1 = ticker1, ticker2 = ticker2,
      tau = tau,
      return_full = FALSE
    )
  }

  max_attempts <- 5
  attempt <- 1
  opt <- NULL

  while (attempt <= max_attempts) {

    cat("Optimization attempt:", attempt, "\n")

    this_maxeval <- if (ticker1 == "^FTSE" && ticker2 == "^FCHI") 10000 else maxeval

    nl_result <- tryCatch(
      nloptr(
        x0 = theta_start,
        eval_f = eval_f_wrapper,
        lb = lower,
        ub = upper,
        opts = list(
          algorithm = "NLOPT_GN_CRS2_LM",
          maxeval = this_maxeval,
          xtol_rel = 1e-7,
          print_level = 0
        ),
        Y1 = Y1_train, Y2 = Y2_train,
        w1 = w1, w2 = w2,
        alpha = alpha,
        tau = -1
      ),
      error = function(e) {
        cat("Error:", e$message, "\n")
        NULL
      }
    )

    if (!is.null(nl_result) && nl_result$status >= 1) {
      opt <- nl_result$solution
      cat("nloptr converged:", nl_result$objective, "\n")
      break
    }

    attempt <- attempt + 1
    if (attempt <= max_attempts) {
      theta_start <- theta_start + runif(length(theta_start), -0.03, 0.03)
      theta_start <- pmin(pmax(theta_start, lower), upper)
    }
  }

  if (is.null(opt)) stop("Optimization failed")

  # Local refinement (BFGS)
  result_train <- optim(
    par = opt,
    fn = BP_GAS_FZ_Loss_BlockB_Final,
    Y1 = Y1_train, Y2 = Y2_train,
    w1 = w1, w2 = w2,
    alpha = alpha,
    ticker1 = ticker1, ticker2 = ticker2,
    tau = -1,
    return_full = FALSE,
    method = "BFGS",
    control = list(maxit = 2000, trace = 1, reltol = 1e-8)
  )

  opt <- result_train$par

  # Extract parameters
  b <- -exp(opt[1]); c <- pnorm(opt[2]); a <- c * b

  b11 <- max(0.7, pnorm(opt[3]) * 2 - 1)
  b12 <- 0.95 * tanh(opt[4])
  b21 <- 0.95 * tanh(opt[5])
  b22 <- max(0.7, pnorm(opt[6]) * 2 - 1)
  b33 <- max(0.7, pnorm(opt[7]) * 2 - 1)

  Bmat <- matrix(
    c(b11, b12, 0,
      b21, b22, 0,
      0,   0,   b33),
    nrow = 3, byrow = TRUE
  )

  A_vec <- -exp(opt[8:16])
  Amat <- matrix(A_vec, nrow = 3, ncol = 3, byrow = TRUE)

  cat("Estimated parameters:\n")
  cat(sprintf(" a = %.4f, b = %.4f (c = %.4f)\n", a, b, c))
  cat(" B matrix:\n"); print(round(Bmat, 4))
  cat(" A matrix:\n"); print(round(Amat, 4))
  cat("Final FZ loss:", result_train$value, "\n")

  # In-sample paths (optional)
  in_sample <- BP_GAS_FZ_Loss_BlockB_Final(
    opt, Y1_train, Y2_train, w1, w2, alpha,
    ticker1, ticker2, tau = -1, return_full = TRUE
  )

  # OOS forecast
  oos_results <- BP_GAS_FZ_BlockB_Forecast(
    opt,
    Y1_train, Y2_train,
    Y1_oos, Y2_oos,
    w1, w2,
    alpha,
    ticker1, ticker2
  )

  Yp_oos <- w1 * Y1_oos + w2 * Y2_oos

  eval_results <- run_evaluation_tests(
    Yp_oos = Yp_oos,
    v_oos = oos_results$v_oos,
    e_oos = oos_results$e_oos,
    alpha = alpha,
    ticker1 = ticker1,
    ticker2 = ticker2
  )

  cat("Volatility:", sd(Y1_train), sd(Y2_train), "\n")
  cat("Correlation:", cor(Y1_train, Y2_train), "\n")

  result <- list(
    parameters = list(a = a, b = b, c = c, B = Bmat, A = Amat),
    fz_loss = result_train$value,
    in_sample = in_sample,
    oos_results = oos_results,
    Yp_oos = Yp_oos,
    Y1_train = Y1_train,
    Y2_train = Y2_train,
    Y1_oos = Y1_oos,
    Y2_oos = Y2_oos,
    Yp_train = w1 * Y1_train + w2 * Y2_train,
    w1 = w1, w2 = w2,
    ticker1 = ticker1, ticker2 = ticker2,
    eval_results = eval_results,
    summary = function() print_summary_blockB(result, ticker1, ticker2)
  )

  result
}
