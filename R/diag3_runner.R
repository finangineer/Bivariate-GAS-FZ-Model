# ------------------------------------------------------------
# 3-Asset Diagonal BP-GAS-FZ runner
# Dependencies (must be sourced/available):
#   - load_asset_data_3()              [from R/data_load.R]
#   - sampleve()                       [from R/utils.R]
#   - BP_GAS_FZ_Loss_Diag_3Asset()      [from R/model_diag3_loss.R]
#   - BP_GAS_FZ_Diag_3Asset_Forecast()  [from R/model_diag3_forecast.R]
# ------------------------------------------------------------

run_evaluation_tests_diag3 <- function(Yp_oos, v_oos, e_oos, alpha, ticker1, ticker2, ticker3 = "") {

  VaR_Kupiec_Test <- function(Yp, VaR, alpha) {
    breaches <- sum(Yp < VaR, na.rm = TRUE)
    N <- length(Yp)
    p <- breaches / N

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

  list(kupiec = kupiec_result, cc_test = cc_test, berkowitz = berkowitz_result, breaches = breaches)
}


print_summary_diag3 <- function(result, ticker1, ticker2, ticker3 = "") {

  cat("----------------------------------------------------\n")
  cat(ticker1, "vs.", ticker2, if (ticker3 != "") paste("vs.", ticker3) else "", "\n")
  cat("----------------------------------------------------\n")

  cat("Estimated parameters:\n")
  cat(sprintf(" a = %.4f, b = %.4f (c = %.4f)\n",
              result$parameters$a, result$parameters$b, result$parameters$c))

  cat(" B matrix:\n"); print(round(result$parameters$B, 4))
  cat(" A matrix:\n"); print(round(result$parameters$A, 4))
  cat("Final FZ loss:", result$fz_loss, "\n")

  cat("\nKupiec Test Results:\n")
  cat(sprintf(" LR = %.4f\n p-value = %.4f\n breaches = %d\n expected_breaches = %.1f\n",
              result$eval_results$kupiec$LR,
              result$eval_results$kupiec$p_value,
              result$eval_results$kupiec$breaches,
              result$eval_results$kupiec$expected_breaches))

  cat("\nConditional Coverage Test:\n")
  if (!is.null(result$eval_results$cc_test) && !inherits(result$eval_results$cc_test, "try-error")) {
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
  if (ticker3 != "") cat(ticker3, ":", sd(result$Y3_train), "\n")

  cat("Correlation:\n")
  cat(" ", ticker1, "-", ticker2, ":", cor(result$Y1_train, result$Y2_train), "\n")
  if (ticker3 != "") {
    cat(" ", ticker1, "-", ticker3, ":", cor(result$Y1_train, result$Y3_train), "\n")
    cat(" ", ticker2, "-", ticker3, ":", cor(result$Y2_train, result$Y3_train), "\n")
  }

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


run_diag3_model <- function(ticker1, ticker2, ticker3,
                            start_date, end_date,
                            w1, w2, w3,
                            alpha,
                            maxeval = 15000) {

  data <- load_asset_data_3(ticker1, ticker2, ticker3, start_date, end_date)
  Y1 <- data$Y1; Y2 <- data$Y2; Y3 <- data$Y3

  T_total <- length(Y1)
  T_train <- floor(0.75 * T_total)

  Y1_train <- Y1[1:T_train]; Y2_train <- Y2[1:T_train]; Y3_train <- Y3[1:T_train]
  Y1_oos <- Y1[(T_train + 1):T_total]
  Y2_oos <- Y2[(T_train + 1):T_total]
  Y3_oos <- Y3[(T_train + 1):T_total]

  cat("Train:", T_train, "OOS:", T_total - T_train, "\n")

  # GARCH initialization for diagonal persistence (optional; you already used tryCatch)
  garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0)),
    distribution.model = "std"
  )

  garch_fit1 <- tryCatch(ugarchfit(spec = garch_spec, data = Y1_train, solver = "hybrid"),
                         error = function(e) { cat("GARCH fit failed for", ticker1, ":", e$message, "\n"); NULL })
  garch_fit2 <- tryCatch(ugarchfit(spec = garch_spec, data = Y2_train, solver = "hybrid"),
                         error = function(e) { cat("GARCH fit failed for", ticker2, ":", e$message, "\n"); NULL })
  garch_fit3 <- tryCatch(ugarchfit(spec = garch_spec, data = Y3_train, solver = "hybrid"),
                         error = function(e) { cat("GARCH fit failed for", ticker3, ":", e$message, "\n"); NULL })

  beta1_init <- if (is.null(garch_fit1)) 0.9 else as.numeric(coef(garch_fit1)["beta1"])
  beta2_init <- if (is.null(garch_fit2)) 0.9 else as.numeric(coef(garch_fit2)["beta1"])
  beta3_init <- if (is.null(garch_fit3)) 0.9 else as.numeric(coef(garch_fit3)["beta1"])

  # Starting parameters (your approach, with minimal guards)
  Yp_train <- w1 * Y1_train + w2 * Y2_train + w3 * Y3_train
  vebar <- sampleve(Yp_train, alpha)

  roll_sd <- zoo::rollapply(Yp_train, width = 20, FUN = sd, fill = NA)
  avg_sqrtQ <- mean(na.omit(roll_sd), na.rm = TRUE)
  if (!is.finite(avg_sqrtQ) || avg_sqrtQ <= 0) avg_sqrtQ <- sd(Yp_train)

  initial_a <- vebar[1, 1] / avg_sqrtQ
  initial_b <- vebar[1, 2] / avg_sqrtQ
  initial_c <- initial_a / initial_b

  var1 <- var(Y1_train); var2 <- var(Y2_train); var3 <- var(Y3_train)

  theta_start <- c(
    log(-initial_b),                 # theta_b
    qnorm(pmin(initial_c, 0.95)),    # theta_c
    qnorm(pmin(beta1_init, 0.95)),   # beta1
    qnorm(pmin(beta2_init, 0.95)),   # beta2
    qnorm(pmin(beta3_init, 0.95)),   # beta3
    qnorm(0.90), qnorm(0.70), qnorm(0.80),   # beta12,beta13,beta23
    log(0.2 * sqrt(var1)),
    log(0.2 * sqrt(var2)),
    log(0.2 * sqrt(var3)),
    log(0.05), log(0.1), log(0.08)           # alpha12,alpha13,alpha23
  )

  lower <- c(-5, -2, rep(1.0364, 6), rep(-7, 3), rep(-8, 3))
  upper <- c(-1.5, 0.524, rep(2, 6), rep(-2, 3), rep(-3, 3))
  theta_start <- pmin(pmax(theta_start, lower), upper)

  eval_f_wrapper <- function(x, Y1, Y2, Y3, w1, w2, w3, alpha, tau) {
    BP_GAS_FZ_Loss_Diag_3Asset(
      x,
      Y1 = Y1, Y2 = Y2, Y3 = Y3,
      w1 = w1, w2 = w2, w3 = w3,
      alpha = alpha,
      ticker1 = ticker1, ticker2 = ticker2, ticker3 = ticker3,
      tau = tau,
      return_full = FALSE
    )
  }

  max_attempts <- 5
  attempt <- 1
  opt <- NULL

  while (attempt <= max_attempts) {

    cat("Optimization attempt:", attempt, "\n")

    nl_result <- tryCatch(
      nloptr(
        x0 = theta_start,
        eval_f = eval_f_wrapper,
        lb = lower,
        ub = upper,
        opts = list(
          algorithm = "NLOPT_GN_CRS2_LM",
          maxeval = maxeval,
          xtol_rel = 1e-7,
          print_level = 0
        ),
        Y1 = Y1_train, Y2 = Y2_train, Y3 = Y3_train,
        w1 = w1, w2 = w2, w3 = w3,
        alpha = alpha,
        tau = -1
      ),
      error = function(e) { cat("Error:", e$message, "\n"); NULL }
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

  result_train <- optim(
    par = opt,
    fn = BP_GAS_FZ_Loss_Diag_3Asset,
    Y1 = Y1_train, Y2 = Y2_train, Y3 = Y3_train,
    w1 = w1, w2 = w2, w3 = w3,
    alpha = alpha,
    ticker1 = ticker1, ticker2 = ticker2, ticker3 = ticker3,
    tau = -1,
    return_full = FALSE,
    method = "BFGS",
    control = list(maxit = 3000, trace = 1, reltol = 1e-9)
  )

  opt <- result_train$par

  b <- -exp(opt[1]); c <- pnorm(opt[2]); a <- c * b
  beta <- pmax(0.7, pnorm(opt[3:8]) * 2 - 1)
  Bmat <- diag(beta)
  alpha_params <- -exp(opt[9:14])
  Amat <- diag(alpha_params)

  cat("Estimated parameters:\n")
  cat(sprintf(" a = %.4f, b = %.4f (c = %.4f)\n", a, b, c))
  cat(" B matrix:\n"); print(round(Bmat, 4))
  cat(" A matrix:\n"); print(round(Amat, 4))
  cat("Final FZ loss:", result_train$value, "\n")

  in_sample <- BP_GAS_FZ_Loss_Diag_3Asset(
    opt, Y1_train, Y2_train, Y3_train, w1, w2, w3, alpha,
    ticker1, ticker2, ticker3, tau = -1, return_full = TRUE
  )

  oos_results <- BP_GAS_FZ_Diag_3Asset_Forecast(
    opt,
    Y1_train, Y2_train, Y3_train,
    Y1_oos, Y2_oos, Y3_oos,
    w1, w2, w3,
    alpha,
    ticker1, ticker2, ticker3
  )

  Yp_oos <- w1 * Y1_oos + w2 * Y2_oos + w3 * Y3_oos

  eval_results <- run_evaluation_tests_diag3(
    Yp_oos, oos_results$v_oos, oos_results$e_oos,
    alpha, ticker1, ticker2, ticker3
  )

  cat("Volatility:", sd(Y1_train), sd(Y2_train), sd(Y3_train), "\n")
  cat("Correlation:\n")
  cat(" ", ticker1, "-", ticker2, ":", cor(Y1_train, Y2_train), "\n")
  cat(" ", ticker1, "-", ticker3, ":", cor(Y1_train, Y3_train), "\n")
  cat(" ", ticker2, "-", ticker3, ":", cor(Y2_train, Y3_train), "\n")

  result <- list(
    parameters = list(a = a, b = b, c = c, B = Bmat, A = Amat),
    fz_loss = result_train$value,
    in_sample = in_sample,
    oos_results = oos_results,
    Yp_oos = Yp_oos,
    Y1_train = Y1_train, Y2_train = Y2_train, Y3_train = Y3_train,
    Y1_oos = Y1_oos, Y2_oos = Y2_oos, Y3_oos = Y3_oos,
    Yp_train = w1 * Y1_train + w2 * Y2_train + w3 * Y3_train,
    w1 = w1, w2 = w2, w3 = w3,
    ticker1 = ticker1, ticker2 = ticker2, ticker3 = ticker3,
    eval_results = eval_results,
    summary = function() print_summary_diag3(result, ticker1, ticker2, ticker3)
  )

  result
}
