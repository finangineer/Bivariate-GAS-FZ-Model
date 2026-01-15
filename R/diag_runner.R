# ------------------------------------------------------------
# Diagonal BP-GAS-FZ: helper routines (printing, evaluation, main runner)
# Notes:
# - Designed for interactive use (readers run blocks as needed).
# - No automatic pipeline assumed.
# ------------------------------------------------------------

print_summary <- function(result, ticker1, ticker2) {

  cat("----------------------------------------------------\n")
  cat(ticker1, "vs.", ticker2, "\n")
  cat("----------------------------------------------------\n")

  cat("Estimated parameters:\n")
  cat(sprintf(" a = %.4f, b = %.4f (c = %.4f)\n",
              result$parameters$a, result$parameters$b, result$parameters$c))

  cat(sprintf(" β = diag(%.4f, %.4f, %.4f)\n",
              result$parameters$beta[1, 1], result$parameters$beta[2, 2], result$parameters$beta[3, 3]))

  cat(sprintf(" α = diag(%.4f, %.4f, %.4f)\n",
              result$parameters$alpha[1, 1], result$parameters$alpha[2, 2], result$parameters$alpha[3, 3]))

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
    cat("Berkowitz test could not be computed due to insufficient breaches.\n")
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

  # Out-of-sample loss CI (simple t-interval)
  losses <- result$oos_results$loss_oos
  mean_loss <- mean(losses, na.rm = TRUE)
  sd_loss <- sd(losses, na.rm = TRUE)
  n_loss <- sum(!is.na(losses))
  error <- qt(0.975, df = n_loss - 1) * sd_loss / sqrt(n_loss)

  cat("\n--- Out-of-Sample Loss ---\n")
  cat(sprintf("Average OOS Loss: %.6f\n", mean_loss))
  cat(sprintf("95%% CI for Loss: [%.6f, %.6f]\n", mean_loss - error, mean_loss + error))
}


run_evaluation_tests <- function(Yp_oos, v_oos, e_oos, alpha, ticker1, ticker2) {

  VaR_Kupiec_Test <- function(Yp, VaR, alpha) {
    breaches <- sum(Yp < VaR)
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
    z <- (Yp[breaches] - e[breaches]) / sd(Yp[breaches], na.rm = TRUE)

    if (length(z) > 0) {
      lm_test <- lm(z ~ 1)
      p_value <- summary(lm_test)$coefficients[1, 4]
      return(list(p_value = p_value))
    }

    list(p_value = NA)
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


run_diagonal_model <- function(ticker1, ticker2,
                               start_date, end_date,
                               w1, w2,
                               alpha,
                               maxeval = 5000) {

  data <- load_asset_data(ticker1, ticker2, start_date, end_date)
  Y1 <- data$Y1
  Y2 <- data$Y2
  dates <- data$dates

  T_total <- length(Y1)
  T_train <- floor(0.75 * T_total)

  Y1_train <- Y1[1:T_train]
  Y2_train <- Y2[1:T_train]
  Y1_oos <- Y1[(T_train + 1):T_total]
  Y2_oos <- Y2[(T_train + 1):T_total]

  cat("Train:", T_train, "OOS:", T_total - T_train, "\n")

  var1 <- var(Y1_train)
  var2 <- var(Y2_train)

  log_var1 <- log(pmax(rollapply(Y1^2, width = 20, FUN = var, fill = NA), 1e-6))
  log_var2 <- log(pmax(rollapply(Y2^2, width = 20, FUN = var, fill = NA), 1e-6))

  ar1_beta1 <- tryCatch(coef(arima(na.omit(log_var1), order = c(1, 0, 0)))[1], error = function(e) 0.95)
  ar1_beta2 <- tryCatch(coef(arima(na.omit(log_var2), order = c(1, 0, 0)))[1], error = function(e) 0.85)

  Yp_train <- w1 * Y1_train + w2 * Y2_train
  vebar <- sampleve(Yp_train, alpha)

  theta_start <- c(
    log(-vebar[1, 2]),
    qnorm(vebar[1, 1] / vebar[1, 2]),
    qnorm((0.95 + 1) / 2),
    qnorm((0.85 + 1) / 2),
    qnorm((0.7 + 1) / 2),
    log(0.05 * sqrt(var1)),
    log(0.05 * sqrt(var2)),
    log(0.02)
  )

  lower <- c(-5, -5, rep(1.0364, 3), rep(-5, 3))
  upper <- c(-1.5, 0.524, rep(2, 3), rep(-3, 3))
  theta_start <- pmin(pmax(theta_start, lower), upper)

  eval_f_wrapper <- function(x, Y1, Y2, w1, w2, alpha, tau) {
    BP_GAS_FZ_Loss_Diag_Final(
      x,
      Y1 = Y1, Y2 = Y2,
      w1 = w1, w2 = w2,
      alpha = alpha,
      ticker1 = ticker1, ticker2 = ticker2,
      tau = tau,
      return_full = FALSE
    )
  }

  max_attempts <- 3
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
        Y1 = Y1_train,
        Y2 = Y2_train,
        w1 = w1,
        w2 = w2,
        alpha = alpha,
        tau = -1
      ),
      error = function(e) {
        cat("Error in nloptr attempt", attempt, ":", e$message, "\n")
        NULL
      }
    )

    if (!is.null(nl_result) && nl_result$status >= 1) {
      opt <- nl_result$solution
      cat("nloptr converged with objective:", nl_result$objective, "\n")
      break
    }

    attempt <- attempt + 1
    if (attempt <= max_attempts) {
      cat("Retrying with perturbed starting point...\n")
      theta_start <- theta_start + runif(length(theta_start), -0.03, 0.03)
      theta_start <- pmin(pmax(theta_start, lower), upper)
    }
  }

  if (is.null(opt)) {
    stop("Optimization failed after ", max_attempts, " attempts. Try increasing maxeval.")
  }

  result_train <- optim(
    par = opt,
    fn = BP_GAS_FZ_Loss_Diag_Final,
    Y1 = Y1_train,
    Y2 = Y2_train,
    w1 = w1,
    w2 = w2,
    alpha = alpha,
    ticker1 = ticker1,
    ticker2 = ticker2,
    tau = -1,
    return_full = FALSE,
    method = "BFGS",
    control = list(maxit = 2000, trace = 1, reltol = 1e-8)
  )

  opt <- result_train$par

  b <- -exp(opt[1])
  c <- pnorm(opt[2])
  a <- c * b

  beta1  <- max(0.7, pnorm(opt[3]) * 2 - 1)
  beta2  <- max(0.7, pnorm(opt[4]) * 2 - 1)
  beta12 <- max(0.7, pnorm(opt[5]) * 2 - 1)

  alpha1  <- -exp(opt[6])
  alpha2  <- -exp(opt[7])
  alpha12 <- -exp(opt[8])

  cat("Estimated parameters:\n")
  cat(sprintf(" a = %.4f, b = %.4f (c = %.4f)\n", a, b, c))
  cat(sprintf(" β = diag(%.4f, %.4f, %.4f)\n", beta1, beta2, beta12))
  cat(sprintf(" α = diag(%.4f, %.4f, %.4f)\n", alpha1, alpha2, alpha12))
  cat("Final FZ loss:", result_train$value, "\n")

  oos_results <- BP_GAS_FZ_Diag_Forecast(
    theta_estimated = opt,
    Y1_train = Y1_train,
    Y2_train = Y2_train,
    Y1_oos = Y1_oos,
    Y2_oos = Y2_oos,
    w1 = w1,
    w2 = w2,
    alpha = alpha,
    ticker1 = ticker1,
    ticker2 = ticker2
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

  cat("Volatility (standard deviation):\n")
  cat(ticker1, ":", sd(Y1_train), ticker2, ":", sd(Y2_train), "\n")
  cat("Correlation:", cor(Y1_train, Y2_train), "\n")

  result <- list(
    parameters = list(
      a = a,
      b = b,
      c = c,
      beta = diag(c(beta1, beta2, beta12)),
      alpha = diag(c(alpha1, alpha2, alpha12))
    ),
    fz_loss = result_train$value,
    oos_results = oos_results,
    Yp_oos = Yp_oos,
    Y1_train = Y1_train,
    Y2_train = Y2_train,
    eval_results = eval_results,
    summary = function() print_summary(result, ticker1, ticker2)
  )

  result
}
