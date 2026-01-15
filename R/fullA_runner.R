# ------------------------------------------------------------
# Full-A BP-GAS-FZ runner
# - B is diagonal (beta1, beta2, beta12)
# - A is full 3x3 (9 parameters)
# Estimation: nloptr global search + BFGS refinement
# ------------------------------------------------------------

run_fullA_model <- function(ticker1, ticker2,
                            start_date, end_date,
                            w1, w2,
                            alpha,
                            maxeval = 7000) {

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

  # GARCH fits (computed for initialization diagnostics; not used directly below)
  garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0)),
    distribution.model = "std"
  )
  garch_fit1 <- ugarchfit(spec = garch_spec, data = Y1_train)
  garch_fit2 <- ugarchfit(spec = garch_spec, data = Y2_train)
  beta1_init <- coef(garch_fit1)["beta1"]
  beta2_init <- coef(garch_fit2)["beta1"]

  # Starting values based on portfolio historical VaR/ES + conservative beta starts
  Yp_train <- w1 * Y1_train + w2 * Y2_train
  vebar <- sampleve(Yp_train, alpha)

  theta_start <- c(
    log(-vebar[1, 2]),             # theta_b
    qnorm(vebar[1, 1] / vebar[1, 2]), # theta_c
    qnorm((0.95 + 1) / 2),         # beta1 start
    qnorm((0.85 + 1) / 2),         # beta2 start
    qnorm((0.75 + 1) / 2),         # beta12 start
    rep(log(0.02), 9)              # A entries: small negative via -exp(.)
  )

  # Bounds (as in your code)
  lower <- c(-5, -5, rep(1.0364, 3), rep(-6, 9))
  upper <- c(-1.5, 0.524, rep(2, 3), rep(-2, 9))
  theta_start <- pmin(pmax(theta_start, lower), upper)

  # Wrapper for nloptr
  eval_f_wrapper <- function(x, Y1, Y2, w1, w2, alpha, tau) {
    BP_GAS_FZ_Loss_FullA_Final(
      x,
      Y1 = Y1, Y2 = Y2,
      w1 = w1, w2 = w2,
      alpha = alpha,
      ticker1 = ticker1, ticker2 = ticker2,
      tau = tau,
      return_full = FALSE
    )
  }

  # Global search (CRS2-LM) with a few retries
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

  if (is.null(opt)) {
    stop("Optimization failed. Consider increasing maxeval or adjusting starting values.")
  }

  # Local refinement (BFGS)
  result_train <- optim(
    par = opt,
    fn = BP_GAS_FZ_Loss_FullA_Final,
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

  # Extract parameters (same mapping)
  b <- -exp(opt[1])
  c <- pnorm(opt[2])
  a <- c * b

  beta1  <- max(0.7, pnorm(opt[3]) * 2 - 1)
  beta2  <- max(0.7, pnorm(opt[4]) * 2 - 1)
  beta12 <- max(0.7, pnorm(opt[5]) * 2 - 1)

  A_vec <- -exp(opt[6:14])
  Amat <- matrix(A_vec, nrow = 3, ncol = 3, byrow = TRUE)

  cat("Estimated parameters:\n")
  cat(sprintf(" a = %.4f, b = %.4f (c = %.4f)\n", a, b, c))
  cat(sprintf(" β = diag(%.4f, %.4f, %.4f)\n", beta1, beta2, beta12))
  cat(" α (A matrix):\n")
  print(round(Amat, 4))
  cat("Final FZ loss:", result_train$value, "\n")

  # Out-of-sample forecasts
  oos_results <- BP_GAS_FZ_FullA_Forecast(
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

  cat("Volatility:", sd(Y1_train), sd(Y2_train), "\n")
  cat("Correlation:", cor(Y1_train, Y2_train), "\n")

  result <- list(
    parameters = list(
      a = a,
      b = b,
      c = c,
      beta = diag(c(beta1, beta2, beta12)),
      alpha = Amat
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
