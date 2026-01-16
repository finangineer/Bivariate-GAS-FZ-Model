# Block-B BP-GAS-FZ: one-step-ahead OOS forecasts for VaR/ES and per-period FZ0 loss
# Returns: v_oos, e_oos, loss_oos, rho_oos
BP_GAS_FZ_BlockB_Forecast <- function(theta_estimated,
                                      Y1_train, Y2_train,
                                      Y1_oos, Y2_oos,
                                      w1, w2,
                                      alpha,
                                      ticker1, ticker2) {

  # -----------------------
  # Parameter transformation
  # -----------------------
  b <- -exp(theta_estimated[1])
  c <- pnorm(theta_estimated[2])
  a <- c * b

  # Basic validity checks
  Tn_oos <- length(Y1_oos)
  if (!is.finite(a) || !is.finite(b) || a == 0 || b >= 0) {
    return(list(
      v_oos = rep(NA_real_, Tn_oos),
      e_oos = rep(NA_real_, Tn_oos),
      loss_oos = rep(1e10, Tn_oos),
      rho_oos = rep(NA_real_, Tn_oos)
    ))
  }

  # Block-B B matrix (3x3): [b11 b12 0; b21 b22 0; 0 0 b33]
  b11 <- max(0.7, pnorm(theta_estimated[3]) * 2 - 1)
  b12 <- 0.95 * tanh(theta_estimated[4])
  b21 <- 0.95 * tanh(theta_estimated[5])
  b22 <- max(0.7, pnorm(theta_estimated[6]) * 2 - 1)
  b33 <- max(0.7, pnorm(theta_estimated[7]) * 2 - 1)

  Bmat <- matrix(
    c(b11, b12, 0,
      b21, b22, 0,
      0,   0,   b33),
    nrow = 3, byrow = TRUE
  )

  # Full A matrix (9 parameters), negative entries
  A_vec <- -exp(theta_estimated[8:16])
  Amat <- matrix(A_vec, nrow = 3, ncol = 3, byrow = TRUE)

  Tn_train <- length(Y1_train)

  # State path includes one extra for t+1 updates
  f_t <- matrix(0, nrow = Tn_train + Tn_oos + 1, ncol = 3)

  # -------------------------
  # In-sample recursion pass
  # -------------------------
  for (t in 1:Tn_train) {

    f1  <- f_t[t, 1]
    f2  <- f_t[t, 2]
    f12 <- f_t[t, 3]

    sigma1_sq <- max(exp(2 * f1), 1e-6)
    sigma2_sq <- max(exp(2 * f2), 1e-6)
    rho <- tanh(f12)

    sigma1  <- sqrt(sigma1_sq)
    sigma2  <- sqrt(sigma2_sq)
    sigma12 <- rho * sigma1 * sigma2

    Q <- w1^2 * sigma1_sq + w2^2 * sigma2_sq + 2 * w1 * w2 * sigma12
    Q <- max(Q, 1e-6)

    sqrtQ <- sqrt(Q)
    v <- a * sqrtQ
    e <- b * sqrtQ

    # If something explodes numerically, stop early with heavy penalty later
    if (!is.finite(v) || !is.finite(e) || e >= 0) {
      f_t[t + 1, ] <- f_t[t, ]
      next
    }

    Yp <- w1 * Y1_train[t] + w2 * Y2_train[t]
    I_t <- (Yp <= v)

    shock <- 1 - (Yp / (alpha * e)) * I_t

    nabla1  <- ((w1^2 * sigma1_sq + w1 * w2 * sigma12) / Q) * shock
    nabla2  <- ((w2^2 * sigma2_sq + w1 * w2 * sigma12) / Q) * shock
    nabla12 <- ((w1 * w2 * (1 - rho^2) * sigma1 * sigma2) / Q) * shock

    nabla <- c(nabla1, nabla2, nabla12)

    f_next <- Bmat %*% f_t[t, ] + Amat %*% nabla
    f_t[t + 1, ] <- pmin(pmax(as.numeric(f_next), -10), 10)
  }

  # -------------------------
  # Out-of-sample forecasting
  # -------------------------
  v_oos <- numeric(Tn_oos)
  e_oos <- numeric(Tn_oos)
  loss_oos <- numeric(Tn_oos)
  rho_oos <- numeric(Tn_oos)

  for (t in 1:Tn_oos) {

    idx <- Tn_train + t

    f1  <- f_t[idx, 1]
    f2  <- f_t[idx, 2]
    f12 <- f_t[idx, 3]

    sigma1_sq <- max(exp(2 * f1), 1e-6)
    sigma2_sq <- max(exp(2 * f2), 1e-6)
    rho <- tanh(f12)
    rho_oos[t] <- rho

    sigma1  <- sqrt(sigma1_sq)
    sigma2  <- sqrt(sigma2_sq)
    sigma12 <- rho * sigma1 * sigma2

    Q <- w1^2 * sigma1_sq + w2^2 * sigma2_sq + 2 * w1 * w2 * sigma12
    Q <- max(Q, 1e-6)

    sqrtQ <- sqrt(Q)

    v_oos[t] <- a * sqrtQ
    e_oos[t] <- b * sqrtQ

    Yp <- w1 * Y1_oos[t] + w2 * Y2_oos[t]
    I_t <- (Yp <= v_oos[t])

    loss_oos[t] <- -1 / (alpha * e_oos[t]) * I_t * (v_oos[t] - Yp) +
      v_oos[t] / e_oos[t] + log(-e_oos[t]) - 1

    # State update for next period
    shock <- 1 - (Yp / (alpha * e_oos[t])) * I_t

    nabla1  <- ((w1^2 * sigma1_sq + w1 * w2 * sigma12) / Q) * shock
    nabla2  <- ((w2^2 * sigma2_sq + w1 * w2 * sigma12) / Q) * shock
    nabla12 <- ((w1 * w2 * (1 - rho^2) * sigma1 * sigma2) / Q) * shock

    nabla <- c(nabla1, nabla2, nabla12)

    f_next <- Bmat %*% f_t[idx, ] + Amat %*% nabla
    f_t[idx + 1, ] <- pmin(pmax(as.numeric(f_next), -10), 10)
  }

  list(
    v_oos = v_oos,
    e_oos = e_oos,
    loss_oos = loss_oos,
    rho_oos = rho_oos
  )
}
