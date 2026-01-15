# Full-A BP-GAS-FZ: one-step-ahead OOS forecasts for VaR/ES and per-period FZ0 loss
BP_GAS_FZ_FullA_Forecast <- function(theta_estimated,
                                     Y1_train, Y2_train,
                                     Y1_oos, Y2_oos,
                                     w1, w2,
                                     alpha,
                                     ticker1, ticker2) {

  # Parameter mapping
  b <- -exp(theta_estimated[1])
  c <- pnorm(theta_estimated[2])
  a <- c * b

  beta1  <- max(0.7, pnorm(theta_estimated[3]) * 2 - 1)
  beta2  <- max(0.7, pnorm(theta_estimated[4]) * 2 - 1)
  beta12 <- max(0.7, pnorm(theta_estimated[5]) * 2 - 1)
  Bmat <- diag(c(beta1, beta2, beta12))

  A_vec <- -exp(theta_estimated[6:14])
  Amat <- matrix(A_vec, nrow = 3, ncol = 3, byrow = TRUE)

  Tn_train <- length(Y1_train)
  Tn_oos <- length(Y1_oos)

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

    sigma1 <- sqrt(sigma1_sq)
    sigma2 <- sqrt(sigma2_sq)
    sigma12 <- rho * sigma1 * sigma2

    Q <- w1^2 * sigma1_sq + w2^2 * sigma2_sq + 2 * w1 * w2 * sigma12
    Q <- max(Q, 1e-6)

    sqrtQ <- sqrt(Q)
    v <- a * sqrtQ
    e <- b * sqrtQ

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

  for (t in 1:Tn_oos) {

    idx <- Tn_train + t

    f1  <- f_t[idx, 1]
    f2  <- f_t[idx, 2]
    f12 <- f_t[idx, 3]

    sigma1_sq <- max(exp(2 * f1), 1e-6)
    sigma2_sq <- max(exp(2 * f2), 1e-6)
    rho <- tanh(f12)

    sigma1 <- sqrt(sigma1_sq)
    sigma2 <- sqrt(sigma2_sq)
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

  list(v_oos = v_oos, e_oos = e_oos, loss_oos = loss_oos)
}
