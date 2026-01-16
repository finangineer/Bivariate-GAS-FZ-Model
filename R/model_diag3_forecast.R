# 3-Asset Diagonal BP-GAS-FZ: one-step-ahead OOS forecasts for VaR/ES and per-period FZ0 loss
# Returns: v_oos, e_oos, loss_oos, rho_oos (list + rho_mat convenience)

BP_GAS_FZ_Diag_3Asset_Forecast <- function(theta_estimated,
                                           Y1_train, Y2_train, Y3_train,
                                           Y1_oos, Y2_oos, Y3_oos,
                                           w1, w2, w3,
                                           alpha,
                                           ticker1, ticker2, ticker3) {

  b <- -exp(theta_estimated[1])
  c <- pnorm(theta_estimated[2])
  a <- c * b

  Tn_oos <- length(Y1_oos)

  if (!is.finite(a) || !is.finite(b) || a == 0 || b >= 0) {
    return(list(
      v_oos = rep(NA_real_, Tn_oos),
      e_oos = rep(NA_real_, Tn_oos),
      loss_oos = rep(1e10, Tn_oos),
      rho_oos = list(
        rho12 = rep(NA_real_, Tn_oos),
        rho13 = rep(NA_real_, Tn_oos),
        rho23 = rep(NA_real_, Tn_oos)
      ),
      rho_mat = cbind(
        rho12 = rep(NA_real_, Tn_oos),
        rho13 = rep(NA_real_, Tn_oos),
        rho23 = rep(NA_real_, Tn_oos)
      )
    ))
  }

  beta <- pmax(0.7, pnorm(theta_estimated[3:8]) * 2 - 1)
  Bmat <- diag(beta)

  alpha_params <- -exp(theta_estimated[9:14])
  Amat <- diag(alpha_params)

  Tn_train <- length(Y1_train)

  # State path: (f1,f2,f3,f12,f13,f23)
  f_t <- matrix(0, nrow = Tn_train + Tn_oos + 1, ncol = 6)

  # -------------------------
  # In-sample recursion pass
  # -------------------------
  for (t in 1:Tn_train) {

    f1  <- f_t[t, 1]; f2  <- f_t[t, 2]; f3  <- f_t[t, 3]
    f12 <- f_t[t, 4]; f13 <- f_t[t, 5]; f23 <- f_t[t, 6]

    sigma1_sq <- max(exp(2 * f1), 1e-6)
    sigma2_sq <- max(exp(2 * f2), 1e-6)
    sigma3_sq <- max(exp(2 * f3), 1e-6)

    rho12 <- tanh(f12); rho13 <- tanh(f13); rho23 <- tanh(f23)

    sigma1 <- sqrt(sigma1_sq); sigma2 <- sqrt(sigma2_sq); sigma3 <- sqrt(sigma3_sq)

    sigma12 <- rho12 * sigma1 * sigma2
    sigma13 <- rho13 * sigma1 * sigma3
    sigma23 <- rho23 * sigma2 * sigma3

    Q <- w1^2 * sigma1_sq + w2^2 * sigma2_sq + w3^2 * sigma3_sq +
      2 * w1 * w2 * sigma12 +
      2 * w1 * w3 * sigma13 +
      2 * w2 * w3 * sigma23
    Q <- max(Q, 1e-6)

    sqrtQ <- sqrt(Q)
    v <- a * sqrtQ
    e <- b * sqrtQ

    if (!is.finite(v) || !is.finite(e) || e >= 0) {
      f_t[t + 1, ] <- f_t[t, ]
      next
    }

    Yp <- w1 * Y1_train[t] + w2 * Y2_train[t] + w3 * Y3_train[t]
    I_t <- as.numeric(Yp <= v)

    shock <- 1 - (Yp / (alpha * e)) * I_t

    nabla1  <- ((w1^2 * sigma1_sq + w1 * w2 * sigma12 + w1 * w3 * sigma13) / Q) * shock
    nabla2  <- ((w2^2 * sigma2_sq + w1 * w2 * sigma12 + w2 * w3 * sigma23) / Q) * shock
    nabla3  <- ((w3^2 * sigma3_sq + w1 * w3 * sigma13 + w2 * w3 * sigma23) / Q) * shock
    nabla12 <- ((w1 * w2 * (1 - rho12^2) * sigma1 * sigma2) / Q) * shock
    nabla13 <- ((w1 * w3 * (1 - rho13^2) * sigma1 * sigma3) / Q) * shock
    nabla23 <- ((w2 * w3 * (1 - rho23^2) * sigma2 * sigma3) / Q) * shock

    nabla <- c(nabla1, nabla2, nabla3, nabla12, nabla13, nabla23)

    f_next <- Bmat %*% f_t[t, ] + Amat %*% nabla
    f_t[t + 1, ] <- pmin(pmax(as.numeric(f_next), -10), 10)
  }

  # -------------------------
  # Out-of-sample forecasting
  # -------------------------
  v_oos <- numeric(Tn_oos)
  e_oos <- numeric(Tn_oos)
  loss_oos <- numeric(Tn_oos)

  rho12_oos <- numeric(Tn_oos)
  rho13_oos <- numeric(Tn_oos)
  rho23_oos <- numeric(Tn_oos)

  for (t in 1:Tn_oos) {

    idx <- Tn_train + t

    f1  <- f_t[idx, 1]; f2  <- f_t[idx, 2]; f3  <- f_t[idx, 3]
    f12 <- f_t[idx, 4]; f13 <- f_t[idx, 5]; f23 <- f_t[idx, 6]

    sigma1_sq <- max(exp(2 * f1), 1e-6)
    sigma2_sq <- max(exp(2 * f2), 1e-6)
    sigma3_sq <- max(exp(2 * f3), 1e-6)

    rho12 <- tanh(f12); rho13 <- tanh(f13); rho23 <- tanh(f23)

    rho12_oos[t] <- rho12
    rho13_oos[t] <- rho13
    rho23_oos[t] <- rho23

    sigma1 <- sqrt(sigma1_sq); sigma2 <- sqrt(sigma2_sq); sigma3 <- sqrt(sigma3_sq)

    sigma12 <- rho12 * sigma1 * sigma2
    sigma13 <- rho13 * sigma1 * sigma3
    sigma23 <- rho23 * sigma2 * sigma3

    Q <- w1^2 * sigma1_sq + w2^2 * sigma2_sq + w3^2 * sigma3_sq +
      2 * w1 * w2 * sigma12 +
      2 * w1 * w3 * sigma13 +
      2 * w2 * w3 * sigma23
    Q <- max(Q, 1e-6)

    sqrtQ <- sqrt(Q)

    v_oos[t] <- a * sqrtQ
    e_oos[t] <- b * sqrtQ

    Yp <- w1 * Y1_oos[t] + w2 * Y2_oos[t] + w3 * Y3_oos[t]
    I_t <- as.numeric(Yp <= v_oos[t])

    loss_oos[t] <- -1 / (alpha * e_oos[t]) * I_t * (v_oos[t] - Yp) +
      v_oos[t] / e_oos[t] + log(-e_oos[t]) - 1

    shock <- 1 - (Yp / (alpha * e_oos[t])) * I_t

    nabla1  <- ((w1^2 * sigma1_sq + w1 * w2 * sigma12 + w1 * w3 * sigma13) / Q) * shock
    nabla2  <- ((w2^2 * sigma2_sq + w1 * w2 * sigma12 + w2 * w3 * sigma23) / Q) * shock
    nabla3  <- ((w3^2 * sigma3_sq + w1 * w3 * sigma13 + w2 * w3 * sigma23) / Q) * shock
    nabla12 <- ((w1 * w2 * (1 - rho12^2) * sigma1 * sigma2) / Q) * shock
    nabla13 <- ((w1 * w3 * (1 - rho13^2) * sigma1 * sigma3) / Q) * shock
    nabla23 <- ((w2 * w3 * (1 - rho23^2) * sigma2 * sigma3) / Q) * shock

    nabla <- c(nabla1, nabla2, nabla3, nabla12, nabla13, nabla23)

    f_next <- Bmat %*% f_t[idx, ] + Amat %*% nabla
    f_t[idx + 1, ] <- pmin(pmax(as.numeric(f_next), -10), 10)
  }

  rho_mat <- cbind(rho12 = rho12_oos, rho13 = rho13_oos, rho23 = rho23_oos)

  list(
    v_oos = v_oos,
    e_oos = e_oos,
    loss_oos = loss_oos,
    rho_oos = list(rho12 = rho12_oos, rho13 = rho13_oos, rho23 = rho23_oos),
    rho_mat = rho_mat
  )
}
