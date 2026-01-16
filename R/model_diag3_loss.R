# 3-Asset Diagonal BP-GAS-FZ objective (14 parameters)
# State vector: (f1, f2, f3, f12, f13, f23)
# B diagonal (6), A diagonal (6). Returns mean FZ0 loss (or full paths).

BP_GAS_FZ_Loss_Diag_3Asset <- function(theta_unconstrained,
                                       Y1, Y2, Y3,
                                       w1, w2, w3,
                                       alpha,
                                       ticker1 = "", ticker2 = "", ticker3 = "",
                                       tau = -1,
                                       return_full = FALSE,
                                       target_rho = c(0.62, 0.26, 0.40),
                                       rho_penalty_weight = 0.02) {

  # -----------------------
  # Parameter transformation
  # -----------------------
  b <- -exp(theta_unconstrained[1])
  c <- pnorm(theta_unconstrained[2])
  a <- c * b

  if (!is.finite(a) || !is.finite(b) || a == 0 || b >= 0) return(1e10)

  # Diagonal B (6): beta1,beta2,beta3,beta12,beta13,beta23 >= 0.7
  beta <- pmax(0.7, pnorm(theta_unconstrained[3:8]) * 2 - 1)
  Bmat <- diag(beta)

  # Diagonal A (6): negative via -exp(.)
  alpha_params <- -exp(theta_unconstrained[9:14])
  Amat <- diag(alpha_params)

  # -----------------------
  # Penalties (static)
  # -----------------------
  es_var_ratio <- abs(b / a)
  if (!is.finite(es_var_ratio)) return(1e10)

  es_var_penalty <- if (es_var_ratio < 1.3) 0.06 * (1.3 - es_var_ratio) else 0

  # Encourage vol betas (1:3) and corr betas (4:6) not to be too small
  beta_penalty <- sum(pmax(0.8 - beta[1:3], 0)^2) * 100 +
    sum(pmax(0.8 - beta[4:6], 0)^2) * 50

  # Your “target” correlation-beta preference (kept as you wrote)
  beta_corr_penalty <- 50 * ((beta[4] - 0.9)^2 + (beta[5] - 0.7)^2 + (beta[6] - 0.8)^2)

  # -----------------------
  # Recursion + per-period loss
  # -----------------------
  Tn <- length(Y1)

  f_t <- matrix(0, nrow = Tn + 1, ncol = 6)
  Q_t <- v_t <- e_t <- loss_t <- numeric(Tn)
  nabla <- matrix(0, nrow = Tn, ncol = 6)

  breaches <- 0
  rho_sum <- 0

  for (t in 1:Tn) {

    f1  <- f_t[t, 1]
    f2  <- f_t[t, 2]
    f3  <- f_t[t, 3]
    f12 <- f_t[t, 4]
    f13 <- f_t[t, 5]
    f23 <- f_t[t, 6]

    sigma1_sq <- max(exp(2 * f1), 1e-6)
    sigma2_sq <- max(exp(2 * f2), 1e-6)
    sigma3_sq <- max(exp(2 * f3), 1e-6)

    rho12 <- tanh(f12)
    rho13 <- tanh(f13)
    rho23 <- tanh(f23)

    sigma1 <- sqrt(sigma1_sq)
    sigma2 <- sqrt(sigma2_sq)
    sigma3 <- sqrt(sigma3_sq)

    sigma12 <- rho12 * sigma1 * sigma2
    sigma13 <- rho13 * sigma1 * sigma3
    sigma23 <- rho23 * sigma2 * sigma3

    # Correlation target penalty accumulator (defaults match your hard-coded values)
    rho_sum <- rho_sum +
      (rho12 - target_rho[1])^2 +
      (rho13 - target_rho[2])^2 +
      (rho23 - target_rho[3])^2

    Q <- w1^2 * sigma1_sq + w2^2 * sigma2_sq + w3^2 * sigma3_sq +
      2 * w1 * w2 * sigma12 +
      2 * w1 * w3 * sigma13 +
      2 * w2 * w3 * sigma23
    Q <- max(Q, 1e-6)

    sqrtQ <- sqrt(Q)
    v <- a * sqrtQ
    e <- b * sqrtQ

    if (!is.finite(v) || !is.finite(e) || e >= 0 || sqrtQ < 1e-6) {
      return(1e10 + beta_penalty + beta_corr_penalty)
    }

    Yp <- w1 * Y1[t] + w2 * Y2[t] + w3 * Y3[t]

    I_t <- if (tau == -1) (Yp <= v) else 1 / (1 + exp(tau * (Yp - v)))
    if (!is.finite(I_t) || length(I_t) != 1) {
      return(1e10 + beta_penalty + beta_corr_penalty)
    }

    loss <- -1 / (alpha * e) * I_t * (v - Yp) + v / e + log(-e) - 1
    loss <- loss + es_var_penalty + beta_penalty + beta_corr_penalty

    if (Yp < v) breaches <- breaches + 1
    if (!is.finite(loss)) return(1e10 + beta_penalty + beta_corr_penalty)

    shock <- 1 - (Yp / (alpha * e)) * I_t

    nabla1  <- ((w1^2 * sigma1_sq + w1 * w2 * sigma12 + w1 * w3 * sigma13) / Q) * shock
    nabla2  <- ((w2^2 * sigma2_sq + w1 * w2 * sigma12 + w2 * w3 * sigma23) / Q) * shock
    nabla3  <- ((w3^2 * sigma3_sq + w1 * w3 * sigma13 + w2 * w3 * sigma23) / Q) * shock
    nabla12 <- ((w1 * w2 * (1 - rho12^2) * sigma1 * sigma2) / Q) * shock
    nabla13 <- ((w1 * w3 * (1 - rho13^2) * sigma1 * sigma3) / Q) * shock
    nabla23 <- ((w2 * w3 * (1 - rho23^2) * sigma2 * sigma3) / Q) * shock

    nabla[t, ] <- c(nabla1, nabla2, nabla3, nabla12, nabla13, nabla23)

    if (any(!is.finite(nabla[t, ])) || any(abs(nabla[t, ]) > 1e5)) {
      return(1e10 + beta_penalty + beta_corr_penalty)
    }

    f_next <- Bmat %*% f_t[t, ] + Amat %*% nabla[t, ]
    f_t[t + 1, ] <- pmin(pmax(as.numeric(f_next), -10), 10)

    Q_t[t] <- Q
    v_t[t] <- v
    e_t[t] <- e
    loss_t[t] <- loss
  }

  # Breach penalty (kept as you wrote)
  target_breaches <- 0.05 * Tn
  if (breaches > target_breaches * 1.1) {
    loss_t <- loss_t + 0.08 * (breaches - target_breaches * 1.1)
  }

  # Correlation penalty (average)
  loss_t <- loss_t + rho_penalty_weight * (rho_sum / Tn)

  avg_loss <- mean(loss_t)
  if (!is.finite(avg_loss)) return(1e10 + beta_penalty + beta_corr_penalty)

  if (!return_full) return(avg_loss)

  rho_ts <- cbind(
    rho12 = tanh(f_t[-1, 4]),
    rho13 = tanh(f_t[-1, 5]),
    rho23 = tanh(f_t[-1, 6])
  )

  list(
    average_loss = avg_loss,
    f_ts = f_t[-1, ],
    Q_ts = Q_t,
    v_ts = v_t,
    e_t  = e_t,
    loss_ts = loss_t,
    nabla_ts = nabla,
    rho_ts = rho_ts,
    breaches = breaches
  )
}
