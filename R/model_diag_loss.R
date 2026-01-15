# Diagonal BP-GAS-FZ: objective (average FZ0 loss) with constraints/penalties
BP_GAS_FZ_Loss_Diag_Final <- function(theta_unconstrained,
                                      Y1, Y2,
                                      w1, w2,
                                      alpha,
                                      ticker1, ticker2,
                                      tau = -1,
                                      return_full = FALSE) {

  # Map unconstrained parameters to constrained region
  b <- -exp(theta_unconstrained[1])
  c <- pnorm(theta_unconstrained[2])
  a <- c * b

  beta1  <- max(0.7, pnorm(theta_unconstrained[3]) * 2 - 1)
  beta2  <- max(0.7, pnorm(theta_unconstrained[4]) * 2 - 1)
  beta12 <- max(0.7, pnorm(theta_unconstrained[5]) * 2 - 1)
  Bmat <- diag(c(beta1, beta2, beta12))

  alpha1  <- -exp(theta_unconstrained[6])
  alpha2  <- -exp(theta_unconstrained[7])
  alpha12 <- -exp(theta_unconstrained[8])
  Amat <- diag(c(alpha1, alpha2, alpha12))

  # Small penalty if persistence parameters are nearly identical
  beta_diff_penalty <- 0
  if (abs(beta1 - beta2) < 0.05 || abs(beta1 - beta12) < 0.05 || abs(beta2 - beta12) < 0.05) {
    beta_diff_penalty <- 0.02
  }

  Tn <- length(Y1)

  f_t <- matrix(0, nrow = Tn + 1, ncol = 3)
  Q_t <- v_t <- e_t <- loss_t <- numeric(Tn)
  nabla <- matrix(0, nrow = Tn, ncol = 3)

  breaches <- 0

  for (t in 1:Tn) {

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

    if (!is.finite(v) || !is.finite(e) || e >= 0 || sqrtQ < 1e-6) {
      return(1e10)
    }

    Yp <- w1 * Y1[t] + w2 * Y2[t]
    I_t <- if (tau == -1) (Yp <= v) else 1 / (1 + exp(tau * (Yp - v)))

    loss <- -1 / (alpha * e) * I_t * (v - Yp) + v / e + log(-e) - 1

    # Penalize unrealistically small ES/VaR ratio (|b/a|)
    es_var_ratio <- abs(b / a)
    if (es_var_ratio < 1.3) {
      loss <- loss + 0.05 * (1.3 - es_var_ratio)
    }

    if (Yp < v) {
      breaches <- breaches + 1
    }

    loss <- loss + beta_diff_penalty

    if (!is.finite(loss)) {
      return(1e10)
    }

    shock <- 1 - (Yp / (alpha * e)) * I_t

    nabla1  <- ((w1^2 * sigma1_sq + w1 * w2 * sigma12) / Q) * shock
    nabla2  <- ((w2^2 * sigma2_sq + w1 * w2 * sigma12) / Q) * shock
    nabla12 <- ((w1 * w2 * (1 - rho^2) * sigma1 * sigma2) / Q) * shock

    nabla[t, ] <- c(nabla1, nabla2, nabla12)

    if (any(!is.finite(nabla[t, ])) || any(abs(nabla[t, ]) > 1e5)) {
      return(1e10)
    }

    f_next <- Bmat %*% f_t[t, ] + Amat %*% nabla[t, ]
    f_t[t + 1, ] <- pmin(pmax(as.numeric(f_next), -10), 10)

    Q_t[t] <- Q
    v_t[t] <- v
    e_t[t] <- e
    loss_t[t] <- loss
  }

  # Pair-specific breach penalty (kept as in original code)
  if (ticker1 == "^GSPC" && ticker2 == "^N225" && breaches > 0.05 * Tn * 1.1) {
    loss_t <- loss_t + 0.01 * (breaches - 0.05 * Tn * 1.1)
  }

  avg_loss <- mean(loss_t)
  if (!is.finite(avg_loss)) {
    return(1e10)
  }

  if (!return_full) {
    return(avg_loss)
  }

  list(
    average_loss = avg_loss,
    f_ts = f_t[-1, ],
    Q_ts = Q_t,
    v_ts = v_t,
    e_t  = e_t,
    loss_ts = loss_t,
    nabla_ts = nabla
  )
}
