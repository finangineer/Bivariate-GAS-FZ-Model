# Full-A BP-GAS-FZ objective:
# - B is diagonal (beta1, beta2, beta12), truncated at 0.7
# - A is full 3x3 (9 parameters), negative via -exp(.)
# Returns mean per-period FZ0 loss (or full paths if return_full=TRUE).
BP_GAS_FZ_Loss_FullA_Final <- function(theta_unconstrained,
                                       Y1, Y2,
                                       w1, w2,
                                       alpha,
                                       ticker1 = "",
                                       ticker2 = "",
                                       tau = -1,
                                       return_full = FALSE) {

  # -----------------------
  # Parameter transformation
  # -----------------------
  b <- -exp(theta_unconstrained[1])
  c <- pnorm(theta_unconstrained[2])
  a <- c * b

  beta1  <- max(0.7, pnorm(theta_unconstrained[3]) * 2 - 1)
  beta2  <- max(0.7, pnorm(theta_unconstrained[4]) * 2 - 1)
  beta12 <- max(0.7, pnorm(theta_unconstrained[5]) * 2 - 1)
  Bmat <- diag(c(beta1, beta2, beta12))

  A_vec <- -exp(theta_unconstrained[6:14])
  Amat <- matrix(A_vec, nrow = 3, ncol = 3, byrow = TRUE)

  # -----------------------
  # Penalties (as in original)
  # -----------------------
  beta_diff_penalty <- 0
  if (abs(beta1 - beta2) < 0.05 || abs(beta1 - beta12) < 0.05 || abs(beta2 - beta12) < 0.05) {
    beta_diff_penalty <- 0.02
  }

  A_offdiag_penalty <- 0.01 * sum(abs(Amat[upper.tri(Amat)]) + abs(Amat[lower.tri(Amat)]))

  # -----------------------
  # Recursion + per-period loss
  # -----------------------
  Tn <- length(Y1)

  f_t <- matrix(0, nrow = Tn + 1, ncol = 3)
  Q_t <- v_t <- e_t <- loss_t <- numeric(Tn)
  nabla <- matrix(0, nrow = Tn, ncol = 3)

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

    # ES/VaR ratio penalty (|b/a|)
    es_var_ratio <- abs(b / a)
    if (es_var_ratio < 1.3) {
      loss <- loss + 0.05 * (1.3 - es_var_ratio)
    }

    loss <- loss + beta_diff_penalty + A_offdiag_penalty

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
