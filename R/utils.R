# Compute sample (historical) VaR/ES for one or more tail levels.
# Returns a matrix with columns: VaR, ES (rows correspond to alpha values).
sampleve <- function(data, alpha) {
  vehat <- matrix(0, nrow = length(alpha), ncol = 2)
  colnames(vehat) <- c("VaR", "ES")

  for (i in seq_along(alpha)) {
    vehat[i, 1] <- quantile(data, alpha[i])
    vehat[i, 2] <- mean(data[data < vehat[i, 1]])
  }

  vehat
}
