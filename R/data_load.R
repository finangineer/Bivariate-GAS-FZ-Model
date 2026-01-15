# Packages used in this script
library(quantmod)
library(nloptr)
library(zoo)
library(rugarch)

# Download two Yahoo Finance series and compute aligned close-to-close log returns.
# Returns numeric vectors Y1, Y2 and the corresponding Date index.
load_asset_data <- function(ticker1, ticker2, start_date, end_date) {

  suppressWarnings(
    getSymbols(
      Symbols = c(ticker1, ticker2),
      src = "yahoo",
      from = start_date,
      to = end_date,
      auto.assign = TRUE
    )
  )

  # quantmod strips '^' from object names (e.g., "^GSPC" -> "GSPC")
  obj1 <- gsub("\\^", "", ticker1)
  obj2 <- gsub("\\^", "", ticker2)

  if (!exists(obj1, inherits = TRUE) || !exists(obj2, inherits = TRUE)) {
    stop("Yahoo download failed (objects not found). Check ticker symbols and date range.")
  }

  p1 <- Cl(get(obj1, inherits = TRUE))
  p2 <- Cl(get(obj2, inherits = TRUE))

  prices <- merge(p1, p2, join = "inner")
  colnames(prices) <- c("Asset1", "Asset2")

  rets <- na.omit(diff(log(prices)))
  Y1 <- as.numeric(rets$Asset1)
  Y2 <- as.numeric(rets$Asset2)
  dates <- index(rets)

  cat("Total observations:", length(Y1), "\n")
  cat("Summary for", ticker1, ":\n")
  print(summary(Y1))
  cat("Summary for", ticker2, ":\n")
  print(summary(Y2))

  list(Y1 = Y1, Y2 = Y2, dates = dates)
}
