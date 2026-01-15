# Example usage: Diagonal BP-GAS-FZ on selected index pairs

result_diag_sp_ftse <- run_diagonal_model(
  ticker1 = "^GSPC",
  ticker2 = "^FTSE",
  start_date = "2015-01-01",
  end_date = "2022-12-31",
  w1 = 0.5,
  w2 = 0.5,
  alpha = 0.05,
  maxeval = 5000
)

result_diag_sp_n225 <- run_diagonal_model(
  ticker1 = "^GSPC",
  ticker2 = "^N225",
  start_date = "2015-01-01",
  end_date = "2022-12-31",
  w1 = 0.5,
  w2 = 0.5,
  alpha = 0.05,
  maxeval = 6000
)

result_diag_ft_fchi <- run_diagonal_model(
  ticker1 = "^FTSE",
  ticker2 = "^FCHI",
  start_date = "2015-01-01",
  end_date = "2022-12-31",
  w1 = 0.5,
  w2 = 0.5,
  alpha = 0.05,
  maxeval = 5000
)

result_diag_ft_n225 <- run_diagonal_model(
  ticker1 = "^FTSE",
  ticker2 = "^N225",
  start_date = "2015-01-01",
  end_date = "2022-12-31",
  w1 = 0.5,
  w2 = 0.5,
  alpha = 0.05,
  maxeval = 5000
)

result_diag_sp_ftse$summary()
result_diag_sp_n225$summary()
result_diag_ft_fchi$summary()
result_diag_ft_n225$summary()
