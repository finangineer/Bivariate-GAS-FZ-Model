# Example usage: Full-A BP-GAS-FZ

result_fullA_sp_ftse <- run_fullA_model("^GSPC", "^FTSE", "2015-01-01", "2022-12-31", 0.5, 0.5, 0.05, maxeval = 7000)
result_fullA_sp_n225 <- run_fullA_model("^GSPC", "^N225", "2015-01-01", "2022-12-31", 0.5, 0.5, 0.05, maxeval = 7000)
result_fullA_ft_fchi <- run_fullA_model("^FTSE", "^FCHI", "2015-01-01", "2022-12-31", 0.5, 0.5, 0.05, maxeval = 7000)
result_fullA_ft_n225 <- run_fullA_model("^FTSE", "^N225", "2015-01-01", "2022-12-31", 0.5, 0.5, 0.05, maxeval = 7000)

result_fullA_sp_ftse$summary()
result_fullA_sp_n225$summary()
result_fullA_ft_fchi$summary()
result_fullA_ft_n225$summary()
