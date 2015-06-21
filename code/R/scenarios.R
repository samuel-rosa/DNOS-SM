# SIMULATE SCENARIOS - EASINESS ################################################

# Create scenario: sat > land > dem > soil > geo
colnames(comb_num)
sim_a <- c(1, 1, 1, 1, 1)
sim_b <- c(1, 1, 1, 2, 1)
sim_c <- c(1, 2, 1, 2, 1)
sim_d <- c(1, 2, 1, 2, 2)
sim_e <- c(2, 2, 1, 2, 2)
sim_f <- c(2, 2, 2, 2, 2)
easy_sim_n <- rbind(sim_a, sim_b, sim_c, sim_d, sim_e, sim_f)
colnames(easy_sim_n) <- colnames(comb_num)
easy_sim <- apply(easy_sim_n, 1, function (X) {
  apply(comb_num, 1, identical, X)
})
easy_sim <- which(easy_sim == TRUE, arr.ind = TRUE)
combinations[easy_sim[, 1]]
easy_clay_formulas   <- clay_formulas[easy_sim[, 1]]
easy_carbon_formulas <- carbon_formulas[easy_sim[, 1]]
easy_ecec_formulas   <- ecec_formulas[easy_sim[, 1]]
# clay
rm(data)
data <- cal_data_1@data
clay_easy <- buildModel(easy_clay_formulas, data, vif = TRUE, aic = TRUE)
clay_easy_stats <- modelStats(clay_easy)
clay_easy_stats <- cbind(easy_sim_n, clay_easy_stats)
rm(data)
# carbon
data <- cal_data_1@data
carbon_easy <- buildModel(easy_carbon_formulas, data, vif = TRUE, aic = TRUE)
carbon_easy_stats <- modelStats(carbon_easy)
carbon_easy_stats <- cbind(easy_sim_n, carbon_easy_stats)
rm(data)
# ecec
data <- cal_data_1@data
ecec_easy <- buildModel(easy_ecec_formulas, data, vif = TRUE, aic = TRUE)
ecec_easy_stats <- modelStats(ecec_easy)
ecec_easy_stats <- cbind(easy_sim_n, ecec_easy_stats)
rm(data)
# latex table
Scenario <- data.frame(clay_easy_stats[, c(6, 7, 9, 12)],
                       carbon_easy_stats[, c(6, 7, 9, 12)],
                       ecec_easy_stats[, c(6, 7, 9, 12)])
colnames(Scenario) <- rep(c("df", "AIC", "RMSE", "Adj R$^2$"), 3)
rownames(Scenario) <- c("\\texttt{start}", "\\texttt{sat}", "\\texttt{land}",
                        "\\texttt{dem}", "\\texttt{soil}", "\\texttt{geo}")
Scenario[, c(2, 6, 10)] <- round(Scenario[, c(2, 6, 10)])
long_cap <- "Performance statistics of linear models built to predict soil \\texttt{clay}, \\texttt{carbon}, and \\texttt{ECEC}, for a scenario in which environmental covariates are updated according to the expected effort needed for the updating."
foot <- "* The simulated scenario is composed by updating, in chronological order, the satellite image (\\texttt{sat}), land use map (\\texttt{land}), digital elevation model (\\texttt{dem}), area-class soil map (\\texttt{soil}), and geological map (\\texttt{geo}).\\newline ** Performance statistics: df - degrees of freedom, AIC - Akaike Information Criteria, RMSE - normalized root mean square error, Adj R$^2$ - adjusted coefficient of determination.\\newline *** RMSE was normalized to the standard deviation of soil attributes (clay - 1.78, carbon - 0.57, ECEC - 0.54), while R$^2$ was adjusted to the initial number of covariates offered to the algorithm to build the model."
short_cap <- "Model performance when covariates are updated according to the expected effort."
latex(Scenario, file =  paste(update_covars_dir, "easy-scenario.tex", sep = ""),
      label = "tab:easy-scenario", table.env = TRUE, longtable = FALSE,
      cgroup = c("Clay","Carbon", "ECEC"),
      n.cgroup = c(4, 4, 4), na.blank = TRUE, ctable = TRUE,
      caption = long_cap, caption.lot = short_cap, where = "ht",
      size = "scriptsize", insert.bottom = foot)

# SIMULATE SCENARIOS - EXPECTED ADDED VALUE ####################################

# Create scenario: clay ~ SOIL > DEM  > GEO  > LAND > SAT
colnames(comb_num)
sim_a <- c(1, 1, 1, 1, 1)
sim_b <- c(2, 1, 1, 1, 1)
sim_c <- c(2, 1, 1, 1, 2)
sim_d <- c(2, 1, 2, 1, 2)
sim_e <- c(2, 2, 2, 1, 2)
sim_f <- c(2, 2, 2, 2, 2)
value_sim_n <- rbind(sim_a, sim_b, sim_c, sim_d, sim_e, sim_f)
colnames(value_sim_n) <- colnames(comb_num)
value_sim <- apply(value_sim_n, 1, function (X) {
  apply(comb_num, 1, identical, X)
})
value_sim <- which(value_sim == TRUE, arr.ind = TRUE)
combinations[value_sim[, 1]]
value_clay_formulas   <- clay_formulas[value_sim[, 1]]
# fit models
data <- cal_data_1@data
clay_value <- buildModel(value_clay_formulas, data, penalize = TRUE, vif = TRUE,
                         aic = TRUE, aic.direction = "both")
clay_value_stats <- modelStats(clay_value)
clay_value_stats <- cbind(value_sim_n, clay_value_stats)
rm(data)
# latex table
Scenario <- data.frame(clay_value_stats[, c(6, 7, 9, 12)])
colnames(Scenario) <- c("df", "AIC", "RMSE", "Adj R$^2$")
rownames(Scenario) <- c("\\texttt{start}", "\\texttt{soil}", "\\texttt{dem}",
                        "\\texttt{geo}", "\\texttt{land}", "\\texttt{sat}")
long <- "Performance statistics of linear models built to predict soil \\texttt{clay} for a scenario in which environmental covariates are updated according to the expected added value of the updating."
foot <- "* The simulated scenario is composed by updating, in chronological order, the area-class soil map (\\texttt{soil}), digital elevation model (\\texttt{dem}), geological map (\\texttt{geo}), land use map (\\texttt{land}), and satellite image (\\texttt{sat}).\\newline ** Performance statistics: df - degrees of freedom, AIC - Akaike Information Criteria, RMSE - normalized root mean square error, Adj R$^2$ - adjusted coefficient of determination.\\newline *** RMSE was normalized to the standard deviation of clay (1.78), while R$^2$ was adjusted to the initial number of covariates offered to the algorithm to build the model."
short_cap <- "Clay model performance when covariates are updated according to the expected added value."
latex(Scenario, file =  paste(update_covars_dir, "clay-value-scenario.tex", sep = ""),
      label = "tab:clay-value-scenario", table.env = TRUE, longtable = FALSE,
      na.blank = TRUE, ctable = TRUE,
      caption = long, caption.lot = short_cap, where = "ht",
      size = "scriptsize", insert.bottom = foot)

# Create scenario: carbon ~ SOIL > LAND > DEM  > SAT  > GEO
colnames(comb_num)
sim_a <- c(1, 1, 1, 1, 1)
sim_b <- c(2, 1, 1, 1, 1)
sim_c <- c(2, 2, 1, 1, 1)
sim_d <- c(2, 2, 1, 1, 2)
sim_e <- c(2, 2, 1, 2, 2)
sim_f <- c(2, 2, 2, 2, 2)
value_sim_n <- rbind(sim_a, sim_b, sim_c, sim_d, sim_e, sim_f)
colnames(value_sim_n) <- colnames(comb_num)
value_sim <- apply(value_sim_n, 1, function (X) {
  apply(comb_num, 1, identical, X)
})
value_sim <- which(value_sim == TRUE, arr.ind = TRUE)
combinations[value_sim[, 1]]
value_carbon_formulas   <- carbon_formulas[value_sim[, 1]]
# fit models
data <- cal_data_1@data
carbon_value <- buildModel(value_carbon_formulas, data, penalize = TRUE, vif = TRUE,
                           aic = TRUE, aic.direction = "both")
carbon_value_stats <- modelStats(carbon_value)
carbon_value_stats <- cbind(value_sim_n, carbon_value_stats)
rm(data)
# latex table
Scenario <- data.frame(carbon_value_stats[, c(6, 7, 9, 12)])
colnames(Scenario) <- c("df", "AIC", "RMSE", "Adj R$^2$")
rownames(Scenario) <- c("\\texttt{start}", "\\texttt{soil}", "\\texttt{land}",
                        "\\texttt{dem}", "\\texttt{sat}", "\\texttt{geo}")
long <- "Performance statistics of linear models built to predict soil \\texttt{carbon} for a scenario in which environmental covariates are updated according to the expected added value of the updating."
foot <- "* The simulated scenario is composed by updating, in chronological order, the area-class soil map (\\texttt{soil}), land use map (\\texttt{land}), digital elevation model (\\texttt{dem}), satellite image (\\texttt{sat}), and geological map (\\texttt{geo}).\\newline ** Performance statistics: df - degrees of freedom, AIC - Akaike Information Criteria, RMSE - normalized root mean square error, Adj R$^2$ - adjusted coefficient of determination.\\newline *** RMSE was normalized to the standard deviation of carbon (0.57), while R$^2$ was adjusted to the initial number of covariates offered to the algorithm to build the model."
short_cap <- "carbon model performance when covariates are updated according to the expected added value."
latex(Scenario, file =  paste(update_covars_dir, "carbon-value-scenario.tex", sep = ""),
      label = "tab:carbon-value-scenario", table.env = TRUE, longtable = FALSE,
      na.blank = TRUE, ctable = TRUE,
      caption = long, caption.lot = short_cap, where = "ht",
      size = "scriptsize", insert.bottom = foot)

# Create scenario: ecec ~ SOIL > GEO  > LAND > DEM  > SAT
colnames(comb_num)
sim_a <- c(1, 1, 1, 1, 1)
sim_b <- c(2, 1, 1, 1, 1)
sim_c <- c(2, 1, 2, 1, 1)
sim_d <- c(2, 2, 2, 1, 1)
sim_e <- c(2, 2, 2, 1, 2)
sim_f <- c(2, 2, 2, 2, 2)
value_sim_n <- rbind(sim_a, sim_b, sim_c, sim_d, sim_e, sim_f)
colnames(value_sim_n) <- colnames(comb_num)
value_sim <- apply(value_sim_n, 1, function (X) {
  apply(comb_num, 1, identical, X)
})
value_sim <- which(value_sim == TRUE, arr.ind = TRUE)
combinations[value_sim[, 1]]
value_ecec_formulas   <- ecec_formulas[value_sim[, 1]]
# fit models
data <- cal_data_1@data
ecec_value <- buildModel(value_ecec_formulas, data, penalize = TRUE, vif = TRUE,
                         aic = TRUE, aic.direction = "both")
ecec_value_stats <- modelStats(ecec_value)
ecec_value_stats <- cbind(value_sim_n, ecec_value_stats)
rm(data)
# latex table
Scenario <- data.frame(ecec_value_stats[, c(6, 7, 9, 12)])
colnames(Scenario) <- c("df", "AIC", "RMSE", "Adj R$^2$")
rownames(Scenario) <- c("\\texttt{start}", "\\texttt{soil}", "\\texttt{geo}",
                        "\\texttt{land}", "\\texttt{dem}", "\\texttt{sat}")
long <- "Performance statistics of linear models built to predict soil \\texttt{ecec} for a scenario in which environmental covariates are updated according to the expected added value of the updating."
foot <- "* The simulated scenario is composed by updating, in chronological order, the area-class soil map (\\texttt{soil}), geological map (\\texttt{geo}), land use map (\\texttt{land}), digital elevation model (\\texttt{dem}), and satellite image (\\texttt{sat}).\\newline ** Performance statistics: df - degrees of freedom, AIC - Akaike Information Criteria, RMSE - normalized root mean square error, Adj R$^2$ - adjusted coefficient of determination.\\newline *** RMSE was normalized to the standard deviation of ecec (0.54), while R$^2$ was adjusted to the initial number of covariates offered to the algorithm to build the model."
short_cap <- "ecec model performance when covariates are updated according to the expected added value."
latex(Scenario, file =  paste(update_covars_dir, "ecec-value-scenario.tex", sep = ""),
      label = "tab:ecec-value-scenario", table.env = TRUE, longtable = FALSE,
      na.blank = TRUE, ctable = TRUE,
      caption = long, caption.lot = short_cap, where = "ht",
      size = "scriptsize", insert.bottom = foot)