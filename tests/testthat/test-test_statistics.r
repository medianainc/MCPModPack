context("test_statistics")

library(DoseFinding)

n_cases = 1

test_that("Derivation of test statistics for trials with a normal endpoint", {

	# Select models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
	models = list(linear = NA, quadratic = -1, exponential = 1, emax = 0.2, logistic = c(0.1, 1), sigemax = c(0.1, 1))

	# One-sided Type I error rate
	alpha = 0.025

	# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
	direction = "increasing"

	# Model selection criterion
	model_selection = "AIC"

	endpoint_type = "Normal"	

	Delta = 0.5

	# List of test statistics derived by MCPModPack
	test_statistic1_list = list()

	# List of test statistics derived by DoseFinding
	test_statistic2_list = list()

	for (j in 1:n_cases) {

		# Total number of patients
		n = round(runif(1, 100, 200))

		# Number of doses
		n_doses = round(runif(1, 5, 8))

		# Generate doses and responses assuming a linear trend
		dose = rep(0, n)
		resp = rep(0, n)
		for (i in 1:n) {

			dose[i] = sample(0:(n_doses - 1))[1]
			resp[i] = rnorm(1, dose[i] / n_doses, 1)

		}

		# Run the MCPMod analysis
		results = MCPModAnalysis(endpoint_type = endpoint_type, 
			                models = models, 
			                dose = dose, 
			                resp = resp, 
			                alpha = alpha, 
			                direction = direction, 
			                model_selection = model_selection, 
			                Delta = Delta)

		dose_levels = 0:(n_doses - 1)

		data_set = data.frame(dose, resp)

		df_models = Mods(linear = NULL, quadratic = models$quadratic[1], exponential = models$exponential[1], emax = models$emax[1], logistic = models$logistic[1:2], sigEmax = models$sigemax[1:2], doses = dose_levels)

	    df_output = MCPMod(dose, resp, data_set, df_models, Delta = Delta)

        test_statistic1_list[[j]] = results$mcp_results$test_statistics

        test_statistic2_list[[j]] = as.numeric(df_output$MCTtest$tStat)


	}

    expect_equivalent(test_statistic1_list, test_statistic2_list, 0.001)

})

test_that("Derivation of optimal dose-response contrasts for trials with a binary endpoint", {

	# Select models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
	models = list(linear = NA, quadratic = -1, exponential = 1, emax = 0.2, logistic = c(0.1, 1), sigemax = c(0.1, 1))

	# One-sided Type I error rate
	alpha = 0.025

	# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
	direction = "increasing"

	# Model selection criterion
	model_selection = "AIC"

	endpoint_type = "Binary"	

	Delta = 0.2

	# List of test statistics derived by MCPModPack
	test_statistic1_list = list()

	# List of test statistics derived by DoseFinding
	test_statistic2_list = list()

	for (j in 1:n_cases) {

		# Total number of patients
		n = round(runif(1, 100, 200))

		# Number of doses
		n_doses = round(runif(1, 5, 8))

		# Generate doses and responses assuming a linear trend
		dose = rep(0, n)
		resp = rep(0, n)
		for (i in 1:n) {

			dose[i] = sample(0:(n_doses - 1))[1]
			resp[i] = rnorm(1, dose[i] / n_doses, 1) 
			resp[i] = as.numeric(resp[i] >= 0.5)			

		}

		# Run the MCPMod analysis
		results = MCPModAnalysis(endpoint_type = endpoint_type, 
			                models = models, 
			                dose = dose, 
			                resp = resp, 
			                alpha = alpha, 
			                direction = direction, 
			                model_selection = model_selection, 
			                Delta = Delta)

		dose_levels = 0:(n_doses - 1)

		data_set = data.frame(dose, resp)

		df_models = Mods(linear = NULL, quadratic = models$quadratic[1], exponential = models$exponential[1], emax = models$emax[1], logistic = models$logistic[1:2], sigEmax = models$sigemax[1:2], doses = dose_levels)

		n_patients = rep(0, length(dose_levels))
		mu_hat = rep(0, length(dose_levels))
		diag_mat = rep(0, length(dose_levels))

		resp_rate = rep(0, length(dose_levels))

		for (i in 1:length(dose_levels)) {
			n_patients[i] = sum(dose == dose_levels[i])
			resp_rate[i] = sum(resp[dose == dose_levels[i]]) / n_patients[i]
			diag_mat[i] = 1 / (n_patients[i] * resp_rate[i] * (1 - resp_rate[i]))
			mu_hat[i] = log(resp_rate[i] / (1 - resp_rate[i]))
		}

		S = diag(diag_mat)
		df_output = MCPMod(dose_levels, mu_hat, S = S, models = df_models, type = "general", Delta = Delta)

        test_statistic1_list[[j]] = results$mcp_results$test_statistics

        test_statistic2_list[[j]] = as.numeric(df_output$MCTtest$tStat)

	}

    expect_equivalent(test_statistic1_list, test_statistic2_list, 0.001)

})

test_that("Derivation of optimal dose-response contrasts for trials with a count endpoint", {

	# Select models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
	models = list(linear = NA, quadratic = -1, exponential = 1, emax = 0.2, logistic = c(0.1, 1), sigemax = c(0.1, 1))

	# One-sided Type I error rate
	alpha = 0.025

	# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
	direction = "increasing"

	# Model selection criterion
	model_selection = "AIC"

	endpoint_type = "Count"	

	Delta = 2

	# List of test statistics derived by MCPModPack
	test_statistic1_list = list()

	# List of test statistics derived by DoseFinding
	test_statistic2_list = list()

	for (j in 1:n_cases) {

		# Total number of patients
		n = round(runif(1, 100, 200))

		# Number of doses
		n_doses = round(runif(1, 5, 8))

		# Generate doses and responses assuming a linear trend
		dose = rep(0, n)
		resp = rep(0, n)
		for (i in 1:n) {

			dose[i] = sample(0:(n_doses - 1))[1]
			resp[i] = ceiling(abs(rnorm(1, 2 + 5 * dose[i] / n_doses, 1))) 		

		}

		# Vector of over dispersion parameters for count endpoints
		theta = rep(2, n_doses)

		# Run the MCPMod analysis
		results = MCPModAnalysis(endpoint_type = endpoint_type, 
			                models = models, 
			                dose = dose, 
			                resp = resp, 
			                alpha = alpha, 
			                direction = direction, 
			                model_selection = model_selection, 
			                Delta = Delta,
			                theta = theta)

		dose_levels = 0:(n_doses - 1)

		data_set = data.frame(dose, resp)

		df_models = Mods(linear = NULL, quadratic = models$quadratic[1], exponential = models$exponential[1], emax = models$emax[1], logistic = models$logistic[1:2], sigEmax = models$sigemax[1:2], doses = dose_levels)

		n_patients = rep(0, length(dose_levels))
		mu_hat = rep(0, length(dose_levels))
		diag_mat = rep(0, length(dose_levels))
		ave_count = rep(0, length(dose_levels))

		for (i in 1:length(dose_levels)) {
			n_patients[i] = sum(dose == dose_levels[i])
			ave_count[i] = sum(resp[dose == dose_levels[i]]) / n_patients[i]
			diag_mat[i] = (theta[i] + ave_count[i]) / (n_patients[i] * theta[i] * ave_count[i])
			mu_hat[i] = log(ave_count[i])
		}

		S = diag(diag_mat)
		df_output = MCPMod(dose_levels, mu_hat, S = S, models = df_models, type = "general", Delta = Delta)

        test_statistic1_list[[j]] = results$mcp_results$test_statistics

        test_statistic2_list[[j]] = as.numeric(df_output$MCTtest$tStat)


	}

    expect_equivalent(test_statistic1_list, test_statistic2_list, 0.001)

})
