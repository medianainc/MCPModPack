context("test-statistics")

library(DoseFinding)

test_that("Derivation of test statistics for trials with a normal endpoint", {

	# Select models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
	models = list(linear = NA, 
		          quadratic = c(-1), 		
		          exponential = c(1), 
		          emax = c(0.2), 
		          logistic = c(0.1, 1),
		          sigemax = c(0.1, 1)) 

	# One-sided Type I error rate
	alpha = 0.025

	# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
	direction = "increasing"

	# Model selection criterion
	model_selection = "AIC"

	endpoint_type = "Normal"	

	Delta = 0.5

    # The dose and resp variables are extracted from the 
	# built-in data set (normal) as shown below

	dose = normal$dose
	resp = normal$resp 

	# Run the MCPMod analysis
	results = MCPModAnalysis(endpoint_type = endpoint_type, 
				                models = models, 
				                dose = dose, 
				                resp = resp, 
				                alpha = alpha, 
				                direction = direction, 
				                model_selection = model_selection, 
				                Delta = Delta)

	data_set = data.frame(dose, resp)

    doses = sort(unique(dose))
	n_doses = length(doses)

	df_models = Mods(linear = NULL, 
		             quadratic = models$quadratic[1], 
		             exponential = models$exponential[1], 
		             emax = models$emax[1], 
		             logistic = models$logistic[1:2], 
		             sigEmax = models$sigemax[1:2], 
		             doses = doses)

    df_output = MCPMod(dose, resp, data_set, df_models, Delta = Delta)

	test_statistic1 = results$mcp_results$test_statistics

    test_statistic2 = as.numeric(df_output$MCTtest$tStat)

    expect_equivalent(test_statistic1, test_statistic2, 0.001)

})

test_that("Derivation of test statistics for trials with a binary endpoint", {

	# Select models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
	models = list(linear = NA, 
		          quadratic = c(-1), 		
		          exponential = c(1), 
		          emax = c(0.2), 
		          logistic = c(0.1, 1),
		          sigemax = c(0.1, 1)) 

	# One-sided Type I error rate
	alpha = 0.025

	# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
	direction = "increasing"

	# Model selection criterion
	model_selection = "AIC"

	endpoint_type = "Binary"	

	Delta = 0.2

    # The dose and resp variables are extracted from the 
	# built-in data set (binary) as shown below

	dose = binary$dose
	resp = binary$resp 

	# Run the MCPMod analysis
	results = MCPModAnalysis(endpoint_type = endpoint_type, 
				                models = models, 
				                dose = dose, 
				                resp = resp, 
				                alpha = alpha, 
				                direction = direction, 
				                model_selection = model_selection, 
				                Delta = Delta)

	data_set = data.frame(dose, resp)

    doses = sort(unique(dose))
	n_doses = length(doses)

	df_models = Mods(linear = NULL, 
		             quadratic = models$quadratic[1], 
		             exponential = models$exponential[1], 
		             emax = models$emax[1], 
		             logistic = models$logistic[1:2], 
		             sigEmax = models$sigemax[1:2], 
		             doses = doses)

	n_patients = rep(0, length(doses))
	mu_hat = rep(0, length(doses))
	diag_mat = rep(0, length(doses))

	resp_rate = rep(0, length(doses))

	for (i in 1:length(doses)) {
		n_patients[i] = sum(dose == doses[i])
		resp_rate[i] = sum(resp[dose == doses[i]]) / n_patients[i]
		diag_mat[i] = 1 / (n_patients[i] * resp_rate[i] * (1 - resp_rate[i]))
		mu_hat[i] = log(resp_rate[i] / (1 - resp_rate[i]))
	}

	S = diag(diag_mat)
	df_output = MCPMod(doses, mu_hat, S = S, models = df_models, type = "general", Delta = Delta)

    test_statistic1 = results$mcp_results$test_statistics

    test_statistic2 = as.numeric(df_output$MCTtest$tStat)

    expect_equivalent(test_statistic1, test_statistic2, 0.001)


})

test_that("Derivation of test statistics for trials with a count endpoint", {

	# Select models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
	models = list(linear = NA, 
		          quadratic = c(-1), 		
		          exponential = c(1), 
		          emax = c(0.2), 
		          logistic = c(0.1, 1),
		          sigemax = c(0.1, 1)) 

	# One-sided Type I error rate
	alpha = 0.025

	# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
	direction = "increasing"

	# Model selection criterion
	model_selection = "AIC"

	endpoint_type = "Count"	

	Delta = 0.2

    # The dose and resp variables are extracted from the 
	# built-in data set (count) as shown below

	dose = count$dose
	resp = count$resp 

	# Vector of over dispersion parameters for count endpoints
    doses = sort(unique(dose))
	n_doses = length(doses)
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

	data_set = data.frame(dose, resp)

	df_models = Mods(linear = NULL, 
		             quadratic = models$quadratic[1], 
		             exponential = models$exponential[1], 
		             emax = models$emax[1], 
		             logistic = models$logistic[1:2], 
		             sigEmax = models$sigemax[1:2], 
		             doses = doses)

	n_patients = rep(0, length(doses))
	mean = rep(0, length(doses))
	mu_hat = rep(0, length(doses))
	diag_mat = rep(0, length(doses))

	resp_rate = rep(0, length(doses))

	for (i in 1:length(doses)) {
		n_patients[i] = sum(dose == doses[i])
		mean[i] = mean(resp[dose == doses[i]]) 
		diag_mat[i] = (theta[i] + mean[i]) / (n_patients[i] * theta[i] * mean[i])
		mu_hat[i] = log(mean[i])
	}

	S = diag(diag_mat)
	df_output = MCPMod(doses, mu_hat, S = S, models = df_models, type = "general", Delta = Delta)

    test_statistic1 = results$mcp_results$test_statistics

    test_statistic2 = as.numeric(df_output$MCTtest$tStat)

    expect_equivalent(test_statistic1, test_statistic2, 0.001)

})