options("scipen" = 100, "digits" = 4, warn = -1)

require(devEMF)
require(officer)
require(flextable)
require(mvtnorm)

DF_endpoint_list = c("Normal", "Binary", "Count")
DF_model_list = c("Linear", "Quadratic", "Exponential", "Emax", "Logistic", "SigEmax")
DF_model_list_short = c("Lin", "Quad", "Exp", "Emax", "Logist", "SigEmax")
DF_model_parameters = list(c("e0", "delta"), c("e0", "delta1", "delta2"), c("e0", "e1", "delta"), c("e0", "eMax", "ed50"), c("e0", "eMax", "ed50", "delta"), c("e0", "eMax", "ed50", "h"))


n_evaluation_points = 100

Logit = function(x) {
    return(log(x /(1 - x)))
}

AntiLogit = function(x) {
    return(1 / (1 + exp(-x)))
}

tocap = function(x) {
  s = strsplit(x, " ")[[1]]
  paste0(toupper(substring(s, 1,1)), substring(s, 2))
}

FormatMatrix = function(mat, format) {
    nrows = dim(mat)[1]
    ncols = dim(mat)[2] 

    for_mat = matrix(0, nrows, ncols)

    for (i in 1:nrows) {
      for (j in 1:ncols) {
        for_mat[i, j] = sprintf(format, mat[i, j]) 
      }
    }

    return(for_mat)

} 

# Check if the number is an integer
is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# Arguments:
# parameter: Parameter's value
# n_values: Required number of values
# lower_values: Lower range
# lower_values_sign: Inequality for evaluating the lower range
# upper_values: Upper range
# upper_values_sign: Inequality for evaluating the upper range
# parameter_name: Parameter's name
# component_name: Names of the individual components
# type: Parameter's type (double or integer)
# default_value: Default value
ContinuousErrorCheck = function(parameter, n_values, lower_values, lower_values_sign, upper_values, upper_values_sign, parameter_name, component_name, type = "double", default_value = NA) {

    if (is.null(parameter)) {

        if (!is.na(default_value)) {
            for (i in 1:n_values) {
                parameter[i] = default_value
            }
            return(parameter)
        } else {
            error_message = paste0(parameter_name, " must be specified.") 
            stop(error_message, call. = FALSE)
        }
    } 

    if (!is.na(n_values)) {

      if (length(parameter) != n_values) {
          error_message = paste0(parameter_name, ": ", n_values, " values must be specified.") 
          stop(error_message, call. = FALSE)
      } 

    } else {

      n_values = length(parameter)

    }
    
    for (i in 1:n_values) {

        if (type == "double") {

            if (!is.numeric(parameter[i])) {
                error_message = paste0(parameter_name, ": ", component_name[i], " must be numeric.") 
                stop(error_message, call. = FALSE)            
            }
        }

        if (type == "integer") {

            if (!is.wholenumber(parameter[i])) {
                error_message = paste0(parameter_name, ": ", component_name[i], " must be an integer.") 
                stop(error_message, call. = FALSE)            
            }
        }

        if (length(lower_values) == 1) {

          if (!is.na(lower_values)) {            
              if (lower_values_sign == ">" & parameter[i] <= lower_values) {
                  error_message = paste0(parameter_name, ": Each value must be > ", lower_values, ".") 
                  stop(error_message, call. = FALSE)
              }
              if (lower_values_sign == ">=" & parameter[i] < lower_values) {
                  error_message = paste0(parameter_name, ": Each value must be >= ", lower_values, ".") 
                  stop(error_message, call. = FALSE)
              }
          }


        } else {

          if (!is.na(lower_values[i])) {            
              if (lower_values_sign[i] == ">" & parameter[i] <= lower_values[i]) {
                  error_message = paste0(parameter_name, ": ", component_name[i], " must be > ", lower_values[i], ".") 
                  stop(error_message, call. = FALSE)
              }
              if (lower_values_sign[i] == ">=" & parameter[i] < lower_values[i]) {
                  error_message = paste0(parameter_name, ": ", component_name[i], " must be >= ", lower_values[i], ".") 
                  stop(error_message, call. = FALSE)
              }
          }

        }

        if (length(upper_values) == 1) {

          if (!is.na(upper_values)) {            
              if (upper_values_sign == "<" & parameter[i] >= upper_values) {
                  error_message = paste0(parameter_name, ": Each value must be < ", upper_values, ".") 
                  stop(error_message, call. = FALSE)
              }
              if (upper_values_sign == "<=" & parameter[i] > upper_values) {
                  error_message = paste0(parameter_name, ": Each value must be <= ", upper_values, ".") 
                  stop(error_message, call. = FALSE)
              }
          }

        } else {

          if (!is.na(upper_values[i])) {
              if (upper_values_sign[i] == "<" & parameter[i] >= upper_values[i]) {
                  error_message = paste0(parameter_name, ": ", component_name[i], " must be < ", upper_values[i], ".") 
                  stop(error_message, call. = FALSE)
              }
              if (upper_values_sign[i] == "<=" & parameter[i] > upper_values[i]) {
                  error_message = paste0(parameter_name, ": ", component_name[i], " must be <= ", upper_values[i], ".") 
                  stop(error_message, call. = FALSE)
              }
            }
        }
    }

    return(parameter)
}

# Dose-response functions
DRFunction = function(model_index, coef, x) {

    # Linear model
    if (model_index == 1) {
        y = coef[1] + coef[2] * x
    }

    # Quadratic model
    if (model_index == 2) {
        y = coef[1] + coef[2] * x + coef[3] * x^2
    }

    # Exponential model
    if (model_index == 3) {
        y = coef[1] + coef[2] * (exp(x / coef[3]) - 1)
    }

    # Emax model
    if (model_index == 4) {
        y = coef[1] + coef[2] * x / (coef[3] + x)
    }

    # Logistic model
    if (model_index == 5) {
            den = 1.0 + exp((coef[3] - x) / coef[4])
            y = coef[1] + coef[2] / den
    }

    # SigEmax model
    if (model_index == 6) {
        den = x^coef[4] + coef[3]^coef[4]
        y = coef[1] + coef[2] * x^coef[4] / den
    }

    return(y)

}

# Compute the model parameters to match the placebo and maximum effects
ComputeDRFunctionParameters = function(model_index, placebo_effect, max_effect, max_dose, parameters) {

    # Linear model
    if (model_index == 1) {
        coef = rep(0, 2)
        coef[1] = placebo_effect
        coef[2] = max_effect / max_dose
    }

    # Quadratic model (maximum is assumed to be achieved at the mid-point of the dose range)
    if (model_index == 2) {
        coef = rep(0, 3)
        coef[1] = placebo_effect
        coef[2] = 4 * max_effect / max_dose
        coef[3] = - coef[2] / max_dose 
    }

    # Exponential model
    if (model_index == 3) {
        coef = rep(0, 3)
        coef[1] = placebo_effect
        coef[2] = max_effect / (exp(max_dose / parameters[1]) - 1)
        coef[3] = parameters[1]
    }

    # Emax model
    if (model_index == 4) {
        coef = rep(0, 3)
        coef[1] = placebo_effect
        coef[2] = max_effect * (parameters[1] + max_dose) / max_dose
        coef[3] = parameters[1]
    }

    # Logistic model
    if (model_index == 5) {
        coef = rep(0, 4)
        temp_coef = c(0, 1, parameters[1], parameters[2])
        temp =  max_effect / (DRFunction(5, temp_coef, max_dose) - DRFunction(5, temp_coef, 0))
        coef[1] = placebo_effect- temp * DRFunction(5, temp_coef, 0)
        coef[2] = temp
        coef[3] = parameters[1]
        coef[4] = parameters[2]
    }

    # SigEmax model
    if (model_index == 6) {
        coef = rep(0, 4)
        coef[1] = placebo_effect
        coef[2] = max_effect * (parameters[1]^parameters[2] + max_dose^parameters[2]) / max_dose^parameters[2]
        coef[3] = parameters[1]
        coef[4] = parameters[2]
    }

    return(coef)

}

# Evaluate a dose-response function over the range of doses
EvaluateDRFunction = function(model_index, endpoint_index, coef, dose) {

    # Evaluate the dose-response function for the given endpoint type

    x = seq(from = min(dose), to = max(dose), length.out = n_evaluation_points)
    y = rep(0, length(x))

    for (i in 1:length(x)) {

        # Normal endpoint
        if (endpoint_index == 1) y[i] = DRFunction(model_index, coef, x[i])

        # Binary endpoint
        if (endpoint_index == 2) y[i] = AntiLogit(DRFunction(model_index, coef, x[i]))

        # Count endpoint
        if (endpoint_index == 3) y[i] = exp(DRFunction(model_index, coef, x[i]))

    }

    return(list(x = x , y = y))

}


# Fit dose-response models
ModStep = function(endpoint_index, selected_models, theta_vector, dose, resp, delta, direction_index) {

    # Maximum number of iterations to find maximum likelihood estimates
    maxit = 300

    withCallingHandlers({          

        model_fit = MCPModFitDRModels(endpoint_index, selected_models, dose, resp, delta, direction_index, maxit, theta_vector)

        },

        warning = function(c) {
          msg <- conditionMessage(c)
          if ( grepl("the line search step became smaller than the minimum value allowed", msg, fixed = TRUE) ) {
            invokeRestart("muffleWarning")
          }
        }

    )

    results = list(model_fit = model_fit,
                   dose = dose, 
                   resp = resp)

    return(results)

}

# Compute the optimal contrasts, contrast correlation matrix and adjusted critical value
ContrastStep = function(endpoint_index, selected_models, user_specified, n_groups, dose_levels, alpha, direction_index, mean_group, theta) { 

    #####################################################

    # Total number of models
    n_models = length(DF_model_list)

    # List of selected models
    model_list = (1:n_models)[selected_models]
    n_selected_models = length(model_list)

    doses = dose_levels

    n_doses = length(doses)
    n_patients = sum(n_groups)

    max_dose = max(doses)
    diag_vec = rep(0, n_doses)

    corr_matrix = 0

    # Normal endpoint
    if (endpoint_index == 1) {

      for (i in 1:n_doses) diag_vec[i] = n_groups[i] 

    }

    # Binary endpoint  
    if (endpoint_index == 2) {

      for (i in 1:n_doses) diag_vec[i] = n_groups[i] * mean_group[i] * (1 - mean_group[i])

    }

    # Count endpoint  
    if (endpoint_index == 3) {

      for (i in 1:n_doses) diag_vec[i] = n_groups[i] * theta[i] * mean_group[i] / (theta[i] + mean_group[i])

    }

    S =  diag(1 / diag_vec)  
    Sinv = diag(diag_vec)

    dr_model = rep(0, n_doses)

    #####################################################

    # Compute the model-specific optimal contrasts

    opt_contrast = matrix(0, n_doses, n_models)

    # Set up dose-response models based on the initial values
    for (i in 1:n_models) {

        # Define the vector of starting values
        n_parameters = length(DF_model_parameters[[i]]) 

        non_linear_parameters = user_specified[[i]][1:n_parameters]

        if (length(non_linear_parameters) >= 3) non_linear_parameters = non_linear_parameters[3:length(non_linear_parameters)] else non_linear_parameters = 0

        # Parameters of a standardized model
        if (i != 2) parameter_values = ComputeDRFunctionParameters(i, 0, 1, max_dose, non_linear_parameters)

        # Alternative standardization for the quadratic model
        if (i == 2) {

          parameter_values = rep(0, 3)

          if (direction_index == 1) {
  
            temp = -0.5 / non_linear_parameters[1]
            parameter_values[1] = 0
            parameter_values[2] = 1 / (temp + non_linear_parameters[1] * temp^2)
            parameter_values[3] = non_linear_parameters[1] * parameter_values[2] 

          }

          if (direction_index == -1) {
  
            temp = 0.5 / non_linear_parameters[1]
            parameter_values[1] = 0
            parameter_values[2] = -1 / (temp - non_linear_parameters[1] * temp^2)
            parameter_values[3] = -non_linear_parameters[1] * parameter_values[2] 

          }

        }

        for (j in 1:n_doses) {

          dr_model[j] = DRFunction(i, parameter_values, doses[j]) 
          if (i == 2 & direction_index == -1) dr_model[j] = -DRFunction(i, parameter_values, doses[j]) 

        }

        # Optimal contrasts  
        dr_expected = sum(dr_model * rowSums(Sinv))/sum(rowSums(Sinv)) 
        contrast = Sinv %*% (dr_model - dr_expected)
        contrast = contrast - sum(contrast)
        opt_contrast[, i] = contrast / sqrt(sum(contrast^2))

    }

    #####################################################

    # Apply the list of selected models
    opt_contrast = opt_contrast[, model_list] 

    #####################################################

    if (n_selected_models >= 2) {

      # Compute the correlation matrix for the model-specific test statistics

      cov_mat = t(opt_contrast) %*% S %*% opt_contrast
      diag_mat = diag(sqrt(diag(cov_mat)))
      corr_matrix = solve(diag_mat) %*% cov_mat %*% solve(diag_mat)

    }

    #####################################################

    # Critical value based on a univariate or multivariate t distribution
    if (n_selected_models >= 2) crit_value = qmvt(p = 1 - alpha, tail = "lower.tail", df = n_patients - n_doses, corr = corr_matrix, maxpts = 30000, abseps = 0.001, releps = 0, algorithm = GenzBretz())$quantile else crit_value = qt(p = 1 - alpha, df = n_patients - n_doses)

    # Account for the direction of the dose-response relationship  
    crit_value = crit_value * direction_index

    results = list(opt_contrast = as.matrix(opt_contrast),
                   corr_matrix = corr_matrix,
                   crit_value = crit_value)

    return(results)


}
# End of ContrastStep

# Compute the test statistics
MCPStep = function(endpoint_index, contrast_results, selected_models, user_specified, dose, resp, alpha, direction_index) {

    #####################################################

    # Total number of models
    n_models = length(DF_model_list)

    # List of selected models
    model_list = (1:n_models)[selected_models]
    n_selected_models = length(model_list)

    dose_levels = sort(unique(dose))
    n_doses = length(dose_levels)
    n_groups = table(dose)
    n_patients = length(resp)

    test_statistics = rep(0, n_selected_models)

    theta = user_specified$theta

    corr_matrix = as.matrix(contrast_results$corr_matrix)
    opt_contrast = as.matrix(contrast_results$opt_contrast)

    #####################################################

    # Compute the model-specific test statistics

    # Normal endpoints
    if (endpoint_index == 1) {

        # Group-specific means
        mean_group = rep(0, n_doses)

        for (j in 1:n_doses) {
          for (i in 1:n_patients) {
            if (dose[i] == dose_levels[j]) {
              mean_group[j] = mean_group[j] + resp[i]
            }
          }
          mean_group[j] = mean_group[j] / n_groups[j]
        }

        # Pooled variance estimate
        pooled_variance = 0

        for (j in 1:n_doses) {
          for (i in 1:n_patients) {
            if (dose[i] == dose_levels[j]) {
              pooled_variance = pooled_variance + (resp[i] - mean_group[j])^2      
            }
          }
        }

        pooled_variance = pooled_variance / (n_patients - n_doses)

        # Model-specific test statistics
        for (i in 1:n_selected_models) { 

          num = 0
          den = 0

          for (j in 1:n_doses) { 

            num = num + opt_contrast[j, i] * mean_group[j]
            den = den + pooled_variance * opt_contrast[j, i]^2 / n_groups[j]

          }

          test_statistics[i] = num / sqrt(den)

        }

    }

    # Binary endpoints
    if (endpoint_index == 2) {

        # Group-specific means, logits and variances
        mean_group = rep(0, n_doses)
        logit_group = rep(0, n_doses)
        variance_group = rep(0, n_doses)

        for (j in 1:n_doses) {
          for (i in 1:n_patients) {
            if (dose[i] == dose_levels[j]) {
              mean_group[j] = mean_group[j] + resp[i]
            }
          }
          if (mean_group[j] == 0) mean_group[j] = 1 / (3 * n_groups[j] + 2)
          if (mean_group[j] == n_groups[j]) mean_group[j] = (3 * n_groups[j] + 1) / (3 * n_groups[j] + 2)
          if (mean_group[j] > 0 & mean_group[j] < n_groups[j]) mean_group[j] = mean_group[j] / n_groups[j]
          logit_group[j] = log(mean_group[j] / (1 - mean_group[j]))
          variance_group[j] = 1 / (n_groups[j] * mean_group[j] * (1 - mean_group[j]))
        }

        # Model-specific test statistics
        for (i in 1:n_selected_models) { 

          num = 0
          den = 0

          for (j in 1:n_doses) { 

            num = num + opt_contrast[j, i] * logit_group[j]
            den = den + variance_group[j] * opt_contrast[j, i]^2

          }

          test_statistics[i] = num / sqrt(den)

        }

    }

    # Count endpoints
    if (endpoint_index == 3) {

        # Group-specific means, logits and variances
        mean_group = rep(0, n_doses)
        logit_group = rep(0, n_doses)
        variance_group = rep(0, n_doses)

        for (j in 1:n_doses) {
          for (i in 1:n_patients) {
            if (dose[i] == dose_levels[j]) {
              mean_group[j] = mean_group[j] + resp[i]
            }
          }
          if (mean_group[j] == 0) mean_group[j] = qgamma(0.5, shape = 1 / 3, scale = 1 / n_groups[j])
          if (mean_group[j] > 0) mean_group[j] = mean_group[j] / n_groups[j]
          variance_group[j] = (theta[j] + mean_group[j]) / (n_groups[j] * theta[j] * mean_group[j])
        }

        # Model-specific test statistics
        for (i in 1:n_selected_models) { 

          num = 0
          den = 0

          for (j in 1:n_doses) { 

            num = num + opt_contrast[j, i] * log(mean_group[j])
            den = den + variance_group[j] * opt_contrast[j, i]^2

          }

          test_statistics[i] = num / sqrt(den)

        }

    }

    # Apply the list of selected models
    # test_statistics = test_statistics[model_list]

    adj_pvalues = rep(0, n_selected_models)

    if (n_selected_models >= 2) {

      # Compute adjusted p-values from a multivariate t distribution
      for (i in 1:n_selected_models) {

        if (direction_index == 1) {

            adj_pvalues[i] = 1 - pmvt(lower = rep(-Inf, n_selected_models), upper = rep(test_statistics[i], n_selected_models), df = n_patients - n_doses, corr = corr_matrix, maxpts = 30000, abseps = 0.001, releps = 0)

        } 

        if (direction_index == -1) {

            adj_pvalues[i] = 1 - pmvt(lower = rep(test_statistics[i], n_selected_models), upper = rep(Inf, n_selected_models), df = n_patients - n_doses, corr = corr_matrix, maxpts = 30000, abseps = 0.001, releps = 0)

        }

      }

    } else {

      # Compute the adjusted p-value from a univariate t distribution
      if (direction_index == 1) adj_pvalues[1] = 1 - pt(test_statistics[1], df = n_patients - n_doses)
      if (direction_index == -1) adj_pvalues[1] = pt(test_statistics[1], df = n_patients - n_doses)

    }

    # Identify significant models
    sign_model = as.numeric(adj_pvalues <= alpha)

    #####################################################


    results = list(test_statistics = test_statistics,
                   adj_pvalues = adj_pvalues,
                   sign_model = sign_model)

    return(results)

}

MCPModSimulation = function(endpoint_type, models, alpha = 0.025, direction = "increasing", model_selection = "AIC", Delta, theta = 0, sim_models, sim_parameters) {

    if (missing(endpoint_type)) stop("Endpoint type (endpoint_type): Value must be specified.", call. = FALSE)      
    if (missing(models)) stop("Candidate dose-response models (models): Value must be specified.", call. = FALSE)

    if (missing(Delta)) stop("Treatment effect for identifying the target dose (Delta): Value must be specified.", call. = FALSE)

    if (missing(sim_models)) stop("Simulation models (sim_models): Value must be specified.", call. = FALSE)

    if (missing(sim_parameters)) stop("Simulation parameters (sim_parameters): Value must be specified.", call. = FALSE)

    # Total number of models
    n_models = length(DF_model_list)

    # Error checks 

    if (!tolower(endpoint_type) %in% tolower(DF_endpoint_list)) stop("MCPModSimulation: Endpoint type (endpoint_type): Value must be Normal, Binary or Count.", call. = FALSE)

    endpoint_index = 1  

    for (i in 1:length(DF_endpoint_list)) {
        if (tolower(DF_endpoint_list[i]) == tolower(endpoint_type)) endpoint_index = i
    }   

    if (!tolower(direction) %in% c("increasing", "decreasing")) stop("MCPModSimulation: Direction of the dose-response relationship (direction): Value must be Increasing or Decreasing.", call. = FALSE)

    if (tolower(direction) == "decreasing") direction_index = -1
    if (tolower(direction) == "increasing") direction_index = 1

    if (length(models) < 1) stop("MCPModSimulation: List of dose-response models and initial parameter values (models): At least one model must be specified.", call. = FALSE)

    selected_models = rep(FALSE, n_models)
    user_specified = list()

    if (endpoint_index == 1) user_specified$linear = c(0, 0, 1) else user_specified$linear = c(0, 0)
    if (!is.null(models$linear)) selected_models[1] = TRUE

    if (endpoint_index == 1) user_specified$quadratic = c(0, 0, 0, 1) else user_specified$quadratic = c(0, 0, 0)
    if (!is.null(models$quadratic)) {
      selected_models[2] = TRUE 

      if (direction_index == 1) {

        user_specified$quadratic[3] =  ContinuousErrorCheck(models$quadratic[1], 
                                                             1, 
                                                             lower_values = NA,
                                                             lower_values_sign = NA,
                                                             upper_values = 0,
                                                             upper_values_sign = "<",
                                                             "Quadratic model (quadratic)",
                                                             c("delta2"),
                                                             "double",
                                                             NA) 

      }

      if (direction_index == -1) {

        user_specified$quadratic[3] =  ContinuousErrorCheck(models$quadratic[1], 
                                                             1, 
                                                             lower_values = 0,
                                                             lower_values_sign = ">",
                                                             upper_values = NA,
                                                             upper_values_sign = NA,
                                                             "Quadratic model (quadratic)",
                                                             c("delta2"),
                                                             "double",
                                                             NA) 

      }    
    }

    if (endpoint_index == 1) user_specified$exponential = c(0, 0, 0, 1) else user_specified$exponential = c(0, 0, 0)
    if (!is.null(models$exponential)) {
      selected_models[3] = TRUE 
      user_specified$exponential[3] =  ContinuousErrorCheck(models$exponential[1], 
                                                           1, 
                                                           lower_values = c(0),
                                                           lower_values_sign = c(">"),
                                                           upper_values = c(NA),
                                                           upper_values_sign = c(NA),
                                                           "Exponential model (exponential)",
                                                           c("delta"),
                                                           "double",
                                                           NA) 
    } 

    if (endpoint_index == 1) user_specified$emax = c(0, 0, 0, 1) else user_specified$emax = c(0, 0, 0)
    if (!is.null(models$emax)) {
      selected_models[4] = TRUE 
      user_specified$emax[3] =  ContinuousErrorCheck(models$emax[1], 
                                                           1, 
                                                           lower_values = c(0),
                                                           lower_values_sign = c(">"),
                                                           upper_values = c(NA),
                                                           upper_values_sign = c(NA),
                                                           "Emax model (emax)",
                                                           c("ED50"),
                                                           "double",
                                                           NA) 

    } 

    if (endpoint_index == 1) user_specified$logistic = c(0, 0, 0, 0, 1) else user_specified$logistic = c(0, 0, 0, 0)
    if (!is.null(models$logistic)) {
      selected_models[5] = TRUE 
      user_specified$logistic[3:4] =  ContinuousErrorCheck(models$logistic[1:2], 
                                                           2, 
                                                           lower_values = c(0, 0),
                                                           lower_values_sign = c(">", ">"),
                                                           upper_values = c(NA, NA),
                                                           upper_values_sign = c(NA, NA),
                                                           "Logistic model (logistic)",
                                                           c("ED50", "delta"),
                                                           "double",
                                                           NA) 

    } 

    if (endpoint_index == 1) user_specified$sigemax = c(0, 0, 0, 0, 1) else user_specified$sigemax = c(0, 0, 0, 0)
    if (!is.null(models$sigemax)) {
      selected_models[6] = TRUE 
      user_specified$sigemax[3:4] =  ContinuousErrorCheck(models$sigemax[1:2], 
                                                           2, 
                                                           lower_values = c(0, 0),
                                                           lower_values_sign = c(">", ">"),
                                                           upper_values = c(NA, NA),
                                                           upper_values_sign = c(NA, NA),
                                                           "SigEmax model (sigemax)",
                                                           c("ED50", "h"),
                                                           "double",
                                                           NA) 
    } 

    n_selected_models = sum(selected_models)

    if (n_selected_models == 0) stop("MCPModSimulation: List of models and initial parameter values (models): At least one model must be specified (linear, quadratic, exponential, emax, logistic or sigemax).", call. = FALSE)

    alpha = 
          ContinuousErrorCheck(alpha, 
                               1, 
                               lower_values = c(0.001),
                               lower_values_sign = c(">"),
                               upper_values = c(0.999),
                               upper_values_sign = c("<"),
                               "One-sided Type I error rate (alpha)",
                               c("Value"),
                               "double",
                               NA) 

    if (!model_selection %in% c("AIC", "maxT", "aveAIC")) stop("MCPModSimulation: Model selection criterion (model_selection): Value must be AIC, maxT or aveAIC.", call. = FALSE)

    if (model_selection == "AIC") model_selection_index = 1  
    if (model_selection == "maxT") model_selection_index = 2  
    if (model_selection == "aveAIC") model_selection_index = 3  

    delta = 
      ContinuousErrorCheck(Delta, 
                           1, 
                           lower_values = c(NA),
                           lower_values_sign = c(NA),
                           upper_values = c(NA),
                           upper_values_sign = c(NA),
                           "Treatment effect for identifying the target dose (Delta)",
                           c("Value"),
                           "double",
                           NA) 

    if (direction_index == 1 & delta <= 0) stop("MCPModSimulation: Treatment effect for identifying the target dose (Delta): Value must be positive if the direction of the dose-response relationship (direction) is Increasing.", call. = FALSE)

    if (direction_index == -1 & delta >= 0) stop("MCPModSimulation: Treatment effect for identifying the target dose (Delta): Value must be negative if the direction of the dose-response relationship (direction) is Decreasing.", call. = FALSE)

      n = ContinuousErrorCheck(sim_parameters$n, 
                         NA, 
                         lower_values = 0,
                         lower_values_sign = ">",
                         upper_values = NA,
                         upper_values_sign = NA,
                         "Number of patients in the simulation model (n)",
                         NA,
                         "double",
                         NA) 

      dose_levels = ContinuousErrorCheck(sim_parameters$doses, 
                                         NA, 
                                         lower_values = 0,
                                         lower_values_sign = ">=",
                                         upper_values = 1000,
                                         upper_values_sign = "<=",
                                         "Dose levels in the simulation model (doses)",
                                         NA,
                                         "double",
                                         NA) 

    if(length(dose_levels) != length(n)) stop("MCPModSimulation: The length of the dose vector (doses) must be equal to the length of the sample size vector (n) in the simulation model.", call. = FALSE)    

    n_doses = length(dose_levels)

    if (!is.null(sim_parameters$dropout_rate)) {

      dropout_rate = 
        ContinuousErrorCheck(sim_parameters$dropout_rate, 
                             1, 
                             lower_values = c(0),
                             lower_values_sign = c(">="),
                             upper_values = c(1),
                             upper_values_sign = c("<"),
                             "Patient dropout rate in the simulation model (dropout_rate)",
                             c("Value"),
                             "double",
                             NA) 

    } else {
      dropout_rate = 0
    }

    go_threshold = 
      ContinuousErrorCheck(sim_parameters$go_threshold, 
                           1, 
                           lower_values = c(NA),
                           lower_values_sign = c(NA),
                           upper_values = c(NA),
                           upper_values_sign = c(NA),
                           "Threshold for computing go probabilities (go_threshold)",
                           c("Value"),
                           "double",
                           NA) 

    if (direction_index == 1 & go_threshold <= 0) stop("MCPModSimulation: Threshold for computing go probabilities (go_threshold): Value must be positive if the direction of the dose-response relationship (direction) is Increasing.", call. = FALSE)

    if (direction_index == -1 & go_threshold >= 0) stop("MCPModSimulation: Threshold for computing go probabilities (go_threshold): Value must be negative if the direction of the dose-response relationship (direction) is Decreasing.", call. = FALSE)

    if (!is.null(sim_parameters$nsims)) {

      nsims = 
        ContinuousErrorCheck(sim_parameters$nsims, 
                             1, 
                             lower_values = c(1),
                             lower_values_sign = c(">="),
                             upper_values = c(10000),
                             upper_values_sign = c("<="),
                             "Number of simulations (nsims)",
                             c("Value"),
                             "int",
                             NA) 

    } else {
      nsims = 1000
    }

    sim_parameter_list = list(n = n,
                              doses = dose_levels,
                              dropout_rate = dropout_rate,
                              nsims = nsims,
                              n_patients = sum(n))

    if (endpoint_index == 3) {

      theta = 
        ContinuousErrorCheck(theta, 
                             n_doses, 
                             lower_values = 0,
                             lower_values_sign = ">",
                             upper_values = NA,
                             upper_values_sign = NA,
                             "Overdispersion parameters (theta)",
                             c("Value"),
                             "double",
                             NA)       

      # Create a long vector of overdispersion parameters in the negative binomial distribution (one value for each patient)
      theta_vector = rep(theta, n)                        

    } else {

      theta = 0
      theta_vector = 0

    }

    user_specified$theta = theta
    user_specified$theta_vector = theta_vector

    ######################################################################

    # Simulation models

    if (!is.null(sim_models$placebo_effect)) {

      placebo_effect = 
        ContinuousErrorCheck(sim_models$placebo_effect, 
                             1, 
                             lower_values = c(NA),
                             lower_values_sign = c(NA),
                             upper_values = c(NA),
                             upper_values_sign = c(NA),
                             "Placebo effect in the simulation model (placebo_effect)",
                             c("Value"),
                             "double",
                             NA) 

    } else {
      stop("MCPModSimulation: Placebo effect in the simulation model (placebo_effect): Value must be specified.", call. = FALSE)
    }

    if (!is.null(sim_models$max_effect)) {

      max_effect = ContinuousErrorCheck(sim_models$max_effect, 
                                         NA, 
                                         lower_values = NA,
                                         lower_values_sign = NA,
                                         upper_values = NA,
                                         upper_values_sign = NA,
                                         "Maximum effect over placebo in the simulation model (max_effect)",
                                         NA,
                                         "double",
                                         NA) 


    } else {
      stop("MCPModSimulation: Maximum effect over placebo in the simulation model (max_effect): Value must be specified.", call. = FALSE)
    }

    if (direction_index == 1 & any(max_effect < 0)) stop("MCPModSimulation: Maximum effect over placebo in the simulation model (max_effect): Value must be positive if the direction of the dose-response relationship (direction) is Increasing.", call. = FALSE)

    if (direction_index == -1 & any(max_effect > 0)) stop("MCPModSimulation: Maximum effect over placebo in the simulation model (max_effect): Value must be negative if the direction of the dose-response relationship (direction) is Decreasing.", call. = FALSE)

    max_dose = max(dose_levels)
    n_scenarios = length(max_effect)
    placebo_effect_temp = placebo_effect
    max_effect_temp = max_effect

    # Normal endpoint
    if (endpoint_index == 1) {

      for (i in 1:n_scenarios) { 

        if (direction_index == 1 & max_effect_temp[i] < 0) stop("MCPModSimulation: Maximum effect over placebo in the simulation model (max_effect): Value must be non-negative.", call. = FALSE)

        if (direction_index == -1 & max_effect_temp[i] > 0) stop("MCPModSimulation: Maximum effect over placebo in the simulation model (max_effect): Value must be non-positive.", call. = FALSE)

      }

    }

    # Binary endpoint
    if (endpoint_index == 2) {

      if (placebo_effect_temp < 0 | placebo_effect_temp > 1) stop("MCPModSimulation: Placebo effect in the simulation model (placebo_effect): Value must be >= 0 and <= 1.", call. = FALSE)

      if (placebo_effect_temp == 0) placebo_effect_temp = 0.001  
      if (placebo_effect_temp == 1) placebo_effect_temp = 0.999  

      for (i in 1:n_scenarios) {

        if (direction_index == 1 & placebo_effect_temp + max_effect_temp[i] > 0.999) stop(paste0("MCPModSimulation: Maximum effect over placebo in the simulation model (max_effect): Value must be less than ", 0.999 - placebo_effect_temp,"."), call. = FALSE)

        if (direction_index == -1 & placebo_effect_temp + max_effect_temp[i] < 0.001) stop(paste0("MCPModSimulation: Maximum effect over placebo in the simulation model (max_effect): Value must be less than ", 0.001 - placebo_effect_temp,"."), call. = FALSE)

        max_effect_temp[i] = Logit(placebo_effect_temp + max_effect_temp[i]) - Logit(placebo_effect_temp)

      }

      placebo_effect_temp = Logit(placebo_effect_temp)

    }

    # Count endpoint
    if (endpoint_index == 3) {

      if (placebo_effect_temp < 0) stop("MCPModSimulation: Placebo effect in the simulation model (placebo_effect): Value must be >= 0.", call. = FALSE)

      if (placebo_effect_temp == 0) placebo_effect_temp = 0.001  

      for (i in 1:n_scenarios) {

        if (direction_index == -1 & placebo_effect_temp + max_effect_temp[i] < 0.001) stop(paste0("MCPModSimulation: Maximum effect over placebo in the simulation model (max_effect): Value must be less than ", 0.001 - placebo_effect_temp,"."), call. = FALSE)

        max_effect_temp[i] = log(placebo_effect_temp + max_effect_temp[i]) - log(placebo_effect_temp)

      }

      placebo_effect_temp = log(placebo_effect_temp)

    }

      # Standard deviations are required for normal endpoints
      if (endpoint_index == 1) {

         if (is.null(sim_models$sd)) stop("MCPModSimulation: Standard deviations of the response variable in the simulation model (sd): Value must be specified.", call. = FALSE)

      sd = ContinuousErrorCheck(sim_models$sd, 
                                         NA, 
                                         lower_values = 0,
                                         lower_values_sign = ">",
                                         upper_values = NA,
                                         upper_values_sign = NA,
                                         "Standard deviations of the response variable in the simulation model (sd)",
                                         NA,
                                         "double",
                                         NA) 

         if(length(sd) != n_doses) stop("MCPModSimulation: The length of the dose vector (doses) must be equal to the length of the standard deviation vector (sd) in the simulation model.", call. = FALSE)    


      } else {

        sd = rep(0, n_doses)

      }

    # Compute parameters of the assumed dose-response model to match the placebo and maximum effects 
    sim_model_index = 0

    if (!is.null(sim_models$linear)) {

      sim_model_index = 1

      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 2)

      parameters = 0

      for (i in 1:n_scenarios) {

          coef = ComputeDRFunctionParameters(sim_model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:2) sim_parameter_values[i, j] = coef[j]
      
      }

    }

    if (!is.null(sim_models$quadratic)) {

      sim_model_index = 2

      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 3)

      parameters = sim_models$quadratic

      if (length(parameters) != 1) stop("One parameter must be specified for the simulation model (quadratic).", call. = FALSE)    

      for (i in 1:n_scenarios) { 

          coef = ComputeDRFunctionParameters(sim_model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:3) sim_parameter_values[i, j] = coef[j]

      }

    }

    if (!is.null(sim_models$exponential)) {

      sim_model_index = 3

      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 3)

      parameters = sim_models$exponential

      if (length(parameters) != 1) stop("One parameter must be specified for the simulation model (exponential).", call. = FALSE)    

      for (i in 1:n_scenarios) {

          coef = ComputeDRFunctionParameters(sim_model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:3) sim_parameter_values[i, j] = coef[j]

      }

    }

    if (!is.null(sim_models$emax)) {

      sim_model_index = 4

      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 3)

      parameters = sim_models$emax

      if (length(parameters) != 1) stop("One parameter must be specified for the simulation model (emax).", call. = FALSE)    

      for (i in 1:n_scenarios) {

          coef = ComputeDRFunctionParameters(sim_model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:3) sim_parameter_values[i, j] = coef[j]

      }

    }

    if (!is.null(sim_models$logistic)) {

      sim_model_index = 5

      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 4)

      parameters = sim_models$logistic

      if (length(parameters) != 2) stop("Two parameters must be specified for the simulation model (logistic).", call. = FALSE)

      for (i in 1:n_scenarios) {

          coef = ComputeDRFunctionParameters(sim_model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:4) sim_parameter_values[i, j] = coef[j]

      }

    }

    if (!is.null(sim_models$sigemax)) {

      sim_model_index = 6

      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 4)

      parameters = sim_models$sigemax

      if (length(parameters) != 2) stop("Two parameters must be specified for the simulation model (sigemax).", call. = FALSE)

      for (i in 1:n_scenarios) {

          coef = ComputeDRFunctionParameters(sim_model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:4) sim_parameter_values[i, j] = coef[j]

      }

    }

    if (sim_model_index == 0) stop("MCPModSimulation: List of simulation dose-response models and parameter values (sim_models): At least one model must be specified  (linear, quadratic, exponential, emax, logistic or sigemax).", call. = FALSE)

    sim_model_list = list(sim_model_index = sim_model_index,
                          sim_parameter_values = sim_parameter_values,
                          sd = sd)

    ######################################################################

    # Save the input parameters

    input_parameters = list(direction_index = direction_index,
                            alpha = alpha,
                            delta = delta,
                            user_specified = user_specified,
                            model_selection = model_selection,
                            endpoint_index = endpoint_index,
                            sim_parameters = sim_parameters,
                            placebo_effect = placebo_effect,
                            max_effect = max_effect,
                            sim_model_list = sim_model_list,
                            theta = theta,
                            go_threshold = go_threshold)

    ######################################################################

    # Compute the optimal contrasts and adjusted critical values
    opt_contrast = list()
    contrast_results = list()

    n_groups = n
    doses = dose_levels    

    # Continuous endpoint: A single set of results
    if (endpoint_index == 1) {

      mean_group = 0

      contrast_info = ContrastStep(endpoint_index, selected_models, user_specified, n_groups, doses, alpha, direction_index, mean_group, theta)
      contrast_results[[1]] = contrast_info
      opt_contrast[[1]] = contrast_info$opt_contrast
      crit_value = contrast_info$crit_value
    
    }

    # Binary or count endpoints: A separate set of results for each max effect scenario
    if (endpoint_index %in% c(2, 3)) {

      mean_group = rep(0, n_doses)
      crit_value = rep(0, n_scenarios)

      for (i in 1:n_scenarios) {

        # Compute the rates or average number of events at each dose under each max effect scenario

        for (j in 1:n_doses) {

          # Binary endpoints
          if (endpoint_index == 2) mean_group[j] = AntiLogit(DRFunction(sim_model_list$sim_model_index, sim_model_list$sim_parameter_values[i, ], doses[j]))

          # Count endpoints
          if (endpoint_index == 3) mean_group[j] = exp(DRFunction(sim_model_list$sim_model_index, sim_model_list$sim_parameter_values[i, ], doses[j]))

        }    

      contrast_info = ContrastStep(endpoint_index, selected_models, user_specified, n_groups, doses, alpha, direction_index, mean_group, theta)
      contrast_results[[i]] = contrast_info      
      opt_contrast[[i]] = contrast_info$opt_contrast
      crit_value[i] = contrast_info$crit_value

      } 
    }

    # Number of points used in dose-response plots
    n_points = 20

    # Maximum number of iterations to find maximum likelihood estimates
    maxit = 50

    # Go threshold is defined relative to the placebo effect
    go_threshold = go_threshold + placebo_effect 

    withCallingHandlers({          

      # Run simulations
      sim_results = MCPModRunSimulations(endpoint_index, selected_models, theta, theta_vector, delta, model_selection_index, opt_contrast, crit_value, sim_parameter_list, sim_model_list, direction_index, go_threshold, n_points, maxit)

              },

      warning = function(c) {
        msg <- conditionMessage(c)
        if ( grepl("the line search step became smaller than the minimum value allowed", msg, fixed = TRUE) ) {
          invokeRestart("muffleWarning")
        }
      }

    )

    results = list(contrast_results = contrast_results,
                   input_parameters = input_parameters,
                   selected_models = selected_models,
                   sim_results = sim_results)

    class(results) = "MCPModSimulationResults"

    return(results)

}
# End of MCPModSimulation

MCPModAnalysis = function(endpoint_type, models, dose, resp, alpha = 0.025, direction = "increasing", model_selection = "AIC", Delta, theta = 0) {

    if (missing(endpoint_type)) stop("Endpoint type (endpoint_type): Value must be specified.", call. = FALSE)      
    if (missing(models)) stop("Candidate dose-response models (models): Value must be specified.", call. = FALSE)

    if (missing(dose)) stop("Dose values (dose): Value must be specified.", call. = FALSE)

    if (missing(dose)) stop("Response values (resp): Value must be specified.", call. = FALSE)

    if (missing(Delta)) stop("Treatment effect for identifying the target dose (Delta): Value must be specified.", call. = FALSE)

    # Total number of models
    n_models = length(DF_model_list)

    # Error checks 

    if (!tolower(endpoint_type) %in% tolower(DF_endpoint_list)) stop("MCPModAnalysis: Endpoint type (endpoint_type): Value must be Normal, Binary or Count.", call. = FALSE)

    endpoint_index = 1  

    for (i in 1:length(DF_endpoint_list)) {
        if (tolower(DF_endpoint_list[i]) == tolower(endpoint_type)) endpoint_index = i
    }   

    if (!tolower(direction) %in% c("increasing", "decreasing")) stop("MCPModAnalysis: Direction of the dose-response relationship (direction): Value must be Increasing or Decreasing.", call. = FALSE)

    if (tolower(direction) == "decreasing") direction_index = -1
    if (tolower(direction) == "increasing") direction_index = 1

    if (length(models) == 0) stop("MCPModAnalysis: List of models and initial parameter values (models): At least one model must be specified.", call. = FALSE)

    selected_models = rep(FALSE, n_models)
    user_specified = list()

    # Find the initial values of the model parameters
    if (endpoint_index == 1) user_specified$linear = c(0, 0, 1) else user_specified$linear = c(0, 0)
    if (!is.null(models$linear)) selected_models[1] = TRUE 

    if (endpoint_index == 1) user_specified$quadratic = c(0, 0, 0, 1) else user_specified$quadratic = c(0, 0, 0)
    if (!is.null(models$quadratic)) {
      selected_models[2] = TRUE  

      if (direction_index == 1) {

        user_specified$quadratic[3] =  ContinuousErrorCheck(models$quadratic[1], 
                                                             1, 
                                                             lower_values = NA,
                                                             lower_values_sign = NA,
                                                             upper_values = 0,
                                                             upper_values_sign = "<",
                                                             "Quadratic model (quadratic)",
                                                             c("delta2"),
                                                             "double",
                                                             NA) 

      }

      if (direction_index == -1) {

        user_specified$quadratic[3] =  ContinuousErrorCheck(models$quadratic[1], 
                                                             1, 
                                                             lower_values = 0,
                                                             lower_values_sign = ">",
                                                             upper_values = NA,
                                                             upper_values_sign = NA,
                                                             "Quadratic model (quadratic)",
                                                             c("delta2"),
                                                             "double",
                                                             NA) 

      }

    } 

    if (endpoint_index == 1) user_specified$exponential = c(0, 0, 0, 1) else user_specified$exponential = c(0, 0, 0)
    if (!is.null(models$exponential)) {
      selected_models[3] = TRUE  
      user_specified$exponential[3] =  ContinuousErrorCheck(models$exponential[1], 
                                                           1, 
                                                           lower_values = c(0),
                                                           lower_values_sign = c(">"),
                                                           upper_values = c(NA),
                                                           upper_values_sign = c(NA),
                                                           "Exponential model (exponential)",
                                                           c("delta"),
                                                           "double",
                                                           NA) 
    } 

    if (endpoint_index == 1) user_specified$emax = c(0, 0, 0, 1) else user_specified$emax = c(0, 0, 0)
    if (!is.null(models$emax)) {
      selected_models[4] = TRUE  
      user_specified$emax[3] =  ContinuousErrorCheck(models$emax[1], 
                                                           1, 
                                                           lower_values = c(0),
                                                           lower_values_sign = c(">"),
                                                           upper_values = c(NA),
                                                           upper_values_sign = c(NA),
                                                           "Emax model (emax)",
                                                           c("ED50"),
                                                           "double",
                                                           NA) 

    } 

    if (endpoint_index == 1) user_specified$logistic = c(0, 0, 0, 0, 1) else user_specified$logistic = c(0, 0, 0, 0)
    if (!is.null(models$logistic)) {
      selected_models[5] = TRUE  
      user_specified$logistic[3:4] =  ContinuousErrorCheck(models$logistic[1:2], 
                                                           2, 
                                                           lower_values = c(0, 0),
                                                           lower_values_sign = c(">", ">"),
                                                           upper_values = c(NA, NA),
                                                           upper_values_sign = c(NA, NA),
                                                           "Logistic model (logistic)",
                                                           c("ED50", "delta"),
                                                           "double",
                                                           NA) 

    } 

    if (endpoint_index == 1) user_specified$sigemax = c(0, 0, 0, 0, 1) else user_specified$sigemax = c(0, 0, 0, 0)
    if (!is.null(models$sigemax)) {
      selected_models[6] = TRUE 
      user_specified$sigemax[3:4] =  ContinuousErrorCheck(models$sigemax[1:2], 
                                                           2, 
                                                           lower_values = c(0, 0),
                                                           lower_values_sign = c(">", ">"),
                                                           upper_values = c(NA, NA),
                                                           upper_values_sign = c(NA, NA),
                                                           "SigEmax model (sigemax)",
                                                           c("ED50", "h"),
                                                           "double",
                                                           NA) 
    } 

    n_selected_models = sum(selected_models)

    if (n_selected_models == 0) stop("MCPModAnalysis: List of models and initial parameter values (models): At least one model must be specified (linear, quadratic, exponential, emax, logistic or sigemax).", call. = FALSE)

    alpha = 
          ContinuousErrorCheck(alpha, 
                               1, 
                               lower_values = c(0.001),
                               lower_values_sign = c(">"),
                               upper_values = c(0.999),
                               upper_values_sign = c("<"),
                               "One-sided Type I error rate (alpha)",
                               c("Value"),
                               "double",
                               NA) 

    if (!model_selection %in% c("AIC", "maxT", "aveAIC")) stop("MCPModAnalysis: Model selection criterion (model_selection): Value must be AIC, maxT or aveAIC.", call. = FALSE)

    delta = 
      ContinuousErrorCheck(Delta, 
                           1, 
                           lower_values = c(NA),
                           lower_values_sign = c(NA),
                           upper_values = c(NA),
                           upper_values_sign = c(NA),
                           "Treatment effect for identifying the target dose (Delta)",
                           c("Value"),
                           "double",
                           NA) 

    if (direction_index == 1 & delta <= 0) stop("MCPModAnalysis: Treatment effect for identifying the target dose (Delta): Value must be positive if the direction of the dose-response relationship (direction) is Increasing.", call. = FALSE)

    if (direction_index == -1 & delta >= 0) stop("MCPModAnalysis: Treatment effect for identifying the target dose (Delta): Value must be negative if the direction of the dose-response relationship (direction) is Decreasing.", call. = FALSE)

    dose = ContinuousErrorCheck(dose, 
                       NA, 
                       lower_values = 0,
                       lower_values_sign = c(">="),
                       upper_values = 1000,
                       upper_values_sign = c("<="),
                       "Dose levels (dose)",
                       NA,
                       "double",
                       NA) 

    dose_levels = sort(unique(dose))
    max_dose = max(dose_levels)
    n_doses = length(dose_levels)
    n_groups = table(dose)

    # Normal endpoints
    if (endpoint_index == 1) {

      resp = ContinuousErrorCheck(resp, 
                         NA, 
                         lower_values = NA,
                         lower_values_sign = c(NA),
                         upper_values = c(NA),
                         upper_values_sign = c(NA),
                         "Responses (resp)",
                         NA,
                         "double",
                         NA) 

    }

    # Binary endpoints
    if (endpoint_index == 2) {

      resp = ContinuousErrorCheck(resp, 
                         NA, 
                         lower_values = 0,
                         lower_values_sign = ">=",
                         upper_values = 1,
                         upper_values_sign = "<=",
                         "Resp variable (resp)",
                         NA,
                         "double",
                         NA) 

    }

    # Count endpoints
    if (endpoint_index == 3) {

      resp = ContinuousErrorCheck(resp, 
                         NA, 
                         lower_values = 0,
                         lower_values_sign = ">=",
                         upper_values = NA,
                         upper_values_sign = NA,
                         "Resp variable (resp)",
                         NA,
                         "double",
                         NA) 

      theta = 
        ContinuousErrorCheck(theta, 
                             n_doses, 
                             lower_values = 0,
                             lower_values_sign = ">",
                             upper_values = NA,
                             upper_values_sign = NA,
                             "Overdispersion parameters (theta)",
                             c("Value"),
                             "double",
                             NA)     

      # Create a long vector of overdispersion parameters in the negative binomial distribution (one value for each patient)
      theta_vector = rep(theta, n_groups)                        

    } else {

      theta = 0
      theta_vector = 0

    }

    user_specified$theta = theta

    if(length(dose) != length(resp)) stop("MCPModAnalysis: The length of the dose vector (dose) must be equal to the length of the response vector (resp).", call. = FALSE)    

    ######################################################################

    # Save the input parameters

    input_parameters = list(direction_index = direction_index,
                            alpha = alpha,
                            delta = delta,
                            user_specified = user_specified,
                            model_selection = model_selection,
                            endpoint_index = endpoint_index,
                            theta = theta)

    ######################################################################

    # Descriptive statistics

    lower_cl = rep(0, n_doses)
    upper_cl = rep(0, n_doses)
    mean_group = rep(0, n_doses)
    sderr = rep(0, n_doses)

    # Normal endpoints
    if (endpoint_index == 1) {

        for (j in 1:n_doses) {

            mean_group[j] = mean(resp[dose == dose_levels[j]])
            sderr[j] = sd(resp) / sqrt(n_groups[j]) 
            lower_cl[j] = mean_group[j] - sderr[j] * qnorm(1 - alpha) 
            upper_cl[j] = mean_group[j] + sderr[j] * qnorm(1 - alpha)  

        }

    }

    # Binary endpoints
    if (endpoint_index == 2) {

        for (j in 1:n_doses) {

            mean_group[j] = mean(resp[dose == dose_levels[j]])
            sderr[j] = sqrt(mean_group[j] * (1 - mean_group[j]) / n_groups[j]) 
            lower_cl[j] = mean_group[j] - sderr[j] * qnorm(1 - alpha)
            upper_cl[j] = mean_group[j] + sderr[j] * qnorm(1 - alpha)

        }

    }

    # Count endpoints
    if (endpoint_index == 3) {

        for (j in 1:n_doses) {

            mean_group[j] = mean(resp[dose == dose_levels[j]])
            sderr[j] = sqrt((theta[j] + mean_group[j]) / (n_groups[j] * theta[j] * mean_group[j]))
            lower_cl[j] = exp(log(mean_group[j]) - sderr[j] * qnorm(1 - alpha)) 
            upper_cl[j] = exp(log(mean_group[j]) + sderr[j] * qnorm(1 - alpha)) 

        }

    }

    descriptive_statistics = list(dose_levels = dose_levels,
                                  n_groups = n_groups,
                                  mean_group = mean_group,
                                  sderr = sderr,
                                  lower_cl = lower_cl,
                                  upper_cl = upper_cl)

    ######################################################################

    # Hypothesis testing step 

    # Compute the optimal contrasts, contrast correlation matrix and adjusted critical value
    contrast_results = ContrastStep(endpoint_index, selected_models, user_specified, n_groups, dose_levels, alpha, direction_index, mean_group, theta)

    # Compute the test statistics
    mcp_results = MCPStep(endpoint_index, contrast_results, selected_models, user_specified, dose, resp, alpha, direction_index)

    ######################################################################

    # Modeling step

    mod_results = ModStep(endpoint_index, selected_models, theta_vector, dose, resp, delta, direction_index)

    # Include selected models only
    selected_models = unlist(selected_models)

    model_fit = list()
    k = 1
    
    for (i in 1:length(mod_results$model_fit)) {

      if (selected_models[i]) {
        model_fit[[k]] = mod_results$model_fit[[i]]
        k = k + 1
      }

    }

    mod_results$model_fit = model_fit

    results = list(contrast_results = contrast_results,
                   mcp_results = mcp_results,
                   mod_results = mod_results,
                   input_parameters = input_parameters,
                   selected_models = selected_models,
                   descriptive_statistics = descriptive_statistics)

    class(results) = "MCPModAnalysisResults"

    return(results)

}
# End of MCPModAnalysis

print.MCPModAnalysisResults = function (x, digits = 3, ...) {

    results = x

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract the optimal contrasts and contrast correlation matrix
    contrast_results = results$contrast_results

    # Extract the test statistics
    mcp_results = results$mcp_results

    # Extract the Mod step results
    mod_results = results$mod_results

    # Extract descriptive statistics
    descriptive_statistics = results$descriptive_statistics

    # Extract the list of model fit parameters
    model_fit = mod_results$model_fit

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    n_doses = length(mcp_results$dose_levels)

    DF_selected_model_list = DF_model_list[selected_models]

    endpoint_index = input_parameters$endpoint_index

    cat("***************************************\n\n")
    cat("Descriptive statistics\n\n")
    cat("***************************************\n\n")

    # Normal endpoint
    if (input_parameters$endpoint_index == 1) {

      x = cbind(descriptive_statistics$dose_levels, 
                descriptive_statistics$n_groups, 
                round(descriptive_statistics$mean_group, digits),
                paste0("(", round(descriptive_statistics$lower_cl, digits), ", ", round(descriptive_statistics$upper_cl, digits), ")"),
                round(descriptive_statistics$sderr, digits))
      x = as.data.frame(x)
      colnames(x) = c("Dose", "n", "Mean", "95% CI", "SE")
      print(x, row.names = FALSE)

    }

    # Binary endpoint
    if (input_parameters$endpoint_index == 2) {

      x = cbind(descriptive_statistics$dose_levels, 
                descriptive_statistics$n_groups, 
                round(descriptive_statistics$mean_group, digits),
                paste0("(", round(descriptive_statistics$lower_cl, digits), ", ", round(descriptive_statistics$upper_cl, digits), ")"))
      x = as.data.frame(x)
      colnames(x) = c("Dose", "n", "Rate", "95% CI")
      print(x, row.names = FALSE)

    }

    # Count endpoint
    if (input_parameters$endpoint_index == 3) {

      x = cbind(descriptive_statistics$dose_levels, 
                descriptive_statistics$n_groups, 
                round(descriptive_statistics$mean_group, digits),
                paste0("(", round(descriptive_statistics$lower_cl, digits), ", ", round(descriptive_statistics$upper_cl, digits), ")"),
                round(descriptive_statistics$sderr, digits),
                round(input_parameters$theta, digits))
      x = as.data.frame(x)
      colnames(x) = c("Dose", "n", "Mean", "95% CI", "SE", "Theta")
      print(x, row.names = FALSE)

    }

    cat("\n***************************************\n\n")
    cat("Hypothesis testing and model selection\n\n")
    cat("***************************************\n\n")

    cat("Model-specific dose-response contrasts\n\n")

    x = cbind(descriptive_statistics$dose_levels, FormatMatrix(as.matrix(contrast_results$opt_contrast), "%0.3f"))
    x = as.data.frame(x)
    colnames(x) = c("Dose", DF_selected_model_list)
    print(x, row.names = FALSE)

    cat("\n Contrast correlation matrix\n\n")

    x = cbind(DF_selected_model_list, FormatMatrix(as.matrix(contrast_results$corr_matrix), "%0.3f"))
    x = as.data.frame(x)
    colnames(x) = c("Models", DF_selected_model_list)
    print(x, row.names = FALSE)

    cat("\n Model-specific contrast tests\n\n")

    sign = rep("No", n_selected_models)
    for (i in 1:n_selected_models) {
      if (mcp_results$sign_model[i] == 1) sign[i] = "Yes"
    }

    x = cbind(sprintf("%0.3f", mcp_results$test_statistics), 
              sprintf("%0.4f",mcp_results$adj_pvalues),
              sign)
    x = cbind(DF_selected_model_list, x)
    x = as.data.frame(x)
    colnames(x) = c("Model", "Test statistic", "Adjusted p-value", "Significant contrast")
    print(x, row.names = FALSE)

    cat("\nAdjusted critical value: ", round(contrast_results$crit_value, 3), sep = "") 

    ##################################################################

    cat("\n\n***************************************\n\n")
    cat("Dose-response modeling\n\n")
    cat("***************************************\n\n")

    for (i in 1:n_selected_models) {

        current_model = model_fit[[i]]

        model = current_model$model

        # Print well-defined models only
        if (current_model$status >= 0) {

            n_parameters = length(DF_model_parameters[[model]]) 

            cat(paste0("Dose-response model: ", DF_model_list[model]), "\n")

            cat("Parameter estimates\n")
            coef = round(current_model$coef[1:n_parameters], digits)
            names(coef) = DF_model_parameters[[model]]
            print(coef)

            cat("\n***************************************\n\n")

        } else {

            cat(paste0("Dose-response model: ", DF_model_list[model]), "\n")
            cat("Parameter could not be estimated\n")

            cat("\n***************************************\n\n")

        }

    }

    ##################################################################

    cat("Dose selection\n\n")
    cat("***************************************\n\n")

    sign = rep("No", n_selected_models)
    criterion = rep(NA, n_selected_models)
    test_statistics = rep(NA, n_selected_models)
    target_dose = rep(NA, n_selected_models)
    for (i in 1:n_selected_models) {
      if (mcp_results$sign_model[i] == 1) {
        sign[i] = "Yes"
        criterion[i] = model_fit[[i]]$criterion
        test_statistics[i] = mcp_results$test_statistics[i]
        if (model_fit[[i]]$target_dose >= 0) target_dose[i] = model_fit[[i]]$target_dose

      }
    }

    model_weight = rep(NA, n_selected_models)
    for (i in 1:n_selected_models) {

      if (model_fit[[i]]$status >= 0) {
        current_criterion = criterion[i]
        denominator = 0
        for (j in 1:n_selected_models) {
          if (mcp_results$sign_model[j] == 1 & model_fit[[j]]$status >= 0) denominator = denominator + exp(- 0.5 * (criterion[j] - current_criterion))
        }
        if (mcp_results$sign_model[i] == 1 & abs(denominator) > 0.0001) model_weight[i] = 1 / denominator
      } else {
        criterion[i] = NA
      }
    }

    cat("Model selection criteria\n\n") 

    x = cbind(DF_selected_model_list, 
              sign, 
              sprintf("%0.2f", criterion),
              sprintf("%0.3f", test_statistics),
              sprintf("%0.3f", model_weight))
    x = as.data.frame(x)
    colnames(x) = c("Model", "Significant contrast", "AIC", "Test statistic", "Model weight")
    print(x, row.names = FALSE)

    if (input_parameters$model_selection == "AIC") {

      criterion_label = "based on the smallest AIC"

      if (sum(mcp_results$sign_model) > 0) {

        index = which.min(criterion)

        cat("\nSelected model (", criterion_label, "): ", DF_selected_model_list[index], "\n\n", sep = "") 

      } else {

        cat("\nSelected model (", criterion_label, "): No model is significant. \n\n", sep = "") 

      }

    }

    if (input_parameters$model_selection == "maxT") {

      criterion_label = "based on the most significant test statistic"

      if (sum(mcp_results$sign_model) > 0) {

        index = which.max(test_statistics)

        cat("\nSelected model (", criterion_label, "): ", DF_selected_model_list[index], "\n\n", sep = "") 

      } else {

        cat("\nSelected model (", criterion_label, "): No model is significant. \n\n", sep = "") 

      }

    }

    if (input_parameters$model_selection == "aveAIC") {
    
      criterion_label = "based on weighted model averaging"
      cat("\n")
    
    }

    cat("***************************************\n\n")

    cat("Model-specific estimated target doses (based on Delta = ", input_parameters$delta, ")\n\n", sep = "") 

    x = cbind(DF_selected_model_list, 
              sprintf("%0.3f", target_dose))
    x = as.data.frame(x)
    colnames(x) = c("Model", "Target dose")
    print(x, row.names = FALSE)

    if (input_parameters$model_selection == "AIC" | input_parameters$model_selection == "maxT") {    

      if (sum(mcp_results$sign_model) > 0) {

        cat("\nSelected target dose (", criterion_label, "): ", round(target_dose[index], 3)," \n\n", sep="") 

      } else {

        cat("\nSelected target dose cannot be determined. \n\n", sep = "") 

      }

    }

    if (input_parameters$model_selection == "aveAIC") {    

      weighted_dose = 0
      for (i in 1:n_selected_models) {

        if (mcp_results$sign_model[i] == 1) weighted_dose = weighted_dose + target_dose[i] * model_weight[i]

      }

      if (sum(mcp_results$sign_model) > 0) {

        cat("\nSelected target dose (", criterion_label, "): ", round(weighted_dose, 3), " \n\n", sep="") 

      } else {

        cat("\nSelected target dose cannot be determined. \n\n", sep = "") 

      }

    }

}

print.MCPModSimulationResults = function (x, digits = 3, ...) {

    results = x
    #############################################################################

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract the user-specified values of the model parameters
    user_specified = input_parameters$user_specified

    # Extract the initial values of the model parameters
    initial_values = input_parameters$initial_values

    # Extract the optimal contrasts and contrast correlation matrix
    contrast_results = results$contrast_results

    # Simulation results
    sim_results = results$sim_results

    # Total number of models
    n_models = length(DF_model_list)

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    DF_selected_model_list = DF_model_list[selected_models]

    endpoint_index = input_parameters$endpoint_index

    dose_levels = input_parameters$sim_parameters$doses

    n_doses = length(dose_levels)

    cat("***************************************\n\n")
    cat("Simulation results\n\n")
    cat("***************************************\n\n")

    cat("Power summary\n\n")

    x = data.frame(input_parameters$max_effect, 
                   round(sim_results$power, digits))
    colnames(x) = c("Maximum effect over placebo", "Power")
    print(x, row.names = FALSE)

    if (input_parameters$model_selection != "aveAIC") {

      cat(paste0("\nGo probability summary (based on the threshold of ", input_parameters$go_threshold, ")\n\n"))

      x = data.frame(input_parameters$max_effect, 
                     round(sim_results$go_prob, digits))
      colnames(x) = c("Maximum effect over placebo", "Pr(Go)")
      print(x, row.names = FALSE)

      cat("\nProbability of selecting a dose-response model\n\n")

      x = data.frame(input_parameters$max_effect,
                     round(sim_results$model_index_summary, digits))

      colnames(x) = c("Maximum effect over placebo", "No model selected", DF_model_list_short[selected_models])
      print(x, row.names = FALSE)

    }

    cat("\nEstimated target doses (based on Delta = ", input_parameters$delta, ")\n\n", sep = "") 

    true_target_dose = round(sim_results$true_target_dose, digits)
    true_target_dose[true_target_dose == -1] = "NA"
    true_target_dose = as.character(true_target_dose)

    lower_bound_target_dose = round(sim_results$target_dose_summary[, 1], digits)
    lower_bound_target_dose[lower_bound_target_dose == -1] = NA

    mean_target_dose = round(sim_results$target_dose_summary[, 2], digits)
    mean_target_dose[mean_target_dose == -1] = NA

    upper_bound_target_dose = round(sim_results$target_dose_summary[, 3], digits)
    upper_bound_target_dose[upper_bound_target_dose == -1] = NA

    x = data.frame(input_parameters$max_effect, 
                   true_target_dose,
                   paste0(mean_target_dose, " (", lower_bound_target_dose, ", ", upper_bound_target_dose, ")"))
    colnames(x) = c("Maximum effect over placebo", "True target dose", "Mean target dose (95% CI)")
    print(x, row.names = FALSE)

    cat("\nProbability of identifying the target dose\n\n")

    x = data.frame(input_parameters$max_effect, 
                   round(sim_results$target_dose_categorical_summary, digits))

    colnames(x) = c("Maximum effect over placebo", "No dose found", dose_levels[2:n_doses], "Greater than max dose")
    print(x, row.names = FALSE)

}

SaveReport = function(report, report_title) {

  # Create a docx object
  doc = officer::read_docx(system.file(package = "MCPModPack", "template/report_template.docx"))

  dim_doc = officer::docx_dim(doc)

  # Report's title
  doc = officer::set_doc_properties(doc, title = report_title)
  doc = officer::body_add_par(doc, value = report_title, style = "heading 1")

  # Text formatting
  my.text.format = officer::fp_text(font.size = 12, font.family = "Arial")

  # Table formatting
  header.cellProperties = officer::fp_cell(border.left = officer::fp_border(width = 0), border.right = officer::fp_border(width = 0), border.bottom = officer::fp_border(width = 2), border.top = officer::fp_border(width = 2), background.color = "#eeeeee")
  data.cellProperties = officer::fp_cell(border.left = officer::fp_border(width = 0), border.right = officer::fp_border(width = 0), border.bottom = officer::fp_border(width = 0), border.top = officer::fp_border(width = 0))

  header.textProperties = officer::fp_text(font.size = 12, bold = TRUE, font.family = "Arial")
  data.textProperties = officer::fp_text(font.size = 12, font.family = "Arial")

  thick_border = fp_border(color = "black", width = 2)

  leftPar = officer::fp_par(text.align = "left")
  rightPar = officer::fp_par(text.align = "right")
  centerPar = officer::fp_par(text.align = "center")

  # Number of sections in the report (the report's title is not counted)
  n_sections = length(report) 

  # Loop over the sections in the report
  for(section_index in 1:n_sections) {

      # Determine the item's type (text by default)
      type = report[[section_index]]$type

      # Determine the item's label 
      label = report[[section_index]]$label

      # Determine the item's label 
      footnote = report[[section_index]]$footnote

      # Determine the item's value 
      value = report[[section_index]]$value

      # Determine column width 
      column_width = report[[section_index]]$column_width

      # Determine the page break status 
      page_break = report[[section_index]]$page_break
      if (is.null(page_break)) page_break = FALSE

      # Determine the figure's location (for figures only)
      filename = report[[section_index]]$filename

      # Determine the figure's dimensions (for figures only)
      dim = report[[section_index]]$dim

      if (!is.null(type)) {

        # Fully formatted data frame 
        if (type == "table") {

            doc = officer::body_add_par(doc, value = label, style = "heading 2")

            summary_table = flextable::regulartable(data = value)
            summary_table = flextable::style(summary_table, pr_p = leftPar, pr_c = header.cellProperties, pr_t = header.textProperties, part = "header")
            summary_table = flextable::style(summary_table, pr_p = leftPar, pr_c = data.cellProperties, pr_t = data.textProperties, part = "body")

            summary_table = flextable::hline_bottom(summary_table, part = "body", border = thick_border )

            summary_table = flextable::width(summary_table, width = column_width)

            doc = flextable::body_add_flextable(doc, summary_table)

            if (!is.null(footnote)) doc = officer::body_add_par(doc, value = footnote, style = "Normal")

            if (page_break) doc = officer::body_add_break(doc, pos = "after")

        }

        # Enhanced metafile graphics produced by package devEMF 
        if (type == "emf_plot") {

            doc = officer::body_add_par(doc, value = label, style = "heading 2")

            doc = officer::body_add_img(doc, src = filename, width = dim[1], height = dim[2]) 

            if (!is.null(footnote)) doc = officer::body_add_par(doc, value = footnote, style = "Normal")

            if (page_break) doc = officer::body_add_break(doc, pos = "after")

            # Delete the figure
            if (file.exists(filename)) file.remove(filename)   

        }

      }    

  }

  return(doc)       

}
# End SaveReport


CreateTable = function(data_frame, column_names, column_width, title, page_break, footnote = NULL) {

    if (is.null(column_width)) {
      column_width = rep(2, dim(data_frame)[2])
    } 
     
    data_frame = as.data.frame(data_frame)

    colnames(data_frame) = column_names

    item_list = list(label = title, 
                     value = data_frame,
                     column_width = column_width,
                     type = "table",
                     footnote = footnote,
                     page_break = page_break)

    return(item_list)

  }
# End of CreateTable         


GenerateAnalysisReport = function(results, report_title) { 

    #############################################################################

    # Error checks

    if (class(results) != "MCPModAnalysisResults") stop("AnalysisReport: The object was not created by the MCPModAnalysis function.", call. = FALSE)

    if (!requireNamespace("officer", quietly = TRUE)) {
     stop("AnalysisReport: The officer package is required to generate this report. Please install it.", call. = FALSE)
    }

    if (!requireNamespace("flextable", quietly = TRUE)) {
     stop("AnalysisReport: The flextable package is required to generate this report. Please install it.", call. = FALSE)
    }

    if (!requireNamespace("devEMF", quietly = TRUE)) {
     stop("AnalysisReport: The devEMF package is required to generate this report. Please install it.", call. = FALSE)
    }

    #############################################################################

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract the user-specified values of the model parameters
    user_specified = input_parameters$user_specified

    # Extract the initial values of the model parameters
    initial_values = input_parameters$initial_values

    # Extract the optimal contrasts and contrast correlation matrix
    contrast_results = results$contrast_results

    # Extract the test statistics
    mcp_results = results$mcp_results

    # Extract the Mod step results
    mod_results = results$mod_results

    # Extract the list of model fit parameters
    model_fit = mod_results$model_fit

    # Extract descriptive statistics
    descriptive_statistics = results$descriptive_statistics

    # Number of doses
    n_doses = length(descriptive_statistics$dose_levels)

    # Total number of models
    n_models = length(DF_model_list)

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    DF_selected_model_list = DF_model_list[selected_models]
    DF_selected_model_index_list = (1:n_models)[selected_models]

    # Extract the dose levels and responses
    dose = mod_results$dose
    resp = mod_results$resp

    endpoint_index = input_parameters$endpoint_index

    digits = 3

    #############################################################################

    # Empty list of tables to be included in the report
    item_list = list()
    item_index = 1
    table_index = 1
    figure_index = 1

    #############################################################################

    column_names = c("Parameter", "Value")

    col1 = NULL
    col2 = NULL

    if (input_parameters$direction_index == 1) direction_label = "Increasing" else direction_label = "Decreasing"

    if (endpoint_index != 3) {

      col1 = c(col1, "Endpoint type", "Candidate models", "One-sided Type I error rate", "Direction of the dose-response relationship")
      col2 = c(col2, DF_endpoint_list[endpoint_index], 
                     paste0(DF_selected_model_list, collapse = ", "), 
                     input_parameters$alpha,
                     direction_label)

    } else {

      col1 = c(col1, "Endpoint type", "Candidate models", "One-sided Type I error rate", "Direction of the dose-response relationship", "Overdispersion parameters")
      col2 = c(col2, DF_endpoint_list[endpoint_index], 
                     paste0(DF_selected_model_list, collapse = ", "), 
                     input_parameters$alpha,
                     direction_label, 
                     paste0(round(input_parameters$theta, digits), collapse = ", "))

    }



    data_frame = cbind(col1, col2)
    title = paste0("Table ", table_index, ". General parameters.")

    column_width = c(3.5, 3)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    column_names = c("Candidate model", "Parameter values")

    col1 = NULL
    col2 = NULL

    col1 = DF_selected_model_list
    for (i in 1:n_models) {
     
        if (selected_models[i] == TRUE) {

          n_parameters = length(DF_model_parameters[[i]])
          model_user_specified = round(user_specified[[i]][1:n_parameters], digits)
          for (j in 1:n_parameters) model_user_specified[j] = paste0(DF_model_parameters[[i]][j], " = ", model_user_specified[j])
          col2 = c(col2, paste0(model_user_specified, collapse = ", "))

        }    
    }

    data_frame = cbind(col1, col2)
    title = paste0("Table ", table_index, ". Initial values of the model parameters.")

    column_width = c(2.5, 4)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    # Descriptive statistics

    # Normal endpoint
    if (input_parameters$endpoint_index == 1) {

      column_names = c("Dose", "Number of patients", "Mean", "95% CI", "SE")
      data_frame = cbind(descriptive_statistics$dose_levels, 
                descriptive_statistics$n_groups, 
                round(descriptive_statistics$mean_group, digits),
                paste0("(", round(descriptive_statistics$lower_cl, digits), ", ", round(descriptive_statistics$upper_cl, digits), ")"),
                round(descriptive_statistics$sderr, digits))

      column_width = c(1, 1.5, 1.5, 1.5, 1)

    }

    # Binary endpoint
    if (input_parameters$endpoint_index == 2) {

      column_names = c("Dose", "Number of patients", "Rate", "95% CI")

      data_frame = cbind(descriptive_statistics$dose_levels, 
                descriptive_statistics$n_groups, 
                round(descriptive_statistics$mean_group, digits),
                paste0("(", round(descriptive_statistics$lower_cl, digits), ", ", round(descriptive_statistics$upper_cl, digits), ")"))

      column_width = c(1, 1.5, 2, 2)

    }

    # Count endpoint
    if (input_parameters$endpoint_index == 3) {

      column_names = c("Dose", "Number of patients", "Mean", "95% CI", "SE", "Theta")

      data_frame = cbind(descriptive_statistics$dose_levels, 
                descriptive_statistics$n_groups, 
                round(descriptive_statistics$mean_group, digits),
                paste0("(", round(descriptive_statistics$lower_cl, digits), ", ", round(descriptive_statistics$upper_cl, digits), ")"),
                round(descriptive_statistics$sderr, digits),
                round(user_specified$theta, digits))

      column_width = c(1, 1, 1, 1.5, 1, 1)

    }


    title = paste0("Table ", table_index, ". Descriptive statistics.")

    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    # Hypothesis testing step

    #############################################################################

    column_names = c("Dose", DF_selected_model_list)

    data_frame = cbind(descriptive_statistics$dose_levels, round(contrast_results$opt_contrast, 3))

    title = paste0("Table ", table_index, ". Model-specific dose-response contrasts.")

    original_width = c(0.9, 0.9, 1.1, 0.8, 0.9, 0.9) # c(6.5 - 0.95 * 6, rep(0.95, 6))
    column_width = c(1, 5.5 * original_width[selected_models] / sum(original_width[selected_models]))
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    column_names = c("Models", DF_selected_model_list)

    data_frame = cbind(DF_selected_model_list, round(contrast_results$corr_matrix, 3))

    title = paste0("Table ", table_index, ". Contrast correlation matrix.")

    original_width = c(0.9, 0.9, 1.1, 0.8, 0.9, 0.9) # c(6.5 - 0.95 * 6, rep(0.95, 6))
    column_width = c(1, 5.5 * original_width[selected_models] / sum(original_width[selected_models]))
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    column_names = c("Model", "Test statistic", "Adjusted p-value", "Significant contrast")

    sign = rep("No", n_selected_models)
    for (i in 1:n_selected_models) {
      if (mcp_results$sign_model[i] == 1) sign[i] = "Yes"
    }

    data_frame = cbind(DF_selected_model_list, 
                       round(mcp_results$test_statistics, 3), 
                       round(mcp_results$adj_pvalues, 4),
                       sign)

    title = paste0("Table ", table_index, ". Model-specific contrast tests.")

    footnote = paste0("Adjusted critical value: ", round(contrast_results$crit_value, 3), ".") 

    column_width = c(1.25, 1.5, 1.5, 2.25)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    # Dose-response modeling step

    #############################################################################

    column_names = c("Model", "Parameter", "Estimate", "Convergence criterion")

    col1 = NULL
    col2 = NULL
    col3 = NULL
    col4 = NULL

    k = 1

    for (i in 1:n_models) {

        if (selected_models[i] == TRUE) {

          # Current dose-response model
          current_model = model_fit[[k]]
          coef = current_model$coef

          test_statistic = mcp_results$test_statistics[k]
          adj_pvalue = mcp_results$adj_pvalues[k]
          if (mcp_results$sign_model[k] == 1) sign = "Yes" else sign = "No"

          n_parameters = length(DF_model_parameters[[i]])
          col1 = c(col1, DF_model_list[i], rep("", n_parameters - 1))
          col2 = c(col2, DF_model_parameters[[i]])
          # Display the estimated parameters if the model converged
          if (current_model$status >= 0) parameter_estimates = round(current_model$coef[1:n_parameters], 3) else parameter_estimates = rep("NA", n_parameters)
          col3 = c(col3, parameter_estimates) 
          if (current_model$status >= 0) convergence_criterion = round(current_model$convergence_criterion, 3) else convergence_criterion = "NA"        
          col4 = c(col4, convergence_criterion, rep("", n_parameters - 1)) 

          k = k + 1

        }

    }

    data_frame = cbind(col1, col2, col3, col4)

    title = paste0("Table ", table_index, ". Estimated parameters of dose-response models.")

    column_width = c(1.5, 1.25, 1.25, 2.5)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE, "The convergence criterion is defined as the length of the gradient vector evaluated at the maximum likelihood estimate and a high value of the convergence criterion suggests lack of convergence.")
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    column_names = c("Parameter", "Value")

    col1 = NULL
    col2 = NULL

    if (input_parameters$model_selection == "AIC") model_selection_label = "Select the model with the smallest AIC"
    if (input_parameters$model_selection == "maxT") model_selection_label = "Select the model corresponding to the largest test statistic"
    if (input_parameters$model_selection == "aveAIC") model_selection_label = "Average the models using AIC-based weights"

    col1 = c(col1, "Model selection criterion", "Delta")
    col2 = c(col2, model_selection_label, input_parameters$delta)

    data_frame = cbind(col1, col2)
    title = paste0("Table ", table_index, ". Model selection parameters.")

    footnote = "Delta is the pre-defined clinically meaningful improvement over placebo."

    column_width = c(3, 3.5)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1  


    #############################################################################

    column_names = c("Model", "Significant contrast", "AIC", "Test statistic", "Model weight")

    sign = rep("No", n_selected_models)
    criterion = rep(NA, n_selected_models)
    test_statistics = rep(NA, n_selected_models)
    target_dose = rep(NA, n_selected_models)
    for (i in 1:n_selected_models) {
      if (mcp_results$sign_model[i] == 1) {
        sign[i] = "Yes"
        if (model_fit[[i]]$status >= 0) criterion[i] = model_fit[[i]]$criterion else criterion[i] = NA
        test_statistics[i] = mcp_results$test_statistics[i]
        if (model_fit[[i]]$target_dose >= 0) target_dose[i] = model_fit[[i]]$target_dose
      }
    }

    model_weight = rep(NA, n_selected_models)
    for (i in 1:n_selected_models) {

      if (model_fit[[i]]$status >= 0) {
        current_criterion = criterion[i]
        denominator = 0
        for (j in 1:n_selected_models) {
          if (mcp_results$sign_model[j] == 1 & model_fit[[j]]$status >= 0) denominator = denominator + exp(- 0.5 * (criterion[j] - current_criterion))
        }
        if (mcp_results$sign_model[i] == 1 & abs(denominator) > 0.0001) model_weight[i] = 1 / denominator
      } else {
        criterion[i] = NA
      }
    }

    data_frame = cbind(DF_selected_model_list, 
                       sign, 
                       sprintf("%0.2f", criterion),
                       sprintf("%0.3f", test_statistics),
                       sprintf("%0.3f", model_weight))

    title = paste0("Table ", table_index, ". Model selection criteria.")

    column_width = c(1.25, 1.5, 1.25, 1.5, 1)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1 

    #############################################################################

    column_names = c("Model", "Target dose")

    data_frame = cbind(DF_selected_model_list, 
                       sprintf("%0.3f", target_dose))

    title = paste0("Table ", table_index, ". Model-specific estimated target doses.")

    column_width = c(2, 4.5)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1 

    #############################################################################

    column_names = c("Parameter", "Value")

    selected_model = "No model is significant" 
    selected_target_dose = "Target dose cannot be determined"

    if (input_parameters$model_selection == "AIC") {
      if (sum(mcp_results$sign_model) > 0) {
        index = which.min(criterion)
        selected_model = DF_selected_model_list[index]
      } 
    }

    if (input_parameters$model_selection == "maxT") {
      if (sum(mcp_results$sign_model) > 0) {
        index = which.max(test_statistics)
        selected_model = DF_selected_model_list[index]
      } 
    }

    if (input_parameters$model_selection == "aveAIC") {
        selected_model = "Dose selection is based on averaging the significant models"
    }

    if (input_parameters$model_selection == "AIC" | input_parameters$model_selection == "maxT") {    
      if (sum(mcp_results$sign_model) > 0) {
        if (!is.na(target_dose[index])) selected_target_dose = round(target_dose[index], 3)
      }
    }

    if (input_parameters$model_selection == "aveAIC") {    
      weighted_dose = 0
      for (i in 1:n_models) {
        if (mcp_results$sign_model[i] == 1) weighted_dose = weighted_dose + target_dose[i] * model_weight[i]
      }
      if (sum(mcp_results$sign_model) > 0) {
        if (!is.na(weighted_dose)) selected_target_dose = round(weighted_dose, 3)

      }
    }      

    data_frame = cbind(c("Selected model", "Selected target dose"), c(selected_model, selected_target_dose))

    title = paste0("Table ", table_index, ". Selected model and target dose.")

    column_width = c(2, 4.5)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE)
    item_index = item_index + 1
    table_index = table_index + 1 

    #############################################################################

    # Plot all models

    width = 6.5
    height = 5

    eval_function_list = list()  

    x_limit = c(min(descriptive_statistics$dose_levels), max(descriptive_statistics$dose_levels))
    y_limit = c(min(descriptive_statistics$lower_cl), max(descriptive_statistics$upper_cl))

    # Evaluate the dose-response functions and determine the axis ranges
    for (i in 1:n_selected_models) {

        current_model = model_fit[[i]]

        if (current_model$status >= 0) {

          # Evaluate the current dose-response model
          model_index = DF_selected_model_index_list[i]
          coef = current_model$coef

          eval_function_list[[i]] = EvaluateDRFunction(model_index, endpoint_index, coef, dose) 

          # Exclude the cases with missing coefficients
          if (!any(is.na(eval_function_list[[i]]$y))) {

            if (min(eval_function_list[[i]]$y) < y_limit[1]) y_limit[1] = min(eval_function_list[[i]]$y)
            if (max(eval_function_list[[i]]$y) > y_limit[2]) y_limit[2] = max(eval_function_list[[i]]$y)

          }

      } else {

        # Parameters cannot be estimated
        eval_function_list[[i]] = list(x = rep(NA, n_evaluation_points), y = rep(NA, n_evaluation_points)) 

      }
    }

    bar_width = (x_limit[2] - x_limit[1]) / 100

    for (i in 1:n_selected_models) {

        if (i < n_selected_models) page_break = TRUE else page_break = FALSE

        current_model = model_fit[[i]]

        filename = paste0("model_fit", i, ".emf")

        xlab = "Dose"    
        if (endpoint_index == 1) ylab = "Mean response"    
        if (endpoint_index == 2) ylab = "Response rate"    
        if (endpoint_index == 3) ylab = "Average number of events"    

        test_statistic = mcp_results$test_statistics[i]
        target_dose = current_model$target_dose       
        if (mcp_results$sign_model[i] == 1) {
          if (current_model$status >= 0) sign = "This model is included in the set of significant models" else sign = "This model is included in the set of significant models but the parameters cannot be estimated" 
          if (current_model$target_dose >= 0) target_dose = paste0("the target dose is ", round(current_model$target_dose, 3), ".") else target_dose = "the target dose cannot be determined."
        } else {
          sign = "This model is not included in the set of significant models"
          target_dose = "the target dose is undefined."
        }

        emf(file = filename, width = width, height = height)
        plot(x = descriptive_statistics$dose_levels, y = descriptive_statistics$mean_group, xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="p", pch = 19) 
        lines(x = eval_function_list[[i]]$x, y = eval_function_list[[i]]$y, col = "black", lwd = 2)
        for (j in 1:n_doses) {
          lines(x = rep(descriptive_statistics$dose_levels[j], 2), y = c(descriptive_statistics$lower_cl[j], descriptive_statistics$upper_cl[j]), col = "black", lwd = 1)
          lines(x = c(descriptive_statistics$dose_levels[j] - bar_width, descriptive_statistics$dose_levels[j] + bar_width), y = rep(descriptive_statistics$lower_cl[j], 2), col = "black", lwd = 1)
          lines(x = c(descriptive_statistics$dose_levels[j] - bar_width, descriptive_statistics$dose_levels[j] + bar_width), y = rep(descriptive_statistics$upper_cl[j], 2), col = "black", lwd = 1)
        }
        if (mcp_results$sign_model[i] == 1 & current_model$target_dose >= x_limit[1] & current_model$target_dose <= x_limit[2]) abline(v = current_model$target_dose, lty = "dashed") 
        dev.off()

        item_list[[item_index]] =  list(label = paste0("Figure ", figure_index, ". ", DF_selected_model_list[i], " model."), 
                                 filename = filename,
                                 dim = c(width, height),
                                 type = "emf_plot",
                                 footnote = paste0(sign, ", ", target_dose),
                                 page_break = page_break)

        item_index = item_index + 1
        figure_index = figure_index + 1

    }

  #############################################################################

  report = item_list

  doc = SaveReport(report, report_title)

  return(doc)

}
# End of GenerateAnalysisReport      

GenerateSimulationReport = function(results, report_title) { 

    #############################################################################

    # Error checks

    if (class(results) != "MCPModSimulationResults") stop("SimulationReport: The object was not created by the MCPModSimulation function.", call. = FALSE)

    if (!requireNamespace("officer", quietly = TRUE)) {
     stop("SimulationReport: The officer package is required to generate this report. Please install it.", call. = FALSE)
    }

    if (!requireNamespace("flextable", quietly = TRUE)) {
     stop("SimulationReport: The flextable package is required to generate this report. Please install it.", call. = FALSE)
    }

    if (!requireNamespace("devEMF", quietly = TRUE)) {
     stop("SimulationReport: The devEMF package is required to generate this report. Please install it.", call. = FALSE)
    }

    #############################################################################

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract the user-specified values of the model parameters
    user_specified = input_parameters$user_specified

    # Extract the initial values of the model parameters
    initial_values = input_parameters$initial_values

    # Simulation results
    sim_results = results$sim_results

    # Total number of models
    n_models = length(DF_model_list)

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    DF_selected_model_list = DF_model_list[selected_models]

    endpoint_index = input_parameters$endpoint_index

    dose_levels = input_parameters$sim_parameters$doses

    n_doses = length(dose_levels)

    digits = 3

    #############################################################################

    # Empty list of tables to be included in the report
    item_list = list()
    item_index = 1
    table_index = 1
    figure_index = 1

    #############################################################################

    column_names = c("Parameter", "Value")

    col1 = NULL
    col2 = NULL

    if (input_parameters$direction_index == 1) direction_label = "Increasing" else direction_label = "Decreasing"

    if (endpoint_index != 3) {

      col1 = c(col1, "Endpoint type", "Candidate models", "One-sided Type I error rate", "Direction of the dose-response relationship")
      col2 = c(col2, DF_endpoint_list[endpoint_index], 
                     paste0(DF_selected_model_list, collapse = ", "), 
                     input_parameters$alpha,
                     direction_label)

    } else {

      col1 = c(col1, "Endpoint type", "Candidate models", "One-sided Type I error rate", "Direction of the dose-response relationship", "Overdispersion parameters")
      col2 = c(col2, DF_endpoint_list[endpoint_index], 
                     paste0(DF_selected_model_list, collapse = ", "), 
                     input_parameters$alpha,
                     direction_label, 
                     paste0(round(input_parameters$theta, digits), collapse = ", "))

    }

    data_frame = cbind(col1, col2)
    title = paste0("Table ", table_index, ". General parameters.")

    column_width = c(2.5, 4)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    column_names = c("Candidate model", "Parameter values")

    col1 = NULL
    col2 = NULL

    col1 = DF_selected_model_list
    for (i in 1:n_models) {
     
        if (selected_models[i] == TRUE) {

          n_parameters = length(DF_model_parameters[[i]])
          model_user_specified = round(user_specified[[i]][1:n_parameters], digits)
          for (j in 1:n_parameters) model_user_specified[j] = paste0(DF_model_parameters[[i]][j], " = ", model_user_specified[j])
          col2 = c(col2, paste0(model_user_specified, collapse = ", "))

        }    
    }

    data_frame = cbind(col1, col2)
    title = paste0("Table ", table_index, ". Initial values of the model parameters.")

    column_width = c(2.5, 4)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    column_names = c("Parameter", "Value")

    col1 = NULL
    col2 = NULL

    if (input_parameters$model_selection == "AIC") model_selection_label = "Select the model with the smallest AIC"
    if (input_parameters$model_selection == "maxT") model_selection_label = "Select the model corresponding to the largest test statistic"
    if (input_parameters$model_selection == "aveAIC") model_selection_label = "Average the models using AIC-based weights"

    col1 = c(col1, "Model selection criterion", "Delta")
    col2 = c(col2, model_selection_label, input_parameters$delta)

    footnote = "Delta is the pre-defined clinically meaningful improvement over placebo."

    data_frame = cbind(col1, col2)
    title = paste0("Table ", table_index, ". Model selection parameters.")

    column_width = c(3, 3.5)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    if (endpoint_index == 1) {

      column_names = c("Dose", "Number of patients", "Standard deviation")  
      data_frame = cbind(dose_levels, 
                         sprintf("%d", input_parameters$sim_parameters$n),
                         input_parameters$sim_model_list$sd)
      column_width = c(2.5, 2, 2)

    } else {

      column_names = c("Dose", "Number of patients")  
      data_frame = cbind(dose_levels, 
                         sprintf("%d", input_parameters$sim_parameters$n))
      column_width = c(2.5, 4)

    } 

    title = paste0("Table ", table_index, ". Trial parameters.")

    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    column_names = c("Parameter", "Value")

    col1 = c("Patient dropout rate", "Number of simulation runs")
    col2 = c(input_parameters$sim_parameters$dropout_rate,
             input_parameters$sim_parameters$nsims)

    data_frame = cbind(col1, col2)

    title = paste0("Table ", table_index, ". Simulation parameters.")

    column_width = c(3.5, 3)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1 

    #############################################################################

    column_names = c("Scenario", "Placebo effect", "Maximum effect over placebo", "Simulation model parameters")

    n_scenarios = length(input_parameters$max_effect)
    scenario_list = 1:n_scenarios

    model_index = input_parameters$sim_model_list$sim_model_index
    n_parameters = length(DF_model_parameters[[model_index]])

    model_parameters = rep("", n_scenarios)
    temp_string = rep("", n_parameters)

    for (i in 1:n_scenarios) {

        for (j in 1:n_parameters) temp_string[j] = paste0(DF_model_parameters[[model_index]][j], " = ", round(input_parameters$sim_model_list$sim_parameter_values[i, j], digits))

        model_parameters[i] = paste0(temp_string, collapse = ", ")

    }

    data_frame = cbind(sprintf("%d", scenario_list), 
                       rep(round(input_parameters$placebo_effect, digits), n_scenarios),
                       round(input_parameters$max_effect, digits),
                       model_parameters)

    column_width = c(1, 1, 1.5, 3)
    title = paste0("Table ", table_index, ". Assumed dose-response scenarios (", DF_model_list[model_index], " model).")

    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE)
    item_index = item_index + 1
    table_index = table_index + 1  

    #############################################################################

    # Plot simulation models

    width = 6.5
    height = 5

    true_target_dose = round(sim_results$true_target_dose, digits)
    true_target_dose[true_target_dose == -1] = NA

    eval_function_list = list()  

    # Determine the axis ranges
    x_limit = c(min(dose_levels), max(dose_levels))
    temp = c(input_parameters$placebo_effect, input_parameters$placebo_effect + input_parameters$max_effect)
    y_limit = c(min(temp), max(temp))

    xlab = "Dose"    
    if (endpoint_index == 1) {
      ylab = "Mean response"    
    }
    if (endpoint_index == 2) {
      ylab = "Response rate"    
      y_limit = c(0, 1)
    }
    if (endpoint_index == 3) {
      ylab = "Average number of events"    
    }

    for (i in 1:n_scenarios) {

        filename = paste0("assumed_model", i, ".emf")

        eval_function_list[[i]] = EvaluateDRFunction(model_index, endpoint_index, input_parameters$sim_model_list$sim_parameter_values[i, ], dose_levels) 

        emf(file = filename, width = width, height = height)
        plot(x = eval_function_list[[i]]$x, y = eval_function_list[[i]]$y, xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="l", lwd = 2) 
        if (!is.na(true_target_dose[i])) abline(v = true_target_dose[i], lty = "dashed") 
        dev.off()

        if (!is.na(true_target_dose[i])) target_dose = paste0("The true target dose is ", true_target_dose[i], ".") else target_dose = "The true target dose is undefined."

        item_list[[item_index]] =  list(label = paste0("Figure ", figure_index, ". Assumed dose-response model (Scenario ", i, ")."), 
                                 filename = filename,
                                 dim = c(width, height),
                                 type = "emf_plot",
                                 footnote = target_dose,
                                 page_break = TRUE)

        item_index = item_index + 1
        figure_index = figure_index + 1

    }

    #############################################################################

    column_names = c("Scenario", "Maximum effect over placebo", "Power")

    data_frame = cbind(sprintf("%d", scenario_list),
                       round(input_parameters$max_effect, digits), 
                       round(sim_results$power, digits))

    title = paste0("Table ", table_index, ". Simulation results: Power.")

    footnote = "Power is the probability that the best dose-response contrast is significant."

    column_width = c(2, 2.5, 2)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1 

    #############################################################################

    # Go probabilities unless the model selection method is aveAIC
    if (input_parameters$model_selection != "aveAIC") {

      column_names = c("Scenario", "Maximum effect over placebo", "Go probability")

      data_frame = cbind(sprintf("%d", scenario_list),
                         round(input_parameters$max_effect, digits), 
                         round(sim_results$go_prob, digits))

      title = paste0("Table ", table_index, ". Simulation results: Go probabilities based on the threshold of ", input_parameters$go_threshold, ".")

      footnote = "The go probability is the probability that the best dose-response contrast is significant and the maximum effect for the corresponding model exceeds the pre-defined go threshold."

      column_width = c(2, 2.5, 2)
      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE, footnote)
      item_index = item_index + 1
      table_index = table_index + 1 

    }

    #############################################################################

    if (input_parameters$model_selection != "aveAIC") {

      column_names = c("Scenario", "No model selected", DF_model_list_short[selected_models])

      data_frame = cbind(sprintf("%d", scenario_list), 
                         round(sim_results$model_index_summary, digits))

      title = paste0("Table ", table_index, ". Simulation results: Probability of selecting a dose-response model.")

      original_width = c(0.5, 0.8, 0.8, 1.1, 0.7, 0.8, 0.8) 
      column_width = c(1, 1, rep(4.5 / n_selected_models, n_selected_models))

      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
      item_index = item_index + 1
      table_index = table_index + 1 

    }

    #############################################################################

    column_names = c("Scenario", "True target dose", "Mean target dose (95% CI)")

    true_target_dose = round(sim_results$true_target_dose, digits)
    true_target_dose[true_target_dose == -1] = "NA"
    true_target_dose = as.character(true_target_dose)

    lower_bound_target_dose = round(sim_results$target_dose_summary[, 1], digits)
    lower_bound_target_dose[lower_bound_target_dose == -1] = NA

    mean_target_dose = round(sim_results$target_dose_summary[, 2], digits)
    mean_target_dose[mean_target_dose == -1] = NA

    upper_bound_target_dose = round(sim_results$target_dose_summary[, 3], digits)
    upper_bound_target_dose[upper_bound_target_dose == -1] = NA

    data_frame = cbind(sprintf("%d", scenario_list), 
                       true_target_dose,
                       paste0(mean_target_dose, " (", lower_bound_target_dose, ", ", upper_bound_target_dose, ")"))

    title = paste0("Table ", table_index, ". Simulation results: Target dose estimates.")

    column_width = c(2, 2, 2.5) 
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1 

    #############################################################################

    column_names = c("Scenario", "No dose found", dose_levels[2:n_doses], "Greater than max dose")

    data_frame = cbind(sprintf("%d", scenario_list), 
                       round(sim_results$target_dose_categorical_summary, digits))

    title = paste0("Table ", table_index, ". Simulation results: Probability of identifying the target dose.")

    footnote = "Each column presents the probability that the estimated target dose is less than or equal to the current dose and is strictly greater than the next lower dose."

    column_width = c(1, rep(5.5 / (n_doses + 1), n_doses + 1)) 
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1 

    #############################################################################

    # Plot power

    width = 6.5
    height = 5

    # Determine the axis ranges
    x_limit = c(min(input_parameters$max_effect), max(input_parameters$max_effect))
    y_limit = c(0, 1)

    filename = paste0("power.emf")

    xlab = "Maximum effect over placebo"    
    ylab = "Power"    

    emf(file = filename, width = width, height = height)
    plot(x = input_parameters$max_effect, y = sim_results$power, xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="l", lwd = 2) 
    dev.off()

    item_list[[item_index]] =  list(label = paste0("Figure ", figure_index, ". Simulation results: Power."), 
                             filename = filename,
                             dim = c(width, height),
                             type = "emf_plot",
                             footnote = "",
                             page_break = FALSE)

    item_index = item_index + 1
    figure_index = figure_index + 1

    #############################################################################

    # Plot go probabilities  unless the model selection method is aveAIC
    if (input_parameters$model_selection != "aveAIC") {

      width = 6.5
      height = 5

      # Determine the axis ranges
      x_limit = c(min(input_parameters$max_effect), max(input_parameters$max_effect))
      y_limit = c(0, 1)

      page_break = FALSE

      filename = paste0("go_prob.emf")

      xlab = "Maximum effect over placebo"    
      ylab = "Probability"    
       
      emf(file = filename, width = width, height = height)
      plot(x = input_parameters$max_effect, y = sim_results$go_prob, xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="l", lwd = 2) 
      dev.off()

      footnote = ""

      item_list[[item_index]] =  list(label = paste0("Figure ", figure_index, ". Simulation results: Go probabilities based on the threshold of ", input_parameters$go_threshold, "."), 
                               filename = filename,
                               dim = c(width, height),
                               type = "emf_plot",
                               footnote = footnote,
                               page_break = page_break)

      item_index = item_index + 1
      figure_index = figure_index + 1

    }


    #############################################################################

    # Plot estimated dose-response curves unless the model selection method is aveAIC
    if (input_parameters$model_selection != "aveAIC") {

      width = 6.5
      height = 5

      # Determine the axis ranges
      x_limit = c(min(dose_levels), max(dose_levels))
      temp = c(input_parameters$placebo_effect, input_parameters$placebo_effect + input_parameters$max_effect, sim_results$dose_response_lower, sim_results$dose_response_upper)    
      y_limit = c(min(temp, na.rm = TRUE), max(temp, na.rm = TRUE))

      if (endpoint_index == 2) {
        y_limit[1] = min(0, y_limit[1])
        y_limit[2] = max(1, y_limit[2])
      }

      for (i in 1:n_scenarios) {

          if (i < n_scenarios) page_break = TRUE else page_break = FALSE

          filename = paste0("dose_response_model", i, ".emf")

          xlab = "Dose"    
          if (endpoint_index == 1) ylab = "Mean response"    
          if (endpoint_index == 2) ylab = "Response rate"    
          if (endpoint_index == 3) ylab = "Average number of events"    
       
          emf(file = filename, width = width, height = height)
          plot(x = sim_results$dosex, y = sim_results$dose_response_mean[, i], xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="l", lwd = 2) 
          polygon(c(rev(sim_results$dosex), sim_results$dosex), c(rev(sim_results$dose_response_upper[, i]), sim_results$dose_response_lower[, i]), col = "grey80", border = NA)
          lines(x = sim_results$dosex, y = sim_results$dose_response_mean[, i], col="black", lwd = 2)  
          lines(x = eval_function_list[[i]]$x, y = eval_function_list[[i]]$y, col="red", lwd = 2)  
          
          if (!is.na(sim_results$target_dose_summary[i, 2])) abline(v = sim_results$target_dose_summary[i, 2], lty = "dashed") 
          dev.off()

          if (!is.na(sim_results$target_dose_summary[i, 2])) footnote = paste0("Red curve: Assumed dose-response curve. Black curve: Estimated dose-response curve with a 95% confidence band. The estimated target dose is ", round(sim_results$target_dose_summary[i, 2], digits), ".") else footnote = "Red curve: Assumed dose-response curve. Black curve: Estimated dose-response curve with a 95% confidence band. The estimated target dose is undefined."

          item_list[[item_index]] =  list(label = paste0("Figure ", figure_index, ". Simulation results: Assumed and estimated dose-response curve based on the best model selected by MCPMod (Scenario ", i, ", Maximum effect over placebo is ", input_parameters$max_effect[i], ")."), 
                                   filename = filename,
                                   dim = c(width, height),
                                   type = "emf_plot",
                                   footnote = footnote,
                                   page_break = page_break)

          item_index = item_index + 1
          figure_index = figure_index + 1
          
      }

    }

    #############################################################################

    report = item_list

    doc = SaveReport(report, report_title)

    return(doc)


}
# End of GenerateSimulationReport

AnalysisReport = function(results, report_title, report_filename) { 

    # Generate the report
    doc = GenerateAnalysisReport(results, report_title)

    # Save the report
    print(doc, target = report_filename)          


}

SimulationReport = function(results, report_title, report_filename) { 

    # Generate the report
    doc = GenerateSimulationReport(results, report_title)

    # Save the report
    print(doc, target = report_filename)          


}