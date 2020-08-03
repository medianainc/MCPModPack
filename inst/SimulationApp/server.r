DF_endpoint_list = c("Normal", "Binary", "Count")
DF_model_list = c("Linear", "Quadratic", "Exponential", "Emax", "Logistic", "SigEmax")
DF_model_list_short = c("Lin", "Quad", "Exp", "Emax", "Logist", "SigEmax")
DF_model_parameters = list(c("e0", "delta"), c("e0", "delta1", "delta2"), c("e0", "e1", "delta"), c("e0", "eMax", "ed50"), c("e0", "eMax", "ed50", "delta"), c("e0", "eMax", "ed50", "h"))
n_models = length(DF_model_list)
digits = 3

Logit = function(x) {
    return(log(x /(1 - x)))
}

shinyServer(function(input, output, session) {

  output$AssumedDoseResponse = renderPlot({
    # Extract input parameters
    select_scenario_model = as.numeric(input$select_scenario_model)
    endpoint_index = as.numeric(input$endpoint_index)
    dose_levels = as.numeric(unlist(strsplit(input$sim_parameters_doses, ",")))
    placebo_effect = input$sim_models_placebo_effect
    placebo_effect_temp = placebo_effect
    max_effect = as.numeric(unlist(strsplit(input$sim_models_max_effect, ",")))
    max_effect_temp = max_effect
    model_index = as.numeric(input$sim_model_index)
    n_scenarios = length(max_effect)
    direction_index = as.numeric(input$direction)
    max_dose = max(dose_levels)
    
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

    # Compute parameters of the assumed dose-response model to match the placebo and maximum effects 
    sim_model = list(linear = NA)
    if (as.numeric(input$sim_model_index) == 2) {
      sim_model = list(quadratic = NA)
    } else
    if (as.numeric(input$sim_model_index) == 3) {
      sim_model = list(exponential = input$sim_model_exponential_delta)
    } else
    if (as.numeric(input$sim_model_index) == 4) {
      sim_model = list(emax = input$sim_model_emax_ed50)
    } else
    if (as.numeric(input$sim_model_index) == 5) {
      sim_model = list(logistic = c(input$sim_model_logistic_ed50, input$sim_model_logistic_delta))
    } else
    if (as.numeric(input$sim_model_index) == 6) {
      sim_model = list(sigemax = c(input$sim_model_sigemax_ed50, input$sim_model_sigemax_h))
    }

    if (model_index == 1) { 
      # linear
      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 2)
      parameters = 0
      for (i in 1:n_scenarios) {
          coef = ComputeDRFunctionParameters(model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:2) sim_parameter_values[i, j] = coef[j]
      }
    }

    if (model_index == 2) {
      # quadratic
      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 3)
      parameters = sim_model$quadratic

      for (i in 1:n_scenarios) {
          coef = ComputeDRFunctionParameters(model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:3) sim_parameter_values[i, j] = coef[j]
      }
    }

    if (model_index == 3) {
      # exponential
      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 3)
      parameters = sim_model$exponential

      for (i in 1:n_scenarios) {
          coef = ComputeDRFunctionParameters(model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:3) sim_parameter_values[i, j] = coef[j]
      }
    }

    if (model_index == 4) {
      # emax
      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 3)
      parameters = sim_model$emax

      for (i in 1:n_scenarios) {
          coef = ComputeDRFunctionParameters(model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:3) sim_parameter_values[i, j] = coef[j]
      }
    }

    if (model_index == 5) {
      # logistic
      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 4)
      parameters = sim_model$logistic

      for (i in 1:n_scenarios) {
          coef = ComputeDRFunctionParameters(model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:4) sim_parameter_values[i, j] = coef[j]
      }
    }

    if (model_index == 6) {
      # sigemax
      # Create a matrix with possible values of model parameters
      sim_parameter_values = matrix(0, n_scenarios, 4)
      parameters = sim_model$sigemax

      for (i in 1:n_scenarios) {
          coef = ComputeDRFunctionParameters(model_index, placebo_effect_temp, max_effect_temp[i], max_dose, parameters)
          for (j in 1:4) sim_parameter_values[i, j] = coef[j]
      }
    }

    # Determine the axis ranges
    x_limit = c(min(dose_levels), max(dose_levels))
    temp = c(placebo_effect, placebo_effect + max_effect)
    y_limit = c(min(temp), max(temp))

    title = paste0("Assumed dose-response model (Scenario ", select_scenario_model, ")")
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

    eval_function = EvaluateDRFunction(model_index, endpoint_index, sim_parameter_values[select_scenario_model, ], dose_levels) 
    plot(x = eval_function$x, y = eval_function$y, xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="l", lwd = 2, main = title)
  })

  RunMCPModSimulation = reactive({
    # List of candidate dose-response models
    models = list()

    if (input$linear_model) models[["linear"]] = NA
    if (input$quadratic_model) models[["quadratic"]] = as.numeric(input$quadratic_model_par1)
    if (input$exponential_model) models[["exponential"]] = as.numeric(input$exponential_model_par1)
    if (input$emax_model) models[["emax"]] = as.numeric(input$emax_model_par1)
    if (input$logistic_model) models[["logistic"]] = c(as.numeric(input$logistic_model_par1), as.numeric(input$logistic_model_par2)) 
    if (input$sigemax_model) models[["sigemax"]] = c(as.numeric(input$sigemax_model_par1), as.numeric(input$sigemax_model_par2)) 

    # One-sided Type I error rate
    alpha = as.numeric(input$alpha)

    # Endpoint type
    if (as.numeric(input$endpoint_index) == 1) endpoint_type = "Normal"
    if (as.numeric(input$endpoint_index) == 2) endpoint_type = "Binary"
    if (as.numeric(input$endpoint_index) == 3) endpoint_type = "Count"

    # Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
    if (as.numeric(input$direction) == 1) direction = "increasing"
    if (as.numeric(input$direction) == -1) direction = "decreasing"

    # Model selection criterion
    if (as.numeric(input$model_selection) == 1) model_selection = "AIC"
    if (as.numeric(input$model_selection) == 2) model_selection = "maxT"
    if (as.numeric(input$model_selection) == 3) model_selection = "aveAIC"

    # The treatment effect for identifying the target dose 
    Delta = as.numeric(input$delta)

    # Vector of overdispersion parameters
    if (endpoint_type == "Count") theta = as.numeric(unlist(strsplit(input$theta, ",")))
    else theta = 0

    # Select the assumed dose-response model, values of the non-linear model parameters and the standard deviation of the outcome variable in each trial arm (required for normally distributed endpoint)
    sim_model = list(linear = NA)
    if (as.numeric(input$sim_model_index) == 2) {
      sim_model = list(quadratic = NA)
    } else
    if (as.numeric(input$sim_model_index) == 3) {
      sim_model = list(exponential = input$sim_model_exponential_delta)
    } else
    if (as.numeric(input$sim_model_index) == 4) {
      sim_model = list(emax = input$sim_model_emax_ed50)
    } else
    if (as.numeric(input$sim_model_index) == 5) {
      sim_model = list(logistic = c(input$sim_model_logistic_ed50, input$sim_model_logistic_delta))
    } else
    if (as.numeric(input$sim_model_index) == 6) {
      sim_model = list(sigemax = c(input$sim_model_sigemax_ed50, input$sim_model_sigemax_h))
    }

    placebo_effect = input$sim_models_placebo_effect
    max_effect = as.numeric(unlist(strsplit(input$sim_models_max_effect, ",")))
    sd = as.numeric(unlist(strsplit(input$sim_models_sd, ",")))

    sim_models = c(
                  sim_model,
                  list( 
                    placebo_effect = placebo_effect, 
                    max_effect = max_effect, 
                    sd = sd))
    
    npatients = as.numeric(unlist(strsplit(input$sim_parameters_n, ",")))
    doses = as.numeric(unlist(strsplit(input$sim_parameters_doses, ",")))
    dropout_rate = as.numeric(input$sim_parameters_dropout_rate)
    nsims = as.numeric(input$sim_parameters_nsims)
    go_threshold = as.numeric(input$sim_parameters_go_threshold)

    # Simulation parameters (number of patients in each trial arm, dose levels, dropout rate and number of simulations)
    sim_parameters = list(n = npatients,
                        doses = doses,
                        dropout_rate = dropout_rate,
                        nsims = nsims,
                        go_threshold = go_threshold)

    withProgress(message = "Running simulations", value = 1, {

      # Perform an MCPMod-based simulation of the trial's data
      results = MCPModSimulation(endpoint_type = endpoint_type, 
                              models = models, 
                              alpha = alpha, 
                              direction = direction, 
                              model_selection = model_selection, 
                              Delta = Delta,
                              theta = theta,
                              sim_models = sim_models,
                              sim_parameters = sim_parameters)
    })

    # Return the list of results
    results
  })

  EstimatedDoseResponseModels = reactive({
    results = RunMCPModSimulation()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    DF_selected_model_list = DF_model_list[selected_models]
    DF_selected_model_index_list = (1:n_models)[selected_models]

    # Simulation results
    sim_results = results$sim_results

    # Extract the Mod step results
    mod_results = results$mod_results

    # Extract the list of model fit parameters
    model_fit = mod_results$model_fit

    endpoint_index = input_parameters$endpoint_index

    dose_levels = input_parameters$sim_parameters$doses

    model_index = input_parameters$sim_model_list$sim_model_index

    # Evaluate the dose-response function under the selected scenario
    i = as.numeric(input$select_scenario)

    eval_function_list = EvaluateDRFunction(model_index, endpoint_index, input_parameters$sim_model_list$sim_parameter_values[i, ], dose_levels) 

    result_list = list(eval_function_list = eval_function_list,
                              dose_response_mean = sim_results$dose_response_mean[, i],
                              dose_response_lower = sim_results$dose_response_lower[, i],
                              dose_response_upper = sim_results$dose_response_upper[, i])

    # Return the list of results
    result_list

  })    

  output$PowerSummary = renderTable({
    results = RunMCPModSimulation()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract simulation results
    sim_results = results$sim_results

    column_names = c("Maximum effect over placebo", "Power")
    data_frame = cbind(input_parameters$max_effect, 
                   round(sim_results$power, digits))
    colnames(data_frame) = column_names

    data_frame
  }, digits = digits)

  output$ProbabilityOfModel = renderTable({
    results = RunMCPModSimulation()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract simulation results
    sim_results = results$sim_results

    # Selected models
    selected_models = results$selected_models

    if (input_parameters$model_selection != "aveAIC") {

      data_frame = cbind(input_parameters$max_effect, 
                    round(sim_results$model_index_summary, digits))
      colnames(data_frame) = c("Maximum effect over placebo", "No model selected", DF_model_list[selected_models])

      data_frame
    }

  }, digits = digits)

  output$TargetDoseEstimates = renderTable({
    results = RunMCPModSimulation()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract simulation results
    sim_results = results$sim_results

    true_target_dose = round(sim_results$true_target_dose, digits)
    true_target_dose[true_target_dose == -1] = "NA"
    true_target_dose = as.character(true_target_dose)

    lower_bound_target_dose = round(sim_results$target_dose_summary[, 1], digits)
    lower_bound_target_dose[lower_bound_target_dose == -1] = NA

    mean_target_dose = round(sim_results$target_dose_summary[, 2], digits)
    mean_target_dose[mean_target_dose == -1] = NA

    upper_bound_target_dose = round(sim_results$target_dose_summary[, 3], digits)
    upper_bound_target_dose[upper_bound_target_dose == -1] = NA

    data_frame = cbind(input_parameters$max_effect, 
                   true_target_dose,
                   paste0(mean_target_dose, " (", lower_bound_target_dose, ", ", upper_bound_target_dose, ")"))
    colnames(data_frame) = c("Maximum effect over placebo", "True target dose", "Mean target dose (95% CI)")

    data_frame
  })

  output$GoProbabilities = renderTable({
    results = RunMCPModSimulation()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract simulation results
    sim_results = results$sim_results

    n_scenarios = length(input_parameters$max_effect)
    scenario_list = 1:n_scenarios

    data_frame = cbind(round(input_parameters$max_effect, digits), 
                        round(sim_results$go_prob, digits))
    colnames(data_frame) = c("Maximum effect over placebo", "Go probability")

    data_frame
  })

  output$ProbabilityOfSelectingADose = renderTable({
    results = RunMCPModSimulation()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract simulation results
    sim_results = results$sim_results

    dose_levels = input_parameters$sim_parameters$doses
    n_doses = length(dose_levels)

    data_frame = cbind(input_parameters$max_effect, 
                   round(sim_results$target_dose_categorical_summary, digits))

    colnames(data_frame) = c("Maximum effect over placebo", "No dose selected", dose_levels[2:n_doses], "Greater than max dose")

    data_frame
  }, digits = digits)

  output$DownloadResults = downloadHandler(
    filename = function() {
      "Report.docx"
    },

    content = function(file) {
      results = RunMCPModSimulation()

      doc = GenerateSimulationReport(results, "MCPMod simulation report")

      # Save the report
      xfile = paste0(file, ".docx")
      print(doc, target = xfile)          
      file.rename(xfile, file)
    }
  )

  output$Power = renderPlot({
    results = RunMCPModSimulation()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract simulation results
    sim_results = results$sim_results

    # Determine the axis ranges
    x_limit = c(min(input_parameters$max_effect), max(input_parameters$max_effect))
    y_limit = c(0, 1)

    xlab = "Maximum effect over placebo"    
    ylab = "Power"    

    plot(x = input_parameters$max_effect, y = sim_results$power, xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="l", lwd = 2) 
  })

  output$GoProbability = renderPlot({
    results = RunMCPModSimulation()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract simulation results
    sim_results = results$sim_results

    # Determine the axis ranges
    x_limit = c(min(input_parameters$max_effect), max(input_parameters$max_effect))
    y_limit = c(0, 1)

    xlab = "Maximum effect over placebo"    
    ylab = "Probability"    

    plot(x = input_parameters$max_effect, y = sim_results$go_prob, xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="l", lwd = 2) 
  })

  output$EstimatedDoseResponse = renderPlot({
    results = RunMCPModSimulation()

    result_list = EstimatedDoseResponseModels()

    # Extract input parameters
    input_parameters = results$input_parameters

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

    # Determine the axis ranges
    x_limit = c(min(dose_levels), max(dose_levels))
    temp = c(input_parameters$placebo_effect, input_parameters$placebo_effect + input_parameters$max_effect, sim_results$dose_response_lower, sim_results$dose_response_upper)    
    y_limit = c(min(temp, na.rm = TRUE), max(temp, na.rm = TRUE))

    if (endpoint_index == 2) {
      y_limit[1] = min(0, y_limit[1])
      y_limit[2] = max(1, y_limit[2])
    }

    xlab = "Dose"    
    if (endpoint_index == 1) ylab = "Mean response"    
    if (endpoint_index == 2) ylab = "Response rate"    
    if (endpoint_index == 3) ylab = "Average number of events"    
  
    plot(x = sim_results$dosex, y = result_list$dose_response_mean, xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="l", lwd = 2) 
    polygon(x = c(rev(sim_results$dosex), sim_results$dosex), y = c(rev(result_list$dose_response_upper), result_list$dose_response_lower), col = "grey80", border = NA)
    lines(x = sim_results$dosex, y = result_list$dose_response_mean, col="black", lwd = 2)  
    lines(x = result_list$eval_function_list$x, y = result_list$eval_function_list$y, col="red", lwd = 2)  
  })

  observeEvent(input$sim_models_max_effect, {
    max_effect = as.numeric(unlist(strsplit(input$sim_models_max_effect, ",")))

    # Update UI radio buttons with plot variants
    max_effect_variants = setNames(1:length(max_effect), max_effect)
    updateRadioButtons(session, inputId = "select_scenario",
      choices = max_effect_variants,
      selected = min(max_effect_variants),
      inline = TRUE
    )

    updateRadioButtons(session, inputId = "select_scenario_model",
      choices = max_effect_variants,
      selected = min(max_effect_variants),
      inline = TRUE
    )
  })

  observeEvent(input$jump_to_panel2, {
        updateTabItems(session, "sidebarMCPModSimulation",
                          selected = "responseModels")
  })

  observeEvent(input$jump_to_panel3, {
        updateTabItems(session, "sidebarMCPModSimulation",
                          selected = "simulation")
  })

  observeEvent(input$jump_to_panel4, {
        updateTabItems(session, "sidebarMCPModSimulation",
                          selected = "report")
  })

})
