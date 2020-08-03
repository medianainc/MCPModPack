library(shiny)

DF_endpoint_list = c("Normal", "Binary", "Count")
DF_model_list = c("Linear", "Quadratic", "Exponential", "Emax", "Logistic", "SigEmax")
DF_model_list_short = c("Lin", "Quad", "Exp", "Emax", "Logist", "SigEmax")
DF_model_parameters = list(c("e0", "delta"), c("e0", "delta1", "delta2"), c("e0", "e1", "delta"), c("e0", "eMax", "ed50"), c("e0", "eMax", "ed50", "delta"), c("e0", "eMax", "ed50", "h"))
n_models = length(DF_model_list)
digits = 3

shinyServer(function(input, output, session) {


  GetDataSet = reactive({

      validate(
        need(input$data_set, message = "Error: Please load the trial's data set.")
      )

      data_set_object = input$data_set

      if (is.null(data_set_object)) return(NULL)

      tryCatch({
          data_set = read.csv(data_set_object$datapath, header = TRUE, na = ".", stringsAsFactors = FALSE)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
              
      return(data_set)
  })

  output$error_message <- renderText({ 
      
      # Open the trial data set  
      validate(
        need(input$data_set, message = "Error: Please load the trial's data set.")
      )

      data_set = GetDataSet()  

      # Check that the dose and resp variables are defined
      validate(
        need(!is.null(data_set$dose), message = "Error: The trial's data set does not include the dose variable.")
      )

      validate(
        need(!is.null(data_set$resp), message = "Error: The trial's data set does not include the resp variable.")
      )

      n_doses = length(unique(data_set$dose))

      # Endpoint type
      if (as.numeric(input$endpoint_index) == 1) endpoint_type = "Normal"
      if (as.numeric(input$endpoint_index) == 2) endpoint_type = "Binary"
      if (as.numeric(input$endpoint_index) == 3) endpoint_type = "Count"

      # Vector of overdispersion parameters
      if (endpoint_type == "Count") {
        theta = as.numeric(unlist(strsplit(input$theta, ",")))

        validate(
          need(length(theta) == n_doses, message = "Error: The number of overdispersion parameters must be equal to the number of dose levels in the trial's data set.")
        )

      }



  })


  RunMCPModAnalysis = reactive({
    # Open the trial data set  
    data_set = GetDataSet()  

    # List of candidate dose-response models
    models = list()

    if (input$linear_model) models[["linear"]] = NA
    if (input$quadratic_model) models[["quadratic"]] = as.numeric(input$quadratic_model_par1)
    if (input$exponential_model) models[["exponential"]] = as.numeric(input$exponential_model_par1)
    if (input$emax_model) models[["emax"]] = as.numeric(input$emax_model_par1)
    if (input$logistic_model) models[["logistic"]] = c(as.numeric(input$logistic_model_par1), as.numeric(input$logistic_model_par2)) 
    if (input$sigemax_model) models[["sigemax"]] = c(as.numeric(input$sigemax_model_par1), as.numeric(input$sigemax_model_par2)) 

    # Update UI choises for plot
    current_models = c()
    if (input$linear_model) current_models[["Linear"]] = 1
    if (input$quadratic_model) current_models[["Quadratic"]] = 2
    if (input$exponential_model) current_models[["Exponential"]] = 3
    if (input$emax_model) current_models[["Emax"]] = 4
    if (input$logistic_model) current_models[["Logistic"]] = 5
    if (input$sigemax_model) current_models[["SigEmax"]] = 6
    current_models = setNames(1:length(current_models), labels(current_models))

    updateRadioButtons(session, inputId = "select_plot_model",
      choices = current_models,
      selected = min(current_models),
      inline = TRUE
    )

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
    if (endpoint_type == "Count") theta = as.numeric(unlist(strsplit(input$theta, ","))) else theta = 0

    withProgress(message = "Running analysis", value = 1, {

      # Perform an MCPMod-based analysis of the trial's data
      results = MCPModAnalysis(endpoint_type = endpoint_type, 
                              models = models, 
                              dose = data_set$dose, 
                              resp = data_set$resp, 
                              alpha = alpha, 
                              direction = direction, 
                              model_selection = model_selection, 
                              Delta = Delta,
                              theta = theta)

    })

    # Return the list of results
    results
  })

  EvaluateDoseResponseModels = reactive({
    results = RunMCPModAnalysis()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    DF_selected_model_list = DF_model_list[selected_models]
    DF_selected_model_index_list = (1:n_models)[selected_models]

    # Extract the Mod step results
    mod_results = results$mod_results

    # Extract the list of model fit parameters
    model_fit = mod_results$model_fit

    # Extract the dose levels and responses
    dose = mod_results$dose

    # Evaluate the selected dose-response function
    i = as.numeric(input$select_plot_model)
    if (i > length(model_fit)) i = 1

    current_model = model_fit[[i]]

    # Evaluate the current dose-response model
    model_index = DF_selected_model_index_list[i]
    coef = current_model$coef

    eval_function_list = EvaluateDRFunction(model_index, input_parameters$endpoint_index, coef, dose)

    # Return the list of results
    eval_function_list
  })    

  output$DescriptiveStatistics = renderTable({
    results = RunMCPModAnalysis()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Extract descriptive statistics
    descriptive_statistics = results$descriptive_statistics

    # Normal endpoint
    if (input_parameters$endpoint_index == 1) {
      column_names = c("Dose", "n", "Mean", "95% CI", "SD")
      data_frame = cbind(descriptive_statistics$dose_levels, 
                descriptive_statistics$n_groups, 
                round(descriptive_statistics$mean_group, digits),
                paste0("(", round(descriptive_statistics$lower_cl, digits), ", ", round(descriptive_statistics$upper_cl, digits), ")"),
                round(descriptive_statistics$sd, digits))
    }

    # Binary endpoint
    if (input_parameters$endpoint_index == 2) {
      column_names = c("Dose", "n", "Rate", "95% CI")
      data_frame = cbind(descriptive_statistics$dose_levels, 
                descriptive_statistics$n_groups, 
                round(descriptive_statistics$mean_group, digits),
                paste0("(", round(descriptive_statistics$lower_cl, digits), ", ", round(descriptive_statistics$upper_cl, digits), ")"))
    }

    # Count endpoint
    if (input_parameters$endpoint_index == 3) {
      column_names = c("Dose", "n", "Mean", "95% CI", "SD", "Theta")

      data_frame = cbind(descriptive_statistics$dose_levels, 
                descriptive_statistics$n_groups, 
                round(descriptive_statistics$mean_group, digits),
                paste0("(", round(descriptive_statistics$lower_cl, digits), ", ", round(descriptive_statistics$upper_cl, digits), ")"),
                round(descriptive_statistics$sd, digits),
                round(input_parameters$theta, digits))
    }

    colnames(data_frame) = column_names

    data_frame
  })

  output$Contrasts = renderTable({
    results = RunMCPModAnalysis()

    # Extract descriptive statistics
    descriptive_statistics = results$descriptive_statistics

    # Extract the optimal contrasts and contrast correlation matrix
    contrast_results = results$contrast_results

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    DF_selected_model_list = DF_model_list[selected_models]

    column_names = c("Dose", DF_selected_model_list)

    data_frame = cbind(descriptive_statistics$dose_levels, round(contrast_results$opt_contrast, digits))

    colnames(data_frame) = column_names

    data_frame
  })

  output$CorrelationMatrix = renderTable({
    results = RunMCPModAnalysis()

    # Extract descriptive statistics
    descriptive_statistics = results$descriptive_statistics

    # Extract the optimal contrasts and contrast correlation matrix
    contrast_results = results$contrast_results

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    DF_selected_model_list = DF_model_list[selected_models]

    column_names = c("Models", DF_selected_model_list)

    data_frame = cbind(DF_selected_model_list, round(contrast_results$corr_matrix, digits))

    colnames(data_frame) = column_names

    data_frame
  })

  output$ContrastTests = renderTable({
    results = RunMCPModAnalysis()

    # Extract descriptive statistics
    descriptive_statistics = results$descriptive_statistics

    # Extract the optimal contrasts and contrast correlation matrix
    contrast_results = results$contrast_results

    # Selected models
    selected_models = results$selected_models

    # Extract the test statistics
    mcp_results = results$mcp_results

    # Number of selected models
    n_selected_models = sum(selected_models)

    DF_selected_model_list = DF_model_list[selected_models]

    column_names = c("Model", "Test statistic", "Adjusted p-value", "Significant contrast")

    sign = rep("No", n_selected_models)
    for (i in 1:n_selected_models) {
      if (mcp_results$sign_model[i] == 1) sign[i] = "Yes"
    }

    data_frame = cbind(DF_selected_model_list, 
                      round(mcp_results$test_statistics, digits), 
                      round(mcp_results$adj_pvalues, 4),
                      sign)

    colnames(data_frame) = column_names

    data_frame
  })

  output$DoseResponseModels = renderTable({
    results = RunMCPModAnalysis()

    # Extract descriptive statistics
    descriptive_statistics = results$descriptive_statistics

    # Extract the optimal contrasts and contrast correlation matrix
    contrast_results = results$contrast_results

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    # Extract the test statistics
    mcp_results = results$mcp_results

    # Extract the Mod step results
    mod_results = results$mod_results

    # Extract the list of model fit parameters
    model_fit = mod_results$model_fit

    DF_selected_model_list = DF_model_list[selected_models]

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

    colnames(data_frame) = column_names

    data_frame
  })

  output$ModelSelectionCriteria = renderTable({
    results = RunMCPModAnalysis()

    # Extract input parameters
    input_parameters = results$input_parameters

    column_names = c("Parameter", "Value")

    if (input_parameters$model_selection == "AIC") model_selection_label = "Select the model with the smallest AIC"
    if (input_parameters$model_selection == "maxT") model_selection_label = "Select the model corresponding to the largest test statistic"
    if (input_parameters$model_selection == "aveAIC") model_selection_label = "Average the models using AIC-based weights"

    col1 = c("Model selection criterion", "Delta")
    col2 = c(model_selection_label, input_parameters$delta)

    data_frame = cbind(col1, col2)

    colnames(data_frame) = column_names

    data_frame
  })

  output$ModelSelection = renderTable({
    results = RunMCPModAnalysis()

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    # Extract the test statistics
    mcp_results = results$mcp_results

    # Extract the Mod step results
    mod_results = results$mod_results

    # Extract the list of model fit parameters
    model_fit = mod_results$model_fit

    DF_selected_model_list = DF_model_list[selected_models]

    column_names = c("Model", "Significant contrast", "AIC", "Test statistic", "Model weight")

    sign = rep("No", n_selected_models)
    criterion = rep(NA, n_selected_models)
    test_statistics = rep(NA, n_selected_models)
    for (i in 1:n_selected_models) {
      if (mcp_results$sign_model[i] == 1) {
        sign[i] = "Yes"
        criterion[i] = model_fit[[i]]$criterion
        test_statistics[i] = mcp_results$test_statistics[i]
      }
    }

    model_weight = rep(NA, n_selected_models)
    for (i in 1:n_selected_models) {
      current_criterion = criterion[i]
      denominator = 0
      for (j in 1:n_selected_models) {
        if (mcp_results$sign_model[j] == 1) denominator = denominator + exp(- 0.5 * (criterion[j] - current_criterion))
      }
      if (mcp_results$sign_model[i] == 1 & abs(denominator) > 0.0001) model_weight[i] = 1 / denominator
    }

    data_frame = cbind(DF_selected_model_list, 
                      sign, 
                      sprintf("%0.2f", criterion),
                      sprintf("%0.3f", test_statistics),
                      sprintf("%0.3f", model_weight))

    colnames(data_frame) = column_names

    data_frame
  })

  output$EstimatedTargetDoses = renderTable({
    results = RunMCPModAnalysis()

    # Selected models
    selected_models = results$selected_models

    # Number of selected models
    n_selected_models = sum(selected_models)

    # Extract the test statistics
    mcp_results = results$mcp_results

    # Extract the Mod step results
    mod_results = results$mod_results

    # Extract the list of model fit parameters
    model_fit = mod_results$model_fit

    DF_selected_model_list = DF_model_list[selected_models]

    target_dose = rep(NA, n_selected_models)
    for (i in 1:n_selected_models) {
      if (mcp_results$sign_model[i] == 1) {
        if (model_fit[[i]]$target_dose >= 0) target_dose[i] = model_fit[[i]]$target_dose
      }
    }

    column_names = c("Model", "Target dose")

    data_frame = cbind(DF_selected_model_list, 
                      sprintf("%0.3f", target_dose))

    colnames(data_frame) = column_names

    data_frame
  })

  output$TargetDose = renderTable({
    results = RunMCPModAnalysis()

    # Extract input parameters
    input_parameters = results$input_parameters

    # Selected models
    selected_models = results$selected_models

    # Extract the test statistics
    mcp_results = results$mcp_results

    # Extract the Mod step results
    mod_results = results$mod_results

    # Extract the list of model fit parameters
    model_fit = mod_results$model_fit

    # Number of selected models
    n_selected_models = sum(selected_models)

    DF_selected_model_list = DF_model_list[selected_models]

    criterion = rep(NA, n_selected_models)
    test_statistics = rep(NA, n_selected_models)
    target_dose = rep(NA, n_selected_models)
    for (i in 1:n_selected_models) {
      if (mcp_results$sign_model[i] == 1) {
        criterion[i] = model_fit[[i]]$criterion
        test_statistics[i] = mcp_results$test_statistics[i]
        if (model_fit[[i]]$target_dose >= 0) target_dose[i] = model_fit[[i]]$target_dose
      }
    }

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

    model_weight = rep(NA, n_selected_models)
    for (i in 1:n_selected_models) {
      current_criterion = criterion[i]
      denominator = 0
      for (j in 1:n_selected_models) {
        if (mcp_results$sign_model[j] == 1) denominator = denominator + exp(- 0.5 * (criterion[j] - current_criterion))
      }
      if (mcp_results$sign_model[i] == 1 & abs(denominator) > 0.0001) model_weight[i] = 1 / denominator
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

    colnames(data_frame) = column_names

    data_frame
  })

output$DoseResponseModel = renderPlot({

  results = RunMCPModAnalysis()

  eval_function_list = EvaluateDoseResponseModels()

  # Extract descriptive statistics
  descriptive_statistics = results$descriptive_statistics

  # Extract input parameters
  input_parameters = results$input_parameters

  # Selected models
  selected_models = results$selected_models

  # Number of selected models
  n_selected_models = sum(selected_models)

  # Extract the Mod step results
  mod_results = results$mod_results

  # Extract the list of model fit parameters
  model_fit = mod_results$model_fit

  # Total number of models
  n_models = length(DF_model_list)

  DF_selected_model_list = DF_model_list[selected_models]
  DF_selected_model_index_list = (1:n_models)[selected_models]

  # Extract the dose levels and responses
  dose = mod_results$dose
  resp = mod_results$resp

  # Number of doses
  n_doses = length(descriptive_statistics$dose_levels)

  endpoint_index = input_parameters$endpoint_index

  xlab = "Dose"
  if (input_parameters$endpoint_index == 1) ylab = "Mean response"
  if (input_parameters$endpoint_index == 2) ylab = "Response rate"
  if (input_parameters$endpoint_index == 3) ylab = "Average number of events"

  x_limit = c(min(descriptive_statistics$dose_levels), max(descriptive_statistics$dose_levels))
  y_limit = c(min(descriptive_statistics$lower_cl), max(descriptive_statistics$upper_cl))

  bar_width = (x_limit[2] - x_limit[1]) / 100

  # Exclude the cases with missing coefficients
  if (!any(is.na(eval_function_list$y))) {

    if (min(eval_function_list$y) < y_limit[1]) y_limit[1] = min(eval_function_list$y)
    if (max(eval_function_list$y) > y_limit[2]) y_limit[2] = max(eval_function_list$y)

  }

  plot(x = descriptive_statistics$dose_levels, y = descriptive_statistics$mean_group, xlab=xlab, ylab=ylab, xlim = x_limit, ylim = y_limit, col="black", type="p", pch = 19)
  lines(x = eval_function_list$x, y = eval_function_list$y, col = "black", lwd = 2)
  for (j in 1:n_doses) {
    lines(x = rep(descriptive_statistics$dose_levels[j], 2), y = c(descriptive_statistics$lower_cl[j], descriptive_statistics$upper_cl[j]), col = "black", lwd = 1)
    lines(x = c(descriptive_statistics$dose_levels[j] - bar_width, descriptive_statistics$dose_levels[j] + bar_width), y = rep(descriptive_statistics$lower_cl[j], 2), col = "black", lwd = 1)
    lines(x = c(descriptive_statistics$dose_levels[j] - bar_width, descriptive_statistics$dose_levels[j] + bar_width), y = rep(descriptive_statistics$upper_cl[j], 2), col = "black", lwd = 1)
  }


  })

  output$DownloadResults = downloadHandler(
    filename = function() {
      "Report.docx"
    },

    content = function(file) {
      results = RunMCPModAnalysis()

      doc = GenerateAnalysisReport(results, "MCPMod analysis report")

      # Save the report
      xfile = paste0(file, ".docx")
      print(doc, target = xfile)          
      file.rename(xfile, file)
  })         

  observeEvent(input$jump_to_panel2, {
        updateTabItems(session, "sidebarMCPModAnalysis",
                          selected = "analysis")
  })

  observeEvent(input$jump_to_panel3, {
        updateTabItems(session, "sidebarMCPModAnalysis",
                          selected = "report")
  })

})