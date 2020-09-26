library(shiny)
library(shinydashboard)

shinyUI(
  dashboardPage(
    dashboardHeader(title = "MCPMod Simulation"),
    dashboardSidebar(
      sidebarMenu(id = "sidebarMCPModSimulation",
        menuItem("Parameters", tabName = "parameters", icon = icon("sliders-h")),
        menuItem("Dose-response models", tabName = "responseModels", icon = icon("chart-area")),
        menuItem("Simulation", tabName = "simulation", icon = icon("table")), # chart-bar
        menuItem("Report", tabName = "report", icon = icon("file")),
        menuItem("References", tabName = "references", icon = icon("link"))
      )
    ),
    dashboardBody(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "main.css?2")
      ),
      tabItems(
        # Parameters tab
        tabItem(tabName = "parameters",

          fluidRow(
            box(
              title = "Simulation of dose-finding trials using MCPMod",
              solidHeader = FALSE,
              collapsible = TRUE,
              width = 12,

              "This web-based tool performs a simulation-based evaluation of dose-finding trial designs using MCPMod."
            )
          ),

          fluidRow(
            column(6, class = "col-md-4",
              box(
                title = "Trial parameters",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                selectInput("endpoint_index", label = "Primary endpoint's type",
                          c("Normal" = 1, "Binary" = 2, "Count" = 3)),

                selectInput("direction", label = "Dose-response relationship",
                            c("Increasing" = 1, "Decreasing" = -1)),
                tags$p(class = "help-block",
                  "With an increasing dose-response relationship, the endpoint's higher value corresponds to a beneficial treatment effect and, with a decreasing dose-response relationship, the endpoint's lower value indicates a beneficial treatment effect.")
              ),

              box(
                title = "Model and target dose selection",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                selectInput("model_selection", label = "Model selection criterion",
                          c("AIC" = 1, "maxT" = 2, "aveAIC" = 3)),

                conditionalPanel(
                  condition = "input.delta <= 0 && input.direction == 1",
                  div(class = "alert alert-danger", "Treatment effect for identifying the target dose must be positive if the direction of the dose-response relationship is Increasing.")
                ),
                conditionalPanel(
                  condition = "input.delta >= 0 && input.direction == -1",
                  div(class = "alert alert-danger", "Treatment effect for identifying the target dose must be negative if the direction of the dose-response relationship is Decreasing.")
                ),
                numericInput(inputId = "delta", label = "Treatment effect for identifying the target dose", value = 0.4),
                tags$p(class = "help-block",
                "The treatment effect for identifying the target dose is defined relative to the placebo effect.")
              ),

              box(
                title = "Other parameters",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                conditionalPanel(
                  condition = "input.endpoint_index == 3",
                  conditionalPanel(
                    condition = "!RegExp('^\\\\s*\\\\d*\\\\.?\\\\d+\\\\s*(,\\\\s*\\\\d*\\\\.?\\\\d+\\\\s*)*$').test(input.theta)",
                    div(class = "alert alert-danger", "Overdispersion parameters must be specified using a list of positive numbers separated by a comma.")
                  ),
                  textInput(inputId = "theta", label = "List of overdispersion parameters in trial arms",
                              value = "2, 2, 2, 2, 2")
                ),

                conditionalPanel(
                  condition = "input.alpha <= 0 || input.alpha >= 1",
                  div(class = "alert alert-danger", "One-sided Type I error rate value must be > 0 and < 1.")
                ),
                numericInput(inputId = "alpha", label = "One-sided Type I error rate (alpha)",
                            value = 0.025, min = 0.001, max = 0.999)
              )
            ),

            column(6, class = "col-md-4",
              box(
                title = "Candidate dose-response models",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tags$p("Select the candidate dose-response models and specify the initial values of model parameters."),
                conditionalPanel(
                  condition = paste("input.linear_model == false && ",
                                    "input.quadratic_model == false && ",
                                    "input.exponential_model == false && ",
                                    "input.emax_model == false && ",
                                    "input.logistic_model == false && ",
                                    "input.sigemax_model == false"),
                  div(class = "alert alert-danger", "At least one model must be specified.")
                ),
                hr(),

                checkboxInput("linear_model", "Linear model", TRUE),
                hr(),

                checkboxInput("quadratic_model", "Quadratic model", TRUE),
                conditionalPanel(
                  condition = "input.quadratic_model == true",
                  numericInput(inputId = "quadratic_model_par1", label = "Delta", value = -0.5)
                ),
                hr(),

                checkboxInput("exponential_model", "Exponential model", TRUE),
                conditionalPanel(
                  condition = "input.exponential_model == true",
                  conditionalPanel(
                    condition = "input.exponential_model_par1 <= 0",
                    div(class = "alert alert-danger", "Delta must be > 0.")
                  ),
                  numericInput(inputId = "exponential_model_par1", label = "Delta", value = 0.3)
                ),
                hr(),

                checkboxInput("emax_model", "Emax model", TRUE),
                conditionalPanel(
                  condition = "input.emax_model == true",
                  conditionalPanel(
                    condition = "input.emax_model_par1 <= 0",
                    div(class = "alert alert-danger", "ED50 must be > 0.")
                  ),
                  numericInput(inputId = "emax_model_par1", label = "ED50", value = 0.3)
                ),
                hr(),
        
                checkboxInput("logistic_model", "Logistic model", TRUE),
                conditionalPanel(
                  condition = "input.logistic_model == true",
                  conditionalPanel(
                    condition = "input.logistic_model_par1 <= 0 || input.logistic_model_par2 <= 0",
                    div(class = "alert alert-danger", "Both ED50 and Delta must be > 0.")
                  ),
                  numericInput(inputId = "logistic_model_par1", label = "ED50", value = 0.5),
                  numericInput(inputId = "logistic_model_par2", label = "Delta", value = 0.1)
                ),
                hr(),
        
                checkboxInput("sigemax_model", "SigEmax model", TRUE),
                conditionalPanel(
                  condition = "input.sigemax_model == true",
                  conditionalPanel(
                    condition = "input.sigemax_model_par1 <= 0 || input.sigemax_model_par2 <= 0",
                    div(class = "alert alert-danger", "Both ED50 and h must be > 0.")
                  ),
                  numericInput(inputId = "sigemax_model_par1", label = "ED50", value = 0.5),
                  numericInput(inputId = "sigemax_model_par2", label = "h", value = 5)
                )
              )
            ),

            column(6, class = "col-md-4",
              box(
                title = "Simulation model",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                selectInput("sim_model_index", label = "Assumed dose-response model",
                            c("Linear" = 1, "Quadratic" = 2, "Exponential" = 3, "Emax" = 4,
                              "Logistic" = 5, "SigEmax" = 6)),

                conditionalPanel(
                  condition = "input.sim_model_index == 3",
                  conditionalPanel(
                    condition = "input.sim_model_exponential_delta <= 0",
                    div(class = "alert alert-danger", "Delta must be > 0.")
                  ),
                  numericInput(inputId = "sim_model_exponential_delta", label = "Delta", value = 1)
                ),
                conditionalPanel(
                  condition = "input.sim_model_index == 4",
                  conditionalPanel(
                    condition = "input.sim_model_emax_ed50 <= 0",
                    div(class = "alert alert-danger", "ED50 must be > 0.")
                  ),
                  numericInput(inputId = "sim_model_emax_ed50", label = "ED50", value = 1)
                ),
                conditionalPanel(
                  condition = "input.sim_model_index == 5",
                  conditionalPanel(
                    condition = "input.sim_model_logistic_ed50 <= 0 || input.sim_model_logistic_delta <= 0",
                    div(class = "alert alert-danger", "Both ED50 and Delta must be > 0.")
                  ),
                  numericInput(inputId = "sim_model_logistic_ed50", label = "ED50", value = 1),
                  numericInput(inputId = "sim_model_logistic_delta", label = "Delta", value = 1)
                ),
                conditionalPanel(
                  condition = "input.sim_model_index == 6",
                  conditionalPanel(
                    condition = "input.sim_model_sigemax_ed50 <= 0 || input.sim_model_sigemax_h <= 0",
                    div(class = "alert alert-danger", "Both ED50 and h must be > 0.")
                  ),
                  numericInput(inputId = "sim_model_sigemax_ed50", label = "ED50", value = 1),
                  numericInput(inputId = "sim_model_sigemax_h", label = "h", value = 1)
                ),

                hr(),

                numericInput(inputId = "sim_models_placebo_effect", label = "Placebo effect",
                            value = 0.4, min = 0.01, max = 0.99, step = 0.01),

                conditionalPanel(
                  condition = "!RegExp('^\\\\s*\\\\-?\\\\d*\\\\.?\\\\d+\\\\s*(,\\\\s*\\\\-?\\\\d*\\\\.?\\\\d+\\\\s*)*$').test(input.sim_models_max_effect)",
                  div(class = "alert alert-danger", "Maximum effects must be specified using a list of numbers separated by a comma.")
                ),
                textInput(inputId = "sim_models_max_effect", label = "List of maximum effects over placebo",
                            value = "0, 0.1, 0.2, 0.3, 0.4, 0.5"),
                tags$p(class = "help-block",
                  "The maximum effect is defined as the difference between the effect at the highest dose and placebo effect."),

                conditionalPanel(
                  condition = "input.endpoint_index == 1",
                  conditionalPanel(
                    condition = "!RegExp('^\\\\s*\\\\d*\\\\.?\\\\d+\\\\s*(,\\\\s*\\\\d*\\\\.?\\\\d+\\\\s*)*$').test(input.sim_models_sd)",
                    div(class = "alert alert-danger", "Standard deviations must be specified using a list of positive numbers separated by a comma.")
                  ),
                  textInput(inputId = "sim_models_sd", label = "List of standard deviations in trial arms",
                              value = "0.5, 0.5, 0.75, 0.75, 0.95")
                )
              )
            ),

            column(6, class = "col-md-4",
              box(
                title = "Simulation parameters",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                conditionalPanel(
                  condition = "!RegExp('^[0-9]+([ ]*,[ ]*[0-9]+)*$').test(input.sim_parameters_n)",
                  div(class = "alert alert-danger", "Sample sizes must be specified using a list of positive integers separated by a comma.")
                ),
                textInput(inputId = "sim_parameters_n", label = "List of sample sizes in trial arms",
                            value = "40, 40, 40, 40, 40"),

                conditionalPanel(
                  condition = "!RegExp('^\\\\s*\\\\d*\\\\.?\\\\d+\\\\s*(,\\\\s*\\\\d*\\\\.?\\\\d+\\\\s*)*$').test(input.sim_parameters_doses)",
                  div(class = "alert alert-danger", "Dose levels must be specified using a list of non-negative numbers separated by a comma.")
                ),
                textInput(inputId = "sim_parameters_doses", label = "List of dose levels",
                            value = "0, 0.05, 0.2, 0.6, 1"),

                conditionalPanel(
                  condition = "input.sim_parameters_dropout_rate < 0 || input.sim_parameters_dropout_rate >= 1",
                  div(class = "alert alert-danger", "Dropout rate must be >= 0 and < 1.")
                ),
                numericInput(inputId = "sim_parameters_dropout_rate", label = "Dropout rate", 
                              value = 0.05, min = 0, max = 1, step = 0.05),

                conditionalPanel(
                  condition = "input.sim_parameters_go_threshold < 0 && input.direction == 1",
                  div(class = "alert alert-danger", "Threshold for computing go probabilities must be positive if the direction of the dose-response relationship is Increasing.")
                ),
                conditionalPanel(
                  condition = "input.sim_parameters_go_threshold > 0 && input.direction == -1",
                  div(class = "alert alert-danger", "Threshold for computing go probabilities must be negative if the direction of the dose-response relationship is Decreasing.")
                ),
                numericInput(inputId = "sim_parameters_go_threshold", label = "Threshold for computing go probabilities", 
                              value = 0.3, min = -100, max = 100, step = 1),
                tags$p(class = "help-block",
                "The threshold for computing go probabilities is defined relative to the placebo effect."),

                conditionalPanel(
                  condition = "!RegExp('^[0-9]+$').test(input.sim_parameters_nsims) || input.sim_parameters_nsims < 1",
                  div(class = "alert alert-danger", "Number of simulation runs must be a positive integer.")
                ),
                numericInput(inputId = "sim_parameters_nsims", label = "Number of simulation runs", 
                              value = 1000, min = 100, max = 10000, step = 100)
              )
            )
          ),

          fluidRow(
            box(
              title = "Next step",
              background = "red",
              width = 12,

              tags$p("Proceed to the next tab to view assumed dose-response models."),
              actionButton("jump_to_panel2", "Next tab")
            )
          )

        ),

        tabItem(tabName = "responseModels",
          fluidRow(
            box(
              title = "Assumed dose-response models",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              width = 12,

              radioButtons("select_scenario_model", "Maximum effect over placebo:",
                          c("0" = 1, "0.1" = 2, "0.2" = 3, "0.3" = 4, "0.4" = 5, "0.5" = 6),
                          1, TRUE),

              plotOutput("AssumedDoseResponse")
            ),

            box(
              title = "Next step",
              background = "red",
              width = 12,

              tags$p("Proceed to the next tab to view simulation summary."),
              actionButton("jump_to_panel3", "Next tab")
            )

          )
        ),

        # Simulation tab
        tabItem(tabName = "simulation",
          fluidRow(
            box(
              title = "Summary of MCPMod simulation results",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = 12
            ),
            
            column(12, class="col-md-4",
              box(
                title = "Power summary",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("PowerSummary"),

                tags$p(class = "help-block",
                  "Power is the probability that the best dose-response contrast is significant.")

              )
            ),
          
            conditionalPanel(
              condition = "input.model_selection == 1",
              column(12, class="col-md-8",
                box(
                  title = "Probability of selecting a dose-response model",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = NULL,

                  tableOutput("ProbabilityOfModel")

                )
              )
            )
          ),

          fluidRow(
            column(12, class = "col-md-7",
              box(
                title = "Target dose estimates",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("TargetDoseEstimates")
              )
            ),

            column(12, class = "col-md-5",
              box(
                title = "Summary of go probabilities",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                conditionalPanel(
                  condition = "input.model_selection != 3",
                  tableOutput("GoProbabilities"),

                tags$p(class = "help-block",
                  "The go probability is the probability that the best dose-response contrast is significant and the maximum effect for the corresponding model exceeds the pre-defined go threshold.")
                
                ),

                conditionalPanel(
                  condition = "input.model_selection == 3",
                  tags$p("Go probabilities cannot be computed if model averaging (aveAIC) is requested.")
                )
              )
            )
          ),

          fluidRow(
            column(12,
              box(
                title = "Probability of selecting a dose",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("ProbabilityOfSelectingADose"),

                tags$p(class = "help-block",
                  "Each column presents the probability that the estimated target dose is less than or equal to the current dose and is strictly greater than the next lower dose.")


              )
            )
          ),

          fluidRow(
            column(12, class = "col-lg-6",
              box(
                title = "Power",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                plotOutput("Power")
              )
            ),

            column(12, class = "col-lg-6",
              box(
                title = "Go probability",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                conditionalPanel(
                  condition = "input.model_selection != 3",
                  plotOutput("GoProbability")
                ),

                conditionalPanel(
                  condition = "input.model_selection == 3",
                  tags$p("Go probabilities cannot be computed if model averaging (aveAIC) is requested.")
                )
              )
            )
          ),

          fluidRow(
            column(12, class = "col-lg-8",
              box(
                title = "Assumed and estimated dose-response models",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                conditionalPanel(
                  condition = "input.model_selection != 3",
                  
                  radioButtons("select_scenario", "Maximum effect over placebo:",
                              c("0" = 1, "0.1" = 2, "0.2" = 3, "0.3" = 4, "0.4" = 5, "0.5" = 6),
                              1, TRUE),

                  plotOutput("EstimatedDoseResponse"),
                  tags$p("Red curve: Assumed dose-response curve. Black curve: Estimated dose-response curve with a 95% confidence band.")
                ),

                conditionalPanel(
                  condition = "input.model_selection == 3",
                  tags$p("Dose-response models cannot be estimated if model averaging (aveAIC) is requested.")
                )
              )
            )
          ),

          fluidRow(
            box(
              title = "Next step",
              background = "red",
              width = 12,

              tags$p("Proceed to the next tab to create a detailed simulation report."),
              actionButton("jump_to_panel4", "Next tab")
            )
          )                     

        ),

        tabItem(tabName = "report",
          fluidRow(
            box(
              title = "Create an MCPMod simulation report",
              background = "red",
              width = 12,

              tags$p("Click the Download button to create and save a detailed MCPMod simulation report in a Microsoft Word format."),
              downloadButton("DownloadResults", "Download")
            )
          )
        ),

        tabItem(tabName = "references",
          fluidRow(
            column(12, class = "col-lg-6",
              wellPanel(
                h4("MCPModSimulation function"),

                tags$p("This Shiny application is based on the MCPModSimulation function that implements the simulation-based analysis of dose-finding clinical trials with normally distributed, binary and count endpoints using the MCPMod methodology. For more information on the underlying methodology, see the technical manual in the package's doc folder."),

                br(),
                h4("References"),

                tags$p("Bornkamp, B., Bezlyak, V., Bretz, F. (2015). Implementing the MCP-Mod procedure for dose-response testing and estimation. Modern Approaches to Clinical Trials Using SAS.  Menon, S., Zink, R. (editors). SAS Press: Cary, NC."),
                tags$p("Bretz, F., Pinheiro, J.C.,  Branson, M. (2005). Combining multiple comparisons and modeling techniques in dose response studies. Biometrics. 61, 738-748."),
                tags$p("Bretz, F., Tamhane, A.C., Pinheiro, J. (2009). Multiple testing in dose response problems. Multiple Testing Problems in Pharmaceutical Statistics. Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). New York: Chapman and Hall/CRC Press."),
                tags$p("Genz, A., Bretz, F. (2002). Methods for the computation of multivariate t-probabilities. Journal of Computational and Graphical Statistics. 11, 950-971."),
                tags$p("Nandakumar, S., Dmitrienko, A., Lipkovich, I. (2017). Dose-finding methods. Analysis of Clinical Trials Using SAS: A Practical Guide (Second Edition). Dmitrienko, A., Koch, G.G. (editors). SAS Press: Cary, NC."),
                tags$p("Pinheiro, J.C., Bornkamp, B., Bretz, F. (2006). Design and analysis of dose finding studies combining multiple comparisons and modeling procedures. Journal of Biopharmaceutical Statistics. 16, 639-656."),
                tags$p("Pinheiro J., Bornkamp B., Glimm E., Bretz F. (2013). Model-based dose finding under model uncertainty using general parametric models. Statistics in Medicine. 33, 1646-1661.")

              )
            )
          )
        )
      )
    )
  )
)
