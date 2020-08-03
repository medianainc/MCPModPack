library(shiny)
library(shinydashboard)

shinyUI(
  dashboardPage(
    dashboardHeader(title = "MCPMod Analysis"),
    dashboardSidebar(
      sidebarMenu(id = "sidebarMCPModAnalysis",
        menuItem("Parameters", tabName = "parameters", icon = icon("sliders-h")),
        menuItem("Analysis", tabName = "analysis", icon = icon("chart-bar")),
        menuItem("Report", tabName = "report", icon = icon("file")),
        menuItem("References", tabName = "references", icon = icon("link"))
      )
    ),
    dashboardBody(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "main.css")
      ),
      tabItems(
        # Parameters tab
        tabItem(tabName = "parameters",
          fluidRow(
            box(
              title = "Analysis of dose-finding trials using MCPMod",
              solidHeader = FALSE,
              collapsible = TRUE,
              width = 12,

              "This web-based tool performs an MCPMod-based analysis of a dose-finding trial."
            )
          ),

          fluidRow(
            column(6,
              box(
                title = "Upload data set",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tags$p("Select a data set with the dose and response information (comma-separated value file):"),
                fileInput('data_set', '', accept = c('.csv') ),
                tags$p(class = "help-block",
                  "The data set is required to include the dose and resp variables with a single record per patient.")

              ),

              box(
                title = "Trial parameters",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                selectInput("endpoint_index", label = "Primary endpoint's type",
                          c("Normal" = 1, "Binary" = 2, "Count" = 3)),

                selectInput("direction", label = "Direction of the dose-response relationship",
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
                numericInput(inputId = "delta", label = "Treatment effect for identifying the target dose",value = 0.2),
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

            column(6,
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
            )
          ),

          fluidRow(
            box(
              title = "Next step",
              background = "red",
              width = 12,

              tags$p("Proceed to the next tab to perform the MCPMod analysis of the trial's data."),
              actionButton("jump_to_panel2", "Next tab")
            )
          )
        ),

        # Analysis tab
        tabItem(tabName = "analysis",
          fluidRow(
            box(
              title = "Summary of MCPMod analysis results",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = 12,

              textOutput("error_message")
            )
          ),

          fluidRow(
            column(12, class = "col-lg-6",
              box(
                title = "Descriptive statistics",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("DescriptiveStatistics")
              )
            ),
            column(12, class = "col-lg-6",
              box(
                title = "Model-specific dose-response contrasts",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("Contrasts")
              )
            )
          ),

          fluidRow(
            column(12, class = "col-lg-6",
              box(
                title = "Contrast correlation matrix",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("CorrelationMatrix")
              )
            ),
            column(12, class = "col-lg-6",
              box(
                title = "Model-specific contrast tests",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("ContrastTests")
              )
            )
          ),

          fluidRow(
            column(12, class = "col-lg-5",
              box(
                title = "Parameters of dose-response models",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("DoseResponseModels"),
                tags$p(class = "help-block", "The convergence criterion is defined as the length of the gradient vector evaluated at the maximum likelihood estimate and a high value of the convergence criterion suggests lack of convergence.")
              )
            ),

            column(12, class = "col-lg-7",
              box(
                title = "Model selection",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                h5("Model selection parameters"),
                tableOutput("ModelSelectionCriteria"),

                h5("Model selection criteria"),
                tableOutput("ModelSelection")
              )
            )
          ),

          fluidRow(
            column(12, class = "col-lg-5",
              box(
                title = "Target dose estimation",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                h5("Model-specific estimated target doses"),
                tableOutput("EstimatedTargetDoses"),

                h5("Selected model and target dose"),
                tableOutput("TargetDose")
              )
            ),

            column(12, class = "col-lg-7",
              box(
                title = "Dose-response models",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                radioButtons("select_plot_model", "",
                            c("Linear" = 1, "Quadratic" = 2, "Exponential" = 3, "Emax" = 4, "Logistic" = 5, "SigEmax" = 6),
                            1, TRUE),

                plotOutput("DoseResponseModel")
              )
            )
          ),

          fluidRow(
            box(
              title = "Next step",
              background = "red",
              width = 12,

              tags$p("Proceed to the next tab to create a detailed analysis report."),
              actionButton("jump_to_panel3", "Next tab")
            )
          )
        ),

        # Report tab
        tabItem(tabName = "report",
          fluidRow(
            box(
              title = "Create an MCPMod analysis report",
              background = "red",
              width = 12,

              tags$p("Click the Download button to create and save a detailed MCPMod analysis report in a Microsoft Word format."),
              downloadButton("DownloadResults", "Download")
            )
          )
        ),

        # References tab
        tabItem(tabName = "references",
          fluidRow(
            column(12, class = "col-lg-6",
              wellPanel(
                h4("MCPModAnalysis function"),

                tags$p("This Shiny application is based on the MCPModAnalysis function that implements the MCPMod-based analysis of dose-finding clinical trials with normally distributed, binary and count endpoints. For more information on the underlying methodology, see the technical manual in the package's doc folder."),

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
