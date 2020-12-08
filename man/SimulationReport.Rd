\name{SimulationReport}
\alias{SimulationReport}
\alias{GenerateSimulationReport}
\title{
Generation of a Word-based summary of MCPMod simulation results  
}
\description{
This function creates a detailed summary of MCPMod simulation results in a Microsoft Word format.
}
\usage{
SimulationReport(results, report_title, report_filename)
}
\arguments{
  \item{results}{Object of class \samp{MCPModSimulationResults} created by the \code{\link{MCPModSimulation}} function.  
  }
  \item{report_title}{Character value defining the report's title.
  }
  \item{report_filename}{Character value defining the report's filename. The report is saved in the current working directory.
  }
}
\author{
Alex Dmitrienko <admitrienko@medianainc.com>
}
\seealso{
  \code{\link{MCPModSimulation}}, \code{\link{AnalysisReport}}
}
\examples{
  \donttest{  
# Simulation-based evaluation of dose-finding trials with a binary endpoint

# Endpoint type
endpoint_type = "Binary"

# Select the candidate dose-response models and initial values 
# of the non-linear model parameters (linear, quadratic, exponential, 
# emax, logistic and sigemax)
models = list(linear = NA, 
              quadratic = -0.5, 
              exponential = 0.3, 
              emax = 0.3, 
              logistic = c(0.5, 0.1), 
              sigemax = c(0.5, 5))

# One-sided Type I error rate
alpha = 0.025

# Direction of the dose-response relationship
direction = "increasing"

# Model selection criterion
model_selection = "AIC"

# The treatment effect for identifying the target dose 
# (this effect is defined relative to the placebo effect)
Delta = 0.3

# Select the assumed dose-response model and values of the non-linear model parameters
sim_models = list(linear = NA, 
                  placebo_effect = 0.2, 
                  max_effect = seq(from = 0, to = 0.5, by = 0.1))

# Simulation parameters 
sim_parameters = list(n = rep(40, 5),
                      doses = c(0, 0.05, 0.2, 0.6, 1),
                      dropout_rate = 0.05,
                      nsims = 1000,
                      go_threshold = 0.3)

# Perform an MCPMod-based simulation
results = MCPModSimulation(endpoint_type = endpoint_type, 
                           models = models, 
                           alpha = alpha, 
                           direction = direction, 
                           model_selection = model_selection, 
                           Delta = Delta,
                           sim_models = sim_models,
                           sim_parameters = sim_parameters)

# Simple summary of the MCPMod simulation results
results

# Detailed summary of the MCPMod simulation results (remove tempfile)
SimulationReport(results, 
  "MCPMod simulation summary (Binary endpoint)", 
  tempfile("MCPMod simulation summary (Binary endpoint).docx", fileext=".docx")) 
  }
}