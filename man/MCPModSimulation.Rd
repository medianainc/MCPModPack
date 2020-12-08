\name{MCPModSimulation}
\alias{MCPModSimulation}
\alias{ComputeDRFunctionParameters}
\alias{EvaluateDRFunction}
\title{
MCPMod-based simulation of dose-finding trial designs
}
\description{
This function implements the simulation-based analysis of dose-finding clinical trials with normally distributed, binary and count endpoints using the MCPMod methodology. For more information, see the technical manual in the package's \code{doc} folder.
}
\usage{
MCPModSimulation(endpoint_type, models, alpha, direction, 
                 model_selection, Delta, theta, sim_models, sim_parameters)
}
\arguments{
  \item{endpoint_type}{Character value defining the primary endpoint's type. Possible values:  
   \itemize{
    \item \code{"Normal"}: Normally distributed primary endpoint.
    \item \code{"Binary"}: Binary primary endpoint.
    \item \code{"Count"}: Count-type primary endpoint.
  }
  }
  \item{models}{List of candidate dose-response models with initial values of the non-linear model parameters. The package supports the following dose-response models: linear, quadratic, exponential, emax, logistic and sigemax. No initial value is required for the linear model, a single initial value is required for the quadratic, exponential and Emax models, and two initial values are required for the logistic and sigEmax models.
  }
  \item{alpha}{Numeric value defining the one-sided significance level (default value is 0.025).
  }
  \item{direction}{Character value defining the direction of the dose-response relationship. Possible values:  
   \itemize{
    \item \code{"Increasing"}: A larger value of the treatment difference corresponds to a beneficial treatment effect.
    \item \code{"Decreasing"}: A smaller value of the treatment difference indicates a beneficial treatment effect.
  }
  }
  \item{model_selection}{Character value defining the criterion for selecting the best dose-response model. Possible values:  
   \itemize{
    \item \code{"AIC"}: Akaike information criterion (AIC).
    \item \code{"maxT"}: Most significant test statistic.
    \item \code{"aveAIC"}: Weighted AIC-based criterion.
  }
  }
  \item{Delta}{Numeric value specifying the treatment effect for identifying the target dose. The treatment effect is defined relative to the placebo effect. A positive value is required if \code{direction = "Increasing"} and a negative value is required otherwise.
  }
  \item{theta}{Numeric vector defining the overdispersion parameter in each trial arm (required only with count-type primary endpoints).
  }    
  \item{sim_models}{List defining the assumed dose-response model and initial values of the non-linear parameters in the simulated trial. The package supports the following dose-response models: linear, quadratic, exponential, Emax, logistic and sigEmax. No initial value is required for the linear model, a single initial value is required for the quadratic, exponential and Emax models, and two initial values are required for the logistic and sigEmax models. The following components are required:
   \itemize{
    \item \code{"placebo_effect"}: Numeric value defining the placebo effect in the assumed dose-response model.
    \item \code{"max_effect"}: Numeric vector defining the effects at the maximum dose in the assumed dose-response model. These effects are defined relative to the placebo effect. Positive values are required if \cr \code{direction = "Increasing"} and negative values are required otherwise.
    \item \code{"sd"}: Numeric vector defining the standard deviations of the response variable in each trial arm (required for normally distributed endpoints).
  }
  }
  \item{sim_parameters}{List defining the design and simulation parameters in the simulated trial. The following components are required:
   \itemize{
    \item \code{"n"}: Integer vector defining the number of patients in each trial arm.
    \item \code{"doses"}: Numeric vector defining the dose levels in each trial arm.
    \item \code{"dropout_rate"}: Numeric value defining the dropout rate in the simulated trial (between 0 and 1).
    \item \code{"nsims"}: Integer value defining the number of simulation runs.
    \item \code{"go_threshold"}: Numeric value specifying the threshold for computing go probabilities, i.e., probabilities that the maximum effect for the best dose-response model corresponding to a significant contrast exceeds a pre-defined value. The threshold is defined relative to the placebo effect. A positive value is required if \code{direction =} \code{"Increasing"} and a negative value is required otherwise.
  }
  }

}
\value{
The function returns an object of class \samp{MCPModSimulationResults}. This object is a list with the following components:
  \item{input_parameters}{a list containing the user-specified parameters, e.g, endpoint type, model selection criteria, etc.}
  \item{selected_models}{a logical vector defining the candidate dose-response models.}
  \item{sim_results}{a list containing the simulation results based on the assumed dose-response model, e.g., power, target dose estimates, etc.}
  A detailed summary of the simulation results can be generated using the SimulationReport function.
}

\author{
Alex Dmitrienko <admitrienko@medianainc.com>
}
\seealso{
  \code{\link{MCPModAnalysis}}
}
\examples{
  \donttest{
# Simulation-based evaluation of dose-finding trials with a count endpoint

# Endpoint type
endpoint_type = "Count"

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
Delta = 2

# Vector of overdispersion parameters
theta = c(2, 2, 2, 2, 2)

# Select the assumed dose-response model and values of the non-linear model parameters
sim_models = list(emax = 1, 
                  placebo_effect = 3, 
                  max_effect = seq(from = 0, to = 3, by = 1))

# Simulation parameters 
# (go threshold is defined relative to the placebo effect)
sim_parameters = list(n = c(50, 50, 50, 50, 50),
                      doses = c(0, 0.05, 0.2, 0.6, 1),
                      dropout_rate = 0.05,
                      nsims = 1000,
                      go_threshold = 2)

# Perform an MCPMod-based simulation
results = MCPModSimulation(endpoint_type = endpoint_type, 
                           models = models, 
                           alpha = alpha, 
                           direction = direction, 
                           model_selection = model_selection, 
                           Delta = Delta,
                           theta = theta,
                           sim_models = sim_models,
                           sim_parameters = sim_parameters)

# Simple summary of the MCPMod simulation results
results

# Detailed summary of the MCPMod simulation results (remove tempfile)
SimulationReport(results, 
  "MCPMod simulation summary (Count endpoint)", 
  tempfile("MCPMod simulation summary (Count endpoint).docx", fileext=".docx"))

  }
}