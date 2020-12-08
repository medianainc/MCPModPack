\name{MCPModAnalysis}
\alias{MCPModAnalysis}
\title{
MCPMod-based analysis of dose-finding clinical trials with normally distributed, binary and count endpoints
}
\description{
This function implements the MCPMod-based analysis of dose-finding clinical trials with normally distributed, binary and count endpoints, including derivation of the optimal contrasts for the candidate dose-response models, evaluation of dose-response tests based on the optimal contrasts, selection of the significant dose-response models and estimation of the target dose. For more information, see the technical manual in the package's \code{doc} folder.
}
\usage{
MCPModAnalysis(endpoint_type, models, dose, resp, alpha, direction, 
               model_selection, Delta, theta)
}
\arguments{
  \item{endpoint_type}{Character value defining the primary endpoint's type. Possible values:  
   \itemize{
    \item \code{"Normal"}: Normally distributed primary endpoint.
    \item \code{"Binary"}: Binary primary endpoint.
    \item \code{"Count"}: Count-type primary endpoint.
  }
  }
  \item{models}{List of candidate dose-response models with initial values of the non-linear model parameters. The package supports the following dose-response models: linear, quadratic, exponential, Emax, logistic and sigEmax. No initial value is required for the linear model, a single initial value is required for the quadratic, exponential and Emax models, and two initial values are required for the logistic and sigEmax models.
  }
  \item{dose, resp}{Numeric vectors of equal length specifying the dose and response values.
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
}
\value{
The function returns an object of class \samp{MCPModAnalysisResults}. This object is a list with the following components:
  \item{input_parameters}{a list containing the user-specified parameters, e.g, endpoint type, model selection criteria, etc.}
  \item{selected_models}{a logical vector defining the candidate dose-response models.}
  \item{descriptive_statistics}{a list containing the descriptive statistics computed from the trial's data set.}
  \item{contrast_results}{a list containing the contrast evaluation results for the candidate dose-response models, e.g., the model-specific optimal dose-response contrasts and contrast correlation matrix.}
  \item{mcp_results}{a list containing the multiplicity adjustment results for the candidate dose-response models, e.g., the model-specific test statistics and adjusted p-values.}
  \item{mod_results}{a list containing the modeling results for the candidate dose-response models, e.g., estimated model parameters, target dose estimate.}
  A detailed summary of the MCPMod analysis results can be generated using the AnalysisReport function.  
}

\author{
Alex Dmitrienko <admitrienko@medianainc.com>
}

\seealso{
  \code{\link{MCPModSimulation}}
}
\examples{
  \donttest{
# MCPMod-based analysis of a dose-finding trial with a binary endpoint

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

# Perform an MCPMod-based analysis of the trial's data
# The data set binary is included in the package
results = MCPModAnalysis(endpoint_type = endpoint_type, 
                     models = models, 
                     dose = binary$dose, 
                     resp = binary$resp, 
                     alpha = alpha, 
                     direction = direction, 
                     model_selection = model_selection, 
                     Delta = Delta)

# Simple summary of the MCPMod analysis results
results

# Detailed summary of the MCPMod analysis results (remove tempfile)
AnalysisReport(results, 
  "MCPMod analysis summary (Binary endpoint)", 
  tempfile("MCPMod analysis summary (Binary endpoint).docx", fileext=".docx")) 
  }
}


