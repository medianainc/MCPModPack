\name{AnalysisReport}
\alias{AnalysisReport}
\alias{GenerateAnalysisReport}
\title{
Generation of a Word-based summary of MCPMod analysis results  
}
\description{
This function creates a detailed summary of MCPMod analysis results in a Microsoft Word format.
}
\usage{
AnalysisReport(results, report_title, report_filename)
}
\arguments{
  \item{results}{Object of class \samp{MCPModAnalysisResults} created by the \code{\link{MCPModAnalysis}} function.  
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
  \code{\link{MCPModAnalysis}}, \code{\link{SimulationReport}}
}
\examples{
  \donttest{  
# MCPMod-based analysis of a dose-finding trial with a 
# binary endpoint

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