\name{MCPModPack-package}
\alias{MCPModPack-package}
\alias{MCPModPack}
\alias{TestDoseResponseFunction}
\alias{TestFindTargetDose}
\docType{package}
\title{
Design and analysis of dose-finding trials
}
\description{
The MCPModPack package facilitates the design and analysis of dose-finding clinical trials with normally distributed, binary and count endpoints using the MCPMod methodology. 
}
\details{
\tabular{ll}{
Package: \tab MCPModPack\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2020-08-01\cr
License: \tab GPL-3
}
Key functions included in the package:
\itemize{
\item \code{\link{MCPModAnalysis}}: Analyze data from a dose-finding trial using MCPMod.
\item \code{\link{AnalysisReport}}: Generate a detailed summary of MCPMod analysis results in a Microsoft Word format.
\item \code{\link{AnalysisApp}}: Launch a Shiny-based graphical user interface to analyze data from a dose-finding trial.
\item \code{\link{MCPModSimulation}}: Perform a simulation-based evaluation of dose-finding trial designs using MCPMod.
\item \code{\link{SimulationReport}}: Generate a detailed summary of MCPMod simulation results in a Microsoft Word format.
\item \code{\link{SimulationApp}}: Launch a Shiny-based graphical user interface to perform a simulation-based evaluation of dose-finding trial designs.
}
The package comes with three example data sets:
\itemize{
\item \code{\link{normal}}: Data set based on a dose-finding trial with a normally distributed endpoint.
\item \code{\link{binary}}: Data set based on a dose-finding trial with a binary endpoint.
\item \code{\link{count}}: Data set based on a dose-finding trial with a count endpoint. 
}
}
\author{
Alex Dmitrienko <admitrienko@medianainc.com>
}
\references{
Bornkamp, B., Bezlyak, V., Bretz, F. (2015). Implementing the MCP-Mod procedure for dose-response testing and estimation. \emph{Modern Approaches to Clinical Trials Using SAS}.  Menon, S., Zink, R. (editors). SAS Press: Cary, NC.

Bretz, F., Pinheiro, J.C.,  Branson, M. (2005). Combining multiple comparisons and modeling techniques in dose response studies. \emph{Biometrics}. 61, 738-748.

Bretz, F., Tamhane, A.C., Pinheiro, J. (2009). Multiple testing in dose response problems. \emph{Multiple Testing Problems in Pharmaceutical Statistics}. Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). New York: Chapman and Hall/CRC Press.

Nandakumar, S., Dmitrienko, A., Lipkovich, I. (2017). Dose-finding methods. \emph{Analysis of Clinical Trials Using SAS: A Practical Guide} (Second Edition). Dmitrienko, A., Koch, G.G. (editors). SAS Press: Cary, NC.

Pinheiro, J. C., Bornkamp, B., Bretz, F. (2006). Design and analysis of dose finding studies combining multiple comparisons and modeling procedures. \emph{Journal of Biopharmaceutical Statistics}. 16, 639-656.

Pinheiro J., Bornkamp B., Glimm E., Bretz F. (2013). Model-based dose finding under model uncertainty using general parametric models. \emph{Statistics in Medicine}. 33, 1646-1661.
}
\keyword{ package }
\examples{
    \donttest{
# MCPMod-based analysis of a dose-finding trial with a normally distributed endpoint

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
Delta = 0.5

# Perform an MCPMod-based analysis of the trial's data
# The data set normal is included in the package
results = MCPModAnalysis(endpoint_type = "Normal", 
                     models = models, 
                     dose = normal$dose, 
                     resp = normal$resp, 
                     alpha = alpha, 
                     direction = direction, 
                     model_selection = model_selection, 
                     Delta = Delta)

# Simple summary of the MCPMod analysis results
results

# Detailed summary of the MCPMod analysis results (remove tempfile)
AnalysisReport(results, 
  "MCPMod analysis summary (Normally distributed endpoint)", 
  tempfile("MCPMod analysis summary (Normally distributed endpoint).docx", fileext=".docx")) 
  }
}
