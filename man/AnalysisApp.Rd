\name{AnalysisApp}
\alias{AnalysisApp}
\title{
Creation of a Shiny-based interface to perform an MCPMod analysis of a dose-finding trial
}
\description{
This function creates a Shiny-based graphical user interface to perform MCPMod-based analysis of a dose-finding trial. The application requires data sets with the dose and response information (comma-separated value file). The data set is required to include the \code{dose} and \code{resp} variables with a single record per patient, see for example the \code{\link{normal}} data set.
}
\usage{
AnalysisApp()
}
\arguments{
  \item{data}{The data set is required to include the dose and resp variables with a single record per patient. }  
}
\author{
Alex Dmitrienko <admitrienko@medianainc.com>
}
\seealso{
  \code{\link{MCPModAnalysis}}, \code{\link{AnalysisReport}}
}
\examples{
  \donttest{  
# Launch a Shiny-based interface to perform an MCPMod-based analysis 
# of a dose-finding trial
AnalysisApp()
	}
}
