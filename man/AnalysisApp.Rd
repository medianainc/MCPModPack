\name{AnalysisApp}
\alias{AnalysisApp}
\title{
Creation of a Shiny-based interface to perform an MCPMod analysis of a dose-finding trial
}
\description{
This function creates a Shiny-based graphical user interface to perform MCPMod-based analysis of a dose-finding trial. The application requires data sets with the dose and response information (comma-separated value file). The data set is required to include the \code{dose} and \code{resp} variables with a single record per patient. The first row of the file should contain the required variable names as follows: "dose","resp" (including the quotation marks and the comma). Subsequent rows should contain comma-separated values corresponding to the dose and response values for each patient (quotation marks are not required for the rows of data).  See for example the \code{\link{normal}} data set. 
}
\usage{
AnalysisApp()
}
\author{
Alex Dmitrienko <admitrienko@medianainc.com>
}
\seealso{
  \code{\link{MCPModAnalysis}}, \code{\link{AnalysisReport}}
}
