# MCPModPack 
## Simulation-Based Design and Analysis of Dose-Finding Trials

[![Downloads from the RStudio CRAN mirror](https://cranlogs.r-pkg.org/badges/MCPModPack)](https://cran.r-project.org/package=MCPModPack)

This package provides an efficient implementation of the MCPMod (Multiple Comparisons and Modeling) method to support a simulation-based design and analysis of dose-finding trials with normally distributed, binary and count endpoints.

## Functions

The following functions are included in the package:

* MCPModAnalysis: Analyze data from a dose-finding trial using MCPMod. 
* MCPModSimulation: Perform a simulation-based evaluation of dose-finding trial designs using MCPMod. 
* AnalysisReport: Generate a detailed summary of MCPMod analysis results in a Microsoft Word format.
* SimulationReport: Generate a detailed summary of MCPMod simulation results in a Microsoft Word format.
* AnalysisApp: Launch an easy-to-use Shiny-based graphical user interface to analyze data from a dose-finding trial.
* SimulationApp: Launch an easy-to-use Shiny-based graphical user interface to perform a simulation-based evaluation of dose-finding trial designs.

## References

Bornkamp, B., Bezlyak, V., Bretz, F. (2015). Implementing the MCP-Mod procedure for dose-response testing and estimation. Modern Approaches to Clinical Trials Using SAS.  Menon, S., Zink, R. (editors). SAS Press: Cary, NC.

Bretz, F., Pinheiro, J.C.,  Branson, M. (2005). Combining multiple comparisons and modeling techniques in dose response studies. Biometrics. 61, 738-748.

Bretz, F., Tamhane, A.C., Pinheiro, J. (2009). Multiple testing in dose response problems. Multiple Testing Problems in Pharmaceutical Statistics. Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). New York: Chapman and Hall/CRC Press.

Nandakumar, S., Dmitrienko, A., Lipkovich, I. (2017). Dose-finding methods. Analysis of Clinical Trials Using SAS: A Practical Guide (Second Edition). Dmitrienko, A., Koch, G.G. (editors). SAS Press: Cary, NC.

Pinheiro, J. C., Bornkamp, B., Bretz, F. (2006). Design and analysis of dose finding studies combining multiple comparisons and modeling procedures. Journal of Biopharmaceutical Statistics. 16, 639-656.

Pinheiro J., Bornkamp B., Glimm E., Bretz F. (2013). Model-based dose finding under model uncertainty using general parametric models. Statistics in Medicine. 33, 1646-1661.