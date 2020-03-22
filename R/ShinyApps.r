AnalysisApp = function() {

  appDir = system.file("AnalysisApp", package = "MCPModPack")
  shiny::runApp(appDir, display.mode = "normal")

}

SimulationApp = function() {

  appDir = system.file("SimulationApp", package = "MCPModPack")
  shiny::runApp(appDir, display.mode = "normal")

}
