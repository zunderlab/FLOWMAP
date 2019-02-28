
#' GUI for launching shiny Application
#'
#' This code was adapted from the cytofkit R package,
#' available here: \url{https://github.com/JinmiaoChenLab/cytofkit/}.
#'
#' Specifically, the code in this file is adapted from cytofkit_GUI.R, but
#' modified to take user-specified inputs for FLOW-MAP runs and
#' then interface with the FLOWMAPR library back-end.
#'
#' This code, authored and maintained by Jinmiao Chen and Hao Chen,
#' was used under the Artistic License 2.0 or Artistic-2.0 available
#' here: \url{https://opensource.org/licenses/Artistic-2.0}.
#'
#' @examples
#' LaunchGUI()
#' @export
LaunchGUI <- function() {
  require(shiny)
  runApp(appDir = file.path(system.file(package = "FLOWMAPR"), "application"))
}
