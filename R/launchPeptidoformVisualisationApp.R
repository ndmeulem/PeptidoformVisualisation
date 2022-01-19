#' Launches the Peptidoform Visualisation Shiny App
#'
#' @param maxSize maximum memory size that input files are allowed to have in Mb
#'
#' @export launchPeptidoformVisualisation
#'
#' @return shiny application object
#'
#' @example
#' \dontrun{launchPeptidoformVisualisation()}
#'
#' @import shiny shinymeta shinyjs BiocManager
#'


# wrapper for shiny::shinyApp()
launchPeptidoformVisualisation <- function(maxSize=500) {
  shinyjs::useShinyjs()
  onStart = shinyhelper::observe_helpers(help_dir = system.file("helpfiles", package="PeptidoformVisualisation"))
  options(shiny.maxRequestSize=maxSize*1024^2)
  shinyApp(ui = ui, server = server)
}
