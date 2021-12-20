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
#' @import shiny shinymeta shinyjs
#'


# wrapper for shiny::shinyApp()
launchPeptidoformVisualisation <- function(maxSize=500) {
  #shinyjs::useShinyjs()
  options(shiny.maxRequestSize=maxSize*1024^2)
  shinyApp(ui = ui, server = server)
}
