#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
reactomeVisualizer <- function(pathId, placeHolder, width = 930, height = 500, elementId = NULL) {

  # forward options using x
  x <- list(pathId, placeHolder)
  # create widget
  htmlwidgets::createWidget(
                   name = 'reactomeVisualizer',
                   x = x,
                   width = width,
                   height = height,
                   package = 'reactomeVisualizer',
                   elementId = elementId
               )
}

#' Shiny bindings for reactomeVisualizer
#'
#' Output and render functions for using reactomeVisualizer within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a reactomeVisualizer
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name reactomeVisualizer-shiny
#'
#' @export
reactomeVisualizerOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'reactomeVisualizer', width, height, package = 'reactomeVisualizer')
}

#' @rdname reactomeVisualizer-shiny
#' @export
renderReactomeVisualizer <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, reactomeVisualizerOutput, env, quoted = TRUE)
}
