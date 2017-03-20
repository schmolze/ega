#' @export
#' @title Plot a Clarke Error Grid
#' @description The function uses \code{\link[ggplot2]{ggplot}} to draw the
#' Clarke error grid lines according to the criteria described in the
#' original publication by Clarke et. al. (see reference below). If zones
#' have not already been assigned via the \code{zones} parameter, the
#' function \code{\link{getClarkeZones}} is called first. The values in
#' \code{referenceVals} and \code{testVals} are then superimposed as a scatter
#' plot. Some basic plot parameters can be specified as arguments, but the
#' return value can also be stored and modified further before plotting
#' (see examples and vignette).
#' @param referenceVals A vector of glucose values obtained via the reference
#' method.
#' @param testVals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{referenceVals}, so the length should be the same.
#' @param title The main plot title. Defaults to "Clarke Error Grid".
#' @param xlab The x-axis label. Defaults to "Reference Glucose
#' Concentration (mg/dL)".
#' @param ylab The y-axis label. Defaults to "Test Glucose Concentration
#' (mg/dL)".
#' @param linesize The size to be used when drawing the zone lines. The
#' acceptable values are the same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is 0.5.
#' @param linetype The type of line to be used when drawing the zone lines. The
#' acceptable values are the same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is "solid".
#' @param linecolor The color of the zone lines. The acceptable values are the
#' same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is "black".
#' @param linealpha The alpha (transparency) level to be used when drawing
#' the zone lines. The acceptable values are the same as for
#' \code{\link[ggplot2]{geom_segment}}. The default is 0.6.
#' @param pointsize The size to be used when plotting the glucose data points.
#' The acceptable values are the same as for \code{\link[ggplot2]{geom_point}}.
#' The default is 2.
#' @param pointalpha The alpha (transparency) level to be used when plotting
#' the glucose data points. The acceptable values are the same as for
#' \code{\link[ggplot2]{geom_point}}. The default is 1.
#' @param zones An optional character vector specifying the Clarke zones
#' for each paired value. If this is not supplied, \code{\link{getClarkeZones}}
#' will be called to generate zone labels.
#' @param unit A string specifying the units of measurement. This should be either
#' \code{"gram"} (the default) for \code{mg/dl} or \code{"mol"} for \code{mmol/l}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the return
#' value is not assigned, a plot is drawn.
#' @examples
#' library(ggplot2)
#'
#' # default
#' plotClarkeGrid(glucose_data$ref, glucose_data$test)
#'
#' # with options
#' plotClarkeGrid(glucose_data$ref, glucose_data$test,
#'               pointsize=1.5,
#'               pointalpha=0.6,
#'               linetype="dashed")
#'
#' # store return value and modify
#' ceg <- plotClarkeGrid(glucose_data$ref, glucose_data$test)
#'
#' ceg + theme_gray() +
#'    theme(plot.title = element_text(size = rel(2), colour = "blue"))
#' @seealso
#' \code{\link{getClarkeZones}} \code{\link[ggplot2]{ggplot}}
#' @references
#' Clarke, W. L., D. Cox, L. A. Gonder-Frederick, W. Carter, and S. L. Pohl.
#' "Evaluating Clinical Accuracy of Systems for Self-Monitoring of Blood
#' Glucose." Diabetes Care 10, no. 5 (September 1, 1987): 622-28.
plotClarkeGrid <- function (referenceVals, testVals, title = "Clarke Error Grid", xlab="", ylab="",
                            linesize = 0.5, linetype = "solid", linecolor = "black",
                            linealpha = 0.6, pointsize = 2, pointalpha = 1, zones = NA, unit='gram')
{
  #unit selection and control
  if (unit != "mol" & unit != "gram") {
    stop("'unit' must be either 'mol' or 'gram'.")
  }
  #axis named depending on unit
  if (unit == "mol") {
    n <- 18 #where n is a scaling factor for unit conversion
    if (xlab==""){
      xlab="Reference Glucose Concentration (mmol/L)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mmol/L)"
    }
  } else {
    n <- 1
    if (xlab==""){
      xlab="Reference Glucose Concentration (mg/dL)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mg/dL)"
    }
  }

  # use default zone assignment if none is provided
  if (is.na (zones)) {
    zones <- getClarkeZones (referenceVals, testVals, unit)
  }
  tolerance <- 0.2

  # create a df for ggplot (NULL to appease CRAN)
  ref <- test <- NULL
  data <- data.frame (ref=referenceVals, test=testVals, zones=zones)

  #better solution for scaling axis automatically with some extra space
  maxX <- max (max (data$ref) + 20 / n, 550 / n)
  maxY <- max ( (data$test + 20 / n), (1 + tolerance) * maxX, 650 / n)

  #labels with coordinats and colors
  labels <- data.frame (x=c (240 / n, 120 / n, 350 / n, 120 / n, 163 / n,
                             35 / n, 350 / n, 35 / n, 350 / n),
                        y=c (230 / n, 200 / n, 230 / n, 300 / n, 20 / n,
                             130 / n, 130 / n, 300 / n, 35 / n),
                        label=c ("A", "B", "B", "C", "C", "D", "D", "E", "E"),
                        color=c ("blue", "blue", "blue", "blue", "blue",
                                 "red", "red", "red", "red"))

  #segment endpoints for borders (NULL to appease CRAN)
  x1 <- y1 <- xend <- yend <- NULL
  border <- data.frame (x1=c (58.3 / n, 70 / n, 70 / n, 70 / n, 0, 240 / n,
                              0, 70 / n, 180 / n, 240 / n, 180 / n, 130 / n),
                        y1=c (70 / n, 56 / n, 180 / n, 83 / n, 180 / n,
                              180 / n, 70 / n, 0, 70 / n, 70 / n, 0, 0),
                        xend=c (maxX, maxX, 550 / n, 70 / n, 70 / n, maxX,
                                58.3 / n, 70 / n, maxX, 240 / n,
                                180 / n, 180 / n),
                        yend=c ((1 + tolerance) * maxX, (1 - tolerance) * maxX,
                                660 / n, maxY, 180 / n, 180 / n, 70 / n,
                                56 / n, 70 / n, 180 / n, 70 / n, 70 / n))

  ceg <- ggplot(data, aes(x = ref, y = test)) +
    scale_x_continuous (breaks = c (round (70 / n, digits=1),
                                    round (100 / n, digits=1),
                                    round (150 / n, digits=1),
                                    round (180 / n, digits=1),
                                    round (240 / n, digits=1),
                                    round (300 / n, digits=1),
                                    round (350 / n, digits=1),
                                    round (400 / n, digits=1),
                                    round (450 / n, digits=1),
                                    round (500 / n, digits=1),
                                    round (550 / n, digits=1),
                                    round (600 / n, digits=1),
                                    round (650 / n, digits=1),
                                    round (700 / n, digits=1),
                                    round (750 / n, digits=1),
                                    round (800 / n, digits=1),
                                    round (850 / n, digits=1),
                                    round (900 / n, digits=1),
                                    round (950 / n, digits=1),
                                    round (1000 / n, digits=1)),
                        expand = c (0, 0)) +
    scale_y_continuous (breaks = c (round (70 / n, digits=1),
                                    round (100 / n, digits=1),
                                    round (150 / n, digits=1),
                                    round (180 / n, digits=1),
                                    round (240 / n, digits=1),
                                    round (300 / n, digits=1),
                                    round (350 / n, digits=1),
                                    round (400 / n, digits=1),
                                    round (450 / n, digits=1),
                                    round (500 / n, digits=1),
                                    round (550 / n, digits=1),
                                    round (600 / n, digits=1),
                                    round (650 / n, digits=1),
                                    round (700 / n, digits=1),
                                    round (750 / n, digits=1),
                                    round (800 / n, digits=1),
                                    round (850 / n, digits=1),
                                    round (900 / n, digits=1),
                                    round (950 / n, digits=1),
                                    round (1000 / n, digits=1)),
                        expand = c (0, 0)) +
    geom_point(aes(color = zones), size = pointsize, alpha = pointalpha) +
    geom_segment (aes (x = x1, y = y1, xend = xend, yend = yend),
                  data=border, linetype=linetype) +
    annotate (geom="text", x = labels$x, y = labels$y, size = 6,
              label = labels$label, color=labels$color) +
    theme_bw () +
    #    theme (legend.position = "none") +
    ggtitle (title) +
    xlab (xlab) +
    ylab (ylab)
  ceg
}


#' @export
#' @title Plot a Parkes (Consensus) Error Grid
#' @description The function uses \code{\link[ggplot2]{ggplot}} to draw the
#' Parkes (consensus) error grid lines according to the criteria described in
#' the publications listed in the References section (see below). If zones
#' have not already been assigned via the \code{zones} parameter, the
#' function \code{\link{getParkesZones}} is called first. The values in
#' \code{referenceVals} and \code{testVals} are then superimposed as a scatter
#' plot. Some basic plot parameters can be specified as arguments, but the
#' return value can also be stored and modified further before plotting
#' (see examples and vignette).
#' @param referenceVals A vector of glucose values obtained via the reference
#' method.
#' @param testVals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{referenceVals}, so the length should be the same.
#' @param type An integer (1 or 2) specifying whether to plot the grid for Type 1
#' or Type 2 diabetes. Defaults to 1.
#' @param title The main plot title. Defaults to "Parkes (Consensus) Error Grid
#' for Type [type] Diabetes".
#' @param xlab The x-axis label. Defaults to "Reference Glucose
#' Concentration (mg/dL)".
#' @param ylab The y-axis label. Defaults to "Test Glucose Concentration
#' (mg/dL)".
#' @param linesize The size to be used when drawing the zone lines. The
#' acceptable values are the same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is 0.5.
#' @param linetype The type of line to be used when drawing the zone lines. The
#' acceptable values are the same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is "solid".
#' @param linecolor The color of the zone lines. The acceptable values are the
#' same as for \code{\link[ggplot2]{geom_segment}}.
#' The default is "black".
#' @param linealpha The alpha (transparency) level to be used when drawing
#' the zone lines. The acceptable values are the same as for
#' \code{\link[ggplot2]{geom_segment}}. The default is 0.6.
#' @param pointsize The size to be used when plotting the glucose data points.
#' The acceptable values are the same as for \code{\link[ggplot2]{geom_point}}.
#' The default is 2.
#' @param pointalpha The alpha (transparency) level to be used when plotting
#' the glucose data points. The acceptable values are the same as for
#' \code{\link[ggplot2]{geom_point}}. The default is 1.
#' @param zones An optional character vector specifying the Clarke zones
#' for each paired value. If this is not supplied, \code{\link{getClarkeZones}}
#' will be called to generate zone labels.
#' @param unit A string specifying the units of measurement. This should be either
#' \code{"gram"} (the default) for \code{mg/dl} or \code{"mol"} for \code{mmol/l}.
#' @return A \code{\link[ggplot2]{ggplot}} object is returned. If the return
#' value is not assigned, a plot is drawn.
#' @examples
#' library(ggplot2)
#'
#' # default
#' plotParkesGrid(glucose_data$ref, glucose_data$test)
#'
#' # with options
#' plotParkesGrid(glucose_data$ref, glucose_data$test,
#'               pointsize=2,
#'               pointalpha=0.5,
#'               linesize=2,
#'               linealpha=0.3,
#'               linetype="dotdash")
#'
#' # store return value and modify
#' peg <- plotParkesGrid(glucose_data$ref, glucose_data$test, type=2)
#'
#' peg + theme_gray() +
#'    theme(plot.title = element_text(size = rel(2), colour = "red"))
#' @seealso
#' \code{\link{getParkesZones}} \code{\link[ggplot2]{ggplot}}
#' @references
#' Parkes, J. L., S. L. Slatin, S. Pardo, and B.H. Ginsberg. "A New Consensus
#' Error Grid to Evaluate the Clinical Significance of Inaccuracies in the
#' Measurement of Blood Glucose." Diabetes Care 23, no. 8 (August 2000):
#' 1143-48
#'
#' Pfutzner, Andreas, David C. Klonoff, Scott Pardo, and Joan L. Parkes.
#' "Technical Aspects of the Parkes Error Grid." Journal of Diabetes Science
#' and Technology 7, no. 5 (September 2013): 1275-81
plotParkesGrid <- function (referenceVals, testVals, type=1, title="", xlab="", ylab="", linesize = 0.5,
                            linetype = "solid", linecolor = "black", linealpha = 0.6, pointsize = 2, pointalpha = 1, zones = NA, unit='gram') {

  if (type != 1 & type != 2) {
    stop("'type' must be 1 or 2.")
  }
  if (unit != "mol" & unit != "gram") {
    stop("'unit' must be either 'mol' or 'gram'.")
  }
  if (title==""){
    title <- paste("Parkes (Consensus) Error Grid for Type ", type, " Diabetes")
  }
  if (unit == "mol") {
    n <- 18
    if (xlab==""){
      xlab="Reference Glucose Concentration (mmol/L)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mmol/L)"
    }
  } else {
    n <- 1
    if (xlab==""){
      xlab="Reference Glucose Concentration (mg/dL)"
    }
    if (ylab==""){
      ylab="Test Glucose Concentration (mg/dL)"
    }
  }


  if (is.na (zones)) {
    zones <- getParkesZones (referenceVals, testVals, type, unit)
  }

  # NULL to appease CRAN
  ref <- test <- NULL
  data <- data.frame (ref=referenceVals, test=testVals, zones=zones)

  maxX <- max (max (data$ref) + 20 / n, 550 / n)#better solution
  maxY <- max (max (data$test) + 20 / n, 550 / n)

  labels <- data.frame (x=c (320 / n, 220 / n, 385 / n, 140 / n, 405 / n, 415 / n, 75 / n, 21 / n),
                        y=c (320 / n, 360 / n, 235 / n, 375 / n, 145 / n,  50 / n, 383 / n, 383 / n),
                        label=c ("A", "B", "B", "C", "C", "D", "D", "E"),
                        color=c ("green", "blue", "blue", "red", "red", "red", "red", "red"))


  if (type==1){#type 1 diabetes
    ce <- .coef (35, 155, 50, 550)
    cdu <- .coef (80, 215, 125, 550)
    cdl <- .coef (250, 40, 550, 150)
    ccu <- .coef (70, 110, 260, 550)
    ccl <- .coef (260, 130, 550, 250)
    cbu <- .coef (280, 380, 430, 550)
    cbl <- .coef (385, 300, 550, 450)

    # NULL to appease CRAN
    x1 <- y1 <- xend <- yend <- NULL

    border <- data.frame (x1=c (0 / n,
                                0 / n, 30 / n, 140 / n, 280 / n,
                                50 / n, 50 / n, 170 / n, 385 / n,
                                0 / n, 30 / n, 50 / n, 70 / n,
                                120 / n, 120 / n, 260 / n,
                                0 / n, 25 / n, 50 / n, 80 / n,
                                250 / n, 250 / n,
                                0 / n, 35 / n),
                          y1=c (0 / n,
                                50 / n, 50 / n, 170 / n, 380 / n,
                                0 / n, 30 / n, 145 / n, 300 / n,
                                60 / n, 60 / n, 80 / n, 110 / n,
                                0 / n, 30 / n, 130 / n,
                                100 / n, 100 / n, 125 / n, 215 / n,
                                0 / n, 40 / n,
                                150 / n, 155 / n),
                          xend=c (min (maxX, maxY),
                                  30 / n, 140 / n, 280 / n, .endx (280 / n, 380 / n, maxY, cbu),
                                  50 / n, 170 / n, 385 / n, maxX,
                                  30 / n, 50 / n, 70 / n, .endx (70 / n, 110 / n, maxY, ccu),
                                  120 / n, 260 / n, maxX,
                                  25 / n, 50 / n, 80 / n, .endx (80 / n, 215 / n, maxY, cdu),
                                  250 / n, maxX,
                                  35 / n, .endx (35 / n, 155 / n, maxY, ce) ),
                          yend=c (min (maxX, maxY),
                                  50 / n, 170 / n, 380 / n, maxY,
                                  30 / n, 145 / n, 300 /n, .endy (385 / n, 300 / n, maxX, cbl),
                                  60 / n, 80 / n, 110 / n, maxY,
                                  30 / n, 130 / n, .endy (260 / n, 130 / n, maxX, ccl),
                                  100 / n, 125 / n, 215 / n, maxY,
                                  40 / n, .endy (410 / n, 110 / n, maxX, cdl),
                                  155 / n, maxY))

  } else {#type 2 diabetes
    ce <- .coef (35, 200, 50, 550)
    cdu <- .coef (35, 90, 125, 550)
    cdl <- .coef (410, 110, 550, 160)
    ccu <- .coef (30, 60, 280, 550)
    ccl <- .coef (260, 130, 550, 250)
    cbu <- .coef (230, 330, 440, 550)
    cbl <- .coef (330, 230, 550, 450)

    border <- data.frame (x1=c (0 / n,
                                0 / n, 30 / n, 230 / n,
                                50 / n, 50 / n, 90 / n, 330 / n,
                                0 / n, 30 / n,
                                90 / n, 260 / n,
                                0 / n, 25 / n, 35 / n,
                                250 / n, 250 / n, 410 / n,
                                0 / n, 35 / n),
                          y1=c (0 / n,
                                50 / n, 50 / n, 330 / n,
                                0 / n, 30 / n, 80 / n, 230 / n,
                                60 / n, 60 / n,
                                0 / n, 130 / n,
                                80 / n, 80 / n, 90 / n,
                                0 / n, 40 / n, 110 / n,
                                200 / n, 200 / n),
                          xend=c (min (maxX, maxY),
                                  30 / n, 230 / n, .endx (230 / n, 330 / n, maxY, cbu),
                                  50 / n, 90 / n, 330 / n, maxX,
                                  30 / n, .endx (30 / n, 60 / n, maxY, ccu),
                                  260 / n, maxX,
                                  25 / n, 35 / n, .endx (35 / n, 90 / n, maxY, cdu),
                                  250 / n, 410 / n, maxX,
                                  35 / n, .endx (35 / n, 200 / n, maxY, ce) ),
                          yend=c (min (maxX, maxY),
                                  50 / n, 330 / n, maxY,
                                  30 / n, 80 / n, 230 / n, .endy (330 / n, 230 / n, maxX, cbl),
                                  60 / n, maxY,
                                  130 / n, .endy (260 / n, 130 / n, maxX, ccl),
                                  80 / n, 90 / n, maxY,
                                  40 / n, 110 / n, .endy (410 / n, 110 / n, maxX, cdl),
                                  200 / n, maxY))

  }

  peg <- ggplot(data, aes(x = ref, y = test)) +
    scale_x_continuous (breaks = c (round (70 / n, digits=1),
                                    round (100 / n, digits=1),
                                    round (150 / n, digits=1),
                                    round (180 / n, digits=1),
                                    round (240 / n, digits=1),
                                    round (300 / n, digits=1),
                                    round (350 / n, digits=1),
                                    round (400 / n, digits=1),
                                    round (450 / n, digits=1),
                                    round (500 / n, digits=1),
                                    round (550 / n, digits=1),
                                    round (600 / n, digits=1),
                                    round (650 / n, digits=1),
                                    round (700 / n, digits=1),
                                    round (750 / n, digits=1),
                                    round (800 / n, digits=1),
                                    round (850 / n, digits=1),
                                    round (900 / n, digits=1),
                                    round (950 / n, digits=1),
                                    round (1000 / n, digits=1)),
                        expand = c (0, 0)) +
    scale_y_continuous (breaks = c (round (70 / n, digits=1),
                                    round (100 / n, digits=1),
                                    round (150 / n, digits=1),
                                    round (180 / n, digits=1),
                                    round (240 / n, digits=1),
                                    round (300 / n, digits=1),
                                    round (350 / n, digits=1),
                                    round (400 / n, digits=1),
                                    round (450 / n, digits=1),
                                    round (500 / n, digits=1),
                                    round (550 / n, digits=1),
                                    round (600 / n, digits=1),
                                    round (650 / n, digits=1),
                                    round (700 / n, digits=1),
                                    round (750 / n, digits=1),
                                    round (800 / n, digits=1),
                                    round (850 / n, digits=1),
                                    round (900 / n, digits=1),
                                    round (950 / n, digits=1),
                                    round (1000 / n, digits=1)),
                        expand = c (0, 0)) +
    geom_point(aes(color = zones), size=pointsize, alpha=pointalpha) +
    geom_segment (aes (x=x1, y=y1, xend=xend, yend=yend), data=border, linetype=linetype) +
    annotate (geom="text", x=labels$x, y=labels$y, size=6, label=labels$label, color=labels$color) +
    theme_bw () +
    #    theme (legend.position = "none") +
    ggtitle (title) +
    xlab (xlab) +
    ylab (ylab)
  peg
}


