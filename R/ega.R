#' Clarke and Parkes (Consensus) error grid analysis
#'
#' @name ega-package
#' @import ggplot2
#' @import mgcv
#' @docType package
NULL


#' 5072 paired reference and test glucose values.
#'
#' A dataset containing 5072 paired reference method and test method
#' glucose values (in mg/dL).
#'
#' @format A data frame with 5072 rows and 2 variables:
#' \describe{
#'   \item{ref}{Reference method glucose value, in mg/dL}
#'   \item{test}{Test method glucose value, in mg/dL}
#' }
#' @source The data is from a modified clinical dataset.
"glucose_data"


#' @export
#' @title Assign Clarke error grid zones to paired glucose values
#' @description \code{referenceVals} and \code{testVals} are assumed to contain
#' paired glucose values from a reference method and a test method,
#' respectively. \code{unit} contains info on the unit of measurement. Two
#' options exist: \code{"gram"} for mg/dL and \code{"mol"} for mmol/l
#' with \code{"gram"} applied by default. The discrepancy between the two values
#' is used to place the pair into a Clarke error grid zone according to the
#' criteria described in the original paper by Clarke et. al. (see reference below).
#' @param referenceVals A vector of glucose values obtained via the reference
#' method.
#' @param testVals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{referenceVals}, so the length should be the same.
#' @param unit A string specifying the units of measurement. This should be either
#' \code{"gram"} (the default) for \code{mg/dl} or \code{"mol"} for \code{mmol/l}.
#' @return A character vector is returned, with each element being one of
#' \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, or \code{"E"}.
#' @examples
#' zones <- getClarkeZones (glucose_data$ref / 18, glucose_data$test / 18,
#' unit="mol")
#'
#' # counts
#' table(zones)
#'
#' # percentages
#' round (table (zones) / length (zones) * 100, digits=2)
#'
#' @references
#' Clarke, W. L., D. Cox, L. A. Gonder-Frederick, W. Carter, and S. L. Pohl.
#' "Evaluating Clinical Accuracy of Systems for Self-Monitoring of Blood
#' Glucose." Diabetes Care 10, no. 5 (September 1, 1987): 622-28.
getClarkeZones <- function (referenceVals, testVals, unit="gram"){
  if (unit != "mol" & unit != "gram") {
    stop("'unit' must be either 'mol' or 'gram'.")
  }
  if (unit == "mol") {
    n <- 18 #scaling factor for mmol/l conversion from mg/dl
  } else {
    n <- 1
  }
  zones <- vector (mode="character", length=length (referenceVals))
  bias <- testVals - referenceVals

  # absolute relative error = abs(bias)/reference*100
  are <- abs (bias) / referenceVals * 100
  eq1 <- (7 / 5) * (referenceVals - 130 / n)
  eq2 <- referenceVals + 110 / n

  # zone D: ref < 70 and (test > 70 and test < 180) or
  #   ref > 240 and (test > 70 and test < 180)
  test_D <- testVals >= 70 / n & testVals < 180 / n#error corrected >=70 instead of >70
  zoneD <- (referenceVals < 70 / n & test_D) |
    (referenceVals > 240 / n & test_D)
  zones[zoneD] <- "D"

  # assign A after D, since part of A will overwrite D

  # zone C: (ref >= 130 and ref <= 180 and test < eq1) or
  #   (ref > 70 and ref > 180 and ref > eq2)
  zoneC <- (referenceVals >= 130 / n & referenceVals <= 180 / n & testVals < eq1) |
    (referenceVals > 70 / n & testVals > 180 / n & testVals > eq2)
  zones[zoneC] <- "C"


  #Assign A after C, since part of C will override A

  # zone A: are <= 20  or (ref < 58.3 and test < 70)
  zoneA <- (are <= 20) |
    (referenceVals < 70 / n & testVals < 70 / n)#error solved
  zones[zoneA] <- "A"

  # zone E: (ref <= 70 and test >= 180) or (ref >=180 and test <=70)
  zoneE <- (referenceVals <= 70 / n & testVals >= 180 / n) |
    (referenceVals >= 180 / n & testVals <= 70 / n)
  zones[zoneE] <- "E"

  # the rest are zone B
  zones <- replace (zones, zones == "", "B")

  return (zones)

}


#' @export
#' @title Assign Parkes (Consensus) error grid zones to paired glucose values
#' @description \code{referenceVals} and \code{testVals} are assumed to contain
#' paired glucose values from a reference method and a test method,
#' respectively. The discrepancy between the two values, as well as the
#' type of error grid desired (Type 1 or Type 2 diabetes), is used to place the
#' pair into a Parkes (Consensus) error grid zone, according to the
#' criteria described in the second reference below. \code{unit} contains info
#' on the unit of measurement. Two options exist: \code{"gram"} for mg/dL and
#' \code{"mol"} for mmol/l with \code{"gram"} applied by default.
#' @param referenceVals A vector of glucose values obtained via the reference
#' method.
#' @param testVals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{referenceVals}, so the length should be the same.
#' @param type An integer (1 or 2) specifying whether to obtain zones for Type 1
#' or Type 2 diabetes. Defaults to 1.
#' @param unit A string specifying the units of measurement. This should be either
#' \code{"gram"} (the default) for \code{mg/dl} or \code{"mol"} for \code{mmol/l}.
#' @return A character vector is returned, with each element being one of
#' \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, or \code{"E"}.
#' @examples
#' zones <- getParkesZones(glucose_data$ref, glucose_data$test)
#'
#' # counts
#' table(zones)
#'
#' # percentages
#' round (table (zones) / length (zones) * 100, digits=2)
#' @references
#' Parkes, J. L., S. L. Slatin, S. Pardo, and B.H. Ginsberg. "A New Consensus
#' Error Grid to Evaluate the Clinical Significance of Inaccuracies in the
#' Measurement of Blood Glucose." Diabetes Care 23, no. 8 (August 2000):
#' 1143-48
#'
#' Pfutzner, Andreas, David C. Klonoff, Scott Pardo, and Joan L. Parkes.
#' "Technical Aspects of the Parkes Error Grid." Journal of Diabetes Science
#' and Technology 7, no. 5 (September 2013): 1275-81
getParkesZones <- function (referenceVals, testVals, type=1, unit="gram"){
  if (type != 1 & type != 2){
    stop("'type' must be 1 or 2.")
  }
  if (unit != "mol" & unit != "gram") {
    stop("'unit' must be either 'mol' or 'gram'.")
  }
  if (unit == "mol") {
    n <- 18 #scaling factor for mmol/l conversion from mg/dl
  } else {
    n <- 1
  }

  #setting the graph limits with some space to accomondate all the datapoints
  maxX <- max (max (referenceVals) + 20 / n, 550 / n)#better solution
  maxY <- max ( (testVals + 20 / n), maxX, 550 / n)

  #common block
  testdf <- as.matrix (cbind (as.numeric (referenceVals), as.numeric (testVals)))
  zones <- vector (mode = "character", length = length (referenceVals))
  zones[1:length (referenceVals)] <- "A" #all datapoints are A by default

  #diffirenciation by diabites type
  if (type==1){
    #determining line coefficients for final segments
    ce <- .coef (35, 155, 50, 550)
    cdu <- .coef (80, 215, 125, 550)
    cdl <- .coef (250, 40, 550, 150)
    ccu <- .coef (70, 110, 260, 550)
    ccl <- .coef (260, 130, 550, 250)
    cbu <- .coef (280, 380, 430, 550)
    cbl <- .coef (385, 300, 550, 450)

    #diabetes type 1 zones - creates polygons of a size dependent on data
    limitE1 <- matrix (data= c (0, 35 / n, .endx (35 / n, 155 / n, maxY, ce), 0, 0, #x limits E upper
                                150 / n, 155 / n, maxY, maxY, 150 / n),#y limits E upper
                       ncol=2, byrow=FALSE)

    limitD1L <- matrix (data= c (250 / n, 250 / n, maxX, maxX, 250 / n,#x limits D lower
                                 0, 40 / n, .endy (410 / n, 110 / n, maxX, cdl), 0, 0),#y limits D lower
                        ncol=2, byrow=FALSE)

    limitD1U <- matrix (data= c (0, 25 / n, 50 / n, 80 / n, .endx (80 / n, 215 / n, maxY, cdu), 0, 0,#x limits D upper
                                 100 / n, 100 / n, 125 / n, 215 / n , maxY, maxY, 100 / n),#y limits D upper
                        ncol=2, byrow=FALSE)

    limitC1L <- matrix (data= c (120 / n, 120 / n, 260 / n, maxX, maxX, 120 / n, #x limits C lower
                                 0, 30 / n, 130 / n, .endy (260 / n, 130 / n, maxX, ccl) , 0, 0),#y limits C lower
                        ncol=2, byrow=FALSE)

    limitC1U <- matrix (data= c (0, 30 / n, 50 / n, 70 / n, .endx (70 / n, 110 / n, maxY, ccu), 0, 0, #x limits C upper
                                 60 / n, 60 / n, 80 / n, 110 / n , maxY, maxY, 60 / n),#y limits C upper
                        ncol=2, byrow=FALSE)

    limitB1L <- matrix (data= c (50 / n, 50 / n, 170 / n, 385 / n, maxX, maxX, 50 / n, #x limits B lower
                                 0, 30 / n, 145 / n, 300 / n , .endy (385 / n, 300 / n, maxX, cbl), 0, 0),#y limits B lower
                        ncol=2, byrow=FALSE)


    limitB1U <- matrix (data= c (0, 30 / n, 140 / n, 280 / n, .endx (280 / n, 380 / n, maxY, cbu), 0, 0, #x limits B upper
                                 50 / n, 50 / n, 170 / n, 380 / n , maxY, maxY, 50 / n),#y limits B upper
                        ncol=2, byrow=FALSE)

    #labelling zones using in.out function from mgcv package
    zones[which (in.out (limitB1L, testdf))] <- "B"
    zones[which (in.out (limitB1U, testdf))] <- "B"
    zones[which (in.out (limitC1L, testdf))] <- "C"
    zones[which (in.out (limitC1U, testdf))] <- "C"
    zones[which (in.out (limitD1L, testdf))] <- "D"
    zones[which (in.out (limitD1U, testdf))] <- "D"
    zones[which (in.out (limitE1, testdf))] <- "E"

  }else{ #type 2 diabetes
    ce <- .coef (35, 200, 50, 550)
    cdu <- .coef (35, 90, 125, 550)
    cdl <- .coef (410, 110, 550, 160)
    ccu <- .coef (30, 60, 280, 550)
    ccl <- .coef (260, 130, 550, 250)
    cbu <- .coef (230, 330, 440, 550)
    cbl <- .coef (330, 230, 550, 450)
    #diabetes type 2 zones - creates polygons of a size dependent on data
    limitE2 <- matrix (data= c (0, 35 / n, .endx (35 / n, 200 / n, maxY, ce), 0, 0, #x limits E upper
                                200 / n, 200 / n, maxY, maxY, 200 / n),#y limits E upper
                       ncol=2, byrow=FALSE)

    limitD2L <- matrix (data= c (250 / n, 250 / n, 410 / n, maxX, maxX, 250 / n,#x limits D lower
                                 0, 40 / n, 110 / n, .endy (410 / n, 110 / n, maxX, cdl), 0, 0),#y limits D lower
                        ncol=2, byrow=FALSE)

    limitD2U <- matrix (data= c (0, 25 / n, 35 / n, .endx (35 / n, 90 / n, maxY, cdu), 0, 0, #x limits D upper
                                 80 / n, 80 / n, 90 / n, maxY, maxY, 80 / n),#y limits D upper
                        ncol=2, byrow=FALSE)

    limitC2L <- matrix (data= c (90 / n, 260 / n, maxX, maxX, 90 / n, #x limits C lower
                                 0, 130 / n, .endy (260 / n, 130 / n, maxX, ccl), 0, 0),#y limits C lower
                        ncol=2, byrow=FALSE)

    limitC2U <- matrix (data= c (0, 30 / n, .endx (30 / n, 60 / n, maxY, ccu), 0, 0, #x limits C upper
                                 60 / n, 60 / n, maxY, maxY, 60 / n),#y limits C upper
                        ncol=2, byrow=FALSE)

    limitB2L <- matrix (data= c (50 / n, 50 / n, 90 / n, 330 / n, maxX, maxX, 50 / n, #x limits B lower
                                 0, 30 / n, 80 / n, 230 / n , .endy (330 / n, 230 / n, maxX, cbl), 0, 0),#y limits B lower
                        ncol=2, byrow=FALSE)

    limitB2U <- matrix (data= c (0, 30 / n, 230 / n, .endx (230 / n, 330 / n, maxY, cbu), 0, 0, #x limits B upper
                                 50 / n, 50 / n, 330 / n , maxY, maxY, 50 / n),#y limits B upper
                        ncol=2, byrow=FALSE)

    #labelling zones using in.out function from mgcv package
    zones[which (in.out (limitB2L, testdf))] <- "B"
    zones[which (in.out (limitB2U, testdf))] <- "B"
    zones[which (in.out (limitC2L, testdf))] <- "C"
    zones[which (in.out (limitC2U, testdf))] <- "C"
    zones[which (in.out (limitD2L, testdf))] <- "D"
    zones[which (in.out (limitD2U, testdf))] <- "D"
    zones[which (in.out (limitE2, testdf))] <- "E"
  }
  return (zones)
}


#finding tan coefficient for a line going trough 2 points
.coef <- function (x, y, xend, yend){
  if (xend==x) {
    stop("Vertical line - function inapplicable")
  }
  return ( (yend - y) / (xend - x))
}

#setting y axis end point with tan coeff and x endpoint known
.endy <- function (startx, starty, maxx, coef){
  return ( (maxx - startx) * coef + starty)
}

#setting x axis end point with tan coeff and y endpoint known
.endx <- function (startx, starty, maxy, coef){
  return ( (maxy - starty)  / coef + startx)
}
