#' Provides Error Grid Analysis to aid clinical interpretation of blood glucose
#' meter results.
#'
#' @name ega-package
#' @import ggplot2
#' @docType package
NULL


#' @export
#' @title Assign Clarke zones to paired glucose readings.
#' @description Description goes here.
#' @param reference_vals A vector of glucose values obtained via the reference
#' method.
#' @param test_vals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{reference_vals}, so the length should be the same.
#' @return A character vector is returned, with each element being one of
#' \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, or \code{"E"}.
#' @references
#' Clarke, W. L., D. Cox, L. A. Gonder-Frederick, W. Carter, and S. L. Pohl.
#' “Evaluating Clinical Accuracy of Systems for Self-Monitoring of Blood Glucose.”
#' Diabetes Care 10, no. 5 (September 1, 1987): 622–28.
#' doi:10.2337/diacare.10.5.622.
getClarkeZones <- function(reference_vals, test_vals) {

  zones <- vector(mode="character", length = length(reference_vals))

  bias <- test_vals-reference_vals

  # absolute relative error = abs(bias)/reference*100
  are <- abs(bias)/reference_vals*100

  eq1 <- (7/5)*(reference_vals-130)
  eq2 <- reference_vals+110

  # zone A: are <= 20  or (ref < 70 and test < 70)
  zoneA <- (are <= 20) | (reference_vals < 70 & test_vals < 70)

  zones[zoneA] <- "A"

  # zone E: (ref <= 70 and test >= 180) or (ref >=180 and test <=70)
  zoneE <- (reference_vals <= 70 & test_vals >= 180) |
    (reference_vals >= 180 & test_vals <= 70)

  zones[zoneE] <- "E"

  # zone D: ref < 70 and (test > 70 and test < 180) or
  #   ref > 240 and (test > 70 and test < 180)
  test_ok <- test_vals > 70 & test_vals < 180
  zoneD <- (reference_vals < 70 & test_ok) | (reference_vals > 240 & test_ok)

  zones[zoneD] <- "D"

  # zone C: (ref >= 130 and ref <= 180 and test < eq1) or
  #   (ref > 70 and ref > 180 and ref > eq2)
  zoneC <- (reference_vals >= 130 & reference_vals <= 180 & test_vals < eq1) |
    (reference_vals > 70 & test_vals > 180 & test_vals > eq2)

  zones[zoneC] <- "C"

  # the rest are zone B
  zones <- replace(zones, zones=="", "B")

  return(zones)

}


#' @export
#' @title Assign Parkes zones to paired glucose readings.
#' @description Description goes here.
#' @param reference_vals A vector of glucose values obtained via the reference
#' method.
#' @param test_vals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{reference_vals}, so the length should be the same.
#' @return A character vector is returned, with each element being one of
#' \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, or \code{"E"}.
#' @references
#' Parkes, J. L., S. L. Slatin, S. Pardo, and B. H. Ginsberg. “A New Consensus
#' Error Grid to Evaluate the Clinical Significance of Inaccuracies in the
#' Measurement of Blood Glucose.” Diabetes Care 23, no. 8 (August 2000): 1143–48.
#'
#' Pfützner, Andreas, David C. Klonoff, Scott Pardo, and
#' Joan L. Parkes. “Technical Aspects of the Parkes Error Grid.” Journal of
#' Diabetes Science and Technology 7, no. 5 (September 2013): 1275–81.
getParkesZones <- function(reference_vals, test_vals) {

  zones <- vector(mode="character", length = length(reference_vals))

  # B lower lines:
  # 50/0->50/30->170/145->385/300->5000/4495.45
  bl_x <- c(50, 50, 170, 385, 5000)
  bl_y <- c(0, 30, 145, 300, 4495.45)

  # B lower line equations
  bl_lines <- getLineEqs(bl_x, bl_y, reference_vals)

  bl1 <- bl_lines[[1]]
  bl2 <- bl_lines[[2]]
  bl3 <- bl_lines[[3]]
  bl4 <- bl_lines[[4]]

  # B upper lines:
  # 0/50->30/50->140/170->280/380->5000/5729.3
  bu_x <- c(0, 30, 140, 280, 5000)
  bu_y <- c(50, 50, 170, 380, 5729.3)

  # B upper line equations
  bu_lines <- getLineEqs(bu_x, bu_y, reference_vals)

  bu1 <- bu_lines[[1]]
  bu2 <- bu_lines[[2]]
  bu3 <- bu_lines[[3]]
  bu4 <- bu_lines[[4]]

  # C lower lines:
  # 120/0->120/30->260/130->5000/2091.38
  cl_x <- c(120, 120, 260, 5000)
  cl_y <- c(0, 30, 130, 2091.38)

  # C lower line equations:
  cl_lines <- getLineEqs(cl_x, cl_y, reference_vals)

  cl1 <- cl_lines[[1]]
  cl2 <- cl_lines[[2]]
  cl3 <- cl_lines[[3]]

  # C upper lines:
  # 0/60->30/60->50/80->70/110->5000/11526.84
  cu_x <- c(0, 30, 50, 70, 5000)
  cu_y <- c(60, 60, 80, 110, 11526.84)

  # C upper line equations:
  cu_lines <- getLineEqs(cu_x, cu_y, reference_vals)

  cu1 <- cu_lines[[1]]
  cu2 <- cu_lines[[2]]
  cu3 <- cu_lines[[3]]
  cu4 <- cu_lines[[4]]

  # D lower lines:
  # 250/0->250/40->5000/1781.67
  dl_x <- c(250, 250, 5000)
  dl_y <- c(0, 40, 1781.67)

  # D lower line equations:
  dl_lines <- getLineEqs(dl_x, dl_y, reference_vals)

  dl1 <- dl_lines[[1]]
  dl2 <- dl_lines[[2]]

  # D upper lines:
  # 0/100->25/100->50/125->80/215->5000/36841.67
  du_x <- c(0, 25, 50, 80, 5000)
  du_y <- c(100, 100, 125, 215, 36841.67)

  # D upper line equations
  du_lines <- getLineEqs(du_x, du_y, reference_vals)

  du1 <- du_lines[[1]]
  du2 <- du_lines[[2]]
  du3 <- du_lines[[3]]
  du4 <- du_lines[[4]]

  # E lines:
  # 0/150->35/155->5000/130900
  eu_x <- c(0, 35, 5000)
  eu_y <- c(150, 155, 130900)

  # E line equations
  eu_lines <- getLineEqs(eu_x, eu_y, reference_vals)

  eu1 <- eu_lines[[1]]
  eu2 <- eu_lines[[2]]

  # B lower lines:
  # 50/0->50/30->170/145->385/300->5000/4495.45
  zoneB_lower <- (test_vals > 0 & test_vals <= 30 & reference_vals >= 50) |
    (test_vals > 30 & test_vals <= 145 & test_vals < bl2) |
    (test_vals > 145 & test_vals <= 300 & test_vals < bl3) |
    (test_vals > 300 & test_vals <= 4495.45 & test_vals < bl4)

  # B upper lines:
  # 0/50->30/50->140/170->280/380->5000/5729.3
  zoneB_upper <- (reference_vals > 0 & reference_vals <= 30 & test_vals > 50) |
    (test_vals > 50 & test_vals <= 170 & test_vals > bu2) |
    (test_vals > 170 & test_vals <= 380 & test_vals > bu3) |
    (test_vals > 380 & test_vals <= 5729.3 & test_vals > bu4)

  # C lower lines:
  # 120/0->120/30->260/130->5000/2091.38
  zoneC_lower <- (test_vals > 0 & test_vals <= 30 & reference_vals >= 120) |
    (test_vals > 30 & test_vals <= 130 & test_vals < cl2) |
    (test_vals > 130 & test_vals <= 2091.38 & test_vals < cl3)

  # C upper lines:
  # 0/60->30/60->50/80->70/110->5000/11526.84
  zoneC_upper <- (reference_vals > 0 & reference_vals <= 30 & test_vals > 60) |
    (test_vals > 60 & test_vals <= 80 & test_vals > cu2) |
    (test_vals > 80 & test_vals <= 110 & test_vals > cu3) |
    (test_vals > 110 & test_vals <= 11526.84 & test_vals > cu4)

  # D lower lines:
  # 250/0->250/40->5000/1781.67
  zoneD_lower <- (test_vals > 0 & test_vals <= 40 & reference_vals >= 250) |
    (test_vals > 40 & test_vals <= 1781.67 & test_vals < dl2)

  # D upper lines:
  # 0/100->25/100->50/125->80/215->5000/36841.67
  zoneD_upper <- (reference_vals > 0 & reference_vals <= 25 & test_vals > 100) |
    (test_vals > 100 & test_vals <= 125 & test_vals > du2) |
    (test_vals > 125 & test_vals <= 215 & test_vals > du3) |
    (test_vals > 215 & test_vals <= 36841.67 & test_vals > du4)

  # E lines:
  # 0/150->35/155->5000/130900
  zoneE_upper <- (reference_vals > 0 & reference_vals <= 35 & test_vals > eu1) |
    (test_vals > 155 & test_vals <= 130900 & test_vals > eu2)


  zones[zoneB_lower] <- "B"
  zones[zoneC_lower] <- "C"
  zones[zoneD_lower] <- "D"

  zones[zoneB_upper] <- "B"
  zones[zoneC_upper] <- "C"
  zones[zoneD_upper] <- "D"

  zones[zoneE_upper] <- "E"

  # the rest are zone A
  zones <- replace(zones, zones=="", "A")

  zones

}


#' @export
#' @title Generate "reference" and "test" glucose data.
#' @description Description goes here.
generateGlucoseData <- function(n=100, precision=0.2, lower=90, upper=130) {

    nl_range <- lower:upper

    reference_vals <- sample(nl_range, n, TRUE)

    generateTestValues <- function(reference_val, precision) {

      # possible percentages by which to alter the test value
      precision_seq <- seq(0, precision, by=0.01)

      dir <- sample(c(1,-1), 1)
      amt <- sample(precision_seq, 1) * dir

      test_val <- reference_val + (amt*reference_val)

      return(test_val)

    }

    test_vals <- sapply(reference_vals, generateTestValues, precision)

    list("reference_vals"=reference_vals, "test_vals"=test_vals)

}


#' @title Get line equations.
#' @description Description goes here.
getLineEqs <- function(xcoords, ycoords, xvals) {

  num_lines <- length(xcoords)-1
  line_eqs <- vector(mode="list", length=num_lines)

  for (i in 1:num_lines) {

    x1 <- xcoords[i]
    y1 <- ycoords[i]
    x2 <- xcoords[i+1]
    y2 <- ycoords[i+1]

    m <- (y2-y1)/(x2-x1)
    b <- (-x1*m)+y1

    line_eqs[[i]] <- m*xvals+b

  }

  line_eqs

}
