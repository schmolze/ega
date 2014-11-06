#' @export
#' @title Plot a Clarke Error Grid for the given data.
#' @description Description goes here.
#' @param reference_vals A vector of glucose values obtained via the reference
#' method.
#' @param test_vals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{reference_vals}, so the length should be the same.
#' @param zones An optional character vector specifying the Clarke zones
#' for each paired value. If this is not supplied, \code{\link{getClarkeZones}}
#' will be called to generate zone labels.
#' @references
#' Clarke, W. L., D. Cox, L. A. Gonder-Frederick, W. Carter, and S. L. Pohl.
#' “Evaluating Clinical Accuracy of Systems for Self-Monitoring of Blood Glucose.”
#' Diabetes Care 10, no. 5 (September 1, 1987): 622–28.
#' doi:10.2337/diacare.10.5.622.
plotClarkeGrid <- function(reference_vals, test_vals, zones=NA) {

  # use default zone assignment if none is provided
  if (is.na(zones))
    zones <- getClarkeZones(reference_vals, test_vals)

  tolerance <- 0.2

  # create a df for ggplot
  data <- data.frame("ref"=reference_vals, "test"=test_vals, "zones"=zones)

  zoneA <- subset(data, zones=="A")
  maxA <- max(zoneA[["ref"]])

  # calculate line for zone A upper
  x1 <- 70
  y1 <- 84
  x2 <- maxA
  y2 <- maxA+tolerance*maxA
  slope_Au <- (y2-y1)/(x2-x1) # rise over run!
  intercept_Au <- (-70*slope_Au)+84

  # choose large x,y end values so line goes off the chart
  xend_Au <- 5000
  yend_Au <- slope_Au*5000+intercept_Au # y=mx+b

  # calculate line for zone A lower
  x1 <- 70
  y1 <- 56
  x2 <- maxA
  y2 <- maxA-tolerance*maxA
  slope_Al <- (y2-y1)/(x2-x1)
  intercept_Al <- (-70*slope_Al)+56

  # choose large x,y end values so line goes off the chart
  xend_Al <- 5000
  yend_Al <- slope_Al*5000+intercept_Al

  # zone C upper
  x1 <- 70
  y1 <- 180
  slope_Cu <- 1
  intercept_Cu <- (-70*slope_Cu)+180

  # intersection of C upper and A upper
  # C upper: y = ax+c
  # A upper: y = bx+d

  # intersection: (d-c)/(a-b), (ad-bc)/(a-b)
  xend_Cu <- (intercept_Au-intercept_Cu)/(slope_Cu-slope_Au)
  yend_Cu <- ((slope_Cu*intercept_Au)-(slope_Au*intercept_Cu))/(slope_Cu-slope_Au)

  ceg <- ggplot(data, aes(x=ref, y=test)) +

    # label the clinially relevant levels
    scale_x_continuous(breaks=c(70, 100, 150, 180, 240, 300, 350, 400, 450,
                                         500, 550, 600, 650, 700, 750, 800,
                                         850, 900, 950, 1000),
                       expand = c(0,0)) +

    scale_y_continuous(breaks=c(70, 100, 150, 180, 250, 300, 350, 400, 450,
                                         500, 550, 600, 650, 700, 750, 800,
                                         850, 900, 950, 1000),
                       expand=c(0, 0)) +

    # color by location type
    geom_point(aes(color=zones)) +

    # zone E upper
    annotate("segment", x=0, y=180, xend=70, yend=180, alpha=0.6) +
    annotate("segment", x=70, y=180, xend=70, yend=Inf, alpha=0.6) +

    # zone E lower
    annotate("segment", x=180, y=0, xend=180, yend=70, alpha=0.6) +
    annotate("segment", x=180, y=70, xend=Inf, yend=70, alpha=0.6) +

    #zone D right
    annotate("segment", x=240, y=70, xend=240, yend=180, alpha=0.6) +
    annotate("segment", x=240, y=180, xend=Inf, yend=180, alpha=0.6) +

    #zone D left
    annotate("segment", x=0, y=180, xend=70, yend=180, alpha=0.6) +
    annotate("segment", x=0, y=70, xend=70, yend=70, alpha=0.6) +
    annotate("segment", x=70, y=70, xend=70, yend=180, alpha=0.6) +

    # zone A upper
    annotate("segment", x=70, y=84, xend=xend_Au, yend=yend_Au, alpha=0.6) +

    # zone A lower
    annotate("segment", x=70, y=56, xend=xend_Al, yend=yend_Al, alpha=0.6) +

    # zone A vertical line
    annotate("segment", x=70, y=0, xend=70, yend=56, alpha=0.6) +

    # zone C lower
    annotate("segment", x=130, y=0, xend=180, yend=70, alpha=0.6) +

    # zone C upper
    annotate("segment", x=70, y=180, xend=xend_Cu, yend=yend_Cu, alpha=0.6) +

    # now add the zone text labels
    annotate("text", x = 240, y = 230, size=6, label = "A") +
    annotate("text", x = 120, y = 185, size=6, label = "B") +
    annotate("text", x = 350, y = 230, size=6, label = "B") +
    annotate("text", x = 120, y = 300, size=6, label = "C") +
    annotate("text", x = 163, y = 20, size=6, label = "C") +
    annotate("text", x = 35, y = 130, size=6, label = "D") +
    annotate("text", x = 350, y = 130, size=6, label = "D") +
    annotate("text", x = 35, y = 300, size=6, label = "E") +
    annotate("text", x = 350, y = 35, size=6, label = "E") +

    # dummy values to expand right and top
    annotate("text", x = 0, y = 350, size=6, label = "") +
    annotate("text", x = 400, y = 0, size=6, label = "") +

    theme_bw() +
    theme(legend.position="none")

  ceg

}


#' @export
#' @title Plot a Parkes (Consensus) Error Grid for the given data.
#' @description Description goes here.
#' @param reference_vals A vector of glucose values obtained via the reference
#' method.
#' @param test_vals A vector of glucose values obtained via a non-reference
#' method (e.g. a new meter). The values in this vector are paired with those
#' in \code{reference_vals}, so the length should be the same.
#' @param zones An optional character vector specifying the Clarke zones
#' for each paired value. If this is not supplied, \code{\link{getClarkeZones}}
#' will be called to generate zone labels.
#' @references
#' Parkes, J. L., S. L. Slatin, S. Pardo, and B. H. Ginsberg. “A New Consensus
#' Error Grid to Evaluate the Clinical Significance of Inaccuracies in the
#' Measurement of Blood Glucose.” Diabetes Care 23, no. 8 (August 2000): 1143–48.
#'
#' Pfützner, Andreas, David C. Klonoff, Scott Pardo, and
#' Joan L. Parkes. “Technical Aspects of the Parkes Error Grid.” Journal of
#' Diabetes Science and Technology 7, no. 5 (September 2013): 1275–81.
plotParkesGrid <- function(reference_vals, test_vals, zones=NA) {

  # use default zone assignment if none is provided
  if (is.na(zones))
    zones <- getParkesZones(reference_vals, test_vals)

  # create a df for ggplot
  data <- data.frame("ref"=reference_vals, "test"=test_vals, "zones"=zones)

  peg <- ggplot(data, aes(x=ref, y=test)) +

    # label the clinially relevant levels
    scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400,
                                         450, 500, 550, 600, 650, 700, 750, 800,
                                         850, 900, 950, 1000),
                                expand = c(0,0)) +

    scale_y_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300, 350, 400,
                                         450, 500, 550, 600, 650, 700, 750, 800,
                                         850, 900, 950, 1000),
                                expand=c(0, 0)) +

    # color by zone
    geom_point(aes(color=zones)) +

    # zone B upper
    # 0/50->30/50->140/170->280/380->5000/5729.3
    annotate("segment", x=0, y=50, xend=30, yend=50, alpha=0.6) +
    annotate("segment", x=30, y=50, xend=140, yend=170, alpha=0.6) +
    annotate("segment", x=140, y=170, xend=280, yend=380, alpha=0.6) +
    annotate("segment", x=280, y=380, xend=5000, yend=5729.3, alpha=0.6) +

    # zone B lower
    # 50/0->50/30->170/145->385/300->5000/4495.45
    annotate("segment", x=50, y=0, xend=50, yend=30, alpha=0.6) +
    annotate("segment", x=50, y=30, xend=170, yend=145, alpha=0.6) +
    annotate("segment", x=170, y=145, xend=385, yend=300, alpha=0.6) +
    annotate("segment", x=385, y=300, xend=5000, yend=4495.45, alpha=0.6) +

    # zone C upper
    # 0/60->30/60->50/80->70/110->5000/11526.84
    annotate("segment", x=0, y=60, xend=30, yend=60, alpha=0.6) +
    annotate("segment", x=30, y=60, xend=50, yend=80, alpha=0.6) +
    annotate("segment", x=50, y=80, xend=70, yend=110, alpha=0.6) +
    annotate("segment", x=70, y=110, xend=5000, yend=11526.84, alpha=0.6) +

    # zone C lower
    # 120/0->120/30->260/130->5000/2091.38
    annotate("segment", x=120, y=0, xend=120, yend=30, alpha=0.6) +
    annotate("segment", x=120, y=30, xend=260, yend=130, alpha=0.6) +
    annotate("segment", x=260, y=130, xend=5000, yend=2091.38, alpha=0.6) +

    # zone D upper
    # 0/100->25/100->50/125->80/215->5000/36841.67
    annotate("segment", x=0, y=100, xend=25, yend=100, alpha=0.6) +
    annotate("segment", x=25, y=100, xend=50, yend=125, alpha=0.6) +
    annotate("segment", x=50, y=125, xend=80, yend=215, alpha=0.6) +
    annotate("segment", x=80, y=215, xend=5000, yend=36841.67, alpha=0.6) +

    # zone D lower
    # 250/0->250/40->5000/1781.67
    annotate("segment", x=250, y=0, xend=250, yend=40, alpha=0.6) +
    annotate("segment", x=250, y=40, xend=5000, yend=1781.67, alpha=0.6) +

    # zone E upper
    # 0/150->35/155->5000/130900
    annotate("segment", x=0, y=150, xend=35, yend=155, alpha=0.6) +
    annotate("segment", x=35, y=155, xend=5000, yend=130900, alpha=0.6) +

    # zone text labels
    annotate("text", x = 320, y = 320, size=6, label = "A") +
    annotate("text", x = 220, y = 360, size=6, label = "B") +
    annotate("text", x = 385, y = 235, size=6, label = "B") +
    annotate("text", x = 140, y = 375, size=6, label = "C") +
    annotate("text", x = 405, y = 145, size=6, label = "C") +
    annotate("text", x = 415, y = 50, size=6, label = "D") +
    annotate("text", x = 75, y = 383, size=6, label = "D") +
    annotate("text", x = 21, y = 383, size=6, label = "E") +

    # dummy values to force minimum size
    annotate("text", x = 0, y = 550, size=6, label = "") +
    annotate("text", x = 550, y = 0, size=6, label = "") +

    theme_bw() +
    theme(legend.position="none")

  peg

}
