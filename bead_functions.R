#' Plot helpful cytograms for estimating the D1, D2 and FSC coordinates of the inflection point (corresponds to location of 1µm beads).
#'
#' @param dataframe containing EVT data.
#' @return D1, D2 and fsc values of presumed 1 µm beads
#' @export
inflection.point <- function(DF) {
  QUANTILES <- c(2.5, 50.0, 97.5)

  def_par <- par(no.readonly = TRUE) # save default, for resetting...
  par(mfrow = c(1,3), pty = "s")

  popcycle::plot_cyt(DF, "fsc_small", "pe")

  poly_beads <- splancs::getpoly(quiet = TRUE)
  points <- DF[, c("fsc_small", "pe")]
  colnames(points) <- c("x", "y")
  b <- subset(DF, splancs::inout(points, poly = poly_beads, bound = TRUE, quiet = TRUE))

  popcycle::plot_cyt(b, "fsc_small", "D1")
  abline(h = 29000, lwd = 1, col = "red3")
  abline(v = 44500, lwd = 1, col = "red3")

  polyd1 <- splancs::getpoly(quiet=TRUE)
  points <- b[, c("fsc_small", "D1")]
  colnames(points) <- c("x", "y")
  opp.d1 <- subset(b,splancs::inout(points, poly = polyd1, bound = TRUE, quiet = TRUE))

  popcycle::plot_cyt(b, "fsc_small", "D2")
  abline(h = 29000, lwd = 1, col = "red3")
  abline(v = 44500, lwd = 1, col = "red3")

  polyd2 <- splancs::getpoly(quiet=TRUE)
  points <- b[, c("fsc_small", "D2")]
  colnames(points) <- c("x", "y")
  opp.d2 <- subset(b, splancs::inout(points, poly = polyd2, bound = TRUE, quiet = TRUE))

  FSC <- round(summary(c(opp.d1$fsc_small, opp.d2$fsc_small)))
  D1 <- round(summary(opp.d1$D1))
  D2 <- round(summary(opp.d2$D2))

  inflection <- data.frame()
  for (quant in QUANTILES) {
    if (quant == 2.5) {
      i <- 2; j <- 5
    } else if (quant == 50.0) {
      i <- j <- 3
    } else if (quant == 97.5) {
      i <- 5; j <- 2
    }
    fsc <- as.vector(FSC[i])
    d1 <- as.vector(D1[j])
    d2 <- as.vector(D2[j])
    newrow <- data.frame(quantile = quant, fsc, d1, d2, stringsAsFactors = FALSE)
    inflection <- rbind(inflection, newrow)
  }
  par(def_par)

  return(inflection)
}


create.filter.params <- function(inst, fsc, d1, d2, min.d1, min.d2, width = width, slope.file = NULL) {
  QUANTILES <- c(2.5, 50.0, 97.5)

  # Rename to get correct dataframe headers
  beads.fsc.small <- as.numeric(fsc)
  beads.D1 <- as.numeric(d1)
  beads.D2 <- as.numeric(d2)
  min.D1 <- as.numeric(min.d1)
  min.D2 <- as.numeric(min.d2)

  width <- as.numeric(width)

  if (is.null(slope.file)) {
    slope.file <- "https://raw.githubusercontent.com/armbrustlab/seaflow-virtualcore/master/1.bead_calibration/seaflow_filter_slopes.csv"
  }
  slopes <- read.csv(slope.file)

  filter.params <- data.frame()
  headers <- c(
    "quantile", "beads.fsc.small",
    "beads.D1", "beads.D2", "width",
    "notch.small.D1", "notch.small.D2",
    "notch.large.D1", "notch.large.D2",
    "offset.small.D1", "offset.small.D2",
    "offset.large.D1", "offset.large.D2"
  )
  for (quant in QUANTILES) {
    if (quant == 2.5) {
      suffix <- "_2.5"
      i <- 1
    } else if (quant == 97.5) {
      suffix <- "_97.5"
      i <- 3
    } else if (quant == 50.0) {
      suffix <- ""
      i <- 2
    }

    # Small particles
    offset.small.D1 <- min.D1
    offset.small.D2 <- min.D2
    notch.small.D1 <- round((beads.D1[i]-min.D1)/beads.fsc.small[i],3)
    notch.small.D2 <- round((beads.D2[i]-min.D2)/beads.fsc.small[i],3)

    # Large particles
    notch.large.D1 <- round(slopes[slopes$ins == inst, paste0('notch.large.D1', suffix)], 3)
    notch.large.D2 <- round(slopes[slopes$ins == inst, paste0('notch.large.D2', suffix)], 3)
    offset.large.D1 <- round(beads.D1[i] - notch.large.D1 * beads.fsc.small[i])
    offset.large.D2 <- round(beads.D2[i] - notch.large.D2 * beads.fsc.small[i])

    newrow <- data.frame(
      quant, beads.fsc.small[i],
      beads.D1[i], beads.D2[i], width,
      notch.small.D1, notch.small.D2,
      notch.large.D1, notch.large.D2,
      offset.small.D1, offset.small.D2,
      offset.large.D1, offset.large.D2,
      stringsAsFactors=FALSE
    )
    names(newrow) <- headers
    filter.params <- rbind(filter.params, newrow)
  }

  return(filter.params)
}