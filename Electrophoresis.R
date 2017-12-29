#!/usr/bin/env Rscript
# Electrophoresis.R
# Creates a scatterplot of the molecular weight ladder with Rf values vs log(MW).
# An interpolated line is drawn with its equation shown. Then for each lane, 
# draw a line for each band with its indicative Rf value until it intersects the
# line of best fit whereby a perpendicular line is drawn until it intersects the
# logarithmic axis.
#
# Currently reads from two files: MolecularLadder.csv with "Rf,MW" as entries
# and Lanes.txt with each lane having a new line with its Rf values separated by
# commas.
#
# You are free to make your own aesthetic changes. Scroll down to line ~180 to
# find the main section.
#
# Author: Matt Preston (website: matthewpreston.github.io)
# Created On: Nov ?? 2016	V1.0
# Revised On: Dec 24 2016	V1.1 - Added dashed, multicoloured interpolation
#								   lines, extrapolation dashes, improved legend,
#								   different drawing order, 'straightened' the
#								   exponential equation
#             Dec 27 2016   V1.2 - Added command line support, added subtitle
#								   feature, titles can now be bold/italics/etc.

# ==== INCLUDES ================================================================

# Used for parsing command line arguments
if (!suppressWarnings(require(optparse, quietly=TRUE))) {
  install.packages("optparse", dependencies=TRUE)
  suppressPackageStartupMessages(library(optparse))
}

# Used for subscripting/superscripting legend entries
if (!suppressWarnings(require(latex2exp, quietly=TRUE))) {
  install.packages("latex2exp", dependencies=TRUE)
  suppressPackageStartupMessages(library(latex2exp))
}

# ==== FUNCTIONS ===============================================================

# Sees if file names exists
#
# @param files		A character vector of file names
# @stop				Stops if a file name does not exist
checkFiles <- function(files) {
  for (file in files) {
    if (file.access(file) == -1) {
	  stop(sprintf("Specified file ( %s ) does not exist", file))
	}
  }
}

# Reads from a text file where each lane is on a separate line. The first entry
# is what substances were run in the lanes with the rest of the entries being Rf
# values. Each entry is separated by a comma. The title for each lanes can be
# made superscripted/subscripted if the syntax for LaTex is correct. (See:
# ftp://cran.r-project.org/pub/R/web/packages/latex2exp/vignettes/using-latex2exp.html)
#
# @param textFile   Name of the file to read from
# @param header     Whether the file has a header line (default=TRUE)
# @return           Returns a list of numerical vectors (each lane is a vector)
loadLanes <- function(textFile, header=TRUE) {
  input <- readLines(textFile)
  if (header) {
    input <- input[-1]
  }
  lanes = list()
  for (lane in input) {
    line = strsplit(lane, ",")[[1]]
    lanes[[line[1]]] = as.numeric(line[-1])
  }
  return(lanes)
}

# Returns either the maximum or minimum y limit when determining the range of y
# values to plot the scatter graph.
#
# @param ladder     Data frame with Rf and MW (use code above to create ladder)
# @param lanes      List of vectors containing Rf values for each band
# @param comparison Find either minimum or maximum (default = min)
# @return           Returns the minimum/maximum y limit
findRfLimit <- function(ladder, lanes, comparison=min) {
  # Create a single vector of Rf values
  RfValues <- ladder$Rf
  for (lane in lanes) {
    RfValues <- c(RfValues, lane)
  }
  # Using the given function, return the min/max y value
  return(comparison(RfValues))
}

# Generates a vector of suggested tick marks when using a log10 scale. 
# Ex. When n=9 there will be 9 tick marks between [ 10^x, 10^(x+1) )
# The parameter n must be either 1, 2, 5, or 9 as these were 'relatively' spaced
#
# @param tickMin    The lower bound
# @param tickMax    The upper bound
# @param n          The amount of ticks between [ 10^x, 10^(x+1) )
#                   (default = 9)
# @stop				Stops if n is not 1, 2, 5, or 9
# @return           Returns a vector of tick marks to be used in axis()
getLogTicks <- function(tickMin, tickMax, n=9) {
  if (!(n %in% c(1,2,5,9))) { stop("n is not 1, 2, 5, or 9") }
  # Find lowest tick mark based off of tickMin (a*10^p >= tickMin)
  lowestPower <- floor(log10(tickMin))
  lowest <- 10 ^ lowestPower
  a = 1
  while (a*lowest < tickMin) { a = a + 1 }
  lowest = a * lowest
  # Find highest tick mark based off of tickMax (a*10^p <= tickMax)
  highestPower <- floor(log10(tickMax))
  highest <- 10 ^ highestPower
  a = 1
  while ((a+1)*highest <= tickMax) { a = a + 1 }
  highest = a * highest
  # Create a vector that has all the tick marks from 10 ^ lowest power to
  # 10 ^ highest power first
  if (n == 9) {
    significands = c(1,2,3,4,5,6,7,8,9)
  } else if (n == 5) {
    significands = c(1,2,4,6,8)
  } else if (n == 2) {
    significands = c(1,5)
  } else {
    significands = c(1)
  }
  exponentials <- 10^(lowestPower:highestPower)
  result <- c()
  for (exponential in exponentials) {
    result <- c(result, significands*exponential)
  }
  # Filter out the tick marks that are not within [lowest,highest] and return
  result <- result[result >= lowest]
  return(result[result <= highest])
}

# Provides prettier axis tick marks (Ex. 2e+05 -> 2x10^5, but with superscript)
#
# @param logTicks   Output from getLogTicks
# @return           Returns labels that can be used in axis()
getLogLabels <- function(logTicks) {
  exponents <- floor(log10(logTicks))
  sigs <- logTicks / (10 ^ exponents)
  result <- c()
  for (i in 1:length(logTicks)) {
    result <- c(result, as.expression(bquote(.(sigs[i])%*%10^.(exponents[i]))))
  }
  return(result)
}

# Turns a number into an exponential string which LaTeX can interpret to
# superscript the power
# Ex. (sigDigits=2): 206971 -> "$2.1\\times10^5$" (2.1x10^5)
#
# @param number     Number to transmute into an exponential
# @param sigDigits  Number of significant digits to include (default=2)
# @return           Returns a exponential string safe for LaTeX
exponentialLaTeX <- function(number, sigDigits=2) {
  exponent <- floor(log10(number))
  number <- number / (10 ^ exponent)
  return(sprintf("$%1.*f\\times 10^%d$", sigDigits-1, number, exponent))
}

# Adds some modifiers (ex bold, italic, etc.) when creating a LaTeX expression.
# Typically, TeX returns an expression(), which is used when displaying
# complicated text expressions. However, one can't easily add 'bold()' and such,
# as used in plot titles. This function mediates this deficiency.
#
# @param string     A string to be used for LaTeX formatting
# @param modifiers  A character vector of modifiers to be added
#                   See the R manual for mathematic expressions
# @return           Returns a modified expression to be used in displaying text
modifiedLaTeXExpression <- function(string, modifiers=c()) {
  s <- paste(deparse(as.list(TeX(string))[[1]]), collapse="")
  for (modifier in modifiers) {
    s <- paste(as.character(modifier), "(", s, ")", sep="")
  }
  s <- parse(text=s)
  substitute(s)
}

# ==== MAIN ====================================================================

# Parse argument options
option_list <- list(
  make_option(c("-o", "--output"), default="Molecular Weight Interpolation.png",
			  help="Name out output file [%default]", metavar="FILE"),
  make_option(c("-m", "--main"),
			  default=paste("Interpolation of Molecular Weights from Standard ",
							"Molecular Weight Ladder", sep=""),
			  help="Title of plot. Supports LaTeX formatting. [%default]",
			  metavar="STR"),
  make_option(c("-s", "--subtitle"),
			  help="Subtitle of plot. Supports LaTeX formatting.",
			  metavar="STR"),
  make_option(c("-x", "--x_label"), default="$R_f$ Value",
			  help="X axis label. Supports LaTeX formatting. [%default]",
			  metavar="STR"),
  make_option(c("-y", "--y_label"), default="Molecular Weight (Da)",
			  help="Y axis label. Supports LaTeX formatting. [%default]",
			  metavar="STR"),
  make_option(c("--width"), type="integer", default=960,
			  help="Width of output file in pixels [%default]", metavar="INT"),
  make_option(c("--height"), type="integer", default=540,
			  help="Height of output file in pixels [%default]", metavar="INT"),
  make_option(c("--resolution"), type="integer", default=72,
			  help="Resolution of output file in pixels per inch [%default]",
			  metavar="INT")
  )
parser <- OptionParser(usage="%prog [options] Ladder.csv Lanes.txt",
					   option_list=option_list)
arguments <- tryCatch({
  parse_args(parser, positional_arguments=2)
}, error=function(err) {
  cat("\nImproper number of positional arguments, 2 required\n\nTraceback:\n")
  stop(err)
})
opts <- arguments$options
args <- arguments$args
checkFiles(args)

# Load molecular weight ladder and data files
ladder = read.csv(args[1], header=TRUE, stringsAsFactors=FALSE)
names(ladder) <- c("Rf", "MW")
lanes = loadLanes(args[2], header=TRUE)

# Output file
png(file=opts$output, width=opts$width, height=opts$height, res=opts$resolution)

# Find linear regression line
model = lm(log(ladder$MW)~ladder$Rf)
slope = model$coefficients[[2]]
yIntercept = model$coefficients[[1]]

# Prepare for plotting
opar <- par()  # Store current state
par(mar=c(5,5,4,1))
plotMargins <- par()$mar
plotWidth <- opts$width - (plotMargins[2]+plotMargins[4])/5*opts$resolution
plotHeight <- opts$height - (plotMargins[1]+plotMargins[3])/5*opts$resolution

# Now to begin plotting. Steps are performed in this order:
#
# 1) Draw plot with labelled axes and a title. Ladder data are not drawn yet.
# 2) Find distinct bands and draw interpolated, multicoloured lines for each
#	 one of these bands
# 3) Add a corresponding Rf value and interpolated MW for each of these bands
# 4) Now plot the ladder points with its regression line
# 5) Print the exponential equation used being parallel to the line
# 6) Finally, add a legend to explain everything
#
# The reason for this order is purely for looks (a tradeoff for potential
# efficiency)

# Draw out plot skeleton
xMin <- findRfLimit(ladder, lanes, min)
xMax <- findRfLimit(ladder, lanes, max)
xAxisMin <- floor(10*xMin)/10
xAxisMax <- ceiling(10*xMax)/10
yMin <- exp(slope*xMax+yIntercept) # Reverse maximum and minimum
yMax <- exp(slope*xMin+yIntercept) # since of negative slope
yAxisMin <- yMin / 1.8 # Add a small buffer room for aesthetics
yAxisMax <- yMax * 1.8
plot(ladder$Rf, ladder$MW, log="y",
     xlim=c(xAxisMin, xAxisMax), xaxs="i", xlab=NA, xaxt="n",
     ylim=c(yAxisMin, yAxisMax), yaxs="i", ylab=NA, yaxt="n",
     cex.main=1.5, type="n")
if (!is.null(opts$subtitle)) {					# Add main title (and subtitle)
  mtext(modifiedLaTeXExpression(opts$main, "bold"), side=3, line=2, cex=1.5)
  mtext(modifiedLaTeXExpression(opts$subtitle, "bold"), side=3,line=0.5,cex=1.5)
} else {
  mtext(modifiedLaTeXExpression(opts$main, "bold"), side=3, line=1.5, cex=1.5)
}
xTicks <- seq(xAxisMin, xAxisMax, by=0.1)
axis(1, at=xTicks, labels=xTicks, las=1, cex.axis=1)			# Add x axis
mtext(TeX(opts$x_label), side=1, line=3, cex=1)
yTicks <- getLogTicks(yAxisMin, yAxisMax, n=9)
yLabels <- getLogLabels(yTicks)
axis(2, at=yTicks, labels=yLabels, las=1, cex.axis=1, hadj=0.85)# Add y axis
mtext(TeX(opts$y_label), side=2, line=3.5, cex=1)

# Find distinct bands
distinctBands <- list()
for (i in 1:length(lanes)) {
  for (band in lanes[[i]]) {
    band <- as.character(band)
    distinctBands[[band]] <- c(distinctBands[[band]], i)
  }
}

# Draw these bands using dashed, multicoloured lines
# Make each 'dash' be about 1/6 of an inch; this was chosen aesthetically
xDivision <- opts$resolution / 6 / plotWidth * (xAxisMax-xAxisMin)
yDivision <- opts$resolution / 6 / plotHeight * log10(yAxisMax/yAxisMin)
colours <- rainbow(length(lanes))
for (band in names(distinctBands)) {
  RfValue <- as.numeric(band)
  yIntersection <- exp(slope*RfValue+yIntercept)
  # Draw a straight horizontal line originating from the intersection
  # point to the y-axis
  x0Segments <- seq(xAxisMin, RfValue, xDivision)
  x1Segments <- c(x0Segments[-1], RfValue)
  segments(x0=x0Segments, y0=yIntersection,
           x1=x1Segments, y1=yIntersection,
           col=colours[distinctBands[[band]]], lwd=2, lty=2)
  # Draw a straight vertical line originating from the Rf value to the
  # regression line
  y0Segments <- 10 ^ seq(log10(yAxisMin), log10(yIntersection), yDivision)
  y1Segments <- c(y0Segments[-1], yIntersection)
  segments(x0=RfValue, y0=y0Segments,
           x1=RfValue, y1=y1Segments,
           col=colours[distinctBands[[band]]], lwd=2, lty=2)
}

# Add a label showing interpolated MW and original Rf value
# I do iterate through this loop again to get the text on top of the
# interpolated lines
for (band in names(distinctBands)) {
  RfValue <- as.numeric(band)
  yIntersection <- exp(slope*RfValue+yIntercept)
  text(xAxisMin, yIntersection,
       TeX(exponentialLaTeX(yIntersection, sigDigits = 3)),
       cex=0.8, adj=c(-0.03, -0.15))
  text(RfValue, yAxisMin, sprintf("%.6f", RfValue),
       cex=0.8, srt=270, adj=c(1, -0.5))
}

# Plot the molecular ladder points and linear regression line
points(ladder$Rf, ladder$MW, pch=19, cex=1.5)
ladderRfMin <- min(ladder$Rf)
ladderRfMax <- max(ladder$Rf)
ladderMWMin <- exp(slope*ladderRfMax+yIntercept) # Reverse maximum and minimum
ladderMWMax <- exp(slope*ladderRfMin+yIntercept) # because of negative slope
point1 <- c(xMin, exp(slope*xMin+yIntercept)) # Leftmost point (lowest Rf)
point2 <- c(xMax, exp(slope*xMax+yIntercept)) # Rightmost point (highest Rf)
lines(c(ladderRfMin, ladderRfMax), c(ladderMWMax, ladderMWMin), lwd=2)
if (point1[1] < ladderRfMin) { # Extrapolate leftward
  lines(c(point1[1], ladderRfMin), c(point1[2], ladderMWMax), lwd=1, lty=2)
}
if (point2[1] > ladderRfMax) { # Extrapolate rightward
  lines(c(point2[1], ladderRfMax), c(point2[2], ladderMWMin), lwd=1, lty=2)
}

# Plot margins can be given in terms of the number of lines. To convert this to
# pixels (to be used in trigonometry...), we need to know how many pixels
# comprise a line. PNG resolution is given in ppi or pixel per inch. And
# conveniently R can tell you how many inches a line takes up. Use:
#
# par()$mai/par()$mar
#
# and you'll see that there are 0.2 inches per line. Then using math, we can
# figure out the width and height of a plot in pixels followed by finding the
# pixels covered by the regression line and finally apply arctangent. All of
# this just to find a convenient unit to apply trigonometry in order to make
# this one silly piece of text be parallel with the regression line!
deltaX <- point2[1] - point1[1]
deltaY <- log10(point2[2]/point1[2])
deltaX <- deltaX / (xAxisMax - xAxisMin) * plotWidth
deltaY <- deltaY / log10(yAxisMax/yAxisMin) * plotHeight

# Plot equation of regression line
meanX <- mean(c(point1[1], point2[1]))
coefficient <- exp(yIntercept)
text(meanX*1.05, coefficient*exp(slope*meanX)*1.05,
     TeX(sprintf("$MW=%de^{%.3fR_f}$", round(coefficient), slope)),
     #as.expression(bquote(MW==.(coefficient)*e^(.(slope)*R[f]))),
     cex=1.5, srt=atan(deltaY/deltaX)/pi*180)

# Put up a legend to explain the 'pretty' colours
legend("topright",
       c("Molecular Weight Ladder", "Regression Line", TeX(names(lanes)),
         "Multicoloured Lines Indicate Lane Overlap"),
       col=c("black", "black", colours, "black"),
       lty=c(NA, rep(1, length(names(lanes))+1), NA),
       pch=c(19, rep(NA, length(names(lanes))+1), NA),
       bty="n", pt.cex=1.1, lwd=2, cex=1)

# Finish off
suppressWarnings(par(opar))
garbage <- dev.off()