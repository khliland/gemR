% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confints.R
\name{confints}
\alias{confints}
\alias{plot.confints}
\title{Confidence Intervals of Effect Differences}
\usage{
confints(X1, X2, confidence = 0.95, df.used = 0)

\method{plot}{confints}(
  x,
  y,
  xlab = "",
  ylab = "normalised log2",
  sorted = TRUE,
  labels = FALSE,
  nonZero = FALSE,
  xlim = NULL,
  ylim = NULL,
  text.pt = 12,
  ...
)
}
\arguments{
\item{X1}{\code{data.frame} containing first effect.}

\item{X2}{\code{data.frame} containing second effect.}

\item{confidence}{Level of confidence, default = 0.95.}

\item{df.used}{Optional argument indicating how many degrees of freedom have been consumed during deflation. Default = 0.}

\item{x}{Object of class \code{confint}.}

\item{y}{Not used.}

\item{xlab}{X label (\code{character})}

\item{ylab}{Y label (\code{character})}

\item{sorted}{Logical indicating if intervals should be sorted according to their mean values, or a vector of indices/labels to sort by.}

\item{labels}{Logical indicating if sample labels should be used on x axis.}

\item{nonZero}{Logical indicating if intervals are required not to include zero.}

\item{xlim}{Limits of the horizontal scale.}

\item{ylim}{Limits of the vertical scale.}

\item{text.pt}{Size scaling of text in the plot (default = 16).}

\item{...}{Further arguments to \code{qplot}.}
}
\value{
An object of class \code{confints}, which holds
the information needed to perform statistics or plot the
confidence intervals is returned from \code{confints}.
The plotting routine returns a ggplot structure for plotting.
}
\description{
Confidence Intervals of Effect Differences
}
\examples{
data(MS)
# Compare MS and non-MS patients within group 1
conf <- with(MS, confints(proteins[MS == "yes" & group == 1,],
                          proteins[MS == "no"  & group == 1,]))
p1 <- plot(conf)
p2 <- plot(conf, nonZero = TRUE) # Only intervals without 0.
grid.arrange(p1,p2)

# Shorter plot with labels
confShort <- conf[1:10,]
p1 <- plot(confShort, labels = TRUE)
p2 <- plot(confShort, labels = TRUE, nonZero = TRUE)
grid.arrange(p1,p2)
}
\seealso{
\code{\link{GEM}}, \code{\link{elastic}} and \code{\link{pls}}.
}
