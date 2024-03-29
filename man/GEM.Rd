% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gem.R
\name{GEM}
\alias{GEM}
\alias{print.GEM}
\alias{plot.GEM}
\alias{tableGEM}
\title{General Effect Modelling}
\usage{
GEM(formula, data)

\method{print}{GEM}(x, ...)

\method{plot}{GEM}(
  x,
  y = 1,
  what = "raw",
  col = NULL,
  pch = NULL,
  model.line = (what \%in\% c("raw")),
  ylim = NULL,
  ylab = "",
  xlab = "",
  main = NULL,
  ...
)

tableGEM(object, variable)
}
\arguments{
\item{formula}{a model formula specifying features and effects.}

\item{data}{a \code{data.frame} containing response variables (features) and design factors or other groupings/continuous variables.}

\item{x}{Object of class \code{GEM}.}

\item{...}{Additional arguments to \code{plot}}

\item{y}{Response name or number.}

\item{what}{What part of GEM to plot; \code{raw} data (default), \code{fits}, \code{residuals} or
a named model effect (can be combined with 'effect', see \code{Examples}).}

\item{col}{Color of points, defaults to grouping. Usually set to a factor name or a column name in the input data with custom colours.}

\item{pch}{Plot character of points, defaults to 1. Usually set to a factor name or a column name in the input data with custom symbols}

\item{model.line}{Include line indicating estimates, default = TRUE. Can be an effect name.}

\item{ylim}{Y axis limits (\code{numeric}, but defaults to NULL)}

\item{ylab}{Y label (\code{character})}

\item{xlab}{X label (\code{character})}

\item{main}{Main title, defaults to \code{y} with description from \code{what}.}

\item{object}{GEM object.}

\item{variable}{Numeric for selecting a variable for extraction.}
}
\value{
\code{GEM} returns an object of class \code{GEM} containing effects, ER values (effect + residuals),
fitted values, residuals, features, coefficients, dummy design, symbolic design, dimensions,
highest level interaction and feature names.
}
\description{
General Effect Modelling
}
\examples{
## Multiple Sclerosis
data(MS, package = "gemR")
gem <- GEM(proteins ~ MS * group, data = MS)
print(gem)
plot(gem)                                       # Raw data, first feature
plot(gem,2)                                     # Raw data, numbered feature
plot(gem,'Q76L83', col='MS', pch='group')     # Selected colour and plot character
plot(gem,'Q76L83', what='effect MS',
     model.line='effect group')              # Comparison of factors (points and lines)
\donttest{
  # Example compound plot
  old.par <- par(c("mfrow", "mar"))
  # on.exit(par(old.par))
  par(mfrow = c(3,3), mar = c(2,4,4,1))
  plot(gem,'Q76L83')                         # Raw data, named feature
  plot(gem,'Q76L83', what='fits')            # Fitted values
  plot(gem,'Q76L83', what='residuals')       # Residuals
  plot(gem,'Q76L83', what='effect MS')       # Effect levels
  plot(gem,'Q76L83', what='effect group')    # ----||----
  plot(gem,'Q76L83', what='effect MS:group') # ----||----
  plot(gem,'Q76L83', what='MS')              # ER values
  plot(gem,'Q76L83', what='group')           # --------||---------
  plot(gem,'Q76L83', what='MS:group')        # --------||---------
  par(old.par)
}

# Complete overview of GEM
tab <- tableGEM(gem, 1)

# In general there can be more than two, effects, more than two levels, and continuous effects:
# MS$three <- factor(c(rep(1:3,33),1:2))
# gem3    <- GEM(proteins ~ MS * group + three, data = MS)


## Lactobacillus
data(Lactobacillus, package = "gemR")
gemLac <- GEM(proteome ~ strain * growthrate, data = Lactobacillus)
print(gemLac)
plot(gemLac)                            # Raw data, first feature
plot(gemLac,2)                          # Raw data, numbered feature
plot(gemLac,'P.LSA0316', col='strain',
    pch='growthrate')                  # Selected colour and plot character
plot(gemLac,'P.LSA0316', what='strain',
    model.line='growthrate')           # Selected model.line


## Diabetes
data(Diabetes, package = "gemR")
gemDia <- GEM(transcriptome ~ surgery * T2D, data = Diabetes)
print(gemDia)
plot(gemDia)                            # Raw data, first feature
plot(gemDia,2)                          # Raw data, numbered feature
plot(gemDia,'ILMN_1720829', col='surgery',
    pch='T2D')                         # Selected colour and plot character

}
\references{
* Mosleth et al. (2021) Cerebrospinal fluid proteome shows disrupted neuronal development in multiple sclerosis. Scientific Report, 11,4087. <doi:10.1038/s41598-021-82388-w>

* E.F. Mosleth et al. (2020). Comprehensive Chemometrics, 2nd edition; Brown, S., Tauler, R., & Walczak, B. (Eds.). Chapter 4.22. Analysis of Megavariate Data in Functional Omics. Elsevier. <doi:10.1016/B978-0-12-409547-2.14882-6>
}
\seealso{
Analyses using \code{GEM}: \code{\link{elastic}} and \code{\link{pls}}. Confidence interval plots: \code{\link{confints}}. Convenience knock-in and knock-out of effects: \code{\link{knock.in}}.
}
