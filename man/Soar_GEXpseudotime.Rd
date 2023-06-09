% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Soar_GEXpseudotime.R
\name{Soar_GEXpseudotime}
\alias{Soar_GEXpseudotime}
\title{View Expression of given genes over pseudotime}
\usage{
Soar_GEXpseudotime(egret, pick.features, smooth.span = 0.75)
}
\arguments{
\item{egret}{An Egret object}

\item{pick.features}{A character vector of feature(s) contained in the gene expression data.}

\item{smooth.span}{Value controlling the smoothness of the loess span for curve fitting.}
}
\value{
a plot
}
\description{
View expression of given genes over pseudotime. Useful for manually picked genes (e.g. known markers or functional complexes).
potential visualization for trajectory-stable reference and trajectory-unstable noise genes as well.
}
