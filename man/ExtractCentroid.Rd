% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Centroids.R
\name{ExtractCentroid}
\alias{ExtractCentroid}
\title{Visualisation of fuzzy centroids:}
\usage{
ExtractCentroid(
  ecotonefinder,
  method = c("fanny", "cmeans", "vegclust"),
  normalized = c("species", "cluster", "none"),
  position = "dodge",
  threshold = 0,
  plot = TRUE,
  col = NULL,
  return.plot = TRUE,
  labels = ggplot2::waiver(),
  main = "Community composition",
  xlab = "species",
  ylab = "Centroid contribution",
  cex.x = 12
)
}
\arguments{
\item{ecotonefinder}{A list containing elements named in the same way
than EcotoneFinder function outcomes. Must contain “cmeans”, “fanny”
or “vegclust” results.}

\item{method}{The fuzzy clustering results from which the centroids
will be extracted.}

\item{normalized}{Method to normalise the centroid values, either
by “species” or “cluster”. If “none”, the centroids are plotted
without transformation. See details.}

\item{position}{Set the positions of the bars for the barchart.
This is passed down to the geom_bar function of ggplot. Default
is set to “dodge”.}

\item{threshold}{Threshold for centroid contribution value under
which the species will not be plotted. Can be used to simplify
plots containing many species. See Details.}

\item{plot}{Logical. Should the plot be displayed. If FALSE,
the centroids matrix is returned without plotting.}

\item{col}{Colour vector for the plot. Should be of the same
length that the number of fuzzy clusters.}

\item{return.plot}{Logical. Should the GGplot object be
stored internally (e.g. for multi-ploting). Default is
TRUE.}

\item{labels}{Character vectors of labels for the legend. Must
be of the same length that the number of fuzzy clusters.}

\item{main}{Main title for the plot. See plot.}

\item{xlab}{A title for the x-axis. See plot.}

\item{ylab}{A title for the y-axis. See plot.}

\item{cex.x}{cex for the x-axis labels.}
}
\value{
A matrix containing the cluster centroids.
}
\description{
Visualisation of fuzzy centroids:
}
\details{
This function extracts and plots the fuzzy centroids species
 contributions, according to user-defined normalisation steps and
 threshold value. The contributions of the different species in the
 fuzzy centroids may be used as a proxy for community compositions.
 The cmeans function (cmeans package) and vegclust function (vegclust
 package) internally compute the centroid compositions and their
 outputs are directly used by the ExtractCentroid function.
 The fanny function (cluster package), however, does not provide
 internal centroids calculation. They are computed here as:
 \deqn{Centroid[cluster j] = \sum[ij] (Membership[ij] x Observation[ij]) / \sum[j] Membership[j]
 }
 Where the centroid of a cluster is the mean of all observations,
 weighted by their degree of belonging to the cluster.
 The obtained species contributions to the centroids of the fuzzy
 clusters can then be plotted as they are, if normalised = “none”.
 To obtain more intuitive units for the interpretation of the
 species contributions, two normalisation methods are proposed.
 If normalised = “cluster”, the species contributions are given
 in percent per clusters (i.e. the sum of all species contributions
 in each cluster centroid equals 100). If normalised = “species”,
 each species has its contributions summed to 100 (i.e. each species
 is in percent per cluster).
 For normalised = “none” and normalised = “cluster”, a threshold
 value can be specified. Species that do not score above this
 threshold will not be displayed on the resulting plot. This can
 be used to simplify the outputs, for dataset containing large
 number of species.
}
\examples{
##### Artificial dataset & analyses:
 SyntheticTrial <- SyntheticData(SpeciesNum = 20, CommunityNum = 4,
                                 SpCo = NULL ,Length = 500,
                                 Parameters=list(a = rep(60, 4),
                                                 b = c(0,150,350,500),
                                                 c = rep(0.015,4)),
                                 dev.c = 0.007,
                                 pal = c("#008585", "#B8CDAE", "#E6C186", "#C7522B"))

 EcoFinder <- EcotoneFinder(SyntheticTrial[,-1],
                            dist = SyntheticTrial$Distance,
                            method = "all", groups=4,
                            standardize = "hellinger",
                            diversity="all")

 ##### Centroid plot without normalisation:
 Centroid <- ExtractCentroid(EcoFinder, method = "fanny",
                             normalized = "none",  threshold = 0,
                             plot = TRUE, position = "dodge",
                             col = colorspace::heat_hcl(4))

 ##### Centroid plot normalised by clusters:
 Centroid <- ExtractCentroid(EcoFinder, method = "fanny",
                             normalized = "cluster",  threshold = 0,
                             plot = TRUE, position = "dodge",
                             col = colorspace::heat_hcl(4))

}
