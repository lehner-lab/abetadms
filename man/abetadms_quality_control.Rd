% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abetadms_quality_control.R
\name{abetadms_quality_control}
\alias{abetadms_quality_control}
\title{abetadms_quality_control}
\usage{
abetadms_quality_control(fitness_list, outpath, colour_scheme,
  execute = TRUE)
}
\arguments{
\item{fitness_list}{named list of folder paths with fitness estimates (required)}

\item{outpath}{output path for plots and saved objects (required)}

\item{colour_scheme}{colour scheme file (required)}

\item{execute}{whether or not to execute the analysis (default: TRUE)}
}
\value{
Nothing
}
\description{
Quality control plots including scatterplots comparing fitness estimates between replicates.
}
