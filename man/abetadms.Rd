% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abetadms.R
\name{abetadms}
\alias{abetadms}
\title{abetadms}
\usage{
abetadms(startStage = 1, stopStage = 0,
  base_dir = "/nfs/users/project/prj004631/afaure/DMS/Results/abetadms_proj_bdf",
  DMS2structure_path = "/users/blehner/afaure/DMS/Code/DMS2structure",
  bayesian_double_fitness = F, rerun_epistasis = F,
  rerun_structure = F, numCores = 10)
}
\arguments{
\item{startStage}{Start at a specified analysis stage (default:1)}

\item{stopStage}{Stop at a specified analysis stage (default:0 i.e. no stop condition)}

\item{base_dir}{Base directory for all output file (default:NB private CRG server path; change accordingly)}

\item{DMS2structure_path}{Path to DMS2structure repository (default:NB private CRG server path; change accordingly)}

\item{bayesian_double_fitness}{Estimate fitness of double mutants using bayesian framework (default:F)}

\item{rerun_epistasis}{re-run epistasis analysis? (default:F)}

\item{rerun_structure}{re-run structure analysis? (default:F)}

\item{numCores}{Number of available CPU cores (default:10)}
}
\value{
Nothing
}
\description{
Main analysis script.
}
