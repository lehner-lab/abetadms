% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abetadms_preprocess_fitness.R
\name{abetadms_preprocess_fitness}
\alias{abetadms_preprocess_fitness}
\title{abetadms_preprocess_fitness}
\usage{
abetadms_preprocess_fitness(fitness_path, outpath, all_reps,
  bayesian_double_fitness = F, min_mean_input_read_count = 10,
  min_input_read_count_doubles = 50, lam_d = 0.025, numCores = 10,
  execute = TRUE)
}
\arguments{
\item{fitness_path}{path to DiMSum output (required)}

\item{outpath}{output path for plots and saved objects (required)}

\item{all_reps}{list of replicates to retain (required)}

\item{bayesian_double_fitness}{Estimate fitness of double mutants using bayesian framework (default:F)}

\item{min_mean_input_read_count}{minimum mean input read count for high confidence variants (default:10)}

\item{min_input_read_count_doubles}{minimum input read count for doubles used to derive prior for Bayesian doubles correction (default:50)}

\item{lam_d}{Poisson distribution for score likelihood (default:0.025)}

\item{numCores}{Number of available CPU cores (default:1)}

\item{execute}{whether or not to execute the analysis (default: TRUE)}
}
\value{
Nothing
}
\description{
Estimate fitness of doubles mutants using bayesian framework from DiMSum fitness output.
}
