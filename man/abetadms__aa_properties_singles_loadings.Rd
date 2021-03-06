% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abetadms__aa_properties_singles_loadings.R
\name{abetadms__aa_properties_singles_loadings}
\alias{abetadms__aa_properties_singles_loadings}
\title{abetadms__aa_properties_singles_loadings}
\usage{
abetadms__aa_properties_singles_loadings(input_dt, aa_properties_file,
  selected_identifiers = NULL)
}
\arguments{
\item{input_dt}{data.table with double mutant codes (required)}

\item{aa_properties_file}{path to amino acid properties file (required)}

\item{selected_identifiers}{path to file with selected subset of identifiers}
}
\value{
input_dt with additional columns for loadings = delta(AA property scale) for each mutant
}
\description{
Add AA properties loadings for each mutant in input data.table
}
