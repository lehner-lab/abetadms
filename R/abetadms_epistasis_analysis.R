
#' abetadms_epistasis_analysis
#'
#' Epistasis analysis.
#'
#' @param fitness_path path to R data object with single and double mutant fitness values (required)
#' @param miscpath path to misc scripts and data directory (required)
#' @param DMS2structure_path Path to DMS2structure repository (required)
#' @param numCores Number of available CPU cores (default:1)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms_epistasis_analysis <- function(
  fitness_path,
  miscpath,
  DMS2structure_path,
  numCores = 1,
  execute = TRUE
  ){

  #Do nothing if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: abetadms_epistasis_analysis (this might take a while; you may want to adjust 'numCores' argument, DEFAULT=10)", "*******\n\n"))

  #Create output directory
  abetadms__create_dir(abetadms_dir = file.path(miscpath, "misc_epistasis_analysis"))

	### Run misc script on command-line
	###########################

  system(paste0(
  	file.path(miscpath, "scripts", "abetadms__epistasis_analysis.R"),
  	" -i ",
  	fitness_path,
  	" -o ",
  	file.path(miscpath, "misc_epistasis_analysis"),
  	" -d ",
  	DMS2structure_path,
  	" -c ",
  	numCores))

}








