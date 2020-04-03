
#' abetadms_AB42_structure_propensities
#'
#' AB42 structure predictions.
#'
#' @param result_dir directory with AB42 structure match results (required)
#' @param outpath output path for plots and saved objects (required)
#' @param miscpath path to misc scripts and data directory (required)
#' @param DMS2structure_path Path to DMS2structure repository (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#' @param rerun_structure re-run structure analysis? (default:F)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms_AB42_structure_propensities <- function(
  result_dir,
  outpath,
  miscpath,
  DMS2structure_path,
  colour_scheme,
  execute = TRUE,
  rerun_structure = F
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: abetadms_AB42_structure_propensities", "*******\n\n"))

	#Create output directory
	abetadms__create_dir(abetadms_dir = outpath)

	### Re-run structure analysis
	###########################

	if(rerun_structure){

	  #Display status
	  message(paste("\n\n*******", "re-running structure propensity calculations (this might take a while)", "*******\n\n"))

	  #Create output directory
	  abetadms__create_dir(abetadms_dir = file.path(miscpath, "misc_AB42_structures"))

		#Run misc script on command-line
	  system(paste0(
	  	file.path(miscpath, "scripts", "abetadms__AB42_structures.R"),
	  	" -o ",
	  	file.path(miscpath, "misc_AB42_structures"),
	  	" -d ",
	  	DMS2structure_path,
	  	" -e ",
	  	file.path(miscpath, "misc_epistasis_analysis"),
	  	" -p ",
	  	file.path(miscpath, "pdb")))
	}

	### Setup
	###########################

	#List of structures
	pdb_list <- list(
	  "1iyt_monomer" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 1),
	  "2beg_filament" = list("aa_seq" = "LVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 17),
	  "2beg_monomer" = list("aa_seq" = "LVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 17),
	  "2lmn_fibril" = list("aa_seq" = "GYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 9),
	  "2lmn_filament" = list("aa_seq" = "GYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 9),
	  "2lmn_monomer" = list("aa_seq" = "GYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 9),
	  "2lmp_fibril" = list("aa_seq" = "GYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 9),
	  "2lmp_filament" = list("aa_seq" = "GYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 9),
	  "2lmp_monomer" = list("aa_seq" = "GYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 9),
	  "2lmq_fibril" = list("aa_seq" = "GYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 9),
	  "2lmq_filament" = list("aa_seq" = "GYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 9),
	  "2lmq_monomer" = list("aa_seq" = "GYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 9),
	  "2lnq_filament" = list("aa_seq" = "QKLVFFAENVGSNKGAIIGLMVGGVV", "idx_start" = 15),
	  "2lnq_monomer" = list("aa_seq" = "QKLVFFAENVGSNKGAIIGLMVGGVV", "idx_start" = 15),
	  "2m4j_fibril" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 1),
	  "2m4j_filament" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 1),
	  "2m4j_monomer" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 1),
	  "2mpz_fibril" = list("aa_seq" = "QKLVFFAENVGSNKGAIIGLMVGGVV", "idx_start" = 15),
	  "2mpz_filament" = list("aa_seq" = "QKLVFFAENVGSNKGAIIGLMVGGVV", "idx_start" = 15),
	  "2mpz_monomer" = list("aa_seq" = "QKLVFFAENVGSNKGAIIGLMVGGVV", "idx_start" = 15),
	  "2mvx_fibril" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFADVGSNKGAIIGLMVGGVV", "idx_start" = 1),
	  "2mvx_filament" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFADVGSNKGAIIGLMVGGVV", "idx_start" = 1),
	  "2mvx_monomer" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFADVGSNKGAIIGLMVGGVV", "idx_start" = 1),
	  "2mxu_filament" = list("aa_seq" = "EVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 11),
	  "2mxu_monomer" = list("aa_seq" = "EVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 11),
	  "2nao_fibril" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 1),
	  "2nao_filament" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 1),
	  "2nao_monomer" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 1),
	  "2onv_monomer" = list("aa_seq" = "GGVVIA", "idx_start" = 37),
	  "2y3j_fibril" = list("aa_seq" = "AIIGLM", "idx_start" = 30),
	  "2y3j_filament" = list("aa_seq" = "AIIGLM", "idx_start" = 30),
	  "2y3j_monomer" = list("aa_seq" = "AIIGLM", "idx_start" = 30),
	  "5aef_fibril" = list("aa_seq" = "QKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 15),
	  "5aef_monomer" = list("aa_seq" = "QKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 15),
	  "5kk3_fibril" = list("aa_seq" = "EVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 11),
	  "5kk3_filament" = list("aa_seq" = "EVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 11),
	  "5kk3_monomer" = list("aa_seq" = "EVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 11),
	  "5oqv_fibril" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 1),
	  "5oqv_filament" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 1),
	  "5oqv_monomer" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA", "idx_start" = 1),
	  "6shs_fibril" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 1),
	  "6shs_filament" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 1),
	  "6shs_monomer" = list("aa_seq" = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV", "idx_start" = 1))

	kernel_scores_kernal_width <- list()
	for(i in names(pdb_list)){
	  kernel_scores_kernal_width[[i]] <- fread(file.path(result_dir, "pdb", i, "processed_data", "rand_strategy_kernal_width_kernel_structure_propensity.txt"))
	}

	### Kernel smooth scores - randomisation strategy: kernal width
	###########################

	#Distinct segments
	pdb_sets <- list(
		"1iyt" = c("1iyt_monomer"),
		"2beg" = c("2beg_filament", "2beg_monomer"),
		"2lmn" = c("2lmn_fibril", "2lmn_filament", "2lmn_monomer"),
		"2lmp" = c("2lmp_fibril", "2lmp_filament", "2lmp_monomer"),
		"2lmq" = c("2lmq_fibril", "2lmq_filament", "2lmq_monomer"),
		"2lnq" = c("2lnq_filament", "2lnq_monomer"),
		"2m4j" = c("2m4j_fibril", "2m4j_filament", "2m4j_monomer"),
		"2mpz" = c("2mpz_fibril", "2mpz_filament", "2mpz_monomer"),
		"2mvx" = c("2mvx_fibril", "2mvx_filament", "2mvx_monomer"),
		"2mxu" = c("2mxu_filament", "2mxu_monomer"),
		"2nao" = c("2nao_fibril", "2nao_filament", "2nao_monomer"),
		"2onv" = c("2onv_monomer"),
		"2y3j" = c("2y3j_fibril", "2y3j_filament", "2y3j_monomer"),
		"5aef" = c("5aef_fibril", "5aef_monomer"),
		"5kk3" = c("5kk3_fibril", "5kk3_filament", "5kk3_monomer"),
		"5oqv" = c("5oqv_fibril", "5oqv_filament", "5oqv_monomer"),
		"6shs" = c("6shs_fibril", "6shs_filament", "6shs_monomer"))
	# #LARKS mutants
	# pdb_LARKS <- c("5whn", "5whp", "5wkb")
	# pdb_sets <- list("segs" = pdb_segs)
	# pdb_sets <- as.list(pdb_segs)
	# names(pdb_sets) <- unlist(pdb_sets)
	for(pdb_set_name in names(pdb_sets)){
	  pdb_names <- pdb_sets[[pdb_set_name]]
	  for(score_name in c("association_score_norm", "epistasis_score_norm")){
	    abetadms__AB42_structure_propensity_plot(
	    	kernel_scores=kernel_scores_kernal_width, 
	    	outpath=file.path(outpath, paste0("kernel_propensity_kernal_width_", pdb_set_name, "_", score_name, ".pdf")),
	    	colour_scheme=colour_scheme,
	    	pdb_names=pdb_names, 
	    	score_type=score_name, 
	    	pdb_list=pdb_list, 
	    	position_offset=1)
	   }
	}

}








