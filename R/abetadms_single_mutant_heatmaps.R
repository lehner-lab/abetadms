
#' abetadms_single_mutant_heatmaps
#'
#' Combine fitness estimates from two different DMS experiments (regions).
#'
#' @param fitness_dt data.table with single mutant fitness values (required)
#' @param outpath output path for plots and saved objects (required)
#' @param disease_mut_file table of human disease mutations and classifications (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms_single_mutant_heatmaps <- function(
  fitness_dt,
  disease_mut_file,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: abetadms_single_mutant_heatmaps", "*******\n\n"))

	#Create output directory
	abetadms__create_dir(abetadms_dir = outpath)

	#Subset to single AA mutants only
	dms_dt <- copy(fitness_dt[Nmut_aa==1])

	#Fitness hotspot positions
	mean_fitness <- dms_dt[STOP==F,mean(abs(fitness))]
	Pos_abs_hotspot <- dms_dt[STOP==F,.(hotspot = mean(abs(fitness))>mean_fitness),by=Pos_abs][hotspot==T,Pos_abs]
	dms_dt[, hotspot := as.numeric(Pos_abs %in% Pos_abs_hotspot)]

	#Disease mutations
	dis_mut <- read.table(disease_mut_file, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
	fAD_muts <- rownames(dis_mut)

	### Heatmap of transformed data
	###########################

	#Fitness
	abetadms__single_mutant_heatmap(
		input_df = as.data.frame(dms_dt), 
		variable_name="fitness", 
		output_file = file.path(outpath, 'heatmap_fitness_center_scale.pdf'),
		x_annotation="sequence", 
		xaxis_angle=0, 
		na_colour="white", 
		na_text="-", 
		colour_low=colour_scheme[["shade 0"]][[3]], 
		colour_high=colour_scheme[["shade 0"]][[1]])

	#Fitness - cluster positions
	abetadms__single_mutant_heatmap(
		input_df = as.data.frame(dms_dt), 
		variable_name="fitness", 
		output_file = file.path(outpath, 'heatmap_fitness_center_scale_colcluster.pdf'),
		x_annotation="both", 
		xaxis_angle=90, 
		na_colour="white", 
		na_text="-", 
		colour_low=colour_scheme[["shade 0"]][[3]], 
		colour_high=colour_scheme[["shade 0"]][[1]],
		cluster="column",
		width=24,
		height=8)

	#Fitness (ALS mutations)
	abetadms__single_mutant_heatmap(
		input_df = as.data.frame(dms_dt), 
		variable_name="fitness", 
		mut_dict=list("F" = fAD_muts),
		output_file = file.path(outpath, 'heatmap_fitness_center_scale_hmuts.pdf'),
		x_annotation="sequence", 
		xaxis_angle=0, 
		na_colour="white", 
		na_text="-", 
		colour_low=colour_scheme[["shade 0"]][[3]], 
		colour_high=colour_scheme[["shade 0"]][[1]])

	#Residual fitness (after controlling for hydrophobicity and position)
	dms_dt[STOP==F, fitness_resid := lm(fitness ~ .*., data = .SD)[['residuals']], .SDcols = c('fitness', 'PC1 (Hydrophobicity)', 'hotspot')]
	# dms_dt[STOP==F, fitness_resid := loess(fitness ~ .SD[[2]] + .SD[[3]] + .SD[[2]]*.SD[[3]], data = .SD)[['residuals']], .SDcols = c('fitness', 'PC1 (Hydrophobicity)', 'hotspot')]
	dms_dt[, fitness_resid_sig := p.adjust(2*pnorm(-abs(fitness_resid/sigma)), method = "bonferroni")<0.01]
	abetadms__single_mutant_heatmap(
		input_df = as.data.frame(dms_dt), 
		variable_name="fitness_resid", 
		mut_dict=list("*" = dms_dt[fitness_resid_sig==T,mut_code]),
		output_file = file.path(outpath, 'heatmap_resid_fitness_center_scale.pdf'),
		x_annotation="sequence", 
		xaxis_angle=0, 
		na_colour="white", 
		na_text="-", 
		colour_low=colour_scheme[["shade 0"]][[3]], 
		colour_high=colour_scheme[["shade 0"]][[1]])

}

