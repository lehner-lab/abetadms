
#' abetadms_fitness_model_summary
#'
#' Dot plots showing explained variance of models to predict variant fitness.
#'
#' @param fitness_dt data.table with single and double mutant fitness values and their mutant effects from AA PCA (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms_fitness_model_summary <- function(
  fitness_dt,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: abetadms_fitness_model_summary", "*******\n\n"))

	#Create output directory
	abetadms__create_dir(abetadms_dir = outpath)

	### Setup
	###########################

	#Single AA mutants only (no STOPs)
	singles_dt <- copy(fitness_dt[Nmut_aa==1 & !STOP])
	#Double AA mutants only (no STOPs)
	doubles_dt <- copy(fitness_dt[Nmut_aa==2 & !STOP])
	doubles_dt[, fitness := fitness_cond]

	#Fitness hotspot positions
	mean_fitness <- singles_dt[STOP==F,mean(abs(fitness))]
	Pos_abs_hotspot <- singles_dt[Nmut_aa==1 & !STOP,.(hotspot = mean(abs(fitness))>mean_fitness),by=Pos_abs][hotspot==T,Pos_abs]
	singles_dt[, hotspot := as.numeric(Pos_abs %in% Pos_abs_hotspot)*2]
	doubles_dt[, hotspot := as.numeric(Pos_abs1 %in% Pos_abs_hotspot) + as.numeric(Pos_abs2 %in% Pos_abs_hotspot)]

	#Aggregation tool column names
	agg_tool_columns <- c(
		"DisEMBL_COILS",
		"DisEMBL_REM465",
		"DisEMBL_HOTLOOPS",
		"IUPred2A_short",
		"IUPred2A_long",
		"PolyPhen",
		"Tango",
		"Waltz",
		"Zyggregator")

	#AA properties column names
	aa_prop_columns <- names(singles_dt)[grep("^PC", names(singles_dt))]

	### Linear models to predict fitness - singles
	###########################

	feat_df_formodel <- cbind(
		data.frame(
			Position = singles_dt[,hotspot], 
			Fitness = singles_dt[,fitness]), 
		singles_dt[,.SD,.SDcols = aa_prop_columns], 
		singles_dt[,.SD,.SDcols = agg_tool_columns[agg_tool_columns %in% names(singles_dt)]])
	abetadms__fitness_lm_summary(
		input_df = feat_df_formodel, 
		output_file = file.path(outpath, 'fitness_lm_alldata_explainedvar_summary_singles.pdf'), 
		colour_scheme = colour_scheme,
		width = 9,
		height = 4)

	### Linear models to predict fitness - doubles
	###########################

	feat_df_formodel <- cbind(
		data.frame(
			Position = doubles_dt[,hotspot], 
			Fitness = doubles_dt[,fitness]), 
		doubles_dt[,.SD,.SDcols = aa_prop_columns], 
		doubles_dt[,.SD,.SDcols = agg_tool_columns[agg_tool_columns %in% names(doubles_dt)]])
	abetadms__fitness_lm_summary(
		input_df = feat_df_formodel, 
		output_file = file.path(outpath, 'fitness_lm_alldata_explainedvar_summary_doubles.pdf'), 
		colour_scheme = colour_scheme,
		width = 9,
		height = 4)

	### Linear models to predict fitness - singles and doubles
	###########################

	combined_dt <- rbind(singles_dt, doubles_dt)
	#Number of variants used in model
	print(paste0("All non-STOP AA variants (single, double AA substitutions): ", combined_dt[,.N], " (", paste(singles_dt[,.N], doubles_dt[,.N], sep = ", "), ")"))
	feat_df_formodel <- cbind(
		data.frame(
			Position = combined_dt[,hotspot], 
			Fitness = combined_dt[,fitness]), 
		combined_dt[,.SD,.SDcols = aa_prop_columns], 
		combined_dt[,.SD,.SDcols = agg_tool_columns[agg_tool_columns %in% names(combined_dt)]])
	lm_results <- abetadms__fitness_lm_summary(
		input_df = feat_df_formodel, 
		output_file = file.path(outpath, 'fitness_lm_alldata_explainedvar_summary_both.pdf'), 
		return_results = T,
		colour_scheme = colour_scheme,
		width = 9,
		height = 4)

	print(paste0("R-squared of Hydrophobicity in hotspot: ", round(lm_results[["PC1 (Hydrophobicity)_inside"]], 2)))
	print(paste0("Number of variants in hotspot: ", combined_dt[hotspot==max(hotspot),.N]))
	print(paste0("R-squared of Hydrophobicity inside/outside hotspot with hotspot variable: ", round(lm_results[["PC1 (Hydrophobicity)_pos"]], 2)))

	### Linear models to predict residual fitness - singles and doubles
	###########################

	combined_dt <- rbind(singles_dt, doubles_dt)
	combined_dt[, fitness_resid := lm(fitness ~ .*., data = .SD)[['residuals']], .SDcols = c('fitness', 'PC1 (Hydrophobicity)', 'hotspot')]
	#Number of variants used in model
	print(paste0("All non-STOP AA variants (single, double AA substitutions): ", combined_dt[,.N], " (", paste(singles_dt[,.N], doubles_dt[,.N], sep = ", "), ")"))
	feat_df_formodel <- cbind(
		data.frame(
			Position = combined_dt[,hotspot], 
			Fitness = combined_dt[,fitness_resid]), 
		combined_dt[,.SD,.SDcols = aa_prop_columns], 
		combined_dt[,.SD,.SDcols = agg_tool_columns[agg_tool_columns %in% names(combined_dt)]])
	lm_results <- abetadms__fitness_lm_summary(
		input_df = feat_df_formodel, 
		output_file = file.path(outpath, 'fitness_resid_lm_alldata_explainedvar_summary_both.pdf'), 
		return_results = T,
		colour_scheme = colour_scheme,
		width = 9,
		height = 4)

	print(paste0("R-squared of Zyggregator in hotspot (after controlling for Hydrophobicity and position): ", round(lm_results[["Zyggregator_inside"]], 2)))

}

