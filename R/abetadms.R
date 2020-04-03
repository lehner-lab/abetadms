
#' abetadms
#'
#' Main analysis script.
#'
#' @param startStage Start at a specified analysis stage (default:1)
#' @param stopStage Stop at a specified analysis stage (default:0 i.e. no stop condition)
#' @param base_dir Base directory for all output file (default:NB private CRG server path; change accordingly)
#' @param DMS2structure_path Path to DMS2structure repository (default:NB private CRG server path; change accordingly)
#' @param bayesian_double_fitness Estimate fitness of double mutants using bayesian framework (default:F)
#' @param rerun_epistasis re-run epistasis analysis? (default:F)
#' @param rerun_structure re-run structure analysis? (default:F)
#' @param numCores Number of available CPU cores (default:10)
#'
#' @return Nothing
#' @export
abetadms <- function(
  startStage=1,
  stopStage=0,
  base_dir = "/nfs/users/project/prj004631/afaure/DMS/Results/abetadms_proj_bdf",
  DMS2structure_path = "/users/blehner/afaure/DMS/Code/DMS2structure",
  bayesian_double_fitness = F,
  rerun_epistasis = F,
  rerun_structure = F,
  numCores = 10
  ){

	colour_scheme <- list(
		"shade 0" = list(
			"#F4270C",
			"#F4AD0C",
			"#1B38A6",
			"#09B636"),
		"shade 1" = list(
			"#FFB0A5",
			"#FFE4A5",
			"#9DACE3",
			"#97E9AD"),
		"shade 2" = list(
			"#FF6A56",
			"#FFCB56",
			"#4C63B7",
			"#43C766"),
		"shade 3" = list(
			"#A31300",
			"#A37200",
			"#0C226F",
			"#007A20"),
		"shade 4" = list(
			"#410800",
			"#412D00",
			"#020B2C",
			"#00300D"))

  #First and last analysis stages
  first_stage <- startStage
  last_stage <- stopStage

	#DiMSum fitness
	stagenum <- 1
	abetadms_preprocess_fitness(
	  fitness_path = file.path(base_dir, "misc", "DiMSum_fitness", "AB42_COMP_ABCD_any50_fitness_replicates.RData"),
	  outpath = abetadms__format_dir(dir_suffix="_abetadms_preprocess_fitness", stagenum=stagenum, base_dir=base_dir),
	  all_reps = 1:4,
	  bayesian_double_fitness = bayesian_double_fitness,
		numCores = numCores,
	  execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Quality control plots
	stagenum <- 2
	abetadms_quality_control(
		fitness_list = list(
			"1" = file.path(base_dir, "001_abetadms_preprocess_fitness")),
		outpath = abetadms__format_dir(dir_suffix="_abetadms_quality_control", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Combine fitness estimates
	stagenum <- 3
	fitness_norm_dt <- abetadms_combine_fitness(
		fitness_list = list(
			"1" = file.path(base_dir, "001_abetadms_preprocess_fitness")),
		disease_mut_file = file.path(base_dir, "misc", "abeta_reported_mutations.txt"),
		outpath = abetadms__format_dir(dir_suffix="_abetadms_combine_fitness", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Calculate single and double mutant effects from AA PCA (using single mutants)
	stagenum <- 4
	fitness_aaprop_dt <- abetadms_aa_properties_mutant_effects(
		fitness_dt = fitness_norm_dt,
		outpath = abetadms__format_dir(dir_suffix="_abetadms_aa_properties_mutant_effects", stagenum=stagenum, base_dir=base_dir),
		aaprop_file = file.path(base_dir, "misc", "amino_acid_properties", "amino_acid_properties_annotated_supplementary.txt"),
		aaprop_file_selected = file.path(base_dir, "misc", "amino_acid_properties", "selected.amino_acid_properties.txt"),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Calculate single and double mutant effects from aggregation tool predictions (using single mutants)
	stagenum <- 5
	fitness_aggtools_dt <- abetadms_agg_tools_mutant_effects(
		fitness_dt = fitness_aaprop_dt,
		aggtool_results_file = file.path(base_dir, "misc", "aggregation_tools_results.RData"),
		# aggscale_file = file.path(base_dir, "misc", "amino_acid_properties", "aggregation_scales.txt"),
		outpath = abetadms__format_dir(dir_suffix="_abetadms_agg_tools_mutant_effects", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Single mutant heatmaps
	stagenum <- 6
	abetadms_single_mutant_heatmaps(
		fitness_dt = fitness_aggtools_dt,
		disease_mut_file = file.path(base_dir, "misc", "abeta_reported_mutations.txt"),
		outpath = abetadms__format_dir(dir_suffix="_abetadms_single_mutant_heatmaps", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Human disease mutations
	stagenum <- 7
	abetadms_human_disease_mutations(
		fitness_dt = fitness_aggtools_dt,
		missense_AF_file = file.path(base_dir, "misc", "app_gnomAD_r2.1.1_missense_allelefreqs.txt"),
		disease_mut_file = file.path(base_dir, "misc", "abeta_reported_mutations.txt"),
		outpath = abetadms__format_dir(dir_suffix="_abetadms_human_disease_mutations", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Dot plots showing explained variance of models to predict variant fitness.
	stagenum <- 8
	abetadms_fitness_model_summary(
		fitness_dt = fitness_aggtools_dt,
		outpath = abetadms__format_dir(dir_suffix="_abetadms_fitness_model_summary", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	# #Violin plots showing (residual) fitness versus number of introduced AAs of various properties
	# stagenum <- 9
	# abetadms_num_introduced_aa_violins(
	# 	fitness_dt = fitness_aggtools_dt,
	# 	outpath = abetadms__format_dir(dir_suffix="_abetadms_num_introduced_aa_violins", stagenum=stagenum, base_dir=base_dir),
	# 	colour_scheme = colour_scheme,
	# 	execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Helix propensity of WT and fitness hotspot line plot
	stagenum <- 9
	abetadms_wt_helix_propensity(
		fitness_dt = fitness_aggtools_dt,
		outpath = abetadms__format_dir(dir_suffix="_abetadms_wt_helix_propensity", stagenum=stagenum, base_dir=base_dir),
		aaprop_file = file.path(base_dir, "misc", "amino_acid_properties", "amino_acid_properties_annotated_supplementary.txt"),
		aaprop_file_selected = file.path(base_dir, "misc", "amino_acid_properties", "selected.amino_acid_properties.txt"),
		wtfasta_path = file.path(base_dir, "misc", "AB42.fa"),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))

	#Epistasis analysis
	stagenum <- 10
	abetadms_epistasis_analysis(
		fitness_path = file.path(base_dir, "003_abetadms_combine_fitness", "combine_fitness_values.RData"),
		miscpath = file.path(base_dir, "misc"),
		DMS2structure_path = DMS2structure_path,
		numCores = numCores,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum) & rerun_epistasis))

	#Secondary structure predictions
	stagenum <- 11
	abetadms_secondary_structure_predictions(
		fitness_dt = fitness_aggtools_dt,
		outpath = abetadms__format_dir(dir_suffix="_abetadms_secondary_structure_predictions", stagenum=stagenum, base_dir=base_dir),
		miscpath = file.path(base_dir, "misc"),
		DMS2structure_path = DMS2structure_path,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)),
		rerun_structure = rerun_structure)

	#Guenther structure propensities
	stagenum <- 12
	abetadms_AB42_structure_propensities(
		result_dir = file.path(base_dir, "misc", "misc_AB42_structures"),
		outpath = abetadms__format_dir(dir_suffix="_abetadms_AB42_structure_propensities", stagenum=stagenum, base_dir=base_dir),
		miscpath = file.path(base_dir, "misc"),
		DMS2structure_path = DMS2structure_path,
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)),
		rerun_structure = rerun_structure)

	#PWI heatmaps
	stagenum <- 13
	abetadms_PWI_heatmaps(
		PWI_dir_list = list(
			"1" = file.path(base_dir, "misc", "misc_epistasis_analysis/doubles_cond_1/processed_data/")),
		outpath = abetadms__format_dir(dir_suffix="_abetadms_PWI_heatmaps", stagenum=stagenum, base_dir=base_dir),
		colour_scheme = colour_scheme,
		execute = (first_stage <= stagenum & (last_stage == 0 | last_stage >= stagenum)))
}
