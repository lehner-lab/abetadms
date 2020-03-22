
#' abetadms_quality_control
#'
#' Quality control plots including scatterplots comparing fitness estimates between replicates.
#'
#' @param fitness_list named list of folder paths with fitness estimates (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms_quality_control <- function(
  fitness_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: abetadms_quality_control", "*******\n\n"))

	#Create output directory
	abetadms__create_dir(abetadms_dir = outpath)
	
	### Load data
	###########################

	#DMS variant data (for individual replicates)
	dms_dt1 <- abetadms__load_fitness(fitness_list[[1]], object=T, read_filter=F)

	#Rename column names consistently
	names(dms_dt1)[grep('fitness[1-4]$', names(dms_dt1))] <- paste0("fitness", 1:length(grep('fitness[1-4]$', names(dms_dt1))))
	names(dms_dt1)[grep('fitness[1-4]_uncorr$', names(dms_dt1))] <- paste0("fitness", 1:length(grep('fitness[1-4]_uncorr$', names(dms_dt1))), "_uncorr")
	names(dms_dt1)[grep('fitness[1-4]_cond$', names(dms_dt1))] <- paste0("fitness", 1:length(grep('fitness[1-4]_cond$', names(dms_dt1))), "_cond")
	names(dms_dt1)[grep('sigma[1-4]$', names(dms_dt1))] <- paste0("sigma", 1:length(grep('sigma[1-4]$', names(dms_dt1))))
	names(dms_dt1)[grep('sigma[1-4]_uncorr$', names(dms_dt1))] <- paste0("sigma", 1:length(grep('sigma[1-4]_uncorr$', names(dms_dt1))), "_uncorr")
	names(dms_dt1)[grep('sigma[1-4]_cond$', names(dms_dt1))] <- paste0("sigma", 1:length(grep('sigma[1-4]_cond$', names(dms_dt1))), "_cond")

	#Combine regions for supplemental table
	dms_dt1[, region := names(fitness_list)[1]]
	dms_dt <- dms_dt1[,.SD,,.SDcols = c("Pos", "WT_AA", "Mut", "Nmut_nt", "Nmut_codons", "sigma", "fitness", "fitness_cond",
			"Pos1", "Pos2", "WT_AA1", "WT_AA2", "Mut1", "Mut2", "Nmut_aa", "region", "STOP", "mean_count", "is.reads0",
			names(dms_dt1)[grep('fitness[1-4]$', names(dms_dt1))], names(dms_dt1)[grep('fitness[1-4]_uncorr$', names(dms_dt1))], names(dms_dt1)[grep('fitness[1-4]_cond$', names(dms_dt1))],
			names(dms_dt1)[grep('sigma[1-4]$', names(dms_dt1))], names(dms_dt1)[grep('sigma[1-4]_uncorr$', names(dms_dt1))], names(dms_dt1)[grep('sigma[1-4]_cond$', names(dms_dt1))])]
	#Remove missing fitness values
	dms_dt <- dms_dt[!is.na(fitness) | !is.na(fitness_cond),]
	#Absolute position (singles)
	dms_dt[, Pos_abs := Pos+as.numeric(region)-1]
	#Absolute position (doubles)
	dms_dt[, Pos_abs1 := Pos1+as.numeric(region)-1]
	dms_dt[, Pos_abs2 := Pos2+as.numeric(region)-1]
	#Mutation code (singles)
	dms_dt[, mut_code := paste0(WT_AA, Pos_abs, Mut)]
	#Mutation code (doubles)
	dms_dt[, mut_code1 := paste0(WT_AA1, Pos_abs1, Mut1)]
	dms_dt[, mut_code2 := paste0(WT_AA2, Pos_abs2, Mut2)]
	#Supplementary data file
	fwrite(dms_dt[is.reads0==T,], file = file.path(outpath, "supplementary_table_raw_fitness_estimates_all.txt"), sep = "\t")

	#Combine regions for QC
	dms_dt <- dms_dt1[,.SD,,.SDcols = c("Nmut_aa", "region", "fitness", "fitness_cond", "fitness_uncorr", "STOP", "mean_count", "is.reads0",
			names(dms_dt1)[grep('fitness[1-4]$', names(dms_dt1))], names(dms_dt1)[grep('fitness[1-4]_cond$', names(dms_dt1))])]

	### Mean counts vs uncorrected fitness
	###########################

	#Set uncorrected fitness to fitness for singles
	dms_dt[Nmut_aa==1, fitness_uncorr := fitness]
	#Set conditional fitness to fitness for singles
	dms_dt[Nmut_aa==1, fitness_cond := fitness]
	plot_df <- reshape2::melt(as.data.frame(dms_dt[,.(mean_count, fitness_uncorr, fitness_cond, region)]), id = c("mean_count", "region"))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(mean_count, value)) +
		ggplot2::scale_x_log10() +
    ggplot2::stat_binhex(bins=30) +
	  ggplot2::xlab("Mean input read count") +
	  ggplot2::ylab("Fitness") +
	  ggplot2::geom_hline(yintercept = 0, linetype=2) +
	  ggplot2::geom_vline(xintercept = 10, linetype=2, colour = colour_scheme[["shade 0"]][[1]]) +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(variable ~ region)
	suppressWarnings(suppressMessages(ggplot2::ggsave(file=file.path(outpath, 'fitness_vs_mean_input_read_count.pdf'), width=5, height=4, useDingbats=FALSE)))

	### Scatterplots comparing fitness estimates between replicates
	###########################

	#Singles, doubles and STOPs mutants with defined fitness
	dms_dt <- dms_dt[is.reads0==T & (Nmut_aa==1 | Nmut_aa==2) & (!is.na(fitness) | !is.na(fitness_cond))]
	# dms_dt <- dms_dt[(Nmut_aa==1 | Nmut_aa==2) & (!is.na(fitness) | !is.na(fitness_cond))]
	dms_dt[Nmut_aa==2, fitness1 := fitness1_cond]
	dms_dt[Nmut_aa==2, fitness2 := fitness2_cond]
	dms_dt[Nmut_aa==2, fitness3 := fitness3_cond]
	dms_dt[Nmut_aa==2, fitness4 := fitness4_cond]

	#Convert infinite fitness values to NA
	dms_dt[is.infinite(fitness1), fitness1 := NA]
	dms_dt[is.infinite(fitness2), fitness2 := NA]
	dms_dt[is.infinite(fitness3), fitness3 := NA]
	dms_dt[is.infinite(fitness4), fitness4 := NA]

	#Number of variants with defined fitness values in each replicate
	print(paste0("All AA variants (single, double AA substitutions): ", dms_dt[,.N], " (", paste0(unlist(dms_dt[,.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep1 AA variants (single, double AA substitutions): ", dms_dt[!is.na(fitness1),.N], " (",paste0(unlist(dms_dt[!is.na(fitness1),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep2 AA variants (single, double AA substitutions): ", dms_dt[!is.na(fitness2),.N], " (",paste0(unlist(dms_dt[!is.na(fitness2),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep3 AA variants (single, double AA substitutions): ", dms_dt[!is.na(fitness3),.N], " (",paste0(unlist(dms_dt[!is.na(fitness3),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	print(paste0("Rep4 AA variants (single, double AA substitutions): ", dms_dt[!is.na(fitness4),.N], " (",paste0(unlist(dms_dt[!is.na(fitness4),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	# dms_dt[!is.na(fitness4),.N]

	#Plot replicate fitness correlation
	#num_rep = 3
	num_rep = 4
	abetadms__ggpairs_binhex(
		input_dt = dms_dt[Nmut_aa %in% c(1,2),.SD,.SDcols = paste0('fitness', c(1:num_rep))], 
		# input_dt_upper = dms_dt[Nmut_aa==1,.SD,.SDcols = paste0('fitness', c(1:num_rep))],
		output_file = file.path(outpath, "replicate_scatter_binhex.pdf"),
		xlab = "Fitness",
		ylab = "Fitness",
		width = 5,
		height = 5,
		label_size = 2,
  	plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]]),
  	colour_limits = c(0, 100))

	#Plot replicate 1 and 2 correlation
	plot_dt <- copy(dms_dt)[Nmut_aa %in% c(1,2) & !is.na(fitness1) & !is.na(fitness2)]
  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  temp_cor <- plyr::ddply(plot_dt, c("region"), plyr::summarize, cor = round(cor(fitness1, fitness2, use = "pairwise.complete.obs"), 2), n = length(fitness1))
  d <- ggplot2::ggplot(plot_dt, ggplot2::aes(fitness1, fitness2)) +
  	ggplot2::stat_binhex(bins=50) +
  	ggplot2::geom_abline(linetype = 2) + 
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 0, y = 0, size = 2) +
    ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 100)) +
    ggplot2::ylab("Fitness (rep2)") +
    ggplot2::xlab("Fitness (rep1)") +
    ggplot2::theme_classic() + 
    ggplot2::facet_grid(region ~.)
  ggplot2::ggsave(file=file.path(outpath, "replicate1_replicate2_scatter.pdf"), width=5, height=5, useDingbats=FALSE)

	# ### Scale fitness by median standard deviation of AA mutants in each library
	# ###########################

	# abetadms__scale_fitness(
	#   fitness_path=fitness_list[[1]],
	#   outpath=file.path(outpath, names(fitness_list)[1]))

	### Load data
	###########################

	#DMS variant data (for individual replicates)
	dms_dt1 <- abetadms__load_fitness(fitness_list[[1]], object=T, read_filter=F)

	#Rename column names consistently
	names(dms_dt1)[grep('fitness[1-4]$', names(dms_dt1))] <- paste0("fitness", 1:length(grep('fitness[1-4]$', names(dms_dt1))))
	names(dms_dt1)[grep('fitness[1-4]_cond$', names(dms_dt1))] <- paste0("fitness", 1:length(grep('fitness[1-4]_cond$', names(dms_dt1))), "_cond")

	#Combine regions
	dms_dt1[, region := names(fitness_list)[1]]
	dms_dt <- dms_dt1[,.SD,,.SDcols = c("Nmut_aa", "region", "fitness", "fitness_cond", "fitness_uncorr", "STOP", "mean_count", "is.reads0",
			names(dms_dt1)[grep('fitness[1-4]$', names(dms_dt1))], names(dms_dt1)[grep('fitness[1-4]_cond$', names(dms_dt1))])]

	### Mean counts vs uncorrected fitness
	###########################

	#Set uncorrected fitness to fitness for singles
	dms_dt[Nmut_aa==1, fitness_uncorr := fitness]
	#Set conditional fitness to fitness for singles
	dms_dt[Nmut_aa==1, fitness_cond := fitness]
	plot_df <- reshape2::melt(as.data.frame(dms_dt[,.(mean_count, fitness_uncorr, fitness_cond, region)]), id = c("mean_count", "region"))
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(mean_count, value)) +
		ggplot2::scale_x_log10() +
    ggplot2::stat_binhex(bins=30) +
	  ggplot2::xlab("Mean input read count") +
	  ggplot2::ylab("Fitness") +
	  ggplot2::geom_hline(yintercept = 0, linetype=2) +
	  ggplot2::geom_vline(xintercept = 10, linetype=2, colour = colour_scheme[["shade 0"]][[1]]) +
	  ggplot2::theme_bw() +
	  ggplot2::facet_grid(variable ~ region)
	suppressWarnings(suppressMessages(ggplot2::ggsave(file=file.path(outpath, 'fitness_vs_mean_input_read_count_scale.pdf'), width=5, height=4, useDingbats=FALSE)))

	### Scatterplots comparing fitness estimates between replicates
	###########################

	#Singles, doubles and STOPs mutants with defined fitness
	dms_dt <- dms_dt[is.reads0==T & (Nmut_aa==1 | Nmut_aa==2) & (!is.na(fitness) | !is.na(fitness_cond))]
	# dms_dt <- dms_dt[(Nmut_aa==1 | Nmut_aa==2) & (!is.na(fitness) | !is.na(fitness_cond))]
	dms_dt[Nmut_aa==2, fitness1 := fitness1_cond]
	dms_dt[Nmut_aa==2, fitness2 := fitness2_cond]
	dms_dt[Nmut_aa==2, fitness3 := fitness3_cond]
	dms_dt[Nmut_aa==2, fitness4 := fitness4_cond]

	#Convert infinite fitness values to NA
	dms_dt[is.infinite(fitness1), fitness1 := NA]
	dms_dt[is.infinite(fitness2), fitness2 := NA]
	dms_dt[is.infinite(fitness3), fitness3 := NA]
	dms_dt[is.infinite(fitness4), fitness4 := NA]

	# #Number of variants with defined fitness values in each replicate
	# print(paste0("All AA variants (single, double AA substitutions): ", dms_dt[,.N], " (", paste0(unlist(dms_dt[,.N,Nmut_aa][,2]), collapse = ", "), ")"))
	# print(paste0("Rep1 AA variants (single, double AA substitutions): ", dms_dt[!is.na(fitness1),.N], " (",paste0(unlist(dms_dt[!is.na(fitness1),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	# print(paste0("Rep2 AA variants (single, double AA substitutions): ", dms_dt[!is.na(fitness2),.N], " (",paste0(unlist(dms_dt[!is.na(fitness2),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	# print(paste0("Rep3 AA variants (single, double AA substitutions): ", dms_dt[!is.na(fitness3),.N], " (",paste0(unlist(dms_dt[!is.na(fitness3),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	# print(paste0("Rep4 AA variants (single, double AA substitutions): ", dms_dt[!is.na(fitness4),.N], " (",paste0(unlist(dms_dt[!is.na(fitness4),.N,Nmut_aa][,2]), collapse = ", "), ")"))
	# dms_dt[!is.na(fitness4),.N]

	#Plot replicate fitness correlation
	# num_rep = 3
	num_rep = 4
	abetadms__ggpairs_binhex(
		input_dt = dms_dt[Nmut_aa %in% c(1,2),.SD,.SDcols = paste0('fitness', c(1:num_rep))], 
		# input_dt_upper = dms_dt[Nmut_aa==1,.SD,.SDcols = paste0('fitness', c(1:num_rep))],
		output_file = file.path(outpath, "replicate_scatter_binhex_scale.pdf"),
		xlab = "Fitness",
		ylab = "Fitness",
		width = 5,
		height = 5,
		label_size = 2,
  	plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]]),
  	colour_limits = c(0, 100))

	#Plot replicate 1 and 2 correlation
	plot_dt <- copy(dms_dt)[Nmut_aa %in% c(1,2) & !is.na(fitness1) & !is.na(fitness2)]
  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]], colour_scheme[["shade 0"]][[3]])
  temp_cor <- plyr::ddply(plot_dt, c("region"), plyr::summarize, cor = round(cor(fitness1, fitness2, use = "pairwise.complete.obs"), 2), n = length(fitness1))
  d <- ggplot2::ggplot(plot_dt, ggplot2::aes(fitness1, fitness2)) +
  	ggplot2::stat_binhex(bins=50) +
  	ggplot2::geom_abline(linetype = 2) + 
	  ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("R = ", cor, "\nn = ", n)), x = 0, y = 0, size = 2) +
    ggplot2::scale_fill_gradientn(colours=plot_colours,name = "Frequency", na.value=plot_colours[length(plot_colours)], limits=c(0, 100)) +
    ggplot2::ylab("Fitness (rep2)") +
    ggplot2::xlab("Fitness (rep1)") +
    ggplot2::theme_classic() + 
    ggplot2::facet_grid(region ~.)
  ggplot2::ggsave(file=file.path(outpath, "replicate1_replicate2_scatter_scale.pdf"), width=5, height=5, useDingbats=FALSE)

}

