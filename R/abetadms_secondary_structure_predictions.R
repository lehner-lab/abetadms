
#' abetadms_secondary_structure_predictions
#'
#' Secondary structure predictions.
#'
#' @param fitness_dt data.table with single and double mutant fitness values (required)
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
abetadms_secondary_structure_predictions <- function(
  fitness_dt,
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
  message(paste("\n\n*******", "running stage: abetadms_secondary_structure_predictions", "*******\n\n"))

	#Create output directory
	abetadms__create_dir(abetadms_dir = outpath)

	### Re-run structure analysis
	###########################

	if(rerun_structure){

	  #Display status
	  message(paste("\n\n*******", "re-running structure predictions (this might take a while)", "*******\n\n"))

	  #Create output directory
	  abetadms__create_dir(abetadms_dir = file.path(miscpath, "misc_secondary_structure"))

		#Run misc script on command-line
	  system(paste0(
	  	file.path(miscpath, "scripts", "abetadms__secondary_structure.R"),
	  	" -o ",
	  	file.path(miscpath, "misc_secondary_structure"),
	  	" -d ",
	  	DMS2structure_path,
	  	" -e ",
	  	file.path(miscpath, "misc_epistasis_analysis")))
	}

	### Setup
	###########################

	#Load secondary structure predictions
	load(file.path(miscpath, "misc_secondary_structure", "result_list.RData"))
	ss_dt <- rbind(
		result_list[["PWI_cond_1"]][["secondary_structure_score"]][["association_score"]][, region := "1"][, score_type := "association_score"],
		result_list[["PWI_cond_1"]][["secondary_structure_score"]][["posE_pcor"]][, region := "1"][, score_type := "posE_pcor"],
		result_list[["PWI_cond_1"]][["secondary_structure_score"]][["negE_pcor"]][, region := "1"][, score_type := "negE_pcor"])
	ss_dt[, Pos_abs := as.numeric(region)+Pos1-1]

	#Single AA mutants only (no STOPs)
	singles_dt <- copy(fitness_dt[Nmut_aa==1 & !STOP])
	#Double AA mutants only (no STOPs)
	doubles_dt <- copy(fitness_dt[Nmut_aa==2 & !STOP])
	doubles_dt[, fitness := fitness_cond]

	#Fitness hotspot positions
	mean_fitness <- singles_dt[,mean(abs(fitness))]
	Pos_abs_hotspot <- singles_dt[Nmut_aa==1 & !STOP,.(hotspot = mean(abs(fitness))>mean_fitness),by=Pos_abs][hotspot==T,Pos_abs]

	#Mutant positions
	Pos_abs_all <- unique(singles_dt[,Pos_abs])

	# #Structure coordinates
	# NTD_coords <- c(1, 78)
	# RRM1_coords <- c(101, 176)
	# RRM2_coords <- c(191, 262)
	# PRLD_coords <- c(274, 414)
	# LARKS_coords <- c(312, 317)
	# LARKS_coords_focus <- c(312, 314)
	# helix_coords <- c(321, 330)
	# helix_coords_focus <- c(326, 330)

	### Secondary structure significance - association score (both positive and negative epistasis)
	###########################

	#Line plot
	plot_dt <- copy(ss_dt)[score_type == "association_score",]
	plot_dt[, alpha_helix_pval := -log10(alpha_p_seed0)]
	plot_dt[, beta_strand_pval := -log10(beta_p_seed0)]
	plot_df <- reshape2::melt(as.data.frame(plot_dt[,.(Pos_abs, alpha_helix_pval, beta_strand_pval)]), id = 'Pos_abs')
	data_min <- min(plot_df[,'value']) + 0.2*min(plot_df[,'value'])
	data_max <- max(plot_df[,'value']) + 0.2*max(plot_df[,'value'])
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(Pos_abs, value, color=variable)) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_hotspot), xmax=max(Pos_abs_hotspot), ymin=data_min, ymax=data_max), fill = "lightgrey", linetype = 0, alpha = 0.1) +
	  # ggplot2::geom_rect(ggplot2::aes(xmin=min(LARKS_coords), xmax=max(LARKS_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	  # ggplot2::geom_rect(ggplot2::aes(xmin=min(helix_coords), xmax=max(helix_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	  ggplot2::geom_vline(xintercept = median(Pos_abs_all), linetype = 2) +
	  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
	  ggplot2::geom_line() +
	  ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(4, 2)]))
  }
	ggplot2::ggsave(file=file.path(outpath, 'tdp43_ss_pred_association_score_all.pdf'), width=8, height=3, useDingbats=FALSE)

	### Secondary structure significance - positive and negative association score
	###########################

	#Line plot
	plot_dt <- copy(ss_dt)[score_type %in% c("posE_pcor", "negE_pcor"),]
	plot_dt[, alpha_helix_pval := -log10(alpha_p_seed0)]
	plot_dt[, beta_strand_pval := -log10(beta_p_seed0)]
	plot_df <- reshape2::melt(as.data.frame(plot_dt[,.(Pos_abs, alpha_helix_pval, beta_strand_pval, score_type)]), id = c('Pos_abs', 'score_type'))
	data_min <- min(plot_df[,'value']) + 0.2*min(plot_df[,'value'])
	data_max <- max(plot_df[,'value']) + 0.2*max(plot_df[,'value'])
	d <- ggplot2::ggplot(plot_df, ggplot2::aes(Pos_abs, value, color=variable)) +
	  ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_hotspot), xmax=max(Pos_abs_hotspot), ymin=data_min, ymax=data_max), fill = "lightgrey", linetype = 0, alpha = 0.1) +
	  # ggplot2::geom_rect(ggplot2::aes(xmin=min(LARKS_coords), xmax=max(LARKS_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	  # ggplot2::geom_rect(ggplot2::aes(xmin=min(helix_coords), xmax=max(helix_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	  ggplot2::geom_vline(xintercept = median(Pos_abs_all), linetype = 2) +
	  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
	  ggplot2::geom_line(ggplot2::aes(linetype = score_type)) +
	  ggplot2::theme_classic()
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(4, 2)]))
  }
	ggplot2::ggsave(file=file.path(outpath, 'tdp43_ss_pred_association_score_pos_neg.pdf'), width=8, height=3, useDingbats=FALSE)

	### Correlation with epistasis in/outside regions - singles
	###########################

	epi_dt <- fread(file.path(miscpath, "misc_epistasis_analysis", "doubles_cond_1", "processed_data", "DMS_epistasis.txt"))[order(mut_code1, mut_code2)]
	doubles_dt <- fitness_dt[Nmut_aa==2,][order(mut_code1, mut_code2)]

	doubles_dt[, epistasis := epi_dt[,epistasis]]
	doubles_dt[, sigmaE := epi_dt[,sigmaE]]
	doubles_dt[, pos_epistasis_sig := epi_dt[,pos_epistasis_sig]]
	doubles_dt[, neg_epistasis_sig := epi_dt[,neg_epistasis_sig]]

  weighted_cor <- function(input_dt){
  	temp_dt <- input_dt[,.SD,,.SDcols = names(input_dt)[names(input_dt)!="wt"]]
  	temp_wt <- input_dt[,wt]
  	return(cov.wt(as.data.frame(temp_dt), wt = unlist(temp_wt), cor = T)$cor[1,2])
  }

	plot_double_epistasis_boxplot <- function(
		input_dt, 
		feature_name, 
		ranges_list,
		width = 5,
		height = 5,
		text_code = F,
		plot_path
		){
		input_dt[, metric := .SD,,.SDcols = feature_name]
		temp_list <- list()
		for(i in names(ranges_list)){
			temp_list[[i]] <- input_dt[Pos_abs1 %in% ranges_list[[i]] & Pos_abs2 %in% ranges_list[[i]],]
			temp_list[[i]][, region := paste0(i, ":", i)]
			temp_list[[i]] <- temp_list[[i]][!is.na(metric),.SD,,.SDcols = c("epistasis", "sigmaE", "metric", "region", "pos_epistasis_sig", "neg_epistasis_sig")]
		}
		plot_dt <- do.call("rbind", temp_list)
		plot_dt[, region := factor(region, levels = paste0(names(ranges_list), ":", names(ranges_list)))]
		plot_dt[, epistasis_sig := "Neither"]
		plot_dt[neg_epistasis_sig==T, epistasis_sig := "Negative"]
		plot_dt[pos_epistasis_sig==T, epistasis_sig := "Positive"]
	  plot_colours = c(colour_scheme[["shade 0"]][[1]], colour_scheme[["shade 0"]][[2]])
	  # plot_dt[, wt := 1/sigmaE^2]
	  # temp_cor <- plot_dt[,.(cor = weighted_cor(.SD), n = .N),region,.SDcols = c("epistasis", "metric", "wt")]
		d <- ggplot2::ggplot(plot_dt, ggplot2::aes(epistasis_sig, metric, colour = epistasis_sig)) +
	    ggplot2::geom_boxplot(notch = T) +
		  ggplot2::xlab("Significant epistasis") +
		  ggplot2::ylab(feature_name) +
		  ggplot2::theme_bw() +
		  ggplot2::geom_smooth(method = "lm", colour = "black") +
		  ggplot2::facet_wrap(region~.) +
		  # ggplot2::coord_cartesian(ylim=c(-0.4, 0.4)) +
		  # ggplot2::geom_text(data = temp_cor, ggplot2::aes(label = paste0("Rw = ", round(cor, 2), "\nn = ", n)), x = 0, y = 0, colour = "black", 
		  # 	vjust = "inward", hjust = "inward")
		# if(text_code){
		# 	d <- d + ggplot2::geom_text(size = 2, color = "black")
		# }
		ggplot2::ggsave(file=plot_path, width=width, height=height, useDingbats=FALSE)
	}

	feature_name <- "AGADIR"
	plot_double_epistasis_boxplot(
		input_dt = copy(doubles_dt)[!is.na(sigmaE)], 
		feature_name = feature_name, 
		ranges_list = list(
	  	Nterm = 1:26,
	  	Cterm = 27:42,
	  	Helix = 22:28),
		width = 10, height = 5,
		plot_path = file.path(outpath, paste0('epistasis_vs_', feature_name, '_boxplot_hotspot_doubles.pdf')))



	# #Dotplots
	# #Adjust p-values
	# plot_dt2 <- rbind(
	# 	copy(plot_dt)[, ss_type := "beta_strand"][, value := beta_strand_pval],
	# 	copy(plot_dt)[, ss_type := "alpha_helix"][, value := alpha_helix_pval])[, value_adjust := -log10(p.adjust(10^-value, method = "BH"))]
	# #Subset to LARKS and helix positions
	# plot_dt2[Pos_abs %in% seq(LARKS_coords[1], LARKS_coords[2]), feature := "LARKS"]
	# plot_dt2[feature=="LARKS" & ss_type=="alpha_helix", feature_cat := "LARKS_alpha"]
	# plot_dt2[feature=="LARKS" & ss_type=="beta_strand", feature_cat := "LARKS_beta"]
	# plot_dt2[Pos_abs %in% seq(helix_coords[1], helix_coords[2]), feature := "helix"]
	# plot_dt2[feature=="helix" & ss_type=="alpha_helix", feature_cat := "helix_alpha"]
	# plot_dt2[feature=="helix" & ss_type=="beta_strand", feature_cat := "helix_beta"]
	# plot_df <- as.data.frame(plot_dt2[feature=="LARKS" | feature=="helix",])
	# plot_df[,"feature_cat"] <- factor(plot_df[,"feature_cat"], levels = c("LARKS_beta", "LARKS_alpha", "helix_beta", "helix_alpha"))
	# plot_df[,"score_type"] <- factor(plot_df[,"score_type"], levels = c("posE_pcor", "negE_pcor"))
	# d <- ggplot2::ggplot(plot_df, ggplot2::aes(feature_cat, value_adjust, fill=score_type)) +
	#   ggplot2::geom_dotplot(binaxis = "y", stackdir = "center") +
	#   ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
	#   ggplot2::theme_classic() +
	#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
 #  if(!is.null(colour_scheme)){
 #    d <- d + ggplot2::scale_fill_manual(values = unlist(colour_scheme[["shade 0"]][c(3,1)]))
 #  }
	# suppressMessages(ggplot2::ggsave(file=file.path(outpath, 'tdp43_ss_pred_association_score_pos_neg_features_dotplot.pdf'), width=4, height=4, useDingbats=FALSE))

	### Secondary structure significance - association score (both positive and negative epistasis) - deepcontact
	###########################

	# ss_dt <- rbind(
	# 	result_list[["PWI_cond_1_deepcontact"]][["secondary_structure_score"]][["association_score"]][, region := "1"][, score_type := "association_score"],
	# 	result_list[["PWI_cond_332_deepcontact"]][["secondary_structure_score"]][["association_score"]][, region := "332"][, score_type := "association_score"])
	# ss_dt[, Pos_abs := as.numeric(region)+Pos1-1]

	# #Line plot
	# plot_dt <- copy(ss_dt)[score_type == "association_score",]
	# plot_dt[, alpha_helix_pval := -log10(alpha_p_seed0)]
	# plot_dt[, beta_strand_pval := -log10(beta_p_seed0)]
	# plot_df <- reshape2::melt(as.data.frame(plot_dt[,.(Pos_abs, alpha_helix_pval, beta_strand_pval)]), id = 'Pos_abs')
	# data_min <- min(plot_df[,'value']) + 0.2*min(plot_df[,'value'])
	# data_max <- max(plot_df[,'value']) + 0.2*max(plot_df[,'value'])
	# d <- ggplot2::ggplot(plot_df, ggplot2::aes(Pos_abs, value, color=variable)) +
	#   ggplot2::geom_rect(ggplot2::aes(xmin=min(Pos_abs_hotspot), xmax=max(Pos_abs_hotspot), ymin=data_min, ymax=data_max), fill = "lightgrey", linetype = 0, alpha = 0.1) +
	#   ggplot2::geom_rect(ggplot2::aes(xmin=min(LARKS_coords), xmax=max(LARKS_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	#   ggplot2::geom_rect(ggplot2::aes(xmin=min(helix_coords), xmax=max(helix_coords), ymin=data_min, ymax=data_max), fill = NA, linetype = 1, color = "black") +
	#   ggplot2::geom_vline(xintercept = median(Pos_abs_all), linetype = 2) +
	#   ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2) +
	#   ggplot2::geom_line() +
	#   ggplot2::theme_classic()
 #  if(!is.null(colour_scheme)){
 #    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][c(4, 2)]))
 #  }
	# ggplot2::ggsave(file=file.path(outpath, 'tdp43_ss_pred_association_score_all_deepcontact.pdf'), width=8, height=3, useDingbats=FALSE)


}








