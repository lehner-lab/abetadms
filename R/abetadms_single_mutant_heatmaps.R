
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
		cluster="column")

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

	#Fitness (ALS mutations)
	abetadms__single_mutant_heatmap(
		input_df = as.data.frame(dms_dt), 
		variable_name="fitness", 
		mut_dict=list("F" = fAD_muts),
		output_file = file.path(outpath, 'heatmap_fitness_center_scale_hmuts_both.pdf'),
		x_annotation="both", 
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

	### K-mediods clustering of residues based on singles
	###########################

	#Determine optimimum number of clusters (sample each fitness value from error distribution) - average silhouette width (for K = 1-10)
	set.seed(1)
	pamk_list <- list()
	for(i in 1:100){
		dms_dt_rand <- copy(fitness_dt[Nmut_aa==1])
		dms_dt_rand[,fitness := rnorm(1, mean = fitness, sd = sigma),mut_code]
		heat_mat <- abetadms__datatable_to_matrix(dms_dt_rand)
		d <- dist(t(heat_mat), method = "euclidean") # distance matrix
		pamk_list[[i]] <- fpc::pamk(d,krange=1:10)
	}
	mean_sil_dt <- as.data.table(do.call("rbind", sapply(pamk_list, '[', "crit")))
	names(mean_sil_dt) <- paste0("K=", 1:10)
	#Plot
	plot_dt <- reshape2::melt(mean_sil_dt, measure.vars = 1:10)
	d <- ggplot2::ggplot(plot_dt, ggplot2::aes(variable, value)) +
	  ggplot2::geom_boxplot() +
	  ggplot2::theme_bw() + 
	  ggplot2::ylab("Average Silhouette width") +
	  ggplot2::xlab("Number of clusters")
	ggplot2::ggsave(file=file.path(outpath, 'kmedoids_average_silhouette_width_boxplot.pdf'), width=5, height=5, useDingbats=FALSE)

	#Silhouette plot (K=2)
	heat_mat <- abetadms__datatable_to_matrix(dms_dt)
	d <- dist(t(heat_mat), method = "euclidean") # distance matrix
	set.seed(1)
	clust_obj <- fpc::pamk(d,krange=1:10)
	#Plot
	plot_dt <- as.data.table(clust_obj[["pamobject"]][["silinfo"]][["widths"]])
	plot_dt[, Pos := rownames(clust_obj[["pamobject"]][["silinfo"]][["widths"]])]
	plot_dt[, Pos := factor(Pos, levels = Pos)]
	plot_dt[, cluster := factor(cluster)]
	pos_order <- plot_dt[,Pos]
	d <- ggplot2::ggplot(plot_dt, ggplot2::aes(Pos, sil_width, fill = cluster)) +
	  ggplot2::geom_bar(stat = "identity") +
	  ggplot2::theme_bw() +
	  ggplot2::ylab("Silhouette width") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
	ggplot2::ggsave(file=file.path(outpath, 'kmedoids_silhouette_width_barplot.pdf'), width=7, height=5, useDingbats=FALSE)

	#Silhouette plot (K=2)
	sil_list <- sapply(sapply(sapply(pamk_list, '[', "pamobject"), '[', "silinfo"), '[', "widths")
	sil_df <- do.call("rbind", sil_list)
	sil_dt <- as.data.table(sil_df)
	sil_dt[, Pos := factor(rownames(sil_df), levels = pos_order)]
	sil_dt[, cluster := factor(cluster)]
	d <- ggplot2::ggplot(sil_dt, ggplot2::aes(Pos, sil_width, fill = cluster)) +
	  ggplot2::geom_boxplot() +
	  ggplot2::theme_bw() + 
	  ggplot2::ylab("Silhouette width") +
	  ggplot2::xlab("Pos") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::facet_grid(cluster~., scales = "free_x")
	ggplot2::ggsave(file=file.path(outpath, 'kmedoids_silhouette_width_boxplots.pdf'), width=7, height=7, useDingbats=FALSE)


	# ### Heatmap of average single effects from doubles
	# ###########################

	# #Subset to double AA mutants only
	# dms_dt1 <- copy(fitness_dt[Nmut_aa==2])
	# dms_dt2 <- copy(fitness_dt[Nmut_aa==2])

	# #Pos and Mut
	# dms_dt1[, Pos_abs := Pos1]
	# dms_dt1[, Mut := Mut1]
	# dms_dt1[, WT_AA := WT_AA1]
	# dms_dt2[, Pos_abs := Pos2]
	# dms_dt2[, Mut := Mut2]
	# dms_dt2[, WT_AA := WT_AA2]

	# #Combine
	# dms_dt <- rbind(dms_dt1, dms_dt2)

	# #Mean effect of singles
	# double_means <- dms_dt[STOP==F,.(fitness = mean(fitness_cond)),.(Pos_abs,Mut,WT_AA)]
	# double_means[,Nmut_aa := 1]

	# #Fitness
	# abetadms__single_mutant_heatmap(
	# 	input_df = as.data.frame(double_means), 
	# 	variable_name="fitness", 
	# 	output_file = file.path(outpath, 'heatmap_fitness_meandoubles.pdf'),
	# 	x_annotation="sequence", 
	# 	xaxis_angle=0, 
	# 	na_colour="white", 
	# 	na_text="-", 
	# 	colour_low=colour_scheme[["shade 0"]][[3]], 
	# 	colour_high=colour_scheme[["shade 0"]][[1]])


}

