
PWI_datatable_to_matrix <- function(input_dt, feature_name = "epistasis_score", score_threshold = 0.1){
  input_df <- as.data.frame(input_dt)
  #Calculate absolute position
  if(!"Pos1_abs" %in% colnames(input_df) | !"Pos2_abs" %in% colnames(input_df)){
    input_df$Pos1_abs <- NA
    input_df$Pos1_abs <- input_df$Pos1+as.numeric(input_df$region)-1
    input_df$Pos2_abs <- NA
    input_df$Pos2_abs <- input_df$Pos2+as.numeric(input_df$region)-1
  }
  #Position range
  pos_range <- range(input_df[,c("Pos1_abs","Pos2_abs")])
  epi_mat <- matrix(NA, nrow = pos_range[2]-pos_range[1]+1, ncol = pos_range[2]-pos_range[1]+1)
  colnames(epi_mat) <- pos_range[1]:pos_range[2]
  rownames(epi_mat) <- pos_range[1]:pos_range[2]
  for(i in pos_range[1]:(pos_range[2])){
    if(i!=pos_range[2]){
      for(j in (i+1):pos_range[2]){
        if(sum(input_df$Pos1_abs==i & input_df$Pos2_abs==j)==1){
          temp_score <- as.numeric(input_df[input_df$Pos1_abs==i & input_df$Pos2_abs==j,feature_name])
          if(temp_score >= score_threshold){
            epi_mat[as.character(i),as.character(j)] <- temp_score
          }else{
            epi_mat[as.character(i),as.character(j)] <- 0
          }
          epi_mat[as.character(j),as.character(i)] <- epi_mat[as.character(i),as.character(j)]
        }
      }
    }
  }
  return(epi_mat)
}

#' abetadms_PWI_heatmaps
#'
#' Plot pair-wise interacton matrix heatmaps.
#'
#' @param PWI_dir_list list of PWI directories (required)
#' @param outpath output path for plots and saved objects (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms_PWI_heatmaps <- function(
  PWI_dir_list,
  outpath,
  colour_scheme,
  execute = TRUE
  ){

	#Return previous results if analysis not executed
	if(!execute){
		return()
	}

  #Display status
  message(paste("\n\n*******", "running stage: abetadms_PWI_heatmaps", "*******\n\n"))

	#Create output directory
	abetadms__create_dir(abetadms_dir = outpath)

	### Setup
	###########################

	#Load PWI data
	PWI_cond_list <- list()
	for(d in names(PWI_dir_list)){
		for(f in list.files(PWI_dir_list[[d]], pattern = "DMS_PWI.txt")){
		  PWI_cond_list[[paste0("PWI_cond_", d, gsub('DMS_PWI|.txt', '', f))]] <- fread(file.path(PWI_dir_list[[d]], f))[,region := d]
		}
	}
	#Merge
	pwi_dt <- do.call('rbind', PWI_cond_list)
	#Absolute position
  # pwi_dt[, Pos_abs1 := as.numeric(region)+Pos1-1]
  # pwi_dt[, Pos_abs2 := as.numeric(region)+Pos2-1]
  pwi_dt[, Pos_abs1 := Pos1]
  pwi_dt[, Pos_abs2 := Pos2]

	# #Structure coordinates
	# LARKS_coords <- c(312, 317)
	# helix_coords <- c(321, 330)

	### PWI heatmaps
	###########################

  for(i in c('Epistasis', 'Association', 'Combined')){
    num_top_ints <- 10
    plot_mat <- PWI_datatable_to_matrix(
      input_dt = pwi_dt, 
      feature_name = paste0(tolower(i), "_score"), 
      score_threshold = 0)
    # top_thresh <- plot_mat[lower.tri(plot_mat)]
    # top_thresh <- top_thresh[order(top_thresh, decreasing = T)][num_top_ints]
    input_matrix_text <- replace(plot_mat, T, '')
    input_matrix_text[c(28, 42), c(28, 42)] <- '*'
    input_matrix_text[c(23, 28), c(23, 28)] <- '*'
    input_matrix_text[c(3, 28), c(3, 28)] <- '*'
    abetadms__tile_heatmap_wrapper(
      input_matrix = plot_mat,
      input_matrix_text = input_matrix_text, text_size = 5, text_colour = colour_scheme[["shade 0"]][1],
      output_file = file.path(outpath, paste0("heatmap_cond_", tolower(i), "_score.pdf")), 
      width=7, height=7, colour_low=colour_scheme[["shade 0"]][1], colour_high=colour_scheme[["shade 0"]][3], colour_clip=FALSE, 
      cluster='none', xaxis_angle=90, xaxis_hjust=0.5, xaxis_size=NULL, 
      plot_title = paste0(i, " score"), 
      xlab = "AB amino acid position", ylab = "AB amino acid position")
  }

  for(i in c('posE_enr', 'negE_enr')){
    num_top_ints <- 10
    plot_mat <- PWI_datatable_to_matrix(
      input_dt = pwi_dt, 
      feature_name = i, 
      score_threshold = 0)
    # top_thresh <- plot_mat[lower.tri(plot_mat)]
    # top_thresh <- top_thresh[order(top_thresh, decreasing = T)][num_top_ints]
    input_matrix_text <- replace(plot_mat, T, '')
    input_matrix_text[c(28, 42), c(28, 42)] <- '*'
    input_matrix_text[c(23, 28), c(23, 28)] <- '*'
    input_matrix_text[c(3, 28), c(3, 28)] <- '*'
    abetadms__tile_heatmap_wrapper(
      input_matrix = plot_mat,
      input_matrix_text = input_matrix_text, text_size = 5, text_colour = colour_scheme[["shade 0"]][1],
      output_file = file.path(outpath, paste0("heatmap_cond_", i, ".pdf")), 
      width=7, height=7, colour_low=colour_scheme[["shade 0"]][1], colour_high=colour_scheme[["shade 0"]][3], colour_clip=FALSE, 
      cluster='none', xaxis_angle=90, xaxis_hjust=0.5, xaxis_size=NULL, 
      plot_title = paste0(i, " score"), 
      xlab = "AB amino acid position", ylab = "AB amino acid position")
  }

	# #Plot - helix positions - association score
 #  num_top_ints <- 7
	# plot_mat <- PWI_datatable_to_matrix(
	# 	input_dt = pwi_dt[Pos_abs1 %in% seq(helix_coords[1],helix_coords[2]) & Pos_abs2 %in% seq(helix_coords[1],helix_coords[2]),], 
	# 	feature_name = "association_score", 
	# 	score_threshold = 0)
 #  top_thresh <- plot_mat[lower.tri(plot_mat)]
 #  top_thresh <- top_thresh[order(top_thresh, decreasing = T)][num_top_ints]
 #  input_matrix_text <- replace(plot_mat, T, '')
 #  input_matrix_text[which(plot_mat>=top_thresh)] <- '*'
	# abetadms__tile_heatmap_wrapper(
 #    input_matrix = plot_mat,
 #    input_matrix_text = input_matrix_text, text_size = 10, text_colour = colour_scheme[["shade 0"]][1],
 #    output_file = file.path(outpath, paste0("heatmap_cond_association_score", "_helix", ".pdf")), 
 #    width=5.5, height=5, colour_low=colour_scheme[["shade 0"]][1], colour_high=colour_scheme[["shade 0"]][3], colour_clip=0.7, 
 #    cluster='none', xaxis_angle=90, xaxis_hjust=0.5, xaxis_size=NULL, 
 #    plot_title = "Association enrichment", 
 #    xlab = "TDP-43 amino acid position", ylab = "TDP-43 amino acid position")

	# #Plot - LARKS positions - association score
 #  num_top_ints <- 7
	# plot_mat <- PWI_datatable_to_matrix(
	# 	input_dt = pwi_dt[Pos_abs1 %in% seq(LARKS_coords[1],LARKS_coords[2]) & Pos_abs2 %in% seq(LARKS_coords[1],LARKS_coords[2]),], 
	# 	feature_name = "association_score", 
	# 	score_threshold = 0)
 #  top_thresh <- plot_mat[lower.tri(plot_mat)]
 #  top_thresh <- top_thresh[order(top_thresh, decreasing = T)][num_top_ints]
 #  input_matrix_text <- replace(plot_mat, T, '')
 #  input_matrix_text[which(plot_mat>=top_thresh)] <- '*'
	# abetadms__tile_heatmap_wrapper(
 #    input_matrix = plot_mat,
 #    input_matrix_text = input_matrix_text, text_size = 10, text_colour = colour_scheme[["shade 0"]][1],
 #    output_file = file.path(outpath, paste0("heatmap_cond_association_score", "_LARKS", ".pdf")), 
 #    width=5.5, height=5, colour_low=colour_scheme[["shade 0"]][1], colour_high=colour_scheme[["shade 0"]][3], colour_clip=0.7, 
 #    cluster='none', xaxis_angle=90, xaxis_hjust=0.5, xaxis_size=NULL, 
 #    plot_title = "Association enrichment", 
 #    xlab = "TDP-43 amino acid position", ylab = "TDP-43 amino acid position")

  # #Plot - helix and LARKS positions - epistasis score - positive
  # num_top_ints <- 7
  # plot_mat <- PWI_datatable_to_matrix(
  #   input_dt = pwi_dt[Pos_abs1 %in% seq(LARKS_coords[1],helix_coords[2]) & Pos_abs2 %in% seq(LARKS_coords[1],helix_coords[2]),], 
  #   feature_name = "posE_enr", 
  #   score_threshold = 0)
  # top_thresh <- plot_mat[lower.tri(plot_mat)]
  # top_thresh <- top_thresh[order(top_thresh, decreasing = T)][num_top_ints]
  # input_matrix_text <- replace(plot_mat, T, '')
  # input_matrix_text[which(plot_mat>=top_thresh)] <- '*'
  # abetadms__tile_heatmap_wrapper(
  #   input_matrix = plot_mat,
  #   input_matrix_text = input_matrix_text, text_size = 10, text_colour = colour_scheme[["shade 0"]][1],
  #   output_file = file.path(outpath, paste0("heatmap_cond_epistasis_score_positive", "_helix_and_LARKS", ".pdf")), 
  #   width=5.5, height=5, colour_low=colour_scheme[["shade 0"]][1], colour_high=colour_scheme[["shade 0"]][3], colour_midpoint=1.5, 
  #   cluster='none', xaxis_angle=90, xaxis_hjust=0.5, xaxis_size=NULL, 
  #   plot_title = "Epistasis enrichment - positive", 
  #   xlab = "TDP-43 amino acid position", ylab = "TDP-43 amino acid position")

  # #Plot - helix and LARKS positions - epistasis score - positive
  # num_top_ints <- 7
  # plot_mat <- PWI_datatable_to_matrix(
  #   input_dt = pwi_dt[Pos_abs1 %in% seq(LARKS_coords[1],helix_coords[2]) & Pos_abs2 %in% seq(LARKS_coords[1],helix_coords[2]),], 
  #   feature_name = "negE_enr", 
  #   score_threshold = 0)
  # top_thresh <- plot_mat[lower.tri(plot_mat)]
  # top_thresh <- top_thresh[order(top_thresh, decreasing = T)][num_top_ints]
  # input_matrix_text <- replace(plot_mat, T, '')
  # input_matrix_text[which(plot_mat>=top_thresh)] <- '*'
  # abetadms__tile_heatmap_wrapper(
  #   input_matrix = plot_mat,
  #   input_matrix_text = input_matrix_text, text_size = 10, text_colour = colour_scheme[["shade 0"]][1],
  #   output_file = file.path(outpath, paste0("heatmap_cond_epistasis_score_negative", "_helix_and_LARKS", ".pdf")), 
  #   width=5.5, height=5, colour_low=colour_scheme[["shade 0"]][1], colour_high=colour_scheme[["shade 0"]][3], colour_midpoint=1.5, 
  #   cluster='none', xaxis_angle=90, xaxis_hjust=0.5, xaxis_size=NULL, 
  #   plot_title = "Epistasis enrichment - negative", 
  #   xlab = "TDP-43 amino acid position", ylab = "TDP-43 amino acid position")


}

