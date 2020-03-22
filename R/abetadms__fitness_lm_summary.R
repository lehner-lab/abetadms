
#' abetadms__fitness_lm_summary
#'
#' Dot plots showing explained variance of models to predict variant fitness.
#'
#' @param input_df data.frame with fitness and predictor variables (required)
#' @param output_file path for output plot file (required)
#' @param return_results whether or not to return rsq_list (default:F)
#' @param colour_scheme colour scheme list
#' @param subset data subset (default:none)
#' @param width plot width (default:6)
#' @param height plot height (default:6)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms__fitness_lm_summary <- function(
  input_df,
  output_file,
  return_results = F,
  colour_scheme = NULL,
  subset = NULL,
  width = 6, 
  height = 6){
	#Subset
  if(!is.null(subset)){
    input_df <- input_df[subset,]
  }
  #Normalise
  input_df <- abetadms__glm_prenorm(input_df)
  #Initialise
  rsq_list <- list()
  for(this_feat in colnames(input_df)[-grep("Position|Fitness", colnames(input_df))]){
    #Both inside/outside with hotspot variable
    temp_feats <- c(colnames(input_df)[grep("Position|Fitness", colnames(input_df))], this_feat)
    fitness_lm <- lm(Fitness~.*., data = input_df[,colnames(input_df) %in% temp_feats])
    rsq_list[[paste0(this_feat, "_pos")]] <- summary(fitness_lm)$adj.r.squared
  	#Both inside/outside 
    temp_feats <- c(colnames(input_df)[grep("Fitness", colnames(input_df))], this_feat)
    fitness_lm <- lm(Fitness~.*., data = input_df[,colnames(input_df) %in% temp_feats])
    rsq_list[[this_feat]] <- summary(fitness_lm)$adj.r.squared
    #Inside hotspot
    fitness_lm <- lm(Fitness~.*., data = input_df[input_df$Position==max(input_df$Position),colnames(input_df) %in% temp_feats])
    rsq_list[[paste0(this_feat, "_inside")]] <- summary(fitness_lm)$adj.r.squared
    #Outside hotspot
    fitness_lm <- lm(Fitness~.*., data = input_df[input_df$Position==min(input_df$Position),colnames(input_df) %in% temp_feats])
    rsq_list[[paste0(this_feat, "_outside")]] <- summary(fitness_lm)$adj.r.squared
  }

  #All physical/chemical properties features and position and all first order interactions
  temp_feats <- c(colnames(input_df)[grep("Fitness|^PC|Position", colnames(input_df))])
  fitness_lm <- lm(Fitness~.*., data = input_df[,colnames(input_df) %in% temp_feats])
  rsq_list[["All PCs_pos"]] <- summary(fitness_lm)$adj.r.squared
  #All physical/chemical properties features and all first order interactions
  temp_feats <- c(colnames(input_df)[grep("Fitness|^PC", colnames(input_df))])
  fitness_lm <- lm(Fitness~.*., data = input_df[,colnames(input_df) %in% temp_feats])
  rsq_list[["All PCs"]] <- summary(fitness_lm)$adj.r.squared
  #All physical/chemical properties features and all first order interactions - inside hotspot
  fitness_lm <- lm(Fitness~.*., data = input_df[input_df$Position==max(input_df$Position),colnames(input_df) %in% temp_feats])
  rsq_list[["All PCs_inside"]] <- summary(fitness_lm)$adj.r.squared
  #All physical/chemical properties features and all first order interactions - outside hotspot
  fitness_lm <- lm(Fitness~.*., data = input_df[input_df$Position==min(input_df$Position),colnames(input_df) %in% temp_feats])
  rsq_list[["All PCs_outside"]] <- summary(fitness_lm)$adj.r.squared

  #All physical/chemical properties features and position and all first order interactions
  temp_feats <- c(colnames(input_df)[-grep("^PC", colnames(input_df))])
  fitness_lm <- lm(Fitness~.*., data = input_df[,colnames(input_df) %in% temp_feats])
  rsq_list[["All Aggregation tools_pos"]] <- summary(fitness_lm)$adj.r.squared
  #All aggregation tool features and all first order interactions
  temp_feats <- c(colnames(input_df)[-grep("^PC|Position", colnames(input_df))])
  fitness_lm <- lm(Fitness~.*., data = input_df[,colnames(input_df) %in% temp_feats])
  rsq_list[["All Aggregation tools"]] <- summary(fitness_lm)$adj.r.squared
  #All physical/chemical properties features and all first order interactions - inside hotspot
  fitness_lm <- lm(Fitness~.*., data = input_df[input_df$Position==max(input_df$Position),colnames(input_df) %in% temp_feats])
  rsq_list[["All Aggregation tools_inside"]] <- summary(fitness_lm)$adj.r.squared
  #All physical/chemical properties features and all first order interactions - outside hotspot
  fitness_lm <- lm(Fitness~.*., data = input_df[input_df$Position==min(input_df$Position),colnames(input_df) %in% temp_feats])
  rsq_list[["All Aggregation tools_outside"]] <- summary(fitness_lm)$adj.r.squared

  #All physical/chemical properties and aggregation tool features and position and all first order interactions
  fitness_lm <- lm(Fitness~.*., data = input_df)
  rsq_list[["All features_pos"]] <- summary(fitness_lm)$adj.r.squared
  #All physical/chemical properties and aggregation tool features and all first order interactions
  temp_feats <- c(colnames(input_df)[-grep("Position", colnames(input_df))])
  fitness_lm <- lm(Fitness~.*., data = input_df[,colnames(input_df) %in% temp_feats])
  rsq_list[["All features"]] <- summary(fitness_lm)$adj.r.squared
  #All physical/chemical properties and aggregation tool features and all first order interactions - inside hotspot
  fitness_lm <- lm(Fitness~.*., data = input_df[input_df$Position==max(input_df$Position),colnames(input_df) %in% temp_feats])
  rsq_list[["All features_inside"]] <- summary(fitness_lm)$adj.r.squared
  #All physical/chemical properties and aggregation tool features and all first order interactions - outside hotspot
  fitness_lm <- lm(Fitness~.*., data = input_df[input_df$Position==min(input_df$Position),colnames(input_df) %in% temp_feats])
  rsq_list[["All features_outside"]] <- summary(fitness_lm)$adj.r.squared

  #Plot
  plot_df <- data.frame(rsquared = unlist(rsq_list), feature = gsub("_pos$|_inside$|_outside$", "", names(rsq_list)), stringsAsFactors = F)
  plot_df$position <- "both"
  plot_df$position[grepl("_pos$", names(rsq_list))] <- "both (with hotspot variable)"
  plot_df$position[grepl("_inside$", names(rsq_list))] <- "inside hotspot"
  plot_df$position[grepl("_outside$", names(rsq_list))] <- "outside hotspot"
  plot_df$position <- factor(plot_df$position, levels = c(
  	"inside hotspot", 
  	"outside hotspot",
  	"both", 
  	"both (with hotspot variable)"))
  plot_df$feature_category <- "Aggregation tools"
  plot_df$feature_category[grep("^PC", plot_df$feature)] <- "Amino acid properties"
  plot_df$feature_category[grep("^All", plot_df$feature)] <- "Combinations"
  plot_df$feature_category <- factor(plot_df$feature_category, levels = c("Combinations", "Aggregation tools", "Amino acid properties"))
  plot_df$feature <- factor(plot_df$feature, levels = unique(plot_df$feature[order(plot_df$feature_category, plot_df$rsquared, decreasing = T)]))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(feature, rsquared, color=position)) +
    ggplot2::geom_point() +
    ggplot2::ylab("% Variance in Fitness Explained") +
    ggplot2::xlab("Feature set") +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::facet_wrap(~ feature_category, scales = "free_x")
  if(!is.null(colour_scheme)){
    d <- d + ggplot2::scale_colour_manual(values = unlist(colour_scheme, use.names = FALSE)[1:4])
  }
  suppressWarnings(suppressMessages(ggplot2::ggsave(file=output_file, width=width, height=height, useDingbats=FALSE)))
  if(return_results){
    return(rsq_list)
  }
}

