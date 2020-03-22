
#' abetadms_human_disease_mutations
#'
#' Check fitness bias of human disease mutations.
#'
#' @param fitness_dt data.table with single mutant fitness values (required)
#' @param outpath output path for plots and saved objects (required)
#' @param missense_AF_file table of missense mutation allele frequencies (required)
#' @param disease_mut_file table of human disease mutations and classifications (required)
#' @param colour_scheme colour scheme file (required)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms_human_disease_mutations <- function(
  fitness_dt,
  missense_AF_file,
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
  message(paste("\n\n*******", "running stage: abetadms_human_disease_mutations", "*******\n\n"))

	#Create output directory
	abetadms__create_dir(abetadms_dir = outpath)

	#Single AA mutants only (no STOPs)
	singles_dt <- copy(fitness_dt[Nmut_aa==1])
	#Double AA mutants only (no STOPs)
	doubles_dt <- copy(fitness_dt[Nmut_aa==2])
	doubles_dt[, fitness := fitness_cond]
	tox_dt <- rbind(singles_dt, doubles_dt, fill = T)

	#Fitness hotspot positions
	mean_fitness <- singles_dt[STOP==F,mean(abs(fitness))]
	Pos_abs_hotspot <- singles_dt[STOP==F,.(hotspot = mean(abs(fitness))>mean_fitness),by=Pos_abs][hotspot==T,Pos_abs]

	#Disease mutations
	dis_mut <- read.table(disease_mut_file, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
	fAD_muts <- rownames(dis_mut)

	#AA code translation dict
	aa_obj <- Biostrings::AAString("GAVLMIFYWKRHDESTCNQP")
	aa_list <- Biostrings::AMINO_ACID_CODE[strsplit(as.character(aa_obj), NULL)[[1]]]
	aa_list["*"] <- "X"
	aa_list_rev <- names(aa_list)
	names(aa_list_rev) <- aa_list

	#Detected human missense mutations
	miss_mut <- as.data.frame(fread(missense_AF_file))
	miss_mut <- miss_mut[nchar(miss_mut[,1])==11,]
	miss_mut[,1] <- gsub("^p.", "", miss_mut[,1])
	rownames(miss_mut) <- paste0(
		aa_list_rev[substr(miss_mut[,1], 1, 3)], 
		as.integer(substr(miss_mut[,1], 4, 6))-672+1, 
		aa_list_rev[substr(miss_mut[,1], 7, 9)])

	#Add mutation information
	tox_dt[, hmut_cat := "Never observed (gnomAD)"]
	tox_dt[!is.na(miss_mut[mut_code,2]), hmut_cat := "Observed (gnomAD)"]
	tox_dt[mut_code %in% fAD_muts, hmut_cat := "fAD"]

	fwrite(tox_dt[mut_code %in% rownames(dis_mut),], file = file.path(outpath, "dis_mut_tab.tsv"), sep = "\t")

	### Fitness bias of TDP-43 mutations
	###########################

	#Z-test individually
	tox_bias_ind <- tox_dt[hmut_cat %in% c("fAD"),.(mut_code, p_value = pnorm(fitness/sigma, lower.tail=FALSE)*2)][order(p_value, decreasing = F)]

	#Number of variants with defined fitness values in each replicate
	print(paste0("Aggregation propensity of all human disease mutations:"))
	print(tox_bias_ind)

	# #t-test on whole sample
	# tox_dt[hmut_cat %in% c("fAD"),t.test(fitness)]

	#Z-test on whole sample
	fitness_merged <- tox_dt[hmut_cat %in% c("fAD"),sum( fitness / (sigma^2) ) / sum( 1 / (sigma^2) )]
	sigma_merged <- tox_dt[hmut_cat %in% c("fAD"),sqrt( 1 / sum( 1 / (sigma^2) ) )]
	z_score <- (fitness_merged-0)/sigma_merged
	p_value <- pnorm(z_score, lower.tail = FALSE)*2
	tox_bias_all <- data.table(pvalue = p_value, n = unlist(tox_bias_ind[,.N]))

	# #Z-test (reference)
	# fitness_merged_ref <- tox_dt[hmut_cat %in% c("Observed (gnomAD)"),sum( fitness / (sigma^2) ) / sum( 1 / (sigma^2) )]
	# sigma_merged_ref <- tox_dt[hmut_cat %in% c("Observed (gnomAD)"),sqrt( 1 / sum( 1 / (sigma^2) ) )]
	# z_score_ref <- ( fitness_merged - fitness_merged_ref ) / sqrt( (sigma_merged^2) + (sigma_merged_ref^2) )
	# pnorm(z_score_ref, lower.tail = FALSE)*2

	#Fitness histogram of all single and double observed versus unobserved human missense mutations
	# tox_bias_all <- tox_dt[hmut_cat %in% c("fAD"),.(pvalue = wilcox.test(fitness)$p.value, n = .N)]
	set.seed(1)
	plot_df <- as.data.frame(tox_dt[,c("fitness", "Nmut_aa", "STOP", "hmut_cat")])
	plot_df[,"hmut_cat"] <- factor(plot_df[,"hmut_cat"], levels = c("fAD", "Never observed (gnomAD)", "Observed (gnomAD)"))
  d <- ggplot2::ggplot(plot_df, ggplot2::aes(fitness, ..density..)) +
    ggplot2::geom_density() +
    ggplot2::geom_jitter(data = plot_df[!grepl("gnomAD", plot_df[,"hmut_cat"]),], ggplot2::aes(x = fitness, y = fitness*0-1, color = hmut_cat)) +
    ggplot2::xlab("Fitness") +
    ggplot2::geom_vline(xintercept = 0, linetype=2) +
    ggplot2::geom_vline(xintercept = tox_dt[Nmut_aa==1 & STOP==T,median(fitness)], linetype=2, colour = "darkgrey") +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_manual(values = unlist(colour_scheme[["shade 0"]][1:4])) +
    ggplot2::scale_shape_manual(values = c(1, 19)) + 
    ggplot2::annotate("text", label = paste0("P-value (all) = ", format(tox_bias_all[,"pvalue"], digits = 2, scientific = T), " (", tox_bias_all[,"n"], ")") , x = 0.2, y = -0.5)
  ggplot2::ggsave(file=file.path(outpath, 'human_disease_mut_fitness.pdf'), width=5, height=3, useDingbats=FALSE)

	#Test human mutation fitness bias
  # wilcox.test(tox_dt[hmut_cat %in% c("fALS", "sALS", "fALS and sALS"),fitness])
  # wilcox.test(tox_dt[hmut_cat=="fALS" & hmut_recurrent==T,fitness])
  # wilcox.test(tox_dt[hmut_cat=="fALS" & hmut_recurrent==F,fitness])
  # wilcox.test(tox_dt[hmut_cat=="sALS" & hmut_recurrent==T,fitness])
  # wilcox.test(tox_dt[hmut_cat=="sALS" & hmut_recurrent==F,fitness])
  # wilcox.test(tox_dt[hmut_cat=="fALS and sALS" & hmut_recurrent==T,fitness])
  # wilcox.test(tox_dt[hmut_cat=="fALS and sALS" & hmut_recurrent==F,fitness])
  # wilcox.test(tox_dt[hmut_cat %in% c("fALS", "sALS", "fALS and sALS") & !(hmut_cat=="fALS" & hmut_recurrent==T),fitness])
}

