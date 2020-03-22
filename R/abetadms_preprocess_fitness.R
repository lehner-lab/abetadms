
#' abetadms_preprocess_fitness
#'
#' Estimate fitness of doubles mutants using bayesian framework from DiMSum fitness output.
#'
#' @param fitness_path path to DiMSum output (required)
#' @param outpath output path for plots and saved objects (required)
#' @param all_reps list of replicates to retain (required)
#' @param bayesian_double_fitness Estimate fitness of double mutants using bayesian framework (default:F)
#' @param min_mean_input_read_count minimum mean input read count for high confidence variants (default:10)
#' @param min_input_read_count_doubles minimum input read count for doubles used to derive prior for Bayesian doubles correction (default:50)
#' @param lam_d Poisson distribution for score likelihood (default:0.025)
#' @param numCores Number of available CPU cores (default:1)
#' @param execute whether or not to execute the analysis (default: TRUE)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms_preprocess_fitness <- function(
  fitness_path,
  outpath,
  all_reps,
  bayesian_double_fitness = F,
  min_mean_input_read_count = 10,
  min_input_read_count_doubles = 50,
  lam_d = 0.025,
  numCores = 10,
  execute = TRUE
  ){

  #Return previous results if analysis not executed
  if(!execute){
    return()
  }

  #Display status
  message(paste("\n\n*******", "running stage: abetadms_preprocess_fitness", "*******\n\n"))

  #Create output directory
  abetadms__create_dir(abetadms_dir = outpath)

  ### Load fitness and reformat data
  ###########################

  #Load data
  load(fitness_path)
  all_reps_str <- paste0(all_reps, collapse = "")

  dataset_syn <- all_variants
  singles_silent <- dataset_syn[Nham_aa<=1 & is.na(WT),]
  singles_silent[, is.reads0 := TRUE]

  #Add position, mutant AA, WT AA and mean input count
  wt_AAseq <- dataset_syn[WT==T,aa_seq]
  wt_AAseq_split <- strsplit(wt_AAseq,"")[[1]]
  singles_silent[,Pos := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split),aa_seq]
  singles_silent[,Mut := strsplit(aa_seq,"")[[1]][Pos],aa_seq]
  singles_silent[,WT_AA := wt_AAseq_split[Pos],aa_seq]

  #Rename fitness columns
  names(singles_silent) <- gsub("_uncorr$", "", names(singles_silent))
  doubles <- dataset_syn[Nham_aa==2,]
  doubles[, fitness_uncorr := fitness]
  doubles[, sigma_uncorr := sigma]

  #Add position, mutant AA, WT AA and mean input count
  doubles[,Pos1 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[1],aa_seq]
  doubles[,Pos2 := which(strsplit(aa_seq,"")[[1]] !=wt_AAseq_split)[2],aa_seq]
  doubles[,Mut1 := strsplit(aa_seq,"")[[1]][Pos1],aa_seq]
  doubles[,Mut2 := strsplit(aa_seq,"")[[1]][Pos2],aa_seq]
  doubles[,WT_AA1 := wt_AAseq_split[Pos1],aa_seq]
  doubles[,WT_AA2 := wt_AAseq_split[Pos2],aa_seq]
  #Min mean input read count
  doubles[, is.reads0 := TRUE]

  ### Run bayesian fitness estimation for double mutants
  ###########################

  if(bayesian_double_fitness){
    doubles <- abetadms__run_bayesian_double_fitness(
      wt_dt = dataset_syn[WT==T],
      singles_dt = singles_silent[Nham_aa==1],
      doubles_dt = doubles,
      outpath = outpath,
      all_reps = all_reps,
      min_mean_input_read_count = min_mean_input_read_count,
      min_input_read_count_doubles = min_input_read_count_doubles,
      numCores = numCores)
  }else{
    doubles[, fitness_cond := fitness]
    doubles[, sigma_cond := sigma]
    doubles[, fitness1_cond := fitness1_uncorr]
    doubles[, fitness2_cond := fitness2_uncorr]
    doubles[, fitness3_cond := fitness3_uncorr]
    doubles[, fitness4_cond := fitness4_uncorr]
    doubles[, sigma1_cond := sigma1_uncorr]
    doubles[, sigma2_cond := sigma2_uncorr]
    doubles[, sigma3_cond := sigma3_uncorr]
    doubles[, sigma4_cond := sigma4_uncorr]
  }

  ### Output replicate data files
  ###########################

  # define which variants have enough reads
  wildtype[,is.reads0 := TRUE]
  singles_silent[mean_count >= min_mean_input_read_count,is.reads0 := TRUE]
  doubles[mean_count >= min_mean_input_read_count,is.reads0 := TRUE]

  ### Output plain text files
  ###########################

  singles_silent[, Nmut_nt := Nham_nt]
  singles_silent[, Nmut_aa := Nham_aa]
  doubles[, Nmut_nt := Nham_nt]
  doubles[, Nmut_aa := Nham_aa]

  #Save objects
  save(dataset_syn, singles_silent, doubles, file = file.path(outpath, "DMS_processed_data.RData"))

  ##### finalize data.tables
  silent = singles_silent[Nmut_aa==0,.(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,fitness,sigma)]
  singles = singles_silent[Nmut_aa==1,.(Pos,WT_AA,Mut,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,fitness,sigma)]

  #for doubles #add single mutant fitness/sigma values to double mutant table
  doubles[,fitness1 := singles[Pos == Pos1 & Mut == Mut1,fitness],.(Pos1,Mut1)]
  doubles[,sigma1 := singles[Pos == Pos1 & Mut == Mut1,sigma],.(Pos1,Mut1)]
  doubles[,fitness2 := singles[Pos == Pos2 & Mut == Mut2,fitness],.(Pos2,Mut2)]
  doubles[,sigma2 := singles[Pos == Pos2 & Mut == Mut2,sigma],.(Pos2,Mut2)]

  doubles = doubles[,.(Pos1,Pos2,WT_AA1,WT_AA2,Mut1,Mut2,Nmut_nt,Nmut_aa,Nmut_codons,STOP,mean_count,is.reads0,
                       fitness1,sigma1,fitness2,sigma2,
                       fitness_uncorr,sigma_uncorr,
                       fitness_cond,sigma_cond)]

  #Exclude variants with STOP codons from downstream fitness analyses
  wildtype[,is.fitness := TRUE]
  silent[,is.fitness := TRUE]
  singles[,is.fitness := !STOP]
  doubles[,is.fitness := !STOP]

  #write data to files
  write.table(x = wildtype, file = file.path(outpath, "DMS_wildtype.txt"),
              quote = F,row.names = F, col.names = T)
  write.table(x = silent, file = file.path(outpath, "DMS_silent.txt"),
              quote = F,row.names = F, col.names = T)
  write.table(x = singles, file = file.path(outpath, "DMS_singles.txt"),
              quote = F,row.names = F, col.names = T)
  write.table(x = doubles, file = file.path(outpath, "DMS_doubles.txt"),
              quote = F,row.names = F, col.names = T)

}
