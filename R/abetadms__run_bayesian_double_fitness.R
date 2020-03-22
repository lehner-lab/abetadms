
#' abetadms__run_bayesian_double_fitness
#'
#' Estimate fitness of double mutants using bayesian framework.
#'
#' @param wt_dt data.table with WT (required)
#' @param singles_dt data.table with double AA mutants (required)
#' @param doubles_dt data.table with double AA mutants (required)
#' @param outpath output path for plots and saved objects (required)
#' @param all_reps list of replicates to retain (required)
#' @param min_mean_input_read_count minimum mean input read count for high confidence variants (default:10)
#' @param min_input_read_count_doubles minimum input read count for doubles used to derive prior for Bayesian doubles correction (default:50)
#' @param lam_d Poisson distribution for score likelihood (default:0.025)
#' @param numCores Number of available CPU cores (default:1)
#'
#' @return Nothing
#' @export
#' @import data.table
abetadms__run_bayesian_double_fitness <- function(
  wt_dt,
  singles_dt,
  doubles_dt,
  outpath,
  all_reps,
  min_mean_input_read_count = 10,
  min_input_read_count_doubles = 50,
  lam_d = 0.025,
  numCores = 10
  ){

  #Display status
  message(paste("\n\n*******", "Estimating fitness of double mutants using bayesian framework (this might take a while; you may want to adjust 'numCores' argument, DEFAULT=10)", "*******\n\n"))

  ### Globals
  ###########################

  all_reps_str <- paste0(all_reps, collapse = "")

  ### Bayesian framework for fitness estimation for double mutants
  ###########################

  #Bin mean counts for replicate 1
  doubles_dt[,counts_for_bins := .SD[[1]],,.SDcols = paste0("count_e",all_reps[1],"_s0")]
  doubles_dt[,bin_count := findInterval(log10(counts_for_bins),seq(0.5,4,0.25))]
  doubles_dt[,.(.N,mean(counts_for_bins)),bin_count][order(bin_count)]

  #Plot fitness densities for different mean count bins (replicate 1)
  doubles_dt[,fitness_for_bins := .SD[[1]],,.SDcols = paste0("fitness",all_reps[1],"_uncorr")]
  d <- ggplot2::ggplot(doubles_dt[between(bin_count,2,8)],ggplot2::aes(fitness_for_bins,..density..,color=factor(bin_count))) +
    ggplot2::geom_density(adjust=1)
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework1.pdf"), d, width = 7, height = 5, useDingbats=FALSE)

  #Plot fitness densities for mean counts greater/less than 50 (replicate 1)
  d <- ggplot2::ggplot(doubles_dt,ggplot2::aes(fitness_for_bins,..density..,color=counts_for_bins >= min_input_read_count_doubles)) +
    ggplot2::geom_density(adjust=1)
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework2.pdf"), d, width = 7, height = 5, useDingbats=FALSE)
  #>> try to estimate what fitness scores are for variants with low sequence coverage
  # use double mutants with variants >= min_input_read_count_doubles counts 

  # save.image(file = file.path(outpath, "Rsession1.RData"))

  ## calculate posterior double mutant fitness based on prior from single mutants
  postpois_conditioned_singleF <- function(i){  
    require(data.table)
    count_in = double_data[i,count_in]
    count_out = double_data[i,count_out]
    lam_in = exp(seq(floor(log(count_in+0.1)-max(c(0.5,1/log10(count_in+1.75)))),(log(count_in+0.1)+max(c(0.5,1/log10(count_in+1.75)))),lam_d))
    lam_out = exp(seq(floor(log(count_out+0.1)-max(c(0.5,1/log10(count_out+1.75)))),(log(count_out+0.1)+max(c(0.5,1/log10(count_out+1.75)))),lam_d))
    lam_low = range(log(lam_out))[1] - range(log(lam_in))[2]
    lam_high = range(log(lam_out))[2] - range(log(lam_in))[1]
    idx = row(matrix(NA,nrow=length(lam_out),ncol=length(lam_in))) - col(matrix(NA,nrow=length(lam_out),ncol=length(lam_in)))
    likelihood = sapply(split(outer(dpois(count_out,lambda = lam_out),dpois(count_in,lambda = lam_in)),idx),sum)
    score_prior = density(score_prior_cond[,.(fdist = sqrt((double_data[i,F1]-F1)^2+(double_data[i,F2]-F2)^2),F)][
      order(fdist)][1:Nneighbours,F],
      from = (lam_low-wt_corr),
      to = (lam_high-wt_corr),
      n = as.integer(as.character(round((lam_high-lam_low)/lam_d + 1)))) #super weird bug
    posterior = score_prior$y*likelihood
    
    moments = list()
    moments[1] = weighted.mean(x = score_prior$x,w = posterior)
    moments[2] = sqrt(sum(( moments[[1]]-score_prior$x)^2 * posterior)/
                        sum(posterior))
    return(moments)
  }

  # Setup cluster
  clust <- parallel::makeCluster(numCores) #This line will take time

  #Calculate conditional fitness and sigma
  for (E in all_reps) {
    #wildtype "correction" to calculate scores
    wt_corr <- wt_dt[,log(unlist(.SD[,1]) / unlist(.SD[,2])),,.SDcols = c(paste0("count_e",E,"_s1"),paste0("count_e",E,"_s0"))]
    #data for prior calculation
    double_data <- doubles_dt[!is.na(get(paste0("fitness",E,"_uncorr"))),.(Pos1,Mut1,Pos2,Mut2,count_in = unlist(.SD[,1]),count_out = unlist(.SD[,2]),
                             F = unlist(.SD[,3])),,
                          .SDcols = c(paste0("count_e",E,"_s0"),paste0("count_e",E,"_s1"),paste0("fitness",E,"_uncorr"))]
    # double_data = merge(double_data,singles_dt[,.(Pos,Mut,F1 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
    double_data <- merge(double_data,singles_dt[!is.na(singles_dt[,paste0("fitness",E)]) & is.reads0==T,.(Pos,Mut,F1 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos1","Mut1"),by.y = c("Pos","Mut"))
    # double_data = merge(double_data,singles_dt[,.(Pos,Mut,F2 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
    double_data <- merge(double_data,singles_dt[!is.na(singles_dt[,paste0("fitness",E)]) & is.reads0==T,.(Pos,Mut,F2 = .SD),,.SDcols = paste0("fitness",E)],by.x = c("Pos2","Mut2"),by.y = c("Pos","Mut"))
    
    Nneighbours <- 500
    score_prior_cond <- double_data[count_in >= min_input_read_count_doubles & F > -Inf & F1 > -Inf & F2 > -Inf]

    # make variables available to each core's workspace
    parallel::clusterExport(clust, list("double_data","lam_d","wt_corr","score_prior_cond","Nneighbours"), envir = environment())

    #posterior fitness conditioned on single fitness
    t=proc.time()
    helper <- parallel::parSapply(clust,X = 1:nrow(double_data), postpois_conditioned_singleF)
    print(proc.time()-t)
    helper1 <- matrix(unlist(helper),nrow=2)
    double_data[,paste0("fitness",E,"_cond") := helper1[1,]]
    double_data[,paste0("sigma",E,"_cond") := helper1[2,]]
    doubles_dt <- merge(doubles_dt, double_data[,.SD,,.SDcols = c("Pos1", "Pos2", "Mut1", "Mut2", paste0("fitness",E,"_cond"), paste0("sigma",E,"_cond"))], by = c("Pos1", "Pos2", "Mut1", "Mut2"), all.x = T)
  }
  parallel::stopCluster(clust)

  #Scatterplot matrix - singles
  d <- GGally::ggpairs(singles_dt[Nham_aa==1,grep(names(singles_dt),pattern="fitness"),with=F],
          upper=list(continuous = "cor"))
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework_scattermatrix_singles.pdf"), d, width = 10, height = 10, useDingbats=FALSE)

  #Scatterplot matrix - doubles, uncorrected
  set.seed(1)
  d <- GGally::ggpairs(doubles_dt[apply(doubles_dt[,.SD,,.SDcols = paste0("fitness",all_reps,"_uncorr")]==(-Inf), 1, sum)==0
                  ][sample(x = .N,1000),grep(names(doubles_dt),pattern=paste0("fitness[", all_reps_str, "]_uncorr")),with=F],
          upper=list(continuous = "cor"))
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework_scattermatrix_doubles_uncorr.pdf"), d, width = 10, height = 10, useDingbats=FALSE)

  #Scatterplot matrix - doubles, conditional
  set.seed(1)
  d <- GGally::ggpairs(doubles_dt[apply(doubles_dt[,.SD,,.SDcols = paste0("fitness",all_reps,"_uncorr")]==(-Inf), 1, sum)==0
                  ][sample(x = .N,1000),grep(names(doubles_dt),pattern=paste0("fitness[", all_reps_str, "]_cond")),with=F],
          upper=list(continuous = "cor"))
  ggplot2::ggsave(file.path(outpath, "4_doubles_bayesian_framework_scattermatrix_doubles_cond.pdf"), d, width = 10, height = 10, useDingbats=FALSE)

  ### Merge fitness values
  ###########################

  #### doubles
  # #uncorrected fitness
  # fitness_rx = doubles_dt[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]_uncorr"),colnames(doubles_dt))]
  # # sigma_rx = sqrt(doubles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_uncorr"),colnames(doubles_dt))]^2 + 
  # #                   matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
  # sigma_rx = doubles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_uncorr"),colnames(doubles_dt))]
  # # sigma_s2 = random_effect_model(fitness_rx,sigma_rx)
  # # doubles_dt[,fitness_uncorr := rowSums(fitness_rx/(sigma_rx^2 + sigma_s2),na.rm=T)/rowSums(1/(sigma_rx^2 + sigma_s2),na.rm=T)]
  # doubles_dt[,fitness_uncorr := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  # doubles_dt[,sigma_uncorr := sqrt(1/rowSums(1/(sigma_rx^2+sigma_s2),na.rm=T))]
  # doubles_dt[,sigma_uncorr := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  d <- ggplot2::ggplot(doubles_dt,ggplot2::aes(fitness_uncorr,sigma_uncorr)) + 
    ggplot2::geom_hex() + 
    ggplot2::scale_y_log10() + 
    ggplot2::coord_cartesian(ylim=c(0.01,10))
  ggplot2::ggsave(file.path(outpath, "5_sigma_vs_fitness_doubles_uncorr.pdf"), d, width = 5, height = 5, useDingbats=FALSE)

  #conditioned fitness
  fitness_rx = doubles_dt[,.SD,.SDcols = grep(paste0("fitness[", all_reps_str, "]_cond"),colnames(doubles_dt))]
  # sigma_rx = sqrt(doubles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_cond"),colnames(doubles_dt))]^2 + 
  #                   matrix(replicate_error^2,nrow = dim(fitness_rx)[1],ncol = dim(fitness_rx)[2]))
  sigma_rx = doubles_dt[,.SD,.SDcols = grep(paste0("sigma[", all_reps_str, "]_uncorr"),colnames(doubles_dt))]
  # sigma_s2 = random_effect_model(fitness_rx,sigma_rx)
  # doubles_dt[,fitness_cond := rowSums(fitness_rx/(sigma_rx^2 + sigma_s2),na.rm=T)/rowSums(1/(sigma_rx^2 + sigma_s2),na.rm=T)]
  doubles_dt[,fitness_cond := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  # doubles_dt[,sigma_cond := sqrt(1/rowSums(1/(sigma_rx^2+sigma_s2),na.rm=T))]
  doubles_dt[,sigma_cond := sqrt(1/rowSums(1/(sigma_rx^2),na.rm=T))]
  d <- ggplot2::ggplot(doubles_dt,ggplot2::aes(fitness_cond,sigma_cond)) + 
    ggplot2::geom_hex() + 
    ggplot2::scale_y_log10() + 
    ggplot2::coord_cartesian(ylim=c(0.01,10))
  ggplot2::ggsave(file.path(outpath, "5_sigma_vs_fitness_doubles_cond.pdf"), d, width = 5, height = 5, useDingbats=FALSE)

  #Plot to compare double fitness estimates
  p1=ggplot2::ggplot(doubles_dt,ggplot2::aes(mean_count,fitness_uncorr)) + 
    ggplot2::geom_hex() + 
    ggplot2::scale_x_log10() +
    ggplot2::scale_fill_continuous(trans="log10")
  p2=ggplot2::ggplot(doubles_dt,ggplot2::aes(mean_count,fitness_cond)) + 
    ggplot2::geom_hex()+ 
    ggplot2::scale_x_log10() +
    ggplot2::scale_fill_continuous(trans="log10")
  p3=ggplot2::ggplot(doubles_dt[between(bin_count,2,8)],ggplot2::aes(fitness_uncorr,..scaled..,color=factor(bin_count))) +
    ggplot2::geom_density(adjust=1)
  p4=ggplot2::ggplot(doubles_dt[between(bin_count,2,8)],ggplot2::aes(fitness_cond,..scaled..,color=factor(bin_count))) +
    ggplot2::geom_density(adjust=1)
  p5=ggplot2::ggplot(doubles_dt,ggplot2::aes(fitness_uncorr,sigma_uncorr)) + 
    ggplot2::geom_hex() + 
    ggplot2::scale_y_log10()# +
    # ggplot2::coord_cartesian(ylim = c(0.05,2))
  p6=ggplot2::ggplot(doubles_dt,ggplot2::aes(fitness_cond,sigma_cond)) + 
    ggplot2::geom_hex()+ 
    ggplot2::scale_y_log10()# +
    # ggplot2::coord_cartesian(ylim = c(0.05,2))
  ggplot2::theme_set(ggplot2::theme_minimal())
  #Plot
  d <- cowplot::plot_grid(plotlist = list(p1,p2,p3,p4,p5,p6),nrow=3)
  rm(p1,p2,p3,p4,p5,p6)
  ggplot2::ggsave(file.path(outpath, "5_doubles_fitness_estimates.pdf"), d, width = 10, height = 10, useDingbats=FALSE)

  #Plot fitness values against each other
  set.seed(1)
  d <- GGally::ggpairs(doubles_dt[sample(.N,1000),.(fitness_uncorr,fitness_cond)])
  ggplot2::ggsave(file.path(outpath, "5_doubles_fitness_estimates_scattermatrix.pdf"), d, width = 10, height = 10, useDingbats=FALSE)


  #Plot sigma values against each other
  d <- GGally::ggpairs(doubles_dt[,.(sigma_uncorr,sigma_cond)])
  ggplot2::ggsave(file.path(outpath, "5_doubles_sigma_estimates_scattermatrix.pdf"), d, width = 10, height = 10, useDingbats=FALSE)

  return(doubles_dt)
}
