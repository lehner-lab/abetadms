
#' abetadms__datatable_to_matrix
#'
#' Convert data.table to matrix for heatmap plotting (single AA mutants only).
#'
#' @param input_dt input data.table (required)
#' @param variable_name name of variable to use for heatmap cells (defaut:"fitness")
#'
#' @return a matrix for heamtmap plotting
#' @export
abetadms__datatable_to_matrix <- function(
  input_dt,
  variable_name="fitness" 
  ){
  aa_obj <- Biostrings::AAString("GAVLMIFYWKRHDESTCNQP")
  aa_list <- Biostrings::AMINO_ACID_CODE[strsplit(as.character(aa_obj), NULL)[[1]]]
  aa_list["*"] <- "X"

  #Only single AA mutants
  input_df <- as.data.frame(input_dt)
  dms_1aa <- input_df[input_df$Nmut_aa==1,]
  #Absolute position
  # dms_1aa$Pos_abs <- as.numeric(substr(dms_1aa$mut_code, 2, 4))
  #WT sequence
  wt_seq <- unique(dms_1aa[order(dms_1aa[,"Pos_abs"]),c("WT_AA", "Pos_abs")])[,"WT_AA"]

  #Construct heatmap matrix
  heat_mat <- matrix(nrow = length(aa_list), ncol = max(dms_1aa$Pos_abs)-min(dms_1aa$Pos_abs)+1)
  rownames(heat_mat) <- names(aa_list)
  colnames(heat_mat) <- min(dms_1aa$Pos_abs):max(dms_1aa$Pos_abs)
  for(aa_pos in min(dms_1aa$Pos_abs):max(dms_1aa$Pos_abs)){
    for(aa_id in names(aa_list)){
      temp_index <- which(dms_1aa$Pos_abs==aa_pos & dms_1aa$Mut==aa_id)
      if(length(temp_index)==1){
        heat_mat[aa_id,as.character(aa_pos)] <- dms_1aa[temp_index,variable_name]
      }
    }
  }
  return(heat_mat)
}

