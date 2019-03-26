#' Import point mutations for lineplots
#' @param mutationSummary a matrix summarising mutations reported by superFreq. Every row is a muations, every column is a sample and every entry is the VAF for that row mutation for that column sample.
#' @param patientID a character vector specifying the patient/s id/s for which mutations have to be imported. This value has to correspond to the names used to save the `mutationSummary` object.
#' @param tidy Logical. Should the ouput be in a tidy or untidy (list of matrices) format? Default is `tidy = TRUE`.
#' @param sample_name_parts vector with constitutive parts of the names for every sample. For example, if the sample names are pf the format: patientID.Time.Status.Repl.within.batch.Outcome, the vector will be sample_name_parts = c("patientID","Time","Status","Repl.Within","Batch","Outcome"). The component `Time` is required.
#' @param time_order vector specifying the order of time points, e.g. the levels of the `Time` variable `sample_name_parts` definition. This will be used to order samples for plotting purposes.
#'
#' @examples
#'
#' mutationSummary <- list()
#' mutationSummary[[1]] <- data.frame(x1 = c(1000,2000),
#' x2 = c(1000,2001),
#' call = c("1001+TCCG","chr22.A"),
#' stories.D1.Screen.Diag.R1.B1.Rel = c(0.33,0.74),
#' stories.D1.Cyc1.Rem.R1.B1.Rel = c(0,0.68),
#' stories.D1.Cyc2.Rem.R1.B1.Rel = c(0,0.1),
#' errors.D1.Screen.Diag.R1.B1.Rel = c(0,0.5),
#' errors.D1.Cyc1.Rem.R1.B1.Rel = c(0.04, 0.2),
#' errors.D1.Cyc2.Rem.R1.B1.Rel = c(0.2,0.1),
#' GoI = c("FLT3","STAG2"),
#' label = c("FLT3 (13) frameshift","1Mb A"),
#' severity = c(5,0),
#' chr = c(13,22),
#' start = c(2029,4554),
#' end = c(2035,4600))
#' rownames(mutationSummary[[1]]) <- mutationSummary[[1]]$call
#' names(mutationSummary) <- "D1"
#'


import_snvs = function(mutationSummary = NA, patientID = NA, tidy = TRUE,
                                    time_order = c("Screen","Cyc1","Cyc2","Cyc3","Cyc4","Cyc9"),
                                    sample_name_parts = c("PID","Time","Status","Repl.Within","Batch","Outcome")) {

  options(warn=-1)
  if(is.na(mutationSummary)){
    stop("mutationSummary not availble in input.")
  }
  options(warn=0)

  if( is.na(patientID) ){
    stop("patientID is not defined.")
  }

  if(sum(sample_name_parts %in% "Time") == 0){
    stop("`Time` is required in sample name components specified in sample_name_parts.")
  }

  #import data and subset on patient
  mutations = mutationSummary[[patientID]]

  #only point mutations
  mutations = mutations[mutations$x1 == mutations$x2,]

  #return null if no mutations
  if ( nrow(mutations) == 0 ) {
    warning(paste0("No SNVs available in input for patient ",patientID))
    return(NULL)
  }

  #remove chromosome from the label
  mutations$label = gsub('( ?)\\(.*\\)', '', mutations$label)

  #convert to line plot format. (whatever that will be)
  y_matrix = mutations$stories

  ret = list(mutations=mutations$label, y_matrix=y_matrix)
  ret$mutations =  ret$mutations[matrixStats::rowMaxs(ret$y_matrix) > 0.15]
  ret$y_matrix =  ret$y_matrix[matrixStats::rowMaxs(ret$y_matrix) > 0.15,,drop=F]

  if ( length(ret$mutations) == 0 ){
    message(paste0("No SNVs found in patient ",patientID))
    return(NULL)
  }

  options(warn=-1)
  if(tidy){

    tidy_snvs <- data.frame(ret$y_matrix) %>%
      dplyr::mutate(mutation_key = rownames(ret$y_matrix)) %>%
      dplyr::mutate(mutation_det = ret$mutations) %>%
      tidyr::gather(SampleName,VAF, 1:ncol(ret$y_matrix)) %>%
      tidyr::separate(SampleName, into = sample_name_parts,sep="[.]",remove=FALSE) %>%
      dplyr::mutate(Time = factor(Time,levels = time_order)) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time))) %>%
      tidyr::separate(mutation_det, into = c("SYMBOL"),sep=" ",remove=FALSE) %>%
      dplyr::mutate(variant_type = "superfreq-SNVs")


  }
  options(warn=0)

  if(tidy){
    return(tidy_snvs)
  }else{
    return(ret)
  }
}
