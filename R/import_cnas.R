#' Import copy number variation detected by superFreq for a specific patient
#' @param mutationSummary a matrix summarising mutations reported by superFreq. Every row is a CN event, every column is a sample and every entry is the VAF for that row mutation for that column sample.
#' @param tidy Logical. Should the ouput be in a tidy or untidy (list of matrices) format? Default is `tidy = TRUE`.
#' @param patientID a character vector specifying the patient/s id/s for which mutations have to be imported. This value has to correspond to the names used to save the `mutationSummary` object.



import_cnas_for_lineplot = function(mutationSummary, patientID, tidy = TRUE) {
  #import data and subset on patient
  mutations = mutationSummary[[patientID]]

  #only CNAs mutations, not from sex chromosomes
  mutations = mutations[mutations$x1 != mutations$x2 & !(superFreq::xToChr(mutations$x1, genome='hg38') %in% c('X', 'Y')),]

  #return null if no CNAs
  if ( nrow(mutations) == 0 ){
    warning(paste0("No CNA available in input for patient ",patientID))
    return(NULL)
  }

  #remove chromosome from the label, and add affected gene
  mutations$label = gsub('( ?)\\(.*\\)', '', mutations$label)
  #mutations$label = gsub('[0-9]*.bp ', '', mutations$label)
  mutations$label = paste0(mutations$GoI, ' ', mutations$label)

  #convert to line plot format. (whatever that will be)
  y_matrix = mutations$stories

  ret = list(mutations=mutations$label, y_matrix=y_matrix)
  ret$mutations =  ret$mutations[matrixStats::rowMaxs(ret$y_matrix) > 0.15]
  ret$y_matrix =  ret$y_matrix[matrixStats::rowMaxs(ret$y_matrix) > 0.15,,drop=F]

  if ( length(ret$mutations) == 0 ) {
    message(paste0("No CNAs found in patient ",patientID))
    return(NULL)
  }

  options(warn=-1)
  if(tidy){

    tidy_cnas <- data.frame(ret$y_matrix) %>%
      dplyr::mutate(mutation_key = rownames(ret$y_matrix)) %>%
      dplyr::mutate(mutation_det = ret$mutations) %>%
      tidyr::gather(SampleName,VAF, 1:ncol(ret$y_matrix)) %>%
      tidyr::separate(SampleName, into = c("PID","Time","Status","Repl.Within","Batch","Outcome"),sep="[.]",remove=FALSE) %>%
      dplyr::mutate(Time = forcats::fct_relevel(Time,"Screen","Cyc1","Cyc2","Cyc3","Cyc4","Cyc9")) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time))) %>%
      dplyr::mutate(variant_type = "CNAs")

  }
  options(warn=0)

  if(tidy){
    return(tidy_snvs)
  }else{
    return(ret)
  }

}

