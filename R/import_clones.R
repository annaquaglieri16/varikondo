#' Import clones detected by superFreq for a specific patient
#' @param file path to where the superFreq stories.Rdata for `patientID` is saved.
#' @param patientID a character vector specifying the patient/s id/s for which stories have to be imported.
#' @param tidy Logical. Should the ouput be in a tidy or untidy (list of matrices) format? Default is `tidy = TRUE`.



import_clones_for_lineplot = function(file, patientID, tidy = TRUE) {
  #import data and subset on patient
  #file = file.path(stories_dir, patientID, 'stories.Rdata')

  if ( !file.exists(file) ){
    message(paste0(file," does not exist for patient ",patientID))
    return(NULL)
  }

  load(file)
  storyMx = stories$stories[[1]]$consistentClusters$cloneStories$stories
  anchors = c(rownames(stories$stories[[1]]$anchorStories$anchorSNVs),
              rownames(stories$stories[[1]]$anchorStories$anchorCNAs))
  Nmut = sapply(stories$stories[[1]]$consistentClusters$storyList, function(muts) sum(muts %in% anchors))

  use = names(Nmut) != 'germline'

  ret = list(mutations=paste0('clone (', Nmut[use], ' anchors)'), y_matrix=storyMx[use,,drop=F])
  ret$mutations =  ret$mutations[matrixStats::rowMaxs(ret$y_matrix) > 0.15]
  ret$y_matrix =  ret$y_matrix[matrixStats::rowMaxs(ret$y_matrix) > 0.15,,drop=F]


  if ( length(ret$mutations) == 0 ) {
    message(paste0("No clones found in patient ",patientID))
    return(NULL)
  }

  options(warn=-1)
  if(tidy){

    tidy_clones <- data.frame(ret$y_matrix) %>%
      dplyr::mutate(mutation_key = rownames(ret$y_matrix)) %>%
      dplyr::mutate(mutation_det = ret$mutations) %>%
      tidyr::gather(SampleName,VAF, 1:ncol(ret$y_matrix)) %>%
      tidyr::separate(SampleName, into = c("PID","Time","Status","Repl.Within","Batch","Outcome"),sep="[.]",remove=FALSE) %>%
      dplyr::mutate(Time = forcats::fct_relevel(Time,"Screen","Cyc1","Cyc2","Cyc3","Cyc4","Cyc9")) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time))) %>%
      dplyr::mutate(variant_type = "superfreq")
  }
  options(warn=0)

  if(tidy){
    return(tidy_clones)
  }else{
    return(ret)
  }

  }
