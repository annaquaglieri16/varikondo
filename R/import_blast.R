#' Import blast contect for a specific patient
#' @param clinicalData data frame with the percentage of `Blast` information for the patients in the cohort.
#' @param tidy Logical. Should the ouput be in a tidy or untidy (list of matrices) format? Default is `tidy = TRUE`.
#' @param patientID a character vector specifying the patient/s id/s for which stories have to be imported.
#'
#' @examples
#'
#' patientID = "D1"
#'
#' clinicalData <- data.frame(SampleName = c("D1.Screen.Diag.R1.B1.Rel","D1.Cyc1.Rem.R1.B1.Rel",
#'                                      "D1.Cyc2.Rem.R1.B1.Rel","D1.Cyc3.Rel.R1.B1.Rel"),
#'                       AgeDiagnosis = 65,
#'                       Sex = "F",
#' BlastPerc = c(80,5,7,40)) %>%
#'     tidyr::separate(SampleName,into=c("PID","Time","Status","Repl.within","batch","Outcome"),sep = "[.]",remove=FALSE)
#'
#'     import_indels <- import_indels_for_lineplot(variants,
#'     patientID = "D1",
#'     studyGenes = "BCL2",
#'     clinicalData = clinicalData)


import_blast_for_lineplot = function(clinicalData, patientID, tidy = TRUE) {

  #import blast fraction
  clinicalData <- clinicalData %>%
    filter(PID %in% patientID & !is.na(Time)) #some samples with NA time that messed things up

  if ( nrow(clinicalData) == 0 ){
    message(paste0("The clinical data has no rows for patient ",patientID))
    return(NULL)
  }


  #convert to line plot format (whatever that will be)
  samples = unique(clinicalData$SampleName)
  y_matrix = matrix(as.numeric(as.character(clinicalData$Blast))/100, ncol=length(samples), dimnames=list('Blast', samples))

  #blast content sometimes missing or not parsed properly. Doesn't have to mean that it's 0, so throw a warning.
  if ( any(is.na(y_matrix)) ) {
    warning('Replacing NA blast content with 0... :o')
    y_matrix[is.na(y_matrix)] = 0
  }

  ret = list(mutations='Blast', y_matrix=y_matrix)

  options(warn=-1)

  if(tidy){

    tidy_blast <- data.frame(ret$y_matrix) %>%
      dplyr::mutate(mutation_key = rownames(ret$y_matrix)) %>%
      dplyr::mutate(mutation_det = ret$mutations) %>%
      tidyr::gather(SampleName,VAF, 1:ncol(ret$y_matrix)) %>%
      tidyr::separate(SampleName, into = c("PID","Time","Status","Repl.Within","Batch","Outcome"),sep="[.]",remove=FALSE) %>%
      dplyr::mutate(Time = forcats::fct_relevel(Time,"Screen","Cyc1","Cyc2","Cyc3","Cyc4","Cyc9")) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time))) %>%
      dplyr::mutate(variant_type = "Blasts")
  }
  options(warn=0)

  if(tidy){
    return(tidy_blast)
  }else{
    return(ret)
  }

}
