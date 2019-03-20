#' Import clinical information for a specific patient
#' @param clinicalData data frame containing the clinical information `extract_column` and `SampleName` for the patients in the cohort.
#' @param tidy Logical. Should the ouput be in a tidy or untidy (list of matrices) format? Default is `tidy = TRUE`.
#' @param patientID a character vector specifying the patient/s id/s for which stories have to be imported.
#' @param sample_name_parts vector with constitutive parts, separated by a ".", of the `SampleName` column. For example, if the `SampleName` entry for a generic patient looks like: P1.Screening.R1.Responder, `sample_name_parts` will look something like sample_name_parts = c("PID","Time","Replicate","Outcome") and it will be used to create separate columns for the relative information. The components `Time` and `PID` (patient ID) are required. Every sequencing sample needs to have the same `SampleName` structure.
#' @param time_order vector specifying the order of time points, e.g. the levels of the `Time` variable extracted from the `SampleName` column of the input `variants`.
#' @param extract_column vector specifying the column to be extracted. Currently only numeric variables wit values between 0 and 1 are allowed to be consistent with plotting together with clonality and VAF estimates.
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


import_blast_for_lineplot = function(clinicalData = NA, patientID = NA, tidy = TRUE,
                                     extract_column = "Blast",
                                     sample_name_parts = c("PID","Time","Status","Repl.Within","Batch","Outcome"),
                                     time_order = c("Screen","Cyc1","Cyc2","Cyc3","Cyc4","Cyc9")) {

  options(warn=-1)
  if(is.na(clinicalData)){
    stop("clinicalData not availble in input.")
  }
  options(warn=0)

  if( is.na(patientID) ){
    stop("patientID is not defined.")
  }

  if(sum(sample_name_parts %in% "PID") == 0){
    stop("`PID` entry, patient ID, is required in sample_name_parts argument.")
  }

  # Check clinical data columns
  need_columns <- c("SampleName",extract_column)
  check_columns <- sum(!(need_columns %in% colnames(clinicalData)))

  if(check_columns > 0){
    missing <- need_columns[!(need_columns %in% colnames(clinicalData))]
    stop(paste0("The following columns are missing for clinicalData: ",
                paste0(missing,collapse=", ")))
  }

  # Modify and filter clinical data
  clinicalData <- clinicalData %>%
    tidyr::separate(SampleName,into = sample_name_parts,remove = FALSE,sep="[.]") %>%
    dplyr::filter(PID %in% patientID)

  if ( nrow(clinicalData) == 0 ){
    message(paste0("The clinical data has no rows for patient ",patientID))
    return(NULL)
  }

  # The variable to be plotted has to be consistent with VAF
   if(min(clinicalData[,extract_column],na.rm=TRUE) < 0 |
      max(clinicalData[,extract_column],na.rm=TRUE) > 1){
     message(paste0(extract_column," needs to be within 0-1 or 0-100"))
     return(NULL)
   }

  # convert to long format
  samples = unique(clinicalData$SampleName)
  y_matrix = matrix(as.numeric(as.character(clinicalData[,extract_column,drop=TRUE])), ncol=length(samples),
                    dimnames=list(extract_column, samples))

  #blast content sometimes missing or not parsed properly. Doesn't have to mean that it's 0, so throw a warning.
  if ( any(is.na(y_matrix)) ) {
    warning(paste0('Replacing NA ', extract_column,' content with 0... :o'))
    y_matrix[is.na(y_matrix)] = 0
  }

  ret = list(mutations=extract_column, y_matrix=y_matrix)

  options(warn=-1)

  if(tidy){

    tidy_blast <- data.frame(ret$y_matrix) %>%
      dplyr::mutate(mutation_key = rownames(ret$y_matrix)) %>%
      dplyr::mutate(mutation_det = ret$mutations) %>%
      tidyr::gather(SampleName,VAF, 1:ncol(ret$y_matrix)) %>%
      tidyr::separate(SampleName, into = sample_name_parts,sep="[.]",remove=FALSE) %>%
      dplyr::mutate(Time = factor(Time,levels = time_order)) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time))) %>%
      dplyr::mutate(variant_type = extract_column)
  }
  options(warn=0)

  if(tidy){
    return(tidy_blast)
  }else{
    return(ret)
  }

}
