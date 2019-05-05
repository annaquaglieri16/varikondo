#' Import clinical information for a patient
#' @param clinicalData data frame containing the clinical information for a patient. This needs to include a column `SampleName` containing unique names for all the RNA-Seq samples available in the cohort; a `PID` columns for patient ID and a `Time` column that should reflect the time of sample collection.
#' @param patientID a character vector specifying the patient id for which clinical data should be imported.
#' @param extract_column vector specifying the column to be extracted. Currently only numeric variables with values between 0 and 1 are allowed to be consistent with plotting together with clonality and VAF estimates.
#'
#' @examples
#'
#'
#' clinicalData <- data.frame(SampleName = c("D1.Screen.Rel","D1.Cyc1.Rel",
#'                                           "D1.Cyc2.Rel","D1.Cyc3.Rel"),
#'                            PID = "D1",
#'                            Time = c("Screen", "Cyc1","Cyc2","Cyc3"),
#'                            Outcome = "Rel",
#'                            AgeDiagnosis = 65,
#'                            Sex = "F",
#'                            Blast = c(80,5,7,40)/100)
#'
#' import_blast <- import_clinical(clinicalData = clinicalData,
#' extract_column = "Blast",
#' patientID = "D1")
#'
#' @export

#' @import dplyr
#' @import tidyr


import_clinical = function(clinicalData = NA, patientID = NA,
                                     extract_column = "Blast") {

  options(warn=-1)
  if(is.na(clinicalData)){
    stop("clinicalData not availble in input.")
  }

  if( is.na(patientID) ){
    stop("patientID is not defined.")
  }

  # Check clinical data columns
  need_columns <- c("SampleName","Time","PID")
  check_columns <- sum(!(need_columns %in% colnames(clinicalData)))

  if(check_columns > 0){
    missing <- need_columns[!(need_columns %in% colnames(clinicalData))]
    stop(paste0("The following columns are missing for clinicalData: ",
                paste0(missing,collapse=", ")))
  }

  # Modify and filter clinical data
  clinicalData <- clinicalData %>%
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

  # Create tidy format
  tidy_blast <- data.frame(ret$y_matrix) %>%
      dplyr::mutate(mutation_key = rownames(ret$y_matrix)) %>%
      dplyr::mutate(mutation_det = ret$mutations) %>%
      tidyr::gather(SampleName,VAF, 1:ncol(ret$y_matrix)) %>%
      dplyr::mutate(variant_type = extract_column)


  return(tidy_blast)

}
