#' Import clinical information for a patient
#' @param clinicalData data frame containing the clinical information for a patient. This needs to include a column `SampleName` containing unique names for all the RNA-Seq samples available in the cohort; a `PID` column for patient ID and a `Time` column that should reflect the time of sample collection.
#' @param patientID a character vector specifying the patient id for which clinical data should be imported.
#' @param extract_column vector specifying the column to be extracted. Currently only numeric variables with values between 0 and 1 are allowed to be consistent with plotting together with clonality and VAF estimates.
#'
#' @return The function returns a data fram where only `SampleName`, `PID`, `Time` and the variable of interest `extract_column` are kept as well as new variables `VAF`, `mutation_key` and `mutation_det` are added. In this way, clinical information, like tumour purity, is stored in the same way as the variants created with `import_any()` and can be easily combined and plotted together. At the moment, the structure expects `extract_column` to be in the 0-1 range to be consistent with te fact that usually the Variant Allele Frequency (VAF) for a variant is plotted over time.
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
                                     extract_column) {

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

  plot_matrix <- clinicalData %>%
    dplyr::select(.vars = c("SampleName","PID","Time",extract_column))
  colnames(plot_matrix) <- c("SampleName","PID","Time","VAF")

  if(sum(is.na(plot_matrix$Blast)) > 0){
    message(paste0(extract_column," contains NAs. Check before plotting!"))
  }

  plot_matrix <- plot_matrix %>%
    dplyr::mutate(mutation_key = extract_column) %>%
    dplyr::mutate(mutation_det = extract_column) %>%
    dplyr::mutate(variant_type = extract_column)


  return(plot_matrix)

}
