#' Import clones detected by superFreq for lineplots
#' @param clinicalData data frame with the `Blast` information (as %) for the patients in the cohort.
#' @param patientID a character vector specifying the patient/s id/s for which stories have to be imported.

import_blast_for_lineplot = function(clinicalData, patientID) {
  #import blast fraction
  clinicalData <- clinicalData %>%
    filter(OurPID %in% patientID & !is.na(Time)) #some samples with NA time that messed things up

  if ( nrow(clinicalData) == 0 ) return(NULL)

  #convert to line plot format (whatever that will be)
  samples = unique(clinicalData$SampleName)
  y_matrix = matrix(as.numeric(clinicalData$Blast)/100, ncol=length(samples), dimnames=list('Blast', samples))

  #blast content sometimes missing or not parsed properly. Doesn't have to mean that it's 0, so throw a warning.
  if ( any(is.na(y_matrix)) ) {
    warning('Replacing NA blast content with 0... :o')
    y_matrix[is.na(y_matrix)] = 0
  }

  ret = list(mutations='Blast', y_matrix=y_matrix)
  return(ret)
}
