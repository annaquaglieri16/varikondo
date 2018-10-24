#' Import INDELs detected by VarDict
#' @param indels a data frame where every row is an INDEL for one sample at a specific time point
#' @param patientID a character vector specifying the patient/s id/s for which INDELs have to be imported.
#' @param studyGenes genes of interest used to subset the `indels` data frame.
#' @param minQual minimum quality for an INDEL to be kept.
#' @param minLength minimum length of the alternative allele required for an INDEL to be kept.
#' @param clinicalData clinical data about the patients in the cohort. It has to contain a column `SampleName` which is made of the following information `OurPID`,`Time`,`Status`,`Repl.within`,`batch` and `Outcome` separated by a '.'.


#' @details  This function will take as input a data frame of INDELs `indels`. It will keep only the INDELs for `patientID` on `studyGenes` and with a `minQual` and `minLength` of the alt allele. It will return a matrix where the rows are all the unique INDELs found for `patientID` across time and the columns are the samples of `patientID` across time. Every entry of this matrix is the VAF for each sample for that row mutation. An INDEL is reported in output only if the max VAF across time is > 0.15


import_indels_for_lineplot = function(indels = indels, patientID, studyGenes, minQual=20, minLength=2,
                                      clinicalData) {

  if ( nrow(indels) == 0 ) {
    warning(paste0("No variants available in input for patient ",patientID))
    return(NULL)
  }

  # No check done for column names.

  indels <- indels %>%
    dplyr::filter(nchar(alt) > minLength) %>%
    tidyr::separate(SampleName,into=c("OurPID","Time","Status","Repl.within","batch","Outcome"),remove = FALSE,sep="[.]") %>%
    dplyr::filter(OurPID %in% patientID) %>% # restrict analysis to OurPID
    dplyr::filter(SYMBOL %in% studyGenes) %>% # restrict analysis to specific genes
    dplyr::filter(qual > minQual) %>%
    dplyr::filter(IMPACT %in% c("HIGH", 'MODERATE')) %>%
    dplyr::group_by(Location,ref,alt) %>% # for every Location called only keep the highest rank
    dplyr::filter(order(IMPACT_rank) == order(IMPACT_rank)[which.min(order(IMPACT_rank))])

  #return empty if no variants passed the filters
  if ( nrow(indels) == 0 ) return(NULL)

  # needed to get all the sample names
  clinicalData <- clinicalData %>%
    dplyr::filter(OurPID %in% patientID & !is.na(Time)) # some samples with NA time that messed things up, this should be corrected now

  # we will add the clinical information for samples without mutation
  indelsAll <- indels[!duplicated(indels[,c("Location","ref","alt")]),]

  # indelsAll: inlcudes all unique indels called for this PID across time
  # indels: include all indels - the same one repeated twice if called at different times for the same patient
  # mutations: includes unique mutation name (chr_pos_alt) for a patient
  # labels: has the same length as mutations and contains the labels (Symbol + consequence) for the unique mutations

  #make a matrix of y-values for the plot, in this case VAFs
  indelsAll$Mutation = paste0(indelsAll$Location, "_", indelsAll$ref, "_", indelsAll$alt)
  indels$Mutation = paste0(indels$Location, "_", indels$ref, "_", indelsAll$alt)
  mutations = indelsAll$Mutation[!duplicated(indelsAll$Mutation)]
  labels = paste0(indelsAll$SYMBOL, ' ', gsub('&.+', '', indelsAll$Consequence))[!duplicated(indelsAll$Mutation)]
  samples = unique(clinicalData$SampleName)
  vafs =
    sapply(samples, function(sample) # for each sample - a patient at a specific time
      sapply(mutations, function(mut) { # for each of the unique mutations detected
        ret = with(indels, pmin(1, 2*VAF[Mutation == mut & SampleName == sample])) # take the minimum between 1 and the 2*VAF of this mutation for this sample
        if ( length(ret) == 0 ) return(0) #return 0 if not found
        return(ret)
      })
    )
  y_matrix = matrix(vafs, ncol=length(samples), dimnames=list(mutations, samples))

  #return the y matrix only if the max VAF in one row is > 0.15 and meta data for the mutations
  ret = list(mutations=labels, y_matrix=y_matrix)
  ret$mutations =  ret$mutations[matrixStats::rowMaxs(ret$y_matrix) > 0.15]
  ret$y_matrix =  ret$y_matrix[matrixStats::rowMaxs(ret$y_matrix) > 0.15,,drop=F]
  if ( length(ret$mutations) == 0 ) return(NULL)
  return(ret)

}


