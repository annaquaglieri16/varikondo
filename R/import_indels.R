#' Import INDELs detected by VarDict
#' @param indels a data frame where every row is an INDEL for one sample at a specific time point
#' @param patientID a character vector specifying the patient/s id/s for which INDELs have to be imported.
#' @param studyGenes genes of interest used to subset the `indels` data frame.
#' @param minQual minimum quality for an INDEL to be kept.
#' @param minLength minimum length of the alternative allele required for an INDEL to be kept.
#' @param clinicalData clinical data about the patients in the cohort. It has to contain a column `SampleName` which is made of the following information `OurPID`,`Time`,`Status`,`Repl.within`,`batch` and `Outcome` separated by a '.'.
#' #' @param tidy Logical. Should the ouput be in a tidy or untidy (list of matrices) format?
#'
#' @description   This function will take as input a data frame of INDELs `indels`. Out of it it will keep only the INDELs for `patientID` found on `studyGenes` and with a `minQual` and `minLength` of the alt allele. The function expects the mutations to be annotated with VEP which assign an IMPACT value (HIGH, MODERATE, LOW, MODIFIER) based on https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html. `import_indels_for_lineplot` only preserves mutations with HIGH or MODERATE impact and for every sample within every patient if a mutation appears twice with different impatc values it will only keep the one with the most damaging IMPACT. The INDELs are then merged with the clinical information in `clinicalData` for `patientID`. This means that if no INDELs are returned for one time point within one `patientID` default entries fro Variant Allele Frequency (VAF), reference and alterative depths will be created. The default value is 0 for all of the above. It will return a matrix where the rows are all the unique INDELs found for `patientID` across time and the columns are the samples of `patientID` across time. Every entry of this matrix is the VAF for each sample for that row mutation. An INDEL is reported in output only if the max VAF across time is > 0.15

#' @export

import_indels_for_lineplot = function(indels = indels, patientID, studyGenes, minQual=20, minLength=2,
                                      clinicalData, tidy = TRUE) {

  if ( nrow(indels) == 0 ) {
    warning(paste0("No variants available in input for patient ",patientID))
    return(NULL)
  }

  # No check done for column names.

  indels <- indels %>%
    dplyr::filter(nchar(as.character(alt)) > minLength) %>%
    tidyr::separate(SampleName,into=c("PID","Time","Status","Repl.within","batch","Outcome"),remove = FALSE,sep="[.]") %>%
    dplyr::filter(PID %in% patientID) %>% # restrict analysis to OurPID
    dplyr::filter(SYMBOL %in% studyGenes) %>% # restrict analysis to specific genes
    dplyr::filter(qual > minQual) %>% # only keep good quality indels
    dplyr::filter(IMPACT %in% c("HIGH", 'MODERATE')) %>%
    dplyr::group_by(SampleName,Location,ref,alt) %>% # for every Sample: for every Location called and annotated multiple times only keep the highest rank
    dplyr::filter(order(IMPACT_rank) == order(IMPACT_rank)[which.min(order(IMPACT_rank))])

  # needed to get all the sample names
  clinicalData <- clinicalData %>%
    dplyr::filter(PID %in% patientID & !is.na(Time)) # some samples with NA time that messed things up, this should be corrected now

  # we will add the clinical information for samples without mutation
  indelsAll <- indels[!duplicated(indels[,c("Location","ref","alt")]),]

  # indelsAll: inlcudes all unique indels called for this PID across time
  # indels: include all indels - the same one repeated twice if called at different times for the same patient
  # mutations: includes unique mutation name (chr_pos_alt) for a patient
  # labels: has the same length as mutations and contains the labels (Symbol + consequence) for the unique mutations

  # Untidy version of the data
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

  # return the y matrix only if the max VAF in one row is > 0.15 and meta data for the mutations
  ret = list(mutations=labels, y_matrix=y_matrix)
  ret$mutations =  ret$mutations[matrixStats::rowMaxs(ret$y_matrix) > 0.15]
  ret$y_matrix =  ret$y_matrix[matrixStats::rowMaxs(ret$y_matrix) > 0.15,,drop=F]

  if ( length(ret$mutations) == 0 ) {
    message(paste0("No indels found in patient ",patientID))
    return(NULL)
  }

  # Tidy version of the data
  if(tidy){

    tidy_indels <- data.frame(ret$y_matrix) %>%
      dplyr::mutate(mutation_key = rownames(ret$y_matrix)) %>%
      dplyr::mutate(mutation_det = ret$mutations) %>%
      tidyr::gather(SampleName, VAF,1:ncol(ret$y_matrix)) %>%
      tidyr::separate(SampleName, into = c("PID","Time","Status","Repl.Within","Batch","Outcome"),sep="[.]",remove=FALSE) %>%
      dplyr::mutate(Time = forcats::fct_relevel(Time,"Screen","Cyc1","Cyc2","Cyc3","Cyc4","Cyc9")) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time)))
  }

  if (tidy) {
    return(tidy_indels)
  } else {
    return(ret)
  }

}


#' @examples
#'
#` studyGenes
#` studyGenes <- "BCL2"

#' minLength = 2
#' patientID = "D1"
#' minQual=20
#
#' indels <- data.frame(PID = rep("D1",9),
#'                      Repl.within = "R1",
#'                      batch = "B1",
#'                      Status = c(rep("Diag",3),rep("Rem",3),rep("Rel",3)),
#'                      Outcome = "Rel",
#'                      Time = c(rep("Screen",3),rep("Cyc1",3),rep("Cyc2",3)),
#'                      Location = c("chr1_10","chr1_20","chr2_30","chr1_10","chr1_20","chr2_30","chr1_10","chr1_20","chr2_30"),
#'                      alt = c("ACT","ACGTCG","AGG","ACT","ACGTCG","AGG","ACT","ACGTCG","AGG"),
#'                      ref = c("A","G","A","A","G","AGG","A","G","A"),
#'                      ref_depth = c(12,11,9,9,12,8,13,14,8),
#'                      VAF = rep(0.7,9),
#'                      Consequence = rep(c("C1","C2","C3"),times=3),
#'                      IMPACT = rep(c("MODERATE","HIGH","HIGH"),times=3),
#'                      IMPACT_rank = rep(c(3,2,1),times=3),
#'                      qual = 49,
#'                      SYMBOL = "BCL2") %>%
#'   dplyr::mutate(variant_ID = c("A","B","C","A","B","C","A","B","C"),
#'                 tot_depth = ref_depth/(1-VAF)) %>%
#'   dplyr::mutate(alt_depth = tot_depth * VAF) %>%
#'   tidyr::unite(key, Location,alt,sep = "_",remove=FALSE) %>%
#'   tidyr::separate(Location,into = c("chrom","pos"),sep = "_",remove=FALSE) %>%
#'   tidyr::unite(SampleName,PID,Time,Status,Repl.within,batch,Outcome,sep = ".") %>%
#'   tidyr::separate(SampleName,into = c("PID"),sep = "[.]",remove=FALSE)
#'
#' # example clinical
#' clinical <- data.frame(PID = rep("D1",4),
#'                        Time = c("Screen","Cyc1","Cyc2","Cyc3"))
#'




