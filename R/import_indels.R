#' Import INDELs detected by VarDict or Mutect2
#'
#' @param variants a data frame where every row is an INDEL for one sample at a specific time point
#' @param patientID a character vector specifying the patient/s id/s for which INDELs have to be imported.
#' @param studyGenes genes of interest used to subset the `indels` data frame.
#' @param minQual minimum quality for an INDEL to be kept.
#' @param minLength minimum length of the alternative allele required for an INDEL to be kept.
#' @param clinicalData clinical data about the patients in the cohort. It has to contain a column `SampleName` which is made of the following information `OurPID`,`Time`,`Status`,`Repl.within`,`batch` and `Outcome` separated by a '.'.
#' @param tidy Logical. Should the ouput be in a tidy or untidy (list of matrices) format? Default is `tidy = TRUE`.
#'
#' @description   This function will take as input a data frame of variants and return a filtered set of INDELs with sample's clinical infrmation. The variants can derive from any caller but the outout should be standardised to have the following columns: `SampleName`, `PID`, `Time`, `Status`, `Repl.within`, `batch`, `Outcome`, `chrom`, `pos`, `alt`, `ref`, `ref_depth`, `VAF`, `Consequence` and `IMPACT` as annotated by Variant Effect Predictor (VEP).

#' @details This function will keep only the INDELs for `patientID` found on `studyGenes` and with a `minQual` and `minLength` of the alt allele. The function expects the mutations to be annotated with VEP which assign an IMPACT value (HIGH, MODERATE, LOW, MODIFIER) based on https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html. The function only preserves mutations with HIGH or MODERATE impact and for every sample within every patient if a mutation appears twice with different impacts values it will only keep the one with the most damaging impact The INDELs are then merged with the clinical information in `clinicalData` for `patientID`. This means that if no INDELs are returned for one time point within one `patientID` default entries fro Variant Allele Frequency (VAF), reference and alterative depths will be created. The default value is 0 for all of the above. An INDEL is reported for a patient only if at any time point its VAF >= 15 and the total depth is larger than 10. The function can return the INDELs in a tidy (long format) or untidy (wide, matrix and list) format. If `tidy = FALSE` the function will return a matrix where the rows are all the unique INDELs found for `patientID` across time and the columns are the samples of `patientID` across time.

#' @examples
#'
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
#'                      alt_depth = c(2,20,10,2,10,15,20,20,100),
#'                      Consequence = rep(c("C1","C2","C3"),times=3),
#'                      IMPACT = rep(c("MODERATE","HIGH","HIGH"),times=3),
#'                      IMPACT_rank = rep(c(3,2,1),times=3),
#'                      qual = 49,
#'                      SYMBOL = "BCL2") %>%
#'      dplyr::mutate(tot_depth = ref_depth + alt_depth) %>%
#'   dplyr::mutate(VAF = alt_depth/tot_depth) %>%
#'   tidyr::separate(Location,into = c("chrom","pos"),sep = "_") %>%
#'   tidyr::unite(SampleName,PID,Time,Status,Repl.within,batch,Outcome,sep = ".") %>%
#'   tidyr::separate(SampleName,into = c("PID"),sep = "[.]",remove=FALSE)
#'
#' clinicalData <- data.frame(SampleName = c("D1.Screen.Diag.R1.B1.Rel","D1.Cyc1.Rem.R1.B1.Rel",
#'                                      "D1.Cyc2.Rem.R1.B1.Rel","D1.Cyc3.Rel.R1.B1.Rel"),
#'                       AgeDiagnosis = 65,
#'                       Sex = "F",
#'                       BlastPerc = c(80,5,7,40)) %>%
#'     tidyr::separate(SampleName,into=c("PID","Time","Status","Repl.within","batch","Outcome"),sep = "[.]",remove=FALSE)
#'
#'     import_indels <- import_indels_for_lineplot(variants,
#'     patientID = "D1",
#'     studyGenes = "BCL2",
#'     clinicalData = clinicalData)




#' @export

import_indels_for_lineplot = function(variants = variants, patientID, studyGenes, minQual=20, minLength=2,
                                      clinicalData, tidy = TRUE) {

  if ( nrow(variants) == 0 ) {
    warning(paste0("No variants available in input for patient ",patientID))
    return(NULL)
  }

  # No check done for column names.

  # Keep only INDELs of interest based on patient, impact, quality and genes
  indels <- variants %>%
    dplyr::filter(nchar(as.character(alt)) > minLength) %>%
    tidyr::unite(Location, chrom, pos ,sep="_",remove=FALSE) %>% # create columns that will be useful later
    tidyr::unite(mutation_det, Location,ref,alt,sep = "_",remove=FALSE) %>% # create columns that will be useful later
    tidyr::separate(SampleName,into=c("PID","Time","Status","Repl.within","batch","Outcome"),remove = FALSE,sep="[.]") %>%
    dplyr::filter(PID %in% patientID) %>% # restrict analysis to OurPID
    dplyr::filter(SYMBOL %in% studyGenes) %>% # restrict analysis to specific genes
    dplyr::filter(qual > minQual) %>% # only keep good quality indels
    dplyr::filter(IMPACT %in% c("HIGH", 'MODERATE')) %>%
    dplyr::group_by(SampleName,Location,ref,alt) %>% # for every Sample: for every Location called and annotated multiple times only keep the highest rank
    dplyr::filter(qual > minQual) %>%  #quality filter
    dplyr::filter(IMPACT %in% c("HIGH", 'MODERATE')) %>%
    dplyr::group_by(SampleName,Location,ref,alt) %>% # for every Location called only keep the highest rank
    dplyr::filter(order(IMPACT_rank) == order(IMPACT_rank)[which.min(order(IMPACT_rank))])

  ################################################
  # Get all the clinical information for patientID
  ################################################
  clinicalData <- clinicalData %>%
    dplyr::filter(PID %in% patientID & !is.na(Time))
  # some samples with NA time that messed things up, this should be corrected now

  # The following part was added due to the fact that hard coded filters (depht/quality)
  # cased that if a mutation was present at diagnosis but did not passed the filters at remission/relapse
  # it was considered as missing and therefore that sample would have VAF = 0  and ref_depth = 0
  # To improve on this representation we are now keeping all the variants that pass with at least some
  # minimum required filters in any of the time points for patientID. We will then keep
  # all such mutations even for samples that do not meet the filters. This will allow to see if the mutation
  # was actually still present but was partially reduced at some stages.

  # 1. Define all the unique indels found for patientID
  unique_indels <- unique(indels[,c("chrom","pos","ref","Location","mutation_det","alt","SYMBOL","Consequence")])

  # 2. Create all the possible combinations between variants found and clinical samples
  clinical_indels_empty <- merge(clinicalData ,unique_indels,all=TRUE)

  # 3. Fill NAs in VAF and ref_depth where a mutation wasn't found in the unfiltered data frames of indels
  clinical_indels_fill <- merge(clinical_indels_empty,indels,all=TRUE) %>%
    dplyr::mutate(ref_depth = ifelse(is.na(ref_depth),0,ref_depth),
           VAF = ifelse(is.na(VAF),0,VAF),
           alt_depth = ifelse(is.na(alt_depth),0,alt_depth),
           tot_depth = ifelse(is.na(tot_depth),0,tot_depth))

  # 4. Filter based on minimal required ref_depth threshold and re-add lost mutations later
  indels_keep <- clinical_indels_fill %>% filter(tot_depth >= 10 & VAF >= 0.15)

  # indels_leave contains some indels not found and some that do not meet the filters
  indels_leave <- clinical_indels_fill %>% filter(tot_depth < 10 | VAF < 0.15)


  if(nrow(indels_keep) == 0){
    message(paste0("No INDELs found in patient ",patientID))
    return(NULL)
  }

  # 5. If some good quality indels are found then
  # re-add indels whose key (chrom_pos_alt) was present at other time points
  indels_saver <-  dplyr::bind_rows(indels_keep,
                                       subset(indels_leave,mutation_det %in% indels_keep$mutation_det))

  # Re-add indels
  options(warn=-1)
  indels_saver <- indels_saver %>%
      dplyr::mutate(Time = forcats::fct_relevel(Time,"Screen","Cyc1","Cyc2","Cyc3","Cyc4","Cyc9")) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time)))
  options(warn=0)

  ###########################
  # Untidy version of the data
  ###########################

  options(warn=-1)
  indels_untidy <- indels_saver %>%
    dplyr::mutate(Time = forcats::fct_relevel(Time,"Screen","Cyc1","Cyc2","Cyc3","Cyc4","Cyc9")) %>%
    dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time))) %>%
    dplyr::select(mutation_det,SYMBOL,Consequence,VAF,SampleName) %>%
    tidyr::spread(key = SampleName,value=VAF,fill=NA)
  options(warn=0)

  ret <- list()
  ret$y_matrix <- as.matrix(indels_untidy[,4:ncol(indels_untidy)])
  rownames(ret$y_matrix) <- indels_untidy$mutation_det

  ret$mutations <- paste(indels_untidy$SYMBOL,indels_untidy$Consequence)

  if (tidy) {
    return(indels_saver)
  } else {
    return(ret)
  }

  # # Untidy version of the data
  # #make a matrix of y-values for the plot, in this case VAFs
  # #indelsAll$Mutation = paste0(indelsAll$Location, "_", indelsAll$ref, "_", indelsAll$alt)
  # indels$Mutation = paste0(indels$Location, "_", indels$ref, "_", indels$alt)
  # uniqueMutations = !duplicated(indels$Mutation)
  # mutations = indels$Mutation[uniqueMutations]
  # labels = paste0(indels$SYMBOL, ' ', gsub('&.+', '', indels$Consequence))[uniqueMutations]
  # samples = unique(clinicalData$SampleName)
  # vafs =
  #   sapply(samples, function(sample) # for each sample - a patient at a specific time
  #     sapply(mutations, function(mut) { # for each of the unique mutations detected
  #       ret = with(indels, pmin(1, 2*VAF[Mutation == mut & SampleName == sample])) # take the minimum between 1 and the 2*VAF of this mutation for this sample
  #       if ( length(ret) == 0 ) return(0) #return 0 if not found
  #       return(ret)
  #     })
  #   )
  # y_matrix = matrix(vafs, ncol=length(samples), dimnames=list(mutations, samples))
  #
  # # return the y matrix only if the max VAF in one row is > 0.15 and meta data for the mutations
  # ret = list(mutations=labels, y_matrix=y_matrix)
  # ret$mutations =  ret$mutations[matrixStats::rowMaxs(ret$y_matrix) > 0.15]
  # ret$y_matrix =  ret$y_matrix[matrixStats::rowMaxs(ret$y_matrix) > 0.15,,drop=F]

}






