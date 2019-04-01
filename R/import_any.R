#' General import for variants
#'
#' @param variants a data frame where every row is a variant for one sample at a specific time point. The variants can derive from any caller but the input should be standardised to have the following columns: 'SampleName','PID','Time','chrom', 'pos', 'alt', 'ref', 'ref_depth','alt_depth' and (gene) 'SYMBOL' (see more information in Details). The columns `Consequence` and `IMPACT` (as annotated by Variant Effect Predictor (VEP) https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html) are filled with default values if not found. If VEP `Consequence` is not available it could populated with any other informations like exon number, INDEL/SNV label etc...  to add details to each mutations (useful for plotting purposes). The `SampleName` columns is unique for every sequencing sample while `PID` for every patient.
#' @param patientID a character vector specifying the patient/s id/s for which variants have to be imported.
#' @param studyGenes genes of interest.
#' @param minQual minimum quality for a variant to be kept.
#' @param clinicalData clinical data about the patients in the cohort. It has to contain the a column `SampleName`.
#' @param tidy Logical. Should the ouput be in a tidy or untidy (list of matrices) format? Default is `tidy = TRUE`.
#' @param keep_impact vector specifying the IMPACT values to select variants. Values allowed are HIGH, MODERATE, LOW, MODIFIER ( https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html). IMPACT should be a columns of `variants`. If it is not found all variants are kept.
#' @param variant_type Label for the type of variants imported, e.g. vardict-indels.
#'
#' @description   This function will take as input a data frame of variants with specific column information and return a filtered set with sample's clinical infrmation and default variants information also for samples without variants.

#' @details This function will keep only the variants for `patientID` found on `studyGenes` and with a `minQual`. If a sample has no variants, then only clinical information will be returned with default values for the variant information.

#' More details about the `variants` input:
#'
#' - The `Time` column can be defined in any way and it should reflect the time of sample collection. For example it could be defined as Time0, Time1, Time2 etc...
#'
#' - If the column `IMPACT` is not found it will be filled with NAs and no variants will be filtered. Otherwise, values of the columns are checked and if they are within the expected values (HIGH, MODERATE, LOW or MODIFIER) only variants with `keep_impact` entries are kept. If a mutation appears twice with different `IMPACT` values only the most damaging will be kept.
#'
#' The variants are then merged with the clinical information of `patientID`. This step is needed so that if no variants are returned for one time point for one `patientID`, default entries for Variant Allele Frequency (VAF), reference and alterative depths will be created. The default value is 0 for all of the above. A variant is reported for a patient only if at any time point its VAF >= 0.15 and the total depth is >= 10. The function can return the variants in a tidy (long format) or untidy (wide, matrix and list) format. If `tidy = FALSE` the function will return a matrix where the rows are all the unique variants found for `patientID` across time and the columns are the samples of `patientID` across time.

#' @export

#' @examples
#'
#' indels <- data.frame(PID = rep("D1",9),
#'                      Outcome = "Rel",
#'                      SampleName = c("D1.Screen.Rel","D1.Screen.Rel","D1.Screen.Rel",
#'                                     "D1.Cyc1.Rel","D1.Cyc1.Rel","D1.Cyc1.Rel",
#'                                     "D1.Cyc2.Rel","D1.Cyc2.Rel","D1.Cyc2.Rel"),
#'                      Time = c(rep("Screen",3),rep("Cyc1",3),rep("Cyc2",3)),
#'                      chrom = c("chr1","chr1","chr2","chr2","chr1","chr1","chr2","chr1","chr2"),
#'                      pos = c(10,20,30,10,20,30,10,20,30),
#'                      alt = c("ACT","ACGTCG","AGG","ACT","ACGTCG","AGG","ACT","ACGTCG","AGG"),
#'                     ref = c("A","G","A","A","G","AGG","A","G","A"),
#'                      ref_depth = c(12,11,9,9,12,8,13,14,8),
#'                      alt_depth = c(2,20,10,2,10,15,20,20,100),
#'                      SYMBOL = "BCL2",
#'                      Consequence = rep(c("exon1","exon1","exon2"),times=3),
#'                      qual = 49)
#'
#' clinicalData <- data.frame(SampleName = c("D1.Screen.Rel","D1.Cyc1.Rel",
#'                                           "D1.Cyc2.Rel","D1.Cyc3.Rel"),
#'                            PID = "D1",
#'                            AgeDiagnosis = 65,
#'                             Time = c("Screen","Cyc1","Cyc2","Cyc3"),
#'                           Sex = "F",
#'                            BlastPerc = c(80,5,7,40))
#'
#' import_indels <- import_any(variants = indels,
#'                             patientID = "D1",
#'                             studyGenes = "BCL2",
#'                             minQual = 20,
#'                             tidy=TRUE,
#'                             clinicalData = clinicalData)



import_any = function(variants = NULL, patientID = NULL, studyGenes = NULL, minQual=20,
                      clinicalData = NULL, tidy = TRUE,
                      keep_impact = c("HIGH","MODERATE"),
                      variant_type = "indels-vardict") {

  options(warn=-1)

  if( is.null(patientID) ){
    stop("patientID is not defined.")
  }

  if ( is.null(clinicalData)) {
    stop(paste0("No clinicalData available."))
  } else {
    if(nrow(clinicalData) == 0){
      stop(paste0("No lines available in clinicalData."))
    }
  }

  if( is.null(studyGenes) ){
    warning("studyGenes is not defined and all genes will be used.")
    studyGenes <- as.character(unique(variants$SYMBOL))
  }

  # Check variants existence and column requirements
  impact <- TRUE
  search_env <- pryr::where("tidy")

  # Variants existence
  # if(!exists("variants",where = search_env)){
  #  warning("No set of true variants provided.")
  #  return(NULL)
  #}

  # Clinical data existence
  if(!exists("clinicalData",where = search_env)){
    warning("No set of true variants provided.")
    return(NULL)
  }

  if ( nrow(clinicalData) == 0 ) {
    warning(paste0("No lines availble in input for clinicalData"))
    return(NULL)
  }

  ##########################################
  ### Check column requirements for variants
  ##########################################

  facultative_columns <- c("Consequence","IMPACT","qual")
  check_columns <- sum(!(facultative_columns %in% colnames(variants)))

  # Check for facultative columns
  if(check_columns > 0){

    missing <- facultative_columns[!(facultative_columns %in% colnames(variants))]
    #cat(paste0("The following column names are missing in the variants set: ",paste0(missing,collapse=","),".\n"))

    if(sum(missing %in% "Consequence") > 0){
      variants$Consequence <- ""
      cat(paste0("Consequence is missing and will be filled with '.'\n"))
    }

    if(sum(missing %in% "IMPACT") > 0){

      variants$IMPACT <- NA
      impact <- FALSE

      cat(paste0("IMPACT is missing and will be filled with NAs.\n"))

    } else {

      lev <- levels(factor(variants$IMPACT))

      if( sum(!(lev %in% c("HIGH","MODERATE","LOW","MODIFIER"))) > 0) {

        warning(paste0("The impact columns contains unknown values and won't be used for filtering.
                       They should only contain HIGH,MODERATE,LOW,MODIFIER.
                       See https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html.\n"))

        impact <- FALSE

      }

    }

    if(sum(missing %in% "qual") > 0){
      variants$qual <- 0
      minQual <- 0
      cat(paste0("qual is missing and will be filled with 0 and minQual = 0.\n"))
    }

  }

  # necessary columns
  need_columns <- c("SampleName","Time","PID","chrom",
                    "pos","alt","ref","alt_depth","ref_depth","SYMBOL")

  check_columns <- sum(!(need_columns %in% colnames(variants)))

  if(check_columns > 0){
    missing <- need_columns[!(need_columns %in% colnames(variants))]
    stop(paste0("The following columns are missing from variants: ",
                paste0(missing,collapse=","), "\n"))
  }


  # Check clinical data
  need_columns <- c("SampleName","Time","PID")

  check_columns <- sum(!(need_columns %in% colnames(clinicalData)))

  if(check_columns > 0){
    missing <- need_columns[!(need_columns %in% colnames(clinicalData))]
    stop(paste0("The following columns are missing for clinicalData: ",
                paste0(missing,collapse=", ")))
  }


  # Check done

  ##################
  ## Filter variants
  ##################

  # Keep only variants of interest based on patient, impact, quality and genes
  var <- variants %>%
    dplyr::filter(PID %in% patientID) %>% # restrict analysis to OurPID
    tidyr::unite(Location, chrom, pos ,sep="_",remove=FALSE) %>% # create columns that will be useful later: do I need this?
    tidyr::unite(mutation_key, chrom, pos, ref, alt ,sep="_",remove=FALSE) %>% #  - unique IDs for a mutation
    dplyr::filter(SYMBOL %in% studyGenes) %>% # restrict analysis to specific genes
    dplyr::filter(qual >= minQual) %>% # only keep good quality indels
    dplyr::mutate(tot_depth = alt_depth + ref_depth) %>%
    dplyr::mutate(VAF = alt_depth/tot_depth)


  # Missing VEP annotation like with km
  if( sum(is.na(var$IMPACT)) == nrow(var) | !impact) {

    # set all to HIGH
    var <- var %>%
      dplyr::mutate(IMPACT_rank = 1)

  } else {

    # assign numeric ranking
    var <- var %>%
      dplyr::mutate(IMPACT_rank = factor(IMPACT, levels = c("HIGH","MODERATE","LOW","MODIFIER"), labels = c(1,2,3,4))) %>%
      dplyr::mutate(IMPACT_rank = as.numeric(IMPACT_rank))

  }

  var <- var %>%
    dplyr::filter(is.na(IMPACT) | IMPACT %in% keep_impact) %>% # keep NAs and annotated
    dplyr::group_by(SampleName,Location,ref,alt) %>%
    tidyr::unite(mutation_det, SYMBOL, Consequence , sep = " ",remove = FALSE) %>%
    dplyr::filter(order(IMPACT_rank) == order(IMPACT_rank)[which.min(order(IMPACT_rank))]) %>%
    dplyr::filter(!stringr::str_detect(Consequence,c("splice_donor"))) %>%
    dplyr::filter(!stringr::str_detect(Consequence,c("splice_acceptor"))) %>%
    dplyr::mutate(mutation_det = stringr::str_replace(mutation_det,"&.+",""))

  ################################################
  # Get all the clinical information for patientID
  ################################################
  clinicalData <- clinicalData %>%
    dplyr::filter(PID %in% patientID)

  # The following part was added due to the fact that hard coded filters (depht/quality)
  # cased that if a mutation was present at diagnosis but did not passed the filters at remission/relapse
  # it was considered as missing and therefore that sample would have VAF = 0  and ref_depth = 0
  # To improve on this representation we are now keeping all the variants that pass with at least some
  # minimum required filters in any of the time points for patientID. We will then keep
  # all such mutations even for samples that do not meet the filters. This will allow to see if the mutation
  # was actually still present but was partially reduced at some stages.

  # 1. Define all the unique indels found for patientID
  unique_var <- unique(var[,c("chrom","pos","ref","alt","Location","mutation_det","mutation_key","SYMBOL","Consequence")])

  # 2. Create all the possible combinations between variants found and clinical samples
  clinical_var_empty <- merge(clinicalData ,unique_var,all=TRUE)

  # 3. Fill NAs in VAF and ref_depth where a mutation wasn't found in the unfiltered data frames of indels
  clinical_var_fill <- merge(clinical_var_empty,var,all=TRUE) %>%
    dplyr::mutate(ref_depth = ifelse(is.na(ref_depth),0,ref_depth),
                  VAF = ifelse(is.na(VAF),0,VAF),
                  alt_depth = ifelse(is.na(alt_depth),0,alt_depth),
                  tot_depth = ifelse(is.na(tot_depth),0,tot_depth))

  # 4. Filter based on minimal required ref_depth threshold and re-add lost mutations later
  var_keep <- clinical_var_fill %>% dplyr::filter(tot_depth >= 10 & VAF >= 0.15)

  # var_leave contains some indels not found and some that do not meet the filters
  var_leave <- clinical_var_fill %>% dplyr::filter(tot_depth < 10 | VAF < 0.15)


  if(nrow(var_keep) == 0){
    message(paste0("No INDELs found in patient ",patientID))
    return(NULL)
  }

  # 5. If some good quality indels are found then
  # re-add indels whose key (chrom_pos_alt) was present at other time points
  var_saver <-  dplyr::bind_rows(var_keep,
                                 subset(var_leave,mutation_key %in% var_keep$mutation_key))

  # Re-add indels
  var_saver <- var_saver %>%
    dplyr::mutate(variant_type = variant_type)


  ###########################
  # Untidy version of the data
  ###########################

  if( !tidy ){
    var_untidy <- var_saver %>%
      dplyr::select(mutation_det,mutation_key,SYMBOL,Consequence,VAF,SampleName) %>%
      tidyr::spread(key = SampleName, value = VAF, fill = 0)

    ret <- list()
    y_matrix <- var_untidy %>%
      dplyr::select(-mutation_det,-mutation_key,-SYMBOL,-Consequence) %>%
      as.matrix()

    rownames(y_matrix) <- var_untidy$mutation_key

    mutations <- var_untidy$mutation_det

    ret$mutations <- mutations
    ret$y_matrix <- y_matrix

  }

  if (tidy) {
    return(var_saver)
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





