#' Extract SNVs, CNVs and clones on genes of interest from superFreq output
#'
#' @param superFreq_R_path Path to superFreq R folder.
#' @param superFreq_meta_path Path to superFreq run cohort metadata (file ending with .tsv). If the sample names in the NAME field of the metadata file start with a number, you need to add an "X" before the name with paste0() since this is done by default by superFReq.
#' @param studyGenes character vector containing genes of interest.
#' @param patientID a character vector specifying the patient/s id/s for which variants have to be imported.
#' @param ref_genome character vector for the reference genome used in the analysis ('hg38' or 'hg19')
#' @param VAFcut Minimum VAF for a variant found in a patient.
#' @description This function imports superFreq's SNVs, CNAs and clones for one patient and outputs them into a tidy format where every row is a variant.
#' @export

# studyGenes <- read_csv("../../../venetoclax_trial/Recurrent-AML-genes-across-studies.csv")
# studyGenes <- as.character(studyGenes$Symbol)
# superFreq_R_path <- "../../../venetoclax_trial/superFreq/R"
# patientID <- "D1"
# VAFcut <- 0.15
# superFreq_meta_path <- file.path("../../../venetoclax_trial/superFreq/metaData.tsv")
# ref_genome = "hg38"
#
# genes0 <- read.csv("../../../venetoclax_trial/Recurrent-AML-genes-across-studies.csv")
# studyGenes <- as.character(genes0$Symbol[genes0$CBF_AML])
# superFreq_R_path = "../../../cbf_aml_agrf/superFreq/R_full_cohort/"
# superFreq_meta_path = "../../../cbf_aml_agrf/superFreq/runFullCohort/metadata_varscan.tsv"
# meta <- read_delim("../../../cbf_aml_agrf/superFreq/runFullCohort/metadata_varscan.tsv",delim = "\t")
# patients <- unique(meta$INDIVIDUAL)
# time_order = c("Dia","Rem","Rel")
# table(meta$TIMEPOINT)
# patientID <- "CAN02JAB"

import_goi_supefreq <- function(superFreq_R_path = superFreq_R_path,
                                 superFreq_meta_path = superFreq_meta_path,
                                 studyGenes,
                                 patientID = "D1",
                                 ref_genome = "hg38",
                                 VAFcut = 0.15){

  options(warn = -1)
  if( is.na(patientID) ){
    stop("patientID is not defined.")
  }

  storyFiles = file.path(superFreq_R_path, patientID, 'stories.Rdata')
  clusterFiles = file.path(superFreq_R_path, patientID, 'clusters.Rdata')

  cat(patientID, '...', sep='')

  # Check existence
  if ( !file.exists(superFreq_meta_path) ){
    warning(paste0("'superFreq_R_meta' does not exist"))
    return(data.frame())
  }

  if ( !file.exists(storyFiles) ){
    warning(paste0("No 'stories.Rdata' available in ",superFreq_R_path))
    return(data.frame())
  }

  if ( !file.exists(clusterFiles) ){
    warning(paste0("No 'clusterFiles' available in ",superFreq_R_path))
    return(data.frame())
  }

  # Load files if they exist
  load(storyFiles) # should load a `stories` object
  load(clusterFiles)
  metadata <- read_delim(superFreq_meta_path,delim = "\t",progress = FALSE)

  # This list contains one data frame per sample - contains all the variants/somatic/germline/noise called in the VCF
  qs = stories$variants$variants
  # This data frame contains both somatic mutation (SNVs,CN) tracked
  allS = stories$stories[[1]]$all  # q for chris - do I need this [[1]] ?

  #
  #all_long <- allS %>% gather(key = SampleName,value = VAF,colnames(allS$stories))
  #a = do.call(rbind, lapply(colnames(allS$stories), function(sample) data.frame(allS[,1:3], allS$stories[,sample], 'SampleName'=sample)))

  # samples available for this patientID
  samples = names(qs)

  # Only somatic SNvs are annotated
  # Keep only unique names SNVs that satisfies certain filters across samples of this patient - add filter on VAF here
  SNVoInames = unique(
    unlist(
      lapply(
        qs, function(q) {

          q <- q %>%
            dplyr::mutate(VAF = var/cov,
                          mut_name = rownames(q)) %>%
            dplyr::filter((somaticP > 0.5) &
                            (!is.na(q$severity) & severity < 12) &
                            inGene %in% studyGenes &
                            VAF >= VAFcut)
          as.character(q$mut_name)

        }
      )
    )
  )

  SNVoIannotation = do.call(rbind, lapply(SNVoInames, function(snv) {
    qs[[which(sapply(qs, function(q) q[snv,]$somaticP > 0))[1]]][snv,]
  }))

  # DEfine what are the SNVs of interest - at this stage the dataframe includes only VAF
  if(length(SNVoInames) == 0) {
    message(paste0("No SNVs available for patient ",patientID))
    SNVoI <- NULL
  } else{
    SNVoI = allS[SNVoInames,]
    SNVoI$SYMBOL = SNVoIannotation$inGene
    SNVoI$variant_type <- "SNV"
  }

  # Get CNVs for genes of interest
  GoIx = sapply(studyGenes, function(gene) mean(as.numeric(clusters[[1]]$CR[gene,c('x1', 'x2')])))
  CNAoI = do.call(rbind, lapply(studyGenes, function(gene) {
    ret = allS[which(allS$x1 < GoIx[gene] & allS$x2 > GoIx[gene]),]
    ret$SYMBOL = rep(gene, nrow(ret))
    return(ret)
  }))

  if(nrow(CNAoI) == 0) {
    message(paste0("No CNAs available for patient ",patientID))
    CNAoI <- NULL
  } else{
    CNAoI$variant_type <- "CNA"
  }


  if(!is.null(SNVoI) | !is.null(CNAoI)){

    MoI = rbind(SNVoI, CNAoI)
    labels = superFreq:::storyToLabel(MoI, stories$variants, genome=ref_genome, mergeCNAs=F)[as.character(seq(along.with=MoI$x1)),]
    MoI$label = labels$label
    MoI$severity = labels$severity
    MoI$chr = superFreq:::xToChr(MoI$x1)
    MoI$start = superFreq:::xToPos(MoI$x1)
    MoI$end = superFreq:::xToPos(MoI$x2)

    # Convert to long format SNVs
    MoI_stories <- data.frame(MoI$stories)
    MoI <- data.frame(pos = MoI$start,end = MoI$end,
                      chrom = paste0("chr",MoI$chr),
                      SYMBOL = MoI$SYMBOL,
                      call = MoI$call,
                      variant_type = MoI$variant_type,
                      label = MoI$label,
                      MoI_stories)
    MoI_long <- MoI %>%
      tidyr::gather(key = NAME, value = VAF, samples) %>%
      dplyr::left_join(metadata[,c("NAME","TIMEPOINT","INDIVIDUAL")])


    #remove chromosome from the label
    MoI_long$label = gsub('( ?)\\(.*\\)', '', MoI_long$label)

    tidy_superfreq <- MoI_long %>%
      dplyr::mutate(mutation_key = MoI_long$call) %>%
      dplyr::mutate(mutation_det = MoI_long$label) %>%
      dplyr::rename(SampleName = NAME,
                    PID = INDIVIDUAL,
                    Time = TIMEPOINT) %>%
      tidyr::separate(mutation_det, into = c("SYMBOL"),sep=" ",remove=FALSE)

    # add total read depth. qs[[sample]][SNVoI,]$cov

    if(length(SNVoInames) > 0){

      recover_infos <- lapply(samples,
                              function(sample) {
                                sample_rec <- data.frame(tot_depth = qs[[sample]][SNVoInames,]$cov,
                                                         ref_depth = qs[[sample]][SNVoInames,]$ref,
                                                         ref = qs[[sample]][SNVoInames,]$reference,
                                                         alt = qs[[sample]][SNVoInames,]$variant,
                                                         SampleName = sample) %>%
                                  dplyr::mutate(alt_depth = tot_depth - ref_depth)
                                return(sample_rec)
                              })

      recover_infos <-  bind_rows(recover_infos) %>%
        dplyr::left_join(tidy_superfreq)

    } else {

      recover_infos <- lapply(samples,
                              function(sample) {
                                sample_rec <- data.frame(tot_depth = NA,
                                                         ref_depth = NA,
                                                         ref = NA,
                                                         alt = NA,
                                                         alt_depth = NA,
                                                         SampleName = sample)
                                return(sample_rec)
                              })

        recover_infos <- bind_rows(recover_infos) %>%
        dplyr::left_join(tidy_superfreq)

    }

    # Extract clones - to be added to extract genes in clones
    #labels = superFreq:::storyToLabel(MoI, stories$variants, genome=ref_genome, mergeCNAs=F)[as.character(seq(along.with=MoI$x1)),]

    storyMx = stories$stories[[1]]$consistentClusters$cloneStories$stories
    anchors = c(rownames(stories$stories[[1]]$anchorStories$anchorSNVs),
                rownames(stories$stories[[1]]$anchorStories$anchorCNAs))
    Nmut = sapply(stories$stories[[1]]$consistentClusters$storyList, function(muts) sum(muts %in% anchors))

    use = names(Nmut) != 'germline'

    ret = list(mutations=paste0('clone (', Nmut[use], ' anchors)'), y_matrix=storyMx[use,,drop=F])
    ret$mutations =  ret$mutations[matrixStats::rowMaxs(ret$y_matrix) > VAFcut]
    ret$y_matrix =  ret$y_matrix[matrixStats::rowMaxs(ret$y_matrix) > VAFcut,,drop=F]

    if ( length(ret$mutations) == 0 ) {
      message(paste0("No clones found in patient ",patientID))
      return(NULL)
    }

    tidy_clones <- data.frame(ret$y_matrix) %>%
      dplyr::mutate(mutation_key = rownames(ret$y_matrix)) %>%
      dplyr::mutate(mutation_det = ret$mutations) %>%
      tidyr::gather(NAME,VAF, 1:ncol(ret$y_matrix)) %>%
      dplyr::left_join(metadata[,c("NAME","TIMEPOINT","INDIVIDUAL")]) %>%
      dplyr::rename(SampleName = NAME,
                    PID = INDIVIDUAL,
                    Time = TIMEPOINT) %>%
      dplyr::mutate(variant_type = "clones")


    # To add - extract infos about stories
    # fix how I am merging clones and snvs

    clones_moi <- bind_rows(tidy_clones,recover_infos)

    # put this MoI in long format
    cat('found ', nrow(clones_moi), ' mutations/clones of interest.\n', sep='')

    return(clones_moi)

  }

}
