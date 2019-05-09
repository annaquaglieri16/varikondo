#' Extract SNVs, CNVs and clones on genes of interest from superFreq output
#'
#' @param superFreq_R_path Path to superFreq R folder.
#' @param superFreq_meta_path Path to superFreq run cohort metadata (file ending with .tsv). If the sample names in the NAME field of the metadata file start with a number, you need to add an "X" before the name with paste0() since this is done by default by superFReq.
#' @param studyGenes character vector containing genes of interest. If none provided all genes will be used.
#' @param patientID a character vector specifying the patient/s id/s for which variants have to be imported.
#' @param ref_genome character vector for the reference genome used in the analysis ('hg38' or 'hg19')
#' @param min_vaf numeric. Minimum variant allele frequency (VAF) for a variant to be kept at one time point.
#' @param min_alt numeric. Minimum number of reads supporting the alt allele at one time points for a patient.
#' @param severity numeric. Only variants with severity value < severity are kept.
#' @description This function imports superFreq's SNVs, CNAs and clones for one patient and outputs them into a tidy format where every row is a variant.
#' @export

# studyGenes <- read_csv("../../../venetoclax_trial/Recurrent-AML-genes-across-studies.csv")
# studyGenes <- "KIT"
# superFreq_R_path <- "../../../venetoclax_trial/superFreq/R"
# patientID <- "D1"
# min_vaf <- 0.15
# superFreq_meta_path <- file.path("../../../venetoclax_trial/superFreq/metaData.tsv")
# ref_genome = "hg38"
#
# genes0 <- read.csv("../../../venetoclax_trial/Recurrent-AML-genes-across-studies.csv")
# studyGenes <- "KIT"
# superFreq_R_path = "../../../cbf_aml_agrf/superFreq/R_full_cohort/"
# superFreq_meta_path = "../../../cbf_aml_agrf/superFreq/runFullCohort/metadata_varscan.tsv"
# patients <- unique(meta$INDIVIDUAL)
# table(meta$TIMEPOINT)
# patientID <- "RMH07PW"
# ref_genome = "hg38"
# min_vaf = 0.15


#' @import dplyr
#' @importFrom readr read_delim
#' @import tidyr

import_goi_superfreq <- function(superFreq_R_path = superFreq_R_path,
                                 superFreq_meta_path = superFreq_meta_path,
                                 studyGenes = NULL,
                                 patientID = "D1",
                                 ref_genome = "hg38",
                                 min_vaf = 0.15,
                                 min_alt = 2,
                                 severity = 12){

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
  metadata <- readr::read_delim(superFreq_meta_path,delim = "\t",progress = FALSE) %>%
    dplyr::mutate(NAME = make.names(NAME))


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

  # Use all genes for this patient if no studyGenes provided
  bind_qs <- dplyr::bind_rows(qs)
  all_genes <- as.character(unique(bind_qs$inGene))
  if(is.null(studyGenes)){
    cat("All genes are used...")
    studyGenes <- all_genes } else {
      studyGenes <- studyGenes
    }

  SNVoInames = unique(
    unlist(
      lapply(
        qs, function(q) {

          q <- q %>%
            dplyr::mutate(VAF = var/cov,
                          mut_name = rownames(q)) %>%
            dplyr::filter((somaticP > 0.5) &
                            (!is.na(q$severity) & severity < severity) &
                            inGene %in% studyGenes &
                            VAF >= min_vaf &
                            (cov - ref) >= min_alt)
          as.character(q$mut_name)

        }
      )
    )
  )

  SNVoIannotation = do.call(rbind, lapply(SNVoInames, function(snv) {
    qs[[which(sapply(qs, function(q) q[snv,]$somaticP > 0))[1]]][snv,]
  }))

  # DEfine what are the SNVs of interest - at this stage the dataframe includes only VAF as variant information
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
    labels = storyToLabel(MoI, stories$variants, genome=ref_genome, mergeCNAs=F)[as.character(seq(along.with=MoI$x1)),]
    MoI$label = labels$label
    MoI$severity = labels$severity
    MoI$chr = xToChr(MoI$x1, genome=ref_genome)
    MoI$start = xToPos(MoI$x1, genome=ref_genome)
    MoI$end = xToPos(MoI$x2, genome=ref_genome)

    # Convert to long format SNVs
    MoI_stories <- data.frame(MoI$stories)
    MoI <- data.frame(pos = MoI$start,end = MoI$end,
                      chrom = paste0("chr",MoI$chr),
                      SYMBOL = as.character(MoI$SYMBOL),
                      call = as.character(MoI$call),
                      variant_type = as.character(MoI$variant_type),
                      label = as.character(MoI$label),
                      MoI_stories)

    # Converting to long format to merge with metadata file
    MoI_long <- MoI %>%
      tidyr::gather(key = NAME, value = VAF, samples) %>%
      dplyr::left_join(metadata[,c("NAME","TIMEPOINT","INDIVIDUAL")])

    # Reshape CNA to long
    # remove chromosome from the label
    MoI_long$label = gsub('( ?)\\(.*\\)', '', MoI_long$label)
    MoI_long$call <- as.character(MoI_long$call)

    tidy_superfreq <- MoI_long %>%
      dplyr::mutate(mutation_det = label) %>%
      dplyr::rename(SampleName = NAME,
                    PID = INDIVIDUAL,
                    Time = TIMEPOINT)

    #

    # At this stahe I need to recover the variant information for the SNVs found. This includes  tot_depth, ref_depth ect..
    if(length(SNVoInames) > 0){

      recover_infos <- lapply(samples,
                              function(sample) {

                                sample_rec <- qs[[sample]][SNVoInames,] %>%
                                  dplyr::rename(tot_depth = cov,
                                                ref_depth = ref,
                                                ref = reference,
                                                alt = variant,
                                                SYMBOL = inGene) %>%
                                  dplyr::mutate(alt_depth = tot_depth - ref_depth,
                                                call = SNVoInames,
                                                SampleName = sample) %>%
                                  dplyr::select(-x)

                                return(sample_rec)
                              })

      recover_infos <-  bind_rows(recover_infos) %>%
        dplyr::right_join(tidy_superfreq) # Mrge by SampleName, Symbol and mutation key

    } else {

      recover_infos <- tidy_superfreq

    }

    # Add mutation key
    recover_infos$SYMBOL <- as.character(recover_infos$SYMBOL)
    recover_infos$variant_type <- as.character(recover_infos$variant_type)
    recover_infos$mutation_key <- ifelse(recover_infos$variant_type %in% "SNV",
                                         paste(recover_infos$SYMBOL, recover_infos$pos,
                                               recover_infos$ref, recover_infos$alt, sep = "-"),NA)

    recover_infos$mutation_key <- ifelse(recover_infos$variant_type %in% "CNA",
                                         paste(recover_infos$SYMBOL, recover_infos$call,sep = "-"),recover_infos$mutation_key)

    recover_infos$mutation_det <- ifelse(recover_infos$variant_type %in% "CNA",
                                         paste(recover_infos$SYMBOL, recover_infos$mutation_det,sep = " "),recover_infos$mutation_det)

    recover_infos <- recover_infos %>%
      dplyr::select(-call,-label)

    # Extract clones
    # - to be added to extract genes in clones
    #labels = superFreq:::storyToLabel(MoI, stories$variants, genome=ref_genome, mergeCNAs=F)[as.character(seq(along.with=MoI$x1)),]

    # contains mutation_det for every mutation in the clone
    # stories$stories[[1]]$consistentClusters$storyList
    # I should add this as an extra list

    # Extract clones
    storyMx = stories$stories[[1]]$consistentClusters$cloneStories$stories
    anchors = c(rownames(stories$stories[[1]]$anchorStories$anchorSNVs),
                rownames(stories$stories[[1]]$anchorStories$anchorCNAs))
    Nmut = sapply(stories$stories[[1]]$consistentClusters$storyList, function(muts) sum(muts %in% anchors))

    # Only keep somatic clone
    use = names(Nmut) != 'germline'

    ret = list(mutations=paste0('clone (', Nmut[use], ' anchors)'), y_matrix=storyMx[use,,drop=F])
    ret$mutations =  ret$mutations[matrixStats::rowMaxs(ret$y_matrix) > min_vaf]
    ret$y_matrix =  ret$y_matrix[matrixStats::rowMaxs(ret$y_matrix) > min_vaf,,drop=F]

    if ( length(ret$mutations) == 0 ) {
      message(paste0("No clones found in patient ",patientID))
      return(NULL)
    }

    # Convert to long format and join to metadata
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
