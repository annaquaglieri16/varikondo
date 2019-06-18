#' Parse a VCF file and return a data frame with standardised fields across callers
#'
#' @param vcf_path path to where the `.vcf` file for one sample is saved.
#' @param sample_name character. Sample name of the current `vcf` file.
#' @param caller character. One of `mutect`, `vardict`, `varscan`, or `freebayes`.
#' @param vep logical. If TRUE, the annotation fields added by the Variant Effect Predictor will be parsed.
#' @param param same as `param` in `VariantAnnotation::readVcf` to subset the VCF file by genomic coordinate and only import specific regions. An instance of ScanVcfParam or GRanges.

#' @description Currently, it only works with VCF files containing germline calls from the following callers: GATK3 MuTect2, VarScan2, VarDict and freebayes. It uses the Bioconductor package `VariantAnnotation` to read `VCF` files into `R`.
#'
#' @details Freebayes can report more than one alternative allele in output. This means that there will be depths and quality information for every alternative allele. Currently, `parse_vcf_output` uses the `VariantAnnotation` package to read `VCF` fields into `R` but, if multiple entries are reported in one field (alt allele, quality, depth etc..), it only reports the first of them. This should be fixed soon within the `VariantAnnotation` package but in the meantime `parse_vcf_output` parses these fields separately and adds them to the final output. Since, not many variants have a second or thirs alternative allele, the `qual` column reported by `parse_vcf_output` is the sum of the reference base qualities and the first (most common) alternative base qualities divided by the sum of reference and alternative depth for the first allele. The same applies for the variant allele frequency (VAF).
#'
#' @return data frame with standardised fields containing all the variants in the input VCF.
#'
#' @export
#' @examples
#' vcf_path <- system.file("extdata", "chr20_mutect.vcf.gz", package = "varikondo")
#' parsed_vcf_mutect <- parse_vcf_output(vcf_path,
#' caller = "mutect",
#' sample_name = "Sample1",
#' vep = TRUE)
#'
#' vcf_path <- system.file("extdata", "chr20_freebayes.vcf.gz", package = "varikondo")
#' parsed_vcf_freebayes <- parse_vcf_output(vcf_path,
#' caller = "freebayes",
#' sample_name = "Sample1",
#' vep = TRUE)

#' @importFrom VariantAnnotation readVcf
#' @importFrom VariantAnnotation geno
#' @importFrom VariantAnnotation filt
#' @importFrom IRanges ranges
#' @import dplyr
#' @import stringr
#' @importFrom data.table fread
#' @import tidyr
#' @importFrom DelayedArray rowRanges
#' @importFrom tidyselect everything

parse_vcf_output <- function(vcf_path, sample_name = basename(vcf_path), caller, vep = TRUE, param = VariantAnnotation::ScanVcfParam()) {

  # The VCF file needs to be tabixed
  if(!str_detect(pattern = ".bgz$|.gz$", string = vcf_path)){
    vcf_path <- Rsamtools::bgzip(vcf_path,overwrite = TRUE)
    Rsamtools::indexTabix(vcf_path,format="vcf")
    message("The VCF is compressed (Rsamtools::bgzip) and indexed (Rsamtools::indexTabix)")
  } else {
    if(!file.exists(paste0(vcf_path,".tbi"))) {
      Rsamtools::indexTabix(vcf_path,format="vcf")
      message("The VCF is indexed (Rsamtools::indexTabix)")
    }
  }

  # Read VCF in
  vcf <- VariantAnnotation::readVcf(vcf_path, param = param)

  # Check if is comes from somatic calls
  if(ncol(VariantAnnotation::geno(vcf)$GT) > 1){

    stop("parse_vcf_output wasn't implemented for somatic calls")

  } else {

    if(caller == "varscan"){

      vcf_df <- data.frame(data.frame(IRanges::ranges(vcf)),
                           genotype= VariantAnnotation::geno(vcf)$GT[,1],
                           filter = VariantAnnotation::filt(vcf),
                           ref_base_quality = VariantAnnotation::geno(vcf)$RBQ[,1],
                           alt_base_quality = VariantAnnotation::geno(vcf)$ABQ[,1],
                           tot_depth = VariantAnnotation::geno(vcf)$DP[,1],
                           freq = VariantAnnotation::geno(vcf)$FREQ[,1],
                           ref_depth = VariantAnnotation::geno(vcf)$RD[,1],
                           alt_depth = VariantAnnotation::geno(vcf)$AD[,1],
                           ref_forw = VariantAnnotation::geno(vcf)$RDF[,1],
                           ref_rev = VariantAnnotation::geno(vcf)$RDR[,1],
                           alt_forw = VariantAnnotation::geno(vcf)$ADF[,1],
                           alt_rev = VariantAnnotation::geno(vcf)$ADR[,1]) %>%

        tidyr::separate(names,into=c("Location","alleles"),sep="_") %>%
        tidyr::separate(Location,into=c("chrom","pos"),sep=":",remove=FALSE) %>%
        tidyr::separate(alleles,into=c("ref","alt"),sep="/") %>%
        dplyr::mutate(caller="varscan",
                      Location = gsub(":","_",Location),
                      SampleName = sample_name) %>%
        dplyr::mutate(VAF = parse_vaf_varscan(freq),
                      qual = (ref_base_quality + alt_base_quality)/2) %>%  # mean of alt/ref base qualitites
        dplyr::select(-freq) %>%
        dplyr::as_tibble()

    }

    if(caller == "mutect"){

        # Extract QSS field which is not parsed correctly by VariantAnnotation
        readt <- data.table::fread(vcf_path) %>%
          dplyr::as_tibble()
        split_format <- strsplit(as.character(readt$FORMAT),split=":")

        which_qss <- sapply(split_format, function(format) which(format == "QSS"))
        split_format_entries <- strsplit(as.character(readt[,ncol(readt),drop=TRUE]),split=":")
        get_qss <- mapply('[[',split_format_entries,which_qss)

        # each row contains depth for the ref and alt allele, in order
        allele_depths <- do.call(rbind,VariantAnnotation::geno(vcf)$AD[,1])

        vcf_df <- data.frame(IRanges::ranges(vcf)) %>%
         dplyr::mutate(genotype= VariantAnnotation::geno(vcf)$GT[,1],
                       filter = VariantAnnotation::filt(vcf),
                       ref_depth = allele_depths[,1],
                       alt_depth = allele_depths[,2],
                       base_quality = get_qss,
                       VAF = VariantAnnotation::geno(vcf)$AF[,1],
                       alt_forw = VariantAnnotation::geno(vcf)$ALT_F1R2[,1],
                       alt_rev = VariantAnnotation::geno(vcf)$ALT_F2R1[,1],
                       ref_forw = VariantAnnotation::geno(vcf)$REF_F1R2[,1],
                       ref_rev = VariantAnnotation::geno(vcf)$REF_F2R1[,1]) %>%

         dplyr::mutate(ref_depth =  as.numeric(as.character(ref_depth)),
                         alt_depth = as.numeric(as.character(alt_depth)),
                         tot_depth = ref_depth + alt_depth,
                         SampleName = sample_name) %>%

         tidyr::separate(base_quality,into = c("ref_base_quality","alt_base_quality"),sep = ",",remove=TRUE) %>%
         dplyr::mutate(ref_base_quality =  as.numeric(as.character(ref_base_quality)),
                         alt_base_quality = as.numeric(as.character(alt_base_quality))) %>%
         dplyr::mutate(qual_ref = ifelse(ref_depth == 0,0,ref_base_quality/ref_depth),
                         qual_alt= ifelse(alt_depth == 0,0,alt_base_quality/alt_depth)) %>%

           dplyr::mutate(qual = qual_ref + qual_alt/2) %>%
           dplyr::select(-qual_ref,-qual_alt) %>%

           tidyr::separate(names,into=c("Location","alleles"),sep="_") %>%
           tidyr::separate(Location,into=c("chrom","pos"),sep=":",remove=FALSE) %>%
           tidyr::separate(alleles,into=c("ref","alt"),sep="/") %>%
           dplyr::mutate(caller="mutect2",
                         Location = gsub(":","_",Location)) %>%
           as_tibble()

    }


    if(caller == "vardict"){

      vcf_df <- data.frame(data.frame(IRanges::ranges(vcf)),
                           genotype= VariantAnnotation::geno(vcf)$GT[,1],
                           qual = VariantAnnotation::info(vcf)$QUAL,
                           #qual = rowRanges(vcf)$QUAL,
                           filter = VariantAnnotation::filt(vcf),
                           DP = VariantAnnotation::info(vcf)$DP,
                           VAF = VariantAnnotation::info(vcf)$AF,
                           ADJVAF_ADJ_indels = VariantAnnotation::info(vcf)$ADJAF,
                           VD = VariantAnnotation::info(vcf)$VD,
                           REFBIAS = VariantAnnotation::info(vcf)$REFBIAS,
                           VARBIAS = VariantAnnotation::info(vcf)$VARBIAS) %>%

        tidyr::separate(names,into=c("Location","alleles"),sep="_") %>%
        tidyr::separate(Location,into=c("chrom","pos"),sep=":",remove=FALSE) %>%
        tidyr::separate(alleles,into=c("ref","alt"),sep="/") %>%
        tidyr::separate(REFBIAS,into=c("ref_forw","ref_rev"),sep=":") %>%
        tidyr::separate(VARBIAS,into=c("alt_forw","alt_rev"),sep=":") %>%
        dplyr::mutate(ref_depth = DP - VD,
                      SampleName = sample_name) %>%
        dplyr::rename(tot_depth = DP,
                      alt_depth = VD) %>%
        dplyr::select(-start) %>%
        dplyr::mutate(caller="vardict",
                      Location = gsub(":","_",Location))%>%
        dplyr::as_tibble()

    }


    if(caller == "freebayes"){

      # Paste together information relative to various alternative alleles
      paste_qual <- lapply(VariantAnnotation::info(vcf)$QA,
                           function(var) paste(as.character(var),collapse=","))
      df_qual <- do.call(c,paste_qual)

      paste_depth <- lapply(VariantAnnotation::info(vcf)$AO,
                            function(var) paste(as.character(var),collapse=","))
      df_depth <- do.call(c,paste_depth)

      paste_AF <- lapply(VariantAnnotation::info(vcf)$SAF,
                         function(var) paste(as.character(var),collapse=","))
      df_AF <- do.call(c,paste_AF)

      paste_AR <- lapply(VariantAnnotation::info(vcf)$SAR,
                         function(var) paste(as.character(var),collapse=","))
      df_AR <- do.call(c,paste_AR)

      paste_alt <- lapply(DelayedArray::rowRanges(vcf)$ALT,
                         function(var) paste(as.character(var),collapse=","))
      df_alt <- do.call(c,paste_alt)

      vcf_df <- data.frame(data.frame(IRanges::ranges(vcf)),
                           ref = DelayedArray::rowRanges(vcf)$REF,
                           alt = as.character(df_alt),
                           genotype = VariantAnnotation::geno(vcf)$GT[,1],
                           qual_ref = VariantAnnotation::info(vcf)$QR,
                           qual_alt = df_qual,
                           ref_depth = VariantAnnotation::info(vcf)$RO,
                           alt_depth = as.character(df_depth),
                           filter = VariantAnnotation::filt(vcf),
                           tot_depth = VariantAnnotation::info(vcf)$DP,
                           ref_forw = VariantAnnotation::info(vcf)$SRF,
                           ref_rev = VariantAnnotation::info(vcf)$SRR,
                           alt_forw = as.character(df_AF),
                           alt_rev = as.character(df_AR)) %>%

        tidyr::separate(names,into=c("Location"),sep="_") %>%
        tidyr::separate(Location,into=c("chrom","pos"),sep=":",remove=FALSE) %>%
        dplyr::mutate(SampleName = sample_name) %>%

        tidyr::separate(alt_depth, into = c("alt_depth1"),sep=",",remove=FALSE) %>%
        dplyr::mutate_at(.vars = c("alt_depth1"),as.numeric) %>%
        dplyr::mutate(VAF = (alt_depth1)/as.numeric(ref_depth)) %>%

        tidyr::separate(qual_alt, into = c("qual_alt1"),sep=",",remove=FALSE) %>%
        dplyr::mutate_at(.vars = c("qual_alt1"),as.numeric) %>%
        dplyr::mutate(qual = (as.numeric(qual_ref) + qual_alt1)/(qual_alt1 + ref_depth)) %>%

        dplyr::select(-start,-qual_alt1,-alt_depth1) %>%
        dplyr::mutate(caller="freebayes",
                      Location = gsub(":","_",Location))%>%
        dplyr::as_tibble()

    }

  }

  vcf_df <- vcf_df %>% dplyr::select(Location,caller,chrom,pos,end,ref,alt,qual,filter,
                                     genotype,tot_depth,VAF,ref_depth,
                                     alt_depth,ref_forw,ref_rev,alt_forw,alt_rev,everything())


  if(vep){

    parsed_vep <- parse_vep_csq(vcf_path = vcf_path, vcf_df = vcf_df) %>%
      dplyr::as_tibble()

    return(parsed_vep)

  } else {

    return(vcf_df)

  }




}


