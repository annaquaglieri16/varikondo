#' Parse a VCF file and return a data frame with standardised fields across callers to use for caller comparison.
#' @param vcf_path path to where the `.vcf` file for one sample is saved.
#' @param sample_name character. Sample name of the current `vcf` file.
#' @param caller character. One of `mutect`, `vardict` or `varscan`.
#' @param vep logical. If TRUE, the annotation fields added by the Variant Effect Predictor will be parsed.

#' @description It only works woth germline calls and for VCF from the following callers: GATK3 MuTect2, VarScan2 and VarDict.

#' @export


# vcf_path <- "../../../cbf_aml_agrf/variant_calling/vardict/regions_deDupl_both_cohorts/annotated_variants/10.R1.B2.M13ADE05RV.BM.Rem_germline_annotated.vcf"
# sample_name = "10.R1.B2.M13ADE05RV.BM.Rem"
# caller = "vardict"
# vep = TRUE

parse_vcf_output <- function(vcf_path, sample_name, caller, vep = TRUE) {

  vcf <- VariantAnnotation::readVcf(vcf_path)

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
        as_tibble()

    }

    if(caller == "mutect"){

      # each row contains depth for the ref and alt allele, in order
      allele_depths <- do.call(rbind,VariantAnnotation::geno(vcf)$AD[,1])
      base_quality = do.call(rbind,VariantAnnotation::geno(vcf)$QSS[,1])

      vcf_df <- data.frame(data.frame(IRanges::ranges(vcf)),
                           genotype= VariantAnnotation::geno(vcf)$GT[,1],
                           filter = VariantAnnotation::filt(vcf),
                           base_quality = base_quality[,1],
                           ref_depth = allele_depths[,1],
                           alt_depth = allele_depths[,2],
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
                           #SampleName = VariantAnnotation::info(vcf)$SAMPLE,
                           qual = VariantAnnotation::info(vcf)$QUAL,
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
        as_tibble()

    }

  }

  vcf_df <- vcf_df %>% dplyr::select(Location,caller,chrom,pos,end,ref,alt,qual,filter,
                                     genotype,tot_depth,VAF,ref_depth,
                                     alt_depth,ref_forw,ref_rev,alt_forw,alt_rev,everything())


  if(vep){

    parsed_vep <- parse_vep_csq(vcf_path = vcf_path, vcf_df = vcf_df) %>%
      as_tibble()

    return(parsed_vep)

  } else {

    return(vcf_df)

  }




}


# Extract VAF from VarScan output
parse_vaf_varscan <- function(freq){
  freq <- gsub("%","",freq)
  freq <- as.numeric(as.character(freq))/100
  return(freq)
}

