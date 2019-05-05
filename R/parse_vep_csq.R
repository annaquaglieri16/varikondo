#' Parse the CSQ field added to the `VCF` file by the Variant Effect Predictor (VEP)
#' @param vcf_path path to where the `.vcf` file for one sample is saved.
#' @param vcf_df `R` object returned by `parse_vcf_output()` after parsing relevant fields from the `VCF` file.

#' @export

#vcf_path <- "../../../cbf_aml_agrf/variant_calling/vardict/regions_deDupl_both_cohorts/annotated_variants/10.R1.B2.M13ADE05RV.BM.Rem_germline_annotated.vcf"
# sample_name = "10.R1.B2.M13ADE05RV.BM.Rem"
# caller = "vardict"
# vep = TRUE

#' @importFrom VariantAnnotation readVcf
#' @importFrom VariantAnnotation header
#' @importFrom VariantAnnotation info
#' @import dplyr
#' @importFrom tidyselect everything



parse_vep_csq <- function(vcf_path,vcf_df){

  # Read VCF file
  vcf <- VariantAnnotation::readVcf(vcf_path)
  # Check presence of VEP CSQ
  des <- VariantAnnotation::info(VariantAnnotation::header(vcf))
  if(sum(rownames(des) == "CSQ") < 1){

    return(vcf_df)
    stop("No VEP annotation present.")

  } else {

    # Parse VEP output in VCF file
    des_vep <- des$Description[rownames(des) == "CSQ"]
    des_vep <- gsub("Consequence annotations from Ensembl VEP. Format: ","",des_vep)
    des_vep_names <- strsplit(des_vep,split="|",fixed=TRUE)[[1]]

    # Duplicate rows of dataframe to allow VEP annotation
    nrepl <- sapply(VariantAnnotation::info(vcf)$CSQ,length)
    df_repl <- vcf_df %>% dplyr::slice(base::rep(1:dplyr::n(), times = nrepl))
    vcf_df <- df_repl %>% dplyr::mutate(INFO_VEP = unlist(VariantAnnotation::info(vcf)$CSQ)) %>%
      tidyr::separate(INFO_VEP,into = des_vep_names,sep="[|]")


    # Order IMPACT based on https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html
    vcf_df <- vcf_df %>% dplyr::mutate(IMPACT = as.character(IMPACT)) %>%
      dplyr::mutate(IMPACT = ifelse(IMPACT == "" , NA, IMPACT)) %>%
      dplyr::mutate(IMPACT = factor(IMPACT, levels = c("MODIFIER","LOW","MODERATE","HIGH"))) %>%
      dplyr::mutate(IMPACT_rank = as.character(as.numeric(IMPACT)))

    vcf_parsed <- vcf_df %>% dplyr::select(Location,caller,chrom,pos,end,ref,alt,qual,filter,
                                           genotype,tot_depth,VAF,ref_depth,
                                           alt_depth,ref_forw,ref_rev,alt_forw,alt_rev,everything())

    return(vcf_parsed)

  }

}
