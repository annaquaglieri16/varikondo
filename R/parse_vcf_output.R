#' Parse a VCF file and return a data frame with standardised fields across callers to use for caller comparison.
#' @param vcf_path path to where the `.vcf` file for one sample is saved.
#' @param sample_name character. Sample name of the current `vcf` file.
#' @param caller character. One of `mutect`, `vardict` or `varscan`.
#' @param vep logical. If TRUE, the annotation fields added by the Variant Effect Predictor will be parsed.

#' @description Currently, it doesn't work for somatic calls and for VCF from the following callers: GATK3 MuTect2, VarScan2 and VarDict. It used the Bioconductor package `VariantAnnotation` to read `VCF` files into `R`.

#' @export
#' @examples
#' vcf_path <- system.file("extdata", "germline_mutect.vcf", package = "varikondo")
#' parsed_vcf_mutect <- varikondo::parse_vcf_output(annot_vcf_mutect,
#' caller = "mutect",
#' sample_name = "Sample1",
#' vep = TRUE)

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

      # Tmp fix to be updated once variantannotation is upadted
      # base_quality = do.call(rbind,VariantAnnotation::geno(vcf)$QSS[,1])
     # read_qss <- read.table(vcf_path,stringsAsFactors = FALSE) %>%
      #  as_tibble() %>%
      #  dplyr::mutate(names = paste0(V1,":",V2,"_",V4,"/",V5))

  #  ranges <- data.frame(IRanges::ranges(vcf))

   #   if(sum(read_qss$names != ranges$names) > 0){

        warning("This is a tmp fix. MuTect2 VCF will be parsed only with ref base quality.")

        # each row contains depth for the ref and alt allele, in order
        allele_depths <- do.call(rbind,VariantAnnotation::geno(vcf)$AD[,1])
        base_qualities <- do.call(rbind,VariantAnnotation::geno(vcf)$QSS[,1])

        vcf_df <- data.frame(IRanges::ranges(vcf)) %>%
          dplyr::mutate(genotype= VariantAnnotation::geno(vcf)$GT[,1],
                        filter = VariantAnnotation::filt(vcf),
                        ref_depth = allele_depths[,1],
                        alt_depth = allele_depths[,2],
                        base_qualities = base_qualities,
                        VAF = VariantAnnotation::geno(vcf)$AF[,1],
                        alt_forw = VariantAnnotation::geno(vcf)$ALT_F1R2[,1],
                        alt_rev = VariantAnnotation::geno(vcf)$ALT_F2R1[,1],
                        ref_forw = VariantAnnotation::geno(vcf)$REF_F1R2[,1],
                        ref_rev = VariantAnnotation::geno(vcf)$REF_F2R1[,1]) %>%

          dplyr::mutate(ref_depth =  as.numeric(as.character(ref_depth)),
                        alt_depth = as.numeric(as.character(alt_depth)),
                        tot_depth = ref_depth + alt_depth,
                        SampleName = sample_name,
                        base_qualities = as.numeric(as.character(base_qualities))) %>%

          dplyr::mutate(qual = ifelse(ref_depth == 0,0,base_qualities/ref_depth)) %>% # tmp fix
          dplyr::select(-base_qualities) %>%

          tidyr::separate(names,into=c("Location","alleles"),sep="_") %>%
          tidyr::separate(Location,into=c("chrom","pos"),sep=":",remove=FALSE) %>%
          tidyr::separate(alleles,into=c("ref","alt"),sep="/") %>%
          dplyr::mutate(caller="mutect2",
                        Location = gsub(":","_",Location)) %>%
          as_tibble()


    }

      # else {
      #
      #   # each row contains depth for the ref and alt allele, in order
      #   allele_depths <- do.call(rbind,VariantAnnotation::geno(vcf)$AD[,1])
      #   base_qualities <- do.call(rbind,VariantAnnotation::geno(vcf)$QSS[,1])
      #
      #   #names_fields <- read_qss$V9
      #   #names_fields <- strsplit(as.character(names_fields),split=":")[[1]]
      #   #read_qss <- read_qss %>%
      #   #  tidyr::separate(V10, into = names_fields,sep="[:]",remove=FALSE) %>%
      #   #  dplyr::rename(base_quality = QSS) %>%
      #   #  dplyr::select(base_quality,names)
      #
      #   # Bummer!
      #   # Some fields have more entries than others
      #   # [318] "GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1"
      #   # [319] "GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:PGT:PID:QSS:REF_F1R2:REF_F2R1"
      #   # [320] "GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:PGT:PID:QSS:REF_F1R2:REF_F2R1"
      #
      #   vcf_df <- data.frame(IRanges::ranges(vcf)) %>%
      #   dplyr::mutate(genotype= VariantAnnotation::geno(vcf)$GT[,1],
      #                 filter = VariantAnnotation::filt(vcf),
      #                 ref_depth = allele_depths[,1],
      #                 alt_depth = allele_depths[,2],
      #                 #base_qualities = base_qualities,
      #                 qual = base_qualities,
      #                 VAF = VariantAnnotation::geno(vcf)$AF[,1],
      #                 alt_forw = VariantAnnotation::geno(vcf)$ALT_F1R2[,1],
      #                 alt_rev = VariantAnnotation::geno(vcf)$ALT_F2R1[,1],
      #                 ref_forw = VariantAnnotation::geno(vcf)$REF_F1R2[,1],
      #                 ref_rev = VariantAnnotation::geno(vcf)$REF_F2R1[,1]) %>%
      #
      #   dplyr::mutate(ref_depth =  as.numeric(as.character(ref_depth)),
      #                   alt_depth = as.numeric(as.character(alt_depth)),
      #                   tot_depth = ref_depth + alt_depth,
      #                   SampleName = sample_name) %>%
      #
      #  # tidyr::separate(base_quality,into = c("ref_base_quality","alt_base_quality"),sep = ",",remove=TRUE) %>%
      # #  dplyr::mutate(ref_base_quality =  as.numeric(as.character(ref_base_quality)),
      # #                  alt_base_quality = as.numeric(as.character(alt_base_quality))) %>%
      # #  dplyr::mutate(qual_ref = ifelse(ref_depth == 0,0,ref_base_quality/ref_depth),
      # #                  qual_alt= ifelse(alt_depth == 0,0,alt_base_quality/alt_depth)) %>%
      #
      #   #  dplyr::mutate(qual = qual_ref + qual_alt/2) %>%
      #   #  dplyr::select(-qual_ref,-qual_alt) %>%
      #
      #     tidyr::separate(names,into=c("Location","alleles"),sep="_") %>%
      #     tidyr::separate(Location,into=c("chrom","pos"),sep=":",remove=FALSE) %>%
      #     tidyr::separate(alleles,into=c("ref","alt"),sep="/") %>%
      #     dplyr::mutate(caller="mutect2",
      #                   Location = gsub(":","_",Location)) %>%
      #     as_tibble()

    #}


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

