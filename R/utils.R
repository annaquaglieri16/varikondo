#' Functions not to be exported

# Extract VAF from VarScan output

#' @param freq Vector of frequencies as recorded in VarScan VCF file, field FREQ.
#' @keywords internal
#' @return numeric vector of variant allele frequecies

parse_vaf_varscan <- function(freq){
  freq <- gsub("%","",freq)
  freq <- as.numeric(as.character(freq))/100
  return(freq)
}


#  Helper Functions from superFreq
#' @keywords internal

storyToLabel <- function(stories, variants, genome, maxLength = 30, mergeCNAs = FALSE,
          annotationMethod = "VariantAnnotation"){
  call = stories$call
  isSNP = grepl("[0-9]", call)
  x1 = stories$x1
  clonename = rownames(stories)
  dist = stories$x2 - stories$x1
  isCNA = !isSNP & !is.na(dist)
  chr = xToChr(stories$x1, genome = genome)
  if (mergeCNAs & any(isCNA)) {
    dups = unique(cbind(call, chr)[isCNA, , drop = FALSE][duplicated(cbind(call,
                                                                       chr)[isCNA, , drop = FALSE]), , drop = FALSE])
    if (nrow(dups) > 0) {
      for (row in 1:nrow(dups)) {
        thisCall = dups[row, "call"]
        thisChr = dups[row, "chr"]
        muts = call == thisCall & chr == thisChr
        thisDist = sum(dist[muts])
        thisX1 = min(x1[muts])
        thisName = clonename[muts][1]
        call = c(thisCall, call[!muts])
        chr = c(thisChr, chr[!muts])
        dist = c(thisDist, dist[!muts])
        x1 = c(thisX1, x1[!muts])
        clonename = c(thisName, clonename[!muts])
      }
      isSNP = grepl("[0-9]", call)
      isCNA = !isSNP & !is.na(dist)
    }
  }
  label = rep("", length(dist))
  font = rep(1, length(dist))
  colour = rep("black", length(dist))
  severity = rep(100, length(dist))
  q = variants$variants[[1]]
  q = q[q$x %in% x1[isSNP], ]
  gene = q[call, ]$inGene[isSNP]
  if ("severity" %in% names(variants$variants[[1]]) & any(isSNP)) {
    severityMx = sapply(variants$variants, function(q) ifelse(is.na(q[call[isSNP],
                                                                      ]$severity), 100, q[call[isSNP], ]$severity))
    if (sum(isSNP) == 1)
      severityMx = matrix(severityMx, nrow = sum(isSNP))
    mostSevere = ceiling(apply(severityMx, 1, min))
    severity[isSNP] = mostSevere
    type = sapply(mostSevere, function(severity) severityToType(severity,
                                                                annotationMethod))
    type[type == "unknown"] = ""
    label[isSNP] = substr(paste0(gene, " (", chr[isSNP],
                                 ") ", type), 1, maxLength)
  }
  else label[isSNP] = substr(paste0(gene, " (", chr[isSNP],
                                    ") "), 1, maxLength)
  if ("isCosmicCensus" %in% names(variants$variants[[1]]) &
      any(isSNP)) {
    cosmicMx = sapply(variants$variants, function(q) {
      if (!("isCosmicCensus" %in% names(q)))
        return(rep(FALSE, sum(isSNP)))
      return(q[call[isSNP], ]$isCosmicCensus & q[call[isSNP],
                                                 ]$severity < 11)
    })
    if (sum(isSNP) == 1)
      cosmicMx = matrix(unlist(cosmicMx), nrow = sum(isSNP))
    isCensus = apply(cosmicMx, 1, any)
    clinvarMx = sapply(variants$variants, function(q) {
      if (!("ClinVar_ClinicalSignificance" %in% names(q)))
        return(rep(FALSE, sum(isSNP)))
      isPathogenic = grepl("[pP]athogenic$", paste0(q[call[isSNP],
                                                      ]$ClinVar_ClinicalSignificance)) | grepl("[pP]athogenic,",
                                                                                               paste0(q[call[isSNP], ]$ClinVar_ClinicalSignificance))
      return(isPathogenic & q[call[isSNP], ]$severity <
               11)
    })
    if (sum(isSNP) == 1)
      clinvarMx = matrix(unlist(clinvarMx), nrow = sum(isSNP))
    isPathogenic = apply(clinvarMx, 1, any)
    isCensus = isCensus | isPathogenic
    colour[isSNP] = ifelse(isCensus & isPathogenic, mcri("orange"),
                           ifelse(isCensus, mcri("green"), "black"))
    font[isSNP] = ifelse(isCensus, 2, 1)
  }
  distText = ifelse(dist >= 1e+06, paste0(round(dist/1e+06),
                                          "Mbp "), ifelse(dist >= 1000, paste0(round(dist/1000),
                                                                               "kbp "), paste0(dist, "bp ")))
  label[!isSNP & !is.na(dist)] = paste0(distText, call, " (",
                                        chr, ")")[!isSNP & !is.na(dist)]
  label[!isSNP & !is.na(dist)] = shortenCalls(label[!isSNP &
                                                      !is.na(dist)])
  font[!isSNP & !is.na(dist)] = ifelse(grepl("[0-9]AB", label[!isSNP &
                                                                !is.na(dist)]), 4, 3)
  severity[!isSNP & !is.na(dist)] = 0
  label[!isSNP & is.na(dist)] = paste0("clone.", clonename[!isSNP &
                                                             is.na(dist)])
  font[!isSNP & is.na(dist)] = 4
  severity[!isSNP & is.na(dist)] = -1
  return(data.frame(label = label, colour = colour, font = font,
                    severity = severity, stringsAsFactors = FALSE)[order(dist,
                                                                     decreasing = TRUE), ])
}

#' @keywords internal
#'
severityToType <- function(severity, annotationMethod = "VEP"){
  if (annotationMethod == "VariantAnnotation")
    return(VAseverityToConsequence(severity))
  if (severity == 1)
    return("transcript_ablation")
  if (severity == 2)
    return("splice_acceptor_variant")
  if (severity == 3)
    return("splice_donor_variant")
  if (severity == 4)
    return("stop_gained")
  if (severity == 5)
    return("frameshift_variant")
  if (severity == 6)
    return("stop_lost")
  if (severity == 7)
    return("start_lost")
  if (severity == 8)
    return("transcript_amplification")
  if (severity == 9)
    return("inframe_insertion")
  if (severity == 10)
    return("inframe_deletion")
  if (severity == 11)
    return("missense_variant")
  if (severity == 12)
    return("protein_altering_variant")
  if (severity == 13)
    return("splice_region_variant")
  if (severity == 14)
    return("incomplete_terminal_codon_variant")
  if (severity == 15)
    return("stop_retained_variant")
  if (severity == 16)
    return("synonymous_variant")
  if (severity == 17)
    return("coding_sequence_variant")
  if (severity == 18)
    return("mature_miRNA_variant")
  if (severity == 19)
    return("5_prime_UTR_variant")
  if (severity == 20)
    return("3_prime_UTR_variant")
  if (severity == 21)
    return("non_coding_transcript_exon_variant")
  if (severity == 22)
    return("intron_variant")
  if (severity == 23)
    return("NMD_transcript_variant")
  if (severity == 24)
    return("non_coding_transcript_variant")
  if (severity == 25)
    return("upstream_gene_variant")
  if (severity == 26)
    return("downstream_gene_variant")
  if (severity == 27)
    return("TFBS_ablation")
  if (severity == 28)
    return("TFBS_amplification")
  if (severity == 29)
    return("TF_binding_site_variant")
  if (severity == 30)
    return("regulatory_region_ablation")
  if (severity == 31)
    return("regulatory_region_amplification")
  if (severity == 32)
    return("regulatory_region_variant")
  if (severity == 33)
    return("feature_elongation")
  if (severity == 34)
    return("feature_truncation")
  if (severity == 35)
    return("intergenic_variant")
  if (severity == 100)
    return("unknown")
  return("unknown")
}


#' @keywords internal
shortenCalls <- function(calls){
  calls = gsub(" AAAAB", " 4AB", calls)
  calls = gsub(" AAAAAB", " 5AB", calls)
  calls = gsub(" AAAAAAB", " 6AB", calls)
  calls = gsub(" AAAAAAAAAB", " 9AB", calls)
  calls = gsub(" AAAAAAAAAAAAAAAAAAAB", " 19AB", calls)
  calls = gsub(" AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB",
               " 39AB", calls)
  calls = gsub(" AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB",
               " 79AB", calls)
  return(calls)
}

#' @keywords internal
humanChrLengths <- function(){
  humanAllChrLengths()[as.character(c(1, 2, 3, 4, 5, 6, 7,
                                      8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
                                      22, "X", "Y", "M"))]
}

#' @keywords internal
humanAllChrLengths <- function(){
  lengths = c(249250621, 106433, 547496, 243199373, 198022430,
              191154276, 590426, 189789, 191469, 180915260, 171115067,
              4622290, 4795371, 4610396, 4683263, 4833398, 4611984,
              4928567, 159138663, 182896, 146364022, 38914, 37175,
              141213431, 90085, 169874, 187035, 36148, 135534747, 135006516,
              40103, 133851895, 115169878, 107349540, 102531392, 90354753,
              81195210, 1680828, 37498, 81310, 174588, 41001, 78077248,
              4262, 59128983, 92689, 159169, 63025520, 48129895, 27682,
              51304566, 155270560, 59373566, 166566, 186858, 164239,
              137718, 172545, 172294, 172149, 161147, 179198, 161802,
              155397, 186861, 180455, 179693, 211173, 15008, 128374,
              129120, 19913, 43691, 27386, 40652, 45941, 40531, 34474,
              41934, 45867, 39939, 33824, 41933, 42152, 43523, 43341,
              39929, 36651, 38154, 36422, 39786, 38502, 16571)
  names(lengths) = c("1", "1_gl000191_random", "1_gl000192_random",
                     "2", "3", "4", "4_ctg9_hap1", "4_gl000193_random", "4_gl000194_random",
                     "5", "6", "6_apd_hap1", "6_cox_hap2", "6_dbb_hap3", "6_mann_hap4",
                     "6_mcf_hap5", "6_qbl_hap6", "6_ssto_hap7", "7", "7_gl000195_random",
                     "8", "8_gl000196_random", "8_gl000197_random", "9", "9_gl000198_random",
                     "9_gl000199_random", "9_gl000200_random", "9_gl000201_random",
                     "10", "11", "11_gl000202_random", "12", "13", "14", "15",
                     "16", "17", "17_ctg5_hap1", "17_gl000203_random", "17_gl000204_random",
                     "17_gl000205_random", "17_gl000206_random", "18", "18_gl000207_random",
                     "19", "19_gl000208_random", "19_gl000209_random", "20",
                     "21", "21_gl000210_random", "22", "X", "Y", "Un_gl000211",
                     "Un_gl000212", "Un_gl000213", "Un_gl000214", "Un_gl000215",
                     "Un_gl000216", "Un_gl000217", "Un_gl000218", "Un_gl000219",
                     "Un_gl000220", "Un_gl000221", "Un_gl000222", "Un_gl000223",
                     "Un_gl000224", "Un_gl000225", "Un_gl000226", "Un_gl000227",
                     "Un_gl000228", "Un_gl000229", "Un_gl000230", "Un_gl000231",
                     "Un_gl000232", "Un_gl000233", "Un_gl000234", "Un_gl000235",
                     "Un_gl000236", "Un_gl000237", "Un_gl000238", "Un_gl000239",
                     "Un_gl000240", "Un_gl000241", "Un_gl000242", "Un_gl000243",
                     "Un_gl000244", "Un_gl000245", "Un_gl000246", "Un_gl000247",
                     "Un_gl000248", "Un_gl000249", "M")
  return(lengths)
}

#' @keywords internal
#'
chrLengths <- function (genome = "hg19"){
  if (genome == "hg19")
    return(humanChrLengths())
  else if (genome == "mm10")
    return(mouseChrLengths())
  else if (genome == "hg38")
    return(hg38ChrLengths())
  else stop("chrLengths doesnt know about genome ", genome,
            "\n")
}

#' @keywords internal
#'
mouseChrLengths <- function(){
  lengths = c(195471971, 182113224, 160039680, 156508116, 151834684,
              149736546, 145441459, 129401213, 124595110, 130694993,
              122082543, 120129022, 120421639, 124902244, 104043685,
              98207768, 94987271, 90702639, 61431566, 171031299, 91744698,
              16299)
  names(lengths) = c("1", "2", "3", "4", "5", "6", "7", "8",
                     "9", "10", "11", "12", "13", "14", "15", "16", "17",
                     "18", "19", "X", "Y", "M")
  return(lengths)
}

#' @keywords internal
#'
hg38ChrLengths <- function(){
  lengths = c(248956422, 242193529, 198295559, 190214555, 181538259,
              170805979, 159345973, 145138636, 138394717, 133797422,
              135086622, 133275309, 114364328, 107043718, 101991189,
              90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
              50818468, 156040895, 57227415)
  names(lengths) = c("1", "2", "3", "4", "5", "6", "7", "8",
                     "9", "10", "11", "12", "13", "14", "15", "16", "17",
                     "18", "19", "20", "21", "22", "X", "Y")
  return(lengths)
}


#' @keywords internal
#'
xToPos <- function(x, genome = "hg19"){
  chr = xToChr(x, genome)
  pos = x - cumsum(chrLengths(genome))[chr] + chrLengths(genome)[chr]
  return(pos)
}

#' @keywords internal
#'
xToChr <- function(x, genome = "hg19"){
  chrL = chrLengths(genome)
  ret = rep("", length(x))
  for (chr in names(chrL)) {
    ret[x > 0 & x < chrL[chr]] = chr
    x = x - chrL[chr]
  }
  return(ret)
}


#' @keywords internal
#'
VAseverityToConsequence <- function(severity){
  consequence = c(`2` = "spliceSite", `4` = "nonsense", `5` = "frameshift",
                  `10` = "nonsynonymous", `16` = "synonymous", `19` = "promoter",
                  `20` = "threeUTR", `21` = "fiveUTR", `22` = "intron",
                  `30` = "intergenic", `100` = "unknown")
  ret = consequence[as.character(severity)]
  ret[is.na(ret)] = "unknown"
  return(ret)
}

#' @keywords internal
#'
mcri <- function(col = "deafult", al = 1){
  if (length(col) == 0 | length(al) == 0)
    return(character(0))
  if (length(col) == 1 && col[1] == "deafult") {
    cat("Use: mcri('colour'), returning an official MCRI colour.\nAvailable MCRI colours are:\n\ndarkblue\nblue\nlightblue\nazure\ngreen\norange\nviolet\ncyan\nred\ndarkred\nmagenta (aka rose).\n\nReturning default blue.\n")
    return(mcri("blue"))
  }
  if (length(col) > 1 & length(al) == 1)
    return(sapply(col, function(c) mcri(c, al)))
  else if (length(col) > 1 & length(al) == length(col))
    return(sapply(1:length(col), function(i) mcri(col[i],
                                                  al[i])))
  else if (length(col) > 1) {
    warning("length of col and al mismatch in mcri()")
    return(sapply(col, function(c) mcri(c, al[1])))
  }
  if (length(col) == 1 & length(al) > 1)
    return(sapply(al, function(alpha) mcri(col, alpha)))
  if (is.numeric(col)) {
    col = (col%%9) + 1
    if (col == 1)
      col = "blue"
    else if (col == 2)
      col = "red"
    else if (col == 3)
      col = "green"
    else if (col == 4)
      col = "magenta"
    else if (col == 5)
      col = "cyan"
    else if (col == 6)
      col = "orange"
    else if (col == 7)
      col = "violet"
    else if (col == 8)
      col = "darkblue"
    else if (col == 9)
      col = "darkred"
    else col = "black"
  }
  ret = 0
  if (col == "darkblue")
    ret = rgb(9/255, 47/255, 94/255, al)
  if (col == "blue")
    ret = rgb(0, 83/255, 161/255, al)
  if (col == "lightblue")
    ret = rgb(0, 165/255, 210/255, al)
  if (col == "azure")
    ret = rgb(0, 173/255, 239/255, al)
  if (col == "green")
    ret = rgb(141/255, 198/255, 63/255, al)
  if (col == "orange")
    ret = rgb(244/255, 121/255, 32/255, al)
  if (col == "violet")
    ret = rgb(122/255, 82/255, 199/255, al)
  if (col == "cyan")
    ret = rgb(0/255, 183/255, 198/255, al)
  if (col == "red")
    ret = rgb(192/255, 80/255, 77/255, al)
  if (col == "darkred")
    ret = rgb(96/255, 40/255, 38/255, al)
  if (col == "magenta" | col == "rose")
    ret = rgb(236/255, 0/255, 140/255, al)
  if (ret == 0)
    ret = do.call(rgb, as.list(c(col2rgb(col)/255, al)))
  return(ret)
}
