# example variants

# P1 initially has 3 variants
# chr1_10 is refractory
# chr1_20 responds and comes back
# chr2_30 goes away

example_variants <- data.frame(PID = rep(c("P1","P2"),c(6,4)),
                       Time = c(rep("Pre.Treat",3),rep("After.Treat1",1),rep("After.Treat2",2),
                                rep("Pre.Treat",1),rep("After.Treat1",1),rep("After.Treat2",2)),


                       Location = c("chr1_10","chr1_20","chr2_30",
                                    "chr1_10",
                                    "chr1_10","chr1_20",

                                    "chr1_10",
                                    "chr1_10",
                                    "chr1_10","chr4_40"),

                       alt = c("ACT","ACGTCG","AGG",
                               "ACT",
                               "ACT","ACGTCG",

                               "A",
                               "A",
                               "A","C"),

                       ref = c("A","A","A",
                               "A",
                               "A","A",

                               "G",
                               "G",
                               "G","T"),

                       ref_depth = c(12,11,40,
                                     9,
                                     13,14,

                                     20,
                                     25,
                                     23,15),


                       VAF = c(0.55,0.3,0.2,
                               0.48,
                               0.45,0.2,

                               0.3,
                               0.55,
                               0.45,0.1),

                       SYMBOL = c("Gene1","Gene2","Gene3",
                                  "Gene1",
                                  "Gene1","Gene2",

                                  "Gene1",
                                  "Gene1",
                                  "Gene1","Gene4"),

                       Consequence = c("Nonsyn","Damaging","Damaging",
                                       "Nonsyn",
                                       "Nonsyn","Damaging",

                                       "Nonsyn",
                                       "Nonsyn",
                                       "Nonsyn","Damaging")) %>%

  unite(mutation_key, SYMBOL,Location,ref,alt,sep="-",remove=FALSE) %>%
  unite(mutation_det, SYMBOL,Consequence,sep="-",remove=FALSE) %>%
  separate(Location,into = c("chrom","pos"),sep="_",remove=TRUE) %>%
  unite(SampleName,PID,Time,sep=".",remove = FALSE) %>%
  mutate(alt_depth = round(ref_depth * VAF),
         Time = factor(Time,levels = c("Pre.Treat","After.Treat1","After.Treat2")),
         variant_type = case_when(nchar(as.character(alt)) > 1 | nchar(as.character(ref)) > 1 ~ "INDEL",
                                  TRUE ~ "SNV"))


# example clinical
example_metadata <- data.frame(PID = c(rep("P1",3),rep("P2",3)),
                       Time = c("Pre.Treat","After.Treat1","After.Treat2",
                                "Pre.Treat","After.Treat1","After.Treat2"),
                       Blast = c(0.78,0.6,0.7,
                                 0.90,0.05,0.2),
                       Outcome = c(rep("Refractory",3),rep("Relapse",3))) %>%
  unite(SampleName,PID,Time,sep=".",remove = FALSE) %>%
  mutate(Time = factor(Time,levels = c("Pre.Treat","After.Treat1","After.Treat2")))


save(example_variants,example_metadata,file = "./data/example_data.rda")

