library(ggplot2)
library(dplyr)

decodeDNM_raw <- read.table("tables/decode_DNMs.tsv", header = TRUE)
decodeDNM <- cbind(aggregate(decodeDNM_raw[c("Fathers_age_at_conception", "Mothers_age_at_conception")], by = list(decodeDNM_raw$Proband_nr), FUN = mean),
                   mutations = aggregate(decodeDNM_raw, by = list(decodeDNM_raw$Proband_nr), FUN = length)$Fathers_age_at_conception)
human_haploid <- 2.63e9
decodeDNM$mu_rate <- decodeDNM$mutations / (2 * human_haploid)

DOB_table <- read.table("tables/DOB.txt")
names(DOB_table) <- c("Individual", "Date")
DOB_table$Date <- as.Date(DOB_table$Date, format = '%m/%d/%Y')
DOB_vector <- DOB_table$Date
names(DOB_vector) <- DOB_table$Individual

hapsize_table <- read.table("tables/samtools_depth.table")[,seq(2,36,2)]
auto_size <- rowSums(hapsize_table)
names(auto_size) <- paste("trio", 1:11, sep = "")

denovo_candidates <- read.table("denovo_table.txt", header = TRUE)

hetcallability_table <- read.table("tables/het_callabilityGQ70.table")
names(hetcallability_table) <- c("denom", "transmit", "dpgq", "bam", "ab30", "ab35", "ab40", "ab45", "ab50", "ab55", "ab60")

homcallability_table <- read.table("tables/hom_callabilityGQ70.table")
names(homcallability_table) <- c("denom", "transmit", "dpgq", "bamad1", "bamad0")

callability_table <- hetcallability_table[["ab35"]] / hetcallability_table[["transmit"]] *
  (homcallability_table[["bamad0"]] / homcallability_table[["transmit"]])**2

mutation_counts <- aggregate(denovo_candidates[,1], by = list(denovo_candidates$PROBAND), FUN = length)
names(mutation_counts) <- c("proband", "mutations")

mutation_counts$F_phase <- aggregate(denovo_candidates$COMBINED_PHASE == "F", by = list(denovo_candidates$PROBAND), FUN = function(x) { sum(x, na.rm = TRUE)} )$x
mutation_counts$M_phase <- aggregate(denovo_candidates$COMBINED_PHASE == "M", by = list(denovo_candidates$PROBAND), FUN = function(x) { sum(x, na.rm = TRUE)} )$x

mutation_counts$rbF_phase <- aggregate(denovo_candidates$READBASED_PHASE == "F", by = list(denovo_candidates$PROBAND), FUN = function(x) { sum(x, na.rm = TRUE)} )$x
mutation_counts$rbM_phase <- aggregate(denovo_candidates$READBASED_PHASE == "M", by = list(denovo_candidates$PROBAND), FUN = function(x) { sum(x, na.rm = TRUE)} )$x

mutation_counts[["trio"]] <- c("trio11", paste("trio", 1:10, sep = ""))
mutation_counts <- mutation_counts[order(as.numeric(substring(mutation_counts$trio, 5))),]

mutation_counts[["callability"]] <- callability_table 
mutation_counts[["pile_size"]] <- auto_size 

mutation_counts[["callable_size"]] <- mutation_counts$callability * mutation_counts$pile_size
mutation_counts[["rate"]] <- mutation_counts$mutations / (2 * mutation_counts$callable_size)

trio_table <- read.csv("tables/trio_tableMFC.csv", header = FALSE)
names(trio_table) <- c("Mother", "Father", "Child")
trio_table$trio <- paste("trio", 1:11, sep = "")

DOB_idx <- match(mutation_counts$trio, trio_table$trio)
mutation_counts$GTMale <- DOB_vector[as.character(trio_table[DOB_idx,]$Child)] - DOB_vector[as.character(trio_table[DOB_idx,]$Father)]
mutation_counts$GTFemale <- DOB_vector[as.character(trio_table[DOB_idx,]$Child)] - DOB_vector[as.character(trio_table[DOB_idx,]$Mother)]
mutation_counts$GTMale <- as.numeric(mutation_counts$GTMale / 365)
mutation_counts$GTFemale <- as.numeric(mutation_counts$GTFemale / 365)

## Writing positions for samtools faidx to check for CpG site
#
# dnc_pos <- denovo_candidates[,1:2]
# dnc_pos <- data.frame(regions = paste(dnc_pos[["CHROM"]], ":",
#                                       dnc_pos[["POS"]]-1, "-",
#                                       dnc_pos[["POS"]]+1, sep = ""))
# write.table(dnc_pos, file = "dncGQ70_pos.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# $ samtools faidx ref/GCA_000181335.4_Felis_catus_9.0_genomic.fna -r dncGQ70_pos.txt > dncGQ70_seq.txt

dnc_seq <- read.table("tables/dncGQ70_seq.txt", comment.char = ">")
denovo_candidates$triplet <- as.vector(sapply(dnc_seq, toupper))
denovo_candidates$CpG <- substr(denovo_candidates$triplet, 1, 2) == "CG" | substr(denovo_candidates$triplet, 2, 3) == "CG"

collapse_spectra <- function(mtable) {
  tmp <- mtable[which(mtable > 0)]
  
  tmp["A>C"] <- tmp["A>C"] + tmp["T>G"]
  tmp["A>G"] <- tmp["A>G"] + tmp["T>C"]
  tmp["A>T"] <- tmp["A>T"] + tmp["T>A"]
  tmp["C>A"] <- tmp["C>A"] + tmp["G>T"]
  tmp["C>G"] <- tmp["C>G"] + tmp["G>C"]
  tmp["C>T"] <- tmp["C>T"] + tmp["G>A"]
  
  tmp[which(is.na(tmp))] <- 0  
  return(tmp[c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T")])
}

mut_table <- table(paste(denovo_candidates$REF, denovo_candidates$ALT, sep = ">")) %>%
  collapse_spectra
cpg_table <- table(paste(denovo_candidates$REF, denovo_candidates$ALT, sep = ">")[denovo_candidates$CpG]) %>%
  collapse_spectra

human_muttable <- table(paste(decodeDNM_raw$Ref, decodeDNM_raw$Alt, sep = ">")) %>% collapse_spectra

cat_spectrum <- as.numeric(mut_table) / sum(mut_table)
human_spectrum <- as.numeric(human_muttable) / sum(human_muttable)

## Preparing for plots on analysis of spectrum

decode_spectral <- decodeDNM_raw
decode_spectral$mut_class <- paste(decode_spectral$Ref, decode_spectral$Alt, sep = "_")
# Remove indels
decode_spectral <- decode_spectral[!(sapply(decode_spectral$mut_class, nchar) > 3),]
phased_decode <- subset(decode_spectral, !is.na(Phase_combined))

mut_classes <- c("A_C", "A_G", "A_T", "C_A", "C_G", "C_T")
complementary_classes <- c("T_G", "T_C", "T_A", "G_T", "G_C", "G_A")
class_map <- setNames(mut_classes, complementary_classes)

collapse_class <- function(mut_class) {
  if(mut_class %in% mut_classes)
    return(mut_class)
  return(class_map[mut_class])
}

phased_decode$mut_class <- sapply(phased_decode$mut_class, collapse_class)

cref_positions <- subset(phased_decode, mut_class == "C_T")[,1:2]

cref_triplet <- data.frame(regions = paste(cref_positions[["Chr"]], ":",
                                           cref_positions[["Pos_hg38"]]-1, "-",
                                           cref_positions[["Pos_hg38"]]+1, sep = ""))
# write.table(cref_triplet, file = "decode_Ctriplet_pos.txt", quote = FALSE, sep ="\t", row.names = FALSE, col.names = FALSE)
# $ samtools faidx GCA_000001405.15_GRCh38_genomic.fna -r decode_Ctriplet_pos.txt > decode_Ctriplet_seq.txt
cref_sequences <- read.table("tables/decode_Ctriplet_seq.txt", comment.char = ">")
phased_decode$cref_triplet <- NA
phased_decode[phased_decode$mut_class == "C_T",]$cref_triplet <- as.vector(sapply(cref_sequences, toupper))

phased_decode$CpG <- NA
phased_decode$CpG <- substr(phased_decode$cref_triplet, 1, 2) == "CG" | substr(phased_decode$cref_triplet, 2, 3) == "CG"

father_phased <- subset(phased_decode, Phase_combined == "father")
mother_phased <- subset(phased_decode, Phase_combined == "mother")

aggregate_spectra <- function(phased_mutations) {
  mutations_byproband <- aggregate(phased_mutations$mut_class, by = list(Proband = phased_mutations$Proband_nr), length)  
  ages_byproband <- aggregate(phased_mutations[,c("Fathers_age_at_conception", "Mothers_age_at_conception")],
                              by = list(Proband = phased_mutations$Proband_nr), unique)
  mutations_byproband[,c("Fathers_age_at_conception", "Mothers_age_at_conception")] <- ages_byproband[,c(2,3)]
  
  for(class_i in mut_classes) {
    class_byproband <- aggregate(phased_mutations$mut_class, by = list(Proband = phased_mutations$Proband_nr),
                                 FUN = function(classes) { return(sum(classes == class_i)) })
    names(class_byproband)[2] <- class_i
    
    mutations_byproband <- cbind(mutations_byproband, class_byproband[2])
  }
  
  return(mutations_byproband)
}

paternal_spectrum <- aggregate_spectra(father_phased)
maternal_spectrum <- aggregate_spectra(mother_phased)


paternal_models <- list()
for(class in mut_classes) {
  paternal_models[[class]] <- glm(as.formula(paste(class, "~ Fathers_age_at_conception")),
                                  data = paternal_spectrum, family = poisson(link = "identity"))
}
maternal_models <- list()
for(class in mut_classes) {
  maternal_models[[class]] <- glm(as.formula(paste(class, "~ Mothers_age_at_conception")),
                                  data = maternal_spectrum, family = poisson(link = "identity"))
}

predict_spectrum <- function(paternal_age, maternal_age) {
  mutations <- sapply(paternal_models, predict, data.frame(Fathers_age_at_conception = paternal_age)) +
    sapply(maternal_models, predict, data.frame(Mothers_age_at_conception = maternal_age))
  spectrum <- mutations / sum(mutations)
  return(setNames(spectrum, mut_classes))
}
