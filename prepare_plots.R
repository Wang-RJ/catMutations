library(ggplot2)
library(dplyr)

decodeDNM_raw <- read.table("tables/decode_DNMs.tsv", header = TRUE)
decodeDNM <- cbind(aggregate(decodeDNM_raw[c("Fathers_age_at_conception", "Mothers_age_at_conception")], by = list(decodeDNM_raw$Proband_nr), FUN = mean),
                   mutations = aggregate(decodeDNM_raw, by = list(decodeDNM_raw$Proband_nr), FUN = length)$Fathers_age_at_conception)

DOB_table <- read.table("DOB.txt")
names(DOB_table) <- c("Individual", "Date")
DOB_table$Date <- as.Date(DOB_table$Date, format = '%m/%d/%Y')
DOB_vector <- DOB_table$Date
names(DOB_vector) <- DOB_table$Individual

construct_mframe <- function(dnc_frame, callability_bytrio) {
  mu_frame <- aggregate(dnc_frame[,1], by = list(dnc_frame$trio), FUN = length)
  names(mu_frame) <- c("trio", "mutations")
  mu_frame <- mu_frame[order(as.numeric(substring(mu_frame$trio, 5))),]
  
  mu_frame[["callability"]] <- callability_bytrio
  mu_frame[["pile_size"]] <- auto_size
  mu_frame[["callable_size"]] <- callability_bytrio * auto_size
  
  mu_frame[["rate"]] <- mu_frame$mutations / (2 * mu_frame$callable_size)
  mu_frame[["rate_SN"]] <- mu_frame$rate * 1e8
  
  DOB_idx <- match(mu_frame$trio, trio_table$trio)
  mu_frame$GTMale <- DOB_vector[as.character(trio_table[DOB_idx,]$Child)] - DOB_vector[as.character(trio_table[DOB_idx,]$Father)]
  mu_frame$GTFemale <- DOB_vector[as.character(trio_table[DOB_idx,]$Child)] - DOB_vector[as.character(trio_table[DOB_idx,]$Mother)]
  mu_frame$GTMale <- as.numeric(mu_frame$GTMale / 365)
  mu_frame$GTFemale <- as.numeric(mu_frame$GTFemale / 365)
  
  return(mu_frame)
}

phase_list <- lapply(split(factor(denovo_candidates$rb_phase), denovo_candidates$trio), table)
phase_frame <- as.data.frame(do.call(rbind, phase_list))
phase_frame$trio <- rownames(phase_frame)
GATK_mucount$F_phase <- NA
GATK_mucount$M_phase <- NA
GATK_mucount[match(phase_frame$trio, GATK_mucount$trio),][c("F_phase", "M_phase")] <- phase_frame[,1:2]

## Writing positions for samtools faidx to check for CpG site
#
# dnc_pos <- denovo_candidates[,1:2]
# dnc_pos <- data.frame(regions = paste(dnc_pos[["CHROM"]], ":",
#                                       dnc_pos[["POS"]]-1, "-",
#                                       dnc_pos[["POS"]]+1, sep = ""))
# write.table(dnc_pos, file = "dncGQ70_pos.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# $ samtools faidx ref/GCA_000181335.4_Felis_catus_9.0_genomic.fna -r dncGQ70_pos.txt > dncGQ70_seq.txt

dnc_seq <- read.table("dncGQ70_seq.txt", comment.char = ">")
denovo_candidates$triplet <- as.vector(sapply(dnc_seq, toupper))
denovo_candidates$CpG <- substr(denovo_candidates$triplet, 1, 2) == "CG" | substr(denovo_candidates$triplet, 2, 3) == "CG"

cpg_mutable <- sapply(split(denovo_candidates$CpG, denovo_candidates$trio), table)
cpg_mucount <- sapply(cpg_mutable, sum) - sapply(cpg_mutable, "[[", "FALSE")

GATK_mucount$CpG_mutations <- NA
GATK_mucount[match(names(cpg_mucount), GATK_mucount$trio),]$CpG_mutations <- cpg_mucount
GATK_mucount$CpG_rate <- GATK_mucount$CpG_mutations / (2 * GATK_mucount$callable_size * 27735107 / 2.675e9)

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
