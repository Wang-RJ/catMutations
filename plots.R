library(binom)
source("prepare_plots.R")

col1 <- "#d95f02"
col2 <- "#7570b3"
col3 <- "#e7298a"

ggplot(mutation_counts, aes(x = GTMale - 0.5, y = rate)) +
  geom_smooth(fullrange = TRUE, method = lm, color = col1, fill = col1) +
  geom_point(bg = col1, size = 4, shape = 21) + theme_classic() +
  theme(text = element_text(size = 18)) + ylab("mutation rate per generation") + xlab("paternal age")

# Figure S1
ggplot(data = mutation_counts) + labs(x = "Parental age", y = "Phased mutation count") +
  geom_smooth(aes(x = GTMale, y = rbM_phase), method = "glm",
              method.args = list(family = poisson(link = "identity")),
              bg = col1, col = col1, cex = 0.7) + 
  geom_point(aes(x = GTMale, y = rbM_phase), bg = col1, size = 4, shape = 21) +
  geom_smooth(aes(x = GTFemale, y = rbF_phase), method = "glm",
              method.args = list(family = poisson(link = "identity")),
              bg = col2, col = col2, cex = 0.7) +
  geom_point(aes(x = GTFemale, y = rbF_phase), bg = col2, size = 4, shape = 21) +
  theme_classic() + coord_cartesian(ylim = c(0,20), xlim = c(0,13))

# Figure S2
mubar <- barplot(t(matrix(c(mut_table - cpg_table, cpg_table), ncol = 2)) / sum(mut_table), names.arg = names(cpg_table),
                 legend.text = c("non-CpG", "CpG"), args.legend = list(x = "topleft"), ylim = c(0,0.55))

binomCI <- binom.confint(as.vector(mut_table), sum(mut_table), methods = "wilson")
segments(mubar, binomCI$lower, mubar, binomCI$upper)

phased_candidates <- subset(denovo_candidates, !is.na(COMBINED_PHASE))
phased_mut_table <- table(paste(phased_candidates$REF, phased_candidates$ALT, sep = ">")) %>% collapse_spectra

phased_male_table <- table(paste(phased_candidates$REF, phased_candidates$ALT, sep = ">")[phased_candidates$COMBINED_PHASE == "M"]) %>% collapse_spectra
phased_female_table <- table(paste(phased_candidates$REF, phased_candidates$ALT, sep = ">")[phased_candidates$COMBINED_PHASE == "F"]) %>% collapse_spectra

mubar <- barplot(t(matrix(c(phased_male_table, phased_female_table), ncol = 2) / sum(phased_mut_table)), names.arg = names(cpg_table),
                 legend.text = c("M", "F"), args.legend = list(x = "topleft"), ylim = c(0,0.55))
binomCI <- binom.confint(as.vector(phased_mut_table), sum(phased_mut_table), methods = "wilson")
segments(mubar, binomCI$lower, mubar, binomCI$upper)

