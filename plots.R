library(binom)
source("prepare_plots.R")

col1 <- "#d95f02"
col2 <- "#7570b3"
col3 <- "#e7298a"

# Figure 2
# Some calculations for Poisson regression with proper adjustment for callable genome size
# 
mutation_counts$diploid_callable <- mutation_counts$callable_size * 2
mutation_counts$Z_male <- mutation_counts$diploid_callable * mutation_counts$GTMale

cat_pois_model <- glm(mutations ~ Z_male + diploid_callable + 0,
                      family = poisson(link = "identity"), data = mutation_counts)
cat_pois_fit <- data.frame(GTMale = seq(0,15,1),
                           predict(cat_pois_model, data.frame(Z_male = seq(0,15,1) * 2 * cat_haploid,
                                                              diploid_callable = 2 * cat_haploid), se = TRUE))

ggplot(mutation_counts, aes(x = GTMale, y = rate)) + geom_point(bg = col1, size = 4, shape = 21) + 
  geom_line(aes(x = GTMale, y = fit / (2 * cat_haploid)), data = cat_pois_fit, color = col1) +
  geom_ribbon(aes(x = GTMale, y = fit / (2 * cat_haploid),
                  ymin = fit / (2 * cat_haploid) - 1.96 * se.fit / (2 * cat_haploid),
                  ymax = fit / (2 * cat_haploid) + 1.96 * se.fit / (2 * cat_haploid)),
              data = cat_pois_fit, fill = col1, alpha = 0.3) +
  geom_smooth(aes(x = Fathers_age_at_conception - 13, y = mutations / (2 * human_haploid)),
              data = decodeDNM, method = glm, method.args = list(family = poisson(link = "identity")),
              fullrange = TRUE, se = FALSE, cex = 1, lty = 1) +
  geom_smooth(aes(x = Fathers_age_at_conception, y = mutations / (2 * human_haploid)),
              data = decodeDNM, method = glm, method.args = list(family = poisson(link = "identity")),
              fullrange = TRUE, se = FALSE, cex = 1, lty = 2) +
  coord_cartesian(xlim = c(0, 13), ylim = c(0,1.75) * 1e-8) +
  scale_x_continuous(limits = c(0,70), breaks = seq(0,30,5)) +
  scale_y_continuous(limits = c(0,10) * 1e-8, breaks = seq(0,2,0.25) * 1e-8) +
  theme_classic() + theme(text = element_text(size = 18)) +
  ylab("mutation rate") + xlab("Paternal age")

# Figure 3
predicted_spectra_0to50 <- sapply(seq(0,50,1), function(mean_age) { return(predict_spectrum(mean_age, mean_age)) })
humanpredict_ggframe <- data.frame(frequency = as.vector(predicted_spectra_0to50 - human_spectrum),
                                   class = rep(rownames(predicted_spectra_0to50), ncol(predicted_spectra_0to50)),
                                   age = rep(seq(0,50), each = nrow(predicted_spectra_0to50)))

ggplot(humanpredict_ggframe, aes(x = age, y = frequency, col = class)) + geom_line(size = 0.75) + theme_bw() +
  scale_color_brewer(type = 'qual', palette = 'Dark2') + geom_vline(xintercept = 16.8) + geom_vline(xintercept = 3.8)

cathuman_ggframe <- rbind(data.frame(frequency = predict_spectrum(4.68, 2.92), class = mut_classes, group = "total longevity"),
                          data.frame(frequency = cat_spectrum, class = mut_classes, group = "cat"),
                          data.frame(frequency = predict_spectrum(17.68, 15.92), class = mut_classes, group = "reproductive longevity"))

ggplot(cathuman_ggframe, aes(x = class, y = frequency, fill = factor(group, levels = c("total longevity", "cat", "reproductive longevity")))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_bw() + scale_fill_manual(name = "group", values = c("#67a9cf", "#4d4d4d", "#ef8a62"))

# Figure S1
# read-based only for this plot because transmission based will add phased mutations to a few specific trios
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
