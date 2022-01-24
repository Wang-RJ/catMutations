col1 <- "#d95f02"
col2 <- "#7570b3"
col3 <- "#e7298a"

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
