library(tidyverse)
theme_set(theme_classic(base_size = 16))
library(phyloseq)
library(microViz)
library(emmeans)

# data ophalen voor schimmels
ps <- read_rds("../S4B/S4B_ITS.rds")
sample_names(ps)

sample_names(ps) <- str_match(
  sample_names(ps), "i5.(.*)-ITS")[, 2]

# Metadata toevoegen
Metadata <- readxl::read_excel("DNA extractie_codering.xlsx") %>% 
  filter(`Nr epje` %in% sample_names(ps))

metadf <- data.frame(Metadata) %>%
  select(Gewas, Ras, Bodemsoort, Teelt)
row.names(metadf) <- Metadata$`Nr epje`

# Metadata omzetten naar ps-format
metaobj <- sample_data(metadf)

# Metadata toevoegen aan ps object
sample_data(ps) <- metaobj

write_rds(ps, "S4B_ITS.rds")


# NB kost teveel tijd

# plot_bar(ps, "Ras", "Abundance", fill = "Family", facet_grid = "gewas")



# Alleen broccoli
ps_b <- ps %>% ps_filter(Gewas == "Broccoli")

# check ps
phyloseq_validate(ps)
# fix NA's
ps <- tax_fix(ps)

# plaatje
ps %>% 
  ps_filter(Gewas == "Broccoli") %>% 
  ps_filter(Bodemsoort == "Klei") %>% 
  comp_barplot("Family", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_grid(rows = vars(Ras), cols = vars(Teelt), scales = "free") +
  coord_flip()

# prei en mycorrhiza's
ps_prei <- ps %>% 
  tax_fix %>% 
  ps_filter(Gewas == "Prei") %>%
  tax_select(tax_list = "c_Glomeromycetes")

df_prei_myc <- psmelt(ps_prei) %>% 
  group_by(Sample, Ras, Teelt) %>% 
  summarise(Abundance = sum(Abundance))

# grafiek met gemiddelde som van abundantie VAM-schimmels
df_prei_myc %>% 
  filter(Teelt != "Biologisch") %>% 
  ggplot(aes(Ras, Abundance, fill = Teelt)) +
  stat_summary(geom = "bar", position = position_dodge(0.9)) +
  stat_summary(geom = "errorbar", position = position_dodge(0.9), width = 0.5) +
  coord_flip()
# bijbehorende statistiek
fit <- glm(Abundance ~ Ras * Teelt, data = df_prei_myc %>% filter(Ras != "Bodem"), family = poisson())
anova(fit)
r <- rstudent(fit)
hist(r) # check normaliteit
shapiro.test(r)
e1071::skewness(r, type = 2)
e1071::kurtosis(r, type = 2)

em <- emmeans(fit, specs = ~ Ras)
em
plot(em)

ps_prei <- prune_samples(sample_sums(ps_prei)>0, ps_prei)

ps_prei %>%
  comp_barplot("Species", label = NULL) +
  facet_grid(rows = vars(Ras), cols = vars(Teelt), scales = "free") +
  coord_flip()
  
ps_prei %>%
  comp_barplot("Class", label = NULL, tax_transform_for_plot = "identity", bar_width = 1) +
  facet_grid(rows = vars(Ras), cols = vars(Teelt), scales = "fixed") +
  coord_flip()

# multivariate
# perform ordination
ps_nmds <- ps %>%
  ps_filter(Gewas == "Prei",
            # Teelt == "Geen middelen"
            ) %>% 
  tax_agg("Class") %>%
  # tax_transform("clr") %>%
  dist_calc() %>% 
  ord_calc(method = 'NMDS')
# ord_calc will automatically infer you want a "PCA" here
# specify explicitly with method = "PCA", or you can pick another method

# create plot
nmds_plot <- ps_nmds %>%
  ord_plot(
    plot_taxa = 1:6, colour = "Teelt", shape = "Ras", size = 1.5,
    tax_vec_length = 0.325,
    tax_lab_style = tax_lab_style(max_angle = 90, aspect_ratio = 1),
    auto_caption = 8
  )
nmds_plot

# Adonis

adonis_prei <- microViz::dist_permanova(ps_nmds, 
                         variables = c("Teelt", "Ras"), 
                         interactions = "Teelt*Ras", 
                         n_perms = 999,
                         by = "terms")

ps_prei <- subset_samples(ps, Gewas == "Prei")

ps.prop_subset <- transform_sample_counts(ps_prei, function(otu) otu/sum(otu))

bray_dist <- phyloseq::distance(ps.prop_subset, method = "bray")

vegan::adonis2(bray_dist ~ Teelt * Ras, 
               # data = as.data.frame(sample_data(ps.prop_subset)),
               data = test2,
               permutations = 999,
               method = "bray",
               by = "terms")

pairwise_adonis_prei <- pairwiseAdonis::pairwise.adonis2(
  bray_dist ~ Teelt *  Ras, 
  data = test2, 
  p.adjust.m = "bonferroni", 
  perm = 99)

pairwise_adonis_prei <- pairwiseAdonis::pairwise.adonis2(bray_dist ~ sample_data(ps_prei)$Teelt  data = ps_prei, p.adjust.m = "bonferroni", perm = 999)
pairwise_adonis_prei


pairwise_adonis_CycTreat <- vegan::pairwise.adonis2(bray_dist ~ CycTreat, data = sample_info_tab_subset, p.adjust.m = "bonferroni", perm = 999)
pairwise_adonis_CycTreat
write.table(pairwise_adonis_CycTreat, "pairwise_adonis_CycTreat.txt", sep = "\t", quote = FALSE, row. Names = FALSE)

# # customise plot
# customised_plot <- nmds_plot +
#   stat_ellipse(aes(colour = Ras), linewidth = 0.3) + # linewidth not size, since ggplot 3.4.0
#   scale_colour_brewer(palette = "Set1") +
#   theme(legend.position = "bottom") +
#   coord_fixed(ratio = 1, clip = "off") # makes rotated labels align correctly
# 
# # show plot
# customised_plot
# 

