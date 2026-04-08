library(tidyverse)
library(phyloseq)
library(microViz)

# data ophalen voor bacteriën
ps <- read_rds("../S4B/S4B_16S.rds")

# Metadata toevoegen
Metadata <- readxl::read_excel("DNA extractie_codering.xlsx")

metadf <- data.frame(Metadata) %>%
  select(Gewas, Ras, Bodemsoort, Teelt)
row.names(metadf) <- Metadata$`Nr epje`

# Metadata omzetten naar ps-format
metaobj <- sample_data(metadf)

# Metadata toevoegen aan ps object
sample_data(ps) <- metaobj




ps <- ps %>%
  phyloseq::subset_taxa(Order != "Chloroplast") %>%
  phyloseq::subset_taxa(Family != "Mitochondria") %>%
  phyloseq::subset_taxa(Kingdom == "Bacteria") %>%
  phyloseq::subset_taxa(Phylum != "Unclassified ")

# fix NA's
ps <- tax_fix(ps)

# check ps
phyloseq_validate(ps)

# plaatje
ps %>% 
  ps_filter(Gewas == "Broccoli") %>% 
  ps_filter(Bodemsoort == "Klei") %>% 
  comp_barplot("Class", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_grid(rows = vars(Ras), cols = vars(Teelt), scales = "free") +
  coord_flip()

# plaatje
ps %>% 
  ps_filter(Gewas == "Prei") %>% 
  comp_barplot("Class", n_taxa = 15, merge_other = FALSE, label = NULL) +
  facet_grid(rows = vars(Ras), cols = vars(Teelt), scales = "free") +
  coord_flip()


# multivariate
# perform ordination
unconstrained_aitchison_pca <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>% 
  ps_filter(Gewas == "Prei" & Teelt != "Geen middelen") %>% 
  tax_agg("Family") %>%
  tax_transform("clr") %>%
  ord_calc()
# ord_calc will automatically infer you want a "PCA" here
# specify explicitly with method = "PCA", or you can pick another method

# create plot
pca_plot <- unconstrained_aitchison_pca %>%
  ord_plot(
    plot_taxa = 1:6, colour = "Ras", size = 1.5,
    tax_vec_length = 0.325,
    tax_lab_style = tax_lab_style(max_angle = 90, aspect_ratio = 1),
    auto_caption = 8
  )

# customise plot
customised_plot <- pca_plot +
  stat_ellipse(aes(linetype = Ras, colour = Ras), linewidth = 0.3) + # linewidth not size, since ggplot 3.4.0
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom") +
  coord_fixed(ratio = 1, clip = "off") # makes rotated labels align correctly

# show plot
customised_plot
