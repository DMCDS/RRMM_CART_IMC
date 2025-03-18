# Save plots function with output directory
save_plot <- function(plot_list, output_dir, width = 5, height = 6) {
  for (name in names(plot_list)) {
    plot <- plot_list[[name]]
    
    # Check if the object is a ggplot
    if (inherits(plot, "ggplot")) {
      ggsave(filename = file.path(output_dir, paste0(name, ".svg")),
             plot = plot, width = width, height = height)
    } else {
      warning(paste("Object", name, "is not a ggplot. Skipping."))
    }
  }
}

## Loading raw counts matrix
emd_raw <- read.csv("full_results/emd/emd_raw.csv")

## Modifying variables to the correct classes and renaming cell types
str(emd_raw)
emd_raw$Cell_ID <- as.character(emd_raw$Cell_ID)
emd_raw$patient_n <- as.character(emd_raw$patient_n)


emd_raw <- emd_raw %>%
  mutate(major_cell_type = str_replace(major_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         major_cell_type = str_replace(major_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"))

emd_raw <- emd_raw %>%
  mutate(medium_cell_type = str_replace(medium_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         medium_cell_type = str_replace(medium_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"),
         medium_cell_type = str_replace(medium_cell_type, "MDSCs", "Myeloid/MDSCs"))

emd_raw <- emd_raw %>%
  mutate(minor_cell_type = str_replace(minor_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         minor_cell_type = str_replace(minor_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"),
         minor_cell_type = str_replace(minor_cell_type, "MDSCs", "Myeloid/MDSCs"))

emd_raw <- emd_raw %>%
  mutate(functional_minor_cell_type = str_replace(functional_minor_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         functional_minor_cell_type = str_replace(functional_minor_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"),
         functional_minor_cell_type = str_replace(functional_minor_cell_type, "MDSCs", "Myeloid/MDSCs"))

emd_maj_labels <- unique(emd_raw$major_cell_type)
emd_med_labels <- unique(emd_raw$medium_cell_type)
emd_min_labels <- unique(emd_raw$minor_cell_type)
emd_func_labels <- unique(emd_raw$functional_minor_cell_type)

emd_raw <- emd_raw %>%
  mutate(sample_id_alt = case_when(
    sample_id == "15K30184" ~ "PMD1",
    sample_id == "22A64389" ~ "EMD3",
    sample_id == "22K77729" ~ "PMD2",
    sample_id == "23K58040" ~ "EMD1",
    sample_id == "22T25546" ~ "EMD2",
    TRUE ~ "Unknown" # Default if no match is found
  ))

emd_raw <- emd_raw %>%
  mutate(functional_minor_cell_type = case_when(
    functional_minor_cell_type == "PDL1+ Ki-67+ Myeloma" ~ "Ki-67+/PD-L1+ Myeloma",
    functional_minor_cell_type == "Ki-67+ Myeloma" ~ "Ki-67+/PD-L1- Myeloma",
    functional_minor_cell_type == "Myeloma" ~ "Ki-67-/PD-L1- Myeloma",
    TRUE ~ functional_minor_cell_type # Default if no match is found
  ))

emd_raw <- emd_raw |>
  mutate(sample_id_alt = factor(sample_id_alt, levels = c("EMD1", "EMD2", "EMD3", "PMD1", "PMD2")))


emd_maj_labels <- unique(emd_raw$major_cell_type)
emd_med_labels <- unique(emd_raw$medium_cell_type)
emd_min_labels <- unique(emd_raw$minor_cell_type)
emd_func_labels <- unique(emd_raw$functional_minor_cell_type)

## Changed color schemes to accomodate new cell types
emd_palette_medium <- c(
  "B cells" = '#8B4513',
  "CD4 Tmem" = '#9400D3',
  "CD8 GZMB+ Tmem" = "#c3c3c3",
  "CD8 Tex" = "#0A0A0A",
  "CD8 Tmem" = "#969696",
  "ECs" = "#ff7f0e",
  "Fibroblasts" = '#FFFDD0',
  "M1-like macrophage" = "#006400",
  "M2-like macrophage" = "#2ca02c",
  "Myeloid/MDSCs" = "#FFD700",
  "MoDC" = '#00CED1',
  "Monocytes" = "#1f77b4",
  "Myeloma" = '#FF0000',
  "NK cells" = '#FF007F',
  "Treg" = '#D8BFD8'
)

color_palette_minor <- c(
  "Classical\nmonocytes" = "#1f77b4",
  "Nonclassical\nmonocytes" = "#87b4d4",
  "Intermediate\nmonocytes" = '#1E90FF',
  "ECs" = "#ff7f0e",
  "M2-like\nmacrophage" = "#2ca02c",
  "M1-like\nmacrophage" = "#006400",
  "Myeloid/MDSCs" = "#FFD700",
  "Non Th17 CD4 Tmem" = '#9400D3',
  "Th17 CD4 Tmem" ='#4B0082',
  "CD4 Tnaive" = '#DA70D6',
  "CD4 Tex" = '#DDA0DD' ,
  "Treg" = '#D8BFD8',
  "B cells" = '#8B4513',
  "MBCs" = '#A0522D',
  "CD56 bright NK" = '#FF007F',
  "CD56 dim NK" = '#FF69B4',
  "CD8 Tnaive" = "#7f7f7f",
  "CD8 Tmem" = "#969696",
  "CD8 GZMB Tmem" = "#c3c3c3",
  "CD8 Tex" = "#0A0A0A",
  "Myeloma" = '#FF0000',
  "MKs" = "#17becf",
  "MoDC" = '#00CED1',
  "Adipocyte" = '#FFFDD0'
)

color_palette_2 <- c(
  "All" = "black",
  "Classical monocytes" = "#1f77b4",
  "Nonclassical monocytes" = "#87b4d4",
  "Intermediate monocytes" = '#1E90FF',
  "ECs" = "#ff7f0e",
  "M2-like macrophage" = "#2ca02c",
  "M1-like macrophage" = "#006400",
  "Myeloid/MDSCs" = "#d62728",
  "pSTAT1neg Non Th17 CD4 Tmem" = '#9400D3',
  "pSTAT1pos Non Th17 CD4 Tmem" = '#8A2BE2' ,
  "pSTAT1neg Th17 CD4 Tmem" ='#4B0082',
  "pSTAT1pos Th17 CD4 Tmem" = '#6A0DAD' ,
  "CD4 Tnaive" = '#DA70D6',
  "CD4 Tex" = '#DDA0DD' ,
  "Treg" = '#D8BFD8',
  "B cells" = '#8B4513',
  "MBCs" = '#A0522D',
  "pSTAT1neg CD56 dim NK" = '#FF007F',
  "pSTAT1neg CD56 bright NK" = '#FF6EB4',
  "pSTAT1pos CD56 dim NK" = '#FF69B4',
  "pSTAT1pos CD56 bright NK" = '#FF1493',
  "CD8 Tnaive" = "#7f7f7f",
  "pSTAT1pos CD8 Tmem" = "#969696",
  "pSTAT1neg CD8 Tmem" = "#acacac",
  "pSTAT1pos CD8 GZMBpos Tmem" = "#c3c3c3",
  "pSTAT1neg CD8 GZMBpos Tmem" = "#dadada",
  "CD8 Tex" = "#0A0A0A",
  "Myeloma" = '#FF0000',
  "Ki-67pos Myeloma" = '#DC143C',
  "PD-L1pos Myeloma" = '#B22222',
  "MKs" = "#17becf",
  "MoDC" = '#00CED1',
  "Adipocyte" = '#FFFDD0'
)

# Specify the output directory
output_dir <- "emd_data"
dir.create(output_dir, showWarnings = FALSE)


# Displaying cell types per sample
emd_raw |>
  group_by(major_cell_type)|>
  count(major_cell_type)|>
  ungroup()|>
  mutate(perc = n/sum(n)) ##1.62% of cells were unclassified

emd_raw |>
  group_by(medium_cell_type)|>
  count(medium_cell_type)|>
  ungroup()|>
  mutate(perc = n/sum(n)) ##0.62% of cells were double positive T cells

p_emd_comp <- emd_raw |>
  filter(medium_cell_type != "Unclassified")|>
  filter(medium_cell_type != "Double positive T cells")|>
  group_by(sample_id_alt)|>
  count(medium_cell_type)|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(sample_id_alt), y = perc*100, fill = medium_cell_type))+
  geom_bar(stat = 'identity', alpha = 0.8) +
  scale_fill_manual(values = emd_palette_medium)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="Cell type [%]", fill = "Cell type")+
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9)
  )

ggsave(filename = "p_emd_comp.svg", plot = p_emd_comp, path = output_dir,
       width = 6, height = 5)

emd_raw |>
  filter(medium_cell_type != "Unclassified")|>
  filter(medium_cell_type != "Double positive T cells")|>
  group_by(sample_id_alt)|>
  count(medium_cell_type)|>
  mutate(perc = n/sum(n)) |>
  print(n=Inf)

emd_raw |>
  filter(medium_cell_type != "Unclassified")|>
  filter(medium_cell_type != "Double positive T cells")|>
  group_by(sample_id_alt)|>
  count(medium_cell_type)|>
  mutate(perc = n/sum(n)) |>
  print(n=Inf)

emd_raw |>
  filter(major_cell_type != "Unclassified")|>
  filter(major_cell_type != "Double positive T cells")|>
  group_by(sample_id_alt)|>
  count(major_cell_type)|>
  mutate(perc = n/sum(n)) |>
  print(n=Inf)

emd_raw |>
  filter(minor_cell_type != "Unclassified")|>
  filter(minor_cell_type != "Double positive T cells")|>
  group_by(sample_id_alt)|>
  count(minor_cell_type)|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(sample_id_alt), y = perc*100, fill = minor_cell_type))+
  geom_bar(stat = 'identity', alpha = 0.8) +
  #scale_fill_manual(values = emd_palette_medium)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="Cell type [%]", fill = "Cell type")+
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9)
  )

p_emd_mm_comp <- emd_raw |>
  filter(functional_minor_cell_type %in% c("Ki-67-/PD-L1- Myeloma", "Ki-67+/PD-L1- Myeloma", "Ki-67+/PD-L1+ Myeloma"))|>
  group_by(sample_id_alt)|>
  count(functional_minor_cell_type)|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(sample_id_alt), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity', alpha = 0.8) +
  scale_fill_manual(values = c("Ki-67-/PD-L1- Myeloma" = "#FF8F83",
                               "Ki-67+/PD-L1- Myeloma" = "#FF5E48",
                               "Ki-67+/PD-L1+ Myeloma" = '#B22222'))+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y="Myeloma subpopulation [%]", fill = "Subpopulation")+
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9)
  )

ggsave(filename = "p_emd_mm_comp.svg", plot = p_emd_mm_comp, path = output_dir,
       width = 6, height = 5)

emd_raw |>
  filter(functional_minor_cell_type %in% c("Ki-67-/PD-L1- Myeloma", "Ki-67+/PD-L1- Myeloma", "Ki-67+/PD-L1+ Myeloma"))|>
  group_by(sample_id_alt)|>
  count(functional_minor_cell_type)|>
  mutate(perc = n/sum(n)) 


## Distance of immune cells in PMD/EMD samples
distance_emd <- read.csv("full_results/emd/emd_dists.csv")
str(distance_emd)

colnames(distance_emd)<- gsub("\\.", " ", colnames(distance_emd))
colnames(distance_emd)<- gsub("M1 like M  phi ", "M1-like macrophage", colnames(distance_emd))
colnames(distance_emd)<- gsub("M2 like M  phi ", "M2-like macrophage", colnames(distance_emd))
colnames(distance_emd)<- gsub("CD8 GZMB  Tmem", "CD8 GZMB+ Tmem", colnames(distance_emd))
colnames(distance_emd)<- gsub("MDSCs", "Myeloid/MDSCs", colnames(distance_emd))
distance_emd <- distance_emd |>
  mutate(
    medium_cell_type = str_replace(medium_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
    medium_cell_type = str_replace(medium_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"),
    medium_cell_type = str_replace(medium_cell_type, "CD8 GZMB  Tmem", "CD8 GZMB+ Tmem"),
    medium_cell_type = str_replace(medium_cell_type, "MDSCs", "Myeloid/MDSCs")
  )

str(distance_emd)

distance_emd |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = c('Myeloid/MDSCs':'NK cells'), names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = emd_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "All preCNs")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

p_emd_ridge_pmd1 <- distance_emd |>
  filter(library_id == "15K30184")|>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = c('Myeloid/MDSCs':'NK cells'), names_to = "cell_partner", values_to = "cell_distance") |>
  filter(!cell_partner %in% c("Myeloma", "Double positive T cells", "Fibroblasts"))|>
  filter(cell_distance < 400)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = emd_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "PMD1")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

ggsave(filename = "p_emd_ridge_pmd1.svg", plot = p_emd_ridge_pmd1, path = output_dir,
       width = 4, height = 4)

p_emd_ridge_pmd2 <- distance_emd |>
  filter(library_id == "22K77729")|>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = c('Myeloid/MDSCs':'NK cells'), names_to = "cell_partner", values_to = "cell_distance") |>
  filter(!cell_partner %in% c("Myeloma", "Double positive T cells", "Fibroblasts"))|>
  filter(cell_distance < 400)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = emd_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "PMD2")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

ggsave(filename = "p_emd_ridge_pmd2.svg", plot = p_emd_ridge_pmd2, path = output_dir,
       width = 4, height = 4)

p_emd_ridge_bm <- distance_new |>
  filter(!is.na(cn_celltypes)) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner %in% c("B cells", "CD4 Tmem", "CD8 GZMB+ Tmem", "CD8 Tex", "CD8 Tmem", "ECs", "M1-like macrophage", 
                             "M2-like macrophage", "MDSCs", "MoDC", "Monocytes", "NK cells", "Treg"))|>
  mutate(medium_cell_type = str_replace(medium_cell_type, "MDSCs", "Myeloid/MDSCs"),
         cell_partner = str_replace(cell_partner, "MDSCs", "Myeloid/MDSCs"))|>
  filter(cell_distance < 400)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = color_palette_medium_new)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "Bone marrow")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

ggsave(filename = "p_emd_ridge_bm.svg", plot = p_emd_ridge_bm, path = output_dir,
       width = 4, height = 4)

p_emd_ridge_emd3 <- distance_emd |>
  filter(library_id == "22A64389")|>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = c('Myeloid/MDSCs':'NK cells'), names_to = "cell_partner", values_to = "cell_distance") |>
  filter(!cell_partner %in% c("Myeloma", "Double positive T cells", "Fibroblasts"))|>
  filter(cell_distance < 400)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = emd_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "EMD3")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

ggsave(filename = "p_emd_ridge_emd3.svg", plot = p_emd_ridge_emd3, path = output_dir,
       width = 4, height = 4)

p_emd_ridge_emd1 <- distance_emd |>
  filter(library_id == "23K58040")|>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = c('Myeloid/MDSCs':'NK cells'), names_to = "cell_partner", values_to = "cell_distance") |>
  filter(!cell_partner %in% c("Myeloma", "Double positive T cells", "Fibroblasts"))|>
  filter(cell_distance < 400)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = emd_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "EMD1")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

ggsave(filename = "p_emd_ridge_emd1.svg", plot = p_emd_ridge_emd1, path = output_dir,
       width = 4, height = 4)

p_emd_ridge_emd2 <- distance_emd |>
  filter(library_id == "22T25546")|>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = c('Myeloid/MDSCs':'NK cells'), names_to = "cell_partner", values_to = "cell_distance") |>
  filter(!cell_partner %in% c("Myeloma", "Double positive T cells", "Fibroblasts"))|>
  filter(cell_distance < 400)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = emd_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "EMD2")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

ggsave(filename = "p_emd_ridge_emd2.svg", plot = p_emd_ridge_emd2, path = output_dir,
       width = 4, height = 4)


# sample_id == "15K30184" ~ "PMD1",
# sample_id == "22A64389" ~ "EMD3",
# sample_id == "22K77729" ~ "PMD2",
# sample_id == "23K58040" ~ "EMD1",
# sample_id == "22T25546" ~ "EMD2",
