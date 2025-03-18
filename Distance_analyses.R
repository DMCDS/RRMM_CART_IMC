
###
### DISTANCE ANALYSIS FOR BASELINE CNs ASSOCIATED WITH SURVIVAL ----
###

## Data wrangling

distance_new <- read.csv("full_results/distances/baseline_distances.csv")

distance_new$cn_celltypes <- as.character(distance_new$cn_celltypes)

colnames(distance_new)<- gsub("\\.", " ", colnames(distance_new))
colnames(distance_new)<- gsub("M1 like M  phi ", "M1-like macrophage", colnames(distance_new))
colnames(distance_new)<- gsub("M2 like M  phi ", "M2-like macrophage", colnames(distance_new))
colnames(distance_new)<- gsub("CD8 GZMB  Tmem", "CD8 GZMB+ Tmem", colnames(distance_new))
distance_new <- distance_new |>
  mutate(
    medium_cell_type = str_replace(medium_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
    medium_cell_type = str_replace(medium_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"),
    medium_cell_type = str_replace(medium_cell_type, "CD8 GZMB  Tmem", "CD8 GZMB+ Tmem")
  )

str(distance_new)

distance_200_summary <- distance_new_combined |>
  filter(!is.na(cn_celltypes)) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)|>
  group_by(cn_celltypes, cell_partner) %>%
  summarize(
    mean_distance   = mean(cell_distance, na.rm = TRUE),
    se_distance     = sd(cell_distance, na.rm = TRUE) / sqrt(n()),
    median_distance = median(cell_distance, na.rm = TRUE),
    n               = n(),
    .groups = "drop"    # This ungroups afterwards
  )

distance_50_summary <- distance_new_combined |>
  filter(!is.na(cn_celltypes)) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 50)|>
  group_by(cn_celltypes, cell_partner) %>%
  summarize(
    mean_distance   = mean(cell_distance, na.rm = TRUE),
    se_distance     = sd(cell_distance, na.rm = TRUE) / sqrt(n()),
    median_distance = median(cell_distance, na.rm = TRUE),
    n               = n(),
    .groups = "drop"    # This ungroups afterwards
  )

df<- distance_new_combined |>
  filter(!is.na(cn_celltypes)) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)

pairwise_results_ref <- df %>%
  group_by(cell_partner) %>%
  pairwise_wilcox_test(
    cell_distance ~ cn_celltypes,  # Numeric ~ Factor
    ref.group      = "All",        # reference level in cn_celltypes
    p.adjust.method = "BH"
  )


## Looking at all all cell types in comparison to CNs on medium level
p_blborder10cn12_myeloma_distance_all <- distance_new |>
  filter(!is.na(cn_celltypes)) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = color_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "All preCNs")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

p_blborder10cn12_myeloma_distance_allgrey <- distance_new |>
  filter(!is.na(cn_celltypes)) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner)))+
  geom_density_ridges(scale = 2, alpha = 0.3, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5, fill = "lightgrey",
                      color = "darkgrey") +
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "All preCNs")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))


p_blborder10cn12_myeloma_distance_CN0 <- distance_new |>
  filter(cn_celltypes == 0) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = color_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "preCN0")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

p_blborder10cn12_myeloma_distance_CN11 <- distance_new |>
  filter(cn_celltypes == 11) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = color_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "preCN11")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

p_blborder10cn12_myeloma_distance_CN8 <- distance_new |>
  filter(cn_celltypes == 8) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = color_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "preCN8")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

p_blborder10cn12_myeloma_distance_CN9 <- distance_new |>
  filter(cn_celltypes == 9) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = color_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "preCN9")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

p_blborder10cn12_myeloma_distance_CN1 <- distance_new |>
  filter(cn_celltypes == 1) |>
  filter(medium_cell_type == "Myeloma") |>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner != "Myeloma")|>
  filter(cell_distance < 200)|>
  ggplot(aes(x = cell_distance, y = fct_rev(cell_partner), fill = cell_partner))+
  geom_density_ridges(scale = 2, alpha = 0.5, bandwidth = 10, 
                      quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = color_palette_medium)+
  guides(fill = "none")+
  labs(x = "Distance [um]", y ="", title = "preCN1")+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9))

## Saving composition  plots

# Retrieve the plots using the pattern
pattern_myeloma_distance_blborder10cn12 <- "p_blborder10cn12_myeloma_distance_[a-zA-Z0-9]"
plots_myeloma_distance <- mget(ls(pattern = pattern_myeloma_distance_blborder10cn12))

# Define the directory to save the plots
output_dir <- "blborder10cn12_distance"
dir.create(output_dir, showWarnings = FALSE)

# Function to save the plots
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

# Save the plots to the specified directory
save_plot(plots_myeloma_distance, output_dir, width = 3.5, height = 4)


## Comparing distances
# Create a copy of the distance_new dataset and label it as "All"
distance_new_all <- distance_new
distance_new_all$cn_celltypes <- "All"  # This changes the CN column to "All" for all rows

# Combine the original and the "All" datasets
distance_new_combined <- dplyr::bind_rows(distance_new, distance_new_all)

## Comparing the distance between immunosuppressive cells and CD4 Tmem
distance_new |>
  filter(medium_cell_type == "CD4 Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "MDSCs") |>
  filter(cell_distance < 200) |>
ggplot(plot_data, aes(x = factor(CN, levels = c("9", "11", "8", "0", "All")), y = !!sym(current_col))) +
  geom_violin(width = 0.6, color = fill_color, fill = fill_color, alpha = 0.5)+
  geom_pwc(ref.group = "All", hide.ns = F, method = "wilcox.test", p.adjust.method = "fdr",
           label = "p.adj", size = 0.2)+
  scale_x_discrete(labels = c("preCN9", "preCN11", "preCN8", "preCN0", "All"))+
  scale_y_continuous(breaks = c(0,25,50,75,100))+
  coord_flip()+
  theme_classic() +
  labs(
    x = "",
    y = paste("Distance to nearest", cell_type, "[µm]")
  )+
  theme(
    axis.text.x = element_text(
      size = 10            # Adjust text size as needed
    ))


# Focus on PFS-associated CNs
p_blborder10cn12_mm_mdsc_cnoi <- distance_new_combined |>
  filter(medium_cell_type == "Myeloma") |>
 filter(cn_celltypes %in% c('All', 9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "MDSCs") |>
  filter(cell_distance < 50) |>
  ggplot(aes(x = factor(cn_celltypes, levels = c("9", "11", "8", "0", "All")), y = cell_distance,
             fill = cell_partner, color = cell_partner)) +
  geom_boxplot(outliers = T, alpha = 0.5)+
  geom_pwc(ref.group = "All", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  scale_x_discrete(labels = c("preCN9", "preCN11", "preCN8", "preCN0", "All"))+
  scale_y_continuous(breaks = c(0,10,20,30,40,50))+
  scale_fill_manual(values = color_palette_medium)+
  scale_color_manual(values = color_palette_medium)+
  guides(color = "none", fill = "none")+
  labs(y ="Nearest myeloid-lineage cell [µm]", x = "")+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10))

p_blborder10cn12_mm_cd4tex_cnoi <- distance_new_combined |>
  filter(medium_cell_type == "Myeloma") |>
  filter(cn_celltypes %in% c('All', 9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "CD4 Tex") |>
  filter(cell_distance < 50) |>
  ggplot(aes(x = factor(cn_celltypes, levels = c("9", "11", "8", "0", "All")), y = cell_distance,
             fill = cell_partner, color = cell_partner)) +
  geom_boxplot(outliers = T, alpha = 0.5)+
  geom_pwc(ref.group = "All", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  scale_x_discrete(labels = c("preCN9", "preCN11", "preCN8", "preCN0", "All"))+
  scale_y_continuous(breaks = c(0,10,20,30,40,50))+
  scale_fill_manual(values = color_palette_medium)+
  scale_color_manual(values = color_palette_medium)+
  guides(color = "none", fill = "none")+
  labs(y ="Distance of closest myeloid-lineage cell [µm]", x = "")+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10))
        
partner_celltypes <- unique(distance_new_combined$medium_cell_type) 

# Create an empty list to store all plots
plot_list_mm_partner <- list()

output_dir <- "blborder10cn12_distance"
dir.create(output_dir, showWarnings = FALSE)

for (partner in partner_celltypes) {
  
  message("Processing partner: ", partner)
  
  # Filter and reshape data
  plot_data <- distance_new_combined %>%
    filter(medium_cell_type == "Myeloma") %>%
    filter(cn_celltypes %in% c('All','9','0','8','11')) %>%
    pivot_longer(
      cols = order_medium,       # <-- adjust if your real col names differ
      names_to = "cell_partner",
      values_to = "cell_distance"
    ) %>%
    filter(cell_partner == partner) %>%
    filter(cell_distance < 50)   # or whatever threshold you need
  
  # Create plot
  p <- ggplot(
    plot_data,
    aes(
      x = factor(cn_celltypes, levels = c("9", "11", "8", "0", "All")),
      y = cell_distance,
      fill = cell_partner,
      color = cell_partner
    )
  ) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1, alpha = 0.5) +
    geom_pwc(
      ref.group = "All",
      method = "wilcox.test",
      label = "p.adj",
      p.adjust.method = "fdr",
      step.increase = 0.09
    ) +
    scale_x_discrete(labels = c("preCN9", "preCN11", "preCN8", "preCN0", "All")) +
    scale_y_continuous(breaks = c(0,10,20,30,40,50)) +
    scale_fill_manual(values = color_palette_medium) +   # make sure your palette has keys for all partner names
    scale_color_manual(values = color_palette_medium) +
    guides(color = "none", fill = "none") +
    labs(
      x = "",
      y = paste("Distance to", partner, "cell [µm]")  # change label if you like
    ) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10)
    )
  
  # Add this plot to the list
  plot_list_mm_partner[[partner]] <- p
  
  # If you want to save to file:
  # Clean up the partner name for file usage (remove spaces, special chars, etc.)
  file_suffix <- gsub("[^A-Za-z0-9]+", "_", partner)
  file_name <- paste0("boxplot_mm_partner_", file_suffix, ".svg")
  file_path <- file.path(output_dir, file_name)
  
  # Save the plot
  ggsave(file_path, p, width = 5, height = 2)
}


# Create an empty list to store all plots
plot_list_crs_mm_partner <- list()

output_dir <- "blborder10cn12_distance"
dir.create(output_dir, showWarnings = FALSE)

for (partner in partner_celltypes) {
  
  message("Processing partner: ", partner)
  
  # Filter and reshape data
  plot_data <- distance_new_combined %>%
    filter(medium_cell_type == "Myeloma") %>%
    filter(cn_celltypes %in% c('All','0','1')) %>%
    pivot_longer(
      cols = order_medium,       # <-- adjust if your real col names differ
      names_to = "cell_partner",
      values_to = "cell_distance"
    ) %>%
    filter(cell_partner == partner) %>%
    filter(cell_distance < 50)   # or whatever threshold you need
  
  # Create plot
  p <- ggplot(
    plot_data,
    aes(
      x = factor(cn_celltypes, levels = c("1", "0", "All")),
      y = cell_distance,
      fill = cell_partner,
      color = cell_partner
    )
  ) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1, alpha = 0.5) +
    geom_pwc(
      ref.group = "All",
      method = "wilcox.test",
      label = "p.adj",
      p.adjust.method = "fdr",
      step.increase = 0.09
    ) +
    scale_x_discrete(labels = c("preCN1", "preCN0", "All")) +
    scale_y_continuous(breaks = c(0,10,20,30,40,50)) +
    scale_fill_manual(values = color_palette_medium) +   # make sure your palette has keys for all partner names
    scale_color_manual(values = color_palette_medium) +
    guides(color = "none", fill = "none") +
    labs(
      x = "",
      y = paste("Distance to", partner, "cell [µm]")  # change label if you like
    ) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10)
    )
  
  # Add this plot to the list
  plot_list_crs_mm_partner[[partner]] <- p
  
  # If you want to save to file:
  # Clean up the partner name for file usage (remove spaces, special chars, etc.)
  file_suffix <- gsub("[^A-Za-z0-9]+", "_", partner)
  file_name <- paste0("boxplot_crs_mm_partner_", file_suffix, ".svg")
  file_path <- file.path(output_dir, file_name)
  
  # Save the plot
  ggsave(file_path, p, width = 5, height = 1.5)
}


p_blborder10cn12_cd4tmem_myeloma_cnoi <- distance_new |>
  filter(medium_cell_type == "CD4 Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Myeloma") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#FF0000', alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
          step.increase = 0.09) +
  labs(y ="Distance of closest myeloma cell\n to each CD4 Tmem [µm]", x = "")+
  theme_classic()
  

p_blborder10cn12_cd4tmem_m2mac_cnoi <- distance_new |>
  filter(medium_cell_type == "CD4 Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "M2-like macrophage") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#2ca02c", alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest M2-like macrophage\n to each CD4 Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd4tmem_treg_cnoi <- distance_new |>
  filter(medium_cell_type == "CD4 Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Treg") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#D8BFD8', alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest Treg\n to each CD4 Tmem [µm]", x = "")+
  theme_classic()

# All CNs
p_blborder10cn12_cd4tmem_mdsc_all <- distance_new |>
  filter(medium_cell_type == "CD4 Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "MDSCs") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#FFD700", alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          bracket.nudge.y = -0.7step.increase = 0.09) +
  labs(y ="Distance of closest MDSC\n to each CD4 Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd4tmem_myeloma_all <- distance_new |>
  filter(medium_cell_type == "CD4 Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Myeloma") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#FF0000', alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest myeloma cell\n to each CD4 Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd4tmem_m2mac_all <- distance_new |>
  filter(medium_cell_type == "CD4 Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "M2-like macrophage") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#2ca02c", alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest M2-like macrophage\n to each CD4 Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd4tmem_treg_all <- distance_new |>
  filter(medium_cell_type == "CD4 Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Treg") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#D8BFD8', alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest Treg\n to each CD4 Tmem [µm]", x = "")+
  theme_classic()


## Comparing the distance between immunosuppressive cells and CD8 Tmem
# Focus on PFS-associated CNs
p_blborder10cn12_cd8tmem_mdsc_cnoi <- distance_new |>
  filter(medium_cell_type == "CD8 Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "MDSCs") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#FFD700", alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest MDSC\n to each CD8 Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8tmem_myeloma_cnoi <- distance_new |>
  filter(medium_cell_type == "CD8 Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Myeloma") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#FF0000', alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest myeloma cell\n to each CD8 Tmem [µm]", x = "")+
  theme_classic()


p_blborder10cn12_cd8tmem_m2mac_cnoi <- distance_new |>
  filter(medium_cell_type == "CD8 Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "M2-like macrophage") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#2ca02c", alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest M2-like macrophage\n to each CD8 Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8tmem_treg_cnoi <- distance_new |>
  filter(medium_cell_type == "CD8 Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Treg") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#D8BFD8', alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest Treg\n to each CD8 Tmem [µm]", x = "")+
  theme_classic()

# All CNs
p_blborder10cn12_cd8tmem_mdsc_all <- distance_new |>
  filter(medium_cell_type == "CD8 Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "MDSCs") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = F, fill = "#FFD700", alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest MDSC\n to each CD8 Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8tmem_myeloma_all <- distance_new |>
  filter(medium_cell_type == "CD8 Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Myeloma") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#FF0000', alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest myeloma cell\n to each CD8 Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8tmem_m2mac_all <- distance_new |>
  filter(medium_cell_type == "CD8 Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "M2-like macrophage") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#2ca02c", alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest M2-like macrophage\n to each CD8 Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8tmem_treg_all <- distance_new |>
  filter(medium_cell_type == "CD8 Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Treg") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#D8BFD8', alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest Treg\n to each CD8 Tmem [µm]", x = "")+
  theme_classic()


## Comparing the distance between immunosuppressive cells and CD8 GZMB+ Tmem
# Focus on PFS-associated CNs
p_blborder10cn12_cd8gzmbtmem_mdsc_cnoi <- distance_new |>
  filter(medium_cell_type == "CD8 GZMB+ Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "MDSCs") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#FFD700", alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
            step.increase = 0.09) +
  labs(y ="Distance of closest MDSC\n to each CD8 GZMB+ Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8gzmbtmem_myeloma_cnoi <- distance_new |>
  filter(medium_cell_type == "CD8 GZMB+ Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Myeloma") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#FF0000', alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest myeloma cell\n to each CD8 GZMB+ Tmem [µm]", x = "")+
  theme_classic()


p_blborder10cn12_cd8gzmbtmem_m2mac_cnoi <- distance_new |>
  filter(medium_cell_type == "CD8 GZMB+ Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "M2-like macrophage") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#2ca02c", alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest M2-like macrophage\n to each CD8 GZMB+ Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8gzmbtmem_treg_cnoi <- distance_new |>
  filter(medium_cell_type == "CD8 GZMB+ Tmem") |>
  filter(cn_celltypes %in% c(9, 0, 8, 11))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Treg") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#D8BFD8', alpha = 0.5)+
  geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
          step.increase = 0.09) +
  labs(y ="Distance of closest Treg\n to each CD8 GZMB+ Tmem [µm]", x = "")+
  theme_classic()

# All CNs
p_blborder10cn12_cd8gzmbtmem_mdsc_all <- distance_new |>
  filter(medium_cell_type == "CD8 GZMB+ Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "MDSCs") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#FFD700", alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          bracket.nudge.y = -0.7step.increase = 0.09) +
  labs(y ="Distance of closest MDSC\n to each CD8 GZMB+ Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8gzmbtmem_myeloma_all <- distance_new |>
  filter(medium_cell_type == "CD8 GZMB+ Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Myeloma") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#FF0000', alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest myeloma cell\n to each CD8 GZMB+ Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8gzmbtmem_m2mac_all <- distance_new |>
  filter(medium_cell_type == "CD8 GZMB+ Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "M2-like macrophage") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#2ca02c", alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest M2-like macrophage\n to each CD8 GZMB+ Tmem [µm]", x = "")+
  theme_classic()

p_blborder10cn12_cd8gzmbtmem_treg_all <- distance_new |>
  filter(medium_cell_type == "CD8 GZMB+ Tmem") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Treg") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#D8BFD8', alpha = 0.5)+
  # geom_pwc(ref.group = 1, method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest Treg\n to each CD8 GZMB+ Tmem [µm]", x = "")+
  theme_classic()


## Saving plots
##
## Saving composition  plots

# Retrieve the plots using the pattern
pattern_distance_cnoi <- "p_blborder10cn12_[a-zA-Z0-9]+_[a-zA-Z0-9]+_cnoi"
pattern_distance_all <- "p_blborder10cn12_[a-zA-Z0-9]+_[a-zA-Z0-9]+_all"
plots_distance_cnoi <- mget(ls(pattern = pattern_distance_cnoi))
plots_distance_all <- mget(ls(pattern = pattern_distance_all))

# Define the directory to save the plots
output_dir <- "blborder10cn12_distance"
dir.create(output_dir, showWarnings = FALSE)

# Function to save the plots
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

# Save the plots to the specified directory
save_plot(plots_distance_cnoi, output_dir, width = 5, height = 2)
save_plot(plots_distance_all, output_dir, width = 6, height = 4)

###
### DISTANCE ANALYSIS FOR BASELINE CNs ASSOCIATED WITH CRS ----
###


## Distances to monocytes
# Across all preCNs
p_blborder10cn12_monocytes_myeloma_all <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Myeloma") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#FF0000', alpha = 0.5)+
  # geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #           step.increase = 0.09) +
  labs(y ="Distance of closest myeloma cell\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_cd4tmem_all <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "CD4 Tmem") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#9400D3', alpha = 0.5)+
  # geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest CD4 Tmem\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_cd8tmem_all <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "CD8 Tmem") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#969696", alpha = 0.5)+
  # geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest CD8 Tmem\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_m2mac_all <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "M2-like macrophage") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#2ca02c", alpha = 0.5)+
  # geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest M2-like macrophage\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_treg_all <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Treg") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#D8BFD8', alpha = 0.5)+
  # geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest Treg\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_mdsc_all <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(!is.na(cn_celltypes))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "MDSCs") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#FFD700", alpha = 0.5)+
  # geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
  #          step.increase = 0.09) +
  labs(y ="Distance of closest MDSC\n to each monocyte [µm]", x = "")+
  theme_classic()


# Across all CNOIs
p_blborder10cn12_monocytes_myeloma_cnoi <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(cn_celltypes %in% c(0,1,3))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Myeloma") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#FF0000', alpha = 0.5)+
  geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest myeloma cell\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_cd4tmem_cnoi <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(cn_celltypes %in% c(0,1,3))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "CD4 Tmem") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill ='#9400D3', alpha = 0.5)+
  geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest CD4 Tmem\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_cd8tmem_cnoi <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(cn_celltypes %in% c(0,1,3))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "CD8 Tmem") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#969696", alpha = 0.5)+
  geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest CD8 Tmem\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_m2mac_cnoi <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(cn_celltypes %in% c(0,1,3))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "M2-like macrophage") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill =  "#2ca02c", alpha = 0.5)+
  geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest M2-like macrophage\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_treg_cnoi <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(cn_celltypes %in% c(0,1,3))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "Treg") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = '#D8BFD8', alpha = 0.5)+
  geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest Treg\n to each monocyte [µm]", x = "")+
  theme_classic()

p_blborder10cn12_monocytes_mdsc_cnoi <- distance_new |>
  filter(medium_cell_type == "Monocytes") |>
  filter(cn_celltypes %in% c(0,1,3))|>
  pivot_longer(cols = order_medium, names_to = "cell_partner", values_to = "cell_distance") |>
  filter(cell_partner == "MDSCs") |>
  filter(cell_distance < 200) |>
  ggplot(aes(x = reorder(cn_celltypes, cell_distance), y = cell_distance)) +
  geom_boxplot(outliers = T, fill = "#FFD700", alpha = 0.5)+
  geom_pwc(ref.group = "3", method = "wilcox.test", label = "p.adj", p.adjust.method = "fdr",
           step.increase = 0.09) +
  labs(y ="Distance of closest MDSC\n to each monocyte [µm]", x = "")+
  theme_classic()



## Saving plots
##
## Saving composition  plots

# Retrieve the plots using the pattern
pattern_distance_monocyte_cnoi <- "p_blborder10cn12_monocytes_[a-zA-Z0-9]+_cnoi"
pattern_distance_monocyte_all <- "p_blborder10cn12_monocytes_[a-zA-Z0-9]+_all"
plots_distance_monocyte_cnoi <- mget(ls(pattern = pattern_distance_monocyte_cnoi))
plots_distance_monocyte_all <- mget(ls(pattern = pattern_distance_monocyte_all))

# Define the directory to save the plots
output_dir <- "blborder10cn12_distance"
dir.create(output_dir, showWarnings = FALSE)

# Function to save the plots
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

# Save the plots to the specified directory
save_plot(plots_distance_monocyte_cnoi, output_dir, width = 2, height = 4)
save_plot(plots_distance_monocyte_all, output_dir, width = 6, height = 4)
