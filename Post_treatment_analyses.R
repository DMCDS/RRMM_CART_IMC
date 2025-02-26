###
# Save plots function with output directory
save_plot <- function(plot_list, output_dir, prefix = "plot", height = 7, width = 7) {
  for (plot_name in names(plot_list)) {
    filename <- paste0(prefix, plot_name, ".svg")
    filepath <- file.path(output_dir, filename)
    svg(filename = filepath, height = height, width = width)
    print(plot_list[[plot_name]])
    dev.off()
  }
}


### Choosing BORDER 10 micron, 12 clusters for further analysis ----
###

survival_bl_border10_cn12_post <- processed_data_border_post$border10_cn12

survival_bl_border10_cn12_post <- survival_bl_border10_cn12_post %>%
  rename(
    CN0 = X0,
    CN1 = X1,
    CN2 = X2,
    CN3 = X3,
    CN4 = X4,
    CN5 = X5,
    CN6 = X6,
    CN7 = X7,
    CN8 = X8,
    CN9 = X9,
    CN10 = X10,
    CN11 = X11
  )

bl_border10_cn12_post_variables <- c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", 
                                "CN6", "CN7", "CN8", "CN9", "CN10", "CN11", "nan")

bl_border10_cn12_post_variables_wonan <- c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", 
                                      "CN6", "CN7", "CN8", "CN9", "CN10", "CN11")

str(survival_bl_border10_cn12_post)


# List to save plots for grid
p_pfs_blborder10cn12_post <- list()
p_os_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_survival"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables) {
  
  # Directly use the formula in survfit2
  fit <- survfit(Surv(PFS_m, PFS_event) ~ ifelse(get(variable, survival_bl_border10_cn12_post) > median(get(variable, survival_bl_border10_cn12_post)), 1, 0), data = survival_bl_border10_cn12_post)
  
  # Generate the survival plot
  p <- ggsurvplot(fit,
                  ylab="PFS", xlab="Months after CAR-T infusion",
                  break.time.by = 3,
                  xlim = c(0,26),
                  censor.size = 3,
                  pval= TRUE, 
                  pval.coord = c(0.5, 0.1), 
                  pval.size = 3,
                  size = 1.15,
                  axes.offset = F,
                  risk.table = F,
                  risk.table.title = "No. at risk",
                  risk.table.height = 0.19,
                  tables.y.text = FALSE,
                  tables.theme = theme_cleantable(base_size = 2)+ theme(plot.title = element_text(size = 10)),
                  conf.int = F,
                  #surv.median.line = "v",
                  ggtheme = theme_classic2(10)+ theme(axis.line = element_line(linewidth = 0.7)),
                  font.title = c(9),
                  font.tickslab = c(9),
                  font.legend = c(9),
                  font.x = c(9, "bold"),
                  font.y = c(9, "bold"),
                  fontsize =3,
                  legend.labs = c("Low", "High"),
                  legend.title = c(paste0("post",variable)),
                  palette = c("#76AAFF", "#000080")
  )
  
  # Extract only the plot component from ggsurvplot
  survival_plot <- p$plot
  
  # Add the plot to the list
  p_pfs_blborder10cn12_post[[variable]] <- survival_plot
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_PFS_CN_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving PFS plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = survival_plot, width = 3, height = 2)
}

# Convert each plot in p_pfs_blborder10cn12 to a grob (grid graphical object)
grobs_pfs_blborder10cn12_post <- lapply(p_pfs_blborder10cn12_post, ggplotGrob)

# Create the file path for saving the grid image
output_grid_file <- file.path(output_dir, "bl_border10cn12_post_pfs_grid.svg")

# Save the grid of survival plots as a PNG file
svg(output_grid_file, width = 16, height = 12)  # Adjust resolution and size
do.call(grid.arrange, c(grobs_pfs_blborder10cn12_post, ncol = 4))  # Arrange plots in a grid, 4 columns
dev.off()

# Prepare a storage object for raw p-values
p_vals <- numeric(length(bl_border10_cn12_post_variables))  # or p_vals <- list() if you prefer

# Create a named vector or list, helpful for reference later:
names(p_vals) <- bl_border10_cn12_post_variables

for (variable in bl_border10_cn12_post_variables) {
  
  # Directly use the formula in survfit2
  fit <- survfit(Surv(OS_m, OS_event) ~ ifelse(get(variable, survival_bl_border10_cn12_post) > median(get(variable, survival_bl_border10_cn12_post)), 1, 0), data = survival_bl_border10_cn12_post)
  
  p_val <- surv_pvalue(
    fit, 
    data = survival_bl_border10_cn12_post
  )$pval
  
  # Store the raw p-value
  p_vals[[variable]] <- p_val
  
  # Generate the survival plot
  p <- ggsurvplot(fit,
                  ylab="OS", xlab="Months after CAR-T infusion",
                  break.time.by = 3,
                  xlim = c(0,26),
                  censor.size = 3,
                  pval= TRUE, 
                  pval.coord = c(0.5, 0.1), 
                  pval.size = 3,
                  size = 1.15,
                  axes.offset = F,
                  risk.table = F,
                  risk.table.title = "No. at risk",
                  risk.table.height = 0.19,
                  tables.y.text = FALSE,
                  tables.theme = theme_cleantable(base_size = 2)+ theme(plot.title = element_text(size = 10)),
                  conf.int = F,
                  #surv.median.line = "v",
                  ggtheme = theme_classic2(10)+ theme(axis.line = element_line(linewidth = 0.7)),
                  font.title = c(9),
                  font.tickslab = c(9),
                  font.legend = c(9),
                  font.x = c(9, "bold"),
                  font.y = c(9, "bold"),
                  fontsize =3,
                  legend.labs = c("Low", "High"),
                  legend.title = c(paste0("post",variable)),
                  palette = c("#88DC88","#009600"))
  
  # Extract only the plot component from ggsurvplot
  survival_plot <- p$plot
  
  # Add the plot to the list
  p_os_blborder10cn12_post[[variable]] <- survival_plot
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_OS_CN_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving OS plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = survival_plot, width = 3, height = 2)
}


# Convert each plot in p_os_blborder10cn12 to a grob (grid graphical object)
grobs_os_blborder10cn12_post <- lapply(p_os_blborder10cn12_post, ggplotGrob)

# Create the file path for saving the grid image
output_grid_file <- file.path(output_dir, "bl_border10cn12_post_os_grid.svg")

# Save the grid of survival plots as a PNG file
svg(output_grid_file, width = 16, height = 12)  # Adjust resolution and size
do.call(grid.arrange, c(grobs_os_blborder10cn12_post, ncol = 4))  # Arrange plots in a grid, 4 columns
dev.off()

p_vals_fdr <- p.adjust(p_vals, method = "fdr")

###
### COMPARISON OF PRE AND POST
###

### Loading of pre percentages

bl_border10cn12_raw_rbind <- bl_border10cn12_raw |>
  ungroup()|>
  mutate(timepoint = "Pre")

bl_border10cn12_post_raw_rbind <- bl_border10cn12_post_raw |>
  ungroup()|>
  mutate(timepoint = "Post")

bl_border10cn12_raw_all <- rbind(bl_border10cn12_raw_rbind, bl_border10cn12_post_raw_rbind)

## Minor cell type changes
p_prepost_all <- bl_border10cn12_raw_all |>
  mutate(minor_cell_type = factor(minor_cell_type, levels = order_minor),
         timepoint = factor(timepoint, levels = c("Pre", "Post"))) |>
  group_by(timepoint, minor_cell_type) |>
  summarize(number = n()) |>
  mutate(perc = number / sum(number)) |>
  ggplot(aes(x = minor_cell_type, y = perc * 100, group = timepoint, 
             fill = minor_cell_type)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1), color = "black") +
  scale_fill_manual(values = color_palette_minor_2) +
  geom_text(aes(label = paste(timepoint, "\n", sprintf("%.1f%%", perc * 100))),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3, alpha = 1) +
  guides(fill = "none", alpha = "none") +
  labs(y = "Frequency of cell type [%]", x = "") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 35)) +
  theme_classic()
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Subset on T cell compartment, immunosuppressive cells and myeloma
p_prepost_8_cd4 <- bl_border10cn12_raw_all |>
  mutate(minor_cell_type = factor(minor_cell_type, levels = order_minor),
         timepoint = factor(timepoint, levels = c("Pre", "Post"))) |>
  group_by(timepoint, minor_cell_type) |>
  summarize(number = n()) |>
  mutate(perc = number / sum(number)) |>
  filter(minor_cell_type %in% c("CD4 Tex", "CD4 Tnaive", "Th17 CD4 Tmem", "Non Th17 CD4 Tmem"))|>
  ggplot(aes(x = minor_cell_type, y = perc * 100, group = timepoint, 
             fill = minor_cell_type)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1), color = "black") +
  scale_fill_manual(values = color_palette_minor_2) +
  geom_text(aes(label = paste(timepoint, "\n", sprintf("%.1f%%", perc * 100))),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3, alpha = 1) +
  guides(fill = "none", alpha = "none") +
  labs(y = "Frequency of cell type [%]", x = "") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10)) +
  scale_x_discrete(labels = c("CD4 Tex", "CD4 Tnaive", "Th17 CD4\nTmem", "CD4 Tmem"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_prepost_8_cd8 <- bl_border10cn12_raw_all |>
  mutate(minor_cell_type = factor(minor_cell_type, levels = order_minor),
         timepoint = factor(timepoint, levels = c("Pre", "Post"))) |>
  group_by(timepoint, minor_cell_type) |>
  summarize(number = n()) |>
  mutate(perc = number / sum(number)) |>
  filter(minor_cell_type %in% c("CD8 Tex", "CD8 Tnaive", "CD8 Tmem", "CD8 GZMB+ Tmem"))|>
  ggplot(aes(x = minor_cell_type, y = perc * 100, group = timepoint, 
             fill = minor_cell_type)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1), color = "black") +
  scale_fill_manual(values = color_palette_minor_2) +
  geom_text(aes(label = paste(timepoint, "\n", sprintf("%.1f%%", perc * 100))),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3, alpha = 1) +
  guides(fill = "none", alpha = "none") +
  labs(y = "Frequency of cell type [%]", x = "") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5), breaks = c(0, 2.5, 5)) +
  scale_x_discrete(labels = c("CD8 Tex", "CD8 Tnaive", "CD8 Tmem", "CD8 GZMB+\nTmem"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_prepost_6_sup <-  bl_border10cn12_raw_all |>
  mutate(minor_cell_type = factor(minor_cell_type, levels = order_minor),
         timepoint = factor(timepoint, levels = c("Pre", "Post"))) |>
  group_by(timepoint, minor_cell_type) |>
  summarize(number = n()) |>
  mutate(perc = number / sum(number)) |>
  filter(minor_cell_type %in% c("MDSCs", "M2-like macrophage","Treg"))|>
  ggplot(aes(x = minor_cell_type, y = perc * 100, group = timepoint, 
             fill = minor_cell_type)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1), color = "black") +
  scale_fill_manual(values = color_palette_minor_2) +
  geom_text(aes(label = paste(timepoint, "\n", sprintf("%.1f%%", perc * 100))),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3, alpha = 1) +
  guides(fill = "none", alpha = "none") +
  labs(y = "Frequency of cell type [%]", x = "") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +
  scale_x_discrete(limits = c("MDSCs", "M2-like macrophage","Treg"), labels = c("MDSCs", "M2-like\nmacrophage","Treg"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p_prepost_6_mm <- bl_border10cn12_raw_all |>
  group_by(timepoint, functional_minor_cell_type) |>
  summarize(number = n()) |>
  mutate(perc = number / sum(number)) |>
  filter(functional_minor_cell_type %in% c("Ki-67pos Myeloma", "PD-L1pos Myeloma", "Myeloma"))|>
  mutate(functional_minor_cell_type = factor(functional_minor_cell_type, levels = c("Myeloma","Ki-67pos Myeloma", "PD-L1pos Myeloma")),
         timepoint = factor(timepoint, levels = c("Pre", "Post"))) |>
  ggplot(aes(x = functional_minor_cell_type, y = perc * 100, group = timepoint, 
             fill = functional_minor_cell_type)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1), color = "black") +
  scale_fill_manual(values = color_palette_2) +
  geom_text(aes(label = paste(timepoint, "\n", sprintf("%.1f%%", perc * 100))),
            position = position_dodge(width = 0.9), vjust = -0.3, size = 3, alpha = 1) +
  guides(fill = "none", alpha = "none") +
  labs(y = "Frequency of cell type [%]", x = "") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  scale_x_discrete(labels = c(c("Ki-67neg\nPD-L1neg\nMyeloma","Ki-67pos\nMyeloma", "PD-L1pos\nMyeloma")))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

bl_border10cn12_raw_all |>
  group_by(timepoint, functional_minor_cell_type, DFCI_id) |>
  summarize(number = n()) |>
  mutate(perc = number / sum(number)) |>
  filter(functional_minor_cell_type %in% c("Ki-67pos Myeloma", "PD-L1pos Myeloma", "Myeloma"))|>
  mutate(functional_minor_cell_type = factor(functional_minor_cell_type, levels = c("Myeloma","Ki-67pos Myeloma", "PD-L1pos Myeloma")),
         timepoint = factor(timepoint, levels = c("Pre", "Post"))) 
  
## Saving composition  plots
  
# Retrieve the plots using the pattern
pattern_composition_8_prepost <- "p_prepost_8_[a-zA-Z]"
plots_composition_8_prepost <- mget(ls(pattern = pattern_composition_8_prepost))

pattern_composition_6_prepost <- "p_prepost_6_[a-zA-Z]"
plots_composition_6_prepost <- mget(ls(pattern = pattern_composition_6_prepost))


# Define the directory to save the plots
output_dir <- "border10cn12_prepost"
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
save_plot(plots_composition_8_prepost, output_dir, width = 3.5, height = 3)
save_plot(plots_composition_6_prepost, output_dir, width = 3, height = 3)


###
### COMPOSITION OF NEIGHBORHOODS IN BORDER10cn12_post ----
###

## Loading raw counts matrix
bl_border10cn12_post_raw <- read.csv("full_results/border/post/border-10_CNs-12_raw-results___01MOpost.csv")

## Displaying cell types per neighborhood
bl_border10cn12_post_raw <- bl_border10cn12_post_raw %>%
  mutate(medium_cell_type = str_replace(medium_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         medium_cell_type = str_replace(medium_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"))

bl_border10cn12_post_raw <- bl_border10cn12_post_raw %>%
  mutate(minor_cell_type = str_replace(minor_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         minor_cell_type = str_replace(minor_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"))

bl_border10cn12_post_raw <- bl_border10cn12_post_raw %>%
  mutate(functional_minor_cell_type = str_replace(functional_minor_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         functional_minor_cell_type = str_replace(functional_minor_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"))

medium_cell_type_labels <- bl_border10cn12_post_raw |> select(medium_cell_type) |> unique() |> as.vector()

functional_minor_cell_type_labels <- bl_border10cn12_post_raw |> select(functional_minor_cell_type) |> unique() |> as.vector()

unique(bl_border10cn12_post_raw$DFCI_id)

bl_border10cn12_post_raw <- left_join(bl_border10cn12_post_raw, dfci_id_mrn, by=join_by("DFCI_id" == "DFCI_number"))

unique(bl_border10cn12_post_raw$MRN)


p_blborder10cn12_post_CNcomposition <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes)|>
  count(medium_cell_type)|>
  # filter(!medium_cell_type %in% c("Myeloma"))|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(cn_celltypes) , y = perc*100, fill = medium_cell_type))+
  geom_bar(stat = 'identity', alpha = 0.8) +
  # scale_fill_manual(values = c(
  #   "black", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
  #   "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
  #   "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
  #   "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
  #   "darkred"# Additional colors from Set1
  # ))+
  scale_fill_manual(values = color_palette_medium)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="postCN", y="Area [%]", fill = "Cell type")+
  ggtitle("postCN cell composition") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9)
  )

# 
# # Filtering significant neighborhoods associated with PFS
# p_blborder10cn12_post_CNpfs <- bl_border10cn12_post_raw |>
#   filter(cn_celltypes != "NaN")|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
#   group_by(cn_celltypes)|>
#   count(medium_cell_type)|>
#   # filter(!medium_cell_type %in% c("Myeloma"))|>
#   mutate(perc = n/sum(n)) |>
#   ggplot(aes(x=as.factor(cn_celltypes) , y = perc*100, fill = medium_cell_type))+
#   geom_bar(stat = 'identity', alpha = 0.8) +
#   # scale_fill_manual(values = c(
#   #   "black", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
#   #   "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
#   #   "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
#   #   "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
#   #   "darkred"# Additional colors from Set1
#   # ))+
#   scale_fill_manual(values = color_palette_medium)+
#   coord_flip()+
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x="Cell neighborhood", y="Area [%]", fill = "Cell type")+
#   ggtitle("BORDER10CN12 - CN associated with survival") +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 9)
#   )
# 
# 
# # Filtering significant neighborhoods associated with CRS
# p_blborder10cn12_post_CNcrs <- bl_border10cn12_post_raw |>
#   filter(cn_celltypes %in% c(0, 1))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 1))) |>
#   group_by(cn_celltypes)|>
#   count(medium_cell_type)|>
#   # filter(!medium_cell_type %in% c("Myeloma"))|>
#   mutate(perc = n/sum(n)) |>
#   ggplot(aes(x=as.factor(cn_celltypes) , y = perc*100, fill = medium_cell_type))+
#   geom_bar(stat = 'identity', alpha = 0.8) +
#   # scale_fill_manual(values = c(
#   #   "black", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
#   #   "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
#   #   "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
#   #   "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
#   #   "darkred"# Additional colors from Set1
#   # ))+
#   scale_fill_manual(values = color_palette_medium)+
#   coord_flip()+
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x="Cell neighborhood", y="Area [%]", fill = "Cell type")+
#   ggtitle("BORDER10CN12 - CN associated with CRS") +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 9)
#   )

## Displaying cell types per neighborhood
bl_border10cn12_post_raw |>
  # filter(cn_celltypes %in% c("Immune", "Myeloma", "Myeloma-MDSCs")) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(cn_celltypes) , y = perc, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(
    "black", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
    "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "darkred"# Additional colors from Set1
  ))+
  coord_flip()+
  theme_classic() 


##
## Distribution of CNs across patients
##

bl_border10cn12_post_raw |>
  group_by(DFCI_id) |>
  count(cn_celltypes) |>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(DFCI_id) , y = perc*100, fill = as.factor(cn_celltypes)))+
  geom_bar(stat = 'identity', alpha = 0.8) +
  scale_fill_manual(values = c(
    "black", "#8DD3C7", "#FFFFB3", "darkred", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
    "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "darkred"# Additional colors from Set1
  ))+
  #scale_fill_manual(values = color_palette_medium)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="DFCI_id", y="Cells included in area [%]", fill = "CN")+
  #ggtitle("BORDER10CN12 - Neighborhood cell composition") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9)
  )


## Displaying cell types per neighborhood
bl_border10cn12_post_raw |>
  filter(cn_celltypes %in% c(0, 11, 8, 9))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(cn_celltypes) , y = perc, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  coord_flip()+
  theme_classic()

###
### ANALYSIS OF FUNCTIONAL CELL TYPES BASED ON NEIGHBORHOODS
###


###
### ANALYSIS OF FUNCTIONAL CELL TYPES BASED ON NEIGHBORHOODS
###

## MYELOMA COMPARTMENT
p_blborder10cn12_post_myeloma <- bl_border10cn12_post_raw |>
  # filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(cn_celltypes) , y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  coord_flip()+
  theme_classic()+
  theme(legend.position = "top")

p_blborder10cn12_post_myeloma_ki67 <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "Ki-67pos Myeloma")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)) , y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_myeloma_pdl1 <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "PD-L1pos Myeloma")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)) , y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

##
## CD4 COMPARTMENT

p_blborder10cn12_post_cd4 <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD4 Tnaive","pSTAT1neg CD4 Tmem", "pSTAT1pos CD4 Tmem", "CD4 Tex", "Treg"))|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(cn_celltypes)  , y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd4_Tnaive <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD4 Tnaive","pSTAT1neg CD4 Tmem", "pSTAT1pos CD4 Tmem", "CD4 Tex", "Treg"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "CD4 Tnaive")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd4_Tex <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD4 Tnaive","pSTAT1neg CD4 Tmem", "pSTAT1pos CD4 Tmem", "CD4 Tex", "Treg"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "CD4 Tex")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd4_Treg <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD4 Tnaive","pSTAT1neg CD4 Tmem", "pSTAT1pos CD4 Tmem", "CD4 Tex", "Treg"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "Treg")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd4_Tmem <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD4 Tnaive","pSTAT1neg CD4 Tmem", "pSTAT1pos CD4 Tmem", "CD4 Tex", "Treg"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type %in% c("pSTAT1neg CD4 Tmem", "pSTAT1pos CD4 Tmem"))|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

##
## CD8 COMPARTMENT

p_blborder10cn12_post_cd8 <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD8 Tnaive","pSTAT1pos CD8 Tmem","pSTAT1neg CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem","pSTAT1pos CD8 GZMBpos Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd8_Tnaive <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD8 Tnaive","pSTAT1pos CD8 Tmem","pSTAT1neg CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem","pSTAT1pos CD8 GZMBpos Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "CD8 Tnaive")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd8_Tex <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD8 Tnaive","pSTAT1pos CD8 Tmem","pSTAT1neg CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem","pSTAT1pos CD8 GZMBpos Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "CD8 Tex")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd8_pos_neg_Tmem <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD8 Tnaive","pSTAT1pos CD8 Tmem","pSTAT1neg CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem","pSTAT1pos CD8 GZMBpos Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "pSTAT1pos CD8 Tmem")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd8_neg_neg_Tmem <-bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD8 Tnaive","pSTAT1pos CD8 Tmem","pSTAT1neg CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem","pSTAT1pos CD8 GZMBpos Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "pSTAT1neg CD8 Tmem")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd8_neg_pos_Tmem <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD8 Tnaive","pSTAT1pos CD8 Tmem","pSTAT1neg CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem","pSTAT1pos CD8 GZMBpos Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "pSTAT1neg CD8 GZMBpos Tmem")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_cd8_pos_pos_Tmem <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("CD8 Tnaive","pSTAT1pos CD8 Tmem","pSTAT1neg CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem","pSTAT1pos CD8 GZMBpos Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "pSTAT1pos CD8 GZMBpos Tmem")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

##
## NK CELL COMPARTMENT

p_blborder10cn12_post_nks <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1neg CD56 dim NK", "pSTAT1neg CD56 bright NK", "pSTAT1pos CD56 dim NK", "pSTAT1pos CD56 bright NK" ))|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip()+
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_neg_dim_nk <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1neg CD56 dim NK", "pSTAT1neg CD56 bright NK", "pSTAT1pos CD56 dim NK", "pSTAT1pos CD56 bright NK"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "pSTAT1neg CD56 dim NK")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_neg_bright_nk <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1neg CD56 dim NK", "pSTAT1neg CD56 bright NK", "pSTAT1pos CD56 dim NK", "pSTAT1pos CD56 bright NK"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "pSTAT1neg CD56 bright NK")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_pos_dim_nk <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1neg CD56 dim NK", "pSTAT1neg CD56 bright NK", "pSTAT1pos CD56 dim NK", "pSTAT1pos CD56 bright NK"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "pSTAT1pos CD56 dim NK")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_pos_bright_nk <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1neg CD56 dim NK", "pSTAT1neg CD56 bright NK", "pSTAT1pos CD56 dim NK", "pSTAT1pos CD56 bright NK"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "pSTAT1pos CD56 bright NK")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

##
## MACROPHAGE COMPARTMENT

p_blborder10cn12_post_macs <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("M1-like macrophage", "M2-like macrophage"))|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip()+
  theme_classic() 

p_blborder10cn12_post_M1mac <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("M1-like macrophage", "M2-like macrophage"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "M1-like macrophage")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_M2mac <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("M1-like macrophage", "M2-like macrophage"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "M2-like macrophage")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

##
## MONOCYTE COMPARTMENT

p_blborder10cn12_post_monos <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes" ))|>
  mutate(perc = n/sum(n)) |>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip()+
  theme_classic()

p_blborder10cn12_post_class_mono <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "Classical monocytes")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_itm_mono <-bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "Intermediate monocytes")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_nonclass_mono <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "Nonclassical monocytes")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

##
## OTHER CELLS COMPARED TO ALL CELLS

p_blborder10cn12_post_adipo <-bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  # filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "Adipocyte")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_MKs <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  # filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "MKs")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_ECs <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  # filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "ECs")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_MDSCs <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  # filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "MDSCs")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

p_blborder10cn12_post_MoDC <- bl_border10cn12_post_raw |>
  filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  # filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"))|>
  mutate(perc = n/sum(n)) |>
  filter(functional_minor_cell_type == "MoDC")|>
  ggplot(aes(x=as.factor(reorder(cn_celltypes, -perc)), y = perc*100, fill = functional_minor_cell_type))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = color_palette_2)+
  labs(x="Cell neighborhood", y="Area [%]", fill = "Functional cell type")+
  guides(fill = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),   
    axis.text.y = element_text(size = 10)    
  )

##
## Saving composition  plots

# Retrieve the plots using the pattern
pattern_composition_blborder10cn12_post <- "p_blborder10cn12_post_[a-zA-Z]"
plots_composition_post <- mget(ls(pattern = pattern_composition_blborder10cn12_post))

# Define the directory to save the plots
output_dir <- "blborder10cn12_post_composition"
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
save_plot(plots_composition_post, output_dir, width = 2.5, height = 2)

# Save composition plots separately
ggsave(filename = "blborder10cn12_post_CNcomposition.svg", plot = p_blborder10cn12_post_CNcomposition, path = output_dir,
       width = 5, height = 5)


ggsave(filename = "border10cn12_post_CNcomposition.svg", plot = p_mva_pfs, path = output_dir,
       width = 7, height = 5)


# Cell compositions per neighborhood saved as tables including percentages

bl_border10cn12_post_cell_composition_table <- bl_border10cn12_post_raw |>
  # filter(cn_celltypes != "NaN")|>
  # mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
  group_by(cn_celltypes)|>
  count(functional_minor_cell_type)|>
  # filter(functional_minor_cell_type %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"))|>
  mutate(perc = n/sum(n)) |>
  select(-n)|>
  pivot_wider(names_from = functional_minor_cell_type, values_from = perc)

write.csv(bl_border10cn12_post_cell_composition_table, "bl_border10cn12_cell_composition_table.csv")

###
###
###


### NEW COMPOSITION DATA for SURVIVAL ----
# MYELOMA
myeloma_fractions_post_all <- bl_border10cn12_post_raw |>
  group_by(DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

myeloma_fractions_post_cn <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(functional_minor_cell_type) %>%
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n / sum(n)) %>%
  ungroup()

myeloma_fractions_post <- rbind(myeloma_fractions_post_all, myeloma_fractions_post_cn)

p_new_ki67mm_all_post <- myeloma_fractions_post|>
  filter(functional_minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes != "NaN")|>
  #mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloma subpopulation [%]", title = "Ki-67pos myeloma")+
  guides(fill = "none", color = "none")+
  #geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p",bracket.nudge.y = 0.05,
  #         step.increase = 0.08)+
  theme_classic()

p_new_ki67mm_cnoi_post <- myeloma_fractions_post|>
  filter(functional_minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloma subpopulation [%]", title = "Ki-67pos myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  theme_classic()

myeloma_fractions_post|>
  filter(functional_minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))


p_new_mmneg_all_post <- myeloma_fractions_post|>
  filter(functional_minor_cell_type == "Myeloma")|>
  filter(cn_celltypes != "NaN")|>
  #mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloma subpopulation [%]", title = "Ki-67neg/PDL1neg myeloma")+
  guides(fill = "none", color = "none")+
  #geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p",bracket.nudge.y = 0.05,
  #         step.increase = 0.08)+
  theme_classic()

p_new_mmneg_cnoi_post <- myeloma_fractions_post|>
  filter(functional_minor_cell_type == "Myeloma")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloma subpopulation [%]", title = "Ki-67neg/PDL1neg myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  theme_classic()

myeloma_fractions_post|>
  filter(functional_minor_cell_type == "Myeloma")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

p_new_pdl1mm_all_post <- myeloma_fractions_post|>
  filter(functional_minor_cell_type == "PD-L1pos Myeloma")|>
  filter(cn_celltypes != "NaN")|>
  #mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloma subpopulation [%]", title = "PD-L1pos myeloma")+
  guides(fill = "none", color = "none")+
  #geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p",bracket.nudge.y = 0.05,
  #         step.increase = 0.08)+
  theme_classic()

p_new_pdl1mm_cnoi_post <- myeloma_fractions_post|>
  filter(functional_minor_cell_type == "PD-L1pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloma subpopulation [%]", title = "PDL1pos myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  theme_classic()

myeloma_fractions_post|>
  filter(functional_minor_cell_type == "PD-L1pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))


## NEW COMPOSITION DATA
# CD4
cd4_fractions_all_post <- bl_border10cn12_post_raw|>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD4 Tnaive","Non Th17 CD4 Tmem", "Th17 CD4 Tmem", "CD4 Tex", "Treg"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

cd4_fractions_cn <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD4 Tnaive","Non Th17 CD4 Tmem", "Th17 CD4 Tmem", "CD4 Tex", "Treg"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

cd4_fractions <- rbind(cd4_fractions_all, cd4_fractions_cn)

cd4_fractions|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloma subpopulation [%]", title = "Ki-67pos myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.05)+
  theme_classic()+
  facet_wrap(vars(minor_cell_type))
# 
# #CD4 Tnaive
# p_new_cd4naive_all <- cd4_fractions|>
#   filter(minor_cell_type == "CD4 Tnaive")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "CD4 Tnaive")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd4naive_cnoi <- cd4_fractions|>
#   filter(minor_cell_type == "CD4 Tnaive")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD4 T cell subpopulation [%]", title = "CD4 Tnaive")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# cd4_fractions|>
#   filter(minor_cell_type == "CD4 Tnaive")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   group_by(cn_celltypes)|>
#   summarise(median = median(perc))
# 
# #CD Tex
# p_new_cd4ex_all <- cd4_fractions|>
#   filter(minor_cell_type == "CD4 Tex")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "CD4 Tex")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd4ex_cnoi <- cd4_fractions|>
#   filter(minor_cell_type == "CD4 Tex")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD4 T cell subpopulation [%]", title = "CD4 Tex")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# cd4_fractions|>
#   filter(minor_cell_type == "CD4 Tex")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   group_by(cn_celltypes)|>
#   summarise(median = median(perc))
# 
# # Non TH17 CD4 Tmem
# p_new_cd4tmem <-cd4_fractions|>
#   filter(minor_cell_type == "Non Th17 CD4 Tmem")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Non TH17 CD4 Tmem")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.05,
#            step.increase = 0.08)+
#   theme_classic()
# 
# p_new_cd4tmem_all <- cd4_fractions|>
#   filter(minor_cell_type == "Non Th17 CD4 Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Non Th17 CD4 Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd4tmem_cnoi <- cd4_fractions|>
#   filter(minor_cell_type == "Non Th17 CD4 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD4 T cell subpopulation [%]", title = "Non Th17 CD4 Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# # Th17 CD4 Tmem
# p_new_cd4th17 <-cd4_fractions|>
#   filter(minor_cell_type == "Th17 CD4 Tmem")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "TH17 CD4 Tmem")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.5,
#            step.increase = 0.05)+
#   theme_classic()
# 
# p_new_cd4th17_all <- cd4_fractions|>
#   filter(minor_cell_type == "Th17 CD4 Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Th17 CD4 Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd4th17_cnoi <- cd4_fractions|>
#   filter(minor_cell_type == "Th17 CD4 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD4 T cell subpopulation [%]", title = "Th17 CD4 Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# cd4_fractions|>
#   filter(minor_cell_type == "Th17 CD4 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   group_by(cn_celltypes)|>
#   summarise(median = median(perc))
# 
# # Treg
# p_new_cd4treg <-cd4_fractions|>
#   filter(minor_cell_type == "Treg")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Treg")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.25,
#            step.increase = 0.05)+
#   theme_classic()
# 
# p_new_cd4treg_all <- cd4_fractions|>
#   filter(minor_cell_type == "Treg")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Treg")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd4treg_cnoi <- cd4_fractions|>
#   filter(minor_cell_type == "Treg")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,100), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD4 T cell subpopulation [%]", title = "Treg")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()


## NEW COMPOSITION DATA
# Effector CD4 FUNCTIONAL LEVEL including pSTAT1
cd4f_fractions_all <- bl_border10cn12_post_raw |>
  group_by(DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos Non Th17 CD4 Tmem", "pSTAT1pos Th17 CD4 Tmem", "pSTAT1neg Non Th17 CD4 Tmem", "pSTAT1neg Th17 CD4 Tmem"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

cd4f_fractions_cn <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos Non Th17 CD4 Tmem", "pSTAT1pos Th17 CD4 Tmem", "pSTAT1neg Non Th17 CD4 Tmem", "pSTAT1neg Th17 CD4 Tmem"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

cd4f_fractions <- rbind(cd4f_fractions_all, cd4f_fractions_cn)

# # pSTAT1pos Non TH17 CD4 Tmem
# p_new_cd4tmempstat1pos_all <- cd4f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos Non Th17 CD4 Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Non Th17 CD4 Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd4tmempstat1pos_cnoi <- cd4f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos Non Th17 CD4 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Non Th17 CD4 Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# # pSTAT1pos Th17 CD4 Tmem
# p_new_cd4th17pstat1pos_all <- cd4f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos Th17 CD4 Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Th17 CD4 Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd4th17pstat1pos_cnoi <- cd4f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos Non Th17 CD4 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Th17 CD4 Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# # pSTAT1pos Non TH17 CD4 Tmem
# p_new_cd4tmempstat1neg_all <- cd4f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg Non Th17 CD4 Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Non Th17 CD4 Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd4tmempstat1neg_cnoi <- cd4f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg Non Th17 CD4 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Non Th17 CD4 Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# # pSTAT1neg Th17 CD4 Tmem
# p_new_cd4th17pstat1neg_all <- cd4f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg Th17 CD4 Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Th17 CD4 Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd4th17pstat1neg_cnoi <- cd4f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg Th17 CD4 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,40), expand = c(0,0), breaks = c(0,10,20,30,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Th17 CD4 Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()


## NEW COMPOSITION DATA - CD8
# CD8
cd8_fractions_all_post <- bl_border10cn12_post_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD8 Tnaive","CD8 GZMB+ Tmem", "CD8 Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

cd8_fractions_cn_post <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD8 Tnaive","CD8 GZMB+ Tmem", "CD8 Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

cd8_fractions_post <- rbind(cd8_fractions_all_post, cd8_fractions_cn_post)

# cd8_fractions_post|>
#   #filter(minor_cell_type == "Ki-67pos Myeloma")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulations [%]")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
#            step.increase = 0.05)+
#   theme_classic()+
#   facet_wrap(vars(minor_cell_type))
# 
# #CD8 Tnaive
# p_new_cd8naive <-cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tnaive")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tnaive")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.05,
#            step.increase = 0.1)+
#   theme_classic()
# 
# p_new_cd8naive_all <- cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tnaive")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tnaive")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd8naive_cnoi <- cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tnaive")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 Tnaive")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #CD8 Tex
# p_new_cd8ex <-cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tex")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,50), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tex")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.6,
#            step.increase = 0.05)+
#   theme_classic()
# 
# p_new_cd8ex_all <- cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tex")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tex")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd8ex_cnoi <- cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tex")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 Tex")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tex")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   group_by(cn_celltypes)|>
#   summarise(median = median(perc))
# 
# #CD8 Tmem
# p_new_cd8tmem <-cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tmem")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tmem")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = -0.08)+
#   theme_classic()
# 
# p_new_cd8tmem_all <- cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd8tmem_cnoi <- cd8_fractions|>
#   filter(minor_cell_type == "CD8 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #CD8 GZMB+ Tmem
# p_new_cd8gzmb <-cd8_fractions|>
#   filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 GZMB+ Tmem")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.08)+
#   theme_classic()
# 
# p_new_cd8gzmb_all <- cd8_fractions|>
#   filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 GZMB+ Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd8gzmb_cnoi <- cd8_fractions|>
#   filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 GZMB+ Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# cd8_fractions|>
#   filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   group_by(cn_celltypes)|>
#   summarise(median = median(perc))

## NEW COMPOSITION DATA
# Effector CD4 FUNCTIONAL LEVEL including pSTAT1
cd8f_fractions_all_post <- bl_border10cn12_post_raw |>
  group_by(DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos CD8 GZMBpos Tmem", "pSTAT1pos CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem", "pSTAT1neg CD8 Tmem"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

cd8f_fractions_cn_post <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos CD8 GZMBpos Tmem", "pSTAT1pos CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem", "pSTAT1neg CD8 Tmem"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

cd8f_fractions_post <- rbind(cd8f_fractions_all_post, cd8f_fractions_cn_post)

# # pSTAT1pos CD8 GZMBpos Tmem
# p_new_cd8tmemgzmbpstat1pos_all <- cd8f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos CD8 GZMBpos Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos CD8 GZMBpos Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd8tmemgzmbpstat1pos_cnoi <- cd8f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos CD8 GZMBpos Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos CD8 GZMBpos Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# # pSTAT1pos CD8 Tmem
# p_new_cd8tmempstat1pos_all <- cd8f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos CD8 Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos GZMBpos CD8 Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd8tmempstat1pos_cnoi <- cd8f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos CD8 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos GZMBpos CD8 Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# # pSTAT1pos CD8 GZMBpos Tmem
# p_new_cd8tmemgzmbpstat1neg_all <- cd8f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg CD8 GZMBpos Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 GZMBpos Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd8tmemgzmbpstat1neg_cnoi <- cd8f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg CD8 GZMBpos Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 GZMBpos Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# # pSTAT1neg CD8 Tmem
# p_new_cd8tmempstat1neg_all <- cd8f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg CD8 Tmem")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 Tmem")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_cd8tmempstat1neg_cnoi <- cd8f_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg CD8 Tmem")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 Tmem")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()


## NEW COMPOSITION DATA - MYELOID
# MYELOID
myeloid_fractions_all_post <- bl_border10cn12_post_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("Classical monocytes", "Intermediate monocytes", "MDSCs", "Nonclassical monocytes",
                                "M1-like macrophage", "M2-like macrophage", "MoDC", "MKs"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

myeloid_fractions_cn_post <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("Classical monocytes", "Intermediate monocytes", "MDSCs", "Nonclassical monocytes",
                                "M1-like macrophage", "M2-like macrophage", "MoDC", "MKs"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

myeloid_fractions_post <- rbind(myeloid_fractions_all_post, myeloid_fractions_cn_post)

myeloid_fractions_post|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0","11"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0","11")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulations [%]")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.05)+
  theme_classic()+
  facet_wrap(vars(minor_cell_type))
# 
# #Classical moncytes
# p_new_classmono <-myeloid_fractions|>
#   filter(minor_cell_type == "Classical monocytes")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Classical monocytes")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.2,,
#            step.increase = 0.08)+
#   theme_classic()
# 
# p_new_classmono_all <- myeloid_fractions|>
#   filter(minor_cell_type == "Classical monocytes")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Classical monocytes")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_classmono_cnoi <- myeloid_fractions|>
#   filter(minor_cell_type == "Classical monocytes")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Classical monocytes")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #Intermediate monocytes
# p_new_itmmono <-myeloid_fractions|>
#   filter(minor_cell_type == "Intermediate monocytes")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,40), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Intermediate monocytes")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.3,
#            step.increase = 0.05)+
#   theme_classic()
# 
# p_new_itmmono_all <- myeloid_fractions|>
#   filter(minor_cell_type == "Intermediate monocytes")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Intermediate monocytes")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_itmmono_cnoi <- myeloid_fractions|>
#   filter(minor_cell_type == "Intermediate monocytes")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Intermediate monocytes")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #Nonclassical monocytes
# p_new_nonclassmono <-myeloid_fractions|>
#   filter(minor_cell_type == "Nonclassical monocytes")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Nonclassical monocytes")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
#            step.increase = 0.08)+
#   theme_classic()
# 
# p_new_nonclassmono_all <- myeloid_fractions|>
#   filter(minor_cell_type == "Nonclassical monocytes")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Nonclassical monocytes")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_nonclassmono_cnoi <- myeloid_fractions|>
#   filter(minor_cell_type == "Nonclassical monocytes")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Nonclassical monocytes")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()

#MDSCs
p_new_mdsc_post <-myeloid_fractions_post|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes %in% c("All", "0","11"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0","11")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Myeloid/MDSCs")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.08)+
  theme_classic()

p_new_mdsc_all_post <- myeloid_fractions_post|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Post-treatment CNs", y = "Myeloid cell subpopulation [%]", title = "Myeloid/MDSCs")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_mdsc_cnoi_post <- myeloid_fractions_post|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes %in% c("All", "0","11"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0","11")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "postCN0", "postCN11"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Myeloid/MDSCs")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj.format",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

myeloid_fractions_post|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))
# 
# #M1-like macrophage
# p_new_m1mac <-myeloid_fractions|>
#   filter(minor_cell_type == "M1-like macrophage")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,5), expand = c(0,0), breaks = c(0,1,2,3,4,5,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M1-like macrophage")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.3)+
#   theme_classic()
# 
# p_new_m1mac_all <- myeloid_fractions|>
#   filter(minor_cell_type == "M1-like macrophage")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M1-like macrophage")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_m1mac_cnoi <- myeloid_fractions|>
#   filter(minor_cell_type == "M1-like macrophage")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,10), expand = c(0,0), breaks = c(0,2,4,6,8,10))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "Myeloid cell subpopulation [%]", title = "M1-like macrophage")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# 
# 
# #M2-like macrophage
# p_new_m2mac <-myeloid_fractions|>
#   filter(minor_cell_type == "M2-like macrophage")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M2-like macrophage")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
#   theme_classic()
# 
# p_new_m2mac_all <- myeloid_fractions|>
#   filter(minor_cell_type == "M2-like macrophage")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M2-like macrophage")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_m2mac_cnoi <- myeloid_fractions|>
#   filter(minor_cell_type == "M2-like macrophage")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,100), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "Myeloid cell subpopulation [%]", title = "M2-like macrophage")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# myeloid_fractions|>
#   filter(minor_cell_type == "M2-like macrophage")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   group_by(cn_celltypes)|>
#   summarise(median = median(perc))
# 
# #MoDC
# p_new_modc <-myeloid_fractions|>
#   filter(minor_cell_type == "MoDC")|>
#   filter(cn_celltypes %in% c("All", "8","9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,20), expand = c(0,0), breaks = c(0,10,20,40,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "MoDC")+
#   guides(fill = "none", color = "none")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
#            step.increase = 0.1)+
#   theme_classic()
# 
# p_new_modc_all <- myeloid_fractions|>
#   filter(minor_cell_type == "MoDC")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "MoDC")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_modc_cnoi <- myeloid_fractions|>
#   filter(minor_cell_type == "MoDC")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,30), expand = c(0,0), breaks = c(0,5,10,15,20,25))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "Myeloid cell subpopulation [%]", title = "MoDC")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #MKs
# p_new_mk_all <- myeloid_fractions|>
#   filter(minor_cell_type == "MKs")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Megakaryocytes")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_mk_cnoi <- myeloid_fractions|>
#   filter(minor_cell_type == "MKs")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Megakaryocytes")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()

## NEW COMPOSITION DATA - MYELOID
# MYELOID
stroma_fractions_all_post <- bl_border10cn12_post_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("ECs", "Adipocyte"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

stroma_fractions_cn_post <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("ECs", "Adipocyte"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

stroma_fractions_post <- rbind(stroma_fractions_all_post, stroma_fractions_cn_post)

stroma_fractions_post|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulations [%]")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.05)+
  theme_classic()+
  facet_wrap(vars(minor_cell_type))
# 
# #Endothelial cells
# p_new_ec_all <- stroma_fractions|>
#   filter(minor_cell_type == "ECs")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Stroma cell subpopulation [%]", title = "ECs")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_ec_cnoi <- stroma_fractions|>
#   filter(minor_cell_type == "ECs")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "Stroma cell subpopulation [%]", title = "ECs")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #Endothelial cells
# p_new_adipocyte_all <- stroma_fractions|>
#   filter(minor_cell_type == "Adipocyte")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "Stroma cell subpopulation [%]", title = "Adipocyte")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_adipocyte_cnoi <- stroma_fractions|>
#   filter(minor_cell_type == "Adipocyte")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "Stroma cell subpopulation [%]", title = "Adipocyte")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()

## NEW COMPOSITION DATA - NK CELLS
# NK
nk_fractions_all_post <- bl_border10cn12_post_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD56 dim NK", "CD56 bright NK"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

nk_fractions_cn_post <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD56 dim NK", "CD56 bright NK"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

nk_fractions_post <- rbind(nk_fractions_all_post, nk_fractions_cn_post)

nk_fractions_post|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0","11"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0","11")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "NK cell subpopulations [%]")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.05)+
  theme_classic()+
  facet_wrap(vars(minor_cell_type))

# #Dim NK cells
# p_new_dimnk_all <- nk_fractions|>
#   filter(minor_cell_type == "CD56 dim NK")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "CD56 dim NK")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_dimnk_cnoi <- nk_fractions|>
#   filter(minor_cell_type == "CD56 dim NK")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "NK cell subpopulation [%]", title = "CD56 dim NK")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #Bright NK cells
# p_new_brightnk_all <- nk_fractions|>
#   filter(minor_cell_type == "CD56 bright NK")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "CD56 bright NK")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_brightnk_cnoi <- nk_fractions|>
#   filter(minor_cell_type == "CD56 bright NK")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "NK cell subpopulation [%]", title = "CD56 bright NK")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()


## NEW COMPOSITION DATA - FUNCTIONAL NK CELLS
# FUNCTIONAL NK
nkf_fractions_all_post <- bl_border10cn12_post_raw |>
  group_by(DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos CD56 dim NK", "pSTAT1pos CD56 bright NK", "pSTAT1neg CD56 dim NK", "pSTAT1neg CD56 bright NK"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

nkf_fractions_cn_post <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos CD56 dim NK", "pSTAT1pos CD56 bright NK", "pSTAT1neg CD56 dim NK", "pSTAT1neg CD56 bright NK"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

nkf_fractions_post <- rbind(nkf_fractions_all_post, nkf_fractions_cn_post)

nkf_fractions_post|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0","11"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0","11")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "NK cell subpopulations [%]")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.05)+
  theme_classic()+
  facet_wrap(vars(functional_minor_cell_type))

# #pSTAT1pos Dim NK cells
# p_new_pstat1posdimnk_all <- nkf_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos CD56 dim NK")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 dim NK")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_pstat1posdimnk_cnoi <- nkf_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos CD56 dim NK")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 dim NK")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #pSTAT1pos Bright NK cells
# p_new_pstat1posbrightnk_all <- nkf_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos CD56 bright NK")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 bright NK")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_pstat1posbrightnk_cnoi <- nkf_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1pos CD56 bright NK")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 bright NK")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #pSTAT1neg Dim NK cells
# p_new_pstat1negdimnk_all <- nkf_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg CD56 dim NK")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 dim NK")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_pstat1negdimnk_cnoi <- nkf_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg CD56 dim NK")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 dim NK")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #pSTAT1neg Bright NK cells
# p_new_pstat1negbrightnk_all <- nkf_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg CD56 bright NK")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 bright NK")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_pstat1negbrightnk_cnoi <- nkf_fractions|>
#   filter(functional_minor_cell_type == "pSTAT1neg CD56 bright NK")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_2)+
#   scale_fill_manual(values = color_palette_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 bright NK")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()

## NEW COMPOSITION DATA - B CELLS
# B
b_fractions_all_post <- bl_border10cn12_post_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("B cells", "MBCs"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

b_fractions_cn_post <- bl_border10cn12_post_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("B cells", "MBCs"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

b_fractions_post <- rbind(b_fractions_all_post, b_fractions_cn_post)

b_fractions_post|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0","11"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0","11")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "B cell subpopulations [%]")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.05)+
  theme_classic()+
  facet_wrap(vars(minor_cell_type))

# #Dim NK cells
# p_new_b_all <- b_fractions|>
#   filter(minor_cell_type == "B cells")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "B cell subpopulation [%]", title = "B cells")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_b_cnoi <- b_fractions|>
#   filter(minor_cell_type == "B cells")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "B cell subpopulation [%]", title = "B cells")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# #Bright NK cells
# p_new_mbcs_all <- b_fractions|>
#   filter(minor_cell_type == "MBCs")|>
#   filter(cn_celltypes != "NaN")|>
#   ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
#   labs(x = "Cellular neighborhood", y = "B cell subpopulation [%]", title = "MBCs")+
#   guides(fill = "none", color = "none")+
#   theme_classic()
# 
# p_new_mbcs_cnoi <- b_fractions|>
#   filter(minor_cell_type == "MBCs")|>
#   filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
#   mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
#   ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
#   geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
#   geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
#   scale_color_manual(values = color_palette_minor_2)+
#   scale_fill_manual(values = color_palette_minor_2)+
#   scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
#   scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
#   labs(x = "", y = "B cell subpopulation [%]", title = "MBCs")+
#   geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
#            step.increase = 0.08)+
#   guides(fill = "none", color = "none")+
#   theme_classic()

# Retrieve the plots using the pattern
pattern_comp_bl_new_all_post <- "p_new_[a-zA-Z0-9]+_all_post"
plots_comp_bl_new_all <- mget(ls(pattern = pattern_comp_bl_new_all_post))
pattern_comp_bl_new_cnoi_post <- "p_new_[a-zA-Z0-9]+_cnoi_post"
plots_comp_bl_new_cnoi<- mget(ls(pattern = pattern_comp_bl_new_cnoi_post))

# Define the directory to save the plots
output_dir <- "blborder10cn12_post_comp_new"
dir.create(output_dir, showWarnings = FALSE)


# Save the plots to the specified directory
save_plot(plots_comp_bl_new_all, output_dir, width = 5, height = 3)
save_plot(plots_comp_bl_new_cnoi, output_dir, width = 3.3, height = 4)

##
## Survival based on myeloma
##

## Preparing cell proportions table
bl_border10cn12_cc_func <- bl_border10cn12_post_raw |>
  group_by(DFCI_id) |>
  count(functional_minor_cell_type) |>
  pivot_wider(names_from = functional_minor_cell_type, values_from = n)

# Transform and process the data
bl_border10cn12_cc_func$DFCI_id <- as.character(bl_border10cn12_cc_func$DFCI_id)
bl_border10cn12_cc_func_variables <- bl_border10cn12_cc_func |> ungroup() |> select(-DFCI_id) |> colnames()
bl_border10cn12_cc_func_id <- bl_border10cn12_cc_func |> select(DFCI_id)

# Calculation cell proportions for subsequent ILR transformation
bl_border10cn12_cp_func <- bl_border10cn12_cc_func
bl_border10cn12_cp_func$total_cells <- rowSums(bl_border10cn12_cp_func[ , -1], na.rm = TRUE)
bl_border10cn12_cp_func[ , -1] <- bl_border10cn12_cp_func[ , -1] / bl_border10cn12_cp_func$total_cells * 100
bl_border10cn12_cp_func <- bl_border10cn12_cp_func[ , -ncol(bl_border10cn12_cp_func)]

# Isometric log ratio transform data
bl_border10cn12_cp_func_ilr <- as.data.frame(bl_border10cn12_cp_func)
bl_border10cn12_cp_func_trans <- as.data.frame(ilr(bl_border10cn12_cp_func))
bl_border10cn12_cp_func_ilr[2:ncol(bl_border10cn12_cp_func_ilr)] <- bl_border10cn12_cp_func_trans
bl_border10cn12_cp_func_ilr$DFCI_id <- as.character(bl_border10cn12_cp_func_ilr$DFCI_id)
bl_border10cn12_cp_func_ilr <- left_join(baseline_master, bl_border10cn12_cp_func_ilr, by = join_by(DFCI_number == DFCI_id))

survival_bl_border10cn12_cp_func_ilr <- bl_border10cn12_cp_func_ilr |>
  select(all_of(survival_baseline_variables), all_of(bl_border10cn12_cc_func_variables))
survival_bl_border10cn12_cp_func_ilr <- survival_bl_border10cn12_cp_func_ilr |>
  filter(Tissue_ID != "BS-21-J45629")

## Subsetting the myeloma fraction
bl_border10cn12_cc_func_myeloma <- bl_border10cn12_cc_func |> select(DFCI_id, Myeloma, `PD-L1pos Myeloma`, `Ki-67pos Myeloma`)
bl_border10cn12_cc_func_myeloma_variables <- colnames(bl_border10cn12_cc_func_myeloma[2:4])

# Calculation cell proportions for subsequent ILR transformation
bl_border10cn12_cp_func_myeloma <- bl_border10cn12_cc_func_myeloma
bl_border10cn12_cp_func_myeloma$total_cells <- rowSums(bl_border10cn12_cp_func_myeloma[ , -1], na.rm = TRUE)
bl_border10cn12_cp_func_myeloma[ , -1] <- bl_border10cn12_cp_func_myeloma[ , -1] / bl_border10cn12_cp_func_myeloma$total_cells * 100
bl_border10cn12_cp_func_myeloma <- bl_border10cn12_cp_func_myeloma[ , -ncol(bl_border10cn12_cp_func_myeloma)]

# Isometric log ratio transform data
bl_border10cn12_cp_func_myeloma_ilr <- as.data.frame(bl_border10cn12_cp_func_myeloma)
bl_border10cn12_cp_func_myeloma_trans <- as.data.frame(ilr(bl_border10cn12_cp_func_myeloma))
bl_border10cn12_cp_func_myeloma_ilr[2:ncol(bl_border10cn12_cp_func_myeloma_ilr)] <- bl_border10cn12_cp_func_myeloma_trans
bl_border10cn12_cp_func_myeloma_ilr$DFCI_id <- as.character(bl_border10cn12_cp_func_myeloma_ilr$DFCI_id)
bl_border10cn12_cp_func_myeloma_ilr <- left_join(baseline_master, bl_border10cn12_cp_func_myeloma_ilr, by = join_by(DFCI_number == DFCI_id))

survival_bl_border10cn12_cp_func_myeloma_ilr <- bl_border10cn12_cp_func_myeloma_ilr |>
  select(all_of(survival_baseline_variables), all_of(bl_border10cn12_cc_func_myeloma_variables))
survival_bl_border10cn12_cp_func_myeloma_ilr <- survival_bl_border10cn12_cp_func_myeloma_ilr |>
  filter(Tissue_ID != "BS-21-J45629")


km_bl_border10cn12_pdl1MM <- survfit(Surv(PFS_m, PFS_event) ~ ifelse(survival_bl_border10cn12_cp_func_myeloma_ilr$`PD-L1pos Myeloma` > median(survival_bl_border10cn12_cp_func_myeloma_ilr$`PD-L1pos Myeloma`), 1, 0), data = survival_bl_border10cn12_cp_func_myeloma_ilr)
summary(km_bl_border10cn12_pdl1MM)
survdiff(Surv(PFS_m, PFS_event) ~ ifelse(survival_bl_border10cn12_cp_func_myeloma_ilr$`PD-L1pos Myeloma` > median(survival_bl_border10cn12_cp_func_myeloma_ilr$`PD-L1pos Myeloma`), 1, 0), data = survival_bl_border10cn12_cp_func_myeloma_ilr)
p_km_bl_border10cn12_pdl1MM <- ggsurvplot(km_bl_border10cn12_pdl1MM,
                                          ylab="PFS", xlab="Months after CAR-T infusion",
                                          break.time.by = 3,
                                          xlim = c(0,26),
                                          censor.size = 3,
                                          pval= TRUE, 
                                          pval.coord = c(0.5, 0.1), 
                                          pval.size = 3,
                                          size = 1.15,
                                          axes.offset = F,
                                          risk.table = F,
                                          risk.table.title = "No. at risk",
                                          risk.table.height = 0.19,
                                          tables.y.text = FALSE,
                                          tables.theme = theme_cleantable(base_size = 2)+ theme(plot.title = element_text(size = 10)),
                                          conf.int = F,
                                          #surv.median.line = "v",
                                          ggtheme = theme_classic2(10)+ theme(axis.line = element_line(linewidth = 0.7)),
                                          font.title = c(9),
                                          font.tickslab = c(9),
                                          font.legend = c(9),
                                          font.x = c(9, "bold"),
                                          font.y = c(9, "bold"),
                                          fontsize =3,
                                          legend.labs = c("MM PD-L1 low", "MM PD-L1 high"),
                                          legend.title = c(paste("")),
                                          palette = c("darkgrey", "#B22222"))

ggsave(filename = "p_km_bl_border10cn12_pdl1MM.svg", plot = p_km_bl_border10cn12_pdl1MM$plot, path = output_dir,
       width = 3, height = 2)


km_bl_border10cn12_ki67MM <- survfit(Surv(PFS_m, PFS_event) ~ ifelse(survival_bl_border10cn12_cp_func_myeloma_ilr$`Ki-67pos Myeloma` > median(survival_bl_border10cn12_cp_func_myeloma_ilr$`Ki-67pos Myeloma`), 1, 0), data = survival_bl_border10cn12_cp_func_myeloma_ilr)
summary(km_bl_border10cn12_ki67MM)
survdiff(Surv(PFS_m, PFS_event) ~ ifelse(survival_bl_border10cn12_cp_func_myeloma_ilr$`Ki-67pos Myeloma` > median(survival_bl_border10cn12_cp_func_myeloma_ilr$`Ki-67pos Myeloma`), 1, 0), data = survival_bl_border10cn12_cp_func_myeloma_ilr)

p_km_bl_border10cn12_ki67MM <- ggsurvplot(km_bl_border10cn12_ki67MM,
                                          ylab="PFS", xlab="Months after CAR-T infusion",
                                          break.time.by = 3,
                                          xlim = c(0,26),
                                          censor.size = 3,
                                          pval= TRUE, 
                                          pval.coord = c(0.5, 0.1), 
                                          pval.size = 3,
                                          size = 1.15,
                                          axes.offset = F,
                                          risk.table = F,
                                          risk.table.title = "No. at risk",
                                          risk.table.height = 0.19,
                                          tables.y.text = FALSE,
                                          tables.theme = theme_cleantable(base_size = 2)+ theme(plot.title = element_text(size = 10)),
                                          conf.int = F,
                                          #surv.median.line = "v",
                                          ggtheme = theme_classic2(10)+ theme(axis.line = element_line(linewidth = 0.7)),
                                          font.title = c(9),
                                          font.tickslab = c(9),
                                          font.legend = c(9),
                                          font.x = c(9, "bold"),
                                          font.y = c(9, "bold"),
                                          fontsize =3,
                                          legend.labs = c("MM Ki67 low", "MM Ki67 high"),
                                          legend.title = c(paste("")),
                                          palette = c("darkgrey", '#DC143C'))

ggsave(filename = "p_km_bl_border10cn12_ki67MM.svg", plot = p_km_bl_border10cn12_ki67MM$plot, path = output_dir,
       width = 3, height = 2)

km_bl_border10cn12_MM <- survfit(Surv(PFS_m, PFS_event) ~ ifelse(survival_bl_border10cn12_cp_func_myeloma_ilr$Myeloma > median(survival_bl_border10cn12_cp_func_myeloma_ilr$Myeloma), 1, 0), data = survival_bl_border10cn12_cp_func_myeloma_ilr)
summary(km_bl_border10cn12_MM)
survdiff(Surv(PFS_m, PFS_event) ~ ifelse(survival_bl_border10cn12_cp_func_myeloma_ilr$Myeloma > median(survival_bl_border10cn12_cp_func_myeloma_ilr$Myeloma), 1, 0), data = survival_bl_border10cn12_cp_func_myeloma_ilr)
p_km_bl_border10cn12_MM <- ggsurvplot(km_bl_border10cn12_MM,
                                      ylab="PFS", xlab="Months after CAR-T infusion",
                                      break.time.by = 3,
                                      xlim = c(0,26),
                                      censor.size = 3,
                                      pval= TRUE, 
                                      pval.coord = c(0.5, 0.1), 
                                      pval.size = 3,
                                      size = 1.15,
                                      axes.offset = F,
                                      risk.table = F,
                                      risk.table.title = "No. at risk",
                                      risk.table.height = 0.19,
                                      tables.y.text = FALSE,
                                      tables.theme = theme_cleantable(base_size = 2)+ theme(plot.title = element_text(size = 10)),
                                      conf.int = F,
                                      #surv.median.line = "v",
                                      ggtheme = theme_classic2(10)+ theme(axis.line = element_line(linewidth = 0.7)),
                                      font.title = c(9),
                                      font.tickslab = c(9),
                                      font.legend = c(9),
                                      font.x = c(9, "bold"),
                                      font.y = c(9, "bold"),
                                      fontsize =3,
                                      legend.labs = c("Myeloma low", "Myeloma high"),
                                      legend.title = c(paste("")),
                                      palette = c("darkgrey", '#FF0000'))

ggsave(filename = "p_km_bl_border10cn12_MM.svg", plot = p_km_bl_border10cn12_MM$plot, path = output_dir,
       width = 3, height = 2)

###
###
###

###
### CRS ANALYSIS ----
###

##
## Preparation of data
##
survival_bl_border10_cn12_post <- processed_data_border_post$border10_cn12

survival_bl_border10_cn12_post$CRS
survival_bl_border10_cn12_post$`CRP(mg/L)`
survival_bl_border10_cn12_post$`LDH(U/L)`
survival_bl_border10_cn12_post$LDH_log
survival_bl_border10_cn12_post$`SerumMSpike(g/dl)`
survival_bl_border10_cn12_post$`PlasmaCellCore(%)`

survival_bl_border10_cn12_post<-survival_bl_border10_cn12_post |>
  mutate(il6_max_log = log(IL6_max),
         ferritin_max_log = log(Ferritin_max))

survival_bl_border10_cn12_post <- survival_bl_border10_cn12_post |>
  mutate(CRS_1 = ifelse(CRS > 0, 1, 0),
         CRS_2 = ifelse(CRS > 1, 1, 0),
         CRS_onset = as.numeric(difftime(CRS_start, infusion_date, units = "days")),
         CRS_event = ifelse(is.na(CRS_onset),0,1),
         CRS_onset_ci = ifelse(is.na(CRS_onset), OS_days, CRS_onset))|>
  ungroup()

summary(glm(CRS_1 ~ `X0`, survival_bl_border10_cn12_post, family = binomial, na.action = na.omit))


# Define the directory to save the CRS plots
output_dir <- "border10cn12_post_crs"
dir.create(output_dir, showWarnings = FALSE)


##
## Log regression models for CRS development, and high CRS development
##

crs1_glm_border10cn12_post <- data.frame(marker = character(),
                                    coefficient = numeric(),
                                    std_error = numeric(),
                                    p_value = numeric(),
                                    lower95 = numeric(),
                                    upper95 = numeric(),
                                    stringsAsFactors = FALSE)

for(variable in bl_border10_cn12_post_variables) {
  # Fit logistic regression model
  formula_str <- paste0("CRS_1 ~ `",variable,"`")
  
  #print(formula_str)
  model <- glm(formula_str, data = survival_bl_border10_cn12_post, family = binomial)
  
  # Extract coefficients, standard errors, and p-values, and confidence intervals
  coef_summary <- summary(model)$coefficients[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
  coinf <- exp(confint(model))
  
  # Create a data frame with results for the current marker
  marker_results <- data.frame(
    marker = variable,
    coefficient = coef_summary["Estimate"],
    std_error = coef_summary["Std. Error"],
    p_value = coef_summary["Pr(>|z|)"],
    lower95 = coinf[2,1],
    upper95 = coinf[2,2])
  
  # Append results to the main data frame
  crs1_glm_border10cn12_post <- rbind(crs1_glm_border10cn12_post, marker_results)
}

crs1_glm_border10cn12_post <- crs1_glm_border10cn12_post|>
  mutate(OR = exp(coefficient), FDR = p.adjust(p_value, method = "fdr"))

crs1_glm_border10cn12_post |>
  ggplot()+
  geom_point(aes(x = reorder(marker, coefficient) , y = OR),
             size = 4, shape = 19, color = "darkred", alpha = 0.7)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = upper95), color = "darkred")+
  geom_text(aes(x = marker, y = -0.7, 
                label = paste("p =", round(p_value,3))))+
  coord_flip(ylim = c(-1, 6))+
  labs(y = "Odds Ratio (95%CI)", x = "CN", title = "Odds of developing CRS >= 1")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the size and style as needed
  )

crs2_glm_border10cn12_post <- data.frame(marker = character(),
                                    coefficient = numeric(),
                                    std_error = numeric(),
                                    p_value = numeric(),
                                    lower95 = numeric(),
                                    upper95 = numeric(),
                                    stringsAsFactors = FALSE)

for(variable in bl_border10_cn12_post_variables) {
  # Fit logistic regression model
  formula_str <- paste0("CRS_2 ~ `",variable,"`")
  
  #print(formula_str)
  model <- glm(formula_str, data = survival_bl_border10_cn12_post, family = binomial)
  
  # Extract coefficients, standard errors, and p-values, and confidence intervals
  coef_summary <- summary(model)$coefficients[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
  coinf <- exp(confint(model))
  
  # Create a data frame with results for the current marker
  marker_results <- data.frame(
    marker = variable,
    coefficient = coef_summary["Estimate"],
    std_error = coef_summary["Std. Error"],
    p_value = coef_summary["Pr(>|z|)"],
    lower95 = coinf[2,1],
    upper95 = coinf[2,2])
  
  # Append results to the main data frame
  crs2_glm_border10cn12_post <- rbind(crs2_glm_border10cn12_post, marker_results)
}

crs2_glm_border10cn12_post <- crs2_glm_border10cn12_post|>
  mutate(OR = exp(coefficient), FDR = p.adjust(p_value, method = "fdr"))

crs2_glm_border10cn12_post |>
  ggplot()+
  geom_point(aes(x = reorder(marker, coefficient) , y = OR),
             size = 4, shape = 19, color = "darkred", alpha = 0.7)+
  geom_linerange(aes(x = marker, ymin = lower95, 
                     ymax = upper95), color = "darkred")+
  geom_text(aes(x = marker, y = -0.7, 
                label = paste("p =", round(p_value,3))))+
  coord_flip(ylim = c(-1, 6))+
  labs(y = "Odds Ratio (95%CI)", x = "CN", title = "Odds of developing CRS >= 2")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_classic()+ 
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the size and style as needed
  )


##
## CI curves
##

# List to save plots for grid
p_crs1_border10cn12_post <- list()
p_crs2_border10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_crs"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables) {
  
  # Directly use the formula in survfit2
  fit <- survfit(Surv(CRS_onset_ci, as.numeric(CRS_1)) ~ ifelse(get(variable, survival_bl_border10_cn12_post) > median(get(variable, survival_bl_border10_cn12_post)), 1, 0), data = survival_bl_border10_cn12_post, stype = 2)
  
  # Generate the plot
  p <- ggsurvplot(fit,
                  fun = "event",
                  ylab="CI CRS", xlim = c(0,15), ylim = c(0,1), xlab="Days after CAR-T infusion",
                  pval= TRUE, pval.coord = c(0.5, 0.75), pval.size = 3,
                  break.time.by = 5,
                  size = 1.15,
                  axes.offset = FALSE,
                  risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                  tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                  conf.int = FALSE,
                  ggtheme = theme_classic2(10),
                  font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                  fontsize=3,
                  legend.labs = c("Low", "High"),
                  legend.title = c(paste(variable)),
                  palette = c("red2", "red4"))

  
  # Extract only the plot component from ggsurvplot
  crs_plot <- p$plot
  
  # Add the plot to the list
  p_crs1_border10cn12_post[[variable]] <- crs_plot
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_crs1_CN_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving PFS plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = crs_plot, width = 3, height = 2)
}

# Convert each plot in p_pfs_blborder10cn12 to a grob (grid graphical object)
grobs_crs1_blborder10cn12_post <- lapply(p_crs1_border10cn12_post, ggplotGrob)

# Create the file path for saving the grid image
output_grid_file <- file.path(output_dir, "bl_border10cn12_post_crs1_grid.svg")

# Save the grid of survival plots as a PNG file
svg(output_grid_file, width = 16, height = 12)  # Adjust resolution and size
do.call(grid.arrange, c(grobs_crs1_blborder10cn12_post, ncol = 4))  # Arrange plots in a grid, 4 columns
dev.off()



for (variable in bl_border10_cn12_post_variables) {
  
  # Directly use the formula in survfit2
  fit <- survfit(Surv(CRS_onset_ci, as.numeric(CRS_2)) ~ ifelse(get(variable, survival_bl_border10_cn12_post) > median(get(variable, survival_bl_border10_cn12_post)), 1, 0), data = survival_bl_border10_cn12_post, stype = 2)
  
  # Generate the plot
  p <- ggsurvplot(fit,
                  fun = "event",
                  ylab="CI CRS", xlim = c(0,15), ylim = c(0,1), xlab="Days after CAR-T infusion",
                  pval= TRUE, pval.coord = c(0.5, 0.75), pval.size = 3,
                  break.time.by = 5,
                  size = 1.15,
                  axes.offset = FALSE,
                  risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                  tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                  conf.int = FALSE,
                  ggtheme = theme_classic2(10),
                  font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                  fontsize=3,
                  legend.labs = c("Low", "High"),
                  legend.title = c(paste(variable)),
                  palette = c("red2", "red4"))
  
  
  # Extract only the plot component from ggsurvplot
  crs_plot <- p$plot
  
  # Add the plot to the list
  p_crs2_border10cn12_post[[variable]] <- crs_plot
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_crs2_CN_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving OS plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = crs_plot, width = 3, height = 2)
}

# Convert each plot in p_os_blborder10cn12 to a grob (grid graphical object)
grobs_crs2_blborder10cn12_post <- lapply(p_crs2_border10cn12_post, ggplotGrob)

# Create the file path for saving the grid image
output_grid_file <- file.path(output_dir, "bl_border10cn12_post_crs2_grid.svg")

# Save the grid of survival plots as a PNG file
svg(output_grid_file, width = 16, height = 12)  # Adjust resolution and size
do.call(grid.arrange, c(grobs_crs2_blborder10cn12_post, ncol = 4))  # Arrange plots in a grid, 4 columns
dev.off()




##
## Distribution of areas between CRS groups
##

p_crs_border10cn12_post_crs1 <- survival_bl_border10_cn12_post|>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc")|>
  filter(area != "nan")|>
  ggplot(aes(x=as.factor(CRS_1), y=area_perc))+
  geom_boxplot()+
  geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkred')+
  geom_pwc(method = "wilcox_test",
           bracket.nudge.y = 0.05, 
           step.increase = 0.15, 
           tip.length = 0.01, label.size = 3)+
  labs(x = "CRS", y = "Area [norm.]", title = "Development of CRS >= 1")+
  facet_wrap(vars(area), scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )


p_crs_border10cn12_post_crs2 <- survival_bl_border10_cn12_post|>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc")|>
  filter(area != "nan")|>
  ggplot(aes(x=as.factor(CRS_2), y=area_perc))+
  geom_boxplot()+
  geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkred')+
  geom_pwc(method = "wilcox_test",
           bracket.nudge.y = 0.05, 
           step.increase = 0.15, 
           tip.length = 0.01, label.size = 3)+
  labs(x = "CRS", y = "Area [norm.]", title = "Development of CRS >= 2")+
  facet_wrap(vars(area), scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )

p_crs_border10cn12_post_toci <- survival_bl_border10_cn12_post|>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc")|>
  filter(area != "nan")|>
  filter(!is.na(Toc_use))|>
  ggplot(aes(x=as.factor(Toc_use), y=area_perc))+
  geom_boxplot()+
  geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkred')+
  geom_pwc(method = "wilcox_test",
           bracket.nudge.y = 0.05, 
           step.increase = 0.15, 
           tip.length = 0.01, label.size = 3)+
  labs(x = "CRS", y = "Area [norm.]", title = "CRS tocilizumab use")+
  facet_wrap(vars(area), scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )


p_crs_border10cn12_post_steroid <- survival_bl_border10_cn12|>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc")|>
  filter(area != "nan")|>
  filter(!is.na(Steroids_use))|>
  ggplot(aes(x=as.factor(Steroids_use), y=area_perc))+
  geom_boxplot()+
  geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkred')+
  geom_pwc(method = "wilcox_test",
           bracket.nudge.y = 0.05, 
           step.increase = 0.15, 
           tip.length = 0.01, label.size = 3)+
  labs(x = "CRS", y = "Area [norm.]", title = "CRS Steroid use")+
  facet_wrap(vars(area), scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )


p_crs_border10cn12_post_il6_max <- survival_bl_border10_cn12_post |>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc") |>
  filter(area != "nan") |>
  ggplot(aes(x = il6_max_log, y = area_perc)) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkred') +
  geom_smooth(method = "lm", color = "darkorange") +  # Change line color to orange
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +  # Add Spearman correlation and p-value
  labs(x = "IL-6 max", y = "Area [norm.]", title = "Correlation with IL-6 max") +
  facet_wrap(vars(area), scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )


p_crs_border10cn12_post_ferritin_max <- survival_bl_border10_cn12_post|>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc")|>
  filter(area != "nan")|>
  ggplot(aes(x=log(ferritin_max_log), y=area_perc))+
  geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkred') +
  geom_smooth(method = "lm", color = "darkorange") +  # Change line color to orange
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +  # Add Spearman correlation and p-value
  labs(x = "IL-6 max", y = "Area [norm.]", title = "Correlation with Ferritin max") +
  facet_wrap(vars(area), scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )


p_crs_border10cn12_post_il6_max <- survival_bl_border10_cn12_post |>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc") |>
  filter(area != "nan") |>
  ggplot(aes(x = il6_max_log, y = area_perc)) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkred') +
  geom_smooth(method = "lm", color = "darkorange") +  # Change line color to orange
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +  # Add Spearman correlation and p-value
  labs(x = "IL-6 max", y = "Area [norm.]", title = "Correlation with IL-6 max") +
  facet_wrap(vars(area), scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )



##
## Saving CRS plots
##

# Example usage
pattern_crs_border10cn12_post <- "p_crs_border10cn12_post_[a-zA-Z]"
pattern2_crs_border10cn12_post <- "p_crs_border10cn12_post_[a-zA-Z]_[a-zA-Z]"
plots_crs_border10cn12_post <- mget(ls(pattern = pattern_crs_border10cn12_post))
plots2_crs_border10cn12_post <- mget(ls(pattern = pattern2_crs_border10cn12_post))

# Define the directory to save the CRS plots
output_dir <- "border10cn12_post_crs"
dir.create(output_dir, showWarnings = FALSE)

# Save the plots to the specified directory
save_plot(plots_crs_border10cn12_post, output_dir, width =8, height =8)
save_plot(plots2_crs_border10cn12_post, output_dir, width =8, height =8)


###

### Interaction analyses ----
###

bl_border10cn12_post_func_int <- read.csv("full_results/border/post/border-10_CNs-12_interactions__functional_01MOpost.csv")

bl_border10cn12_post_func_int <- bl_border10cn12_post_func_int |> select(-X)
colnames(bl_border10cn12_post_func_int)<- gsub("\\.", " ", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int)<- gsub("Ki 67pos Myeloma", "Ki-67pos\nMyeloma", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int)<- gsub("PD L1pos Myeloma", "PD-L1pos\nMyeloma", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int)<- gsub("M1 like M  phi ", "M1-like\nmacrophage", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int)<- gsub("M2 like M  phi ", "M2-like\nmacrophage", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1pos CD8 Tmem", "pSTAT1pos\nCD8 Tmem", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1neg CD8 Tmem", "pSTAT1neg\nCD8 Tmem", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1pos CD8 GZMBpos Tmem", "pSTAT1pos\nCD8 GZMBpos Tmem", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1neg CD8 GZMBpos Tmem", "pSTAT1neg\nCD8 GZMBpos Tmem", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1neg CD4 Tmem", "pSTAT1neg\nCD4 Tmem", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1pos CD4 Tmem", "pSTAT1pos\nCD4 Tmem", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1neg CD56 dim NK", "pSTAT1neg\nCD56 dim NK", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1neg CD56 bright NK", "pSTAT1neg\nCD56 bright NK", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1pos CD56 dim NK", "pSTAT1pos\nCD56 dim NK", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("pSTAT1pos CD56 bright NK", "pSTAT1pos\nCD56 bright NK", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("Classical monocytes", "Classical\nmonocytes", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub("Nonclassical monocytes", "Nonclassical\nmonocytes", colnames(bl_border10cn12_post_func_int))
colnames(bl_border10cn12_post_func_int) <- gsub( "Intermediate monocytes",  "Intermediate\nmonocytes", colnames(bl_border10cn12_post_func_int))

border10cn12_post_cn_values <- unique(bl_border10cn12_post_func_int$cn_celltypes)


for (cn in border10cn12_post_cn_values) {
  # Create a new data frame for each CN value
  new_df <- bl_border10cn12_post_func_int |>
    filter(cn_celltypes == cn) |>
    select(Adipocyte_Adipocyte:'pSTAT1pos\nCD56 dim NK_pSTAT1pos\nCD56 dim NK') |>
    ungroup()|>
    pivot_longer(cols = everything(), names_to = "interaction", values_to = "value")|>
    group_by(interaction)|>
    summarize(mean = mean(value, na.rm = TRUE)) |>
    separate(interaction, into = c("partner1", "partner2"), sep = "_")|>
    pivot_wider(names_from = partner2, values_from = mean) |>
    tibble::column_to_rownames(var = "partner1")
  
  # Reorder rows and columns to match the order in the item_names
  matching_names <- intersect(item_names, rownames(new_df))
  new_df <- new_df[matching_names, matching_names]
  
  # Dynamically create a new variable name for the data frame
  assign(paste0("bl_border10cn12_post_func_int_CN", cn), new_df)
}

# # Ensure diagonal elements are set to zero
# for (cn in border10cn12_post_cn_values) {
#   df_name <- paste0("bl_border10cn12_post_func_int_CN", cn)
#   df <- get(df_name)
#   for (i in intersect(rownames(df), colnames(df))) {
#     df[i, i] <- 0
#   }
#   assign(df_name, df)
# }


# # Plot the chord diagram with the reordered matrix
# chordDiagram(as.matrix(bl_border10cn12_post_func_int_CN2), scale = TRUE, grid.col = color_palette, link.sort = F)
# circos.clear()
# 
# chordDiagram(as.matrix(bl_border10cn12_post_func_int_CN1), scale = TRUE, grid.col = color_palette, link.sort = F)
# circos.clear()
# 



# chordDiagram(as.matrix(bl_border10cn12_post_func_int_CN1), scale = T, grid.col = color_palette)
# circos.clear() 


# Create a list to hold all the data frames
cn_list_bor <- list()
plots_list_bor <- list()

# Assuming the data frames are named bl_border10cn12_post_func_int_CN0 to bl_border10cn12_post_func_int_CN19
for (i in 0:11) {
  # Dynamically assign the data frames to the list
  cn_list_bor[[i+1]] <- get(paste0("bl_border10cn12_post_func_int_CN", i))
}

# Specify the output directory
output_dir <- "border10cn12_post_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# Iterate over the list of data frames and apply the operations
for (i in 1:length(cn_list_bor)) {
  cn_df <- cn_list_bor[[i]]
  
  # Update the diagonal elements to 0 (if needed)
  # for (cn in intersect(rownames(cn_df), colnames(cn_df))) {
  #   cn_df[cn, cn] <- 0
  # }
  
  # Save each plot as a PNG file
  file_name <- paste0("chord_plot_", i, ".png")
  file_path <- file.path(output_dir, file_name)
  print(paste("Saving plot to:", file_path))  # Debugging line
  
  # Create the PNG
  png(file_path, width = 800, height = 800)
  chordDiagram(as.matrix(cn_df), scale = TRUE, grid.col = color_palette)
  circos.clear()
  dev.off()
  
  # Check if the file was saved correctly
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Load the saved image back as a grob using the correct path
  img <- rasterGrob(png::readPNG(file_path), interpolate = TRUE)
  
  # Add a title above each plot
  title <- textGrob(paste("CN", i-1), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and plot into one grob
  plot_with_title <- arrangeGrob(title, img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the grob in the list
  plots_list_bor[[i]] <- plot_with_title
}


# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_bor, ncol = 5)
dev.off()

chordDiagram(as.matrix(bl_border10cn12_post_func_int_CN0), scale = TRUE, grid.col = color_palette)
circos.clear()

MM_partner <- c("Myeloma", "Ki-67pos\nMyeloma", "PD-L1pos\nMyeloma", "Treg", "CD4 Tnaive", "CD4 Tex", "pSTAT1neg\nCD4 Tmem", "pSTAT1pos\nCD4 Tmem",
                "CD8 Tnaive", "CD8 Tex", "pSTAT1pos\nCD8 Tmem", "pSTAT1neg\nCD8 Tmem", "pSTAT1pos\nCD8 GZMBpos Tmem", "pSTAT1neg\nCD8 GZMBpos Tmem",
                "M1-like\nmacrophage", "M2-like\nmacrophage","MDSCs")

# for (cn in border10cn12_post_cn_values) {
#   # Create a new data frame for each CN value
#   new_df <- bl_border10cn12_post_func_int |>
#     filter(cn_celltypes == cn) |>
#     select(Adipocyte_Adipocyte:'pSTAT1pos\nCD56 dim NK_pSTAT1pos\nCD56 dim NK') |>
#     ungroup()|>
#     pivot_longer(cols = everything(), names_to = "interaction", values_to = "value")|>
#     group_by(interaction)|>
#     summarize(mean = mean(value, na.rm = TRUE)) |>
#     separate(interaction, into = c("partner1", "partner2"), sep = "_")|>
#     filter(partner1 %in% MM_partner & partner2 %in% MM_partner ) |>
#     pivot_wider(names_from = partner2, values_from = mean) |>
#     tibble::column_to_rownames(var = "partner1")
#   
#   # Dynamically create a new variable name for the data frame
#   assign(paste0("bl_border10cn12_post_func_int_MM_CN", cn), new_df)
# }
# 
# for(cn in intersect(rownames(bl_border10cn12_post_func_int_MM_CN1), colnames(bl_border10cn12_post_func_int_MM_CN1))) {
#   bl_border10cn12_post_func_int_MM_CN1[cn, cn] = 0
# }

for (cn in border10cn12_post_cn_values) {
  # Create a new data frame for each CN value
  new_df <- bl_border10cn12_post_func_int |>
    filter(cn_celltypes == cn) |>
    select(Adipocyte_Adipocyte:'pSTAT1pos\nCD56 dim NK_pSTAT1pos\nCD56 dim NK') |>
    ungroup()|>
    pivot_longer(cols = everything(), names_to = "interaction", values_to = "value")|>
    group_by(interaction)|>
    summarize(mean = mean(value, na.rm = TRUE)) |>
    separate(interaction, into = c("partner1", "partner2"), sep = "_")|>
    filter(partner1 %in% MM_partner & partner2 %in% MM_partner ) |>
    pivot_wider(names_from = partner2, values_from = mean) |>
    tibble::column_to_rownames(var = "partner1")
  
  # Reorder rows and columns to match the order in the item_names
  matching_names <- intersect(item_names, rownames(new_df))
  new_df <- new_df[matching_names, matching_names]
  
  # Dynamically create a new variable name for the data frame
  assign(paste0("bl_border10cn12_post_func_int_MM_CN", cn), new_df)
}

# # Ensure diagonal elements are set to zero
# for (cn in border10cn12_post_cn_values) {
#   df_name <- paste0("bl_border10cn12_post_func_int_MM_CN", cn)
#   df <- get(df_name)
#   for (i in intersect(rownames(df), colnames(df))) {
#     df[i, i] <- 0
#   }
#   assign(df_name, df)
# }

# Create a list to hold all the data frames
cn_list_MM_bor <- list()
plots_list_MM_bor <- list()

# Assuming the data frames are named bl_border10cn12_post_func_int_CN0 to bl_border10cn12_post_func_int_CN19
for (i in 0:11) {
  # Dynamically assign the data frames to the list
  cn_list_MM_bor[[i+1]] <- get(paste0("bl_border10cn12_post_func_int_MM_CN", i))
}

# Specify the output directory
output_dir <- "border10cn12_post_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# Iterate over the list of data frames and apply the operations
for (i in 1:length(cn_list_MM_bor)) {
  cn_df <- cn_list_MM_bor[[i]]
  
  # Update the diagonal elements to 0 (if needed)
  # for (cn in intersect(rownames(cn_df), colnames(cn_df))) {
  #   cn_df[cn, cn] <- 0
  # }
  
  # Save each plot as a PNG file
  file_name <- paste0("chord_plot_MM_", i, ".png")
  file_path <- file.path(output_dir, file_name)
  print(paste("Saving plot to:", file_path))  # Debugging line
  
  # Create the PNG
  png(file_path, width = 800, height = 800)
  chordDiagram(as.matrix(cn_df), scale = TRUE, grid.col = color_palette)
  circos.clear()
  dev.off()
  
  # Check if the file was saved correctly
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Load the saved image back as a grob using the correct path
  img <- rasterGrob(png::readPNG(file_path), interpolate = TRUE)
  
  # Add a title above each plot
  title <- textGrob(paste("CN", i-1), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and plot into one grob
  plot_with_title <- arrangeGrob(title, img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the grob in the list
  plots_list_MM_bor[[i]] <- plot_with_title
}


# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_MM_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_MM_bor, ncol = 5)
dev.off()


chordDiagram(as.matrix(bl_border10cn12_post_func_int_MM_CN0), scale = T, grid.col = color_palette)
circos.clear() 

chordDiagram(as.matrix(bl_border10cn12_post_func_int_MM_CN11), scale = T, grid.col = color_palette)
circos.clear() 

chordDiagram(as.matrix(bl_border10cn12_post_func_int_MM_CN8), scale = T, grid.col = color_palette)
circos.clear() 

chordDiagram(as.matrix(bl_border10cn12_post_func_int_MM_CN9), scale = T, grid.col = color_palette)
circos.clear() 

###

### ANALYSIS OF HEMATOTOXICITY ----
###

### NEW ANALYSES

### Analysis of platelets for each neighborhood over time

PLT_time_labels <- c(
  "PLT_d90_grade_grouped" = "Day 90",
  "PLT_d60_grade_grouped" = "Day 60",
  "PLT_d30_grade_grouped" = "Day 30"
)

PLT_time_value_labels <- c(
  "0-1" = "0-1",
  ">=2" = "≥2"
)

# List to save plots for grid
p_plt_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  # Generate the survival plot
  p <- survival_bl_border10_cn12_post|>
    mutate(PLT_d90_grade_grouped = case_when(
      PLT_d90_grade %in% 0:1 ~ "0-1",
      PLT_d90_grade > 1 ~ ">=2",
      TRUE ~ NA_character_  ),
      PLT_d60_grade_grouped = case_when(
        PLT_d60_grade %in% 0:1 ~ "0-1",
        PLT_d60_grade > 1 ~ ">=2",
        TRUE ~ NA_character_  ),
      PLT_d30_grade_grouped = case_when(
        PLT_d30_grade %in% 0:1 ~ "0-1",
        PLT_d30_grade > 1 ~ ">=2",
        TRUE ~ NA_character_  )) |>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area != "nan")|>
    mutate(PLT_d90_grade_grouped = factor(PLT_d90_grade_grouped, levels = c("0-1", ">=2")),
           PLT_d60_grade_grouped = factor(PLT_d60_grade_grouped, levels = c("0-1", ">=2")),
           PLT_d30_grade_grouped = factor(PLT_d30_grade_grouped, levels = c("0-1", ">=2"))) |>
    pivot_longer(cols = c(PLT_d90_grade_grouped, PLT_d60_grade_grouped, PLT_d30_grade_grouped), names_to = "PLT_time", values_to = "PLT_time_value")|>
    filter(area == variable) |>
    filter(!is.na(PLT_time_value))|>
    ggplot(aes(x=PLT_time_value, y=area_perc, fill=PLT_time_value, color = PLT_time_value)) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    facet_wrap(~PLT_time, labeller = labeller(PLT_time = PLT_time_labels)) + 
    geom_pwc(aes(group = PLT_time_value), 
             ref.group = "0-1",
             label = "p.adj.format", 
             method = "wilcox.test",
             p.adjust.method = "fdr")+
    labs(x = "Thrombocytopenia [CTCAE grade}", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
    scale_color_manual(values = c("lightblue", "#17becf"))+
    scale_fill_manual(values = c("lightblue", "#17becf"))+
    scale_x_discrete(labels = PLT_time_value_labels) +
    guides(color = "none", fill = "none")+
    labs(x = "Thrombocytopenia [CTCAE]", y = "Area [norm.]")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )

  
  # Add the plot to the list
  p_plt_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_PLT_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_plt_blborder10cn12_post[[variable]], width = 4, height = 4)
}



survival_bl_border10_cn12_post|>
  mutate(PLT_d90_grade_grouped = case_when(
    PLT_d90_grade %in% 0:1 ~ "0-1",
    PLT_d90_grade > 1 ~ ">=2",
    TRUE ~ NA_character_  ),
    PLT_d60_grade_grouped = case_when(
      PLT_d60_grade %in% 0:1 ~ "0-1",
      PLT_d60_grade > 1 ~ ">=2",
      TRUE ~ NA_character_  ),
    PLT_d30_grade_grouped = case_when(
      PLT_d30_grade %in% 0:1 ~ "0-1",
      PLT_d30_grade > 1 ~ ">=2",
      TRUE ~ NA_character_  )) |>
  pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
  filter(area != "nan")|>
  mutate(PLT_d90_grade_grouped = factor(PLT_d90_grade_grouped, levels = c("0-1", ">=2")),
         PLT_d60_grade_grouped = factor(PLT_d60_grade_grouped, levels = c("0-1", ">=2")),
         PLT_d30_grade_grouped = factor(PLT_d30_grade_grouped, levels = c("0-1", ">=2"))) |>
  pivot_longer(cols = c(PLT_d90_grade_grouped, PLT_d60_grade_grouped, PLT_d30_grade_grouped), names_to = "PLT_time", values_to = "PLT_time_value")|>
  filter(!is.na(PLT_time_value)) |>
  group_by(PLT_time, PLT_time_value, area) |>
  summarise(median=median(area_perc))|>
  filter(area == "CN11")
  print(n=Inf)



### Analysis of neutrophiles for each neighborhood over time

ANC_time_labels <- c(
  "ANC_d90_grade_grouped" = "Day 90",
  "ANC_d60_grade_grouped" = "Day 60",
  "ANC_d30_grade_grouped" = "Day 30"
)

ANC_time_value_labels <- c(
  "0-1" = "0-1",
  "1-2" = "1-2",  
  ">=3" = "≥3"
)


# List to save plots for grid
p_anc_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  # Generate the survival plot
  p <- survival_bl_border10_cn12_post|>
    mutate(ANC_d30_grade_grouped = case_when(
      ANC_d30_grade == 0 ~ "0",
      ANC_d30_grade %in% 1:2 ~ "1-2",
      ANC_d30_grade > 2 ~ "≥3",
      TRUE ~ NA_character_  ),
      ANC_d60_grade_grouped = case_when(
        ANC_d60_grade == 0 ~ "0",
        ANC_d60_grade %in% 1:2 ~ "1-2",
        ANC_d60_grade > 2 ~ "≥3",
        TRUE ~ NA_character_  ),
      ANC_d90_grade_grouped = case_when(
        ANC_d90_grade == 0 ~ "0",
        ANC_d90_grade %in% 1:2 ~ "1-2",
        ANC_d90_grade > 2 ~ "≥3",
        TRUE ~ NA_character_  )) |>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area != "nan")|>
    mutate(ANC_d90_grade_grouped = factor(ANC_d90_grade_grouped, levels = c("0", "1-2", "≥3")),
           ANC_d60_grade_grouped = factor(ANC_d60_grade_grouped, levels = c("0", "1-2", "≥3")),
           ANC_d30_grade_grouped = factor(ANC_d30_grade_grouped, levels = c("0", "1-2", "≥3"))) |>
    pivot_longer(cols = c(ANC_d90_grade_grouped, ANC_d60_grade_grouped, ANC_d30_grade_grouped), names_to = "ANC_time", values_to = "ANC_time_value")|>
    filter(area == variable) |>
    filter(!is.na(ANC_time_value))|>
    ggplot(aes(x=ANC_time_value, y=area_perc, fill=ANC_time_value, color = ANC_time_value)) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    facet_wrap(~ANC_time, labeller = labeller(ANC_time = ANC_time_labels)) + 
    geom_pwc(aes(group = ANC_time_value), 
                       ref.group = "0",
                       label = "p.adj.format", 
                       method = "wilcox.test",
                       p.adjust.method = "fdr")+
    labs(x = "Neutropenia [CTCAE grade]", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
    scale_color_manual(values = c("#ADD8E6", "#89CFF0", "#2A52BE"))+
    scale_fill_manual(values = c("#ADD8E6", "#89CFF0", "#2A52BE"))+
    scale_x_discrete(labels = ANC_time_value_labels) +
    guides(color = "none", fill = "none")+
    labs(x = "Neutropenia [CTCAE]", y = "Area [norm.]")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )
  
  #Adjust multiple testing
 # p <- ggadjust_pvalue(p, p.adjust.method = "fdr")
  
  # Add the plot to the list
  p_anc_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_ANC_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_anc_blborder10cn12_post[[variable]], width = 4, height = 4)
}


### Platelet transfusion 

# List to save plots for grid
p_plttx_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  p <- survival_bl_border10_cn12_post|>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area == variable)|>
    ggplot(aes(x=as.factor(plt_d0_d100), y=area_perc, fill=as.factor(plt_d0_d100), color = as.factor(plt_d0_d100))) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    geom_pwc(ref.group = "0",
             label = "p.format", 
             method = "wilcox.test",
             remove.bracket = F)+
    labs(x = "Platelet transfusion", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
    scale_color_manual(values = c("lightblue", "#17becf"))+
    scale_fill_manual(values = c("lightblue", "#17becf"))+
    scale_x_discrete(labels = c("No", "Yes")) +
    guides(color = "none", fill = "none")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )

  
  # Add the plot to the list
  p_plttx_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_PLTTx_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_plttx_blborder10cn12_post[[variable]], width = 2, height = 4)
}


survival_bl_border10_cn12_post|>
  pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
  group_by(plt_d0_d100, area)|>
  summarise(median = median(area_perc)) |>
  filter(area == "CN11")
  
  
### Severe thrombocytopenia

# List to save plots for grid
p_sevplt_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  p <- survival_bl_border10_cn12_post|>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area == variable)|>
    ggplot(aes(x=as.factor(Severe.Thrombocytopenia...50.G.L..Day.0.100), y=area_perc, fill=as.factor(Severe.Thrombocytopenia...50.G.L..Day.0.100), color = as.factor(Severe.Thrombocytopenia...50.G.L..Day.0.100))) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    geom_pwc(ref.group = "0",
             label = "p.format", 
             method = "wilcox.test",
             remove.bracket = T)+
    labs(x = "Severe Thrombopenia", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
    scale_color_manual(values = c("lightblue", "#17becf"))+
    scale_fill_manual(values = c("lightblue", "#17becf"))+
    scale_x_discrete(labels = c("No", "Yes")) +
    guides(color = "none", fill = "none")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )
  
  
  # Add the plot to the list
  p_sevplt_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_SevPLT_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_sevplt_blborder10cn12_post[[variable]], width = 2, height = 4)
}



### G-CSF administration


# List to save plots for grid
p_gcsf_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  p <- survival_bl_border10_cn12_post|>
    mutate(gcsf_group = ifelse(is.na(First.Day.of.G.CSF), "No", "Yes"))|>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area == variable)|>
    ggplot(aes(x=as.factor(gcsf_group), y=area_perc, fill=as.factor(gcsf_group), color = as.factor(gcsf_group))) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    geom_pwc(ref.group = "No",
             label = "p.format", 
             method = "wilcox.test",
             remove.bracket = F)+
    labs(x = "G-CSF administration", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
    scale_color_manual(values = c("#ADD8E6", "#2A52BE"))+
    scale_fill_manual(values = c("#ADD8E6", "#2A52BE"))+
    guides(color = "none", fill = "none")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )
  
  
  # Add the plot to the list
  p_gcsf_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_GCSF_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_gcsf_blborder10cn12_post[[variable]], width = 2, height = 4)
}


### Hemoglobin


### Analysis of hemoglobin for each neighborhood over time

HgB_time_labels <- c(
  "HgB_d90_grade_grouped" = "Day 90",
  "HgB_d60_grade_grouped" = "Day 60",
  "HgB_d30_grade_grouped" = "Day 30"
)

HgB_time_value_labels <- c(
  "0-1" = "0-1",
  ">=2" = "≥2")


# List to save plots for grid
p_hgb_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  # Generate the survival plot
  p <- survival_bl_border10_cn12_post|>
    mutate(HgB_d30_grade_grouped = case_when(
      HgB_d30_grade %in% 0:1 ~ "0-1",
      HgB_d30_grade > 1 ~ ">=2",
      TRUE ~ NA_character_  ),
      HgB_d60_grade_grouped = case_when(
        HgB_d60_grade %in% 0:1 ~ "0-1",
        HgB_d60_grade > 1 ~ ">=2",
        TRUE ~ NA_character_  ),
      HgB_d90_grade_grouped = case_when(
        HgB_d90_grade %in% 0:1 ~ "0-1",
        HgB_d90_grade > 1 ~ ">=2",
        TRUE ~ NA_character_  )) |>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area != "nan")|>
    mutate(HgB_d90_grade_grouped = factor(HgB_d90_grade_grouped, levels = c("0-1", ">=2")),
           HgB_d60_grade_grouped = factor(HgB_d60_grade_grouped, levels = c("0-1", ">=2")),
           HgB_d30_grade_grouped = factor(HgB_d30_grade_grouped, levels = c("0-1", ">=2"))) |>
    pivot_longer(cols = c(HgB_d90_grade_grouped, HgB_d60_grade_grouped, HgB_d30_grade_grouped), names_to = "HgB_time", values_to = "HgB_time_value")|>
    filter(area == variable) |>
    filter(!is.na(HgB_time_value))|>
    ggplot(aes(x=HgB_time_value, y=area_perc, fill=HgB_time_value, color = HgB_time_value)) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    facet_wrap(~HgB_time, labeller = labeller(HgB_time = HgB_time_labels)) + 
    geom_pwc(ref.group = "0-1",
             label = "p.adj.format", 
             method = "wilcox.test",
             p.adjust.method = "fdr")+
    labs(x = "Anemia [CTCAE grade]", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
    scale_color_manual(values = c('#DDA0DD', "#BA55D3"))+
    scale_fill_manual(values = c('#DDA0DD', "#BA55D3"))+
    scale_x_discrete(labels = HgB_time_value_labels) +
    guides(color = "none", fill = "none")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )
  
  #Adjust multiple testing
  # p <- ggadjust_pvalue(p, p.adjust.method = "fdr")
  
  # Add the plot to the list
  p_hgb_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_HgB_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_hgb_blborder10cn12_post[[variable]], width = 4, height = 4)
}


survival_bl_border10_cn12_post|>
  mutate(HgB_d30_grade_grouped = case_when(
    HgB_d30_grade %in% 0:1 ~ "0-1",
    HgB_d30_grade > 1 ~ ">=2",
    TRUE ~ NA_character_  ),
    HgB_d60_grade_grouped = case_when(
      HgB_d60_grade %in% 0:1 ~ "0-1",
      HgB_d60_grade > 1 ~ ">=2",
      TRUE ~ NA_character_  ),
    HgB_d90_grade_grouped = case_when(
      HgB_d90_grade %in% 0:1 ~ "0-1",
      HgB_d90_grade > 1 ~ ">=2",
      TRUE ~ NA_character_  )) |>
  pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
  filter(area != "nan")|>
  mutate(HgB_d90_grade_grouped = factor(HgB_d90_grade_grouped, levels = c("0-1", ">=2")),
         HgB_d60_grade_grouped = factor(HgB_d60_grade_grouped, levels = c("0-1", ">=2")),
         HgB_d30_grade_grouped = factor(HgB_d30_grade_grouped, levels = c("0-1", ">=2"))) |>
  pivot_longer(cols = c(HgB_d90_grade_grouped, HgB_d60_grade_grouped, HgB_d30_grade_grouped), names_to = "HgB_time", values_to = "HgB_time_value")|>
  filter(!is.na(HgB_time_value))|>
group_by(HgB_time, HgB_time_value, area) |>
  summarise(median=median(area_perc))|>
  filter(area == "CN11")


### Severe Anemia

# List to save plots for grid
p_sevhgb_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  p <- survival_bl_border10_cn12_post|>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area == variable)|>
    ggplot(aes(x=as.factor(Severe.Anemia..Hb...8.g.dL.or.requiring.pRBC.tx..Day.0.100), y=area_perc, fill=as.factor(Severe.Anemia..Hb...8.g.dL.or.requiring.pRBC.tx..Day.0.100), color = as.factor(Severe.Anemia..Hb...8.g.dL.or.requiring.pRBC.tx..Day.0.100))) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    geom_pwc(ref.group = "0",
             label = "p.format", 
             method = "wilcox.test",
             remove.bracket = F)+
    labs(x = "Severe Anemia", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
    scale_color_manual(values = c('#DDA0DD', "#BA55D3"))+
    scale_fill_manual(values = c('#DDA0DD', "#BA55D3"))+
    scale_x_discrete(labels = c("No", "Yes")) +
    guides(color = "none", fill = "none")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )
  
  
  # Add the plot to the list
  p_sevhgb_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_SevHgB_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_sevhgb_blborder10cn12_post[[variable]], width = 2, height = 4)
}



### Severe Anemia

# List to save plots for grid
p_hgbtx_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  p <- survival_bl_border10_cn12_post|>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area == variable)|>
    ggplot(aes(x=as.factor(prbc_d0_d100), y=area_perc, fill=as.factor(prbc_d0_d100), color = as.factor(prbc_d0_d100))) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    geom_pwc(ref.group = "0",
             label = "p.format", 
             method = "wilcox.test",
             remove.bracket = F)+
    labs(x = "RBC transfusion", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
    scale_color_manual(values = c('#DDA0DD', "#BA55D3"))+
    scale_fill_manual(values = c('#DDA0DD', "#BA55D3"))+
    scale_x_discrete(labels = c("No", "Yes")) +
    guides(color = "none", fill = "none")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )
  
  
  # Add the plot to the list
  p_hgbtx_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_HgBTx_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_hgbtx_blborder10cn12_post[[variable]], width = 2, height = 4)
}





survival_bl_border10_cn12_post|>
  pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
  filter(area == "CN4")|>
  ggplot(aes(x=as.factor(late_icaht), y=area_perc)) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(label = "p.format",
           method = "t.test")+
  labs(x = "RBC transfusion", y = "Area [norm.]") +
  guides(fill = "none", color = "none")+
  scale_color_manual(values = c('#DDA0DD', "#BA55D3"))+
  scale_fill_manual(values = c('#DDA0DD', "#BA55D3"))+
  scale_x_discrete(labels = ANC_time_value_labels) +
  guides(color = "none", fill = "none")+
  theme_classic()+
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
  )


##
## EARLY ICAHT


# List to save plots for grid
p_earlyicaht_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  p <- survival_bl_border10_cn12_post|>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area == variable)|>
    ggplot(aes(x=as.factor(early_icaht), y=area_perc)) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    geom_pwc(ref.group = "0",
             label = "p.format", 
             method = "wilcox.test",
             remove.bracket = T)+
    labs(x = "Early ICAHT", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
    scale_x_discrete(labels = c("No", "Yes")) +
    guides(color = "none", fill = "none")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )
  
  
  # Add the plot to the list
  p_earlyicaht_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_earlyICAHT_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_earlyicaht_blborder10cn12_post[[variable]], width = 2, height = 4)
}

# List to save plots for grid
p_lateicaht_blborder10cn12_post <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_post_hemato"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_post_variables_wonan) {
  
  p <- survival_bl_border10_cn12_post|>
    mutate(late_icaht_group = ifelse(late_icaht %in% c(0,1), "low", "high"))|>
    pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
    filter(area == variable)|>
    ggplot(aes(x=as.factor(late_icaht_group), y=area_perc)) +
    geom_boxplot(alpha = 0.2) +
    geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
    geom_pwc(label = "p.format", 
             method = "wilcox.test")+
    labs(x = "Early ICAHT", y = "Area [norm.]") +
    guides(fill = "none", color = "none")+
   # scale_x_discrete(limits = c("low", "high")) +
    guides(color = "none", fill = "none")+
    theme_classic()+
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
    )
  
  
  # Add the plot to the list
  p_lateicaht_blborder10cn12_post[[variable]] <- p
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_lateICAHT_post", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = p_lateicaht_blborder10cn12_post[[variable]], width = 2, height = 4)
}



# 
# 
# 
# survival_bl_border10_cn12_post|>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   ggplot(aes(x=as.factor(plt_d0_d100), y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc()+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# 
# p_hemato_border10cn12_post_HGBd30_1 <-survival_bl_border10_cn12_post|>
#   mutate(HgB_d30_grade_grouped = case_when(
#     HgB_d30_grade == 0 ~ "0",
#     HgB_d30_grade %in% 1:2 ~ "1-2",
#     HgB_d30_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(HgB_d30_grade_grouped = factor(HgB_d30_grade_grouped, levels = c("0", "1-2", ">3"))) |>
#   ggplot(aes(x=HgB_d30_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc(method = "wilcox_test", ref.group = "0", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_HGBd30_2 <-survival_bl_border10_cn12_post|>
#   mutate(HgB_d30_grade_grouped = case_when(
#     HgB_d30_grade %in% 0:1 ~ "0-1",
#     HgB_d30_grade %in% 2 ~ "2",
#     HgB_d30_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(HgB_d30_grade_grouped = factor(HgB_d30_grade_grouped, levels = c("0-1", "2", ">3"))) |>
#   ggplot(aes(x=HgB_d30_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc(method = "wilcox_test", ref.group = "0-1", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_HGBd30_3 <-survival_bl_border10_cn12_post|>
#   mutate(HgB_d30_grade_grouped = case_when(
#     HgB_d30_grade %in% 0:1 ~ "0-1",
#     HgB_d30_grade > 1 ~ ">=2",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(HgB_d30_grade_grouped = factor(HgB_d30_grade_grouped, levels = c("0-1", ">=2"))) |>
#   ggplot(aes(x=HgB_d30_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# 
# p_hemato_border10cn12_post_HGBd60_1 <-survival_bl_border10_cn12_post|>
#   mutate(HgB_d60_grade_grouped = case_when(
#     HgB_d60_grade == 0 ~ "0",
#     HgB_d60_grade %in% 1:2 ~ "1-2",
#     HgB_d60_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(HgB_d60_grade_grouped = factor(HgB_d60_grade_grouped, levels = c("0", "1-2", ">3"))) |>
#   ggplot(aes(x=HgB_d60_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc(method = "wilcox_test", ref.group = "0", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D60")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_HGBd60_2 <-survival_bl_border10_cn12_post|>
#   mutate(HgB_d60_grade_grouped = case_when(
#     HgB_d60_grade %in% 0:1 ~ "0-1",
#     HgB_d60_grade %in% 2 ~ "2",
#     HgB_d60_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(HgB_d60_grade_grouped = factor(HgB_d60_grade_grouped, levels = c("0-1", "2", ">3"))) |>
#   ggplot(aes(x=HgB_d60_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc(method = "wilcox_test", ref.group = "0-1", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D60")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_HGBd60_3 <-survival_bl_border10_cn12_post|>
#   mutate(HgB_d60_grade_grouped = case_when(
#     HgB_d60_grade %in% 0:1 ~ "0-1",
#     HgB_d60_grade > 1 ~ ">=2",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(HgB_d60_grade_grouped = factor(HgB_d60_grade_grouped, levels = c("0-1", ">=2"))) |>
#   ggplot(aes(x=HgB_d60_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D60")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# 
# 
# p_hemato_border10cn12_post_HGBd90_1 <-survival_bl_border10_cn12_post|>
#   mutate(HgB_d90_grade_grouped = case_when(
#     HgB_d90_grade == 0 ~ "0",
#     HgB_d90_grade %in% 1:2 ~ "1-2",
#     HgB_d90_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(HgB_d90_grade_grouped = factor(HgB_d90_grade_grouped, levels = c("0", "1-2", ">3"))) |>
#   ggplot(aes(x=HgB_d90_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc(method = "wilcox_test", ref.group = "0", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D90")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_HGBd90_2 <-survival_bl_border10_cn12_post|>
#   mutate(HgB_d90_grade_grouped = case_when(
#     HgB_d90_grade %in% 0:1 ~ "0-1",
#     HgB_d90_grade %in% 2 ~ "2",
#     HgB_d90_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(HgB_d90_grade_grouped = factor(HgB_d90_grade_grouped, levels = c("0-1", "2", ">3"))) |>
#   ggplot(aes(x=HgB_d90_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc(method = "wilcox_test", ref.group = "0-1", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D90")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_HGBd90_3 <-survival_bl_border10_cn12_post|>
#   mutate(HgB_d90_grade_grouped = case_when(
#     HgB_d90_grade %in% 0:1 ~ "0-1",
#     HgB_d90_grade > 1 ~ ">=2",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(HgB_d90_grade_grouped = factor(HgB_d90_grade_grouped, levels = c("0-1", ">=2"))) |>
#   ggplot(aes(x=HgB_d90_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = '#DDA0DD')+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Anemia CTCAE", y = "Area [norm.]", title = "D90")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# 
# 
# 
# p_hemato_border10cn12_post_PLTd30_1 <-survival_bl_border10_cn12_post|>
#   mutate(PLT_d30_grade_grouped = case_when(
#     PLT_d30_grade == 0 ~ "0",
#     PLT_d30_grade %in% 1:2 ~ "1-2",
#     PLT_d30_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(PLT_d30_grade_grouped = factor(PLT_d30_grade_grouped, levels = c("0", "1-2", ">3"))) |>
#   ggplot(aes(x=PLT_d30_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = "#17becf")+
#   geom_pwc(method = "wilcox_test", ref.group = "0", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Thrombocytopenia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_PLTd30_2 <-survival_bl_border10_cn12_post|>
#   mutate(PLT_d30_grade_grouped = case_when(
#     PLT_d30_grade %in% 0:1 ~ "0-1",
#     PLT_d30_grade %in% 2 ~ "2",
#     PLT_d30_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(PLT_d30_grade_grouped = factor(PLT_d30_grade_grouped, levels = c("0-1", "2", ">3"))) |>
#   ggplot(aes(x=PLT_d30_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = "#17becf")+
#   geom_pwc(method = "wilcox_test", ref.group = "0-1", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Thrombocytopenia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_PLTd30_3 <- survival_bl_border10_cn12_post|>
#   mutate(PLT_d30_grade_grouped = case_when(
#     PLT_d30_grade %in% 0:1 ~ "0-1",
#     PLT_d30_grade > 1 ~ ">=2",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(PLT_d30_grade_grouped = factor(PLT_d30_grade_grouped, levels = c("0-1", ">=2"))) |>
#   ggplot(aes(x=PLT_d30_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = "#17becf")+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Thrombocytopenia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# 
# p_hemato_border10cn12_post_PLTd60_1 <- survival_bl_border10_cn12_post|>
#   mutate(PLT_d60_grade_grouped = case_when(
#     PLT_d60_grade == 0 ~ "0",
#     PLT_d60_grade %in% 1:2 ~ "1-2",
#     PLT_d60_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(PLT_d60_grade_grouped = factor(PLT_d60_grade_grouped, levels = c("0", "1-2", ">3"))) |>
#   ggplot(aes(x=PLT_d60_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = "#17becf")+
#   geom_pwc(method = "wilcox_test", ref.group = "0", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Thrombocytopenia CTCAE", y = "Area [norm.]", title = "D60")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_PLTd60_2 <- survival_bl_border10_cn12_post|>
#   mutate(PLT_d60_grade_grouped = case_when(
#     PLT_d60_grade %in% 0:1 ~ "0-1",
#     PLT_d60_grade %in% 2 ~ "2",
#     PLT_d60_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(PLT_d60_grade_grouped = factor(PLT_d60_grade_grouped, levels = c("0-1", "2", ">3"))) |>
#   ggplot(aes(x=PLT_d60_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = "#17becf")+
#   geom_pwc(method = "wilcox_test", ref.group = "0-1", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Thrombocytopenia CTCAE", y = "Area [norm.]", title = "D60")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_PLTd60_3 <- survival_bl_border10_cn12_post|>
#   mutate(PLT_d60_grade_grouped = case_when(
#     PLT_d60_grade %in% 0:1 ~ "0-1",
#     PLT_d60_grade > 1 ~ ">=2",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(PLT_d60_grade_grouped = factor(PLT_d60_grade_grouped, levels = c("0-1", ">=2"))) |>
#   ggplot(aes(x=PLT_d60_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = "#17becf")+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Thrombocytopenia CTCAE", y = "Area [norm.]", title = "D60")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# 
# 
# p_hemato_border10cn12_post_PLTd90_1 <-survival_bl_border10_cn12_post|>
#   mutate(PLT_d90_grade_grouped = case_when(
#     PLT_d90_grade == 0 ~ "0",
#     PLT_d90_grade %in% 1:2 ~ "1-2",
#     PLT_d90_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(PLT_d90_grade_grouped = factor(PLT_d90_grade_grouped, levels = c("0", "1-2", ">3"))) |>
#   ggplot(aes(x=PLT_d90_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = "#17becf")+
#   geom_pwc(method = "wilcox_test", ref.group = "0", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Thrombocytopenia CTCAE", y = "Area [norm.]", title = "D90")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_PLTd90_2 <-survival_bl_border10_cn12_post|>
#   mutate(PLT_d90_grade_grouped = case_when(
#     PLT_d90_grade %in% 0:1 ~ "0-1",
#     PLT_d90_grade %in% 2 ~ "2",
#     PLT_d90_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(PLT_d90_grade_grouped = factor(PLT_d90_grade_grouped, levels = c("0-1", "2", ">3"))) |>
#   ggplot(aes(x=PLT_d90_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = "#17becf")+
#   geom_pwc(method = "wilcox_test", ref.group = "0-1", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Thrombocytopenia CTCAE", y = "Area [norm.]", title = "D90")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_PLTd90_3 <- survival_bl_border10_cn12_post|>
#   mutate(PLT_d90_grade_grouped = case_when(
#     PLT_d90_grade %in% 0:1 ~ "0-1",
#     PLT_d90_grade > 1 ~ ">=2",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(PLT_d90_grade_grouped = factor(PLT_d90_grade_grouped, levels = c("0-1", ">=2"))) |>
#   ggplot(aes(x=PLT_d90_grade_grouped, y=area_perc))+
#   geom_boxplot(alpha = 0.2) +
#   geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Thrombocytopenia CTCAE", y = "Area [norm.]", title = "D90")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# 
# 
# 
# 
# p_hemato_border10cn12_post_ANCd30_1 <- survival_bl_border10_cn12_post|>
#   mutate(ANC_d30_grade_grouped = case_when(
#     ANC_d30_grade == 0 ~ "0",
#     ANC_d30_grade %in% 1:2 ~ "1-2",
#     ANC_d30_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(ANC_d30_grade_grouped = factor(ANC_d30_grade_grouped, levels = c("0", "1-2", ">3"))) |>
#   ggplot(aes(x=ANC_d30_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkblue')+
#   geom_pwc(method = "wilcox_test", ref.group = "0", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Neutropenia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_ANCd30_2 <- survival_bl_border10_cn12_post|>
#   mutate(ANC_d30_grade_grouped = case_when(
#     ANC_d30_grade %in% 0:1 ~ "0-1",
#     ANC_d30_grade %in% 2 ~ "2",
#     ANC_d30_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(ANC_d30_grade_grouped = factor(ANC_d30_grade_grouped, levels = c("0-1", "2", ">3"))) |>
#   ggplot(aes(x=ANC_d30_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkblue')+
#   geom_pwc(method = "wilcox_test", ref.group = "0-1", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Neutropenia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_ANCd30_3 <- survival_bl_border10_cn12_post|>
#   mutate(ANC_d30_grade_grouped = case_when(
#     ANC_d30_grade %in% 0:1 ~ "0-1",
#     ANC_d30_grade > 1 ~ ">=2",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(ANC_d30_grade_grouped = factor(ANC_d30_grade_grouped, levels = c("0-1", ">=2"))) |>
#   ggplot(aes(x=ANC_d30_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkblue')+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Neutropenia CTCAE", y = "Area [norm.]", title = "D30")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# 
# p_hemato_border10cn12_post_ANCd60_1 <- survival_bl_border10_cn12_post|>
#   mutate(ANC_d60_grade_grouped = case_when(
#     ANC_d60_grade == 0 ~ "0",
#     ANC_d60_grade %in% 1:2 ~ "1-2",
#     ANC_d60_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(ANC_d60_grade_grouped = factor(ANC_d60_grade_grouped, levels = c("0", "1-2", ">3"))) |>
#   ggplot(aes(x=ANC_d60_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkblue')+
#   geom_pwc(method = "wilcox_test", ref.group = "0", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Neutropenia CTCAE", y = "Area [norm.]", title = "D60")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_ANCd60_2 <- survival_bl_border10_cn12_post|>
#   mutate(ANC_d60_grade_grouped = case_when(
#     ANC_d60_grade %in% 0:1 ~ "0-1",
#     ANC_d60_grade %in% 2 ~ "2",
#     ANC_d60_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(ANC_d60_grade_grouped = factor(ANC_d60_grade_grouped, levels = c("0-1", "2", ">3"))) |>
#   ggplot(aes(x=ANC_d60_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkblue')+
#   geom_pwc(method = "wilcox_test", ref.group = "0-1", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Neutropenia CTCAE", y = "Area [norm.]", title = "D60")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_ANCd60_3 <- survival_bl_border10_cn12_post|>
#   mutate(ANC_d60_grade_grouped = case_when(
#     ANC_d60_grade %in% 0:1 ~ "0-1",
#     ANC_d60_grade > 1 ~ ">=2",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(ANC_d60_grade_grouped = factor(ANC_d60_grade_grouped, levels = c("0-1", ">=2"))) |>
#   ggplot(aes(x=ANC_d60_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkblue')+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Neutropenia CTCAE", y = "Area [norm.]", title = "D60")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# 
# 
# p_hemato_border10cn12_post_ANCd90_1 <- survival_bl_border10_cn12_post|>
#   mutate(ANC_d90_grade_grouped = case_when(
#     ANC_d90_grade == 0 ~ "0",
#     ANC_d90_grade %in% 1:2 ~ "1-2",
#     ANC_d90_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(ANC_d90_grade_grouped = factor(ANC_d90_grade_grouped, levels = c("0", "1-2", ">3"))) |>
#   ggplot(aes(x=ANC_d90_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkblue')+
#   geom_pwc(method = "wilcox_test", ref.group = "0", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Neutropenia CTCAE", y = "Area [norm.]", title = "D90")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_ANCd90_2 <- survival_bl_border10_cn12_post|>
#   mutate(ANC_d90_grade_grouped = case_when(
#     ANC_d90_grade %in% 0:1 ~ "0-1",
#     ANC_d90_grade %in% 2 ~ "2",
#     ANC_d90_grade > 2 ~ ">3",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(ANC_d90_grade_grouped = factor(ANC_d90_grade_grouped, levels = c("0-1", "2", ">3"))) |>
#   ggplot(aes(x=ANC_d90_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkblue')+
#   geom_pwc(method = "wilcox_test", ref.group = "0-1", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Neutropenia CTCAE", y = "Area [norm.]", title = "D90")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# p_hemato_border10cn12_post_ANCd90_3 <- survival_bl_border10_cn12_post|>
#   mutate(ANC_d90_grade_grouped = case_when(
#     ANC_d90_grade %in% 0:1 ~ "0-1",
#     ANC_d90_grade > 1 ~ ">=2",
#     TRUE ~ NA_character_  )) |>
#   pivot_longer(cols = all_of(bl_border10_cn12_post_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   mutate(ANC_d90_grade_grouped = factor(ANC_d90_grade_grouped, levels = c("0-1", ">=2"))) |>
#   ggplot(aes(x=ANC_d90_grade_grouped, y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkblue')+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "Neutropenia CTCAE", y = "Area [norm.]", title = "D90")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic() +
#   theme(
#     strip.background = element_blank(),  # Remove the box around facet labels
#     strip.text = element_text(size = 10),    plot.title = element_text(size = 12, hjust = 0.5)  # Adjust the font size of facet labels if needed
#   )
# 
# ##
# ## Saving HEMATOTOX plots
# ##
# 
# # Example usage
# pattern_hemato_border10cn12_post <- "p_hemato_border10cn12_post_[a-zA-Z0-9]+_[0-9]"
# plots_hemato_border10cn12_post <- mget(ls(pattern = pattern_hemato_border10cn12_post))
# 
# # Define the directory to save the CRS plots
# output_dir <- "border10cn12_post_hemato"
# dir.create(output_dir, showWarnings = FALSE)
# 
# # Save the plots to the specified directory
# save_plot(plots_hemato_border10cn12_post, output_dir, width =6, height =8)
# 
# 
# 
# ##
# ##
# ##
# 
# survival_bl_border10_cn12|>
#   pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc")|>
#   filter(area != "nan")|>
#   ggplot(aes(x=factor(hrcyto, levels = c("SR", "HR")), y=area_perc))+
#   geom_boxplot()+
#   geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkgrey')+
#   geom_pwc(method = "wilcox_test", 
#            bracket.nudge.y = 0.05, 
#            step.increase = 0.15, 
#            tip.length = 0.01, label.size = 3)+
#   labs(x = "CAR-HEMATOTOX", y = "Area [norm.]", title = "Hematotoxicity")+
#   facet_wrap(vars(area), scales = "free")+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
#   theme_classic()+
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_text(size = 10),    
#     plot.title = element_text(size = 12, hjust = 0.5)
#   )
# 
# # CN0, CN11, CN12, CN15, CN7 associated with hematotoxicity score
# 
# survival_bl_border10_cn12_post<-survival_bl_border10_cn12_post |>
#   dplyr::mutate(X0_class = ifelse('X0' > median('X0'), 'high', 'low'))|>
#   pivot_longer(cols = c("PLT_LD", "PLT_d0", "PLT_d7", "PLT_d30",
#                         "PLT_d60", "PLT_d90"), names_to = "PLT_time", values_to = "PLT") |>
#   mutate(PLT_time = factor(PLT_time, levels = c("PLT_LD", "PLT_d0", "PLT_d7", 
#                                                 "PLT_d30", "PLT_d60", "PLT_d90"))) |>
#   ggplot(aes(x = PLT_time, y = PLT, color = X0_class, group = X0_class)) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +  # Jittered points
#   geom_smooth(aes(fill = X0_class), method = "loess", se = TRUE, alpha = 0.2) +  # Loess smoothing with ribbons
#   scale_color_manual(values = c("#FF4500", "darkgrey")) +
#   scale_fill_manual(values = c("#FF4500", "lightgrey")) +  # Matching fill colors to point colors
#   theme_classic() +
#   labs(x = "PLT Time", y = "PLT", title = "PLT Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# survival_bl_border10_cn12_post |>
#   mutate(X0_class = ifelse(X0 >= median(X0, na.rm = TRUE), 'high', 'low'))|>
#   pivot_longer(cols = c("ANC_LD", "ANC_d0", "ANC_d7", "ANC_d30",
#                         "ANC_d60", "ANC_d90"), names_to = "ANC_time", values_to = "ANC") |>
#   mutate(ANC_time = factor(ANC_time, levels = c("ANC_LD", "ANC_d0", "ANC_d7", 
#                                                 "ANC_d30", "ANC_d60", "ANC_d90")))|>
#   ggplot( aes(x = ANC_time, y = ANC, color = X0_class, group = X0_class)) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +  # Jittered points
#   geom_smooth(aes(fill = X0_class), method = "loess", se = TRUE, alpha = 0.2) +  # Loess smoothing with ribbons
#   scale_color_manual(values = c("#8A2BE2", "darkgrey")) +  # Shades of violet
#   scale_fill_manual(values = c("#8A2BE2", "lightgrey")) +  # Matching fill colors to point colors
#   ylim(0,5)+
#   theme_classic() +
#   labs(x = "ANC Time", y = "ANC", title = "ANC Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# survival_bl_border10_cn12 |>
#   mutate(X0_class = ifelse(X0 > median(X0, na.rm = TRUE), 'high', 'low'))|>
#   pivot_longer(cols = c("HgB_LD", "HgB_d0", "HgB_d7", "HgB_d30",
#                         "HgB_d60", "HgB_d90"), names_to = "HgB_time", values_to = "HgB") |>
#   mutate(HgB_time = factor(HgB_time, levels = c("HgB_LD", "HgB_d0", "HgB_d7", 
#                                                 "HgB_d30", "HgB_d60", "HgB_d90")))|>
#   ggplot(aes(x = HgB_time, y = HgB, color = X0_class, group = X0_class)) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +  # Jittered points
#   geom_smooth(aes(fill = X0_class), method = "loess", se = TRUE, alpha = 0.3) +  # Loess smoothing with ribbons
#   scale_color_manual(values = c("#20B2AA", "darkgrey")) +  # Shades of turquoise
#   scale_fill_manual(values = c("#20B2AA", "lightgrey")) +  # Matching fill colors to point colors
#   theme_classic() +
#   labs(x = "HgB Time", y = "HgB", title = "HgB Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# survival_bl_border10_cn12 |>
#   mutate(X11_class = ifelse(X11 > median(X11, na.rm = TRUE), 'high', 'low'))|>
#   pivot_longer(cols = c("PLT_LD", "PLT_d0", "PLT_d7", "PLT_d30",
#                         "PLT_d60", "PLT_d90"), names_to = "PLT_time", values_to = "PLT") |>
#   mutate(PLT_time = factor(PLT_time, levels = c("PLT_LD", "PLT_d0", "PLT_d7", 
#                                                 "PLT_d30", "PLT_d60", "PLT_d90"))) |>
#   ggplot(aes(x = PLT_time, y = PLT, color = X11_class, group = X11_class)) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +  # Jittered points
#   geom_smooth(aes(fill = X11_class), method = "loess", se = TRUE, alpha = 0.2) +  # Loess smoothing with ribbons
#   scale_color_manual(values = c("#FF4500", "darkgrey")) +
#   scale_fill_manual(values = c("#FF4500", "lightgrey")) +  # Matching fill colors to point colors
#   theme_classic() +
#   labs(x = "PLT Time", y = "PLT", title = "PLT Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# survival_bl_border10_cn12 |>
#   mutate(X11_class = ifelse(X11 > median(X11, na.rm = TRUE), 'high', 'low'))|>
#   pivot_longer(cols = c("ANC_LD", "ANC_d0", "ANC_d7", "ANC_d30",
#                         "ANC_d60", "ANC_d90"), names_to = "ANC_time", values_to = "ANC") |>
#   mutate(ANC_time = factor(ANC_time, levels = c("ANC_LD", "ANC_d0", "ANC_d7", 
#                                                 "ANC_d30", "ANC_d60", "ANC_d90")))|>
#   ggplot( aes(x = ANC_time, y = ANC, color = X11_class, group = X11_class)) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +  # Jittered points
#   geom_smooth(aes(fill = X11_class), method = "loess", se = TRUE, alpha = 0.2) +  # Loess smoothing with ribbons
#   scale_color_manual(values = c("#8A2BE2", "darkgrey")) +  # Shades of violet
#   scale_fill_manual(values = c("#8A2BE2", "lightgrey")) +  # Matching fill colors to point colors
#   ylim(0,5)+
#   theme_classic() +
#   labs(x = "ANC Time", y = "ANC", title = "ANC Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# survival_bl_border10_cn12 |>
#   mutate(X11_class = ifelse(X11 > median(X11, na.rm = TRUE), 'high', 'low'))|>
#   pivot_longer(cols = c("HgB_LD", "HgB_d0", "HgB_d7", "HgB_d30",
#                         "HgB_d60", "HgB_d90"), names_to = "HgB_time", values_to = "HgB") |>
#   mutate(HgB_time = factor(HgB_time, levels = c("HgB_LD", "HgB_d0", "HgB_d7", 
#                                                 "HgB_d30", "HgB_d60", "HgB_d90")))|>
#   ggplot(aes(x = HgB_time, y = HgB, color = X11_class, group = X11_class)) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +  # Jittered points
#   geom_smooth(aes(fill = X11_class), method = "loess", se = TRUE, alpha = 0.3) +  # Loess smoothing with ribbons
#   scale_color_manual(values = c("#20B2AA", "darkgrey")) +  # Shades of turquoise
#   scale_fill_manual(values = c("#20B2AA", "lightgrey")) +  # Matching fill colors to point colors
#   theme_classic() +
#   labs(x = "HgB Time", y = "HgB", title = "HgB Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# 
# survival_bl_border10_cn12 |>
#   mutate(X12_class = ifelse(X12 > median(X12, na.rm = TRUE), 'high', 'low'))|>
#   pivot_longer(cols = c("PLT_LD", "PLT_d0", "PLT_d7", "PLT_d30",
#                         "PLT_d60", "PLT_d90"), names_to = "PLT_time", values_to = "PLT") |>
#   mutate(PLT_time = factor(PLT_time, levels = c("PLT_LD", "PLT_d0", "PLT_d7", 
#                                                 "PLT_d30", "PLT_d60", "PLT_d90"))) |>
#   ggplot(aes(x = as.factor(PLT_time), y = PLT, color = X12_class, group = X12_class)) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +  # Jittered points
#   geom_smooth(aes(fill = X12_class), method = "loess", se = TRUE, alpha = 0.2) +  # Loess smoothing with ribbons
#   scale_color_manual(values = c("#FF4500", "darkgrey")) +
#   scale_fill_manual(values = c("#FF4500", "lightgrey")) +  # Matching fill colors to point colors
#   theme_classic() +
#   labs(x = "PLT Time", y = "PLT", title = "PLT Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# survival_bl_border10_cn12 |>
#   mutate(X11_class = ifelse(X11 > median(X11, na.rm = TRUE), 'high', 'low')) |>
#   pivot_longer(cols = c("PLT_LD", "PLT_d0", "PLT_d7", "PLT_d30",
#                         "PLT_d60", "PLT_d90"), names_to = "PLT_time", values_to = "PLT") |>
#   mutate(PLT_time = factor(PLT_time, levels = c("PLT_LD", "PLT_d0", "PLT_d7", 
#                                                 "PLT_d30", "PLT_d60", "PLT_d90"))) |>
#   ggplot(aes(x = PLT_time, y = PLT, fill = X11_class)) +
#   geom_boxplot(position = position_dodge(0.75), alpha = 0.7) +  # Grouped boxplots
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.3, size = 1) +  # Jittered points
#   scale_fill_manual(values = c("#FF4500", "darkgrey")) +  # Matching fill colors to point colors
#   theme_classic() +
#   labs(x = "PLT Time", y = "PLT", title = "PLT Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   geom_pwc(aes(group = X11_class), method = "wilcox_test",
#            step.increase = 0.1,
#            label.size = 3,
#            tip.length = 0.01,
#            bracket.nudge.y = 0.05)
# 
# 
# 
# survival_bl_border10_cn12 |>
#   mutate(X12_class = ifelse(X12 > median(X12, na.rm = TRUE), 'high', 'low'))|>
#   pivot_longer(cols = c("ANC_LD", "ANC_d0", "ANC_d7", "ANC_d30",
#                         "ANC_d60", "ANC_d90"), names_to = "ANC_time", values_to = "ANC") |>
#   mutate(ANC_time = factor(ANC_time, levels = c("ANC_LD", "ANC_d0", "ANC_d7", 
#                                                 "ANC_d30", "ANC_d60", "ANC_d90")))|>
#   ggplot( aes(x = ANC_time, y = ANC, color = X12_class, group = X12_class)) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +  # Jittered points
#   geom_smooth(aes(fill = X12_class), method = "loess", se = TRUE, alpha = 0.2) +  # Loess smoothing with ribbons
#   scale_color_manual(values = c("#8A2BE2", "darkgrey")) +  # Shades of violet
#   scale_fill_manual(values = c("#8A2BE2", "lightgrey")) +  # Matching fill colors to point colors
#   ylim(0,5)+
#   theme_classic() +
#   labs(x = "ANC Time", y = "ANC", title = "ANC Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# survival_bl_border10_cn12 |>
#   mutate(X12_class = ifelse(X12 > median(X12, na.rm = TRUE), 'high', 'low'))|>
#   pivot_longer(cols = c("HgB_LD", "HgB_d0", "HgB_d7", "HgB_d30",
#                         "HgB_d60", "HgB_d90"), names_to = "HgB_time", values_to = "HgB") |>
#   mutate(HgB_time = factor(HgB_time, levels = c("HgB_LD", "HgB_d0", "HgB_d7", 
#                                                 "HgB_d30", "HgB_d60", "HgB_d90")))|>
#   ggplot(aes(x = HgB_time, y = HgB, color = X12_class, group = X12_class)) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +  # Jittered points
#   geom_smooth(aes(fill = X12_class), method = "loess", se = TRUE, alpha = 0.3) +  # Loess smoothing with ribbons
#   scale_color_manual(values = c("#20B2AA", "darkgrey")) +  # Shades of turquoise
#   scale_fill_manual(values = c("#20B2AA", "lightgrey")) +  # Matching fill colors to point colors
#   theme_classic() +
#   labs(x = "HgB Time", y = "HgB", title = "HgB Levels Over Time by Hematotox Status") +
#   theme(plot.title = element_text(hjust = 0.5))
# 

