### Function to save plots ----
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

###
### Preparation of baseline data for further analysis ----
###

## Changing neighborhood annotation
survival_bl_border10_cn12 <- processed_data_border$border10_cn12

survival_bl_border10_cn12 <- survival_bl_border10_cn12 %>%
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

bl_border10_cn12_variables <- c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", 
                                "CN6", "CN7", "CN8", "CN9", "CN10", "CN11", "nan")

bl_border10_cn12_variables_wonan <- c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", 
                                     "CN6", "CN7", "CN8", "CN9", "CN10", "CN11")


## Loading and correctly annotating raw counts matrix
## Loading raw counts matrix
bl_border10cn12_raw <- read.csv("full_results/border/baseline/border-10_CNs-12_raw-results___baseline.csv")

## Displaying cell types per neighborhood
bl_border10cn12_raw <- bl_border10cn12_raw %>%
  mutate(medium_cell_type = str_replace(medium_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         medium_cell_type = str_replace(medium_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"))

bl_border10cn12_raw <- bl_border10cn12_raw %>%
  mutate(minor_cell_type = str_replace(minor_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         minor_cell_type = str_replace(minor_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"))

bl_border10cn12_raw <- bl_border10cn12_raw %>%
  mutate(functional_minor_cell_type = str_replace(functional_minor_cell_type, "M2-like M\\$\\\\phi\\$", "M2-like macrophage"),
         functional_minor_cell_type = str_replace(functional_minor_cell_type, "M1-like M\\$\\\\phi\\$", "M1-like macrophage"))


medium_cell_type_labels <- bl_border10cn12_raw |> select(medium_cell_type) |> unique() |> as.vector()

minor_cell_type_labels <- bl_border10cn12_raw |> select(minor_cell_type) |> unique() |> as.vector()

functional_minor_cell_type_labels <- bl_border10cn12_raw |> select(functional_minor_cell_type) |> unique() |> as.vector()

bl_border10cn12_raw$cn_celltypes <- as.character(bl_border10cn12_raw$cn_celltypes)
#bl_border10_cn12_variables <- data_variables_list_border$border10_cn12
#bl_border10_cn12_variables_wonan <- setdiff(bl_border10_cn12_variables, "nan")

###
### SURVIVAL ANALYSES ----
###

# List to save plots for grid
p_pfs_blborder10cn12 <- list()
p_os_blborder10cn12 <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_survival"
dir.create(output_dir, showWarnings = FALSE)

for (variable in bl_border10_cn12_variables) {
  
  # Directly use the formula in survfit2
  fit <- survfit(Surv(PFS_m, PFS_event) ~ ifelse(get(variable, survival_bl_border10_cn12) > median(get(variable, survival_bl_border10_cn12)), 1, 0), data = survival_bl_border10_cn12)
  
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
                  legend.title = c(paste0("pre",variable)),
                  palette = c("#76AAFF", "#000080")
                  )
  
  # Extract only the plot component from ggsurvplot
  survival_plot <- p$plot
  
  # Add the plot to the list
  p_pfs_blborder10cn12[[variable]] <- survival_plot
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_PFS_CN", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving PFS plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = survival_plot, width = 3, height = 2)
}

# Convert each plot in p_pfs_blborder10cn12 to a grob
grobs_pfs_blborder10cn12 <- lapply(p_pfs_blborder10cn12, ggplotGrob)

# Create the file path for saving the grid image
output_grid_file <- file.path(output_dir, "bl_border10cn12_pfs_grid.svg")

# Save the grid of survival plots as a PNG file
svg(output_grid_file, width = 16, height = 8)  # Adjust resolution and size
do.call(grid.arrange, c(grobs_pfs_blborder10cn12, ncol = 4))  # Arrange plots in a grid, 4 columns
dev.off()


for (variable in bl_border10_cn12_variables) {
  # Directly use the formula in survfit2
  fit <- survfit(Surv(OS_m, OS_event) ~ ifelse(get(variable, survival_bl_border10_cn12) > median(get(variable, survival_bl_border10_cn12)), 1, 0), data = survival_bl_border10_cn12)
  
  # Generate the plot
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
                  legend.title = c(paste0("pre",variable)),
                  palette = c("#88DC88","#009600"))

  
  # Extract only the plot component from ggsurvplot
  survival_plot <- p$plot
  
  # Add the plot to the list
  p_os_blborder10cn12[[variable]] <- survival_plot
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_OS_CN", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving OS plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = survival_plot, width = 3, height = 2)
}

# Convert each plot in p_pfs_blborder10cn12 to a grob (grid graphical object)
grobs_os_blborder10cn12 <- lapply(p_os_blborder10cn12, ggplotGrob)

# Create the file path for saving the grid image
output_grid_file <- file.path(output_dir, "bl_border10cn12_os_grid.svg")

# Save the grid of survival plots as a PNG file
svg(output_grid_file, width = 16, height = 10)  # Adjust resolution and size
do.call(grid.arrange, c(grobs_os_blborder10cn12, ncol = 4))  # Arrange plots in a grid, 4 columns
dev.off()

###
### Calculating survival numbers from survfit functions (median and 95%CI)
###

PFS_CN0 <- survfit(Surv(PFS_days, PFS_event) ~ ifelse(CN0 > median(CN0), 'CN0 high', 'CN0 low'), data = survival_bl_border10_cn12)
PFS_CN0 #high 199 146-363 low:575 (414-NA)

PFS_CN11 <- survfit(Surv(PFS_days, PFS_event) ~ ifelse(CN11 > median(CN11), 'CN11 high', 'CN11 low'), data = survival_bl_border10_cn12)
PFS_CN11 #high 241 154-445 low:575 339-NA

PFS_CN8 <- survfit(Surv(PFS_days, PFS_event) ~ ifelse(CN8 > median(CN8), 'CN8 high', 'CN8 low'), data = survival_bl_border10_cn12)
PFS_CN8 #high 193 126-NA low:469 363-NA

PFS_CN9 <- survfit(Surv(PFS_days, PFS_event) ~ ifelse(CN9 > median(CN9), 'CN9 high', 'CN9 low'), data = survival_bl_border10_cn12)
PFS_CN9 #high 445 363-NA low:199 146-NA

##
## Multivariate COX regression analysis
##

# Convert hrcyto (= High-risk cytogenetics) to a factor
survival_bl_border10_cn12$hrcyto <- as.factor(survival_bl_border10_cn12$hrcyto)

# Releveling hrcyto so that "HR" is the reference level
survival_bl_border10_cn12$hrcyto <- relevel(survival_bl_border10_cn12$hrcyto, ref = "SR")

# Convert 'mycare' to a factor
survival_bl_border10_cn12$mycare <- factor(survival_bl_border10_cn12$mycare)

# Change the reference level for 'mycare'
survival_bl_border10_cn12$mycare <- relevel(survival_bl_border10_cn12$mycare, ref = "low")

summary(coxph(Surv(PFS_days, PFS_event) ~ survival_bl_border10_cn12$CN0
      + survival_bl_border10_cn12$CN11
      + survival_bl_border10_cn12$CN8
      + survival_bl_border10_cn12$CN9 
      + survival_bl_border10_cn12$mycare
      + survival_bl_border10_cn12$CAR.T.Type,
      data = survival_bl_border10_cn12))


border10cn12_pfs_mva_allCN <- coxph(Surv(PFS_days, PFS_event) ~ ifelse(survival_bl_border10_cn12$CN0 > median(survival_bl_border10_cn12$CN0), 1, 0)
                                                              + ifelse(survival_bl_border10_cn12$CN11 > median(survival_bl_border10_cn12$CN11), 1, 0)
                                                              + ifelse(survival_bl_border10_cn12$CN8 > median(survival_bl_border10_cn12$CN8), 1, 0)
                                                              + ifelse(survival_bl_border10_cn12$CN9 > median(survival_bl_border10_cn12$CN9), 1, 0) 
                                                              + survival_bl_border10_cn12$mycare
                                                              + survival_bl_border10_cn12$CAR.T.Type,
                                                              data = survival_bl_border10_cn12)
# Used for single models
# ~ ifelse(survival_bl_border10_cn12$X0 > median(survival_bl_border10_cn12$X0), 1, 0)
# + ifelse(survival_bl_border10_cn12$X11 > median(survival_bl_border10_cn12$X11), 1, 0)
# + ifelse(survival_bl_border10_cn12$X8 > median(survival_bl_border10_cn12$X8), 1, 0)
# + ifelse(survival_bl_border10_cn12$X9 > median(survival_bl_border10_cn12$X9), 1, 0) 

summary_pfs_mva_model_allCN <- summary(border10cn12_pfs_mva_allCN)

# Create a data frame with results from the multivariate COX regression model
results_pfs_mva_model <- data.frame(marker = character(),
                                       HR = numeric(),
                                       lower95 = numeric(),
                                       higher95 = numeric(),
                                       p_value = numeric(),
                                       stringsAsFactors = FALSE)

rownames(summary_pfs_mva_model$conf.int)

results_pfs_mva_model <- data.frame(
  marker = rownames(summary_pfs_mva_model$conf.int),
  HR = summary_pfs_mva_model$conf.int[,"exp(coef)"],
  lower95 = summary_pfs_mva_model$conf.int[, "lower .95"],
  higher95 = summary_pfs_mva_model$conf.int[, "upper .95"],
  p_value = summary_pfs_mva_model$coefficients[, "Pr(>|z|)"]
)

results_pfs_mva_model$marker[1] <- "CN0"
results_pfs_mva_model$marker[2] <- "CN11"
results_pfs_mva_model$marker[3] <- "CN8"
results_pfs_mva_model$marker[4] <- "CN9"
results_pfs_mva_model$marker[5] <- "EMD present"
results_pfs_mva_model$marker[6] <- "logFerritin"
results_pfs_mva_model$marker[7] <- "High-risk cytogenetic"

p_mva_pfs <- ggplot(results_pfs_mva_model)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey")+
  geom_errorbar(aes(x = marker, ymin = lower95, ymax = higher95), width = 0.2, size = 0.6, color = "darkgrey")+
  geom_point(aes(x = marker, y = HR), size = 4, shape = 15, color = "darkblue")+
  coord_flip()+
  scale_y_continuous(limits = c (0, 10), expand = c(0.01,0.5))+
  #scale_x_discrete(limits = c("CN0", "CN11", "CN8", "CN9", "EMD present", "High-risk cytogenetic", "logFerritin"))+
  scale_x_discrete(limits = c("logFerritin", "dummy lab",
                              "High-risk cytogenetic", "dummy cyto",
                              "EMD present", "dummy EMD",
                              "CN9", "CN8", "CN11", "CN0", "dummy CN"))+
  labs(x = "", y = "HR [95% CI]")+
  guides(fill = "none", color = "none")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size=11),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 11, margin = margin(r = 10)),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
  )

ggsave(filename = "border10cn12_pfs_mva.svg", plot = p_mva_pfs, path = output_dir,
       width = 7, height = 5)


###
### COMPOSITION OF NEIGHBORHOODS IN BASELINE BORDER10CN12 ----
###

### Testing if myeloma infiltration and pathologist annotation are overlapping
unique(bl_border10cn12_raw$DFCI_id)

dfci_id_mrn <- clinic_master|>
  select(DFCI_number, MRN)

dfci_id_mrn$DFCI_number <- as.integer(dfci_id_mrn$DFCI_number)

bl_border10cn12_raw <- left_join(bl_border10cn12_raw, dfci_id_mrn, by=join_by("DFCI_id" == "DFCI_number"))

unique(bl_border10cn12_raw$MRN)

baseline_major_cell_composition_id <- bl_border10cn12_raw |>
  group_by(DFCI_id) |>
  count(major_cell_type) |>
  mutate(cell_type_perc = n/sum(n)) |>
  select(-n)|>
  pivot_wider(names_from = major_cell_type, values_from = cell_type_perc)

baseline_major_cell_composition$DFCI_id <- as.character(baseline_major_cell_composition$DFCI_id)

id_infiltration <- survival_bl_border10_cn12 |>
  select(DFCI_number, `PlasmaCellCore(%)`, `Plasma.Cells.(%)`, `PlasmaCellAspirate(%)`)

id_infiltration <- left_join(id_infiltration, baseline_major_cell_composition, by=join_by("DFCI_number" == "DFCI_id"))

p_pc_infiltration_cor <- id_infiltration |>
  ggplot(aes(x=`PlasmaCellCore(%)`, y=Myeloma*100))+
  geom_point(size = 3, color = "black", alpha = 0.5)+
  geom_smooth(method = "lm", color = "darkorange", alpha= 0.2) +
  labs(x="Pathologist-annotated PCs [%]", y="IMC-measured PCs [%]")+
  theme_classic()

cor.test(id_infiltration$`PlasmaCellCore(%)`, id_infiltration$Myeloma)

# Specify the output directory for subsetted plots
output_dir <- "blborder10cn12_comp_new"
dir.create(output_dir, showWarnings = FALSE)

ggsave(filename = "p_pc_infiltration_cor.svg", plot = p_pc_infiltration_cor, path = output_dir,
       width = 3, height = 3)


###
p_blborder10cn12_CNcomposition <- bl_border10cn12_raw |>
  # filter(cn_celltypes %in% c("Immune", "Myeloma", "Myeloma-MDSCs")) |>
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
  labs(x="preCN", y="Area [%]", fill = "Cell type")+
  ggtitle("preCN cell composition") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9)
  )

# Filtering significant neighborhoods associated with PFS
p_blborder10cn12_CNpfs <- bl_border10cn12_raw |>
  filter(cn_celltypes != "NaN")|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 11, 8, 9))) |>
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
  labs(x="Cell neighborhood", y="Area [%]", fill = "Cell type")+
  ggtitle("BORDER10CN12 - CN associated with survival") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9)
  )


# Filtering significant neighborhoods associated with CRS
p_blborder10cn12_CNcrs <- bl_border10cn12_raw |>
  filter(cn_celltypes %in% c(0, 1))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c(0, 1))) |>
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
  labs(x="Cell neighborhood", y="Area [%]", fill = "Cell type")+
  ggtitle("BORDER10CN12 - CN associated with CRS") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9)
  )



## Displaying cell types per neighborhood
bl_border10cn12_raw |>
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

bl_border10cn12_raw |>
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
  labs(x="Patient ID", y="BM area covered by CN [%]", fill = "CN")+
  #ggtitle("BORDER10CN12 - Neighborhood cell composition") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9)
  )

## Displaying cell types per neighborhood
bl_border10cn12_raw |>
  filter(cn_celltypes != "NaN")|>
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


##
## NEW MYELOMA ANALYSIS INCLUDING Chi2 TEST
##

# Original data processing
cn_myeloma_fractions <- bl_border10cn12_raw %>%
  filter(cn_celltypes != "NaN") %>%
  group_by(cn_celltypes) %>%
  count(functional_minor_cell_type) %>%
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma")) %>%
  mutate(perc = n / sum(n))

# Create combined data
all_myeloma_fractions <- bl_border10cn12_raw %>%
  filter(cn_celltypes != "NaN") %>%
  count(functional_minor_cell_type) %>%
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma")) %>%
  mutate(
    cn_celltypes = "All",
    perc = n / sum(n)
  )

# Combine original and combined data
myeloma_fractions <- bind_rows(cn_myeloma_fractions, all_myeloma_fractions)

# Get total counts across all cn_celltypes
total_myeloma_counts <- bl_border10cn12_raw %>%
  filter(cn_celltypes != "NaN") %>%
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma")) %>%
  count(functional_minor_cell_type) %>%
  rename(total_n = n)

# Get counts per cn_celltypes
cn_myeloma_counts <- bl_border10cn12_raw %>%
  filter(cn_celltypes != "NaN") %>%
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma")) %>%
  group_by(cn_celltypes) %>%
  count(functional_minor_cell_type) %>%
  rename(cn_n = n)

# Perform Fisher's exact test for CNS of interest
cn_list <- c("0","8","9","11")
myeloma_fisher_results <- list()

for (cn in cn_list) {
  # Get counts for the current cn_celltypes
  cn_data <- cn_myeloma_counts %>%
    filter(cn_celltypes == cn) %>%
    select(functional_minor_cell_type, cn_n)
  
  # Get counts for the rest (combined data)
  rest_data <- total_myeloma_counts %>%
    left_join(cn_data, by = "functional_minor_cell_type") %>%
    mutate(
      cn_n = replace_na(cn_n, 0),
      rest_n = total_n - cn_n
    )
  
  # Create contingency table
  contingency_table <- matrix(
    c(rest_data$cn_n, rest_data$rest_n),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Group = c(paste0("CN_", cn), "Combined"),
      Type = rest_data$functional_minor_cell_type
    )
  )
  
  # Perform Fisher's exact test
  test_result <- chisq.test(contingency_table)
  
  # Store the result
  myeloma_fisher_results[[as.character(cn)]] <- test_result
  
  # Print the result
  cat("Fisher's Exact Test for cn_celltypes =", cn, "\n")
  print(test_result)
  cat("\n")
}

# Collect p-values into a data frame
myeloma_fisher_p_values <- data.frame(
  cn_celltypes = names(myeloma_fisher_results),
  p_value = sapply(fisher_results, function(x) x$p.value)
)

# Convert 'cn_celltypes' to character for consistency
myeloma_fisher_p_values$cn_celltypes <- as.character(myeloma_fisher_p_values$cn_celltypes)

# Adjust p-values using FDR
myeloma_fisher_p_values$adjusted_p_value <- p.adjust(myeloma_fisher_p_values$p_value, method = "fdr")

myeloma_fisher_p_values$labely <- 103

# Create the plot
p_myeloma_fractions <- myeloma_fractions |>
  filter(cn_celltypes %in% c("All", "0", "8", "9", "11")) |>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "9", "11"))) |>
  ggplot(aes(x = cn_celltypes, y = perc * 100)) +
  geom_bar(aes(fill = functional_minor_cell_type), stat = 'identity') +
  scale_fill_manual(values = color_palette_2) +
  labs(x = "preCN", y = "Myeloma\nSubpopulation [%]", fill = "Functional Cell Type") +
  guides(fill = guide_legend(title = NULL)) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_text(
    data = myeloma_fisher_p_values,
    aes(x = cn_celltypes, y = labely, label = paste0("p=", signif(adjusted_p_value, 2))),
    size = 3
  ) +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0))) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.line.x.top = element_line(color = "black", linewidth = 0.5),
    axis.line.y.left = element_line(color = "black", linewidth = 0.5))


###
### Survival based on myeloma fraction
###
bl_border10cn12_func_cp_mm <- bl_border10cn12_raw |>
  select(DFCI_id, functional_minor_cell_type, area) |>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma")) |>
  group_by(DFCI_id, functional_minor_cell_type) |>
  summarise(total_area = sum(area)) |>
  group_by(DFCI_id) |>
  mutate(perc = total_area / sum(total_area)) |>
  ungroup() |>
  select(-total_area)|>
  pivot_wider(values_from = perc, names_from = functional_minor_cell_type) |>
  replace_na(list(Myeloma = 0, `Ki-67pos Myeloma` = 0, `PD-L1pos Myeloma` = 0)) 

## Subsetting myeloma based on functional markers
bl_border10cn12_func_cp_mm_variables <- colnames(bl_border10cn12_func_cp_mm[2:4])

#Isometric log ratio transform data
bl_border10cn12_func_cp_mmilr <- as.data.frame(bl_border10cn12_func_cp_mm)
bl_border10cn12_func_cp_mm_trans <- as.data.frame(ilr(bl_border10cn12_func_cp_mm))
bl_border10cn12_func_cp_mmilr[2:ncol(bl_border10cn12_func_cp_mmilr)] <- bl_border10cn12_func_cp_mm_trans
bl_border10cn12_func_cp_mmilr$DFCI_id <- as.character(bl_border10cn12_func_cp_mmilr$DFCI_id)
bl_border10cn12_func_cp_mmilr <- left_join(baseline_master, bl_border10cn12_func_cp_mmilr, by = join_by(DFCI_number == DFCI_id))
survival_bl_border10cn12_func_cp_mmilr <- bl_border10cn12_func_cp_mmilr |>
  select(all_of(survival_baseline_variables), all_of(bl_border10cn12_func_cp_mm_variables))
survival_bl_border10cn12_func_cp_mmilr <- survival_bl_border10cn12_func_cp_mmilr |>
  filter(Tissue_ID != "BS-21-J45629")

survival_bl_border10cn12_func_cp_mmilr_final <- survival_bl_border10cn12_func_cp_mmilr |>
  filter(!is.na(PFS_days) & !is.na(PFS_event) & !is.na(Myeloma))

negMM_PFS <- survfit(Surv(PFS_days, PFS_event) ~ ifelse(Myeloma > median(Myeloma), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
negMM_PFS
survdiff(Surv(PFS_days, PFS_event) ~ ifelse(Myeloma > median(Myeloma), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
p_neg_MM_PFS <- ggsurvplot(negMM_PFS,
                ylab="PFS", xlab="Days after CAR-T infusion",
                pval= TRUE, pval.coord = c(0.5, 0.2), pval.size = 3,
                size = 1.15,
                axes.offset = TRUE,
                risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
                tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
                conf.int = FALSE,
                ggtheme = theme_classic2(10),
                font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
                fontsize=3,
                legend.title = "Myeloma", 
                legend.labs = c("Low", "High"),
                palette = c("darkgrey", '#FF0000'))

negMM_OS <- survfit(Surv(OS_days, OS_event) ~ ifelse(Myeloma > median(Myeloma), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
negMM_OS
survdiff(Surv(OS_days, OS_event) ~ ifelse(Myeloma > median(Myeloma), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
p_neg_MM_OS <- ggsurvplot(negMM_OS,
           ylab="OS", xlab="Days after CAR-T infusion",
           pval= TRUE, pval.coord = c(0.5, 0.2), pval.size = 3,
           size = 1.15,
           axes.offset = TRUE,
           risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
           tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
           conf.int = FALSE,
           ggtheme = theme_classic2(10),
           font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
           fontsize=3,
           legend.title = "Myeloma", 
           legend.labs = c("Low", "High"),
           palette = c("darkgrey", '#FF0000'))

Ki67MM_PFS <- survfit(Surv(PFS_days, PFS_event) ~ ifelse(`Ki-67pos Myeloma` > median(`Ki-67pos Myeloma`), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
Ki67MM_PFS
survdiff(Surv(PFS_days, PFS_event) ~ ifelse(`Ki-67pos Myeloma` > median(`Ki-67pos Myeloma`), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
p_Ki67_MM_PFS <- ggsurvplot(Ki67MM_PFS,
           ylab="PFS", xlab="Days after CAR-T infusion",
           pval= TRUE, pval.coord = c(0.5, 0.2), pval.size = 3,
           size = 1.15,
           axes.offset = TRUE,
           risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
           tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
           conf.int = FALSE,
           ggtheme = theme_classic2(10),
           font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
           fontsize=3,
          legend.title = "MM Ki-67pos", 
           legend.labs = c("Low", "High"),
           palette = c("darkgrey", '#DC143C'))

Ki67MM_OS <- survfit(Surv(OS_days, OS_event) ~ ifelse(`Ki-67pos Myeloma` > median(`Ki-67pos Myeloma`), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
Ki67MM_OS
survdiff(Surv(OS_days, OS_event) ~ ifelse(`Ki-67pos Myeloma` > median(`Ki-67pos Myeloma`), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
p_Ki67_MM_OS <- ggsurvplot(Ki67MM_OS,
           ylab="OS", xlab="Days after CAR-T infusion",
           pval= TRUE, pval.coord = c(0.5, 0.2), pval.size = 3,
           size = 1.15,
           axes.offset = TRUE,
           risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
           tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
           conf.int = FALSE,
           ggtheme = theme_classic2(10),
           font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
           fontsize=3,
           legend.title = "MM Ki-67pos",
           legend.labs = c("Low", "High"),
           palette = c("darkgrey", '#DC143C'))


PDL1MM_PFS <- survfit(Surv(PFS_days, PFS_event) ~ ifelse(`PD-L1pos Myeloma` > median(`PD-L1pos Myeloma`), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
PDL1MM_PFS
survdiff(Surv(PFS_days, PFS_event) ~ ifelse(`PD-L1pos Myeloma` > median(`PD-L1pos Myeloma`), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
p_PDL1_MM_PFS <- ggsurvplot(PDL1MM_PFS,
           ylab="PFS", xlab="Days after CAR-T infusion",
           pval= TRUE, pval.coord = c(0.5, 0.2), pval.size = 3,
           size = 1.15,
           axes.offset = TRUE,
           risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
           tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
           conf.int = FALSE,
           ggtheme = theme_classic2(10),
           font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
           fontsize=3,
           legend.title = "MM PD-L1pos",
           legend.labs = c("Low", "High"),
           palette = c("darkgrey", '#B22222'))

PDL1MM_OS <- survfit(Surv(OS_days, OS_event) ~ ifelse(`PD-L1pos Myeloma` > median(`PD-L1pos Myeloma`), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
PDL1MM_OS
survdiff(Surv(OS_days, OS_event) ~ ifelse(`PD-L1pos Myeloma` > median(`PD-L1pos Myeloma`), 1, 0), data = survival_bl_border10cn12_func_cp_mmilr_final)
p_PDL1_MM_OS <-ggsurvplot(PDL1MM_OS,
           ylab="OS", xlab="Days after CAR-T infusion",
           pval= TRUE, pval.coord = c(0.5, 0.2), pval.size = 3,
           size = 1.15,
           axes.offset = TRUE,
           risk.table=FALSE, risk.table.title="No. at risk", risk.table.height=.19,
           tables.y.text = FALSE, tables.theme = theme_cleantable(base_size = 2),
           conf.int = FALSE,
           ggtheme = theme_classic2(10),
           font.title=c(9, "bold"), font.tickslab = c(9), font.legend.labs=c(9), font.x = c(9, "bold"), font.y = c(9, "bold"),
           fontsize=3,
           legend.title = "MM PD-L1pos",
           legend.labs = c("Low", "High"),
           palette = c("darkgrey", '#B22222'))


## Saving composition  plots

# Retrieve the plots using the pattern
pattern_survival_MM <- "p_[a-zA-Z0-9]+_MM_[a-zA-Z]"
plots_survival_MM <- mget(ls(pattern = pattern_survival_MM))

# Define the directory to save the plots
output_dir <- "blborder10cn12_survival"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(plots_survival_MM)) {
  
  # Extract the ggplot component from the ggsurvplot object
  MM_survival_plot <- plots_survival_MM[[name]]$plot
  
  # Define the file name and path
  subset_file_name <- paste0("border10cn12_", name, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = MM_survival_plot, width = 3, height = 2)
}


###
### NEW COMPOSITION DATA for SURVIVAL ----
###

# MYELOMA
myeloma_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

myeloma_fractions_cn <- bl_border10cn12_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(functional_minor_cell_type) %>%
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n / sum(n)) %>%
  ungroup()

myeloma_fractions <- rbind(myeloma_fractions_all, myeloma_fractions_cn)

p_new_ki67mm_all <- myeloma_fractions|>
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

p_new_ki67mm_cnoi <- myeloma_fractions|>
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

myeloma_fractions|>
  filter(functional_minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))


p_new_mmneg_all <- myeloma_fractions|>
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

p_new_mmneg_cnoi <- myeloma_fractions|>
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

myeloma_fractions|>
  filter(functional_minor_cell_type == "Myeloma")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

p_new_pdl1mm_all <- myeloma_fractions|>
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

p_new_pdl1mm_cnoi <- myeloma_fractions|>
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

myeloma_fractions|>
  filter(functional_minor_cell_type == "PD-L1pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))
  

## NEW COMPOSITION DATA
# CD4
cd4_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD4 Tnaive","Non Th17 CD4 Tmem", "Th17 CD4 Tmem", "CD4 Tex", "Treg"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

cd4_fractions_cn <- bl_border10cn12_raw |>
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

#CD4 Tnaive
p_new_cd4naive_all <- cd4_fractions|>
  filter(minor_cell_type == "CD4 Tnaive")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "CD4 Tnaive")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd4naive_cnoi <- cd4_fractions|>
  filter(minor_cell_type == "CD4 Tnaive")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "CD4 Tnaive")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

cd4_fractions|>
  filter(minor_cell_type == "CD4 Tnaive")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

#CD Tex
p_new_cd4ex_all <- cd4_fractions|>
  filter(minor_cell_type == "CD4 Tex")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "CD4 Tex")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd4ex_cnoi <- cd4_fractions|>
  filter(minor_cell_type == "CD4 Tex")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "CD4 Tex")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

cd4_fractions|>
  filter(minor_cell_type == "CD4 Tex")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

# Non TH17 CD4 Tmem
p_new_cd4tmem <-cd4_fractions|>
  filter(minor_cell_type == "Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Non TH17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.05,
           step.increase = 0.08)+
  theme_classic()

p_new_cd4tmem_all <- cd4_fractions|>
  filter(minor_cell_type == "Non Th17 CD4 Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Non Th17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd4tmem_cnoi <- cd4_fractions|>
  filter(minor_cell_type == "Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "Non Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# Th17 CD4 Tmem
p_new_cd4th17 <-cd4_fractions|>
  filter(minor_cell_type == "Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "TH17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.5,
           step.increase = 0.05)+
  theme_classic()

p_new_cd4th17_all <- cd4_fractions|>
  filter(minor_cell_type == "Th17 CD4 Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Th17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd4th17_cnoi <- cd4_fractions|>
  filter(minor_cell_type == "Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

cd4_fractions|>
  filter(minor_cell_type == "Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

# Treg
p_new_cd4treg <-cd4_fractions|>
  filter(minor_cell_type == "Treg")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Treg")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.25,
           step.increase = 0.05)+
  theme_classic()

p_new_cd4treg_all <- cd4_fractions|>
  filter(minor_cell_type == "Treg")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Treg")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd4treg_cnoi <- cd4_fractions|>
  filter(minor_cell_type == "Treg")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,100), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "Treg")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


## NEW COMPOSITION DATA
# Effector CD4 FUNCTIONAL LEVEL including pSTAT1
cd4f_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos Non Th17 CD4 Tmem", "pSTAT1pos Th17 CD4 Tmem", "pSTAT1neg Non Th17 CD4 Tmem", "pSTAT1neg Th17 CD4 Tmem"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

cd4f_fractions_cn <- bl_border10cn12_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos Non Th17 CD4 Tmem", "pSTAT1pos Th17 CD4 Tmem", "pSTAT1neg Non Th17 CD4 Tmem", "pSTAT1neg Th17 CD4 Tmem"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

cd4f_fractions <- rbind(cd4f_fractions_all, cd4f_fractions_cn)

# pSTAT1pos Non TH17 CD4 Tmem
p_new_cd4tmempstat1pos_all <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos Non Th17 CD4 Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Non Th17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd4tmempstat1pos_cnoi <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Non Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1pos Th17 CD4 Tmem
p_new_cd4th17pstat1pos_all <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos Th17 CD4 Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Th17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd4th17pstat1pos_cnoi <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1pos Non TH17 CD4 Tmem
p_new_cd4tmempstat1neg_all <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg Non Th17 CD4 Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Non Th17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd4tmempstat1neg_cnoi <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Non Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1neg Th17 CD4 Tmem
p_new_cd4th17pstat1neg_all <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg Th17 CD4 Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Th17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd4th17pstat1neg_cnoi <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,40), expand = c(0,0), breaks = c(0,10,20,30,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


## NEW COMPOSITION DATA - CD8
# CD8
cd8_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD8 Tnaive","CD8 GZMB+ Tmem", "CD8 Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

cd8_fractions_cn <- bl_border10cn12_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD8 Tnaive","CD8 GZMB+ Tmem", "CD8 Tmem", "CD8 Tex"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

cd8_fractions <- rbind(cd8_fractions_all, cd8_fractions_cn)

cd8_fractions|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulations [%]")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.05)+
  theme_classic()+
  facet_wrap(vars(minor_cell_type))

#CD8 Tnaive
p_new_cd8naive <-cd8_fractions|>
  filter(minor_cell_type == "CD8 Tnaive")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tnaive")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.05,
           step.increase = 0.1)+
  theme_classic()

p_new_cd8naive_all <- cd8_fractions|>
  filter(minor_cell_type == "CD8 Tnaive")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tnaive")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd8naive_cnoi <- cd8_fractions|>
  filter(minor_cell_type == "CD8 Tnaive")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 Tnaive")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#CD8 Tex
p_new_cd8ex <-cd8_fractions|>
  filter(minor_cell_type == "CD8 Tex")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,50), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tex")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.6,
           step.increase = 0.05)+
  theme_classic()

p_new_cd8ex_all <- cd8_fractions|>
  filter(minor_cell_type == "CD8 Tex")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tex")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd8ex_cnoi <- cd8_fractions|>
  filter(minor_cell_type == "CD8 Tex")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 Tex")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

cd8_fractions|>
 filter(minor_cell_type == "CD8 Tex")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

#CD8 Tmem
p_new_cd8tmem <-cd8_fractions|>
  filter(minor_cell_type == "CD8 Tmem")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tmem")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = -0.08)+
  theme_classic()

p_new_cd8tmem_all <- cd8_fractions|>
  filter(minor_cell_type == "CD8 Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd8tmem_cnoi <- cd8_fractions|>
  filter(minor_cell_type == "CD8 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#CD8 GZMB+ Tmem
p_new_cd8gzmb <-cd8_fractions|>
  filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 GZMB+ Tmem")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.08)+
  theme_classic()

p_new_cd8gzmb_all <- cd8_fractions|>
  filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 GZMB+ Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd8gzmb_cnoi <- cd8_fractions|>
  filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 GZMB+ Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

cd8_fractions|>
  filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))
  
## NEW COMPOSITION DATA
# Effector CD4 FUNCTIONAL LEVEL including pSTAT1
cd8f_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos CD8 GZMBpos Tmem", "pSTAT1pos CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem", "pSTAT1neg CD8 Tmem"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

cd8f_fractions_cn <- bl_border10cn12_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos CD8 GZMBpos Tmem", "pSTAT1pos CD8 Tmem", "pSTAT1neg CD8 GZMBpos Tmem", "pSTAT1neg CD8 Tmem"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

cd8f_fractions <- rbind(cd8f_fractions_all, cd8f_fractions_cn)

# pSTAT1pos CD8 GZMBpos Tmem
p_new_cd8tmemgzmbpstat1pos_all <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD8 GZMBpos Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos CD8 GZMBpos Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd8tmemgzmbpstat1pos_cnoi <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD8 GZMBpos Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos CD8 GZMBpos Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1pos CD8 Tmem
p_new_cd8tmempstat1pos_all <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD8 Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos GZMBpos CD8 Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd8tmempstat1pos_cnoi <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD8 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos GZMBpos CD8 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1pos CD8 GZMBpos Tmem
p_new_cd8tmemgzmbpstat1neg_all <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD8 GZMBpos Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 GZMBpos Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd8tmemgzmbpstat1neg_cnoi <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD8 GZMBpos Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 GZMBpos Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1neg CD8 Tmem
p_new_cd8tmempstat1neg_all <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD8 Tmem")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 Tmem")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_cd8tmempstat1neg_cnoi <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD8 Tmem")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


## NEW COMPOSITION DATA - MYELOID
# MYELOID
myeloid_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("Classical monocytes", "Intermediate monocytes", "MDSCs", "Nonclassical monocytes",
                                "M1-like macrophage", "M2-like macrophage", "MoDC", "MKs"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

myeloid_fractions_cn <- bl_border10cn12_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("Classical monocytes", "Intermediate monocytes", "MDSCs", "Nonclassical monocytes",
                                "M1-like macrophage", "M2-like macrophage", "MoDC", "MKs"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

myeloid_fractions <- rbind(myeloid_fractions_all, myeloid_fractions_cn)

myeloid_fractions|>
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

#Classical moncytes
p_new_classmono <-myeloid_fractions|>
  filter(minor_cell_type == "Classical monocytes")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Classical monocytes")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.2,,
           step.increase = 0.08)+
  theme_classic()

p_new_classmono_all <- myeloid_fractions|>
  filter(minor_cell_type == "Classical monocytes")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Classical monocytes")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_classmono_cnoi <- myeloid_fractions|>
  filter(minor_cell_type == "Classical monocytes")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Classical monocytes")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#Intermediate monocytes
p_new_itmmono <-myeloid_fractions|>
  filter(minor_cell_type == "Intermediate monocytes")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,40), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Intermediate monocytes")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.3,
           step.increase = 0.05)+
  theme_classic()

p_new_itmmono_all <- myeloid_fractions|>
  filter(minor_cell_type == "Intermediate monocytes")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Intermediate monocytes")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_itmmono_cnoi <- myeloid_fractions|>
  filter(minor_cell_type == "Intermediate monocytes")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Intermediate monocytes")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#Nonclassical monocytes
p_new_nonclassmono <-myeloid_fractions|>
  filter(minor_cell_type == "Nonclassical monocytes")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Nonclassical monocytes")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.08)+
  theme_classic()

p_new_nonclassmono_all <- myeloid_fractions|>
  filter(minor_cell_type == "Nonclassical monocytes")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Nonclassical monocytes")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_nonclassmono_cnoi <- myeloid_fractions|>
  filter(minor_cell_type == "Nonclassical monocytes")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Nonclassical monocytes")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#MDSCs
p_new_mdsc <-myeloid_fractions|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "MDSCs")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.08)+
  theme_classic()

p_new_mdsc_all <- myeloid_fractions|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "MDSCs")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_mdsc_cnoi <- myeloid_fractions|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "MDSCs")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

myeloid_fractions|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

#M1-like macrophage
p_new_m1mac <-myeloid_fractions|>
  filter(minor_cell_type == "M1-like macrophage")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,5), expand = c(0,0), breaks = c(0,1,2,3,4,5,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M1-like macrophage")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.3)+
  theme_classic()

p_new_m1mac_all <- myeloid_fractions|>
  filter(minor_cell_type == "M1-like macrophage")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M1-like macrophage")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_m1mac_cnoi <- myeloid_fractions|>
  filter(minor_cell_type == "M1-like macrophage")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,10), expand = c(0,0), breaks = c(0,2,4,6,8,10))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "M1-like macrophage")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()



#M2-like macrophage
p_new_m2mac <-myeloid_fractions|>
  filter(minor_cell_type == "M2-like macrophage")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M2-like macrophage")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
  theme_classic()

p_new_m2mac_all <- myeloid_fractions|>
  filter(minor_cell_type == "M2-like macrophage")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M2-like macrophage")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_m2mac_cnoi <- myeloid_fractions|>
  filter(minor_cell_type == "M2-like macrophage")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,100), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "M2-like macrophage")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

myeloid_fractions|>
  filter(minor_cell_type == "M2-like macrophage")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))
  
#MoDC
p_new_modc <-myeloid_fractions|>
  filter(minor_cell_type == "MoDC")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,20), expand = c(0,0), breaks = c(0,10,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "MoDC")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.1)+
  theme_classic()

p_new_modc_all <- myeloid_fractions|>
  filter(minor_cell_type == "MoDC")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "MoDC")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_modc_cnoi <- myeloid_fractions|>
  filter(minor_cell_type == "MoDC")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,30), expand = c(0,0), breaks = c(0,5,10,15,20,25))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "MoDC")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#MKs
p_new_mk_all <- myeloid_fractions|>
  filter(minor_cell_type == "MKs")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Megakaryocytes")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_mk_cnoi <- myeloid_fractions|>
  filter(minor_cell_type == "MKs")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Megakaryocytes")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

## NEW COMPOSITION DATA - MYELOID
# MYELOID
stroma_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("ECs", "Adipocyte"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

stroma_fractions_cn <- bl_border10cn12_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("ECs", "Adipocyte"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

stroma_fractions <- rbind(stroma_fractions_all, stroma_fractions_cn)

stroma_fractions|>
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

#Endothelial cells
p_new_ec_all <- stroma_fractions|>
  filter(minor_cell_type == "ECs")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Stroma cell subpopulation [%]", title = "ECs")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_ec_cnoi <- stroma_fractions|>
  filter(minor_cell_type == "ECs")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Stroma cell subpopulation [%]", title = "ECs")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#Endothelial cells
p_new_adipocyte_all <- stroma_fractions|>
  filter(minor_cell_type == "Adipocyte")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Stroma cell subpopulation [%]", title = "Adipocyte")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_adipocyte_cnoi <- stroma_fractions|>
  filter(minor_cell_type == "Adipocyte")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "Stroma cell subpopulation [%]", title = "Adipocyte")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

## NEW COMPOSITION DATA - NK CELLS
# NK
nk_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD56 dim NK", "CD56 bright NK"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

nk_fractions_cn <- bl_border10cn12_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("CD56 dim NK", "CD56 bright NK"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

nk_fractions <- rbind(nk_fractions_all, nk_fractions_cn)

nk_fractions|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
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

#Dim NK cells
p_new_dimnk_all <- nk_fractions|>
  filter(minor_cell_type == "CD56 dim NK")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "CD56 dim NK")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_dimnk_cnoi <- nk_fractions|>
  filter(minor_cell_type == "CD56 dim NK")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "CD56 dim NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#Bright NK cells
p_new_brightnk_all <- nk_fractions|>
  filter(minor_cell_type == "CD56 bright NK")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "CD56 bright NK")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_brightnk_cnoi <- nk_fractions|>
  filter(minor_cell_type == "CD56 bright NK")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "CD56 bright NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


## NEW COMPOSITION DATA - FUNCTIONAL NK CELLS
# FUNCTIONAL NK
nkf_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos CD56 dim NK", "pSTAT1pos CD56 bright NK", "pSTAT1neg CD56 dim NK", "pSTAT1neg CD56 bright NK"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

nkf_fractions_cn <- bl_border10cn12_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(functional_minor_cell_type)|>
  filter(functional_minor_cell_type %in% c("pSTAT1pos CD56 dim NK", "pSTAT1pos CD56 bright NK", "pSTAT1neg CD56 dim NK", "pSTAT1neg CD56 bright NK"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

nkf_fractions <- rbind(nkf_fractions_all, nkf_fractions_cn)

nkf_fractions|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
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

#pSTAT1pos Dim NK cells
p_new_pstat1posdimnk_all <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD56 dim NK")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 dim NK")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_pstat1posdimnk_cnoi <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD56 dim NK")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 dim NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#pSTAT1pos Bright NK cells
p_new_pstat1posbrightnk_all <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD56 bright NK")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 bright NK")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_pstat1posbrightnk_cnoi <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD56 bright NK")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 bright NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#pSTAT1neg Dim NK cells
p_new_pstat1negdimnk_all <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD56 dim NK")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 dim NK")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_pstat1negdimnk_cnoi <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD56 dim NK")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 dim NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#pSTAT1neg Bright NK cells
p_new_pstat1negbrightnk_all <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD56 bright NK")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 bright NK")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_pstat1negbrightnk_cnoi <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD56 bright NK")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 bright NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

## NEW COMPOSITION DATA - B CELLS
# B
b_fractions_all <- bl_border10cn12_raw |>
  group_by(DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("B cells", "MBCs"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

b_fractions_cn <- bl_border10cn12_raw |>
  group_by(cn_celltypes, DFCI_id)|>
  count(minor_cell_type)|>
  filter(minor_cell_type %in% c("B cells", "MBCs"))|>
  mutate(perc = n/sum(n))|>
  ungroup()

b_fractions <- rbind(b_fractions_all, b_fractions_cn)

b_fractions|>
  #filter(minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "8","9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "8","9")))|>
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

#Dim NK cells
p_new_b_all <- b_fractions|>
  filter(minor_cell_type == "B cells")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "B cell subpopulation [%]", title = "B cells")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_b_cnoi <- b_fractions|>
  filter(minor_cell_type == "B cells")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "B cell subpopulation [%]", title = "B cells")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#Bright NK cells
p_new_mbcs_all <- b_fractions|>
  filter(minor_cell_type == "MBCs")|>
  filter(cn_celltypes != "NaN")|>
  ggplot(aes(x=reorder(cn_celltypes, -perc, FUN = median), y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  #scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "B cell subpopulation [%]", title = "MBCs")+
  guides(fill = "none", color = "none")+
  theme_classic()

p_new_mbcs_cnoi <- b_fractions|>
  filter(minor_cell_type == "MBCs")|>
  filter(cn_celltypes %in% c("All", "0", "8", "11", "9"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "8", "11", "9")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN8", "preCN11", "preCN9"))+
  labs(x = "", y = "B cell subpopulation [%]", title = "MBCs")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# Retrieve the plots using the pattern
pattern_comp_bl_new_all <- "p_new_[a-zA-Z0-9]+_all"
plots_comp_bl_new_all <- mget(ls(pattern = pattern_comp_bl_new_all))
pattern_comp_bl_new_cnoi <- "p_new_[a-zA-Z0-9]+_cnoi"
plots_comp_bl_new_cnoi<- mget(ls(pattern = pattern_comp_bl_new_cnoi))

# Define the directory to save the plots
output_dir <- "blborder10cn12_comp_new"
dir.create(output_dir, showWarnings = FALSE)


# Save the plots to the specified directory
save_plot(plots_comp_bl_new_all, output_dir, width = 5, height = 3)
save_plot(plots_comp_bl_new_cnoi, output_dir, width = 3.3, height = 4)

##
## Survival based on myeloma
##

## Preparing cell proportions table
bl_border10cn12_cc_func <- bl_border10cn12_raw |>
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
### Interaction analyses ----
###

bl_border10cn12_func_int <- read.csv("full_results/border/baseline/border-10_CNs-12_interactions__functional_baseline.csv")

bl_border10cn12_func_int <- bl_border10cn12_func_int |> select(-X)
colnames(bl_border10cn12_func_int)<- gsub("\\.", " ", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int)<- gsub("Ki 67pos Myeloma", "Ki-67pos\nMyeloma", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int)<- gsub("PD L1pos Myeloma", "PD-L1pos\nMyeloma", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int)<- gsub("M1 like M  phi ", "M1-like\nmacrophage", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int)<- gsub("M2 like M  phi ", "M2-like\nmacrophage", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1pos CD8 Tmem", "pSTAT1pos\nCD8 Tmem", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1neg CD8 Tmem", "pSTAT1neg\nCD8 Tmem", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1pos CD8 GZMBpos Tmem", "pSTAT1pos\nCD8 GZMBpos Tmem", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1neg CD8 GZMBpos Tmem", "pSTAT1neg\nCD8 GZMBpos Tmem", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1neg Non Th17 CD4 Tmem", "pSTAT1neg\nNon Th17 CD4 Tmem", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1pos Non Th17 CD4 Tmem", "pSTAT1pos\nNon Th17 CD4 Tmem", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1neg Th17 CD4 Tmem", "pSTAT1neg\nTh17 CD4 Tmem", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1pos Th17 CD4 Tmem", "pSTAT1pos\nTh17 CD4 Tmem", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1neg CD56 dim NK", "pSTAT1neg\nCD56 dim NK", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1neg CD56 bright NK", "pSTAT1neg\nCD56 bright NK", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1pos CD56 dim NK", "pSTAT1pos\nCD56 dim NK", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("pSTAT1pos CD56 bright NK", "pSTAT1pos\nCD56 bright NK", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("Classical monocytes", "Classical\nmonocytes", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub("Nonclassical monocytes", "Nonclassical\nmonocytes", colnames(bl_border10cn12_func_int))
colnames(bl_border10cn12_func_int) <- gsub( "Intermediate monocytes",  "Intermediate\nmonocytes", colnames(bl_border10cn12_func_int))

border10cn12_cn_values <- unique(bl_border10cn12_func_int$cn_celltypes)

for (cn in border10cn12_cn_values) {
  # Create a new data frame for each CN value
  new_df <- bl_border10cn12_func_int |>
    filter(cn_celltypes == cn) |>
    select(-cn_celltypes) |>
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
  assign(paste0("bl_border10cn12_func_int_CN", cn), new_df)
}

for (cn in border10cn12_cn_values) {
  # Create a new data frame for each CN value
  new_df <- bl_border10cn12_func_int |>
    filter(cn_celltypes == cn) |>
    select(-cn_celltypes) |>
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
  
  # Fill NA symmetrically
  for (i in seq_len(nrow(new_df))) {
    for (j in seq_len(ncol(new_df))) {
      if (is.na(new_df[i, j]) && !is.na(new_df[j, i])) {
        new_df[i, j] <- new_df[j, i]
      } else if (!is.na(new_df[i, j]) && is.na(new_df[j, i])) {
        new_df[j, i] <- new_df[i, j]
      }
    }
  }
  
  # Dynamically create a new variable name for the data frame
  assign(paste0("bl_border10cn12_func_int_CN", cn), new_df)
}

# Create a list to hold all the data frames
cn_list_bor <- list()
plots_list_bor <- list()

# Assuming the data frames are named bl_border10cn12_func_int_CN0 to bl_border10cn12_func_int_CN19
for (i in 0:11) {
  # Dynamically assign the data frames to the list
  cn_list_bor[[i+1]] <- get(paste0("bl_border10cn12_func_int_CN", i))
}

# Specify the output directory
output_dir <- "border10cn12_chord_plots"
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

chordDiagram(as.matrix(bl_border10cn12_func_int_CN0), scale = TRUE, grid.col = color_palette)
circos.clear()


### Subsetting cell-cell contacts based on specific questions

# Test example
# bl_border10cn12_MM_CD4_int_CN0 <- bl_border10cn12_func_int_CN0 |>
#   select(Myeloma, 'Ki-67pos\nMyeloma', 'PD-L1pos\nMyeloma') |>
#   filter(rownames(bl_border10cn12_func_int_CN0) %in% c("CD4 Tnaive", "CD4 Tex", "pSTAT1neg\nCD4 Tmem", "pSTAT1pos\nCD4 Tmem", "Treg"))
# 
# chordDiagram(as.matrix(MM_CD4_int_CN0), scale = T, grid.col = color_palette)
# circos.clear()


## Subsetting contacts between MM and CD4 T cells
# Create a list to hold the subsetted data frames and plots
MM_CD4_int_cn_list_bor <- list()
plots_list_MM_CD4_int_bor <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_MM_CD4_int_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# List of rows to subset for Myeloma - CD4 T cell contacts
MM_CD4_int_row_names <- c("CD4 Tnaive", "CD4 Tex", "pSTAT1neg\nTh17 CD4 Tmem", "pSTAT1pos\nTh17 CD4 Tmem",
                          "pSTAT1neg\nNon Th17 CD4 Tmem", "pSTAT1pos\nNon Th17 CD4 Tmem", "Treg")

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_list_bor)) {
  cn_df <- cn_list_bor[[i]]
  
  ### Subset the data frame (for Myeloma - CD4 T cell contacts)
  subset_cn_df <- cn_df |>
    select(Myeloma, 'Ki-67pos\nMyeloma', 'PD-L1pos\nMyeloma') |>
    filter(rownames(cn_df) %in% MM_CD4_int_row_names)
  
  # Add the subsetted data frame to the subset list
  MM_CD4_int_cn_list_bor[[i]] <- subset_cn_df
  
  # Save the subset plot as a PNG file
  subset_file_name <- paste0("MM_CD4_int_chord_plot_", i, ".png")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving MM_CD4_int plot to:", subset_file_path))  # Debugging line
  
  # Create the PNG for the subset
  png(subset_file_path, width = 800, height = 800)
  chordDiagram(as.matrix(subset_cn_df), scale = TRUE, grid.col = color_palette)
  circos.clear()
  dev.off()
  
  # Load the saved image back as a grob using the correct path
  subset_img <- rasterGrob(png::readPNG(subset_file_path), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "MM-CD4"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_MM_CD4_int_bor[[i]] <- plot_with_subset_title
}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_MM_CD4_int_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_MM_CD4_int_bor, ncol = 5)
dev.off()



## Subsetting contacts between MM and CD8 T cells
# Create a list to hold the subsetted data frames and plots
MM_CD8_int_cn_list_bor <- list()
plots_list_MM_CD8_int_bor <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_MM_CD8_int_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# List of rows to subset for Myeloma - CD8 T cell contacts
MM_CD8_int_row_names <- c("CD8 Tnaive", "pSTAT1pos\nCD8 Tmem", "pSTAT1neg\nCD8 Tmem", "pSTAT1pos\nCD8 GZMBpos Tmem", "pSTAT1neg\nCD8 GZMBpos Tmem", "CD8 Tex")

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_list_bor)) {
  cn_df <- cn_list_bor[[i]]
  
  ### Subset the data frame (for Myeloma - CD8 T cell contacts)
  subset_cn_df <- cn_df |>
    select(Myeloma, 'Ki-67pos\nMyeloma', 'PD-L1pos\nMyeloma') |>
    filter(rownames(cn_df) %in% MM_CD8_int_row_names)
  
  # Add the subsetted data frame to the subset list
  MM_CD8_int_cn_list_bor[[i]] <- subset_cn_df
  
  # Save the subset plot as a PNG file
  subset_file_name <- paste0("MM_CD8_int_chord_plot_", i, ".png")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving MM_CD8_int plot to:", subset_file_path))  # Debugging line
  
  # Create the PNG for the subset
  png(subset_file_path, width = 800, height = 800)
  chordDiagram(as.matrix(subset_cn_df), scale = TRUE, grid.col = color_palette)
  circos.clear()
  dev.off()
  
  # Load the saved image back as a grob using the correct path
  subset_img <- rasterGrob(png::readPNG(subset_file_path), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "MM-CD8"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_MM_CD8_int_bor[[i]] <- plot_with_subset_title
}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_MM_CD8_int_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_MM_CD8_int_bor, ncol = 5)
dev.off()


## Subsetting contacts between MM and immunosuppressive cells
# Create a list to hold the subsetted data frames and plots
MM_SUP_int_cn_list_bor <- list()
plots_list_MM_SUP_int_bor <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_MM_SUP_int_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# List of rows to subset for Myeloma - SUP T cell contacts
MM_SUP_int_row_names <- c("M2-like\nmacrophage", "M1-like\nmacrophage", "MDSCs")

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_list_bor)) {
  cn_df <- cn_list_bor[[i]]
  
  ### Subset the data frame (for Myeloma - SUP T cell contacts)
  subset_cn_df <- cn_df |>
    select(Myeloma, 'Ki-67pos\nMyeloma', 'PD-L1pos\nMyeloma') |>
    filter(rownames(cn_df) %in% MM_SUP_int_row_names)
  
  # Add the subsetted data frame to the subset list
  MM_SUP_int_cn_list_bor[[i]] <- subset_cn_df
  
  # Save the subset plot as a PNG file
  subset_file_name <- paste0("MM_SUP_int_chord_plot_", i, ".png")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving MM_SUP_int plot to:", subset_file_path))  # Debugging line
  
  # Create the PNG for the subset
  png(subset_file_path, width = 800, height = 800)
  chordDiagram(as.matrix(subset_cn_df), scale = TRUE, grid.col = color_palette)
  circos.clear()
  dev.off()
  
  # Load the saved image back as a grob using the correct path
  subset_img <- rasterGrob(png::readPNG(subset_file_path), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "MM-SUP"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_MM_SUP_int_bor[[i]] <- plot_with_subset_title
}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_MM_SUP_int_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_MM_SUP_int_bor, ncol = 5)
dev.off()



## Subsetting contacts between MM and immunosuppressive cells
# Create a list to hold the subsetted data frames and plots
MM_CD4_SUP_int_cn_list_bor <- list()
plots_list_MM_CD4_SUP_int_bor <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_MM_CD4_SUP_int_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# List of rows to subset for Myeloma - SUP T cell contacts
MM_CD4_SUP_int_row_names <- c("CD4 Tnaive", "CD4 Tex", "pSTAT1neg\nTh17 CD4 Tmem", "pSTAT1pos\nTh17 CD4 Tmem",
                              "pSTAT1neg\nNon Th17 CD4 Tmem", "pSTAT1pos\nNon Th17 CD4 Tmem",  
                              "Treg", "M2-like\nmacrophage", "M1-like\nmacrophage", "MDSCs")

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_list_bor)) {
  cn_df <- cn_list_bor[[i]]
  
  ### Subset the data frame (for Myeloma - SUP T cell contacts)
  subset_cn_df <- cn_df |>
    select(Myeloma, 'Ki-67pos\nMyeloma', 'PD-L1pos\nMyeloma', "CD4 Tnaive", "CD4 Tex", "pSTAT1neg\nTh17 CD4 Tmem", "pSTAT1pos\nTh17 CD4 Tmem",
           "pSTAT1neg\nNon Th17 CD4 Tmem", "pSTAT1pos\nNon Th17 CD4 Tmem",  "Treg") |>
    filter(rownames(cn_df) %in% MM_CD4_SUP_int_row_names)
  
  # Remove links
  subset_cn_df[(4:10),(4:10)] <- NA
  
  # Add the subsetted data frame to the subset list
  MM_CD4_SUP_int_cn_list_bor[[i]] <- subset_cn_df
  
  # Save the subset plot as a PNG file
  subset_file_name_svg <- paste0("MM_CD4_SUP_int_chord_plot_", i, ".svg")
  subset_file_path_svg <- file.path(output_dir, subset_file_name_svg)
  print(paste("Saving MM_CD4_SUP_int plot to:", subset_file_path_svg))  # Debugging line
  
  # Create the SVG for the subset with a pre-allocated track for customizing text
  svg(subset_file_path_svg, width = 9, height = 9)
  chordDiagram(as.matrix(subset_cn_df), scale = TRUE, grid.col = color_palette, 
               annotationTrack = c("grid"), preAllocateTracks = 1)
  
  # Customize sector labels (adjust font size)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name, 
                facing = "clockwise", niceFacing = T, adj = c(0, 0.7), cex = 1.0)  # Adjust 'cex' for font size
  }, bg.border = NA)
  
  # # Add the numbers (axis ticks) to the sectors
  # circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  #   circos.axis(h = "top", major.tick.length = 0.2, labels.cex = 0.6, direction = "outside")
  # }, bg.border = NA)
  
  circos.clear()
  dev.off()
  
  # Convert the SVG to a PNG using rsvg
  subset_file_name_png <- paste0("MM_CD4_SUP_int_chord_plot_", i, ".png")
  subset_file_path_png <- file.path(output_dir, subset_file_name_png)
  
  rsvg_png(subset_file_path_svg, subset_file_path_png)  # Convert SVG to PNG
  
  # Load the saved PNG back as a grob
  subset_img <- rasterGrob(png::readPNG(subset_file_path_png), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "MM-SUP"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_MM_CD4_SUP_int_bor[[i]] <- plot_with_subset_title

}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_MM_CD4_SUP_int_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_MM_CD4_SUP_int_bor, ncol = 5)
dev.off()


## Subsetting contacts between MM and immunosuppressive cells
# Create a list to hold the subsetted data frames and plots
MM_CD8_SUP_int_cn_list_bor <- list()
plots_list_MM_CD8_SUP_int_bor <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_MM_CD8_SUP_int_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# List of rows to subset for Myeloma - SUP T cell contacts
MM_CD8_SUP_int_row_names <- c("CD8 Tnaive", "pSTAT1pos\nCD8 Tmem", "pSTAT1neg\nCD8 Tmem", "pSTAT1pos\nCD8 GZMBpos Tmem", "pSTAT1neg\nCD8 GZMBpos Tmem", "CD8 Tex", "M2-like\nmacrophage", "M1-like\nmacrophage", "MDSCs")

# Iterate over the list of full data frames and apply the subsetting

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_list_bor)) {
  cn_df <- cn_list_bor[[i]]
  
  ### Subset the data frame (for Myeloma - SUP T cell contacts)
  subset_cn_df <- cn_df |>
    select(Myeloma, 'Ki-67pos\nMyeloma', 'PD-L1pos\nMyeloma', "CD8 Tnaive", "pSTAT1pos\nCD8 Tmem", "pSTAT1neg\nCD8 Tmem", "pSTAT1pos\nCD8 GZMBpos Tmem", "pSTAT1neg\nCD8 GZMBpos Tmem", "CD8 Tex") |>
    filter(rownames(cn_df) %in% MM_CD8_SUP_int_row_names)
  
  # Remove links
  subset_cn_df[(4:8), (4:8)] <- NA
  
  # Add the subsetted data frame to the subset list
  MM_CD8_SUP_int_cn_list_bor[[i]] <- subset_cn_df
  
  # Save the subset plot as an SVG file
  subset_file_name_svg <- paste0("MM_CD8_SUP_int_chord_plot_", i, ".svg")
  subset_file_path_svg <- file.path(output_dir, subset_file_name_svg)
  print(paste("Saving MM_CD8_SUP_int plot to:", subset_file_path_svg))  # Debugging line
  
  # Create the SVG for the subset with a pre-allocated track for customizing text
  svg(subset_file_path_svg, width = 9, height = 9)
  chordDiagram(as.matrix(subset_cn_df), scale = TRUE, grid.col = color_palette, 
               annotationTrack = c("grid"), preAllocateTracks = 1)
  
  # Customize sector labels (adjust font size)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name, 
                facing = "clockwise", niceFacing = T, adj = c(0, 0.7), cex = 1.0)  # Adjust 'cex' for font size
  }, bg.border = NA)
  
  # # Add the numbers (axis ticks) to the sectors
  # circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  #   circos.axis(h = "top", major.tick.length = 0.2, labels.cex = 0.6, direction = "outside")
  # }, bg.border = NA)
  
  circos.clear()
  dev.off()
  
  # Convert the SVG to a PNG using rsvg
  subset_file_name_png <- paste0("MM_CD8_SUP_int_chord_plot_", i, ".png")
  subset_file_path_png <- file.path(output_dir, subset_file_name_png)
  
  rsvg_png(subset_file_path_svg, subset_file_path_png)  # Convert SVG to PNG
  
  # Load the saved PNG back as a grob
  subset_img <- rasterGrob(png::readPNG(subset_file_path_png), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "MM-SUP"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_MM_CD8_SUP_int_bor[[i]] <- plot_with_subset_title
}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_MM_CD8_SUP_int_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_MM_CD8_SUP_int_bor, ncol = 5)
dev.off()




## Subsetting contacts between MM and immunosuppressive cells
# Create a list to hold the subsetted data frames and plots
CRS_int_cn_list_bor <- list()
plots_list_CRS_int_bor <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_CRS_int_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# List of rows to subset for Myeloma - SUP T cell contacts
CRS_int_row_names <- c("Myeloma", "Ki-67pos\nMyeloma", "PD-L1pos\nMyeloma", "M2-like\nmacrophage", "M1-like\nmacrophage", "ECs", "Classical\nmonocytes", "Nonclassical\nmonocytes", "Intermediate\nmonocytes")

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_list_bor)) {
  cn_df <- cn_list_bor[[i]]
  
  ### Subset the data frame (for Myeloma - SUP T cell contacts)
  subset_cn_df <- cn_df |>
    select("Myeloma", "Ki-67pos\nMyeloma", "PD-L1pos\nMyeloma", "M2-like\nmacrophage", "M1-like\nmacrophage", "ECs", "Classical\nmonocytes", "Nonclassical\nmonocytes", "Intermediate\nmonocytes") |>
    filter(rownames(cn_df) %in% CRS_int_row_names)
  
  # # # Remove links
  subset_cn_df[(1:3),(1:3)] <- NA
  subset_cn_df[(1:3),(7:9)] <- NA
  subset_cn_df[4,(1:9)] <- NA
  subset_cn_df[(7:9),(1:3)] <- NA
  subset_cn_df[(7:9),(7:9)] <- NA
  subset_cn_df[(5:6),(1:3)] <- NA
  
  # Add the subsetted data frame to the subset list
  CRS_int_cn_list_bor[[i]] <- subset_cn_df
  
  # Save the subset plot as a PNG file
  subset_file_name_svg <- paste0("CRS_int_chord_plot_", i, ".svg")
  subset_file_path_svg <- file.path(output_dir, subset_file_name_svg)
  print(paste("Saving CRS_int plot to:", subset_file_path_svg))  # Debugging line
  
  # Create the SVG for the subset with a pre-allocated track for customizing text
  svg(subset_file_path_svg, width = 9, height = 9)
  chordDiagram(as.matrix(subset_cn_df), scale = TRUE, grid.col = color_palette, 
               annotationTrack = c("grid"), preAllocateTracks = 1)
  
  # Customize sector labels (adjust font size)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name, 
                facing = "clockwise", niceFacing = T, adj = c(0, 0.7), cex = 1.2)  # Adjust 'cex' for font size
  }, bg.border = NA)
  
  # # Add the numbers (axis ticks) to the sectors
  # circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  #   circos.axis(h = "top", major.tick.length = 0.2, labels.cex = 0.6, direction = "outside")
  # }, bg.border = NA)
  
  circos.clear()
  dev.off()
  
  # Convert the SVG to a PNG using rsvg
  subset_file_name_png <- paste0("CRS_int_chord_plot_", i, ".png")
  subset_file_path_png <- file.path(output_dir, subset_file_name_png)
  
  rsvg_png(subset_file_path_svg, subset_file_path_png)  # Convert SVG to PNG
  
  # Load the saved image back as a grob using the correct path
  subset_img <- rasterGrob(png::readPNG(subset_file_path_png), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "CRS"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_CRS_int_bor[[i]] <- plot_with_subset_title
}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_CRS_int_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_CRS_int_bor, ncol = 5)
dev.off()



## Subsetting contacts between MM and immunosuppressive cells
# Create a list to hold the subsetted data frames and plots
CRS_int_cn_list_bor <- list()
plots_list_CRS_int_bor <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_CRS_int_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# List of rows to subset for Myeloma - SUP T cell contacts
CRS_int_row_names <- c("Myeloma", "M2-like\nmacrophage", "M1-like\nmacrophage", "ECs", "Classical\nmonocytes", "Nonclassical\nmonocytes", "Intermediate\nmonocytes")

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_list_bor)) {
  cn_df <- cn_list_bor[[i]]
  
  ### Subset the data frame (for Myeloma - SUP T cell contacts)
  subset_cn_df <- cn_df |>
    select("Myeloma", "M2-like\nmacrophage", "M1-like\nmacrophage", "ECs", "Classical\nmonocytes", "Nonclassical\nmonocytes", "Intermediate\nmonocytes") |>
    filter(rownames(cn_df) %in% CRS_int_row_names)
  
  # # # Remove links
  subset_cn_df[(1:3),(1:3)] <- NA
  subset_cn_df[(1:3),(5:7)] <- NA
  subset_cn_df[4,(1:7)] <- NA
  subset_cn_df[7,1] <- NA
  subset_cn_df[7,(4:7)] <- NA
  subset_cn_df[(5:6),(2:3)] <- NA
  
  # Add the subsetted data frame to the subset list
  CRS_int_cn_list_bor[[i]] <- subset_cn_df
  
  # Save the subset plot as a PNG file
  subset_file_name <- paste0("CRS_int_chord_plot_", i, ".png")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving CRS_int plot to:", subset_file_path))  # Debugging line
  
  # Create the PNG for the subset
  png(subset_file_path, width = 600, height = 600)
  chordDiagram(as.matrix(subset_cn_df), scale = TRUE, grid.col = color_palette)
  circos.clear()
  dev.off()
  
  # Load the saved image back as a grob using the correct path
  subset_img <- rasterGrob(png::readPNG(subset_file_path), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "CRS"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_CRS_int_bor[[i]] <- plot_with_subset_title
}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_CRS_int_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_CRS_int_bor, ncol = 5)
dev.off()


# # Test example - multiple effects of SUP cells
# test <- bl_border10cn12_func_int_CN0 |>
#   select(Myeloma, 'Ki-67pos\nMyeloma', 'PD-L1pos\nMyeloma', "CD4 Tnaive", "CD4 Tex", "pSTAT1neg\nCD4 Tmem", "pSTAT1pos\nCD4 Tmem", "Treg") |>
#   filter(rownames(bl_border10cn12_func_int_CN0) %in% c("CD4 Tnaive", "CD4 Tex", "pSTAT1neg\nCD4 Tmem", "pSTAT1pos\nCD4 Tmem", "Treg", "M2-like\nmacrophage", "M1-like\nmacrophage", "MDSCs"))
# 
# test[(4:8),(4:8)] <- NA
# 
# chordDiagram(as.matrix(test), scale = T, grid.col = color_palette)
# circos.clear()
# 
# test2 <- bl_border10cn12_func_int_CN0 |>
#   select(Myeloma, 'Ki-67pos\nMyeloma', 'PD-L1pos\nMyeloma')
#   filter(rownames(bl_border10cn12_func_int_CN0) %in% c("CD4 Tnaive", "CD4 Tex", "pSTAT1neg\nCD4 Tmem", "pSTAT1pos\nCD4 Tmem", "Treg", "M2-like\nmacrophage", "M1-like\nmacrophage", "MDSCs"))
# 
# test[(4:8),(4:8)] <- NA
# 
# chordDiagram(as.matrix(test2), scale =T, grid.col = color_palette)
# circos.clear()

# Kreisdiagrame / Balkendiagram Kontakte mit MM insgesamt aus MM perspektive
# Damit ranking
# Dann genaue Einteilung mit chord plots

bl_border10cn12_raw |>
  select(medium_cell_type)

medium_cell_type_labels

# # Ensure diagonal elements are set to zero
# for (cn in border10cn12_cn_values) {
#   df_name <- paste0("bl_border10cn12_func_int_CN", cn)
#   df <- get(df_name)
#   for (i in intersect(rownames(df), colnames(df))) {
#     df[i, i] <- 0
#   }
#   assign(df_name, df)
# }


# # Plot the chord diagram with the reordered matrix
# chordDiagram(as.matrix(bl_border10cn12_func_int_CN2), scale = TRUE, grid.col = color_palette, link.sort = F)
# circos.clear()
# 
# chordDiagram(as.matrix(bl_border10cn12_func_int_CN1), scale = TRUE, grid.col = color_palette, link.sort = F)
# circos.clear()
# 

## Contacts of myeloma cells in medium type

bl_border10cn12_med_int <- read.csv("full_results/border/baseline/border-10_CNs-12_interactions__medium_baseline.csv")

bl_border10cn12_med_int <- bl_border10cn12_med_int |> select(-X)
colnames(bl_border10cn12_med_int)<- gsub("\\.", " ", colnames(bl_border10cn12_med_int))
#colnames(bl_border10cn12_med_int)<- gsub("Ki 67pos Myeloma", "Ki-67pos\nMyeloma", colnames(bl_border10cn12_med_int))
#colnames(bl_border10cn12_med_int)<- gsub("PD L1pos Myeloma", "PD-L1pos\nMyeloma", colnames(bl_border10cn12_med_int))
colnames(bl_border10cn12_med_int)<- gsub("M1 like M  phi ", "M1-like\nmacrophage", colnames(bl_border10cn12_med_int))
colnames(bl_border10cn12_med_int)<- gsub("M2 like M  phi ", "M2-like\nmacrophage", colnames(bl_border10cn12_med_int))
colnames(bl_border10cn12_med_int) <- gsub("CD8 GZMB  Tmem", "GZMB+ CD8 Tmem", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("pSTAT1neg CD8 Tmem", "pSTAT1neg\nCD8 Tmem", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("pSTAT1pos CD8 GZMBpos Tmem", "pSTAT1pos\nCD8 GZMBpos Tmem", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("pSTAT1neg CD8 GZMBpos Tmem", "pSTAT1neg\nCD8 GZMBpos Tmem", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("pSTAT1neg CD4 Tmem", "pSTAT1neg\nCD4 Tmem", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("pSTAT1pos CD4 Tmem", "pSTAT1pos\nCD4 Tmem", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("pSTAT1neg CD56 dim NK", "pSTAT1neg\nCD56 dim NK", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("pSTAT1neg CD56 bright NK", "pSTAT1neg\nCD56 bright NK", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("pSTAT1pos CD56 dim NK", "pSTAT1pos\nCD56 dim NK", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("pSTAT1pos CD56 bright NK", "pSTAT1pos\nCD56 bright NK", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("Classical monocytes", "Classical\nmonocytes", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub("Nonclassical monocytes", "Nonclassical\nmonocytes", colnames(bl_border10cn12_med_int))
# colnames(bl_border10cn12_med_int) <- gsub( "Intermediate monocytes",  "Intermediate\nmonocytes", colnames(bl_border10cn12_med_int))

border10cn12_cn_values <- unique(bl_border10cn12_med_int$cn_celltypes)


for (cn in border10cn12_cn_values) {
  # Create a new data frame for each CN value
  new_df <- bl_border10cn12_med_int |>
    filter(cn_celltypes == cn) |>
    select(Adipocyte_Adipocyte:Treg_Treg) |>
    ungroup()|>
    pivot_longer(cols = everything(), names_to = "interaction", values_to = "value")|>
    group_by(interaction)|>
    summarize(mean = mean(value, na.rm = TRUE)) |>
    separate(interaction, into = c("partner1", "partner2"), sep = "_")|>
    pivot_wider(names_from = partner2, values_from = mean) |>
    tibble::column_to_rownames(var = "partner1")
  
  # Reorder rows and columns to match the order in the item_names
  matching_names <- intersect(item_med_int, rownames(new_df))
  new_df <- new_df[matching_names, matching_names]
  
  # Fill NA symmetrically
  for (i in seq_len(nrow(new_df))) {
    for (j in seq_len(ncol(new_df))) {
      if (is.na(new_df[i, j]) && !is.na(new_df[j, i])) {
        new_df[i, j] <- new_df[j, i]
      } else if (!is.na(new_df[i, j]) && is.na(new_df[j, i])) {
        new_df[j, i] <- new_df[i, j]
      }
    }
  }
  
  # Dynamically create a new variable name for the data frame
  assign(paste0("bl_border10cn12_med_int_CN", cn), new_df)
}


# Assuming the data frames are named bl_border10cn12_func_int_CN0 to bl_border10cn12_func_int_CN19
cn_med_list <- list()
for (i in 0:11) {
  # Dynamically assign the data frames to the list
  cn_med_list[[i+1]] <- get(paste0("bl_border10cn12_med_int_CN", i))
}


# bl_border10cn12_med_int_CN0
# 
bl_border10cn12_med_int_CN11 |>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_med_int_CN0)) |>
  mutate(celltype = reorder(celltype, -Myeloma)) |>
  ggplot(aes(x = "", y = Myeloma, fill = celltype)) +
  geom_bar(stat = "identity", width = 1, alpha = 0.5) +
  coord_polar("y", start = 0, direction = -1) +  # Start at 12 o'clock and go clockwise
  theme_void() +
  labs(fill = "Cell Type") +
  scale_fill_manual(values = color_palette_medium2)

## Code for calculating the numbers for the manuscript

bl_border10cn12_med_int_CN8|>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_med_int_CN0)) |>
  mutate(perc = Myeloma / sum(Myeloma, na.rm = TRUE)) |>
  arrange(desc(perc))


## Subsetting contacts between MM and immunosuppressive cells
# Create a list to hold the subsetted data frames and plots
plots_list_pie <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_MM_pie_int"
dir.create(output_dir, showWarnings = FALSE)

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_med_list)) {
  cn_df <- cn_med_list[[i]]
  
  # Generate the pie chart for each CN
  p <- cn_df |>
    select(Myeloma) |>
    mutate(celltype = rownames(bl_border10cn12_med_int_CN0)) |>
    mutate(perc = Myeloma / sum(Myeloma)) |>
    arrange(desc(perc)) |>
    mutate(celltype = reorder(celltype, perc)) |>
    slice_head(n = 10) |>
    ggplot(aes(x = celltype, y = perc, fill = celltype)) +
    geom_bar(stat = "identity", width = 1, alpha = 0.6) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.6)) +
    scale_x_discrete(labels = function(x) gsub("\n", " ", x)) +
    guides(fill = "none") +
    theme_classic() +
    labs(fill = "Cell Type", y = "Myeloma cell interactions [%]", x = "") +
    scale_fill_manual(values = color_palette_medium2) +
    coord_flip() +  # Flip coordinates to make bars horizontal
    theme(
      axis.text.x = element_text(
        angle = 45,          # Rotate text by 45 degrees
        hjust = 1,           # Horizontal justification
        vjust = 1,           # Vertical justification
        size = 10            # Adjust text size as needed
      )
    )
  
  
  # Save the subset plot as a PNG file
  subset_file_name <- paste0("MM_int_pie_chart_CN", i, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving CRS_int plot to:", subset_file_path))  # Debugging line
  
  # Save each pie chart as a PNG file
  ggsave(filename = subset_file_path, plot = p, width = 4, height = 4)
  
  # Modify the plot to remove the legend for the grid arrangement
  p_no_legend <- p + theme(legend.position = "none")
  
  # Add the plot to the list for later display in a grid
  plots_list_pie[[i]] <- p_no_legend
}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "MM_int_pie_charts_grid.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Adjust resolution as needed
grid.arrange(grobs = plots_list_pie, ncol = 4)  # Adjust number of columns based on preference
dev.off()


## Making a NEW GRID PLOT to show connections with MM cells on MEDIUM level
# Creating a table that contains all contacts from all CNs

MM_int_CN0 <- bl_border10cn12_med_int_CN0 |>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_med_int_CN0)) |>
  mutate(perc = Myeloma / sum(Myeloma)) |>
  arrange(desc(perc)) |>
  mutate(celltype = reorder(celltype, perc)) |>
  mutate(CN = "preCN0")

MM_int_CN8 <- bl_border10cn12_med_int_CN8 |>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_med_int_CN8)) |>
  mutate(perc = Myeloma / sum(Myeloma)) |>
  arrange(desc(perc)) |>
  mutate(celltype = reorder(celltype, perc)) |>
  mutate(CN = "preCN8")

MM_int_CN11 <- bl_border10cn12_med_int_CN11 |>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_med_int_CN11)) |>
  mutate(perc = Myeloma / sum(Myeloma)) |>
  arrange(desc(perc)) |>
  mutate(celltype = reorder(celltype, perc)) |>
  mutate(CN = "preCN11")

MM_int_CN9 <- bl_border10cn12_med_int_CN9 |>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_med_int_CN9)) |>
  mutate(perc = Myeloma / sum(Myeloma)) |>
  arrange(desc(perc)) |>
  mutate(celltype = reorder(celltype, perc)) |>
  mutate(CN = "preCN9")

MM_int_cnoi <- rbind(MM_int_CN0, MM_int_CN8, MM_int_CN11, MM_int_CN9)

MM_int_cnoi <- MM_int_cnoi|>
  colnames(bl_border10cn12_med_int)<- gsub("M1 like M  phi ", "M1-like\nmacrophage", colnames(bl_border10cn12_med_int))
colnames(bl_border10cn12_med_int)<- gsub("M2 like M  phi ", "M2-like\nmacrophage", colnames(bl_border10cn12_med_int))
colnames(bl_border10cn12_med_int) <- gsub("CD8 GZMB  Tmem", "GZMB+ CD8 Tmem", colnames(bl_border10cn12_med_int))

# Ensure celltype is a factor sorted alphabetically
MM_int_cnoi$celltype <- factor(MM_int_cnoi$celltype, 
                               levels = sort(unique(MM_int_cnoi$celltype)))

ggplot(MM_int_cnoi, aes(x = celltype, y = CN, color = celltype, size = perc)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = color_palette_medium2) +
  scale_x_discrete(position = "top", limits = order_medium2) +
  scale_size(range = c(1, 8)) +
  # Add evenly spaced y-breaks to create a grid.
  # Adjust these breaks as per your data range.
  guides(color = "none") +
  labs(x = "", y = "") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
 #   panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = -.05)
  )

## Making a NEW GRID PLOT to show connections with MM cells on MINOR level
# Creating a table that contains all contacts from all CNs

MM_int_minor_CN0 <- bl_border10cn12_minor_int_CN0 |>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_minor_int_CN0)) |>
  mutate(perc = Myeloma / sum(Myeloma)) |>
  arrange(desc(perc)) |>
  mutate(celltype = reorder(celltype, perc)) |>
  mutate(CN = "preCN0")

MM_int_minor_CN8 <- bl_border10cn12_minor_int_CN8 |>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_minor_int_CN8)) |>
  mutate(perc = Myeloma / sum(Myeloma)) |>
  arrange(desc(perc)) |>
  mutate(celltype = reorder(celltype, perc)) |>
  mutate(CN = "preCN8")

MM_int_minor_CN11 <- bl_border10cn12_minor_int_CN11 |>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_minor_int_CN11)) |>
  mutate(perc = Myeloma / sum(Myeloma)) |>
  arrange(desc(perc)) |>
  mutate(celltype = reorder(celltype, perc)) |>
  mutate(CN = "preCN11")

MM_int_minor_CN9 <- bl_border10cn12_minor_int_CN9 |>
  select(Myeloma) |>
  mutate(celltype = rownames(bl_border10cn12_minor_int_CN9)) |>
  mutate(perc = Myeloma / sum(Myeloma)) |>
  arrange(desc(perc)) |>
  mutate(celltype = reorder(celltype, perc)) |>
  mutate(CN = "preCN9")

MM_int_minor_cnoi <- rbind(MM_int_minor_CN0, MM_int_minor_CN8, MM_int_minor_CN11, MM_int_minor_CN9)

MM_int_cnoi <- MM_int_cnoi|>
  colnames(bl_border10cn12_med_int)<- gsub("M1 like M  phi ", "M1-like\nmacrophage", colnames(bl_border10cn12_med_int))
colnames(bl_border10cn12_med_int)<- gsub("M2 like M  phi ", "M2-like\nmacrophage", colnames(bl_border10cn12_med_int))
colnames(bl_border10cn12_med_int) <- gsub("CD8 GZMB  Tmem", "GZMB+ CD8 Tmem", colnames(bl_border10cn12_med_int))

# Ensure celltype is a factor sorted alphabetically
MM_int_cnoi$celltype <- factor(MM_int_cnoi$celltype, 
                               levels = sort(unique(MM_int_cnoi$celltype)))

p_MM_int_minor_grid <- ggplot(MM_int_minor_cnoi, aes(x = factor(celltype, levels = item_minor_int_grid), 
                              y = factor(CN, levels = c("preCN9","preCN11","preCN8", "preCN0")),
                              color = celltype, size = perc)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = color_palette_minor) +
  scale_x_discrete(position = "top") +
  scale_size(range = c(2, 10)) +
  # Add evenly spaced y-breaks to create a grid.
  # Adjust these breaks as per your data range.
  guides(color = "none") +
  labs(x = "", y = "", size = "Interactions [%]") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    #   panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = -.05)
  )

ggsave("MM_int_minor_grid.svg", p_MM_int_minor_grid, path = output_dir, width =12, height =2.8)
ggsave("MM_int_minor_grid_legend.svg", p_MM_int_minor_grid, path = output_dir, width =12, height =4)

###
### MINOR CELL TYPE BASED CHORD DIAGRAMS
###


bl_border10cn12_minor_int <- read.csv("full_results/border/baseline/border-10_CNs-12_interactions__minor_baseline.csv")

bl_border10cn12_minor_int <- bl_border10cn12_minor_int |> select(-X)
colnames(bl_border10cn12_minor_int)<- gsub("\\.", " ", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int)<- gsub("Ki 67pos Myeloma", "Ki-67pos\nMyeloma", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int)<- gsub("PD L1pos Myeloma", "PD-L1pos\nMyeloma", colnames(bl_border10cn12_minor_int))
colnames(bl_border10cn12_minor_int)<- gsub("M1 like M  phi ", "M1-like\nmacrophage", colnames(bl_border10cn12_minor_int))
colnames(bl_border10cn12_minor_int)<- gsub("M2 like M  phi ", "M2-like\nmacrophage", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1pos CD8 Tmem", "pSTAT1pos\nCD8 Tmem", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1neg CD8 Tmem", "pSTAT1neg\nCD8 Tmem", colnames(bl_border10cn12_minor_int))
colnames(bl_border10cn12_minor_int) <- gsub("CD8 GZMB  Tmem", "CD8 GZMB Tmem", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1neg CD8 GZMBpos Tmem", "pSTAT1neg\nCD8 GZMBpos Tmem", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1neg Non Th17 CD4 Tmem", "pSTAT1neg\nNon Th17 CD4 Tmem", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1pos Non Th17 CD4 Tmem", "pSTAT1pos\nNon Th17 CD4 Tmem", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1neg Th17 CD4 Tmem", "pSTAT1neg\nTh17 CD4 Tmem", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1pos Th17 CD4 Tmem", "pSTAT1pos\nTh17 CD4 Tmem", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1neg CD56 dim NK", "pSTAT1neg\nCD56 dim NK", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1neg CD56 bright NK", "pSTAT1neg\nCD56 bright NK", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1pos CD56 dim NK", "pSTAT1pos\nCD56 dim NK", colnames(bl_border10cn12_minor_int))
# colnames(bl_border10cn12_minor_int) <- gsub("pSTAT1pos CD56 bright NK", "pSTAT1pos\nCD56 bright NK", colnames(bl_border10cn12_minor_int))
colnames(bl_border10cn12_minor_int) <- gsub("Classical monocytes", "Classical\nmonocytes", colnames(bl_border10cn12_minor_int))
colnames(bl_border10cn12_minor_int) <- gsub("Nonclassical monocytes", "Nonclassical\nmonocytes", colnames(bl_border10cn12_minor_int))
colnames(bl_border10cn12_minor_int) <- gsub( "Intermediate monocytes",  "Intermediate\nmonocytes", colnames(bl_border10cn12_minor_int))

border10cn12_cn_values <- unique(bl_border10cn12_minor_int$cn_celltypes)

for (cn in border10cn12_cn_values) {
  # Create a new data frame for each CN value
  new_df <- bl_border10cn12_minor_int |>
    filter(cn_celltypes == cn) |>
    select(-cn_celltypes) |>
    ungroup()|>
    pivot_longer(cols = everything(), names_to = "interaction", values_to = "value")|>
    group_by(interaction)|>
    summarize(mean = mean(value, na.rm = TRUE)) |>
    separate(interaction, into = c("partner1", "partner2"), sep = "_")|>
    pivot_wider(names_from = partner2, values_from = mean) |>
    tibble::column_to_rownames(var = "partner1")
  
  # Reorder rows and columns to match the order in the item_names
  matching_names <- intersect(item_minor_int, rownames(new_df))
  new_df <- new_df[matching_names, matching_names]
  
  # Fill NA symmetrically
  for (i in seq_len(nrow(new_df))) {
    for (j in seq_len(ncol(new_df))) {
      if (is.na(new_df[i, j]) && !is.na(new_df[j, i])) {
        new_df[i, j] <- new_df[j, i]
      } else if (!is.na(new_df[i, j]) && is.na(new_df[j, i])) {
        new_df[j, i] <- new_df[i, j]
      }
    }
  }
  
  # Dynamically create a new variable name for the data frame
  assign(paste0("bl_border10cn12_minor_int_CN", cn), new_df)
}


# Create a list to hold all the data frames
cn_list_bor_minor <- list()
plots_list_bor_minor <- list()

# Assuming the data frames are named bl_border10cn12_func_int_CN0 to bl_border10cn12_func_int_CN19
for (i in 0:11) {
  # Dynamically assign the data frames to the list
  cn_list_bor_minor[[i+1]] <- get(paste0("bl_border10cn12_minor_int_CN", i))
}

# Specify the output directory
output_dir <- "border10cn12_minor_chord_plots"
dir.create(output_dir, showWarnings = FALSE)

# Iterate over the list of data frames and apply the operations
for (i in 1:length(cn_list_bor_minor)) {
  cn_df <- cn_list_bor_minor[[i]]
  
  # Update the diagonal elements to 0 (if needed)
  # for (cn in intersect(rownames(cn_df), colnames(cn_df))) {
  #   cn_df[cn, cn] <- 0
  # }
  
  # Convert to matrix
  cn_matrix <- as.matrix(cn_df)
  
  # Set diagonal to 0
  diag(cn_matrix) <- 0
  
  # Keep only lower triangle (set upper triangle to 0)
  cn_matrix[upper.tri(cn_matrix)] <- 0
  
  # Alternatively, keep only upper triangle (set lower triangle to 0)
  # cn_matrix[lower.tri(cn_matrix)] <- 0
  
  # # Convert back to tibble with row names as a column
  # cn_df_modified <- cn_matrix %>%
  #   as.data.frame() %>%
  #   rownames_to_column("RowName") %>%
  #   as_tibble()
  
  # Add the subsetted data frame to the subset list
  plots_list_bor_minor[[i]] <- cn_matrix
  
  # Save the subset plot as a PNG file
  subset_file_name_svg <- paste0("CN_int_minor", i, ".svg")
  subset_file_path_svg <- file.path(output_dir, subset_file_name_svg)
  print(paste("Saving Int plot to:", subset_file_path_svg))  # Debugging line
  
  # Create the SVG for the subset with a pre-allocated track for customizing text
  svg(subset_file_path_svg, width = 9, height = 9)
  chordDiagram(as.matrix(cn_matrix), scale = TRUE, grid.col = color_palette_minor, 
               annotationTrack = c("grid"), preAllocateTracks = 1)
  
  # Customize sector labels (adjust font size)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name, 
                facing = "clockwise", niceFacing = T, adj = c(0, 0.7), cex = 1.0)  # Adjust 'cex' for font size
  }, bg.border = NA)
  
  # # Add the numbers (axis ticks) to the sectors
  # circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  #   circos.axis(h = "top", major.tick.length = 0.2, labels.cex = 0.6, direction = "outside")
  # }, bg.border = NA)
  
  circos.clear()
  dev.off()
  
  # Convert the SVG to a PNG using rsvg
  subset_file_name_png <- paste0("CN_int_minor", i, ".png")
  subset_file_path_png <- file.path(output_dir, subset_file_name_png)
  
  rsvg_png(subset_file_path_svg, subset_file_path_png)  # Convert SVG to PNG
  
  # Load the saved PNG back as a grob
  subset_img <- rasterGrob(png::readPNG(subset_file_path_png), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "Int"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_bor_minor[[i]] <- plot_with_subset_title
  
}


# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_grid_minor.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_bor_minor, ncol = 5)
dev.off()



## Subsetting contacts between MM and immunosuppressive cells
# Create a list to hold the subsetted data frames and plots
MM_CD4_SUP_int_cn_list_bor_minor <- list()
plots_list_MM_CD4_SUP_int_bor_minor <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_MM_CD4_SUP_int_chord_plots_minor"
dir.create(output_dir, showWarnings = FALSE)

# List of rows to subset for Myeloma - SUP T cell contacts
MM_CD4_SUP_int_row_names_minor <- c("CD4 Tnaive", "CD4 Tex", "Non Th17 CD4 Tmem", "Th17 CD4 Tmem", 
                              "Treg", "M2-like\nmacrophage", "MDSCs")

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_list_bor_minor)) {
  cn_df <- cn_list_bor_minor[[i]]
  
  ### Subset the data frame (for Myeloma - SUP T cell contacts)
  subset_cn_df <- cn_df |>
    rownames_to_column("RowName") |>
    select("RowName", "Myeloma", "CD4 Tnaive", "CD4 Tex", 
           "Non Th17 CD4 Tmem", "Th17 CD4 Tmem", 
           "Treg", "M2-like\nmacrophage", "MDSCs") |>
    filter(RowName %in% MM_CD4_SUP_int_row_names_minor) |>
    arrange(factor(RowName, levels = c(
       "Non Th17 CD4 Tmem", "Th17 CD4 Tmem","CD4 Tnaive", "CD4 Tex", "Treg", "M2-like\nmacrophage", 
      "MDSCs", "Myeloma"
    ))) |>
    column_to_rownames("RowName")
  
  # Remove links
  subset_cn_df[(1:3),2] <- NA
  subset_cn_df[(1:3),(4:8)] <- NA
  subset_cn_df[4,(2:8)] <- NA
  subset_cn_df[(5:7),(6:8)] <- NA
  
  # Add the subsetted data frame to the subset list
  MM_CD4_SUP_int_cn_list_bor_minor[[i]] <- subset_cn_df
  
  # Save the subset plot as a PNG file
  subset_file_name_svg <- paste0("MM_CD4_SUP_int_chord_plot_minor_", i, ".svg")
  subset_file_path_svg <- file.path(output_dir, subset_file_name_svg)
  print(paste("Saving MM_CD4_SUP_int plot to:", subset_file_path_svg))  # Debugging line
  
  # Create the SVG for the subset with a pre-allocated track for customizing text
  svg(subset_file_path_svg, width = 9, height = 9)
  chordDiagram(as.matrix(subset_cn_df), scale = TRUE, grid.col = color_palette_minor, 
               annotationTrack = c("grid"), preAllocateTracks = 1)
  
  # Customize sector labels (adjust font size)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name, 
                facing = "clockwise", niceFacing = T, adj = c(0, 0.7), cex = 1.0)  # Adjust 'cex' for font size
  }, bg.border = NA)
  
  # # Add the numbers (axis ticks) to the sectors
  # circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  #   circos.axis(h = "top", major.tick.length = 0.2, labels.cex = 0.6, direction = "outside")
  # }, bg.border = NA)
  
  circos.clear()
  dev.off()
  
  # Convert the SVG to a PNG using rsvg
  subset_file_name_png <- paste0("MM_CD4_SUP_int_chord_plot_minor_", i, ".png")
  subset_file_path_png <- file.path(output_dir, subset_file_name_png)
  
  rsvg_png(subset_file_path_svg, subset_file_path_png)  # Convert SVG to PNG
  
  # Load the saved PNG back as a grob
  subset_img <- rasterGrob(png::readPNG(subset_file_path_png), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "MM-SUP-CD4"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_MM_CD4_SUP_int_bor_minor[[i]] <- plot_with_subset_title
  
}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_MM_CD4_SUP_int_grid_minor.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_MM_CD4_SUP_int_bor_minor, ncol = 5)
dev.off()


## Subsetting contacts between MM and immunosuppressive cells
# Create a list to hold the subsetted data frames and plots
MM_CD8_SUP_int_cn_list_bor_minor <- list()
plots_list_MM_CD8_SUP_int_bor_minor <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_MM_CD8_SUP_int_chord_plots_minor"
dir.create(output_dir, showWarnings = FALSE)

# List of rows to subset for Myeloma - SUP T cell contacts
MM_CD8_SUP_int_row_names_minor <- c("CD8 Tnaive", "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tex", "M2-like\nmacrophage", "MDSCs")

# Iterate over the list of full data frames and apply the subsetting
for (i in 1:length(cn_list_bor_minor)) {
  cn_df <- cn_list_bor_minor[[i]]
  
  ### Subset the data frame (for Myeloma - SUP T cell contacts)
  subset_cn_df <- cn_df |>
    rownames_to_column("RowName") |>
    select("RowName", "Myeloma", "CD8 Tnaive", "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tex", "M2-like\nmacrophage", "MDSCs") |>
    filter(RowName %in% MM_CD8_SUP_int_row_names_minor) |>
    arrange(factor(RowName, levels = c(
       "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tnaive", "CD8 Tex", "M2-like\nmacrophage", "MDSCs", "Myeloma"
    ))) |>
    column_to_rownames("RowName")
  
  #Setting order for chord plot
  desired_order <- c(
    "CD8 Tmem", 
    "CD8 GZMB Tmem", 
    "CD8 Tnaive", 
    "CD8 Tex", 
    "M2-like\nmacrophage", 
    "MDSCs", 
    "Myeloma"
  )
  
  # Remove links
  subset_cn_df[(1:3),(2:4)] <- 0
  subset_cn_df[(1:3),(6:7)] <- 0
  subset_cn_df[4,(2:7)] <- 0
  subset_cn_df[(5:6),(6:7)] <- 0
  
  
  # Add the subsetted data frame to the subset list
  MM_CD8_SUP_int_cn_list_bor_minor[[i]] <- subset_cn_df
  
  # Save the subset plot as an SVG file
  subset_file_name_svg <- paste0("MM_CD8_SUP_int_chord_plot_minor_", i, ".svg")
  subset_file_path_svg <- file.path(output_dir, subset_file_name_svg)
  print(paste("Saving MM_CD8_SUP_int plot to:", subset_file_path_svg))  # Debugging line
  
  # Create the SVG for the subset with a pre-allocated track for customizing text
  svg(subset_file_path_svg, width = 9, height = 9)
  chordDiagram(as.matrix(subset_cn_df), scale = TRUE, grid.col = color_palette_minor, 
               annotationTrack = c("grid"), preAllocateTracks = 1, order = desired_order)
  
  # Customize sector labels (adjust font size)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name, 
                facing = "clockwise", niceFacing = T, adj = c(0, 0.7), cex = 1.0)  # Adjust 'cex' for font size
  }, bg.border = NA)
  
  # # Add the numbers (axis ticks) to the sectors
  # circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  #   circos.axis(h = "top", major.tick.length = 0.2, labels.cex = 0.6, direction = "outside")
  # }, bg.border = NA)
  
  circos.clear()
  dev.off()
  
  # Convert the SVG to a PNG using rsvg
  subset_file_name_png <- paste0("MM_CD8_SUP_int_chord_plot_minor_", i, ".png")
  subset_file_path_png <- file.path(output_dir, subset_file_name_png)
  
  rsvg_png(subset_file_path_svg, subset_file_path_png)  # Convert SVG to PNG
  
  # Load the saved PNG back as a grob
  subset_img <- rasterGrob(png::readPNG(subset_file_path_png), interpolate = TRUE)
  
  # Add a title above each subset plot
  subset_title <- textGrob(paste("CN", i-1, "MM-SUP-CD8"), gp = gpar(fontsize = 10, fontface = "bold"))
  
  # Combine the title and subset plot into one grob
  plot_with_subset_title <- arrangeGrob(subset_title, subset_img, ncol = 1, heights = c(0.1, 0.9))
  
  # Store the subset plot in the list
  plots_list_MM_CD8_SUP_int_bor_minor[[i]] <- plot_with_subset_title
}

# Display all plots in a grid
output_grid_file <- file.path(output_dir, "chord_plots_MM_CD8_SUP_int_grid_minor.png")
png(output_grid_file, width = 4000, height = 3000, res = 600)  # Use a higher resolution
grid.arrange(grobs = plots_list_MM_CD8_SUP_int_bor_minor, ncol = 5)
dev.off()

###
### Calculating numbers for MANUSCRIPT TEXT
###

int_manuscript <- data.frame()
int_manuscript <- cn_list_bor_minor[[10]]

# write.csv(int_manuscript, "test_interaction_matrix.csv")


# # Loop through each element in cn_list_bor_minor
# for (i in seq_along(cn_list_bor_minor)) {
#   # Get the matrix from the list
#   int_manuscript <- cn_list_bor_minor[[i]]
#   
#   # Define the file name using the index
#   filename <- paste0("interaction_matrix_", i, ".csv")
#   
#   # Save the matrix as a CSV file
#   write.csv(int_manuscript, filename)
# }

## Calculation of MM-SUP-CD4 interaction numbers
int_manuscript <- int_manuscript|>
  rownames_to_column("RowName") |>
  select("RowName", "Myeloma", "CD4 Tnaive", "CD4 Tex", "Non Th17 CD4 Tmem", "Th17 CD4 Tmem", "Treg", "M2-like\nmacrophage", "MDSCs") |>
  filter(RowName %in% c("CD4 Tnaive", "CD4 Tex", "Non Th17 CD4 Tmem", "Th17 CD4 Tmem", "Treg", "M2-like\nmacrophage", "MDSCs", "Myeloma")) |>
  arrange(factor(RowName, levels = c(
    "Non Th17 CD4 Tmem", "Th17 CD4 Tmem","CD4 Tnaive", "CD4 Tex", "Treg", "M2-like\nmacrophage", 
    "MDSCs", "Myeloma"
  ))) |>
  column_to_rownames("RowName")

int_manuscript <- int_manuscript|>
rownames_to_column("RowName") |>
  select("RowName", "Myeloma", "CD8 Tnaive", "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tex", "Treg", "M2-like\nmacrophage", "MDSCs", "Myeloma") |>
  filter(RowName %in% c("Myeloma", "CD8 Tnaive", "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tex", "Treg", "M2-like\nmacrophage", "MDSCs", "Myeloma")) |>
  arrange(factor(RowName, levels = c(
    "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tnaive", "CD8 Tex", "Treg", "M2-like\nmacrophage", "MDSCs", "Myeloma"
  ))) |>
  column_to_rownames("RowName")

# 1) Identify the common names in both row and column names:
common_names <- intersect(rownames(int_manuscript), colnames(int_manuscript))

# 2) Reorder both the rows and columns of 'df' to have the same order
#    (so that row i corresponds to the same label as column i).
int_manuscript <- int_manuscript[common_names, common_names]

# 3) Convert to a matrix if not already (in case df has non-numeric columns)
df_int_manuscript <- as.matrix(int_manuscript)

# 4) Zero out the diagonal:
diag(df_int_manuscript) <- 0

# 5) Compute column sums and get column-wise percentages:
col_totals <- colSums(df_int_manuscript)
df_percent <- sweep(df_int_manuscript, 2, col_totals, FUN = "/") * 100

# 6) Look at the final percentage matrix
df_percent


# Remove links
int_manuscript[(1:3),2] <- 0
int_manuscript[(1:3),(4:5)] <- 0
int_manuscript[4,(2:8)] <- 0
int_manuscript[(5:7),(6:8)] <- 0

int_manuscript |>
  rownames_to_column("CellType") |>  # Convert rownames to a column
  #rowwise() |>  # Ensures row-wise operation
  mutate(across(-CellType, ~ .x / sum(c_across(-CellType)) * 100)) |>  # Calculate percentage per row, excluding the 'CellType' column
  ungroup() |>  # Removes row-wise grouping
  column_to_rownames("CellType")



int_manuscript <- data.frame()
int_manuscript <- cn_list_bor_minor[[10]]

## Calculation of MM-SUP-CD8 interaction numbers
int_manuscript <- int_manuscript|>
  rownames_to_column("RowName") |>
  select("RowName", "Myeloma", "CD8 Tnaive", "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tex", "M2-like\nmacrophage", "MDSCs") |>
  filter(RowName %in% c("CD8 Tnaive", "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tex", "M2-like\nmacrophage", "MDSCs")) |>
  arrange(factor(RowName, levels = c(
    "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tnaive", "CD8 Tex", "M2-like\nmacrophage", "MDSCs", "Myeloma"
  ))) |>
  column_to_rownames("RowName")

int_manuscript[(1:3),(2:4)] <- 0
int_manuscript[4,(2:7)] <- 0
int_manuscript[(5:6),(6:7)] <- 0

int_manuscript |>
  rownames_to_column("CellType") |>  # Convert rownames to a column
  rowwise() |>  # Ensures row-wise operation
  mutate(across(-CellType, ~ .x / sum(c_across(-CellType)) * 100)) |>  # Calculate percentage per row, excluding the 'CellType' column
  ungroup() |>  # Removes row-wise grouping
  column_to_rownames("CellType")

### CRS ANALYSIS ----
###

##
## Preparation of data
##
# survival_bl_border10_cn12 <- processed_data_border$border10_cn12

survival_bl_border10_cn12$CRS
survival_bl_border10_cn12$`CRP(mg/L)`
survival_bl_border10_cn12$`LDH(U/L)`
survival_bl_border10_cn12$LDH_log
survival_bl_border10_cn12$`SerumMSpike(g/dl)`
survival_bl_border10_cn12$`PlasmaCellCore(%)`

survival_bl_border10_cn12<-survival_bl_border10_cn12 |>
  mutate(il6_max_log = log(IL6_max),
         ferritin_max_log = log(Ferritin_max))

survival_bl_border10_cn12 <- survival_bl_border10_cn12 |>
  mutate(CRS_1 = ifelse(CRS > 0, 1, 0),
         CRS_2 = ifelse(CRS > 1, 1, 0),
         CRS_onset = as.numeric(difftime(CRS_start, infusion_date, units = "days")),
         CRS_event = ifelse(is.na(CRS_onset),0,1),
         CRS_onset_ci = ifelse(is.na(CRS_onset), OS_days, CRS_onset))|>
  ungroup()

summary(glm(CRS_1 ~ `CN0`, survival_bl_border10_cn12, family = binomial, na.action = na.omit))


# Define the directory to save the CRS plots
output_dir <- "border10cn12_crs"
dir.create(output_dir, showWarnings = FALSE)


##
## Log regression models for CRS development, and high CRS development
##

crs1_glm_border10cn12 <- data.frame(marker = character(),
                                 coefficient = numeric(),
                                 std_error = numeric(),
                                 p_value = numeric(),
                                 lower95 = numeric(),
                                 upper95 = numeric(),
                                 stringsAsFactors = FALSE)

for(variable in bl_border10_cn12_variables) {
  # Fit logistic regression model
  formula_str <- paste0("CRS_1 ~ `",variable,"`")
  
  #print(formula_str)
  model <- glm(formula_str, data = survival_bl_border10_cn12, family = binomial)
  
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
  crs1_glm_border10cn12 <- rbind(crs1_glm_border10cn12, marker_results)
}

crs1_glm_border10cn12 <- crs1_glm_border10cn12|>
  mutate(OR = exp(coefficient), FDR = p.adjust(p_value, method = "fdr"))

p_crs_border10cn12_glm_crs1 <- crs1_glm_border10cn12 |>
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

crs2_glm_border10cn12 <- data.frame(marker = character(),
                                 coefficient = numeric(),
                                 std_error = numeric(),
                                 p_value = numeric(),
                                 lower95 = numeric(),
                                 upper95 = numeric(),
                                 stringsAsFactors = FALSE)

for(variable in bl_border10_cn12_variables) {
  # Fit logistic regression model
  formula_str <- paste0("CRS_2 ~ `",variable,"`")
  
  #print(formula_str)
  model <- glm(formula_str, data = survival_bl_border10_cn12, family = binomial)
  
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
  crs2_glm_border10cn12 <- rbind(crs2_glm_border10cn12, marker_results)
}

crs2_glm_border10cn12 <- crs2_glm_border10cn12|>
  mutate(OR = exp(coefficient), FDR = p.adjust(p_value, method = "fdr"))

p_crs_border10cn12_glm_crs2 <- crs2_glm_border10cn12 |>
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
## CI curves for CRS development and high CRS development
##

# List to save plots for grid
p_crs1_border10cn12 <- list()
p_crs2_border10cn12 <- list()

# Specify the output directory for subsetted plots
output_dir <- "border10cn12_crs"
dir.create(output_dir, showWarnings = FALSE)


for (variable in bl_border10_cn12_variables) {
  
  # Directly use the formula in survfit2
  fit <- survfit(Surv(CRS_onset_ci, as.numeric(CRS_1)) ~ ifelse(get(variable, survival_bl_border10_cn12) > median(get(variable, survival_bl_border10_cn12)), 1, 0), data = survival_bl_border10_cn12, stype = 2)
  
  # Generate the plot
  p <- ggsurvplot(fit,
                  fun = "event",
                  ylab="CI CRS grade 1", xlim = c(0,20), ylim = c(0,1), xlab="Days after CAR-T infusion",
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
                  legend.title = c(paste0("pre",variable)),
                  palette = c("red2", "red4"))
  
  
  # Extract only the plot component from ggsurvplot
  crs_plot <- p$plot
  
  # Add the plot to the list
  p_crs1_border10cn12[[variable]] <- crs_plot
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_crs1_CN", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving PFS plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = crs_plot, width = 3, height = 2)
}

# Convert each plot in p_pfs_blborder10cn12 to a grob (grid graphical object)
grobs_crs1_blborder10cn12 <- lapply(p_crs1_border10cn12, ggplotGrob)

# Create the file path for saving the grid image
output_grid_file <- file.path(output_dir, "bl_border10cn12_crs1_grid.svg")

# Save the grid of survival plots as a PNG file
svg(output_grid_file, width = 16, height = 12)  # Adjust resolution and size
do.call(grid.arrange, c(grobs_crs1_blborder10cn12, ncol = 4))  # Arrange plots in a grid, 4 columns
dev.off()



for (variable in bl_border10_cn12_variables) {
  
  # Directly use the formula in survfit2
  fit <- survfit(Surv(CRS_onset_ci, as.numeric(CRS_2)) ~ ifelse(get(variable, survival_bl_border10_cn12) > median(get(variable, survival_bl_border10_cn12)), 1, 0), data = survival_bl_border10_cn12, stype = 2)
  
  # Generate the plot
  p <- ggsurvplot(fit,
                  fun = "event",
                  ylab="CI CRS grade 2", xlim = c(0,20), ylim = c(0,1), xlab="Days after CAR-T infusion",
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
                  legend.title = c(paste0("pre",variable)),
                  palette = c("red2", "red4"))
  
  
  # Extract only the plot component from ggsurvplot
  crs_plot <- p$plot
  
  # Add the plot to the list
  p_crs2_border10cn12[[variable]] <- crs_plot
  
  # Save the plot as a SVG file
  subset_file_name <- paste0("border10cn12_crs2_CN", variable, ".svg")
  subset_file_path <- file.path(output_dir, subset_file_name)
  print(paste("Saving OS plot to:", subset_file_path))  # Debugging line
  
  # Save the plot as a file
  ggsave(filename = subset_file_path, plot = crs_plot, width = 3, height = 2)
}

# Convert each plot in p_os_blborder10cn12 to a grob (grid graphical object)
grobs_crs2_blborder10cn12 <- lapply(p_crs2_border10cn12, ggplotGrob)

# Create the file path for saving the grid image
output_grid_file <- file.path(output_dir, "bl_border10cn12_crs2_grid.svg")

# Save the grid of survival plots as a PNG file
svg(output_grid_file, width = 16, height = 12)  # Adjust resolution and size
do.call(grid.arrange, c(grobs_crs2_blborder10cn12, ncol = 4))  # Arrange plots in a grid, 4 columns
dev.off()



##
## Distribution of areas between CRS groups
##

p_crs_border10cn12_crs1 <- survival_bl_border10_cn12|>
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


p_crs_border10cn12_crs2 <- survival_bl_border10_cn12|>
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

p_crs_border10cn12_toci <- survival_bl_border10_cn12|>
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


p_crs_border10cn12_steroid <- survival_bl_border10_cn12|>
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


p_crs_border10cn12_il6_max <- survival_bl_border10_cn12 |>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc") |>
  filter(area != "nan") |>
  ggplot(aes(x = il6_max_log, y = area_perc)) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkred') +
  geom_smooth(method = "lm", color = "darkorange") +  # Change line color to orange
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
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

p_crs_border10cn12_CN0_il6 <- survival_bl_border10_cn12 |>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc") |>
  filter(area == "X0") |>
  ggplot(aes(x = il6_max_log, y = area_perc)) +
  geom_jitter(size = 3, alpha = 0.3, width = 0.2, color = 'darkred') +
  geom_smooth(method = "lm", color = "darkorange") +  # Change line color to orange
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +  # Add Spearman correlation and p-value
  labs(x = "Max. IL-6 [log]", y = "CN0 area [norm.]", title = "CN0 correlation with IL6 max") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )

p_crs_border10cn12_CN1_il6 <- survival_bl_border10_cn12 |>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc") |>
  filter(area == "X1") |>
  ggplot(aes(x = il6_max_log, y = area_perc)) +
  geom_jitter(size = 3, alpha = 0.3, width = 0.2, color = 'darkred') +
  geom_smooth(method = "lm", color = "darkorange") +  # Change line color to orange
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "bottom", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +  # Add Spearman correlation and p-value
  labs(x = "Max. IL-6 [log]", y = "CN1 area [norm.]", title = "CN1 correlation with IL6 max") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )

ggsave("CN0_IL6max.svg", p_crs_border10cn12_CN0_il6, path = output_dir, width =3, height =3)
ggsave("CN1_IL6max.svg", p_crs_border10cn12_CN1_il6, path = output_dir, width =3, height =3)

p_crs_border10cn12_ferritin_max <- survival_bl_border10_cn12|>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc")|>
  filter(area != "nan")|>
  ggplot(aes(x=ferritin_max_log, y=area_perc))+
  geom_jitter(size = 2, alpha = 0.3, width = 0.2, color = 'darkred') +
  geom_smooth(method = "lm", color = "darkorange") +  # Change line color to orange
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
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

p_crs_border10cn12_CN0_ferritin <- survival_bl_border10_cn12|>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc")|>
  filter(area == "X0")|>
  ggplot(aes(x=ferritin_max_log, y=area_perc))+
  geom_jitter(size = 3, alpha = 0.3, width = 0.2, color = 'darkred') +
  geom_smooth(method = "lm", color = "darkorange") +  # Change line color to orange
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +  # Add Spearman correlation and p-value
  labs(x = "Max. ferritin [log]", y = "CN0 area [norm.]", title = "Correlation with Ferritin max") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )

p_crs_border10cn12_CN1_ferritin <- survival_bl_border10_cn12|>
  pivot_longer(cols = all_of(bl_border10_cn12_variables), names_to = "area", values_to = "area_perc")|>
  filter(area == "X1")|>
  ggplot(aes(x=ferritin_max_log, y=area_perc))+
  geom_jitter(size = 3, alpha = 0.3, width = 0.2, color = 'darkred') +
  geom_smooth(method = "lm", color = "darkorange") +  # Change line color to orange
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "bottom", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +  # Add Spearman correlation and p-value
  labs(x = "Max. ferritin [log]", y = "CN1 area [norm.]", title = "Correlation with Ferritin max") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),  
    strip.text = element_text(size = 10),     
    plot.title = element_text(size = 12, hjust = 0.5)  
  )

ggsave("CN0_ferritinmax.svg", p_crs_border10cn12_CN0_ferritin, path = output_dir, width =3, height =3)
ggsave("CN1_ferritinmax.svg", p_crs_border10cn12_CN1_ferritin, path = output_dir, width =3, height =3)


##
## Saving CRS plots
##

# Example usage
pattern_crs_border10cn12 <- "p_crs_border10cn12_[a-zA-Z]"
pattern2_crs_border10cn12 <- "p_crs_border10cn12_[a-zA-Z]_[a-zA-Z]"
plots_crs_border10cn12 <- mget(ls(pattern = pattern_crs_border10cn12))
plots2_crs_border10cn12 <- mget(ls(pattern = pattern2_crs_border10cn12))

# Define the directory to save the CRS plots
output_dir <- "border10cn12_crs"
dir.create(output_dir, showWarnings = FALSE)

# Save the plots to the specified directory
save_plot(plots_crs_border10cn12, output_dir, width =8, height =8)
save_plot(plots2_crs_border10cn12, output_dir, width =8, height =8)


##
## Clinical plots for CRS
##

### Analysis of Cytokine Release Syndrome
survival_bl_border10_cn12 |>
  group_by(CRS) |>
  count()

survival_bl_border10_cn12$CRS

# Calculate percentages
crs_counts <- survival_bl_border10_cn12 %>%
  group_by(CRS) %>%
  summarize(count = n()) %>%
  mutate(percentage = paste0(round((count / sum(count)) * 100, 1), "%"))

p_CRS_distritbution <- ggplot(survival_bl_border10_cn12, aes(x = "", y =  as.factor(CRS), fill = as.factor(CRS))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y", start = 0, direction = -1) +
  scale_fill_manual(values = c("darkgrey", "red3", "darkred")) +
  scale_y_discrete(breaks = c("0","1","2"))+
  guides(fill = "none")+
  # geom_text(data = crs_counts, aes(x = 1, y = count / 2 + cumsum(c(0, head(count, -1))), label = percentage), 
  #           color = "white", size = 5) +
  ggtitle("Maximum CRS grades") +
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("CRS_distribution.svg", p_CRS_distritbution, path = output_dir, height = 3, width = 3)

p_CRS_duration <- survival_bl_border10_cn12 |>
  filter(CRS != 0) |>
  ggplot(aes(y = CRS_duration, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc()+
  scale_color_manual(values = c("red3", "darkred"))+
  scale_fill_manual(values = c("red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "CRS duration [days]", x = "CRS grade")+
  coord_flip()+
  theme_classic()

ggsave("CRS_duration.svg", p_CRS_duration, path = output_dir, height = 2, width = 4.5)

p_CRS_il6max <- survival_bl_border10_cn12 |>
  ggplot(aes(y = `IL6_max`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(ref.group = "0", method = "wilcox_test", p.adjust.method = "fdr")+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "Max. IL-6 post CAR-T [pg/ul]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_il6max.svg", p_CRS_il6max, path = output_dir, height = 3.5, width = 3)

p_CRS_ferritinmax <- survival_bl_border10_cn12 |>
  ggplot(aes(y = `Ferritin_max`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(ref.group = "0", method = "wilcox_test", p.adjust.method = "fdr", y.position = 6000,
           tip.length = 0.01, step.increase = 0.02)+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  ylim(0,10000)+
  guides(color = "none", fill = "none")+
  labs(y = "Max. ferritin post CAR-T [pg/ul]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_ferritinmax.svg", p_CRS_ferritinmax, path = output_dir, height = 3.5, width = 3)

p_CRS_pcc <- survival_bl_border10_cn12 |>
  ggplot(aes(y = `PlasmaCellCore(%)`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(ref.group = "0", method = "wilcox_test", p.adjust.method = "fdr")+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "PC infiltration pre CAR-T [%]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_pcc.svg", p_CRS_pcc, path = output_dir, height = 3.5, width = 3)

p_CRS_mspike <- survival_bl_border10_cn12 |>
  ggplot(aes(y = `SerumMSpike(g/dl)`,x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(ref.group = "0", method = "wilcox_test", p.adjust.method = "fdr")+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "M-spike pre CAR-T [g/dl]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_mspike.svg", p_CRS_mspike, path = output_dir, height = 3.5, width = 3)

p_CRS_cn0 <-survival_bl_border10_cn12 |>
  ggplot(aes(y = `X0`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(ref.group = "0", method = "wilcox_test", p.adjust.method = "fdr")+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "Area CN0 [norm.]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_cn0.svg", p_CRS_cn0, path = output_dir, height = 3.5, width = 3)


p_CRS_cn1 <- survival_bl_border10_cn12 |>
  ggplot(aes(y = `X1`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc()+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "Area CN1 [norm.]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_cn1.svg", p_CRS_cn1, path = output_dir, height = 3.5, width = 3)

p_CRS_cn3 <-survival_bl_border10_cn12 |>
  ggplot(aes(y = `X3`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc()+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "Area CN3 [norm.]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_cn3.svg", p_CRS_cn3, path = output_dir, height = 3.5, width = 3)


p_CRS_cn9 <-survival_bl_border10_cn12 |>
  ggplot(aes(y = `X9`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc()+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "Area CN9 [norm.]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_cn9.svg", p_CRS_cn9, path = output_dir, height = 3.5, width = 3)

p_CRS_pdl1mm <- survival_bl_border10cn12_cp_func_ilr|>
  ggplot(aes(y = `PD-L1pos Myeloma`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(ref.group = "0", method = "wilcox_test", p.adjust.method = "fdr")+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "PD-L1pos myeloma [norm.]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_pdl1mm.svg", p_CRS_pdl1mm, path = output_dir, height = 3.5, width = 3)


p_CRS_ki67mm <-survival_bl_border10cn12_cp_func_ilr|>
  ggplot(aes(y = `Ki-67pos Myeloma`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(ref.group = "0", method = "wilcox_test", p.adjust.method = "fdr")+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "Ki-67pos myeloma  [norm.]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_ki67mm.svg", p_CRS_ki67mm, path = output_dir, height = 3.5, width = 3)

p_CRS_mm <-survival_bl_border10cn12_cp_func_ilr|>
  ggplot(aes(y = `Myeloma`, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(ref.group = "0", method = "wilcox_test", p.adjust.method = "fdr")+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "PD-L1neg/Ki-67neg myeloma [norm.]", x = "CRS grade")+
  theme_classic()

ggsave("CRS_mm.svg", p_CRS_mm, path = output_dir, height = 3.5, width = 3)

survival_bl_border10cn12_cp_func_ilr|>
  ggplot(aes(y = ``, x = as.factor(CRS), color = as.factor(CRS), fill = as.factor(CRS))) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2)+
  geom_pwc(ref.group = "0", method = "wilcox_test", p.adjust.method = "fdr")+
  scale_color_manual(values = c("darkgrey", "red3", "darkred"))+
  scale_fill_manual(values = c("darkgrey", "red3", "darkred"))+
  guides(color = "none", fill = "none")+
  labs(y = "PD-L1neg/Ki-67neg myeloma [norm.]", x = "CRS grade")+
  theme_classic()

##
## Logistic regression model for CRS2
##

#print(formula_str)

survival_bl_border10_cn12<-survival_bl_border10_cn12 |>
  mutate(mEasix = `LDH(U/L)`*`Creatinine(mg/dl)`/PLT_d0) 

x0_crs_model <- glm(CRS_1 ~ X1 + mEasix, data = survival_bl_border10_cn12, family = binomial)
summary(x0_crs_model)

summary(glm(CRS_1 ~ `LDH(U/L)`, data = survival_bl_border10_cn12, family = binomial))
summary(glm(CRS_1 ~ PLT_d0, data = survival_bl_border10_cn12, family = binomial))
summary(glm(CRS_1 ~ Ferritin_d0, data = survival_bl_border10_cn12, family = binomial))
summary(glm(CRS_1 ~ CRP_d0, data = survival_bl_border10_cn12, family = binomial))
summary(glm(CRS_1 ~ `Creatinine(mg/dl)`, data = survival_bl_border10_cn12, family = binomial))
summary(glm(CRS_1 ~ mEasix, data = survival_bl_border10_cn12, family = binomial))


# Extract coefficients, standard errors, and p-values, and confidence intervals
coef_summary <- summary(model)$coefficients[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
coinf <- exp(confint(model))


survival_bl_border10_cn12$CRS_1


## COMPOSITION OF CNs ASSOCIATED WITH CRS ----
# MYELOMA
p_new_ki67mm_crs <- myeloma_fractions|>
  filter(functional_minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "Myeloma subpopulation [%]", title = "Ki-67pos myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  theme_classic()

p_new_mmneg_crs <- myeloma_fractions|>
  filter(functional_minor_cell_type == "Myeloma")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "Myeloma subpopulation [%]", title = "Ki-67neg/PDL1neg myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  theme_classic()

p_new_pdl1mm_crs <- myeloma_fractions|>
  filter(functional_minor_cell_type == "PD-L1pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "Myeloma subpopulation [%]", title = "PDL1pos myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  theme_classic()

myeloma_fractions|>
  filter(functional_minor_cell_type == "PD-L1pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

## NEW COMPOSITION DATA
# CD4

#CD4 Tnaive
p_new_cd4naive_crs <- cd4_fractions|>
  filter(minor_cell_type == "CD4 Tnaive")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "CD4 Tnaive")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#CD Tex
p_new_cd4ex_crs <- cd4_fractions|>
  filter(minor_cell_type == "CD4 Tex")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "CD4 Tex")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# Non TH17 CD4 Tmem
p_new_cd4tmem_crs <- cd4_fractions|>
  filter(minor_cell_type == "Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "Non Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# Th17 CD4 Tmem
p_new_cd4th17_crs <- cd4_fractions|>
  filter(minor_cell_type == "Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


# Treg
p_new_cd4treg_crs <- cd4_fractions|>
  filter(minor_cell_type == "Treg")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,70), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "Treg")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


## NEW COMPOSITION DATA
# Effector CD4 FUNCTIONAL LEVEL including pSTAT1

# pSTAT1pos Non TH17 CD4 Tmem
p_new_cd4tmempstat1pos_crs <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Non Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1pos Th17 CD4 Tmem
p_new_cd4th17pstat1pos_crs <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1pos Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1pos Non TH17 CD4 Tmem
p_new_cd4tmempstat1neg_crs <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Non Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1neg Th17 CD4 Tmem
p_new_cd4th17pstat1neg_crs <- cd4f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,20), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD4 T cell subpopulation [%]", title = "pSTAT1neg Th17 CD4 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


## NEW COMPOSITION DATA - CD8
# CD8

#CD8 Tnaive
p_new_cd8naive_crs <- cd8_fractions|>
  filter(minor_cell_type == "CD8 Tnaive")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 Tnaive")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#CD8 Tex
p_new_cd8ex_crs <- cd8_fractions|>
  filter(minor_cell_type == "CD8 Tex")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 Tex")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

cd8_fractions|>
  filter(minor_cell_type == "CD8 Tex")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))
  
#CD8 Tmem
p_new_cd8tmem_crs <- cd8_fractions|>
  filter(minor_cell_type == "CD8 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#CD8 GZMB+ Tmem
p_new_cd8gzmb_crs <- cd8_fractions|>
  filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "CD8 GZMB+ Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

## NEW COMPOSITION DATA
# Effector CD8 FUNCTIONAL LEVEL including pSTAT1

# pSTAT1pos CD8 GZMBpos Tmem
p_new_cd8tmemgzmbpstat1pos_crs <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD8 GZMBpos Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos CD8 GZMBpos Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1pos CD8 Tmem
p_new_cd8tmempstat1pos_crs <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD8 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1pos GZMBpos CD8 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1pos CD8 GZMBpos Tmem
p_new_cd8tmemgzmbpstat1neg_crs <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD8 GZMBpos Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 GZMBpos Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# pSTAT1neg CD8 Tmem
p_new_cd8tmempstat1neg_crs <- cd8f_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD8 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN1"))+
  labs(x = "", y = "CD8 T cell subpopulation [%]", title = "pSTAT1neg CD8 Tmem")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


## NEW COMPOSITION DATA - MYELOID
# MYELOID

#Classical moncytes
p_new_classmono_crs_2 <- myeloid_fractions|>
  filter(minor_cell_type == "Classical monocytes")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Classical monocytes")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

myeloid_fractions|>
  filter(minor_cell_type == "Classical monocytes")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))


#Intermediate monocytes
p_new_itmmono_crs_2 <- myeloid_fractions|>
  filter(minor_cell_type == "Intermediate monocytes")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,30), expand = c(0,0), breaks = c(0,10,20,30,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Intermediate monocytes")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


myeloid_fractions|>
  filter(minor_cell_type == "Intermediate monocytes")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

#Nonclassical monocytes
p_new_nonclassmono_crs_2 <- myeloid_fractions|>
  filter(minor_cell_type == "Nonclassical monocytes")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Nonclassical monocytes")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

myeloid_fractions|>
  filter(minor_cell_type == "Nonclassical monocytes")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))

#MDSCs
p_new_mdsc_crs <- myeloid_fractions|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Myeloid/MDSCs")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.01,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

myeloid_fractions|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  group_by(cn_celltypes)|>
  summarise(median = median(perc))
  
#M1-like macrophage
p_new_m1mac_crs_2 <- myeloid_fractions|>
  filter(minor_cell_type == "M1-like macrophage")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,5), expand = c(0,0), breaks = c(0,1,2,3,4,5,6,8,10))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "M1-like macrophage")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#M2-like macrophage
p_new_m2mac_crs_2 <- myeloid_fractions|>
  filter(minor_cell_type == "M2-like macrophage")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "M2-like macrophage")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#MoDC
p_new_modc_crs_2 <- myeloid_fractions|>
  filter(minor_cell_type == "MoDC")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,25), expand = c(0,0), breaks = c(0,5,10,15,20,25))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "MoDC")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#MKs
p_new_mk_crs_2 <- myeloid_fractions|>
  filter(minor_cell_type == "MKs")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Myeloid cell subpopulation [%]", title = "Megakaryocytes")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

## NEW COMPOSITION DATA - STROMAL CELLS
# STROMAL CELLS

#Endothelial cells
p_new_ec_crs_2 <- stroma_fractions|>
  filter(minor_cell_type == "ECs")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Stroma cell subpopulation [%]", title = "ECs")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#Endothelial cells
p_new_adipocyte_crs_2 <- stroma_fractions|>
  filter(minor_cell_type == "Adipocyte")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "Stroma cell subpopulation [%]", title = "Adipocyte")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

## NEW COMPOSITION DATA - NK CELLS
# NK

#Dim NK cells
p_new_dimnk_crs_2 <- nk_fractions|>
  filter(minor_cell_type == "CD56 dim NK")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "CD56 dim NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#Bright NK cells
p_new_brightnk_crs_2 <- nk_fractions|>
  filter(minor_cell_type == "CD56 bright NK")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,140), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "CD56 bright NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()


## NEW COMPOSITION DATA - FUNCTIONAL NK CELLS
# FUNCTIONAL NK

#pSTAT1pos Dim NK cells
p_new_pstat1posdimnk_crs_2 <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD56 dim NK")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 dim NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#pSTAT1pos Bright NK cells
p_new_pstat1posbrightnk_crs_2 <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1pos CD56 bright NK")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1pos CD56 bright NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#pSTAT1neg Dim NK cells
p_new_pstat1negdimnk_crs <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD56 dim NK")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 dim NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#pSTAT1neg Bright NK cells
p_new_pstat1negbrightnk_crs <- nkf_fractions|>
  filter(functional_minor_cell_type == "pSTAT1neg CD56 bright NK")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "NK cell subpopulation [%]", title = "pSTAT1neg CD56 bright NK")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

## NEW COMPOSITION DATA - B CELLS
# B

#B cells
p_new_b_crs_2 <- b_fractions|>
  filter(minor_cell_type == "B cells")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "B cell subpopulation [%]", title = "B cells")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

#MBCs
p_new_mbcs_crs_2 <- b_fractions|>
  filter(minor_cell_type == "MBCs")|>
  filter(cn_celltypes %in% c("All", "0", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "0", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  scale_x_discrete(labels = c("All", "preCN0", "preCN1"))+
  labs(x = "", y = "B cell subpopulation [%]", title = "MBCs")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  guides(fill = "none", color = "none")+
  theme_classic()

# Retrieve the plots using the pattern
pattern_comp_bl_new_crs <- "p_new_[a-zA-Z0-9]+_crs"
plots_comp_bl_new_crs <- mget(ls(pattern = pattern_comp_bl_new_crs))
pattern_comp_bl_new_crs_2 <- "p_new_[a-zA-Z0-9]+_crs_2"
plots_comp_bl_new_crs_2 <- mget(ls(pattern = pattern_comp_bl_new_crs_2))

# Define the directory to save the plots
output_dir <- "blborder10cn12_comp_crs"
dir.create(output_dir, showWarnings = FALSE)


# Save the plots to the specified directory
save_plot(plots_comp_bl_new_crs, output_dir, width = 2.5, height = 3.5)
save_plot(plots_comp_bl_new_crs_2, output_dir, width = 3, height = 3.5)
