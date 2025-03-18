### PRE-POST COMPOSITION COMPARISON
###

myeloma_fractions_all_prepost_minor <- bl_border10cn12_raw_all |>
  group_by(DFCI_id, timepoint)|>
  count(minor_cell_type)|>
  #filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

myeloma_fractions_all_matched_minor <- myeloma_fractions_all_prepost_minor %>%
  group_by(DFCI_id) %>%
  filter(all(c("Pre", "Post") %in% timepoint)) %>%
  ungroup()

unique(myeloma_fractions_all_matched$DFCI_id)

myeloma_fractions_all_prepost_functional <- bl_border10cn12_raw_all |>
  group_by(DFCI_id, timepoint)|>
  count(functional_minor_cell_type)|>
  mutate(perc = n/sum(n)) |>
  mutate(cn_celltypes = "All") |>
  select(cn_celltypes, everything())|>
  ungroup()

myeloma_fractions_all_matched_functional <- myeloma_fractions_all_prepost_functional %>%
  group_by(DFCI_id) %>%
  filter(all(c("Pre", "Post") %in% timepoint)) %>%
  ungroup()
 
p_prepost_comp_minor <- myeloma_fractions_all_prepost_minor |>
  #filter(minor_cell_type == "M1-like macrophage") |>
  filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = minor_cell_type, fill = minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
  geom_point(aes(color = minor_cell_type, fill = minor_cell_type), size = 3, alpha = 0.5) +
  geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) +
  geom_pwc(method = "t.test",
           label = "p",
           step.increase = 0.1,
           bracket.nudge.y= -0.2
  )+
  scale_color_manual(values = color_palette_minor_2) +
  scale_fill_manual(values = color_palette_minor_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  facet_wrap(vars(minor_cell_type), scales = "free_y")


ggadjust_pvalue(p_prepost_comp_minor, p.adjust.method = "fdr", label = "p.format")

p_matched_comp_minor <- myeloma_fractions_all_matched_minor |>
  #filter(minor_cell_type == "M1-like macrophage") |>
  filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = minor_cell_type, fill = minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
  geom_point(aes(color = minor_cell_type, fill = minor_cell_type), size = 3, alpha = 0.5) +
  geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) +
  geom_pwc(method = "t.test",
           label = "p",
           step.increase = 0.1,
           bracket.nudge.y= -0.2
  )+
  scale_color_manual(values = color_palette_minor_2) +
  scale_fill_manual(values = color_palette_minor_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  facet_wrap(vars(minor_cell_type), scales = "free_y")


ggadjust_pvalue(p_matched_comp_minor, p.adjust.method = "fdr", label = "p.format")


p_prepost_comp_functional <- myeloma_fractions_all_prepost_functional |>
  complete(DFCI_id, timepoint, functional_minor_cell_type, fill = list(perc = 0))|>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = functional_minor_cell_type, fill = functional_minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
  geom_point(aes(color = functional_minor_cell_type, fill = functional_minor_cell_type), size = 3, alpha = 0.5) +
  geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) +
  geom_pwc(method = "t.test",
           label = "p",
           step.increase = 0.1,
           bracket.nudge.y= -0.2
  )+
  scale_color_manual(values = color_palette_2) +
  scale_fill_manual(values = color_palette_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  facet_wrap(vars(functional_minor_cell_type), scales = "free_y")


ggadjust_pvalue(p_prepost_comp_functional, p.adjust.method = "fdr", label = "p.adj.format")

# Define the directory to save the CRS plots
output_dir <- "border10cn12_prepost_comparisons"
dir.create(output_dir, showWarnings = FALSE)

p_matched_comp_functional_mm <- myeloma_fractions_all_matched_functional |>
  complete(DFCI_id, timepoint, functional_minor_cell_type, fill = list(perc = 0))|>
  filter(functional_minor_cell_type %in% c("Myeloma", "Ki-67pos Myeloma", "PD-L1pos Myeloma"))|>
  #filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = functional_minor_cell_type, fill = functional_minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
   geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) + 
  geom_point(aes(color = functional_minor_cell_type, fill = functional_minor_cell_type), size = 3, alpha = 0.5) +
  geom_pwc(method = "wilcoxon",
           label = "p",
           label.size = 3,
           step.increase = 0.1,
           bracket.nudge.y= -0.1
  )+
  scale_color_manual(values = color_palette_2) +
  scale_fill_manual(values = color_palette_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  theme(strip.background = element_blank(),strip.text = element_text(size = 10, face = "bold"))+
  facet_wrap(vars(functional_minor_cell_type), scales = "free_y")


p_matched_comp_functional_mm_adj <- ggadjust_pvalue(p_matched_comp_functional_mm, p.adjust.method = "fdr", label = "p.adj.format")

ggsave(filename = "p_matched_comp_functional_mm_adj.svg", plot = p_matched_comp_functional_mm_adj, path = output_dir,
       width = 5, height = 3)


p_matched_comp_minor_cd8 <- myeloma_fractions_all_matched_minor |>
  complete(DFCI_id, timepoint, minor_cell_type, fill = list(perc = 0))|>
  filter(minor_cell_type %in% c("CD8 Tnaive","CD8 GZMB+ Tmem", "CD8 Tmem", "CD8 Tex"))|>
  #filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = minor_cell_type, fill = minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
  geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) + 
  geom_point(aes(color = minor_cell_type, fill = minor_cell_type), size = 3, alpha = 0.5) +
  geom_pwc(method = "wilcoxon",
           label = "p",
           label.size = 3,
           step.increase = 0.1,
           bracket.nudge.y= -0.1
  )+
  scale_color_manual(values = color_palette_minor_2) +
  scale_fill_manual(values = color_palette_minor_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  theme(strip.background = element_blank(),strip.text = element_text(size = 10, face = "bold"))+
  facet_wrap(vars(minor_cell_type), scales = "free_y", nrow = 1)

p_matched_comp_minor_cd8_adj <- ggadjust_pvalue(p_matched_comp_minor_cd8, p.adjust.method = "fdr", label = "p.adj.format")

ggsave(filename = "p_matched_comp_minor_cd8_adj.svg", plot = p_matched_comp_minor_cd8_adj, path = output_dir,
       width = 6.5, height = 3)


p_matched_comp_minor_cd4 <- myeloma_fractions_all_matched_minor |>
  complete(DFCI_id, timepoint, minor_cell_type, fill = list(perc = 0))|>
  filter(minor_cell_type %in% c("CD4 Tnaive","Non Th17 CD4 Tmem", "Th17 CD4 Tmem", "CD4 Tex", "Treg"))|>
  #filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = minor_cell_type, fill = minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
  geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) +
  geom_point(aes(color = minor_cell_type, fill = minor_cell_type), size = 3, alpha = 0.5) +
  geom_pwc(method = "wilcoxon",
           label = "p",
           label.size = 3,
           step.increase = 0.1,
           bracket.nudge.y= -0.1
  )+
  scale_color_manual(values = color_palette_minor_2) +
  scale_fill_manual(values = color_palette_minor_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  theme(strip.background = element_blank(),strip.text = element_text(size = 10, face = "bold"))+
  facet_wrap(vars(minor_cell_type), scales = "free_y", nrow = 1)


p_matched_comp_minor_cd4_adj <- ggadjust_pvalue(p_matched_comp_minor_cd4, p.adjust.method = "fdr", label = "p.adj.format")

ggsave(filename = "p_matched_comp_minor_cd4_adj.svg", plot = p_matched_comp_minor_cd4_adj, path = output_dir,
       width = 7.5, height = 3)

p_matched_comp_minor_myeloid <- myeloma_fractions_all_matched_minor |>
  complete(DFCI_id, timepoint, minor_cell_type, fill = list(perc = 0))|>
  filter(minor_cell_type %in% c("Classical monocytes", "Intermediate monocytes", "MDSCs", "Nonclassical monocytes",
                                "M1-like macrophage", "M2-like macrophage", "MoDC", "MKs"))|>
  #filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = minor_cell_type, fill = minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
  geom_point(aes(color = minor_cell_type, fill = minor_cell_type), size = 3, alpha = 0.5) +
  geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) +
  geom_pwc(method = "wilcoxon",
           label = "p",
           label.size = 3,
           step.increase = 0.1,
           bracket.nudge.y= -0.1
  )+
  scale_color_manual(values = color_palette_minor_2) +
  scale_fill_manual(values = color_palette_minor_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  theme(strip.background = element_blank(),strip.text = element_text(size = 10, face = "bold"))+
  facet_wrap(vars(minor_cell_type), scales = "free_y", nrow = 2)


p_matched_comp_minor_myeloid_adj <- ggadjust_pvalue(p_matched_comp_minor_myeloid, p.adjust.method = "fdr", label = "p.adj.format")

ggsave(filename = "p_matched_comp_minor_myeloid_adj.svg", plot = p_matched_comp_minor_myeloid_adj, path = output_dir,
       width = 7.5, height = 5)


p_matched_comp_minor_stroma <- myeloma_fractions_all_matched_minor |>
  complete(DFCI_id, timepoint, minor_cell_type, fill = list(perc = 0))|>
  filter(minor_cell_type %in% c("ECs", "Adipocyte"))|>
  #filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = minor_cell_type, fill = minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
  geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) +
  geom_point(aes(color = minor_cell_type, fill = minor_cell_type), size = 3, alpha = 0.5) +
  geom_pwc(method = "wilcoxon",
           label = "p",
           label.size = 3,
           step.increase = 0.1,
           bracket.nudge.y= -0.1
  )+
  scale_color_manual(values = color_palette_minor_2) +
  scale_fill_manual(values = color_palette_minor_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  theme(strip.background = element_blank(),strip.text = element_text(size = 10, face = "bold"))+
  facet_wrap(vars(minor_cell_type), scales = "free_y")

p_matched_comp_minor_stroma_adj <- ggadjust_pvalue(p_matched_comp_minor_stroma, p.adjust.method = "fdr", label = "p.adj.format")

ggsave(filename = "p_matched_comp_minor_stroma_adj.svg", plot = p_matched_comp_minor_stroma_adj, path = output_dir,
       width = 4, height = 3)

p_matched_comp_minor_B <- myeloma_fractions_all_matched_minor |>
  complete(DFCI_id, timepoint, minor_cell_type, fill = list(perc = 0))|>
  filter(minor_cell_type %in% c("B cells", "MBCs"))|>
  #filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = minor_cell_type, fill = minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
  geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) +
  geom_point(aes(color = minor_cell_type, fill = minor_cell_type), size = 3, alpha = 0.5) +
  geom_pwc(method = "wilcoxon",
           label = "p",
           label.size = 3,
           step.increase = 0.1,
           bracket.nudge.y= -0.1
  )+
  scale_color_manual(values = color_palette_minor_2) +
  scale_fill_manual(values = color_palette_minor_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  theme(strip.background = element_blank(),strip.text = element_text(size = 10, face = "bold"))+
  facet_wrap(vars(minor_cell_type), scales = "free_y")


p_matched_comp_minor_B_adj <- ggadjust_pvalue(p_matched_comp_minor_B, p.adjust.method = "fdr", label = "p.adj.format")

ggsave(filename = "p_matched_comp_minor_B_adj.svg", plot = p_matched_comp_minor_B_adj, path = output_dir,
       width = 4, height = 3)

p_matched_comp_minor_nk <- myeloma_fractions_all_matched_minor |>
  complete(DFCI_id, timepoint, minor_cell_type, fill = list(perc = 0))|>
  filter(minor_cell_type %in% c("CD56 dim NK", "CD56 bright NK"))|>
  #filter(cn_celltypes != "NaN") |>
  mutate(timepoint = factor(timepoint, levels = c("Pre","Post")))|>
  group_by(timepoint)|>
  ggplot(aes(x = as.factor(timepoint), y = perc * 100)) +
  geom_boxplot(aes(color = minor_cell_type, fill = minor_cell_type), alpha = 0.4, outliers = FALSE, width = 0.4) +
  geom_line(aes(group = DFCI_id), color = "gray60", linetype = "dashed", size = 0.7) +
  geom_point(aes(color = minor_cell_type, fill = minor_cell_type), size = 3, alpha = 0.5) +
  geom_pwc(method = "wilcoxon",
           label = "p",
           label.size = 3,
           step.increase = 0.1,
           bracket.nudge.y= -0.1
  )+
  scale_color_manual(values = color_palette_minor_2) +
  scale_fill_manual(values = color_palette_minor_2) +
  #scale_y_continuous(expand = c(0,5))+
  labs(x = "", 
       y = "Cell population [%]", 
       title = "") +
  guides(fill = "none", color = "none") +
  theme_classic()+
  theme(strip.background = element_blank(),strip.text = element_text(size = 10, face = "bold"))+
  facet_wrap(vars(minor_cell_type), scales = "free_y")


p_matched_comp_minor_nk_adj <- ggadjust_pvalue(p_matched_comp_minor_nk, p.adjust.method = "fdr", label = "p.adj.format")

ggsave(filename = "p_matched_comp_minor_nk_adj.svg", plot = p_matched_comp_minor_nk_adj, path = output_dir,
       width = 4, height = 3)

# Example usage
pattern_matched_comp_adj <- "p_matched_comp_minor_[a-zA-Z]+_adj"
plots_matched_comp_adj <- mget(ls(pattern = pattern_matched_comp_adj))


# Define the directory to save the CRS plots
output_dir <- "border10cn12_prepost_comparisons"
dir.create(output_dir, showWarnings = FALSE)

# Save the plots to the specified directory
save_plot(plots_matched_comp_adj, output_dir, width =3, height =3)
