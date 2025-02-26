
###
### CN1 for CRS
###

## NEW COMPOSITION DATA
# MYELOMA
p_new_ki67mm_cn1 <- myeloma_fractions|>
  filter(functional_minor_cell_type == "Ki-67pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloma subpopulation [%]", title = "Ki-67pos myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  theme_classic()

p_new_mm_cn1 <- myeloma_fractions|>
  filter(functional_minor_cell_type == "Myeloma")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,120), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloma subpopulation [%]", title = "Myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = 0.05,
           step.increase = 0.08)+
  theme_classic()

p_new_pdl1mm_cn1 <- myeloma_fractions|>
  filter(functional_minor_cell_type == "PD-L1pos Myeloma")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=functional_minor_cell_type, fill=functional_minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_2)+
  scale_fill_manual(values = color_palette_2)+
  scale_y_continuous(limits = c(0,80), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloma subpopulation [%]", title = "PD-L1pos myeloma")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.25)+
  theme_classic()



## NEW COMPOSITION DATA
# CD4

#CD4 Tnaive
p_new_cd4naive_cn1 <- cd4_fractions|>
  filter(minor_cell_type == "CD4 Tnaive")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "CD4 Tnaive")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",bracket.nudge.y = -0.45)+
  theme_classic()

#CD Tex
p_new_cd4ex_cn1 <- cd4_fractions|>
  filter(minor_cell_type == "CD4 Tex")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "CD4 Tex")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
  theme_classic()

# Non TH17 CD4 Tmem
p_new_cd4tmem_cn1 <-cd4_fractions|>
  filter(minor_cell_type == "Non Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Non TH17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
  theme_classic()

# Th17 CD4 Tmem
p_new_cd4th17_cn1 <-cd4_fractions|>
  filter(minor_cell_type == "Th17 CD4 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,30), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "TH17 CD4 Tmem")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.7)+
  theme_classic()

# Treg
p_new_cd4treg_cn1 <-cd4_fractions|>
  filter(minor_cell_type == "Treg")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD4 T cell subpopulation [%]", title = "Treg")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.15)+
  theme_classic()


## NEW COMPOSITION DATA - CD8
#CD8 Tnaive
p_new_cd8naive_cn1 <-cd8_fractions|>
  filter(minor_cell_type == "CD8 Tnaive")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tnaive")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.05)+
  theme_classic()

#CD8 Tex
p_new_cd8ex_cn1 <-cd8_fractions|>
  filter(minor_cell_type == "CD8 Tex")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,50), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tex")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.1)+
  theme_classic()

#CD8 Tmem
p_new_cd8tmem_cn1 <-cd8_fractions|>
  filter(minor_cell_type == "CD8 Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 Tmem")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
  theme_classic()

#CD8 GZMB+ Tmem
p_new_cd8gzmb_cn1 <-cd8_fractions|>
  filter(minor_cell_type == "CD8 GZMB+ Tmem")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,70), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "CD8 T cell subpopulation [%]", title = "CD8 GZMB+ Tmem")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
  theme_classic()

## NEW COMPOSITION DATA - MYELOID
#Classical moncytes
p_new_classmono_cn1 <-myeloid_fractions|>
  filter(minor_cell_type == "Classical monocytes")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,60), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Classical monocytes")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.4)+
  theme_classic()

#Intermediate monocytes
p_new_itmmono_cn1 <-myeloid_fractions|>
  filter(minor_cell_type == "Intermediate monocytes")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,20), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Intermediate monocytes")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj", bracket.nudge.y = -0.1)+
  theme_classic()

#Nonclassical monocytes
p_new_nonclassmono_cn1 <-myeloid_fractions|>
  filter(minor_cell_type == "Nonclassical monocytes")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,50), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "Nonclassical monocytes")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
  theme_classic()

#MDSCs
p_new_mdsc_cn1 <-myeloid_fractions|>
  filter(minor_cell_type == "MDSCs")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,110), expand = c(0,0), breaks = c(0,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "MDSCs")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
  theme_classic()

#M1-like macrophage
p_new_m1mac_cn1 <-myeloid_fractions|>
  filter(minor_cell_type == "M1-like macrophage")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,3.5), expand = c(0,0), breaks = c(0,1,2,3,4,5,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M1-like macrophage")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
  theme_classic()


#M2-like macrophage
p_new_m2mac_cn1 <-myeloid_fractions|>
  filter(minor_cell_type == "M2-like macrophage")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,50), expand = c(0,0), breaks = c(0,10,20,30,40,50,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "M2-like macrophage")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj")+
  theme_classic()


#MoDC
p_new_modc_cn1 <-myeloid_fractions|>
  filter(minor_cell_type == "MoDC")|>
  filter(cn_celltypes %in% c("All", "1"))|>
  mutate(cn_celltypes = factor(cn_celltypes, levels = c("All", "1")))|>
  ggplot(aes(x=cn_celltypes, y = perc*100, color=minor_cell_type, fill=minor_cell_type))+
  geom_boxplot(alpha = 0.4, outliers = F, width = .4)+
  geom_jitter(size = 3, alpha = 0.5, width = 0.1)+
  scale_color_manual(values = color_palette_minor_2)+
  scale_fill_manual(values = color_palette_minor_2)+
  scale_y_continuous(limits = c(0,15), expand = c(0,0), breaks = c(0,5,10,15,20,40,60,80,100))+
  labs(x = "Cellular neighborhood", y = "Myeloid cell subpopulation [%]", title = "MoDC")+
  guides(fill = "none", color = "none")+
  geom_pwc(ref.group = "All", method = "wilcox.test", p.adjust.method = "fdr", label = "p.adj",
           step.increase = 0.1)+
  theme_classic()



# Retrieve the plots using the pattern
pattern_comp_bl_new_cn1 <- "p_new_[a-zA-Z0-9]+_cn1"
plots_comp_bl_new_cn1 <- mget(ls(pattern = pattern_comp_bl_new_cn1))

# Define the directory to save the plots
output_dir <- "blborder10cn12_comp_ASH"
dir.create(output_dir, showWarnings = FALSE)


# Save the plots to the specified directory
save_plot(plots_comp_bl_new_cn1, output_dir, width = 2, height = 4)


## Survival plot for ASH presentation


pfs_CN8 <- survfit(Surv(pfs_m, PFS_event) ~ ifelse(survival_bl_border10_cn12$CN8 > median(survival_bl_border10_cn12$CN8),"High", "Low"), data = survival_bl_border10_cn12)
pfs_CN8
summary(pfs_CN8, times = 12)

ggsurvplot(pfs_CN8,
           ylab = "Progression-free survival",
           xlab = "Months after CAR-T infusion",
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
           legend.title = "Baseline CN8",
           legend.labs = c("High", "Low"),
           palette = c("#000080","#76AAFF")
           )

pfs_CN9 <- survfit(Surv(pfs_m, PFS_event) ~ ifelse(survival_bl_border10_cn12$CN9 > median(survival_bl_border10_cn12$CN9),"High", "Low"), data = survival_bl_border10_cn12)
pfs_CN9
summary(pfs_CN9, times = 12)

ggsurvplot(pfs_CN9,
           ylab = "Progression-free survival",
           xlab = "Months after CAR-T infusion",
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
           legend.title = "Baseline CN9",
           legend.labs = c("High", "Low"),
           palette = c("#000080","#76AAFF")
)

survival_bl_border10_cn12 |> group_by(CRS) |> count()


CRS_CI_CN1 <- survfit(Surv(CRS_onset_ci, as.numeric(CRS_2)) ~ ifelse(survival_bl_border10_cn12$CN1 > median(survival_bl_border10_cn12$CN1),"High", "Low"), data = survival_bl_border10_cn12, stype = 2)

# Generate the plot
ggsurvplot(CRS_CI_CN1,
                fun = "event",
                ylab="Cumulative incidence of CRS grade 2", xlim = c(0,20), ylim = c(0,1), xlab="Days after CAR-T infusion",
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
                legend.labs = c("High", "Low"),
                legend.title = c("Baseline CN1"),
                palette = c("red2", "red4")
           )


## POST SURVIVAL

survival_bl_border10_cn12_post <- survival_bl_border10_cn12_post|>
  mutate(pfs_m = PFS_days/30.44)

pfs_CN11_post <- survfit(Surv(pfs_m, PFS_event) ~ ifelse(survival_bl_border10_cn12_post$CN11 > median(survival_bl_border10_cn12_post$CN11),"High", "Low"), data = survival_bl_border10_cn12_post)
pfs_CN11_post
summary(pfs_CN11_post, times = 12)

ggsurvplot(pfs_CN11_post,
           ylab = "Progression-free survival",
           xlab = "Months after CAR-T infusion",
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
           legend.title = "Baseline CN8",
           legend.labs = c("High", "Low"),
           palette = c("#000080","#76AAFF")
)


survival_bl_border10_cn12_post <- survival_bl_border10_cn12_post|>
  mutate(os_m = OS_days/30.44)

os_CN11_post <- survfit(Surv(os_m, OS_event) ~ ifelse(survival_bl_border10_cn12_post$CN11 > median(survival_bl_border10_cn12_post$CN11),"High", "Low"), data = survival_bl_border10_cn12_post)
os_CN11_post
summary(os_CN11_post, times = 12)

ggsurvplot(os_CN11_post,
           ylab = "Overall survival",
           xlab = "Months after CAR-T infusion",
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
           legend.title = "Baseline CN8",
           legend.labs = c("High", "Low"),
           palette = c("#009600","#88DC88")
)

### Survival based on product
pfs_product <- survfit(Surv(PFS_m, PFS_event) ~ survival_bl_border10_cn12$CAR.T.Type, data = survival_bl_border10_cn12)
pfs_product
summary(pfs_product, times = 12)

ggsurvplot(pfs_product,
           ylab = "Progression-free survival",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
           censor.size = 3,
           pval= TRUE, 
           pval.coord = c(0.5, 0.1), 
           pval.size = 3,
           size = 1.15,
           axes.offset = F,
           risk.table = T,
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
           legend.title = "CAR-T product",
           legend.labs = c("Cilta-cel", "Ide-cel"),
           palette = c("#000080","#76AAFF")
)


os_product <- survfit(Surv(OS_m, OS_event) ~ survival_bl_border10_cn12$CAR.T.Type, data = survival_bl_border10_cn12)
os_product
summary(os_product, times = 12)

ggsurvplot(os_product,
           ylab = "Overall survival",
           xlab = "Months after CAR-T infusion",
           break.time.by = 3,
           xlim = c(0,26),
           censor.size = 3,
           pval= TRUE, 
           pval.coord = c(0.5, 0.1), 
           pval.size = 3,
           size = 1.15,
           axes.offset = F,
           risk.table = T,
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
           legend.title = "CAR-T product",
           legend.labs = c("Cilta-cel", "Ide-cel"),
           palette = c("#009600","#88DC88")
)



###
### Test for pre-post composition comparison
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
