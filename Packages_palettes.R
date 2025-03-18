### Packages ----
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(survival)
library(data.table)
library(ggrepel)
library(ggpubr)
library(survminer)
library(openxlsx)
library(gridExtra)
library(compositions)
library(stringr)
library(circlize)
library(grid)
library(cmprsk)
library(svglite)
library(rsvg)
library(tibble)
library(lubridate)
library(prodlim)
library(ggpattern)
library(rstatix)
library(ggridges)
library(forcats)
library(tidytext)

### Color palettes ----
## Color palette for functional cell type annotation
color_palette <- c(
  "Classical\nmonocytes" = "#1f77b4",
  "Nonclassical\nmonocytes" = "#87b4d4",
  "Intermediate\nmonocytes" = '#1E90FF',
  "ECs" = "#ff7f0e",
  "M2-like\nmacrophage" = "#2ca02c",
  "M1-like\nmacrophage" = "#006400",
  "MDSCs" = "#FFD700",
  "pSTAT1neg\nNon Th17 CD4 Tmem" = '#9400D3',
  "pSTAT1pos\nNon Th17 CD4 Tmem" = '#8A2BE2' ,
  "pSTAT1neg\nTh17 CD4 Tmem" ='#4B0082',
  "pSTAT1pos\nTh17 CD4 Tmem" = '#6A0DAD' ,
  "CD4 Tnaive" = '#DA70D6',
  "CD4 Tex" = '#DDA0DD' ,
  "Treg" = '#D8BFD8',
  "B cells" = '#8B4513',
  "MBCs" = '#A0522D',
  "pSTAT1neg\nCD56 dim NK" = '#FF007F',
  "pSTAT1neg\nCD56 bright NK" = '#FF6EB4',
  "pSTAT1pos\nCD56 dim NK" = '#FF69B4',
  "pSTAT1pos\nCD56 bright NK" = '#FF1493',
  "CD8 Tnaive" = "#7f7f7f",
  "pSTAT1pos\nCD8 Tmem" = "#969696",
  "pSTAT1neg\nCD8 Tmem" = "#acacac",
  "pSTAT1pos\nCD8 GZMBpos Tmem" = "#c3c3c3",
  "pSTAT1neg\nCD8 GZMBpos Tmem" = "#dadada",
  "CD8 Tex" = "#0A0A0A",
  "Myeloma" = '#FF0000',
  "Ki-67pos\nMyeloma" = '#DC143C',
  "PD-L1pos\nMyeloma" = '#B22222',
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
  "MDSCs" = "#d62728",
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

## Color palette for medium cell type annotation
color_palette_medium <- c(
  "Adipocyte" = '#FFFDD0',
  "B cells" = '#8B4513',
  "CD4 Tex" = '#DDA0DD' ,
  "CD4 Tmem" = '#9400D3',
  "CD8 GZMB+ Tmem" = "#c3c3c3",
  "CD8 Tex" = "#0A0A0A",
  "CD8 Tmem" = "#969696",
  "ECs" = "#ff7f0e",
  "M1-like macrophage" = "#006400",
  "M2-like macrophage" = "#2ca02c",
  "MDSCs" = "#FFD700",
  "MoDC" = '#00CED1',
  "Monocytes" = "#1f77b4",
  "MKs" = "#17becf",
  "Myeloma" = '#FF0000',
  "NK cells" = '#FF007F',
  "Tnaive" = "#7f7f7f",
  "Treg" = '#D8BFD8'
)

color_palette_medium_new <- c(
  "Adipocyte" = '#FFFDD0',
  "B cells" = '#8B4513',
  "CD4 Tex" = '#DDA0DD' ,
  "CD4 Tmem" = '#9400D3',
  "CD8 GZMB+ Tmem" = "#c3c3c3",
  "CD8 Tex" = "#0A0A0A",
  "CD8 Tmem" = "#969696",
  "ECs" = "#ff7f0e",
  "M1-like macrophage" = "#006400",
  "M2-like macrophage" = "#2ca02c",
  "Myeloid/MDSCs" = "#FFD700",
  "MoDC" = '#00CED1',
  "Monocytes" = "#1f77b4",
  "MKs" = "#17becf",
  "Myeloma" = '#FF0000',
  "NK cells" = '#FF007F',
  "Tnaive" = "#7f7f7f",
  "Treg" = '#D8BFD8'
)

color_palette_medium2 <- c(
  "Monocytes" = "#1f77b4",
  "ECs" = "#ff7f0e",
  "M2-like\nmacrophage" = "#2ca02c",
  "M1-like\nmacrophage" = "#006400",
  "MDSCs" = "#FFD700",
  "CD4 Tmem" = '#9400D3',
  "CD4 Tex" = '#DDA0DD' ,
  "Treg" = '#D8BFD8',
  "B cells" = '#8B4513',
  "NK cells" = '#FF007F',
  "Tnaive" = "#7f7f7f",
  "CD8 Tmem" = "#969696",
  "CD8 GZMB+ Tmem" = "#c3c3c3",
  "CD8 Tex" = "#0A0A0A",
  "Myeloma" = '#FF0000',
  "MKs" = "#17becf",
  "MoDC" = '#00CED1',
  "Adipocyte" = '#FFFDD0'
)

order_medium2 <- c(
  "Adipocyte",
  "B cells",
  "CD4 Tex",
  "CD4 Tmem",
  "CD8 GZMB+ Tmem",
  "CD8 Tex",
  "CD8 Tmem",
  "ECs",
  "M1-like\nmacrophage",
  "M2-like\nmacrophage",
  "MDSCs",
  "MoDC",
  "Monocytes",
  "MKs",
  "Myeloma",
  "NK cells",
  "Tnaive",
  "Treg"
)


order_medium <- c(
  "Monocytes",
  "ECs",
  "M2-like macrophage",
  "M1-like macrophage",
  "MDSCs",
  "CD4 Tmem",
  "CD4 Tex",
  "Treg" ,
  "B cells",
  "NK cells",
  "Tnaive",
  "CD8 Tmem",
  "CD8 GZMB+ Tmem",
  "CD8 Tex",
  "Myeloma",
  "MKs",
  "MoDC",
  "Adipocyte"
)

item_med_int <- c(
  "Monocytes",
  "ECs",
  "M2-like\nmacrophage",
  "M1-like\nmacrophage",
  "MDSCs",
  "CD4 Tmem",
  "CD4 Tex",
  "Treg" ,
  "B cells",
  "NK cells",
  "Tnaive",
  "CD8 Tmem",
  "GZMB+ CD8 Tmem",
  "CD8 Tex",
  "Myeloma",
  "MKs",
  "MoDC",
  "Adipocyte"
)

# Create a vector of names from color_palette
item_names <- names(color_palette)

# Vector with med int names
item_med_int <- c(
  "Monocytes",
  "ECs",
  "M2-like\nmacrophage",
  "M1-like\nmacrophage",
  "MDSCs",
  "CD4 Tmem",
  "CD4 Tex",
  "Treg" ,
  "B cells",
  "NK cells",
  "Tnaive",
  "CD8 Tmem",
  "GZMB+ CD8 Tmem",
  "CD8 Tex",
  "Myeloma",
  "MKs",
  "MoDC",
  "Adipocyte"
)


# Vector with minor int names
item_minor <- c("Adipocyte", "B cells", "CD4 Tex", "CD4 Tnaive", "CD56 bright NK",
                "CD56 dim NK", "CD8 GZMB Tmem", "CD8 Tex", "CD8 Tmem", "CD8 Tnaive",
                "Classical monocytes", "ECs", "Intermediate monocytes", "M1-like macrophage",
                "M2-like macrophage", "MBCs", "MDSCs", "MKs", "MoDC", "Myeloma",
                "Non Th17 CD4 Tmem", "Nonclassical monocytes", "Th17 CD4 Tmem", "Treg")

item_minor_int <- c("Adipocyte", "B cells", "MBCs", "Non Th17 CD4 Tmem", "Th17 CD4 Tmem", "CD4 Tex",  "CD4 Tnaive", "Treg", 
                    "CD56 bright NK", "CD56 dim NK", "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tex",  "CD8 Tnaive",
                    "Classical\nmonocytes", "Intermediate\nmonocytes", "Nonclassical\nmonocytes", "ECs", 
                    "M1-like\nmacrophage", "M2-like\nmacrophage", "MDSCs", "MKs", "MoDC", "Myeloma")

item_minor_int_grid <- c("Myeloma", "B cells", "MBCs", "Non Th17 CD4 Tmem", "Th17 CD4 Tmem", "CD4 Tex",  "CD4 Tnaive", "Treg", 
                     "CD8 Tmem", "CD8 GZMB Tmem", "CD8 Tex",  "CD8 Tnaive", "CD56 bright NK", "CD56 dim NK",
                    "Classical\nmonocytes", "Intermediate\nmonocytes", "Nonclassical\nmonocytes",
                    "M1-like\nmacrophage", "M2-like\nmacrophage", "MDSCs", "MKs", "MoDC", "Adipocyte", "ECs" )

order_minor <- c("Adipocyte", "B cells", "MBCs","CD4 Tex", "CD4 Tnaive", "Th17 CD4 Tmem", "Non Th17 CD4 Tmem", "CD8 Tex", "CD8 Tnaive",
                   "CD8 Tmem", "CD8 GZMB+ Tmem", "Treg", "CD56 bright NK", "CD56 dim NK", 
                "Classical monocytes", "Intermediate monocytes","Nonclassical monocytes","ECs",  "M1-like macrophage",
                "M2-like macrophage",  "MDSCs", "MKs", "MoDC", "Myeloma")

## Color palette for minor cell type annotation
color_palette_minor <- c(
  "Classical\nmonocytes" = "#1f77b4",
  "Nonclassical\nmonocytes" = "#87b4d4",
  "Intermediate\nmonocytes" = '#1E90FF',
  "ECs" = "#ff7f0e",
  "M2-like\nmacrophage" = "#2ca02c",
  "M1-like\nmacrophage" = "#006400",
  "MDSCs" = "#FFD700",
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

color_palette_minor_2 <- c(
  "All" = "black",
  "Classical monocytes" = "#1f77b4",
  "Nonclassical monocytes" = "#87b4d4",
  "Intermediate monocytes" = '#1E90FF',
  "ECs" = "#ff7f0e",
  "M2-like macrophage" = "#2ca02c",
  "M1-like macrophage" = "#006400",
  "MDSCs" = "#FFD700",
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
  "CD8 GZMB+ Tmem" = "#c3c3c3",
  "CD8 Tex" = "#0A0A0A",
  "Myeloma" = '#FF0000',
  "MKs" = "#17becf",
  "MoDC" = '#00CED1',
  "Adipocyte" = '#FFFDD0'
)
