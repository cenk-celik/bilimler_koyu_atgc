rm(list = ls())

library(tidyverse)
library(ComplexHeatmap)

proje_klasoru <- "~/bilimler_koyu_2024"
setwd(proje_klasoru)

kopya_sayisi <- read.csv("brca_tcga_pan_can_atlas_2018/data_cna.txt", header = T, sep = "\t")
dim(kopya_sayisi)
kopya_sayisi[1:5, 1:5]

# gereksiz sutunlari kaldir:
kopya_sayisi$Entrez_Gene_Id <- NULL

# kopya sayisi matrisi olustur:
kopya_sayisi_matris <- kopya_sayisi %>%
  dplyr::select(-c(1)) %>% # Hugo_Symbol haric
  as.matrix()

# gen isimlerini matris satir isimleri olarak degistir:
rownames(kopya_sayisi_matris) <- kopya_sayisi$Hugo_Symbol

# 100'den fazla hastada kazanc (gain) gosteren genler:
kazanc <- which(rowSums(kopya_sayisi_matris == 1) > 100)[1:20]

# 100'den fazla hastada kayip (loss) gosteren genler:
kayip <- which(rowSums(kopya_sayisi_matris == -1) > 100)[1:20]

# iki listeyi birlestir:
genler <- c(kazanc, kayip)

# isi haritasi:
setwd("figurler")
pdf("kopya_sayisi_isi_haritasi.pdf", width = 10, height = 15)
Heatmap(kopya_sayisi_matris[genler, ],
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        col = circlize::colorRamp2(c(-1,0,1), colors = c("orange", "white", "purple")))
dev.off()
