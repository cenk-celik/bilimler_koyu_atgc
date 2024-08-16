rm(list = ls())
# gerekli kutuphaneleri cagir
library(ELMER)
library(ComplexHeatmap)
library(ggpubr)

proje_klasoru <- "~/bilimler_koyu_2024"
setwd(proje_klasoru)

# Metilasyon verisini ice aktar:
load("ham_veri/metilasyon.RData")

# Hizlandirmak icin yalnizca chr8 uzerinde calis:
metilasyon_verisi_chr8 <- subset(metilasyon_verisi, subset = as.character(seqnames(metilasyon_verisi)) %in% c("chr8"))

# NA olan problari kaldir:
metilasyon_verisi_chr8 <- metilasyon_verisi_chr8[rowSums(is.na(assay(metilasyon_verisi_chr8))) == 0, ]

# metilasyon veri cercevesi olustur:
metilasyon_veri_cercevesi <- data.frame(
  "Sample.mean" = colMeans(assay(metilasyon_verisi_chr8), na.rm = TRUE),
  "groups" = metilasyon_verisi_chr8$definition
)

ggboxplot(
  data = metilasyon_veri_cercevesi,
  y = "Sample.mean",
  x = "groups",
  color = "groups",
  add = "jitter",
  ylab = "Ortalama DNA metilasyonu",
  xlab = ""
) + stat_compare_means()

## Diferansiyal metilasyon analizi----------------------------------------------
diferansiyel_metilasyon_sonuc <- TCGAanalyze_DMC(
  data = metilasyon_verisi_chr8,
  groupCol = "definition",
  group1 = "Primary solid Tumor",
  group2 = "Solid Tissue Normal",
  p.cut = 0.05,
  diffmean.cut = 0.15,
  save = FALSE,
  legend = "State",
  plot.filename = "metilasyon_yanardag.png",
  cores = 1
)

## Metilasyon isi haritasi------------------------------------------------------
klinik_veri <- as.data.frame(colData(metilasyon_verisi_chr8))

# Hiper ya da hipo metilasyon bolgeleri:
metilasyon_durumu <- "status"
problar <- rownames(diferansiyel_metilasyon_sonuc)[grep("hypo|hyper", diferansiyel_metilasyon_sonuc$status, ignore.case = TRUE)]
degisken_metilasyon <- metilasyon_verisi_chr8[problar, ]

# klinik veriyi metilasyon verisi ile ayni sekilde siralayacagiz:
sirali_klinik_veri <- klinik_veri[match(substr(colnames(degisken_metilasyon), 1, 12), klinik_veri$patient), ]

# ust lejant bilgisi:
ust_lejant <- HeatmapAnnotation(
  df = sirali_klinik_veri[, c("sample_type", "race")]
)

# sag lejant bilgisi:
sag_lejant <- rowAnnotation(
  df = diferansiyel_metilasyon_sonuc[problar, metilasyon_durumu],
  col = list("status.Primary.solid.Tumor.Solid.Tissue.Normal" = c( "Hypomethylated" = "orange", "Hypermethylated" = "darkgreen"))
)

# isi haritasi:
heatmap  <- Heatmap(
  matrix = assay(degisken_metilasyon),
  name = "DNA methylation",
  show_row_names = FALSE,
  cluster_rows = TRUE,
  right_annotation = sag_lejant,
  top_annotation = ust_lejant,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  column_title = "DNA Metilasyonu"
) 

draw(heatmap, annotation_legend_side =  "bottom")

## Motif analizi----------------------------------------------------------------
library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg19)
options(install.packages.compile.from.source = "always")
library(motifStack)
library(SummarizedExperiment)
library(dplyr)

problar <- rowRanges(metilasyon_verisi_chr8)[rownames(diferansiyel_metilasyon_sonuc)[grep("hypo|hyper", diferansiyel_metilasyon_sonuc$status, ignore.case = TRUE)],]

sekanslar <- GRanges(
  seqnames = as.character(seqnames(problar)),
  IRanges(
    start = ranges(problar) %>% as.data.frame() %>% dplyr::pull("start") - 100,
    end = ranges(problar) %>% as.data.frame() %>% dplyr::pull("end") + 100),
  strand = "*"
)

# motif arayisi:
gadem <- GADEM(sekanslar, verbose = TRUE, genome = Hsapiens)

# bulunan motif sayisi:
nMotifs(gadem)

# bulunan motiflerin gorulme sikligi:
nOccurrences(gadem)

# tum sekans motiflerini goruntule:
consensus(gadem)

# motifi goruntule:
motif <- getPWM(gadem)
pfm <- new("pfm", mat = motif[[1]], name = "Novel Site 1")
plotMotifLogo(pfm)

sessionInfo()