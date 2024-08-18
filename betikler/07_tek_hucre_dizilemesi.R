rm(list = ls())

# Tek hucre dizilemesi analizine giris------------------------------------------
# kutuphaneleri cagir:
library(Seurat)
library(tidyverse)

# veriyi kaynagindan indir:
link <- "https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Breast_Cancer_3p/Breast_Cancer_3p_filtered_feature_bc_matrix.h5"
curl::curl_download(link, destfile = "ham_veri/Breast_Cancer_3p_filtered_feature_bc_matrix.h5")

# veriyi ice aktar:
meme_kanseri_verisi <- Read10X_h5("ham_veri/Breast_Cancer_3p_filtered_feature_bc_matrix.h5")
# Seurat nesnesi olustur:
meme_kanseri <- CreateSeuratObject(meme_kanseri_verisi, project = "meme_kanseri", 
                                   min.cells = 3, 
                                   min.features = 200)

meme_kanseri
# meta veri:
head(meme_kanseri@meta.data)

# gereksiz verileri kaldir:
rm(meme_kanseri_verisi)

## Kalite-kontrol---------------------------------------------------------------
# mitakondrial genler:
meme_kanseri[["percent.mt"]] <- PercentageFeatureSet(meme_kanseri, pattern = "^MT-")

# keman grafigi:
VlnPlot(meme_kanseri, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0.1,
        ncol = 3) &
  xlab("") & theme(aspect.ratio = 1) -> keman_once

# sacma grafigi:
FeatureScatter(meme_kanseri, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  FeatureScatter(meme_kanseri, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") &
  NoLegend() & theme(aspect.ratio = 1) -> sacma_once

# alt kumeleme:
meme_kanseri <- subset(meme_kanseri, 
                       subset = nFeature_RNA > 200 & 
                         nFeature_RNA < 7000 & 
                         nCount_RNA < 50000 &
                         percent.mt < 15 &
                         percent.mt > 2)

VlnPlot(meme_kanseri, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0.1,
        ncol = 3) &
  xlab("") & theme(aspect.ratio = 1) -> keman_sonra

FeatureScatter(meme_kanseri, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  FeatureScatter(meme_kanseri, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") &
  NoLegend() & theme(aspect.ratio = 1) -> sacma_sonra

# butun grafikleri tek panelde goster
cowplot::plot_grid(keman_once, sacma_once, keman_sonra, sacma_sonra, align = "h", nrow = 2, ncol = 2)

# normalizasyon:
meme_kanseri <- NormalizeData(meme_kanseri)

# yuksek degiskenlik gosteren genler
meme_kanseri <- FindVariableFeatures(meme_kanseri, 
                                     selection.method = "vst",
                                     nfeatures = 2000)

# en cok degiskenlik gosteren 10 gen
edg_genler <- head(VariableFeatures(meme_kanseri), 10)

# plot variable features with and without labels
LabelPoints(plot = VariableFeaturePlot(meme_kanseri), points = edg_genler, repel = TRUE) + 
  NoLegend() + 
  theme(aspect.ratio = 1)

## Hucre dongusu fazlarini hesapla:
s_fazi_genleri <- cc.genes.updated.2019$s.genes
g2_m_fazi_genleri <- cc.genes.updated.2019$g2m.genes
# hucre dongusu skoru:
meme_kanseri <- CellCycleScoring(meme_kanseri, 
                                 s.features = s_fazi_genleri,
                                 g2m.features = g2_m_fazi_genleri)

colnames(meme_kanseri@meta.data)
meme_kanseri@meta.data[1:5, c("S.Score", "G2M.Score", "Phase")]
meme_kanseri$CC.Difference <- meme_kanseri$S.Score - meme_kanseri$G2M.Score

# genleri olcekle:
## Her genin ifadesini kaydırarak, hücreler arasındaki ortalama ifade değerini 0 yapar.
## Her genin ifadesini ölçekleyerek, hücreler arasındaki varyansı 1 yapar.
## Bu adım, yüksek ifadesi olan genlerin analizlerde baskın hale gelmemesi için tüm genlere eşit ağırlık verir.
meme_kanseri <- ScaleData(meme_kanseri)

## Boyut indirgeme--------------------------------------------------------------
# Temel bilesen analizi (PCA):
meme_kanseri <- RunPCA(meme_kanseri)

# ilk iki bileseni grafik olarak incele:
DimPlot(meme_kanseri, reduction = "pca", group.by = "Phase") + theme(aspect.ratio = 1)

## istenmeyen varyans kaynaklarini kaldir:
meme_kanseri <- ScaleData(meme_kanseri, vars.to.regress = c("percent.mt", "CC.Difference"))

# Temel bilesen analizi
meme_kanseri <- RunPCA(meme_kanseri)

# ilk iki bileseni grafik olarak incele:
DimPlot(meme_kanseri, reduction = "pca", group.by = "Phase") + theme(aspect.ratio = 1)

# Birinci bileseni yuruten genler:
DimHeatmap(meme_kanseri, dims = 1, cells = 500, balanced = TRUE)

# ilk 15 bileseni yuruten genler:
DimHeatmap(meme_kanseri, dims = 1:15, cells = 500, balanced = TRUE)

# Kac bilesen kullanilmali?
ElbowPlot(meme_kanseri)
bilesen_sayisi = 10

## Hucreleri kumele-------------------------------------------------------------
meme_kanseri <- FindNeighbors(meme_kanseri, dims = 1:bilesen_sayisi)
meme_kanseri <- FindClusters(meme_kanseri, resolution = 0.5)

## Dogrusal olmayan boyut indirgemesi-------------------------------------------
meme_kanseri <- RunUMAP(meme_kanseri, dims = 1:bilesen_sayisi)

# indirgenmis kartezyen duzleminde verinin goruntulenmesi:
DimPlot(meme_kanseri, reduction = "umap")

# renkkorlugu dostu renk paleti:
renk_paleti <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
# indirgenmis kartezyen duzleminde verinin goruntulenmesi:
DimPlot(meme_kanseri, reduction = "umap", label = TRUE, pt.size = 0.2, cols = renk_paleti) + 
  theme(aspect.ratio = 1) + 
  NoLegend() +
  NoAxes()

## Diferansiyel ifade genler----------------------------------------------------
meme_kanseri_belirtecleri <- FindAllMarkers(meme_kanseri,
                                            only.pos = TRUE,
                                            min.pct = 0.25)

# her kume basina, ifadesi 1 log2 kat degisiminden buyuk olan genleri al:
meme_kanseri_belirtecleri %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

meme_kanseri_belirtecleri %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 1) -> genler

# keman grafigi:
VlnPlot(meme_kanseri, features = genler$gene)
# bilinen belirteclerle keman grafigi:
VlnPlot(meme_kanseri, features = c("PIP", "CDH1", "MUC1", "PECAM1", "MKI67")) & 
  NoLegend() & 
  xlab("") & 
  theme(aspect.ratio = 0.8)

# isi haritasi:
meme_kanseri_belirtecleri %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> genler

meme_kanseri <- ScaleData(meme_kanseri, features = genler$gene)
DoHeatmap(meme_kanseri, features = genler$gene, size = 2) + 
  NoLegend() + 
  scale_fill_gradient2()

## Hucre tipi isimlendirilmesi--------------------------------------------------
# https://www.immunesinglecell.org

# genlerin ortalama ifadesi:
meme_kanseri_ortalama <- AggregateExpression(meme_kanseri)
meme_kanseri_ortalama <- round(meme_kanseri_ortalama$RNA, 2)

write.table(meme_kanseri_ortalama, 
            "CELLiD_input.txt", 
            quote = F, 
            col.names = F, 
            row.names = T, 
            sep="\t")

# hucre tipi isimlendirme:
Idents(meme_kanseri) %>% unique()

meme_kanseri <- RenameIdents(meme_kanseri,
                             "0" = "Meme luminal hucresi",
                             "1" = "Meme luminal hucresi",
                             "2" = "Meme luminal hucresi",
                             "3" = "Meme luminal hucresi",
                             "4" = "Makrofaj",
                             "5" = "Fibroblast",
                             "6" = "Bolunen meme luminal oncul hucresi",
                             "7" = "DC")

meme_kanseri@meta.data$hucre_tipleri <- Idents(meme_kanseri)
head(meme_kanseri@meta.data)

# UMAP grafiginde genler:
genler <- c("MYC", "PIP", "LYZ", "PECAM1", "MKI67", "TOP2A")

FeaturePlot(meme_kanseri, features = genler, pt.size = 0.1, ncol = 3) & 
  NoAxes() & 
  theme(aspect.ratio = 1) &
  NoLegend()

DimPlot(meme_kanseri, group.by = "hucre_tipleri", pt.size = 0.1) + 
  theme(aspect.ratio = 1, legend.position = "right") +
  ggtitle("Hücre tipleri") +
  NoAxes()

# Diger analizler---------------------------------------------------------------
# Baska sefere 😊

sessionInfo()