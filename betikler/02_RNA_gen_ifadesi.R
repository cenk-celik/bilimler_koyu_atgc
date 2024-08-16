rm(list = ls())

# gerekli kutuphaneleri cagir
library(tidyverse)
library(limma)
library(edgeR)
library(DGEobj.utils)
library(survival)
library(survminer)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)

proje_klasoru <- "~/bilimler_koyu_2024"
setwd(proje_klasoru)

# RNA dizileme veri analizi-----------------------------------------------------
## RNA dizileme verisi manipulasyonu--------------------------------------------
# Veriyi ice aktar:
# Kanser
RNA_ekspresyonu_kanser <- read.csv("brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt", 
                                   header = T, 
                                   sep = "\t")

# Gereksiz sutunlari kaldir:
RNA_ekspresyonu_kanser$Entrez_Gene_Id <- NULL

# Normal
RNA_ekspresyonu_normal <- read.csv("brca_tcga_pan_can_atlas_2018/normals/data_mrna_seq_v2_rsem_normal_samples.txt", 
                                   header = T, 
                                   sep = "\t")

# Gereksiz sutunlari kaldir:
RNA_ekspresyonu_normal$Entrez_Gene_Id <- NULL

# Normal ve Kanser gen ifadesi veri cercevelerini birlestir:
RNA_ekspresyonu <- merge(RNA_ekspresyonu_kanser, RNA_ekspresyonu_normal, by = "Hugo_Symbol", all = TRUE)

dim(RNA_ekspresyonu)

# ayni genlerin ortalamasini al:
# RNA_ekspresyonu <- summarise_all(group_by(RNA_ekspresyonu, Hugo_Symbol))
RNA_ekspresyonu <- RNA_ekspresyonu %>% group_by(Hugo_Symbol) %>% summarise_all(mean)

# birinci satiri sil (bos satirlar):
RNA_ekspresyonu <- RNA_ekspresyonu[!(RNA_ekspresyonu$Hugo_Symbol == ""), ]

# veri cercevesine donustur:
RNA_ekspresyonu <- as.data.frame(RNA_ekspresyonu)

# Gen isimlerini satir isimlerine ekle:
rownames(RNA_ekspresyonu) <- RNA_ekspresyonu$Hugo_Symbol
RNA_ekspresyonu$Hugo_Symbol <- NULL

# NA degerlerini sifir ile degistir:
RNA_ekspresyonu[is.na(RNA_ekspresyonu)] <- 0

# matrise cevir:
RNA_ekspresyonu <- as.matrix(RNA_ekspresyonu)
RNA_ekspresyonu[1:5, 1:7]

# hasta ID manipulasyonu:
colnames(RNA_ekspresyonu) <- gsub(".01", "", colnames(RNA_ekspresyonu))
colnames(RNA_ekspresyonu) <- gsub(".11", "", colnames(RNA_ekspresyonu))
colnames(RNA_ekspresyonu) <- gsub("\\.", "-", colnames(RNA_ekspresyonu))

# klinik verisini ice aktar:
klinik_verisi <- read.csv("brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", 
                          skip = 4,
                          sep = "\t")

# herhangi bir kategori verisine erisim:
klinik_verisi$SUBTYPE %>% table()

# normal ID'ler:
colnames(RNA_ekspresyonu_normal) <- gsub(".11", "", colnames(RNA_ekspresyonu_normal))
colnames(RNA_ekspresyonu_normal) <- gsub("\\.", "-", colnames(RNA_ekspresyonu_normal))
normal_hastalar <- colnames(RNA_ekspresyonu_normal)

normal_hasta_veri_cercevesi <- data.frame(
  PATIENT_ID = normal_hastalar,
  SUBTYPE = "Normal",
  stringsAsFactors = FALSE
)

klinik_verisi <- bind_rows(klinik_verisi, normal_hasta_veri_cercevesi)
klinik_verisi$SUBTYPE %>% table()

# normal ve kanser kategorilerini ekle:
klinik_verisi$definition <- "Primary Solid Tumour"
klinik_verisi$definition[klinik_verisi$SUBTYPE == "Normal"] <- "Solid Tissue Normal"

klinik_verisi$definition %>% table()

# her iki veri setinde var olan hastalari kullan
klinik_verisi <- klinik_verisi[colnames(RNA_ekspresyonu) %in% klinik_verisi$PATIENT_ID, ]

# SummarisedExperiment nesnesi olustur
meme_kanseri_transkriptomik <- SummarizedExperiment(assays = SimpleList(counts = RNA_ekspresyonu),
                                                    colData = klinik_verisi)
# veri seti boyutu
dim(meme_kanseri_transkriptomik)

# mevcut klinik veri:
colnames(colData(meme_kanseri_transkriptomik))

# kategorik verilere erisim ve ozet tablosu:
table(meme_kanseri_transkriptomik@colData$OS_STATUS)
table(meme_kanseri_transkriptomik@colData$SEX)
table(meme_kanseri_transkriptomik@colData$SUBTYPE)

# ilk 5 hasta gen ekspresyonu verisine erisim:
head(assay(meme_kanseri_transkriptomik)[,1:5])

# veriyi diske kaydet:
dir.create("veri")
# saveRDS(meme_kanseri_transkriptomik, file = "veri/meme_kanseri_transkriptomik.rds")
# save(meme_kanseri_transkriptomik, file = "veri/meme_kanseri_transkriptomik.RData")

## RNA dizileme verisi normalizasyonu-------------------------------------------
# limma_pipeline fonksiyonunu tanımla
# Linear Models for Microarray Data
# Mikroarray ve RNA dizileme verisi icin veri analizi paketi
limma_akisi <- function(islenecek_veri, test_edilecek_degisken, referans=NULL) {
  
  # tcga_data'dan belirtilen test_edilecek_degisken kullanarak tasarım faktörünü çıkar
  tasarim_faktoru = colData(islenecek_veri)[, test_edilecek_degisken, drop=T]
  
  # Tasarım faktörünü bir faktöre dönüştür ve referans belirtilmişse yeniden seviyelendir
  grup = factor(tasarim_faktoru)
  # referans grubu tanimliysa referansi tanimla:
  if (!is.null(referans)) {
    grup = relevel(grup, ref=referans)
  }
  
  # Lineer model için tasarım matrisini oluştur
  tasarim = model.matrix(~ grup)
  
  # Sayımlar, örnek bilgileri ve gen bilgileri ile bir DGEList nesnesi oluştur
  diferansiyel_gen_ifadesi = DGEList(
    counts = assay(islenecek_veri),
    samples = colData(islenecek_veri),
    genes = as.data.frame(rowData(islenecek_veri))
  )
  
  # Düşük ifadelendirilmiş genleri filtrele
  yuksek_ifadelenmis_genler = filterByExpr(diferansiyel_gen_ifadesi, tasarim)
  diferansiyel_gen_ifadesi = diferansiyel_gen_ifadesi[yuksek_ifadelenmis_genler, , keep.lib.sizes=FALSE]
  rm(yuksek_ifadelenmis_genler)
  
  # Veriyi TMM normalizasyonu ve ardından voom dönüşümü kullanarak normalize et
  diferansiyel_gen_ifadesi = calcNormFactors(diferansiyel_gen_ifadesi)
  voom_donusumu = voom(diferansiyel_gen_ifadesi, tasarim, plot=TRUE)
  
  # Voom dönüşümlü veriye lineer bir model uygula
  fit = lmFit(voom_donusumu, tasarim)
  fit = eBayes(fit)
  
  # En fazla farklı şekilde ifade edilen 100 geni p-değerine göre sıralayarak çıkar
  diferansiyel_ifade_genler = topTable(fit, coef = ncol(tasarim), number = 100, sort.by = "p")
  
  # Sonuçları voom nesnesi, fit nesnesi ve en iyi genleri içeren bir liste olarak döndür
  return(list(
    voomObj = voom_donusumu, # Normalizasyon yapılmış veri
    fit = fit, # Uygulanan model ve istatistikler
    topGenes = diferansiyel_ifade_genler # En fazla farklı şekilde ifade edilen 100 gen
  ))
}

# limma akisi uygulanmasi:
limma_sonuclar <- limma_akisi(
  islenecek_veri = meme_kanseri_transkriptomik,
  test_edilecek_degisken = "definition",
  referans="Solid Tissue Normal"
)

# baska degiskenler de test edilebilir:
# limma_sonuclar <- limma_akisi(
#   islenecek_veri = meme_kanseri_transkriptomik,
#   test_edilecek_degisken = "definition",
#   referans="Solid Tissue Normal"
# )

save(limma_sonuclar, file = "veri/limma_sonuclar.RData")

## Gorsellestirme---------------------------------------------------------------
# Temel bilesenler analizi:
plot_PCA <- function(voom_nesnesi, test_edilecek_degisken){
  grup = factor(voom_nesnesi$targets[, test_edilecek_degisken])
  # temel bilesen analizi uygula:
  pca = prcomp(t(voom_nesnesi$E))
  # ilk iki temel bileseni grafik olarak goster:
  plot(pca$x[,1:2], col = grup, pch = 19)
  # lejant ekle:
  legend("topright", inset = .01, levels(grup), pch = 19, col = 1:length(levels(grup)))
  return(pca)
}

pca_sonuc <- plot_PCA(limma_sonuclar$voomObj, "definition")

# Diferansiyel genlerin isi haritasi:
# limma_sonuclar nesnesinden gen ekspresfonu matrisini al:
ekspresyon_matrisi <- as.matrix(t(limma_sonuclar$voomObj$E))

# 20 diferansiyel ifade genler:
head(limma_sonuclar$topGenes, 20)

# 20 en farkli ifade edilmis gen:
genler = rownames(limma_sonuclar$topGenes)[1:20]

# kategori bilgisini al:
kategoriler <- factor(limma_sonuclar$voomObj$targets$definition)

# renk paleti tanimla:
renk_paleti <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)

# hiyerarsik kumeleme fonksiyonu tanimla:
kumeleme_fonksiyonu <- function(x) hclust(x, method = "complete")

# korelasyon matrisinin tersini uzaklik matrisi olarak kullan:
uzaklik_fonksiyonu <- function(x) as.dist((1-cor(t(x)))/2)

# isi haritasi:
isi_haritasi <- gplots::heatmap.2(
  t(ekspresyon_matrisi[, genler]),
  scale = "row", # olceklendirme
  density.info = "none", # lejant yogunluk cizgisi
  trace = "none", # isi haritasi cizgileri
  col = renk_paleti, # renk paleti
  labCol = FALSE, # sutun isimlerini gosterme
  ColSideColors = as.character(as.numeric(kategoriler)), # Show colors for each response class
  dendrogram = "both", # dendrogram
  hclust = kumeleme_fonksiyonu, # hiyerarsik kumeleme
  distfun = uzaklik_fonksiyonu, # uzaklik korelasyon matrisi
  cexRow = 1.0, # gen isimleri punto
  keysize = 1.25, # lejant punto
  margins = c(1, 6) # mizanpaj belirle
)

## Hayatta kalma egrisi (Kaplan-Meier)------------------------------------------
klinik_verisi <- meme_kanseri_transkriptomik@colData

# istenen kategoriler:
kategoriler <- c("PATIENT_ID",
                 "definition",
                 "RADIATION_THERAPY",
                 "OS_STATUS",
                 "OS_MONTHS",
                 "DAYS_LAST_FOLLOWUP",
                 "SEX",
                 "AJCC_PATHOLOGIC_TUMOR_STAGE")

# Yalnizca kanser orneklerini al:
klinik_veri_cercevesi <- klinik_verisi[klinik_verisi$definition == "Primary Solid Tumour", kategoriler]

# ikili (boolean) degisken bilgisi olustur:
klinik_veri_cercevesi$deceased <- klinik_veri_cercevesi$OS_STATUS == "1:DECEASED"

# show first 10 samples
head(klinik_veri_cercevesi)

# hayatta kalma nesnesi olustur:
fit = survfit(Surv(OS_MONTHS, deceased) ~ RADIATION_THERAPY, data = klinik_veri_cercevesi)
print(fit)

# hayatta kalma grafigi:
ggsurvplot(fit, data = klinik_veri_cercevesi)

# estetik:
ggsurvplot(fit,
           data = klinik_veri_cercevesi,
           conf.int = FALSE,
           pval = TRUE,
           palette = "npg",
           risk.table = TRUE,
           ggtheme = theme_survminer() + theme(legend.position = "none"),
           surv.median.line = "hv")

sessionInfo()