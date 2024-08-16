rm(list = ls())
# gerekli kutuphaneleri cagir
library(tidyverse)
library(limma)
library(edgeR)
library(caret)
library(glmnet)

## MÖ ile onkogen tahmini-------------------------------------------------------
load(file = "veri/limma_sonuclar.RData")
# transpoz ve matris donusumu:
matris = as.matrix(t(limma_sonuclar$voomObj$E))

# degiskeni faktor haline getir:
faktor = as.factor(limma_sonuclar$voomObj$targets$definition)

# veriyi test ve egitim seti olarak ikiye bol:
# tekrar edilebilirlik icin rastgele sayi uretimini sabitle:
set.seed(42)
egitim_seti = createDataPartition(faktor, p = 0.8, list = FALSE)

x_egitim = matris[egitim_seti, ]
x_test  = matris[-egitim_seti, ]

y_egitim = faktor[egitim_seti] 
y_test  = faktor[-egitim_seti]

# Elastic-Net modeli:
# Elastic Net modeli eğiteceğiz, bu model LASSO ve Ridge Regresyonu'nun en iyi 
# özelliklerini birleştiren genelleştirilmiş bir doğrusal modeldir.
# Ridge Regresyonu genellikle iyi tahminlerde bulunur, ancak sonuçları pek 
# yorumlanabilir değildir. LASSO, çok fazla gürültüden küçük sinyalleri almakta 
# iyidir, ancak fazlalığı minimize etme eğilimindedir, yani iki gen eşit derecede 
# iyi tahmincilerse (yüksek korelasyona sahip özellikler), genellikle sadece 
# birini seçer. Elastic Net, her iki yöntemin dengeli bir birleşimidir; 
# her bir durumu en iyi tahmin eden genleri veya gen gruplarını (eğer bunlar 
# korelasyona sahipse) seçer ve bu genleri sınıflandırma için model oluşturmakta
# kullanır. Daha sonra, bu genlere bireysel olarak bakarak, sınıflandırma problemi
# için biyolojik olarak önemli olup olmadığını görebiliriz. Elastic Net kullanırken, 
# özellikle alpha gibi ayarlamamız gereken başka parametreler de vardır. 
# Bu parametre, Elastic Net'in LASSO (alpha = 1) gibi mi yoksa Ridge Regresyonu 
# (alpha = 0) gibi mi davranacağını belirler. Basitlik adına alpha değerini 0.5 
# olarak ayarlayacağız, ancak gerçek bir senaryoda hatayı minimize etmek için bu 
# değeri değiştirerek en iyi modeli bulmaya çalışırdık.

# Capraz-dogrulama ile egitim setini kullanip modeli egit:
tahmin_sonuclari <- cv.glmnet(
  x = x_egitim,
  y = y_egitim,
  alpha = 0.5,
  family = "binomial" # yalnizca iki kategori var ise "binomial"; > 2 "multinomial"
)

# Model degerlendirmesi:
# test/tahmin
y_tahmin = predict(tahmin_sonuclari, newx = x_test, type = "class", s = "lambda.min")

# hata/karisiklik matrisi (confusion matrix)
hata_matrisi = table(y_tahmin, y_test)

# degerlendirme istatistikleri
print(hata_matrisi)

print(paste0("Sensitivity: ", sensitivity(hata_matrisi)))
print(paste0("Specificity: ", specificity(hata_matrisi)))
print(paste0("Precision: ",   precision(hata_matrisi)))

# tahminde kullanilan genlerin katsayilari:
katsayilar <- coef(tahmin_sonuclari, s = "lambda.min")
head(katsayilar)

# sifir olmayan katsayilari al:
katsayilar = katsayilar[katsayilar[,1] != 0,]
head(katsayilar)

# birinci satir intercept oldugundan bunu sil:
katsayilar = katsayilar[-1]

# modelde kullanilan genler:
ilgili_genler = names(katsayilar)
length(ilgili_genler)
head(ilgili_genler)

# limma ile ElasticNet karsilastirmasi:
print(intersect(rownames(limma_sonuclar$topGenes), ilgili_genler))

# Hiyerarsik kumeleme-----------------------------------------------------------
# Hem diferansiyel ifade genler hem de 
limma_ile_bulunan_genler = ifelse(
  (ilgili_genler %in% rownames(limma_sonuclar$topGenes)),
  "maroon",
  "lightblue"
)

# renk paleti tanimla:
renk_paleti <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)

# hiyerarsik kumeleme fonksiyonu tanimla:
kumeleme_fonksiyonu <- function(x) hclust(x, method = "complete")

# korelasyon matrisinin tersini uzaklik matrisi olarak kullan:
uzaklik_fonksiyonu <- function(x) as.dist((1-cor(t(x)))/2)

# isi haritasi:
isi_haritasi_matrisi <- t(matris[,ilgili_genler])
rownames(isi_haritasi_matrisi) <- ilgili_genler

gene_heatmap <- gplots::heatmap.2(
  isi_haritasi_matrisi,
  scale = "row",
  density.info = "none",
  trace = "none",
  col = renk_paleti,
  labRow = ilgili_genler,
  RowSideColors = limma_ile_bulunan_genler,
  labCol = FALSE,
  ColSideColors = as.character(as.numeric(faktor)),
  dendrogram = "both",
  hclust = kumeleme_fonksiyonu,
  distfun = uzaklik_fonksiyonu,
  cexRow = .6,
  margins = c(1,5)
)

# Gen Ontolojisi Analizi----------------------------------------------------------
# Hiyerarsik kumeleme bilgisini al:
hiyerarsik_kumeleme = as.hclust(gene_heatmap$rowDendrogram)

# dendrogrami, tumor ve normal olmak uzere ikiye bol:
kumeler <- cutree(hiyerarsik_kumeleme, k = 2)
table(kumeler)
# kanser Gen Ontolojisi analizi:
GO_kanser = enrichGO(names(kumeler[kumeler == 1]),
                     OrgDb = "org.Hs.eg.db", 
                     keyType = "SYMBOL", 
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1, 
                     ont = "BP",
                     pAdjustMethod = "BH")

head(GO_kanser@result)

# normal Gen Ontolojisi analizi:
GO_normal = enrichGO(names(kumeler[kumeler == 2]),
                     OrgDb = "org.Hs.eg.db", 
                     keyType = "SYMBOL", 
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1, 
                     ont = "BP",
                     pAdjustMethod = "BH")

head(GO_normal@result)
