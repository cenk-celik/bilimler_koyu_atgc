library(tidyverse)
# ggplot2 paketi ile gelen "tibble" verisini "elmas" olarak tanımlayalım
elmas <- diamonds
# veri çerçevesine dönüșüm:
elmas <- as.data.frame(elmas)

head(elmas)
tail(elmas)

# ilk üç satır ve diğer tanımlanmış bilgiler:
elmas[1:3, c("carat", "cut", "color", "clarity", "depth", "table", "price")]
# bir diğer yöntem:
elmas[1:3, !(colnames(elmas) %in% c("x", "y", "z"))]

# alt kumeleyelim:
elmas <- elmas[, !(colnames(elmas) %in% c("x", "y", "z"))]
head(elmas)

# yalnızca "carat", "clarity" ve "price" bilgileri:
elmas[1:3, c("carat", "clarity", "price")]

# sütun isimlerini "Turkcelestirelim":
colnames(elmas) <- c("karat", "kesim", "renk", "berraklik", "derinlik", "tablo", "fiyat")
head(elmas)

# "kesim" verilerini Turkcelestirelim:
unique(levels(elmas$kesim))
# "kesim" veri tipi
typeof(elmas$kesim)
elmas$kesim <- as.character(elmas$kesim)
typeof(elmas$kesim)
elmas$kesim[elmas$kesim == "Ideal"] <- "Cok iyi"
elmas$kesim[elmas$kesim == "Premium"] <- "Iyi"
elmas$kesim[elmas$kesim == "Very Good"] <- "Orta"
elmas$kesim[elmas$kesim == "Good"] <- "Kotu"
elmas$kesim[elmas$kesim == "Fair"] <- "Cok Kotu"

head(elmas)

# kalite fark etmeksizin ortalama fiyat
mean(elmas$fiyat)
median(elmas$fiyat)

# histogram
hist(elmas$fiyat)
hist(elmas$fiyat, breaks = 50)

# aynı histogramı güzelleştirelim:
ggplot(elmas, aes(x = fiyat)) +
  geom_histogram()

ggplot(elmas, aes(x = fiyat)) +
  geom_histogram() +
  theme_classic()

ggplot(elmas, aes(x = fiyat)) +
  geom_histogram(bins = 50, colour = "grey", fill = "whitesmoke", alpha = 0.5) +
  theme_classic() +
  theme(aspect.ratio = 1) + 
  xlab("Fiyat") +
  ylab("Adet") +
  ggtitle("Elmas fiyat dağılımı") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = mean(elmas$fiyat), linetype = "dashed", colour = "salmon") +
  geom_vline(xintercept = median(elmas$fiyat), linetype = "dashed", colour = "skyblue")

dir.create("veri_manipulasyonu_figurler")
setwd("veri_manipulasyonu_figurler")

pdf("elmas_fiyat_adet_histogram.pdf", width = 3, height = 3)
ggplot(elmas, aes(x = fiyat)) +
  geom_histogram(bins = 50, colour = "grey", fill = "whitesmoke", alpha = 0.5) +
  theme_classic() +
  theme(aspect.ratio = 1) + 
  xlab("Fiyat") +
  ylab("Adet") +
  ggtitle("Elmas fiyat dagilimi") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = mean(elmas$fiyat), linetype = "dashed", colour = "salmon") +
  geom_vline(xintercept = median(elmas$fiyat), linetype = "dashed", colour = "skyblue")
dev.off()

setwd("..")

# istatistiksel karsilastirma
head(elmas)

# kotu kesim ile iyi kesim arasinda fiyat fark var midir?
elmas_kotu_ve_iyi_kesim <- elmas[elmas$kesim == "Kotu" | elmas$kesim == "Iyi", ]

library(ggpubr)

setwd("veri_manipulasyonu_figurler")

pdf("elmas_kotu_ve_iyi_kesim_karsilastirma.pdf", width = 3, height = 3)
ggplot(elmas_kotu_ve_iyi_kesim, aes(x = kesim, y = fiyat, fill = kesim)) +
  geom_violin() +
  scale_fill_manual(values = c("Kotu" = "royalblue", "Iyi" = "red4")) +
  geom_jitter(alpha = 0.1, size = 0.1) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  xlab("Fiyat") +
  ylab("Adet") +
  ggtitle("Elmas fiyat karsilastirma") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", aspect.ratio = 1) +
  stat_compare_means(method = "t.test")
dev.off()

setwd("..")

# egzersiz: berrakligi "SI2" olan "Orta" elmaslar ile "VVS2" "Cok iyi" elmaslar
# arasinda fiyat farki var mi?
elmas$berraklik_kesim <- paste0(elmas$berraklik, "_", elmas$kesim)
head(elmas)

elmas_SI2_Orta_VVS2_Cok_iyi <- elmas[elmas$berraklik_kesim == "SI2_Orta" | elmas$berraklik_kesim == "VVS2_Cok iyi", ]

setwd("veri_manipulasyonu_figurler")

pdf("elmas_SI2_orta_ve_VVS2_cok_iyi_karsilastirma.pdf", width = 3, height = 3)
ggplot(elmas_SI2_Orta_VVS2_Cok_iyi, aes(x = berraklik_kesim, y = fiyat, fill = berraklik_kesim)) +
  geom_violin() +
  scale_fill_manual(values = c("SI2_Orta" = "royalblue", "VVS2_Cok iyi" = "red4")) +
  geom_jitter(alpha = 0.1, size = 0.1) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  xlab("Fiyat") +
  ylab("Adet") +
  ggtitle("Elmas fiyat karsilastirma") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", aspect.ratio = 1) +
  stat_compare_means(method = "t.test")
dev.off()

setwd("..")

# Gen ekpresyonu veri manipulasyonu---------------------------------------------
gen_ifadesi_verisi <- diff_express
head(gen_ifadesi_verisi)

# gereksiz sutunlardan kurtul:
gen_ifadesi_verisi$detection_call <- NULL

# sutun isimlerini turkcelestir:
colnames(gen_ifadesi_verisi) <- c("sembol", "ortalama", "log2KatDegisimi", "pduz")
head(gen_ifadesi_verisi)
gen_ifadesi_verisi <- na.omit(gen_ifadesi_verisi)

# azalan ya da artan genleri gostermek icin yanardag grafigi:
genler <- gen_ifadesi_verisi$sembol
log2_kat_degisimi <- gen_ifadesi_verisi$log2KatDegisimi
duzeltilmis_p <- gen_ifadesi_verisi$pduz

# yanlis kesif orani:
log10_duzeltilmis_p <- -log10(duzeltilmis_p)

yanardag_grafigi_verisi <- data.frame(gen = genler, log2KD = log2_kat_degisimi, log10pDuz = log10_duzeltilmis_p)
head(yanardag_grafigi_verisi)

# renk verisi:
kat_degisimi_limit = 1
yanlis_kesif_orani_limit = -log(0.05)
yanardag_grafigi_verisi$renk <- "onemsiz"  # ns for not significant
yanardag_grafigi_verisi$renk[yanardag_grafigi_verisi$log2KD > kat_degisimi_limit & yanardag_grafigi_verisi$log10pDuz > yanlis_kesif_orani_limit] <- "artan"
yanardag_grafigi_verisi$renk[yanardag_grafigi_verisi$log2KD < -kat_degisimi_limit & yanardag_grafigi_verisi$log10pDuz > yanlis_kesif_orani_limit] <- "azalan"

renkler <- c("onemsiz" = "grey", "artan" = "red3", "azalan" = "royalblue")

# en cok degisim gosteren 50 gen:
degisen_50_gen <- yanardag_grafigi_verisi[order(yanardag_grafigi_verisi$log10pDuz, decreasing = TRUE),][1:50, ]

library(ggrepel)
# hepsini veri cercevesine topla:
ggplot(yanardag_grafigi_verisi, aes(x = log2KD, y = log10pDuz, color = renk)) +
  geom_point(alpha = 0.2, size = 1) +
  scale_color_manual(values = renkler) +
  geom_hline(yintercept = 1.35, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  ggtitle("Yanardag grafigi") +
  xlab("log2 kat degisimi") +
  ylab("YKO") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1, legend.position = "bottom") +
  geom_text_repel(data = degisen_50_gen, aes(label = gen), size = 3, max.overlaps = 10)
