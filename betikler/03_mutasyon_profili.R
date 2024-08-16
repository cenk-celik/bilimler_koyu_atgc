rm(list = ls())

library(maftools)
library(dplyr)
library(RColorBrewer)

proje_klasoru <- "~/bilimler_koyu_2024"
setwd(proje_klasoru)

# mutasyon verisini getir:
mutasyon_verisi <- read.csv("brca_tcga_pan_can_atlas_2018/data_mutations.txt", 
                            header = T, 
                            sep = "\t")
mutasyon_verisi[1:5, 1:7]

# klinik verisini getir:
klinik_verisi <- read.csv("brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", 
                          skip = 4, 
                          sep = "\t")

klinik_verisi[1:5, 1:7]
colnames(klinik_verisi)

# Verileri yeniden formatla:
# Kaplan-Meier grafikleri cizerken kullanilacak.
klinik_verisi$OS_STATUS[klinik_verisi$OS_STATUS == "0:LIVING"] <- 0
klinik_verisi$OS_STATUS[klinik_verisi$OS_STATUS == "1:DECEASED"] <- 1

klinik_verisi$DFS_STATUS[klinik_verisi$DFS_STATUS == "0:DiseaseFree"] <- 0
klinik_verisi$DFS_STATUS[klinik_verisi$DFS_STATUS == "1:Recurred/Progressed"] <- 1

klinik_verisi$PFS_STATUS[klinik_verisi$PFS_STATUS == "0:CENSORED"] <- 0
klinik_verisi$PFS_STATUS[klinik_verisi$PFS_STATUS == "1:PROGRESSION"] <- 1

# maftools kutuphanesi, klinik veride hasta barkodlari icin ozel
# Tumor_Sample_Barcode sutununa ihtiyac duyar. Bu yuzden yeniden isimlendirebilir,
# ya da yeni bir sutun olusturabiliriz.
klinik_verisi$Tumor_Sample_Barcode <- klinik_verisi$PATIENT_ID

# mutasyon verisi objesi olustur:
# butun takip eden analizler (downstream) bu obje uzerinden yurutulecek
mutasyon <- read.maf(mutasyon_verisi, 
                     clinicalData = klinik_verisi, 
                     rmFlags = TRUE, 
                     isTCGA = TRUE)
mutasyon

# hasta mutasyon verisi ozeti
getSampleSummary(mutasyon)

# hasta gen verisi ozeti
getGeneSummary(mutasyon)

# hasta klinik verisi ozeti
getClinicalData(mutasyon)

# (istege bagli) diske kaydet
# write.mafSummary(maf = mutasyon, basename = "brca")

# Gorsellestirme----------------------------------------------------------------
# mutasyon verisi ozet grafikleri
plotmafSummary(maf = mutasyon,
               rmOutlier = TRUE,
               addStat = "median",
               dashboard = TRUE,
               titvRaw = FALSE)

# minimalist mutasyon cubuk grafigi
mafbarplot(mutasyon)

# onkografik (oncoplot/waterfall plot)
oncoplot(maf = mutasyon, top = 10)

# renkleri degistirmek icin yeni palet olusturma:
renk_paleti = RColorBrewer::brewer.pal(n = 8, name = "Paired")
renk_paleti

# renklere sutun isimleri ekle
names(renk_paleti) = c(
  "Frame_Shift_Del",
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Multi_Hit",
  "Frame_Shift_Ins",
  "In_Frame_Ins",
  "Splice_Site",
  "In_Frame_Del"
)

print(renk_paleti)
oncoplot(maf = mutasyon, colors = renk_paleti, top = 10)

# istege bagli siniflandirma verili onkografik
oncoplot(maf = mutasyon,
         colors = renk_paleti,
         top = 10,
         clinicalFeatures = "SUBTYPE")

# amino asit degisiklikleri ve lolipop grafigi
lollipopPlot(maf = mutasyon,
             gene = "GATA3",
             AACol = "HGVSp_Short",
             showMutationRate = TRUE)

## Transisyon ve transversiyonlar-----------------------------------------------
# Transisyonlar, iki halkalı pürinlerin (A-G) veya tek halkalı pirimidinlerin (C-T) 
# yer değiştirmeleridir: dolayısıyla benzer şekilli bazları içerirler.
# Transversiyonlar, pürin bazlarının pirimidin bazlarıyla yer değiştirmesidir ve
# bu nedenle tek halkalı ve iki halkalı yapıların değişimini içerirler.

transisyon_transversion = titv(maf = mutasyon, plot = FALSE, useSyn = TRUE)

#plot titv summary
plotTiTv(res = transisyon_transversion)

# Somatik etkilesimler----------------------------------------------------------
# Birbirini dislayan (mutually exclusive)
# Birlikte gorulen (co-occuring)
somaticInteractions(maf = mutasyon,
                    top = 10,
                    pvalue = c(0.05, 0.1),
                    fontSize = 1,
                    colPal = "RdBu",
                    leftMar = 10,
                    topMar = 10)

# Kanseri tetikleyen genler-----------------------------------------------------
onemli_genler <- oncodrive(maf = mutasyon,
                           AACol = "HGVSp_Short",
                           minMut = 5,
                           pvalMethod = "zscore")
head(onemli_genler)
dev.off()
plotOncodrive(res = onemli_genler, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

# Kaplan-Meier hayatta kalma egrisi:
mafSurvival(maf = mutasyon,
            genes = "GATA3",
            time = "OS_MONTHS",
            Status = "OS_STATUS",
            showConfInt = TRUE,
            col = c("royalblue", "maroon"),
            isTCGA = TRUE)

mafSurvival(maf = mutasyon, 
            genes = "GATA3",
            time = "DFS_MONTHS",
            Status = "DFS_STATUS",
            showConfInt = TRUE,
            col = c("royalblue", "maroon"),
            isTCGA = TRUE)

mafSurvival(maf = mutasyon, 
            genes = "GATA3",
            time = "PFS_MONTHS",
            Status = "PFS_STATUS",
            showConfInt = TRUE,
            col = c("royalblue", "maroon"),
            isTCGA = TRUE)

# Hangi gen ikilisi hayatta kalmada onemli rol oynuyor?
kanser_ilerlemesi <- survGroup(maf = mutasyon,
                               top = 20, 
                               geneSetSize = 2, 
                               time = "OS_MONTHS", 
                               Status = "OS_STATUS", verbose = FALSE)

print(kanser_ilerlemesi, n = 136)

genler <- c("PIK3CA", "CDH1")

mafSurvival(maf = mutasyon,
            genes = genler,
            time = "OS_MONTHS",
            Status = "OS_STATUS",
            showConfInt = TRUE,
            col = c("royalblue", "maroon"),
            isTCGA = TRUE)

mafSurvival(maf = mutasyon,
            genes = genler,
            time = "DFS_MONTHS",
            Status = "DFS_STATUS",
            showConfInt = TRUE,
            col = c("royalblue", "maroon"),
            isTCGA = TRUE)

# Ilac ile hedeflenebilir genler------------------------------------------------
drugInteractions(maf = mutasyon)
# Spesifik bir gen icin ilac olasiliklari
drugInteractions(genes = "CDH1", drugs = TRUE)


yolaklar <- pathways(maf = mutasyon, plotType = "bar")
plotPathways(maf = mutasyon, pathlist = yolaklar)

# APOBEC zenginlesmesi----------------------------------------------------------
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)

# APOBEC kaynaklı mutasyonlar, katı tümörlerde daha sık görülür ve genellikle 
# TCW motifinde meydana gelen C>T geçiş olaylarıyla ilişkilidir.
# Kısaca, bir örnekteki tüm C>T mutasyonları üzerinde TCW motifinde meydana 
# gelen C>T mutasyonlarının zenginleşmesi, mutasyona uğramış bazların 20 bp 
# çevresinde bulunan arka plan Sitozinler ve TCW"lerle karşılaştırılır.

trinukleotitmatrisi <- trinucleotideMatrix(maf = mutasyon,
                                           prefix = "chr",
                                           add = TRUE, 
                                           ref_genome = "BSgenome.Hsapiens.UCSC.hg19")

plotApobecDiff(tnm = trinukleotitmatrisi, maf = mutasyon, pVal = 0.2)

# imza (Signature) analizi------------------------------------------------------

library("NMF")
imzalar <- estimateSignatures(mat = trinukleotitmatrisi, nTry = 6)

plotCophenetic(res = imzalar)

imzalar = extractSignatures(mat = trinukleotitmatrisi, n = 6)

cosmic_veritabani_ile_karsilastir <- compareSignatures(nmfRes = imzalar, sig_db = "SBS")
pheatmap::pheatmap(mat = cosmic_veritabani_ile_karsilastir$cosine_similarities, 
                   cluster_rows = FALSE, 
                   main = "Dogrulanmis imzalara gore benzerlik")

maftools::plotSignatures(nmfRes = imzalar, title_size = 1.2, sig_db = "SBS")

sessionInfo()