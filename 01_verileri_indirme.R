library(curl)

proje_klasoru <- "~/bilimler_koyu_2024"
dir.create("~/bilimler_koyu_2024")
print(paste0("Proje klasoru olusturuldu: ", proje_klasoru))
setwd(proje_klasoru)
print(paste0("Calisma klasoru ", proje_klasoru, " olarak degistirildi."))

dir.create("ham_veri")
print(paste0("Ham veri klasoru olusturuldu: ", "ham_veri"))

print(paste0("Kanser Genom Atlasi meme kanseri verileri indiriliyor..."))

curl_download("https://cbioportal-datahub.s3.amazonaws.com/brca_tcga_pan_can_atlas_2018.tar.gz",
              destfile = "ham_veri/brca_tcga_pan_can_atlas_2018.tar.gz")
print(paste0("Veriler indirildi... Dosya aciliyor!"))
untar("ham_veri/brca_tcga_pan_can_atlas_2018.tar.gz")

print("Metilasyon verisine erisiliyor...")
library(TCGAbiolinks)
metilasyon <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "DNA Methylation", 
                       platform = "Illumina Human Methylation 27", 
                       data.type = "DNA Methylation")

print("Metilasyon verisi indiriliyor...")
GDCdownload(metilasyon)
metilasyon_verisi <- GDCprepare(query = metilasyon, summarizedExperiment = TRUE)
print("Metilasyon verisi kaydediliyor...")
save(metilasyon_verisi, file = "ham_veri/metilasyon.RData")
print("Metilasyon verisi kaydedildi!")

print("Butun gerekli islemler tamamlandi! Kursta gorusmek uzere!")
print("Cenk")
