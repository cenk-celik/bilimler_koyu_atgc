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

print("Butun gerekli islemler tamamlandi! Kursta gorusmek uzere!")
print("Cenk")
