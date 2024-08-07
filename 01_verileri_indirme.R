library(curl)

proje_klasoru <- "~/bilimler_koyu_2024"
dir.create("~/bilimler_koyu_2024")
setwd(proje_klasoru)

dir.create("ham_veri")
curl_download("https://cbioportal-datahub.s3.amazonaws.com/brca_tcga_pan_can_atlas_2018.tar.gz",
              destfile = "ham_veri/brca_tcga_pan_can_atlas_2018.tar.gz")

untar("ham_veri/brca_tcga_pan_can_atlas_2018.tar.gz")
