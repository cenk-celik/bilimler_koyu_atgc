# Install the sendmailR package
if (!require("sendmailR", quietly = TRUE))
  install.packages("sendmailR")

# Function to send an email notification
send_error_email <- function(error_message) {
  library(sendmailR)
  
  from <- "<your_email@example.com>"
  to <- "<cenk.celik@ucl.ac.uk>"
  subject <- "Kutuphane yukleme hatasi"
  body <- paste("Kutuphane yuklenirken bir hata olustu:", error_message)
  
  sendmail(from, to, subject, body)
}

# Function to install packages with error handling
safe_install <- function(pkg) {
  tryCatch(
    {
      install.packages(pkg, ask = FALSE)
    },
    error = function(e) {
      send_error_email(e$message)
      stop(e)  # rethrow the error after sending the email
    }
  )
}

# Function to install Bioconductor packages with error handling
safe_bioconductor_install <- function(pkg) {
  tryCatch(
    {
      BiocManager::install(pkg, ask = FALSE)
    },
    error = function(e) {
      send_error_email(e$message)
      stop(e)  # rethrow the error after sending the email
    }
  )
}

print("CRAN paketleri indiriliyor...")

# Install CRAN packages
cran_packages <- c("tidyverse", "ggpubr", "ggrepel", "RColorBrewer", "NMF", 
                   "pheatmap", "DGEobj.utils", "hdf5r", "survminer", "caret", 
                   "glmnet", "devtools", "Seurat")

for (pkg in cran_packages) {
  safe_install(pkg)
}

# Install Bioconductor packages
print("Bioconductor indiriliyor...")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19", ask = FALSE)

bioconductor_packages <- c("maftools", "BSgenome.Hsapiens.UCSC.hg19", "edgeR", 
                           "clusterProfiler", "org.Hs.eg.db", "MultiAssayExperiment", 
                           "ComplexHeatmap")

print("Bioconductor paketleri indiriliyor...")

for (pkg in bioconductor_packages) {
  safe_bioconductor_install(pkg)
}

if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("BioinformaticsFMRP/TCGAbiolinks", ask = FALSE)

# Install GitHub packages
print("Yardimci paketler Github'dan indiriliyor...")
tryCatch(
  {
    devtools::install_github('immunogenomics/presto')
  },
  error = function(e) {
    send_error_email(e$message)
    stop(e)  # rethrow the error after sending the email
  }
)

print("Butun paketlerin indirilmesi tamamlandi!")
