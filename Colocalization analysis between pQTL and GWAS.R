# Load packages
library(gwasglue)
library(coloc)
library(locuscomparer)
library(data.table)
library(openxlsx)
library(tidyverse)

# Read and process GWAS data
data1 <- fread('pcad_gwas.txt')  
colnames(data1)

data1$ncase.outcome <- 8496 
data1$ncontrol.outcome <- 26111
data1$samplesize.outcome <- 34607

data1 <- data1 %>% dplyr::select("SNP", "CHR", "BP", 'EA', 'NEA',
                                 "BETA", "SE", "P",
                                 "eaf", "N")
colnames(data1) <- c('SNP', 'chrom', "pos", 'effect_allele', 'other_allele', 
                     "beta", "se", "P", "eaf", "samplesize")

data1$number_cases <- 8496
data1$varbeta <- data1$se^2
data1$MAF <- ifelse(data1$eaf < 0.5, data1$eaf, 1 - data1$eaf)
data1 <- subset(data1, !duplicated(SNP))
data1$s <- data1$number_cases / data1$samplesize
data1$z <- data1$beta / data1$se

# Function for coloc analysis
perform_coloc <- function(pqtl_file, output_prefix, protein_name) {
  # Read pQTL data
  data2 <- fread(pqtl_file)
  
  data2$end <- data2$Pos
  data2 <- data2 %>% dplyr::select("rsids", "Chrom", "Pos", "end",
                                   "effectAllele", "otherAllele",
                                   "ImpMAF", "Beta", "SE", "Pval", "N")
  
  colnames(data2) <- c('SNP', 'chrom', "start", "end", 'effect_allele', 'other_allele',
                       "eaf", "beta", "se", "P", "samplesize")
  
  data2$varbeta <- data2$se^2
  data2$MAF <- ifelse(data2$eaf < 0.5, data2$eaf, 1 - data2$eaf)
  data2$z <- data2$beta / data2$se
  
  lead <- data2 %>% dplyr::arrange(P)
  leadchr <- lead$chrom[1]
  leadstart <- as.numeric(lead$start[1])
  leadend <- as.numeric(lead$end[1])
  
  QTLdata <- data2[data2$chrom == leadchr, ]
  QTLdata <- QTLdata[QTLdata$start > leadstart - 1000000 & QTLdata$end < leadend + 1000000, ]
  QTLdata <- subset(QTLdata, !duplicated(SNP))
  QTLdata <- na.omit(QTLdata)
  
  common_snps <- intersect(data1$SNP, QTLdata$SNP)
  data1_common <- data1[data1$SNP %in% common_snps, ]
  QTLdata_common <- QTLdata[QTLdata$SNP %in% common_snps, ]
  
  result <- coloc.abf(
    dataset1 = list(
      pvalues = data1_common$P, 
      snp = data1_common$SNP, 
      type = "cc", 
      s = data1_common$s[1], 
      N = data1_common$samplesize[1]
    ),
    dataset2 = list(
      pvalues = QTLdata_common$P, 
      snp = QTLdata_common$SNP, 
      type = "quant", 
      N = QTLdata_common$samplesize[1]
    ), 
    MAF = QTLdata_common$MAF
  )
  
  need_result <- result$results %>% dplyr::arrange(desc(SNP.PP.H4))
  write.csv(result$summary, paste0(output_prefix, '_coloc_result.csv'))
  
  gwas_fn <- data1[, c('SNP', 'P')] %>% dplyr::rename(rsid = SNP, pval = P)
  pqtl_fn <- QTLdata[, c('SNP', 'P')] %>% dplyr::rename(rsid = SNP, pval = P)
  
  p <- print(locuscompare(in_fn1 = gwas_fn,
                          in_fn2 = pqtl_fn,
                          title1 = 'GWAS', 
                          title2 = 'pQTL'))
  
  ggsave(
    filename = paste0(protein_name, ".tiff"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  return(result)
}

# Perform coloc analysis for each protein
proteins <- list(
  list("2418_55_APOE_Apo_E.txt.gz", "APOE", "APOE"),
  list("11516_7_FABP1_FABPL.txt.gz", "FABP1", "FABP1"),
  list("15364_101_APOC1_Apo_C_I.txt.gz", "APOC1", "APOC1"),
  list("16057_6_IGF2R_IGF_II_receptor.txt.gz", "IGF2R", "IGF2R"),
  list("3396_54_REN_Renin.txt.gz", "REN", "REN"),
  list("2837_3_MET_Met.txt.gz", "MET", "MET"),
  list("2982_82_LGALS4_Galectin_4.txt.gz", "LGALS4", "LGALS4"),
  list("3216_2_PIGR_PIGR.txt.gz", "PIGR", "PIGR"),
  list("10391_1_ANGPTL3_ANGL3.txt.gz", "ANGPTL3", "ANGPTL3"),
  list("8832_55_BST2_BST_2.txt.gz", "BST2", "BST2"),
  list("2516_57_CCL21_6Ckine.txt.gz", "CCL21", "CCL21")
)

for (protein in proteins) {
  pqtl_file <- protein[[1]]
  output_prefix <- protein[[2]]
  protein_name <- protein[[3]]
  
  cat("Processing", protein_name, "...\n")
  tryCatch({
    perform_coloc(pqtl_file, output_prefix, protein_name)
  }, error = function(e) {
    cat("Error processing", protein_name, ":", e$message, "\n")
  })
}
```