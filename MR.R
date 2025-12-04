library(TwoSampleMR)
library(data.table)
library(openxlsx)
library(dplyr)

# Fatty Acid Binding Protein 1
exp_dat1 = extract_instruments(outcomes = 'prot-a-1011', p1 = 5e-06, clump = TRUE, r2 = 0.001, kb = 10000,
                               force_server = FALSE)

# Calculate F-statistics
MR_F <- function(sample_size, num_IVs, r_square){
  numberator <- r_square * (sample_size - 1 - num_IVs)
  denominator <- (1 - r_square) * num_IVs
  f <- numberator / denominator
  return(f)
}

my_PVE <- function(beta, se_beta, maf, N){
  numberator <- 2 * beta^2
  denominator <- 2 * beta^2 + se_beta^2 * 2 * N
  return(numberator/denominator)
}

# Calculate R^2
exp_dat1$pve <- my_PVE(beta = exp_dat1$beta.exposure,
                       se_beta = exp_dat1$se.exposure,
                       N = exp_dat1$samplesize.exposure)

IV_infor <- exp_dat1 %>%
  group_by(SNP) %>%
  summarise(samplesize = mean(samplesize.exposure),
            nIV = n(),
            R2 = sum(pve))

IV_infor$F_stat <- MR_F(sample_size = IV_infor$samplesize,
                        num_IVs = IV_infor$nIV,
                        r_square = IV_infor$R2)

write.csv(IV_infor, file = "dat_Discoidin, CUB and LCCL Domain Containing 2_PCAD_F.csv")

# Extract outcome data PCAD
out_dat1 <- fread("pcad_mendelian.txt")
colnames(out_dat1) = c("SNP", "chr", "pos", "effect_allele", "other_allele",
                       "beta", "se", "pval", "samplesize")

out_dat1$outcome <- "PCAD"
colnames(out_dat1) = paste0(colnames(out_dat1), '.outcome')
colnames(out_dat1) = c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome",
                       "other_allele.outcome", "beta.outcome", "se.outcome",
                       "pval.outcome", "samplesize.outcome", "outcome")
out_dat1$id.outcome = "pcad"

# Harmonize data
dat1 = harmonise_data(exp_dat1, out_dat1, action = 3)

# MR analysis
res <- mr(dat1)
generate_odds_ratios(res)

# Heterogeneity analysis
mr_heterogeneity(dat1)

# Pleiotropy test
mr_pleiotropy_test(dat1)

# Single SNP analysis
res_single <- mr_singlesnp(dat1)
res_loo <- mr_leaveoneout(dat1)

# Scatter plot
p1 <- mr_scatter_plot(res, dat1)

tiff("Fatty Acid Binding Protein 1_PCAD_plot.tiff",
     width = 10,
     height = 8,
     units = "in",
     res = 300,
     compression = "lzw")

p1[[1]]
dev.off()

p2 <- mr_forest_plot(res_single)

tiff("Fatty Acid Binding Protein 1_PCAD_forest_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")

p2[[1]]
dev.off()

p3 <- mr_leaveoneout_plot(res_loo)

tiff("Fatty Acid Binding Protein 1_PCAD_leaveoneout_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")

p3[[1]]
dev.off()

res_single <- mr_singlesnp(dat1)
p4 <- mr_funnel_plot(res_single)

tiff("Fatty Acid Binding Protein 1_PCAD_funnel_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")

p4[[1]]
dev.off()

#Read exposure data: Average diameter for LDL particles
exp_dat1 = extract_instruments(outcomes = 'met-d-LDL_size', p1 = 5e-05, clump = TRUE, r2 = 0.001, kb = 10000,
                               force_server = FALSE)

# Calculate F-statistics
MR_F <- function(sample_size, num_IVs, r_square){
  numerator <- r_square * (sample_size - 1 - num_IVs)
  denominator <- (1 - r_square) * num_IVs
  f <- numerator / denominator
  return(f)
}

my_PVE <- function(beta, se_beta, maf, N){
  numerator <- 2 * beta^2
  denominator <- 2 * beta^2 + se_beta^2 * 2 * N
  return(numerator/denominator)
}

# Calculate R^2
exp_dat1$pve <- my_PVE(beta = exp_dat1$beta.exposure,
                       se_beta = exp_dat1$se.exposure,
                       N = exp_dat1$samplesize.exposure)

# Calculate F-statistics
library(dplyr)
IV_infor <- exp_dat1 %>%
  group_by(SNP) %>%
  summarise(samplesize = mean(samplesize.exposure),
            nIV = n(),
            R2 = sum(pve))

IV_infor$F_stat <- MR_F(sample_size = IV_infor$samplesize,
                        num_IVs = IV_infor$nIV,
                        r_square = IV_infor$R2)

write.csv(IV_infor, file = "Average diameter for LDL particles_PCAD_F.csv")

# Harmonize data
dat1 = harmonise_data(exp_dat1, out_dat1, action = 3)

# MR analysis
res <- mr(dat1)
generate_odds_ratios(res)

# Heterogeneity analysis
mr_heterogeneity(dat1)

# Pleiotropy test
mr_pleiotropy_test(dat1)

# Single SNP analysis
res_single <- mr_singlesnp(dat1)
res_loo <- mr_leaveoneout(dat1)

# Scatter plot
p1 <- mr_scatter_plot(res, dat1)
tiff("Average diameter for LDL particles_PCAD_scatter_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p1[[1]]
dev.off()

# Forest plot
p2 <- mr_forest_plot(res_single)
tiff("Average diameter for LDL particles_PCAD_forest_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p2[[1]]
dev.off()

# Leave-one-out plot
p3 <- mr_leaveoneout_plot(res_loo)
tiff("Average diameter for LDL particles_PCAD_leaveoneout_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p3[[1]]
dev.off()

# Funnel plot
res_single <- mr_singlesnp(dat1)
p4 <- mr_funnel_plot(res_single)
tiff("Average diameter for LDL particles_PCAD_funnel_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p4[[1]]
dev.off()
```r
# Read exposure data: HPV18
exp_dat1 = extract_instruments(outcomes = 'prot-c-2624_31_2', p1 = 5e-05, clump = TRUE, r2 = 0.001, kb = 10000,
                               force_server = FALSE)

# Calculate F-statistics
MR_F <- function(sample_size, num_IVs, r_square){
  numerator <- r_square * (sample_size - 1 - num_IVs)
  denominator <- (1 - r_square) * num_IVs
  f <- numerator / denominator
  return(f)
}

my_PVE <- function(beta, se_beta, maf, N){
  numerator <- 2 * beta^2
  denominator <- 2 * beta^2 + se_beta^2 * 2 * N
  return(numerator/denominator)
}

# Calculate R^2
exp_dat1$pve <- my_PVE(beta = exp_dat1$beta.exposure,
                       se_beta = exp_dat1$se.exposure,
                       N = exp_dat1$samplesize.exposure)

# Calculate F-statistics
library(dplyr)
IV_infor <- exp_dat1 %>%
  group_by(SNP) %>%
  summarise(samplesize = mean(samplesize.exposure),
            nIV = n(),
            R2 = sum(pve))

IV_infor$F_stat <- MR_F(sample_size = IV_infor$samplesize,
                        num_IVs = IV_infor$nIV,
                        r_square = IV_infor$R2)

# Extract outcome data: Phospholipids to total lipids ratio in chylomicrons and extremely large VLDL
out_dat1 <- extract_outcome_data(
  snps = exp_dat1$SNP,
  outcomes = "ebi-a-GCST90093049"
)

# Harmonize data
dat1 = harmonise_data(exp_dat1, out_dat1, action = 2)

# MR analysis
res <- mr(dat1)
generate_odds_ratios(res)

# Heterogeneity analysis
mr_heterogeneity(dat1)

# Pleiotropy test
mr_pleiotropy_test(dat1)

# Single SNP analysis
res_single <- mr_singlesnp(dat1)
res_loo <- mr_leaveoneout(dat1)

# Scatter plot
p1 <- mr_scatter_plot(res, dat1)
tiff("HPV18_Phospholipids to total lipids ratio in chylomicrons and extremely large VLDL_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p1[[1]]
dev.off()

# Forest plot
p2 <- mr_forest_plot(res_single)
tiff("HPV18_Phospholipids to total lipids ratio in chylomicrons and extremely large VLDL_forest_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p2[[1]]
dev.off()

# Leave-one-out plot
p3 <- mr_leaveoneout_plot(res_loo)
tiff("HPV18_Phospholipids to total lipids ratio in chylomicrons and extremely large VLDL_leaveoneout_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p3[[1]]
dev.off()

# Funnel plot
res_single <- mr_singlesnp(dat1)
p4 <- mr_funnel_plot(res_single)
tiff("HPV18_Phospholipids to total lipids ratio in chylomicrons and extremely large VLDL_funnel_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p4[[1]]
dev.off()
```r
# Read exposure data: HPV16
exp_dat1 = extract_instruments(outcomes = 'prot-c-2623_54_4', p1 = 5e-05, clump = TRUE, r2 = 0.001, kb = 10000,
                               force_server = FALSE)

# Calculate F-statistics
MR_F <- function(sample_size, num_IVs, r_square){
  numerator <- r_square * (sample_size - 1 - num_IVs)
  denominator <- (1 - r_square) * num_IVs
  f <- numerator / denominator
  return(f)
}

my_PVE <- function(beta, se_beta, maf, N){
  numerator <- 2 * beta^2
  denominator <- 2 * beta^2 + se_beta^2 * 2 * N
  return(numerator/denominator)
}

# Calculate R^2
exp_dat1$pve <- my_PVE(beta = exp_dat1$beta.exposure,
                       se_beta = exp_dat1$se.exposure,
                       N = exp_dat1$samplesize.exposure)

# Calculate F-statistics
library(dplyr)
IV_infor <- exp_dat1 %>%
  group_by(SNP) %>%
  summarise(samplesize = mean(samplesize.exposure),
            nIV = n(),
            R2 = sum(pve))

IV_infor$F_stat <- MR_F(sample_size = IV_infor$samplesize,
                        num_IVs = IV_infor$nIV,
                        r_square = IV_infor$R2)

# Extract outcome data: Cholesteryl Esters to Total Lipids in Very Small VLDL percentage
out_dat1 <- extract_outcome_data(
  snps = exp_dat1$SNP,
  outcomes = "met-d-M_VLDL_CE_pct"
)

# Harmonize data
dat1 = harmonise_data(exp_dat1, out_dat1, action = 3)

# MR analysis
res <- mr(dat1)
generate_odds_ratios(res)

# Heterogeneity analysis
mr_heterogeneity(dat1)

# Pleiotropy test
mr_pleiotropy_test(dat1)

# Single SNP analysis
res_single <- mr_singlesnp(dat1)
res_loo <- mr_leaveoneout(dat1)

# Scatter plot
p1 <- mr_scatter_plot(res, dat1)
tiff("HPV16_Cholesteryl Esters to Total Lipids in Very Small VLDL percentage_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p1[[1]]
dev.off()

# Forest plot
p2 <- mr_forest_plot(res_single)
tiff("HPV16_Cholesteryl Esters to Total Lipids in Very Small VLDL percentage_forest_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p2[[1]]
dev.off()

# Leave-one-out plot
p3 <- mr_leaveoneout_plot(res_loo)
tiff("HPV16_Cholesteryl Esters to Total Lipids in Very Small VLDL percentage_leaveoneout_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p3[[1]]
dev.off()

# Funnel plot
res_single <- mr_singlesnp(dat1)
p4 <- mr_funnel_plot(res_single)
tiff("HPV16_Cholesteryl Esters to Total Lipids in Very Small VLDL percentage_funnel_plot.tiff", 
     width = 10, 
     height = 8, 
     units = "in", 
     res = 300,
     compression = "lzw")
p4[[1]]
dev.off()
