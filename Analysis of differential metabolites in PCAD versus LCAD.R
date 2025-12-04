```r
library(data.table)

data1 = fread("LCAD.csv")
data2 = fread("pcad_metabolomics_imputed_pmm.csv")
colnames(data2)[2] = "beta_Hydroxybutyrate"

data = merge(data1, data2, by = "eid")
data = na.omit(data)

clean_feature_names <- function(df) {
  colnames(df) <- gsub(" ", "_", colnames(df))
  colnames(df) <- gsub("-", "_", colnames(df))
  return(df)
}

data <- clean_feature_names(data)

id = which(data$LCAD == 1)
length(id)

str(data)
rm(data1)
rm(data2)

set.seed(123)
train_sub = sample(nrow(data), 7/10 * nrow(data))
Train_data = data[train_sub,]
rownames(Train_data) = Train_data[[1]]
Train_data = Train_data[,-1]

Test_data = data[-train_sub,]
rownames(Test_data) = Test_data[[1]]
Test_data = Test_data[,-1]

x_train = data.matrix(Train_data[,29:279])
x_test = data.matrix(Test_data[,29:279])

x_train = scale(x_train, center = TRUE, scale = TRUE)
x_test = scale(x_test, center = TRUE, scale = TRUE)

library(glmnet)
y_train <- Train_data$LCAD
y_test = Test_data$LCAD

cvfit = cv.glmnet(x = x_train, y = y_train,
                  family = 'binomial',
                  alpha = 1,
                  nfolds = 10)

coef_mat <- coef(cvfit, s = "lambda.1se")
coef_mat <- as.matrix(coef_mat)

lasso_features <- rownames(coef_mat)[which(as.numeric(coef_mat) != 0)]
lasso_features <- setdiff(lasso_features, c("(Intercept)", "Intercept"))
lasso_features

# Metabolite lists
LCAD_metabolites <- c(
  "beta_Hydroxybutyrate",
  "Acetoacetate",
  "Acetone",
  "Alanine",
  "Albumin",
  "Apolipoprotein_B_to_Apolipoprotein_A1_ratio",
  "Average_Diameter_for_HDL_Particles",
  "Average_Diameter_for_LDL_Particles",
  "Average_Diameter_for_VLDL_Particles",
  "Cholesterol_in_IDL",
  "Cholesterol_in_Very_Small_VLDL",
  "Cholesterol_to_Total_Lipids_in_Very_Large_VLDL_percentage",
  "Cholesteryl_Esters_in_Chylomicrons_and_Extremely_Large_VLDL",
  "Cholesteryl_Esters_in_IDL",
  "Cholesteryl_Esters_in_Large_HDL",
  "Cholesteryl_Esters_in_Very_Large_HDL",
  "Cholesteryl_Esters_in_Very_Small_VLDL",
  "Cholesteryl_Esters_to_Total_Lipids_in_Large_VLDL_percentage",
  "Cholesteryl_Esters_to_Total_Lipids_in_Medium_VLDL_percentage",
  "Cholesteryl_Esters_to_Total_Lipids_in_Small_LDL_percentage",
  "Cholesteryl_Esters_to_Total_Lipids_in_Small_VLDL_percentage",
  "Cholesteryl_Esters_to_Total_Lipids_in_Very_Large_HDL_percentage",
  "Cholesteryl_Esters_to_Total_Lipids_in_Very_Small_VLDL_percentage",
  "Concentration_of_Medium_LDL_Particles",
  "Concentration_of_Small_HDL_Particles",
  "Creatinine",
  "Degree_of_Unsaturation",
  "Docosahexaenoic_Acid",
  "Free_Cholesterol_in_IDL",
  "Free_Cholesterol_in_Large_LDL",
  "Free_Cholesterol_in_Very_Small_VLDL",
  "Free_Cholesterol_to_Total_Lipids_in_Chylomicrons_and_Extremely_Large_VLDL_percentage",
  "Free_Cholesterol_to_Total_Lipids_in_Large_HDL_percentage",
  "Free_Cholesterol_to_Total_Lipids_in_Large_VLDL_percentage",
  "Free_Cholesterol_to_Total_Lipids_in_Medium_VLDL_percentage",
  "Free_Cholesterol_to_Total_Lipids_in_Very_Small_VLDL_percentage",
  "Glucose",
  "Glucose_lactate",
  "Glutamine",
  "Glycine",
  "Glycoprotein_Acetyls",
  "Histidine",
  "Leucine",
  "Linoleic_Acid_to_Total_Fatty_Acids_percentage",
  "Monounsaturated_Fatty_Acids",
  "Omega_3_Fatty_Acids",
  "Omega_6_Fatty_Acids",
  "Phenylalanine",
  "Phospholipids_to_Total_Lipids_in_Chylomicrons_and_Extremely_Large_VLDL_percentage",
  "Phospholipids_to_Total_Lipids_in_Large_HDL_percentage",
  "Phospholipids_to_Total_Lipids_in_Large_VLDL_percentage",
  "Phospholipids_to_Total_Lipids_in_Medium_HDL_percentage",
  "Phospholipids_to_Total_Lipids_in_Medium_LDL_percentage",
  "Phospholipids_to_Total_Lipids_in_Very_Large_VLDL_percentage",
  "Pyruvate",
  "Saturated_Fatty_Acids_to_Total_Fatty_Acids_percentage",
  "Spectrometer_corrected_alanine",
  "Total_Cholines",
  "Total_Fatty_Acids",
  "Triglycerides_in_Large_HDL",
  "Triglycerides_in_Medium_VLDL",
  "Triglycerides_in_Small_LDL",
  "Triglycerides_in_Small_VLDL",
  "Triglycerides_in_Very_Small_VLDL",
  "Triglycerides_to_Total_Lipids_in_IDL_percentage",
  "Triglycerides_to_Total_Lipids_in_Large_HDL_percentage",
  "Triglycerides_to_Total_Lipids_in_Large_LDL_percentage",
  "Triglycerides_to_Total_Lipids_in_Medium_LDL_percentage",
  "Triglycerides_to_Total_Lipids_in_Small_LDL_percentage",
  "Triglycerides_to_Total_Lipids_in_Very_Small_VLDL_percentage",
  "Tyrosine",
  "Valine"
)

PCAD_metabolites = c(
  "beta_Hydroxybutyrate",
  "Acetate",
  "Acetoacetate",
  "Average_Diameter_for_HDL_Particles",
  "Average_Diameter_for_LDL_Particles",
  "Average_Diameter_for_VLDL_Particles",
  "Cholesterol_to_Total_Lipids_in_Very_Large_HDL_percentage",
  "Cholesteryl_Esters_in_Medium_VLDL",
  "Cholesteryl_Esters_in_Very_Large_VLDL",
  "Cholesteryl_Esters_to_Total_Lipids_in_IDL_percentage",
  "Cholesteryl_Esters_to_Total_Lipids_in_Small_LDL_percentage",
  "Cholesteryl_Esters_to_Total_Lipids_in_Very_Small_VLDL_percentage",
  "Citrate",
  "Creatinine",
  "Degree_of_Unsaturation",
  "Free_Cholesterol_to_Total_Lipids_in_Chylomicrons_and_Extremely_Large_VLDL_percentage",
  "Glucose_lactate",
  "Glutamine",
  "Glycine",
  "Glycoprotein_Acetyls",
  "Histidine",
  "Linoleic_Acid_to_Total_Fatty_Acids_percentage",
  "Omega_3_Fatty_Acids",
  "Omega_3_Fatty_Acids_to_Total_Fatty_Acids_percentage",
  "Phenylalanine",
  "Phospholipids_to_Total_Lipids_in_Chylomicrons_and_Extremely_Large_VLDL_percentage",
  "Phospholipids_to_Total_Lipids_in_Small_HDL_percentage",
  "Pyruvate",
  "Saturated_Fatty_Acids_to_Total_Fatty_Acids_percentage",
  "Spectrometer_corrected_alanine",
  "Triglycerides_to_Total_Lipids_in_Medium_LDL_percentage",
  "Triglycerides_to_Total_Lipids_in_Very_Large_HDL_percentage",
  "Triglycerides_to_Total_Lipids_in_Very_Small_VLDL_percentage",
  "Valine"
)

inter_list <- intersect(PCAD_metabolites, LCAD_metabolites)
write.table(inter_list, file = 'inter_metabolomics.tsv', sep = '\t')

library(ggplot2)
library(ggvenn)

venn_list <- list(
  PCAD_metabolites = PCAD_metabolites,
  LCAD_metabolites = LCAD_metabolites
)

p = ggvenn(
  venn_list,
  fill_color = c("#00AFBB", "#E7B800"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 5
)

ggsave("venn_plot.tiff", p, width = 10, height = 8, dpi = 300)

PCAD_metabolites_only = setdiff(PCAD_metabolites, inter_list)
write.table(PCAD_metabolites_only, file = 'PCAD_metabolites_only.tsv', sep = '\t')

LCAD_metabolites_only = setdiff(LCAD_metabolites, inter_list)
write.table(LCAD_metabolites_only, file = 'LCAD_metabolites_only.tsv', sep = '\t')
```