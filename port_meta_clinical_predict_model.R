```r
# 1. Data preprocessing
# Load packages
library(data.table)

# Load data
PCAD = fread(file = 'PCAD.csv')
meta = fread(file = 'pcad_metabolomics_imputed_pmm.csv')
colnames(meta)[2] = 'beta-Hydroxybutyrate'
port = fread(file = 'proteomics_imputed_pmm.csv')

# Merge datasets
A = merge(PCAD, meta, by = 'eid')
B = merge(A, port, by = 'eid')

# Define protein and metabolite lists
proteins <- c("ren", "ntprobnp", "angptl3", "met", 
              "lgals4", "pigr", "ccl21", "gast", "bst2", "fabp1")

metabolites <- c(
  "Average_Diameter_for_HDL_Particles",
  "Cholesteryl_Esters_to_Total_Lipids_in_Small_LDL_percentage",
  "Cholesteryl_Esters_to_Total_Lipids_in_Very_Small_VLDL_percentage",
  "Glucose_lactate",
  "Linoleic_Acid_to_Total_Fatty_Acids_percentage",
  "Creatinine",
  "Phospholipids_to_Total_Lipids_in_Chylomicrons_and_Extremely_Large_VLDL_percentage",
  "Average_Diameter_for_VLDL_Particles",
  "Glycine",
  "Acetoacetate"
)

# Prepare data
B = as.data.frame(B)
colnames(B) = gsub(" ", "_", colnames(B))
colnames(B) = gsub("-", "_", colnames(B))

C1 = B[, c(2:29)]
C2 = B[, proteins]
C3 = B[, metabolites]
data_new = cbind(cbind(C1, C2), C3)

# Factor conversion
colnames(data_new)[11] = "Annual_household_income"
str(data_new)

data_new$PCAD = factor(data_new$PCAD)
data_new$Sex = factor(data_new$Sex)
data_new$Ethnicity = factor(data_new$Ethnicity)
data_new$Education = factor(data_new$Education)
data_new$Annual_household_income = factor(data_new$Annual_household_income)
data_new$Drinking = factor(data_new$Drinking)
data_new$Smoking = factor(data_new$Smoking)
data_new$Duration_of_moderate_activity = factor(data_new$Duration_of_moderate_activity)
data_new$Hypertension = factor(data_new$Hypertension)
data_new$Diabetes = factor(data_new$Diabetes)
data_new$Hyperlipidemia = factor(data_new$Hyperlipidemia)

# 2. Model building
# Load packages
library(rms)
library(tableone)
library(broom)
library(pROC)
library(regplot)
library(ResourceSelection)
library(nricens)
library(PredictABEL)
library(ggDCA)
library(caret)
library(survivalROC)

# Split data
set.seed(20251115)
trainIndex <- sample(nrow(data_new), 0.7 * nrow(data_new))
trainset <- data_new[trainIndex, ]
testset <- data_new[-trainIndex, ]

# Rename columns
colnames(data_new)[39:48] = c('AD_for_HDL',
                              'CE_to_TL_in_small_LDL',
                              'CE_to_TL_in_very_small_VLDL',
                              'Glucose_lactate',
                              'LA_to_TFA',
                              'Creatinine',
                              'Phospholipids_to_TL_in_CEL_VLDL',
                              'AD_for_VLDL',
                              'Glycine',
                              'Acetoacetate')

colnames(trainset)[39:48] = c('AD_for_HDL',
                              'CE_to_TL_in_small_LDL',
                              'CE_to_TL_in_very_small_VLDL',
                              'Glucose_lactate',
                              'LA_to_TFA',
                              'Creatinine',
                              'Phospholipids_to_TL_in_CEL_VLDL',
                              'AD_for_VLDL',
                              'Glycine',
                              'Acetoacetate')

colnames(testset)[39:48] = c('AD_for_HDL',
                             'CE_to_TL_in_small_LDL',
                             'CE_to_TL_in_very_small_VLDL',
                             'Glucose_lactate',
                             'LA_to_TFA',
                             'Creatinine',
                             'Phospholipids_to_TL_in_CEL_VLDL',
                             'AD_for_VLDL',
                             'Glycine',
                             'Acetoacetate')

# Set datadist
ddist <- datadist(trainset)
options(datadist = 'ddist')

# Build logistic regression model
model1 <- glm(PCAD ~ Age + Sex + Ethnicity + Smoking + BMI +
                Hyperlipidemia + Hypertension + Diabetes +
                angptl3 + gast + lgals4 + ntprobnp + pigr + ren + fabp1 + met +
                bst2 + ccl21 + AD_for_HDL + 
                CE_to_TL_in_small_LDL + CE_to_TL_in_very_small_VLDL + 
                Glucose_lactate + Creatinine + 
                Phospholipids_to_TL_in_CEL_VLDL + AD_for_VLDL + 
                Glycine + Acetoacetate,
              data = trainset, family = "binomial")

step_model1 <- step(model1, direction = "backward")
summary(step_model1)

# 3. Model visualization and evaluation in trainset
# ROC curve
trainset$prob <- predict(step_model1, newdata = trainset, type = "response")
gmodel1 <- roc(PCAD ~ prob, data = trainset, smooth = FALSE)

tiff(filename = 'ROC_trainset.tiff', 
     width = 5, height = 5, units = 'in', res = 300)

plot(gmodel1, print.auc = FALSE, main = "ROC Curve in Trainset", 
     col = "blue", print.thres.col = "blue", identity.col = "red", 
     identity.lty = 2, identity.lwd = 2, xlim = c(1, 0))

auc_obj <- roc(response = trainset$PCAD, predictor = trainset$prob)
ci_obj <- ci(auc_obj)

text(0.5, 0.6, paste0("AUC: ", round(auc_obj$auc, 3)), col = "black")
text(0.5, 0.5, paste0("95% CI: [", round(ci_obj[1], 3), ", ", 
                      round(ci_obj[3], 3), "]"), col = "black")

dev.off()

# Nomogram
trainset$PCAD <- as.numeric(trainset$PCAD) - 1
ddist <- datadist(trainset)
options(datadist = 'ddist')

nomo1 <- lrm(PCAD ~ Age + Sex + Ethnicity + Smoking + BMI +
               Hyperlipidemia + Hypertension + Diabetes +
               angptl3 + gast + lgals4 + ntprobnp + pigr + ren + fabp1 + met +
               bst2 + ccl21 + AD_for_HDL + 
               CE_to_TL_in_small_LDL + CE_to_TL_in_very_small_VLDL + 
               Glucose_lactate + Creatinine + 
               Phospholipids_to_TL_in_CEL_VLDL + AD_for_VLDL + 
               Glycine + Acetoacetate,
             data = trainset, x = TRUE, y = TRUE)

tiff(filename = 'Nomo for PCAD.tiff', 
     height = 15, width = 12, units = 'in', res = 300)

regplot(nomo1, plot = c('density', 'boxes'), center = TRUE,
        title = 'Nomogram to predict PCAD',
        odds = FALSE, showP = FALSE, observation = TRUE,
        points = TRUE, droplines = FALSE, interval = 'confidence',
        rank = NULL, dencol = 'skyblue', boxcol = 'red')

dev.off()

# Calibration plot
tiff('calibration_plot_trainset.tiff', width = 5, height = 5, units = 'in', res = 300)
val.prob(trainset$prob, trainset$PCAD, xlab = "Predicted Probability in Trainset", 
         ylab = "Actual Probability in Trainset")
dev.off()

# DCA curve
library(rmda)
DCA.model1 <- decision_curve(PCAD ~ Age + Sex + Ethnicity + Smoking + BMI +
                               Hyperlipidemia + Hypertension + Diabetes +
                               angptl3 + gast + lgals4 + ntprobnp + pigr + ren + fabp1 + met +
                               bst2 + ccl21 + AD_for_HDL + 
                               CE_to_TL_in_small_LDL + CE_to_TL_in_very_small_VLDL + 
                               Glucose_lactate + Creatinine + 
                               Phospholipids_to_TL_in_CEL_VLDL + AD_for_VLDL + 
                               Glycine + Acetoacetate, 
                             data = trainset, bootstraps = 100)

tiff('decision_plot_trainset.tiff', width = 5, height = 5, units = 'in', res = 300)
plot_decision_curve(DCA.model1,  
                    curve.names = "DCA Curves in Trainset",
                    cost.benefit.axis = FALSE,
                    col = c('red', 'blue', 'black'),
                    confidence.intervals = FALSE,
                    standardize = FALSE)
dev.off()

# 4. Model evaluation in testset
# ROC curve
testset$prob <- predict(step_model1, newdata = testset, type = "response")
gmodel2 <- roc(PCAD ~ prob, data = testset, smooth = FALSE)

tiff(filename = 'ROC_testset.tiff', 
     width = 5, height = 5, units = 'in', res = 300)

plot(gmodel2, print.auc = FALSE, main = "ROC Curve in Testset", 
     col = "blue", print.thres.col = "blue", identity.col = "red", 
     identity.lty = 2, identity.lwd = 2, xlim = c(1, 0))

auc_obj <- roc(response = testset$PCAD, predictor = testset$prob)
ci_obj <- ci(auc_obj)

text(0.5, 0.6, paste0("AUC: ", round(auc_obj$auc, 3)), col = "black")
text(0.5, 0.5, paste0("95% CI: [", round(ci_obj[1], 3), ", ", 
                      round(ci_obj[3], 3), "]"), col = "black")

dev.off()

# Calibration plot
testset$PCAD <- as.numeric(testset$PCAD) - 1
tiff('calibration_plot_testset.tiff', width = 5, height = 5, 
     units = 'in', res = 300)
val.prob(testset$prob, testset$PCAD)
dev.off()

# DCA curve
DCA.model2 <- decision_curve(PCAD ~ Age + Sex + Ethnicity + Smoking + BMI +
                               Hyperlipidemia + Hypertension + Diabetes +
                               angptl3 + gast + lgals4 + ntprobnp + pigr + ren + fabp1 + met +
                               bst2 + ccl21 + AD_for_HDL + 
                               CE_to_TL_in_small_LDL + CE_to_TL_in_very_small_VLDL + 
                               Glucose_lactate + Creatinine + 
                               Phospholipids_to_TL_in_CEL_VLDL + AD_for_VLDL + 
                               Glycine + Acetoacetate, 
                             data = testset, bootstraps = 100)

tiff('decision_plot_testset.tiff', width = 5, height = 5, units = 'in', res = 300)
plot_decision_curve(DCA.model2,  
                    curve.names = "DCA Curves in Testset",
                    cost.benefit.axis = FALSE,
                    col = c('red', 'blue', 'black'),
                    confidence.intervals = FALSE,
                    standardize = FALSE)
dev.off()

# 5. Pearson correlation analysis between metabolites and proteins
proteins <- c("ren", "ntprobnp", "angptl3", "met", "lgals4", "pigr", 
              "ccl21", "gast", "bst2", "fabp1")

metabolites <- c('AD_for_HDL',
                 'CE_to_TL_in_small_LDL',
                 'CE_to_TL_in_very_small_VLDL',
                 'Glucose_lactate',
                 'LA_to_TFA',
                 'Creatinine',
                 'Phospholipids_to_TL_in_CEL_VLDL',
                 'AD_for_VLDL',
                 'Glycine',
                 'Acetoacetate')

# Prepare correlation data
data3 <- data_new[, c(proteins, metabolites)]
data3 <- as.matrix(data3)

# Load correlation packages
if (!require(corrplot)) install.packages("corrplot")
if (!require(Hmisc)) install.packages("Hmisc")
if (!require(ggcorrplot)) install.packages("ggcorrplot")

library(corrplot)
library(Hmisc)
library(ggcorrplot)

# Calculate correlation matrix and p-values
cor_result <- rcorr(as.matrix(data3), type = "pearson")
cor_matrix <- cor_result$r
p_matrix <- cor_result$P

# Create heatmap
tiff("metabolomics_proteomics_heatmap.tiff", 
     width = 10,
     height = 10,
     units = "in", 
     res = 300,
     compression = "lzw")

p <- ggcorrplot(cor_matrix,
                method = "circle",
                type = "upper",
                hc.order = FALSE,
                lab = TRUE,
                lab_size = 2.5,
                p.mat = p_matrix,
                insig = "blank",
                sig.level = 0.05,
                colors = c("blue", "white", "red"),
                ggtheme = ggplot2::theme_minimal,
                tl.cex = 8,
                tl.srt = 45) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 15
    ),
    axis.text.y = element_text(
      size = 15
    )
  ) +
  labs(title = "Metabolomics-Proteomics Correlation Heatmap")

print(p)
dev.off()
```