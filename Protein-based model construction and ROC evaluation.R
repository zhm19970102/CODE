```r
# 1. Model building using top 10 proteins
# Data preprocessing
data_new = as.data.frame(data)[,top_genes]
data_new = cbind(data[,c(2:29)],data_new)
data_new = as.data.frame(data_new)
Train_data2 = data_new[train_sub,]
Test_data2 = data_new[-train_sub,]

colnames(data_new) = gsub(" ","_",colnames(data_new))
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

# Build logistic regression models
# Protein-only model
model1 <- glm(PCAD ~ .,
              data = Train_data2[,-c(2:28)], family = "binomial")

l <- summary(model1)
result1 <- l$coefficients
result1 <- data.frame(result1)

conf_int <- confint(model1)
result1$OR <- exp(result1$Estimate)
result1$low_95 <- exp(conf_int[, 1])
result1$upper_95 <- exp(conf_int[, 2])

colnames(result1) <- c("Estimate", "Std_Error", "z_value", "p_value",
                       "OR", "OR_low_95", "OR_upper_95")
result1$Variable <- rownames(result1)
result1 <- result1[, c("Variable", "Estimate", "Std_Error", "z_value", "p_value",
                       "OR", "OR_low_95", "OR_upper_95")]
write.csv(result1, "pcad_proteomics_results_only_prot.csv", row.names = FALSE)

# Protein model with clinical variables
colnames(Train_data2)
model2 <- glm(PCAD ~ Age + Sex + Ethnicity + Smoking + BMI +
                Hyperlipidemia + Hypertension + Diabetes +
                angptl3 + gast + lgals4 + ntprobnp + pigr + ren + fabp1 + met +
                bst2 + ccl21,
              data = Train_data2, family = "binomial")

l <- summary(model2)
result <- l$coefficients
result <- data.frame(result)

conf_int <- confint(model2)
result$OR <- exp(result$Estimate)
result$low_95 <- exp(conf_int[, 1])
result$upper_95 <- exp(conf_int[, 2])

colnames(result) <- c("Estimate", "Std_Error", "z_value", "p_value",
                      "OR", "OR_low_95", "OR_upper_95")
result$Variable <- rownames(result)
result <- result[, c("Variable", "Estimate", "Std_Error", "z_value", "p_value",
                     "OR", "OR_low_95", "OR_upper_95")]
write.csv(result, "pcad_proteomics_results_port_add_clinical.csv", row.names = FALSE)

# 2. ROC analysis on training set
library(pROC)
library(RColorBrewer)

pred_train1 <- predict(model1, newdata = Train_data2, type="response")
roc_train1 <- roc(Train_data2$PCAD, pred_train1)

pred_train2 <- predict(model2, newdata = Train_data2, type = "response")
roc_train2 <- roc(Train_data2$PCAD, pred_train2)

ci1 <- ci.auc(roc_train1)
ci2 <- ci.auc(roc_train2)

auc_text1 <- paste0("Model 1 AUC = ", round(auc(roc_train1), 3),
                    " (95% CI: ", round(ci1[1], 3), "-", round(ci1[3], 3), ")")
auc_text2 <- paste0("Model 2 AUC = ", round(auc(roc_train2), 3),
                    " (95% CI: ", round(ci2[1], 3), "-", round(ci2[3], 3), ")")

a <- brewer.pal(n = 3, name = "Set1")[1:2]

tiff("ROC_train_plot.tiff",
     width = 10,
     height = 8,
     units = "in",
     res = 300,
     compression = "lzw")

par(mar = c(5, 5, 4, 2) + 0.1,
    cex.axis = 1.2,
    cex.lab = 1.3,
    cex.main = 1.5)

plot(roc_train1, col = a[1], lwd = 3, main = "ROC Curve Comparison in Training Set", legacy.axes = TRUE)
plot(roc_train2, col = a[2], lwd = 3, add = TRUE)

legend("bottomright",
       legend = c(auc_text1, auc_text2),
       col = a,
       lwd = 3,
       cex = 0.9,
       bg = "white")
dev.off()

# 3. ROC analysis on testing set
model1 <- glm(PCAD ~ .,
              data = Test_data2[,-c(2:28)], family = "binomial")

l <- summary(model1)
result1 <- l$coefficients
result1 <- data.frame(result1)

conf_int <- confint(model1)
result1$OR <- exp(result1$Estimate)
result1$low_95 <- exp(conf_int[, 1])
result1$upper_95 <- exp(conf_int[, 2])

colnames(result1) <- c("Estimate", "Std_Error", "z_value", "p_value",
                       "OR", "OR_low_95", "OR_upper_95")
result1$Variable <- rownames(result1)
result1 <- result1[, c("Variable", "Estimate", "Std_Error", "z_value", "p_value",
                       "OR", "OR_low_95", "OR_upper_95")]
write.csv(result1, "pcad_proteomics_results_only_prot_testset.csv", row.names = FALSE)

model2 <- glm(PCAD ~ Age + Sex + Ethnicity + Smoking + BMI +
                Hyperlipidemia + Hypertension + Diabetes +
                angptl3 + gast + lgals4 + ntprobnp + pigr + ren + fabp1 + met +
                bst2 + ccl21,
              data = Test_data2, family = "binomial")

l <- summary(model2)
result <- l$coefficients
result <- data.frame(result)

conf_int <- confint(model2)
result$OR <- exp(result$Estimate)
result$low_95 <- exp(conf_int[, 1])
result$upper_95 <- exp(conf_int[, 2])

colnames(result) <- c("Estimate", "Std_Error", "z_value", "p_value",
                      "OR", "OR_low_95", "OR_upper_95")
result$Variable <- rownames(result)
result <- result[, c("Variable", "Estimate", "Std_Error", "z_value", "p_value",
                     "OR", "OR_low_95", "OR_upper_95")]
write.csv(result, "pcad_proteomics_results_port_add_clinical_testset.csv", row.names = FALSE)

pred_train1 <- predict(model1, newdata = Test_data2, type="response")
roc_train1 <- roc(Test_data2$PCAD, pred_train1)

pred_train2 <- predict(model2, newdata = Test_data2, type = "response")
roc_train2 <- roc(Test_data2$PCAD, pred_train2)

ci1 <- ci.auc(roc_train1)
ci2 <- ci.auc(roc_train2)

auc_text1 <- paste0("Model 1 AUC = ", round(auc(roc_train1), 3),
                    " (95% CI: ", round(ci1[1], 3), "-", round(ci1[3], 3), ")")
auc_text2 <- paste0("Model 2 AUC = ", round(auc(roc_train2), 3),
                    " (95% CI: ", round(ci2[1], 3), "-", round(ci2[3], 3), ")")

a <- brewer.pal(n = 3, name = "Set1")[1:2]

tiff("ROC_test_plot.tiff",
     width = 10,
     height = 8,
     units = "in",
     res = 300,
     compression = "lzw")

par(mar = c(5, 5, 4, 2) + 0.1,
    cex.axis = 1.2,
    cex.lab = 1.3,
    cex.main = 1.5)

plot(roc_train1, col = a[1], lwd = 3, main = "ROC Curve Comparison in Testing Set", legacy.axes = TRUE)
plot(roc_train2, col = a[2], lwd = 3, add = TRUE)

legend("bottomright",
       legend = c(auc_text1, auc_text2),
       col = a,
       lwd = 3,
       cex = 0.9,
       bg = "white")
dev.off()
```