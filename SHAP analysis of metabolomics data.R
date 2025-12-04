```r
# 1. Load packages
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)
library(RColorBrewer)
library(pheatmap)

# 2. Data preprocessing
top_genes = c(
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

data_new2 = as.data.frame(data)[, top_genes]
data_new2 = cbind(data[, 2], data_new2)
data_new2 = as.data.frame(data_new2)
colnames(data_new2)[1] = 'Group'
data_new2$Group = ifelse(data_new2$Group == 0, "Control", "Treatment")
table(data_new2$Group)
data_new2$Group <- as.factor(data_new2$Group)
data_new2 = na.omit(data_new2)

train_set = data_new2[train_sub, ]
train_set = scale(train_set[, -1], center = TRUE, scale = TRUE)
train_set = as.data.frame(train_set)
train_set$Group = data_new2[train_sub, ][, 1]

test_set = data_new2[-train_sub, ]
test_labels <- test_set$Group
test_features <- test_set[, -1]
test_features = as.data.frame(scale(test_features, center = TRUE, scale = TRUE))

# 3. Multi-model training and ROC comparison
ml_methods <- data.frame(
  ModelName = c("RF", "SVM", "XGB", "GBM", "KNN"),
  MethodID = c("rf", "svmLinear", "xgboost", "gbm", "knn")
)

color_palette = c('#cb4936', '#5ac1b3', '#1616c8', 
                  '#A12568', '#F499C1', '#F7C394', '#B2A157', '#ade87c',
                  '#000000', '#03C4A1', '#bb4316', '#7382BC', '#F0E442',
                  '#3B185F', '#0d6c0d', '#FEC260', '#FD7014', '#a67063',
                  '#1B9B8A', '#D0EBE7', '#713045', '#F6E0EA', '#AD6D28',
                  '#EAB67D', '#5ac1b3', '#EE4590')

auc_results <- c()
model_auc_map <- list()
roc_colors <- color_palette[1:nrow(ml_methods)]
roc_ci_list <- list()

pdf(file = "model_roc_curve_port.pdf", width = 10, height = 8)
par(mar = c(5, 5, 4, 2) + 0.1)

for (i in 1:nrow(ml_methods)) {
  cv_folds <- 3
  mdl_name <- ml_methods$ModelName[i]
  mdl_id <- ml_methods$MethodID[i]
  
  if (nrow(train_set) < 10) stop("Training set too small!")
  
  message("Training model: ", mdl_name)
  
  if (mdl_id == "xgboost") {
    dtrain <- xgb.DMatrix(data = as.matrix(train_set[, -ncol(train_set)]),
                          label = as.numeric(train_set$Group) - 1)
    dtest <- xgb.DMatrix(data = as.matrix(test_features),
                         label = as.numeric(test_labels) - 1)
    
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = 4,
      eta = 0.1,
      nthread = 2,
      verbosity = 0
    )
    
    mdl <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = 100,
      verbose = 0
    )
    
    pred_prob <- predict(mdl, newdata = dtest)
    pred_prob <- data.frame(Treatment = pred_prob,
                            Control = 1 - pred_prob)
  }
  else if (mdl_id == "svmLinear") {
    mdl <- train(Group ~ ., data = train_set, method = mdl_id, prob.model = TRUE,
                 trControl = trainControl(method = "repeatedcv", number = cv_folds, savePredictions = FALSE))
    pred_prob <- predict(mdl, newdata = test_features, type = "prob")
  } else {
    mdl <- train(Group ~ ., data = train_set, method = mdl_id,
                 trControl = trainControl(method = "repeatedcv", number = cv_folds, savePredictions = FALSE))
    pred_prob <- predict(mdl, newdata = test_features, type = "prob")
  }
  
  if (!"Treatment" %in% colnames(pred_prob))
    stop("Prediction probabilities missing 'Treatment' column!")
  
  roc_obj <- roc(as.numeric(test_labels) - 1, as.numeric(pred_prob[, "Treatment"]))
  auc_val <- as.numeric(roc_obj$auc)
  auc_ci <- ci.auc(roc_obj, conf.level = 0.95)
  
  roc_ci_list[[mdl_id]] <- auc_ci
  auc_results <- c(auc_results,
                   sprintf("%s: %.03f [%.03f, %.03f]", mdl_name, auc_val, auc_ci[1], auc_ci[3]))
  model_auc_map[[mdl_id]] <- auc_val
  
  if (i == 1) {
    plot(roc_obj, print.auc = FALSE, legacy.axes = TRUE, main = "", col = roc_colors[i], lwd = 3)
  } else {
    plot(roc_obj, print.auc = FALSE, legacy.axes = TRUE, main = "", col = roc_colors[i], lwd = 3, add = TRUE)
  }
}

legend('bottomright', auc_results, col = roc_colors, lwd = 3, bty = 'n', cex = 0.9)
dev.off()

# 4. Select best model by AUC
auc_vec <- unlist(model_auc_map)
if (length(auc_vec) == 0) stop("AUC results empty!")
best_method <- names(which.max(auc_vec))
cat("Best model by AUC: ", best_method, "\n")

stopifnot(ncol(train_set) >= 2)
train_set$Group <- factor(train_set$Group)
if (!all(c("Treatment", "Control") %in% levels(train_set$Group))) {
  stop("Group must contain 'Treatment' and 'Control' levels")
}

train_set$Group <- relevel(train_set$Group, ref = "Control")
X_train <- train_set[, -ncol(train_set), drop = FALSE]

if (!all(vapply(X_train, is.numeric, logical(1)))) {
  stop("xgboost requires all numeric features")
}
stopifnot(!anyNA(X_train))
stopifnot(all(is.finite(as.matrix(X_train))))

# 5. SHAP value calculation
colnames(train_set) = c('AD_for_HDL',
                        'CE_to_TL_in_small_LDL',
                        'CE_to_TL_in_very_small_VLDL',
                        'Glucose_lactate',
                        'LA_to_TFA',
                        'Creatinine',
                        'Phospholipids_to_TL_in_CEL_VLDL',
                        'AD_for_VLDL',
                        'Glycine',
                        'Acetoacetate',
                        'Group')

colnames(X_train) = c('AD_for_HDL',
                      'CE_to_TL_in_small_LDL',
                      'CE_to_TL_in_very_small_VLDL',
                      'Glucose_lactate',
                      'LA_to_TFA',
                      'Creatinine',
                      'Phospholipids_to_TL_in_CEL_VLDL',
                      'AD_for_VLDL',
                      'Glycine',
                      'Acetoacetate')

dtrain <- xgb.DMatrix(as.matrix(X_train), label = as.numeric(train_set$Group) - 1)

params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.1,
  max_depth = 6
)

final_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 200,
  watchlist = list(train = dtrain),
  verbose = 1
)

shap_obj <- shapviz(final_model, X_pred = as.matrix(train_set[, -ncol(train_set)]), 
                    X = as.matrix(train_set[, -ncol(train_set)]), interactions = TRUE)

shap_importance <- sort(colMeans(abs(shap_obj$S)), decreasing = TRUE)
top_features <- names(shap_importance)

custom_theme <- theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        legend.position = "right")

options(shapviz.colors = color_palette)

print(sv_importance(shap_obj, kind = "bar", number_size = 7,
                    show_numbers = TRUE) + custom_theme)
ggsave("variables_importance.tiff", width = 11, height = 6, dpi = 300)

print(sv_importance(shap_obj, kind = "bee", number_size = 7,
                    show_numbers = TRUE) + custom_theme)
ggsave("variables_importance_bee.tiff", width = 11, 
       height = 6, dpi = 300)

p <- sv_dependence(shap_obj, v = top_features)  
ggsave("variables_importance_dependence.tiff", 
       plot = p,
       width = 22, height = 10, dpi = 400)

print(sv_force(shap_obj, row_id = 15, bar_label_size = 5, annotation_size = 5) + 
        custom_theme)
ggsave("variables_importance_force_plot.tiff", width = 11, 
       height = 6, dpi = 300)
```