```r
# Set language environment
Sys.setenv(LANGUAGE = "en")

# Set global options
options(stringsAsFactors = FALSE)

# Clear environment
rm(list = ls())

# Load required packages
library(dplyr)
library(tidyverse)
library(data.table)
library(tgp)
library(snowfall)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(xgboost)
library(caret)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(circlize)

# 1. Lasso algorithm
RunLasso <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

# 2. Ridge algorithm
RunRidge <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

# 3. Elastic Net regression
RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  cv.fit = cv.glmnet(x = Train_set, 
                     y = Train_label[[classVar]],
                     family = "binomial", alpha = alpha, nfolds = 10) 
  fit = glmnet(x = Train_set, y = Train_label[[classVar]],
               family = "binomial", alpha = alpha, 
               lambda = cv.fit$lambda.1se)
  fit$subFeature = colnames(Train_set)
  if(mode == "Model") return(fit)
  if(mode == "Variable") return(ExtractVar(fit))
}

# 4. Stepwise GLM
RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
  fit <- step(glm(formula = Train_label[[classVar]] ~ ., family = "binomial",
                  data = as.data.frame(Train_set)), direction = direction, 
              trace = 0)
  fit$subFeature = colnames(Train_set)
  if(mode == "Model") return(fit)
  if(mode == "Variable") return(ExtractVar(fit))
}

# 5. PLSRGLM
RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  cv.plsRglm.res = cv.plsRglm(formula = Train_label[[classVar]] ~ .,
                              data = as.data.frame(Train_set),
                              nt = 10, verbose = FALSE)
  fit <- plsRglm(Train_label[[classVar]], as.data.frame(Train_set), 
                 modele = "pls-glm-logistic", verbose = FALSE, sparse = TRUE) 
  fit$subFeature = colnames(Train_set)
  if(mode == "Model") return(fit)
  if(mode == "Variable") return(ExtractVar(fit))
}

# 6. Generalized Linear Model Boosting
RunglmBoost <- function(Train_set, Train_label, mode, classVar) {
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  
  formula <- as.formula(paste(classVar, "~ ."))
  fit <- glmboost(formula, 
                  data = data, 
                  family = Binomial())
  
  cvm <- cvrisk(fit, 
                papply = lapply,
                folds = cv(model.weights(fit), type = "kfold"))
  
  optimal_mstop <- max(mstop(cvm), 40)
  fit <- glmboost(formula,
                  data = data,
                  family = Binomial(),
                  control = boost_control(mstop = optimal_mstop))
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 7. Linear Discriminant Analysis
RunLDA <- function(Train_set, Train_label, mode, classVar) {
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  
  formula <- as.formula(paste(classVar, "~ ."))
  fit <- train(formula, 
               data = data, 
               method = "lda",
               trControl = trainControl(method = "cv"))
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 8. Naive Bayes
RunNaiveBayes <- function(Train_set, Train_label, mode, classVar) {
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  
  formula <- as.formula(paste(classVar, "~ ."))
  fit <- naiveBayes(formula, data = data)
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 9. Random Forest
RunRF <- function(Train_set, Train_label, mode, classVar) {
  rf_nodesize <- 30
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  
  fit <- rfsrc(
    formula = as.formula(paste(classVar, "~.")),
    data = cbind(Train_set, Train_label[classVar]),
    ntree = 300,
    nodesize = rf_nodesize,
    importance = TRUE,
    proximity = FALSE,
    forest = TRUE
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 10. Gradient Boosting Machine
RunGBM <- function(Train_set, Train_label, mode, classVar) {
  fit <- gbm(formula = Train_label[[classVar]] ~ .,
             data = as.data.frame(Train_set),
             distribution = 'bernoulli',
             n.trees = 1000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,
             n.cores = 6)
  
  best <- which.min(fit$cv.error)
  
  fit <- gbm(formula = Train_label[[classVar]] ~ .,
             data = as.data.frame(Train_set),
             distribution = 'bernoulli',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             n.cores = 8)
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 11. XGBoost
RunXGBoost <- function(Train_set, Train_label, mode, classVar) {
  indexes <- createFolds(Train_label[[classVar]], k = 5, list = TRUE)
  
  CV <- unlist(lapply(indexes, function(pt) {
    dtrain <- xgb.DMatrix(data = Train_set[-pt, ], 
                          label = Train_label[-pt, ])
    dtest <- xgb.DMatrix(data = Train_set[pt, ], 
                         label = Train_label[pt, ])
    watchlist <- list(train = dtrain, test = dtest)
    
    bst <- xgb.train(data = dtrain,
                     max.depth = 2,
                     eta = 1,
                     nthread = 4,
                     nrounds = 10,
                     watchlist = watchlist,
                     objective = "binary:logistic",
                     verbose = FALSE)
    which.min(bst$evaluation_log$test_logloss)
  }))
  
  nround <- as.numeric(names(which.max(table(CV))))
  fit <- xgboost(data = Train_set,
                 label = Train_label[[classVar]],
                 max.depth = 2,
                 eta = 1,
                 nthread = 4,
                 nrounds = nround,
                 objective = "binary:logistic",
                 verbose = FALSE)
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# 12. Support Vector Machine
RunSVM <- function(Train_set, Train_label, mode, classVar) {
  library(LiblineaR)
  
  X <- as.matrix(Train_set)
  y <- as.factor(Train_label[[classVar]])
  y_num <- as.numeric(y) - 1
  
  fit <- LiblineaR(
    data = X,
    target = y_num,
    type = 0,
    cost = 1,
    bias = TRUE
  )
  
  fit$subFeature <- colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(fit$W)
}

# 13. Run machine learning algorithm
RunML <- function(method, Train_set, Train_label, mode = "Model", classVar) {
  method <- gsub(" ", "", method)
  
  method_name <- gsub("(\\w+)\\[(.+)\\]", "\\1", method)
  method_param <- gsub("(\\w+)\\[(.+)\\]", "\\2", method)
  
  if (!is.matrix(Train_set)) {
    Train_set <- as.matrix(Train_set)
  }
  if (!is.numeric(Train_set)) {
    Train_set <- apply(Train_set, 2, as.numeric)
  }
  
  method_param <- switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL
  )
  
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, "; using ", ncol(Train_set), " Variables")
  
  args <- list("Train_set" = Train_set,
               "Train_label" = Train_label,
               "mode" = mode,
               "classVar" = classVar)
  args <- c(args, method_param)
  
  obj <- do.call(what = paste0("Run", method_name), args = args)
  
  if (mode == "Variable") {
    message(length(obj), " Variables retained;\n")
  } else {
    message("\n")
  }
  
  return(obj)
}

# 14. Data standardization function
standarize.fun <- function(indata, centerFlag, scaleFlag) { 
  scale(indata, center = centerFlag, scale = scaleFlag)
}

# Batch data scaling
scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL) {
  samplename <- rownames(data)
  
  if (is.null(cohort)) {
    data <- list(data)
    names(data) <- "training"
  } else {
    data <- split(as.data.frame(data), cohort)
  }
  
  if (is.null(centerFlags)) {
    centerFlags <- FALSE
    message("No centerFlags found, set as FALSE")
  }
  if (length(centerFlags) == 1) {
    centerFlags <- rep(centerFlags, length(data))
    message("set centerFlags for all cohort as ", unique(centerFlags))
  }
  if (is.null(names(centerFlags))) {
    names(centerFlags) <- names(data)
    message("match centerFlags with cohort by order\n")
  }
  
  if (is.null(scaleFlags)) {
    scaleFlags <- FALSE
    message("No scaleFlags found, set as FALSE")
  }
  if (length(scaleFlags) == 1) {
    scaleFlags <- rep(scaleFlags, length(data))
    message("set scaleFlags for all cohort as ", unique(scaleFlags))
  }
  if (is.null(names(scaleFlags))) {
    names(scaleFlags) <- names(data)
    message("match scaleFlags with cohort by order\n")
  }
  
  centerFlags <- centerFlags[names(data)]
  scaleFlags <- scaleFlags[names(data)]
  
  outdata <- mapply(standarize.fun, 
                    indata = data, 
                    centerFlag = centerFlags, 
                    scaleFlag = scaleFlags, 
                    SIMPLIFY = FALSE)
  
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  
  return(outdata)
}

# Extract important variables from models
ExtractVar <- function(fit) {
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"      = rownames(coef(fit))[which(coef(fit)[, 1] != 0)],
    "glm"         = names(coef(fit)),
    "svm.formula" = fit$subFeature,
    "train"       = fit$coefnames,
    "glmboost"    = names(coef(fit)[abs(coef(fit)) > 0]),
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs != 0],
    "rfsrc"       = names(which(fit$importance[, 1] > 0.01)),
    "gbm"         = rownames(summary.gbm(fit, plotit = FALSE))[summary.gbm(fit, plotit = FALSE)$rel.inf > 0],
    "xgb.Booster" = fit$subFeature,
    "naiveBayes"  = fit$subFeature
  ))
  
  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  
  return(Feature)
}

# Calculate prediction scores
CalPredictScore <- function(fit, new_data, type = "lp") {
  new_data <- new_data[, fit$subFeature, drop = FALSE]
  
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"       = predict(fit, type = "response", as.matrix(new_data)),
    "glm"          = predict(fit, type = "response", as.data.frame(new_data)),
    "svm.formula"  = {
      prob <- attr(predict(fit, as.data.frame(new_data), probability = TRUE), "probabilities")
      if ("1" %in% colnames(prob)) prob[, "1"] else prob[, 1]
    },
    "train"        = {
      prob <- predict(fit, new_data, type = "prob")
      if ("1" %in% colnames(prob)) prob[, "1"] else prob[, 1]
    },
    "glmboost"     = predict(fit, type = "response", as.data.frame(new_data)),
    "plsRglmmodel" = predict(fit, type = "response", as.data.frame(new_data)),
    "rfsrc"        = {
      pred <- predict(fit, as.data.frame(new_data))$predicted
      if (is.matrix(pred)) pred[, "1"] else pred
    },
    "gbm"          = predict(fit, type = "response", as.data.frame(new_data)),
    "xgb.Booster"  = predict(fit, as.matrix(new_data)),
    "naiveBayes"   = {
      prob <- predict(object = fit, type = "raw", newdata = new_data)
      if ("1" %in% colnames(prob)) prob[, "1"] else prob[, 1]
    }
  ))
  
  RS <- as.numeric(RS)
  
  if (length(RS) != nrow(new_data)) {
    warning(paste0("Prediction length mismatch in model [", class(fit)[1], "] → ",
                   "RS length: ", length(RS), " vs samples: ", nrow(new_data)))
    RS <- rep(NA, nrow(new_data))
  }
  
  names(RS) <- rownames(new_data)
  return(RS)
}

# Predict class
PredictClass <- function(fit, new_data) {
  new_data <- new_data[, fit$subFeature, drop = FALSE]
  
  label <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"       = predict(fit, type = "class", as.matrix(new_data)),
    "glm"          = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 0.5, "1", "0"),
    "svm.formula"  = predict(fit, as.data.frame(new_data), decision.values = TRUE),
    "train"        = predict(fit, new_data, type = "raw"),
    "glmboost"     = predict(fit, type = "class", as.data.frame(new_data)),
    "plsRglmmodel" = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 0.5, "1", "0"),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$class,
    "gbm"          = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 0.5, "1", "0"),
    "xgb.Booster"  = ifelse(predict(fit, as.matrix(new_data)) > 0.5, "1", "0"),
    "naiveBayes"   = predict(object = fit, type = "class", newdata = new_data)
  ))
  
  label <- as.character(label)
  
  if (length(label) != nrow(new_data)) {
    warning(paste0("PredictClass: model [", class(fit)[1], "] length mismatch → ",
                   "predicted: ", length(label), " vs actual: ", nrow(new_data)))
    label <- rep(NA, nrow(new_data))
  }
  
  names(label) <- rownames(new_data)
  return(label)
}

# Evaluate model performance
RunEval <- function(fit,
                    Test_set = NULL,
                    Test_label = NULL,
                    Train_set = NULL,
                    Train_label = NULL,
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    classVar) {
  
  if (!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  }
  
  if ((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(Train_set[, fit$subFeature],
                                 Test_set[, fit$subFeature])
    
    if (!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Training"
    }
    
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    Test_label <- rbind.data.frame(Train_label[, c(cohortVar, classVar)],
                                   Test_label[, c(cohortVar, classVar)])
    Test_label[, 1] <- factor(Test_label[, 1],
                              levels = c(unique(Train_label[, cohortVar]),
                                         setdiff(unique(Test_label[, cohortVar]),
                                                 unique(Train_label[, cohortVar]))))
  } else {
    new_data <- Test_set[, fit$subFeature]
  }
  
  RS <- suppressWarnings(as.numeric(CalPredictScore(fit = fit, new_data = new_data)))
  
  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  
  Predict.out <- split(x = Predict.out, f = Predict.out[, cohortVar])
  
  unlist(lapply(Predict.out, function(data) {
    if (!is.numeric(data$RS) || all(is.na(data$RS))) {
      warning("AUC skipped: prediction RS invalid.")
      return(NA)
    }
    
    if (length(unique(na.omit(data[[classVar]]))) < 2) {
      warning("AUC skipped: only one class in this cohort.")
      return(NA)
    }
    
    tryCatch({
      auc(suppressMessages(roc(response = data[[classVar]], predictor = as.numeric(data$RS))))
    }, error = function(e) {
      warning(paste0("AUC failed: ", e$message))
      return(NA)
    })
  }))
}

# Heatmap plotting function
SimpleHeatmap <- function(Cindex_mat, avg_Cindex,
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns, cluster_rows) {
  
  col_ha <- columnAnnotation(
    "Cohort" = colnames(Cindex_mat),
    col = list("Cohort" = CohortCol),
    show_annotation_name = FALSE
  )
  
  row_ha <- rowAnnotation(
    bar = anno_barplot(
      avg_Cindex,
      bar_width = 0.8,
      border = FALSE,
      gp = gpar(fill = barCol, col = NA),
      add_numbers = TRUE,
      numbers_offset = unit(-10, "mm"),
      axis_param = list("labels_rot" = 0),
      numbers_gp = gpar(fontsize = 9, col = "white"),
      width = unit(3, "cm")
    ),
    show_annotation_name = FALSE
  )
  
  Heatmap(
    as.matrix(Cindex_mat),
    name = "AUC",
    right_annotation = row_ha,
    top_annotation = col_ha,
    col = c("#4195C1", "#FFFFFF", "#FFBC90"),
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
    height = unit(cellheight * nrow(Cindex_mat), "cm"),
    column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
    column_title = NULL,
    cell_fun = function(j, i, x, y, w, h, col) {
      grid.text(
        label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
        x, y, gp = gpar(fontsize = 10)
      )
    }
  )
}

# Suppress output during execution
quiet <- function(..., messages = FALSE, cat = FALSE) {
  if (!cat) {
    sink(tempfile())
    on.exit(sink())
  }
  
  out <- if (messages) eval(...) else suppressMessages(eval(...))
  return(out)
}

########################################################################
########################################################################
# Data loading and splitting
data1 = fread("PCAD.csv")
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

id = which(data$PCAD == 1)
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
y_train <- Train_data$PCAD
y_test = Test_data$PCAD

cvfit = cv.glmnet(x = x_train, y = y_train,
                  family = 'binomial',
                  alpha = 1,
                  nfolds = 10)

coef_mat <- coef(cvfit, s = "lambda.1se")
coef_mat <- as.matrix(coef_mat)

lasso_features <- rownames(coef_mat)[which(as.numeric(coef_mat) != 0)]
lasso_features <- setdiff(lasso_features, c("(Intercept)", "Intercept"))

Train_set <- as.data.frame(Train_data)[, lasso_features, drop = FALSE]
Test_set <- as.data.frame(Test_data)[, lasso_features, drop = FALSE]

Train_class = as.data.frame(y_train)
row.names(Train_class) = row.names(Train_data)
colnames(Train_class) = 'Type'

Test_class = as.data.frame(y_test)
row.names(Test_class) = row.names(Test_data)
colnames(Test_class) = 'Type'

Test_class$Cohort = 'Testing'
Test_class = Test_class[,c(2,1)]

methodRT <- read.table("refer.txt", header = TRUE, sep = "\t", check.names = FALSE)
methods = methodRT$Model
methods <- gsub("-| ", "", methods)

classVar = "Type"
Variable = colnames(Train_set)
preTrain.method = strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))

# Variable selection using machine learning methods
preTrain.var <- list()
set.seed(seed = 123)

if (is.list(Train_set) && !is.data.frame(Train_set)) {
  Train_set <- as.data.frame(Train_set)
}

for (method in preTrain.method) {
  preTrain.var[[method]] <- RunML(
    method = method,
    Train_set = Train_set,
    Train_label = Train_class,
    mode = "Variable",
    classVar = classVar
  )
} 

preTrain.var[["simple"]] <- colnames(Train_set)

# Build machine learning models with selected variables
model <- list()
set.seed(seed = 123)
Train_set_bk <- Train_set
min.selected.var <- 2

for (method in methods) {
  cat(match(method, methods), ":", method, "\n")
  
  parts <- strsplit(method, "\\+")[[1]]
  if (length(parts) == 1) parts <- c("simple", parts)
  
  vars <- preTrain.var[[parts[1]]]
  
  if (length(vars) <= min.selected.var) {
    message("  SKIP ", parts[1], " → only ", length(vars), " variables\n")
    next
  }
  
  ts <- Train_set_bk[, vars, drop = FALSE]
  fit <- RunML(
    method = parts[2],
    Train_set = ts,
    Train_label = Train_class,
    mode = "Model",
    classVar = classVar
  )
  
  if (length(ExtractVar(fit)) <= min.selected.var) {
    message("  DROP ", method, " → only ", length(ExtractVar(fit)), " vars after modelling\n")
  } else {
    model[[method]] <- fit
  }
}

Train_set <- Train_set_bk
saveRDS(model, "model.MLmodel.rds")

# Load machine learning model
model <- readRDS("model.MLmodel.rds")
methodsValid <- names(model)

for (method in methodsValid) {
  cat("Checking:", method, "\n")
  fit <- model[[method]]
  
  train_check <- all(fit$subFeature %in% colnames(Train_set))
  test_check <- all(fit$subFeature %in% colnames(Test_set))
  
  cat("Train features:", train_check, "Test features:", test_check, "\n")
}

RS_list <- list()
for (method in methodsValid) {
  cat("Processing:", method, "\n")
  fit <- model[[method]]
  
  combined_data <- rbind.data.frame(
    Train_set[, fit$subFeature, drop = FALSE],
    Test_set[, fit$subFeature, drop = FALSE]
  )
  
  tryCatch({
    RS_list[[method]] <- CalPredictScore(fit = fit, new_data = combined_data)
    cat("✅ Success:", method, "\n")
  }, error = function(e) {
    cat("❌ ERROR in", method, "→", e$message, "\n")
    RS_list[[method]] <- NULL
  })
}

riskTab <- as.data.frame(t(do.call(rbind, RS_list)))
riskTab <- cbind(id = row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

Class_list <- list()
for (method in methodsValid) {
  fit <- model[[method]]
  combined_data <- rbind.data.frame(
    Train_set[, fit$subFeature, drop = FALSE],
    Test_set[, fit$subFeature, drop = FALSE]
  )
  
  Class_list[[method]] <- PredictClass(fit = fit, new_data = combined_data)
}

Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
classTab <- cbind(id = row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

fea_df <- lapply(methodsValid, function(method) {
  features <- ExtractVar(model[[method]])
  if (length(features) > 0) {
    data.frame(
      features = features,
      algorithm = method,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
})

fea_df <- do.call(rbind, fea_df)
write.table(fea_df, file = "model.genes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

AUC_list <- list()
for (method in methodsValid) {
  AUC_list[[method]] <- RunEval(
    fit = model[[method]],
    Test_set = Test_set,
    Test_label = Test_class,
    Train_set = Train_set,
    Train_label = Train_class,
    Train_name = "Training",
    cohortVar = "Cohort",
    classVar = classVar
  )
}

AUC_mat <- do.call(rbind, AUC_list)
aucTab <- cbind(Method = row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Plot AUC heatmap
library(RColorBrewer)
library(ComplexHeatmap)

AUC_mat <- read.table("model.AUCmatrix.txt", header = TRUE, sep = "\t", 
                      check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)

avg_AUC <- apply(AUC_mat, 1, mean, na.rm = TRUE)
avg_AUC <- sort(avg_AUC, decreasing = TRUE)
AUC_mat <- AUC_mat[names(avg_AUC), , drop = FALSE]

fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")[1:2]
names(CohortCol) <- colnames(AUC_mat)

cellwidth <- 1
cellheight <- 0.5

hm <- SimpleHeatmap(
  Cindex_mat = AUC_mat,
  avg_Cindex = avg_AUC,
  CohortCol = CohortCol,
  barCol = "steelblue",
  cellwidth = cellwidth,
  cellheight = cellheight,
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

pdf(file = "model.AUCheatmap.pdf", 
    width = cellwidth * ncol(AUC_mat) + 8, 
    height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Top 10 important genes heatmap
all_features <- unlist(fea_list, use.names = FALSE)
top_genes_tab <- sort(table(all_features), decreasing = TRUE)
top_genes <- names(top_genes_tab)[1:min(10, length(top_genes_tab))]

importance_mat <- matrix(0, nrow = length(methodsValid), ncol = length(top_genes))
rownames(importance_mat) <- methodsValid
colnames(importance_mat) <- top_genes

for (m in methodsValid) {
  fit <- model[[m]]
  cls <- class(fit)[1]
  imp <- rep(0, length(top_genes))
  names(imp) <- top_genes
  
  if (cls == "lognet") {
    coefs <- coef(fit)
    coef_val <- as.numeric(coefs)
    names(coef_val) <- rownames(coefs)
    imp_val <- coef_val[top_genes]
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
  } else if (cls == "glm") {
    coefs <- coef(fit)
    coef_val <- as.numeric(coefs)
    names(coef_val) <- names(coefs)
    imp_val <- coef_val[top_genes]
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
  } else if (cls == "glmboost") {
    coefs <- coef(fit)
    imp_val <- as.numeric(coefs[top_genes])
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
  } else if (cls == "rfsrc") {
    if (!is.null(fit$importance)) {
      use_genes <- intersect(top_genes, rownames(fit$importance))
      if (length(use_genes) > 0) {
        imp_val <- fit$importance[use_genes, 1]
        imp[use_genes] <- imp_val
      }
    }
  } else if (cls == "gbm") {
    gbm_imp <- summary(fit, plotit = FALSE)
    vals <- gbm_imp$rel.inf
    names(vals) <- gbm_imp$var
    imp_val <- vals[top_genes]
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
  } else if (cls == "xgb.Booster") {
    library(xgboost)
    imp_tab <- xgb.importance(feature_names = fit$subFeature, model = fit)
    vals <- imp_tab$Gain
    names(vals) <- imp_tab$Feature
    imp_val <- vals[top_genes]
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
  }
  
  importance_mat[m, ] <- imp
}

importance_mat[is.na(importance_mat)] <- 0
write.table(importance_mat, "model.importanceMatrix.txt", sep = "\t", quote = FALSE)

cat("All analyses completed successfully!\n")
```