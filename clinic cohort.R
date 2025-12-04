####clinic cohort####
####UKB cohort####
library(gtsummary)
library(survey)
library(haven)
library(tableone)
library(dplyr)
library(plyr)
library(tidyverse)
library(arsenal)
library(readr)
library(readxl)
library(mice)

mydata<- fread("pcad_hpv_end.csv")
str(mydata)
colnames(mydata)
cols_to_convert <- c("Race","Income","Smoking","Drinking" ,"Activity",                       
                     "Status" ,"History_of_hypertension","History_of_hyperlipidemia",      
                     "History_of_diabetes","BMI" ,"Waist","Haemoglobin_concentration","Lymphocyte_count","Monocyte_count",                
                     "Neutrophill_count","Platelet_count" ,"White_blood_cell_count" ,        
                     "CRP","Apolipoprotein_A","Apolipoprotein_B","Cholesterol","HDL","Triglycerides"                  
                     ,"LDL","Lipoprotein_A","IGF1","Cystatin_C")  

cols_to_convert <- gsub(" ", "_",cols_to_convert)
mydata=as.data.frame(mydata)
mydata[,cols_to_convert] <- lapply(mydata[,cols_to_convert], as.numeric)
na_proportion <- colSums(is.na(mydata)) / nrow(mydata)
print(na_proportion)

print(mydata)
imp <- mice(mydata[,(6:28)], m =1, method = 'pmm', maxit = 50, seed = 123)
summary(imp)
completed_data <- complete(imp, 1) 
print(completed_data)


mydata<- read_csv("HPV_UKB_END.csv")
mydata=HPV
str(mydata)
cols_to_convert <- c("Race",                    
                     "Education","Annual_household_income",  
                     "Smoking_status","Drinking_status","Status",                     
                     "HPV","History_of_hyperlipidemia","History_of_hypertension",  
                     "History_of_diabetes","Sex")  
mydata[cols_to_convert] <- lapply(mydata[cols_to_convert], as.character)
colnames(mydata)
mydata$HPV= factor(mydata$HPV)
mydata$Status= factor(mydata$Status)
mydata$Sex = factor(mydata$Sex)
mydata$Race = factor(mydata$Race)
mydata$Education = factor(mydata$Education)
mydata$Annual_household_income = factor(mydata$Annual_household_income)
mydata$Drinking_status= factor(mydata$Drinking_status)
mydata$Smoking_status = factor(mydata$Smoking_status)
mydata$History_of_hyperlipidemia= factor(mydata$History_of_hyperlipidemia)
mydata$History_of_hypertension = factor(mydata$History_of_hypertension)
mydata$History_of_diabetes = factor(mydata$History_of_diabetes)

vars1 <- c( "Status","Age","Sex","Race",                    
            "Education","Annual_household_income",  
            "Smoking_status","Drinking_status","BMI",                     
            "Waist","History_of_hyperlipidemia","History_of_hypertension",  
            "History_of_diabetes","Diastolic_pressure","Systolic_pressure",        
            "Cholesterol","HDL","LDL","TG","CRP","Hemoglobin","Lymphocyte_number",
            "Monocyte_number","Nneutrophils_number","Platelet_count", 
            "White_blood_cell_count")
table1 <- CreateTableOne(vars = vars1, strata = c("HPV"), data = mydata)

vars2 <- c("Age","BMI",                     
           "Waist","Diastolic_pressure","Systolic_pressure",        
           "Cholesterol","HDL","LDL","TG","CRP","Hemoglobin","Lymphocyte_number",
           "Monocyte_number","Nneutrophils_number","Platelet_count", 
           "White_blood_cell_count")

normality_tests <- lapply(mydata[,vars2], shapiro.test)

for (i in seq_along(normality_tests)) {
  cat("Variable:", vars2[i], "\n")
  cat("Test statistic:", normality_tests[[i]]$statistic, "\n")
  cat("p-value:", normality_tests[[i]]$p.value, "\n")
  cat("\n")
}
table2<-print(table1,nonnormal = c("Age","BMI",                     
                                   "Waist","Diastolic_pressure","Systolic_pressure",        
                                   "Cholesterol","HDL","LDL","TG","CRP","Hemoglobin","Lymphocyte_number",
                                   "Monocyte_number","Nneutrophils_number","Platelet_count", 
                                   "White_blood_cell_count"
),showAllLevels = TRUE)

#logistics回归
library(tidyverse)
library(gtsummary)
library(tidyr)
library(survey)
library(readr)
library(broom)
#单因素分析
mydata$Status <- relevel(mydata$Status, ref = "0")

variables <- c("HPV", "Age", "Sex", "Race",                    
               "Education", "Annual_household_income",  
               "Smoking_status", "Drinking_status", "BMI",                     
               "Waist", "History_of_hyperlipidemia", "History_of_hypertension",  
               "History_of_diabetes", "Diastolic_pressure", "Systolic_pressure",        
               "Cholesterol", "HDL", "LDL", "TG", "CRP", "Hemoglobin", "Lymphocyte_number",
               "Monocyte_number", "Nneutrophils_number", "Platelet_count", 
               "White_blood_cell_count"
)

results <- data.frame(
  Variable = character(),
  Level = character(),
  P_value = numeric(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  stringsAsFactors = FALSE
)

for (var in variables) {
  formula_str <- paste("Status ~", var)
  formula <- as.formula(formula_str)
  model2 <- glm(formula, data = mydata, family = binomial())
  
model_summary <- summary(model2)
  
  
  if (is.factor(mydata[[var]]) | is.character(mydata[[var]])) {
   
    n_levels <- length(levels(as.factor(mydata[[var]])))
    
    if (n_levels > 1) {
      for (i in 2:n_levels) {  
      
        P_value <- model_summary$coefficients[i, "Pr(>|z|)"]
        
       
        OR_value <- exp(coef(model2))[i]
        
   
        CI_values <- exp(confint(model2))[i, ]
        CI_lower <- CI_values[1]
        CI_upper <- CI_values[2]
        
       
        level_name <- levels(as.factor(mydata[[var]]))[i]
        
     
        results <- rbind(results, data.frame(
          Variable = var,
          Level = level_name,
          P_value = P_value,
          OR = OR_value,
          CI_lower = CI_lower,
          CI_upper = CI_upper
        ))
      }
    }
  } else {
  
    P_value <- model_summary$coefficients[2, "Pr(>|z|)"]
    
    
    OR_value <- exp(coef(model2))[2]
    
    CI_values <- exp(confint(model2))[2, ]
    CI_lower <- CI_values[1]
    CI_upper <- CI_values[2]
    
    
    results <- rbind(results, data.frame(
      Variable = var,
      Level = "Continuous",
      P_value = P_value,
      OR = OR_value,
      CI_lower = CI_lower,
      CI_upper = CI_upper
    ))
  }
}

#多因素logistic回归
model3 <- glm(Status ~ HPV + Age + Sex +Smoking_status +Drinking_status + BMI + History_of_hyperlipidemia +
                History_of_hypertension + History_of_diabetes + White_blood_cell_count,
              data = mydata, family = binomial())

model_summary2 <- summary(model3)  


P <- model_summary2$coefficients[, "Pr(>|z|)"] 
OR <- exp(coef(model3))
CI <- exp(confint(model3))

results <- data.frame(
  Variable = names(OR),
  P_value = P,
  OR = OR,
  CI_lower = CI[, 1],
  CI_upper = CI[, 2],
  row.names = NULL
)


#多模型
#HPV暴露
pcad=mydata
model1 <- glm(PCAD ~ HPV + Age + Sex + Ethnicity,
              data = pcad, family = binomial())

model_summary1 <- summary(model1)  


P <- model_summary1$coefficients[, "Pr(>|z|)"] 
OR <- exp(coef(model1))
CI <- exp(confint(model1))


model2 <- glm(PCAD ~ HPV + Age + Sex + Ethnicity+
                Education+Average_total_household_income_before_tax,
              data = pcad, family = binomial())

model_summary2 <- summary(model2)  


P <- model_summary2$coefficients[, "Pr(>|z|)"] 
OR <- exp(coef(model2))
CI <- exp(confint(model2))


model3 <- glm(PCAD ~ HPV + Age + Sex + Ethnicity+
                Education+Average_total_household_income_before_tax+
                Smoking+BMI+Hypertension+Diabetes+
                Hyperlipidemia,
              data = pcad, family = binomial())

model_summary3 <- summary(model3)  


P <- model_summary3$coefficients[, "Pr(>|z|)"] 
OR <- exp(coef(model3))
CI <- exp(confint(model3))
#高脂血症病史为暴露
model1 <- glm(PCAD ~ Hyperlipidemia + Age + Sex + Ethnicity,
              data = pcad, family = binomial())

model_summary1 <- summary(model1)  


P <- model_summary1$coefficients[, "Pr(>|z|)"] 
OR <- exp(coef(model1))
CI <- exp(confint(model1))


model2 <- glm(PCAD ~ Hyperlipidemia + Age + Sex + Ethnicity+
                Education+Average_total_household_income_before_tax,
              data = pcad, family = binomial())

model_summary2 <- summary(model2)  


P <- model_summary2$coefficients[, "Pr(>|z|)"] 
OR <- exp(coef(model2))
CI <- exp(confint(model2))


model3 <- glm(PCAD ~ Hyperlipidemia + Age + Sex + Ethnicity+
                Education+Average_total_household_income_before_tax+
                Smoking+BMI+Hypertension+Diabetes+
                Hyperlipidemia,
              data = pcad, family = binomial())

model_summary3 <- summary(model3)  

果
P <- model_summary3$coefficients[, "Pr(>|z|)"] 
OR <- exp(coef(model3))
CI <- exp(confint(model3))
#亚组分析
library(jstable)   
library(survival)  
library(forestploter)  
library(grid)  
library(foreign)
library(tidyverse)  

str(pcad)

pcad$PCAD=as.numeric(pcad$PCAD)
pcad$PCAD = pcad$PCAD-1
 
pcad <- pcad %>%  
  rename('Annual_household_income' = 'Average_total_household_income_before_tax')  


res <- TableSubgroupMultiGLM(
  formula = PCAD ~ HPV,  
  var_subgroups = c("Sex", "Education", "Annual_household_income", 
                    "Smoking", "Hypertension", "Diabetes", "Hyperlipidemia"), 
  data = pcad
)
res

plot_df <- res  
colnames(plot_df)
plot_df[, c(2, 3, 7, 8)][is.na(plot_df[, c(2, 3, 7, 8)])] <- " "  
plot_df$'' <- paste(rep(" ", nrow(plot_df)), collapse = "")  
plot_df[, 4:6] <- apply(plot_df[, 4:6], 2, as.numeric)  
plot_df$"OR (95% CI)" <- ifelse(is.na(plot_df$OR), "",
                                sprintf("%.2f (%.2f to %.2f)",
                                        plot_df$OR, plot_df$Lower, plot_df$Upper))  

plot_df  
colnames(plot_df)
p <- forest(
  data = plot_df[, c(1, 2, 3, 9, 10, 7, 8)],  
  lower = plot_df$Lower,  
  upper = plot_df$Upper,  
  est = plot_df$OR,  
  ci_column = 4,  
  ref_line = 1,  
  x_trans = c("log"),
  ticks_at = c(1, 7.4, 54.6, 403.4)
)
plot(p)  
ggsave("pcad_hpv_ukb_forest_plot.tiff", p, width = 10, height = 8, dpi = 300)
####NHANES cohort####
#1.数据处理及因子化
#加载包
library(gtsummary)
library(survey)
library(haven)
library(tableone)
library(dplyr)
library(plyr)
library(tidyverse)
library(arsenal)
library(readr)
#导入HPV数据
hpv_eocad <- read_csv("hpv_eocad_end.csv")
mydata=hpv_eocad
str(mydata)
colnames(mydata)
mydata$HPV= factor(mydata$HPV)
mydata$Eocad = factor(mydata$Eocad)
mydata$Sex = factor(mydata$Sex)
mydata$Race = factor(mydata$Race)
mydata$Education = factor(mydata$Education)
mydata$Marital_status = factor(mydata$Marital_status)
mydata$Annual_household_income = factor(mydata$Annual_household_income)
mydata$Drinking_status= factor(mydata$Drinking_status)
mydata$Smoking_status = factor(mydata$Smoking_status)
mydata$History_of_hyperlipidemia= factor(mydata$History_of_hyperlipidemia)
mydata$History_of_hypertension = factor(mydata$History_of_hypertension)
mydata$History_of_diabetes = factor(mydata$History_of_diabetes)
#数据加权
NHANES_design = survey::svydesign(data = as.data.frame(mydata),
                                  ids = ~ SDMVPSU.x,
                                  strata = ~ SDMVSTRA.x,
                                  nest = T,
                                  weights = ~ weight.new,
                                  survey.lonely.psu = 'adjust')
colnames(mydata)
tbl = tbl_svysummary(NHANES_design,
                     by = 'Eocad',
                     include = c(HPV,Sex,Age,Race,Education,Marital_status,
                                 Annual_household_income,Drinking_status,
                                 Smoking_status,BMI,Waist,History_of_hyperlipidemia,
                                 History_of_hypertension,History_of_diabetes,
                                 Systolic_pressure,Diastolic_pressure,HDL,TC,TG,LDL,
                                 Cotinine,White_blood_cell_count,Lymphocyte_percent,           
                                 Monocyte_percent,Segmented_neutrophils_percent,Lymphocyte_number,            
                                 Monocyte_number,Segmented_neutrophils_number,Hemoglobin,                  
                                 Platelet_count),
                     statistic = list(all_continuous() ~ '{median} ({p25}, {p75})',
                                      all_categorical() ~ '{n_unweighted} ({p}%)'))%>%
  add_overall()%>%
  add_p(pvalue_fun = function(x) style_number(x, digits = 3))%>%
  add_n(statistic = '{N_nonmiss_unweighted}',
        col_label = '**N**',
        footnote = T)
tbl %>% 
  as_flex_table() %>%
  flextable::save_as_docx(path = 'hpv_eocad_baseline.docx')
#logistic回归分析
library(tidyverse)
library(gtsummary)
library(tidyr)
library(survey)
library(readr)
hpv_eocad_end <- read_csv("hpv_eocad_end.csv")
mydata=hpv_eocad_end
NHANES_design = svydesign(data = mydata,
                          ids = ~ SDMVPSU.x,
                          strata = ~ SDMVSTRA.x,
                          nest = T,
                          weights = ~ weight.new,
                          survey.lonely.psu = 'adjust')

#单因素logistic
colnames(mydata)
vector1= names(mydata)[c(27:36)]
vector= names(mydata)[-c(1,3,4,5,12,13,14:26,37)]

for (i in 1:length(vector)) {
  formula <- as.formula(paste("Eocad ~", vector[i]))
  m <- svyglm(formula, design = NHANES_design, family = binomial)
  tb <- tbl_regression(m, exponentiate = TRUE,
                       pvalue_fun = function(x) style_number(x, digits = 3)) 
  print(tb)
}

#多因素logistic回归
m = svyglm(Eocad ~ HPV+Sex+Age+Education+Race+Marital_status+
             Annual_household_income+
             Smoking_status+BMI+
             History_of_hypertension+History_of_diabetes+
             History_of_hyperlipidemia+
             Monocyte_number,
           design = NHANES_design, family = binomial)
summary(m)
tb <- tbl_regression(m, exponentiate = TRUE,
                     pvalue_fun = function(x) style_number(x, digits = 3))
tb
#多模型
#HPV为暴露
#未调整
colnames(mydata)
m = svyglm(Eocad ~ HPV, design = NHANES_design, family = binomial)
summary(m)
tb <- tbl_regression(m, exponentiate = TRUE,
                     pvalue_fun = function(x) style_number(x, digits = 3))
tb
#Model1
m1 = svyglm(Eocad ~ HPV+Sex+Age+Race
            ,design = NHANES_design, family = binomial)
tb1 <- tbl_regression(m1, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb1
#Model2
m2 = svyglm(Eocad ~ HPV+Sex+Age+Race+Education+Annual_household_income
            ,design = NHANES_design, family = binomial)
tb2 <- tbl_regression(m2, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb2
#Model3
m3 = svyglm(Eocad ~ HPV+Sex+Age+Race+Education+Annual_household_income+
              Smoking_status+BMI+History_of_diabetes+History_of_hypertension+History_of_hyperlipidemia,
            design = NHANES_design, family = binomial)
tb3 <- tbl_regression(m3, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb3

#HDL为暴露
colnames(mydata)
m = svyglm(Eocad ~ HDL, design = NHANES_design, family = binomial)
summary(m)
tb <- tbl_regression(m, exponentiate = TRUE,
                     pvalue_fun = function(x) style_number(x, digits = 3))
tb
#Model1
m1 = svyglm(Eocad ~ HDL+Sex+Age+Race
            ,design = NHANES_design, family = binomial)
tb1 <- tbl_regression(m1, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb1
#Model2
m2 = svyglm(Eocad ~ HDL+Sex+Age+Race+Education+Annual_household_income
            ,design = NHANES_design, family = binomial)
tb2 <- tbl_regression(m2, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb2
#Model3
m3 = svyglm(Eocad ~ HDL+Sex+Age+Race+Education+Annual_household_income+
              Smoking_status+BMI+
              History_of_diabetes+History_of_hypertension,
            design = NHANES_design, family = binomial)
tb3 <- tbl_regression(m3, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb3

#高脂血症病史
colnames(mydata)
m = svyglm(Eocad ~ History_of_hyperlipidemia, design = NHANES_design, family = binomial)
summary(m)
tb <- tbl_regression(m, exponentiate = TRUE,
                     pvalue_fun = function(x) style_number(x, digits = 3))
tb
#Model1
m1 = svyglm(Eocad ~ History_of_hyperlipidemia+Sex+Age+Race
            ,design = NHANES_design, family = binomial)
tb1 <- tbl_regression(m1, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb1
#Model2
m2 = svyglm(Eocad ~ History_of_hyperlipidemia+Sex+Age+Race+Education+Annual_household_income
            ,design = NHANES_design, family = binomial)
tb2 <- tbl_regression(m2, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb2
#Model3
m3 = svyglm(Eocad ~ History_of_hyperlipidemia+Sex+Age+Race+Education+Annual_household_income+
              Smoking_status+BMI+
              History_of_diabetes+History_of_hypertension,
            design = NHANES_design, family = binomial)
tb3 <- tbl_regression(m3, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb3

#LDL为暴露
colnames(mydata)
m = svyglm(Eocad ~ LDL, design = NHANES_design, family = binomial)
summary(m)
tb <- tbl_regression(m, exponentiate = TRUE,
                     pvalue_fun = function(x) style_number(x, digits = 3))
tb
#Model1
m1 = svyglm(Eocad ~ LDL+Sex+Age+Race
            ,design = NHANES_design, family = binomial)
tb1 <- tbl_regression(m1, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb1
#Model2
m2 = svyglm(Eocad ~ LDL+Sex+Age+Race+Education+Annual_household_income
            ,design = NHANES_design, family = binomial)
tb2 <- tbl_regression(m2, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb2
#Model3
m3 = svyglm(Eocad ~ LDL+Sex+Age+Race+Education+Annual_household_income+
              Smoking_status+BMI+
              History_of_diabetes+History_of_hypertension,
            design = NHANES_design, family = binomial)
tb3 <- tbl_regression(m3, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb3

#TG为暴露
colnames(mydata)
m = svyglm(Eocad ~ TG, design = NHANES_design, family = binomial)
summary(m)
tb <- tbl_regression(m, exponentiate = TRUE,
                     pvalue_fun = function(x) style_number(x, digits = 3))
tb
#Model1
m1 = svyglm(Eocad ~ TG+Sex+Age+Race
            ,design = NHANES_design, family = binomial)
tb1 <- tbl_regression(m1, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb1
#Model2
m2 = svyglm(Eocad ~ TG+Sex+Age+Race+Education+Annual_household_income
            ,design = NHANES_design, family = binomial)
tb2 <- tbl_regression(m2, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb2
#Model3
m3 = svyglm(Eocad ~ TG+Sex+Age+Race+Education+Annual_household_income+
              Smoking_status+BMI+
              History_of_diabetes+History_of_hypertension,
            design = NHANES_design, family = binomial)
tb3 <- tbl_regression(m3, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb3

#TC为暴露
colnames(mydata)
m = svyglm(Eocad ~ TC, design = NHANES_design, family = binomial)
summary(m)
tb <- tbl_regression(m, exponentiate = TRUE,
                     pvalue_fun = function(x) style_number(x, digits = 3))
tb
#Model1
m1 = svyglm(Eocad ~ TC+Sex+Age+Race
            ,design = NHANES_design, family = binomial)
tb1 <- tbl_regression(m1, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb1
#Model2
m2 = svyglm(Eocad ~ TC+Sex+Age+Race+Education+Annual_household_income
            ,design = NHANES_design, family = binomial)
tb2 <- tbl_regression(m2, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb2
#Model3
m3 = svyglm(Eocad ~ TC+Sex+Age+Race+Education+Annual_household_income+
              Smoking_status+BMI+
              History_of_diabetes+History_of_hypertension,
            design = NHANES_design, family = binomial)
tb3 <- tbl_regression(m3, exponentiate = TRUE,
                      pvalue_fun = function(x) style_number(x, digits = 3))
tb3

#亚组分析
library(survey)
library(forestplot)
library(dplyr)

NHANES_design <- svydesign(data = mydata,
                           ids = ~ SDMVPSU.x,
                           strata = ~ SDMVSTRA.x,
                           nest = TRUE,
                           weights = ~ weight.new,
                           survey.lonely.psu = 'adjust')


calc_odds_ratio <- function(design, exposure, outcome) {
  svyglm_obj <- svyglm(as.formula(paste(outcome, "~", exposure)), 
                       design = design, 
                       family = quasibinomial())
  coefs <- coef(summary(svyglm_obj))
  or <- exp(coefs[2,1])
  ci_lower <- exp(coefs[2,1] - 1.96 * coefs[2,2])
  ci_upper <- exp(coefs[2,1] + 1.96 * coefs[2,2])
  p_value <- 2 * (1 - pnorm(abs(coefs[2,1] / coefs[2,2])))
  data.frame(OR = or, lower = ci_lower, upper = ci_upper, P = p_value)
}


results <- data.frame(Variable = character(),
                      Level = character(),
                      OR = numeric(),
                      or_lci95 = numeric(),
                      or_uci95 = numeric(),
                      P = numeric())


factor_vars <- c("Sex", "Race", "Education", "Marital_status", "Annual_household_income","Smoking_status","History_of_hypertension"      
                 ,"History_of_hyperlipidemia","History_of_diabetes")

for (var in factor_vars) {
  levels <- levels(mydata[[var]])
  for (level in levels) {
    subset_design <- subset(NHANES_design, mydata[[var]] == level)
    res <- calc_odds_ratio(subset_design, "HPV", "Eocad")
    res$Variable <- var
    res$Level <- level
    results <- rbind(results, res)
  }
}



interaction_p_values <- c()
variables = c("Sex", "Race", "Education", "Marital_status", "Annual_household_income","Smoking_status","History_of_hypertension"      
              ,"History_of_hyperlipidemia","History_of_diabetes")
for (var in variables) {
  formula_int <- as.formula(paste("Eocad ~ HPV *", var))
  svy_glm_int <- svyglm(formula_int, design = NHANES_design, family = quasibinomial())
  p_int <- summary(svy_glm_int)$coefficients[4,4]
  interaction_p_values <- c(interaction_p_values, p_int)
}
results$p_for_interaction <- c(interaction_p_values, rep(NA,13)) 

dt = read.csv(file = 'HPV_results.csv')

names(dt)[6] = "P for interaction"
dt$P = round(dt$P, digits = 3)
dt$`P for interaction` = round(dt$`P for interaction`, digits = 3)
dt$' ' = paste(rep(' ',22),collapse = ' ')
dt$'OR(95%CI)' = ifelse(is.na(dt$OR),'',
                        sprintf('%.2f(%.2f to %.2f)',
                                dt$OR,dt$lower,dt$upper))
dt$P[is.na(dt$P)] = ' '
dt$`P for interaction`[is.na(dt$`P for interaction`)] = ' '
library(grid)
library(forestplot)
library(readxl)
library(forestploter)
tm <- forest_theme(base_size = 10,
                   ci_pch = 20,
                   ci_col = 'orange',
                   ci_lty = 1,
                   ci_lwd = 2.3,
                   ci_Theight = 0.2,
                   refline_lwd = 1.5,
                   refline_lty = 'dashed',
                   refline_col = 'red',
                   summary_fill = '#4575b4',
                   summary_col = '#4575b4',
                   footnote_cex = 1.1,
                   footnote_fontface = 'italic',
                   footnote_col = 'blue')

p <- forest(dt[, c(5,8,7,4,6)],
            est = dt$OR,
            lower = dt$lower,
            upper = dt$upper,
            sizes = 0.6,
            ci_column = 3,
            ref_line = 1,
            xlim = c(0.5, 2.0),  
            ticks_at = c(0.5, 1, 1.5, 2.0),
            arrow_lab = c('protective factor', 'risk factor'),
            footnote = '',
            theme = tm)

print(p)
tiff(filename = 'HPV_Eocad.tif',width = 4500, height = 2800,
     res = 350)
print(p)
dev.off()
####China cohort####
library(foreign) 
library(rms)
library(tableone)
library(broom)
library(pROC)
library(regplot)
library(ResourceSelection)
library(nricens)
library(PredictABEL)
library(ggDCA)
library(devtools)
library(survivalROC) 
library(caret)
library(openxlsx)

PCAD <- read.xlsx("PCAD2.xlsx") 

str(PCAD)

PCAD$Sex <- factor(PCAD$Sex,labels = c("Female","Male"))
PCAD$Smoking <- factor(PCAD$Smoking,labels = c("No","Yes"))
PCAD$Drinking <- factor(PCAD$Drinking,labels = c("No","Yes"))
PCAD$Family_history <- factor(PCAD$Family_history,labels = c("No","Yes"))
PCAD$Hypertension <- factor(PCAD$Hypertension,labels = c("No","Yes"))
PCAD$Diabetes <- factor(PCAD$Diabetes,labels = c("No","Yes"))
PCAD$HBV <- factor(PCAD$HBV,labels = c("No","Yes"))
PCAD$Group <- factor(PCAD$Group,labels = c("Control","Case"))
PCAD$CAG <- factor(PCAD$CAG)
str(PCAD)


PCAD_impute = PCAD[,-(19:22)]
PCAD_case = PCAD[PCAD$Group == "Case",]


library(mice)
imp1 <- mice(PCAD_impute, m = 1, method = "pmm", maxit = 50, print = FALSE)
imputed_data <- complete(imp1) 
sum(is.na(imputed_data))


vars1 <- colnames(imputed_data)[-18]
table1 <- CreateTableOne(vars = vars1, strata = c("Group"), data = imputed_data)


vars2 <- c("Age","BMI","TC","TG","HDL","LDL","apoA","apoB",
           "eFGR","hsCRP")

normality_tests <- lapply(imputed_data[,vars2], shapiro.test)

for (i in seq_along(normality_tests)) {
  cat("Variable:", vars2[i], "\n")
  cat("Test statistic:", normality_tests[[i]]$statistic, "\n")
  cat("p-value:", normality_tests[[i]]$p.value, "\n")
  cat("\n")
}

table2 <- print(table1,nonnormal = vars2,showAllLevels = TRUE)
write.csv(table2,'PCAD_baseline.csv')


imp2 <- mice(PCAD_case, m = 1, method = "pmm", maxit = 50, print = FALSE)
imputed_data2 <- complete(imp2) 
sum(is.na(imputed_data2))


vars3 <- colnames(imputed_data2)[19:22]
vars4 <- c("LVEDd","LAD","LVEF")

normality_tests <- lapply(imputed_data2[,vars4], shapiro.test)

for (i in seq_along(normality_tests)) {
  cat("Variable:", vars4[i], "\n")
  cat("Test statistic:", normality_tests[[i]]$statistic, "\n")
  cat("p-value:", normality_tests[[i]]$p.value, "\n")
  cat("\n")
}

table3 <- CreateTableOne(vars = vars3, data = imputed_data2)
table4 = print(table3, nonnormal = c("LVEDd","LAD","LVEF",
                                     showAllLevels = TRUE))
write.csv(table4,file = 'PACD_features.csv')


colnames(imputed_data)
  
variables <- colnames(imputed_data)[-18]    

#install.packages("logistf")
library(logistf)

results_df <- data.frame(
  Variable = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (var in variables) {  
  formula_str <- paste("Group ~", var)  
  formula <- as.formula(formula_str)  
  model <- logistf(formula, data = imputed_data)  # logistf默认就是二项分布
  

  OR <- exp(coef(model)[2])
  CI <- exp(confint(model)[2, ])
  P <- model$prob[2]  # 直接从模型对象获取p值
  
 
  results_df <- rbind(results_df, data.frame(
    Variable = var,
    OR = round(OR, 3),
    CI_lower = round(CI[1], 3),
    CI_upper = round(CI[2], 3),
    P_value = ifelse(P < 0.001, "<0.001", round(P, 3))
  ))
}


print(results_df)
write.csv(results_df,file = 'Univariate_logistic_regression.csv')



colnames(imputed_data)
model2 <- logistf(Group ~ ., data = imputed_data)


results_df2 <- data.frame(
  Variable = names(coef(model2))[-1],  # 去掉截距项
  OR = exp(coef(model2)[-1]),
  CI_lower = exp(confint(model2)[-1, 1]),
  CI_upper = exp(confint(model2)[-1, 2]),
  P_value = model2$prob[-1]
)


results_df2$OR <- round(results_df2$OR, 3)
results_df2$CI_lower <- round(results_df2$CI_lower, 3)
results_df2$CI_upper <- round(results_df2$CI_upper, 3)
results_df2$P_value <- ifelse(results_df2$P_value < 0.001, "<0.001", 
                              round(results_df2$P_value, 3))

A = print(results_df2)
write.csv(A,file = 'multivariate_logistic_regression.csv')


#Model1
vars3 = c('TC','TG','HDL','LDL','apoA','apoB')
results_df_model1 <- data.frame(
  Variable = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (var in vars3) {  

  formula_str <- paste("Group ~", var, "+ Sex + Age")
  formula <- as.formula(formula_str)
  
  model <- logistf(formula, data = imputed_data)
  
  
  OR = exp(coef(model)[2])
  CI_lower = exp(confint(model)[2, 1])
  CI_upper = exp(confint(model)[2, 2])
  P_value = model$prob[2]
  
  
  results_df_model1 <- rbind(results_df_model1, data.frame(
    Variable = var,
    OR = round(OR, 3),
    CI_lower = round(CI_lower, 3),
    CI_upper = round(CI_upper, 3),
    P_value = ifelse(P < 0.001, "<0.001", round(P, 3))
  ))
}
write.csv(results_df_model1,file = 'model1_logistic_regression.csv')

#Model2
results_df_model2 <- data.frame(
  Variable = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (var in vars3) {  
  
  formula_str <- paste("Group ~", var, "+ Sex + Age + BMI + 
                       Smoking + Family_history + Hypertension + 
                       Diabetes")
  formula <- as.formula(formula_str)
  
  model <- logistf(formula, data = imputed_data)
  
  
  OR = exp(coef(model)[2])
  CI_lower = exp(confint(model)[2, 1])
  CI_upper = exp(confint(model)[2, 2])
  P_value = model$prob[2]
  
  
  results_df_model2 <- rbind(results_df_model2, data.frame(
    Variable = var,
    OR = round(OR, 3),
    CI_lower = round(CI_lower, 3),
    CI_upper = round(CI_upper, 3),
    P_value = ifelse(P < 0.001, "<0.001", round(P, 3))
  ))
}
write.csv(results_df_model2,file = 'model2_logistic_regression.csv')

