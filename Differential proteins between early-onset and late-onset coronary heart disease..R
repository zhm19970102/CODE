```r
# Load data
data1 = fread("LCAD.csv")
sum(is.na(data1))

data2 = fread("proteomics_imputed_pmm.csv")
sum(duplicated(data1$eid))
sum(is.na(data2))

# Merge datasets
data = merge(data1, data2, by = "eid")
sum(duplicated(data$eid))
data = na.omit(data)

# Count PCAD patients
id = which(data$PCAD == 1)
length(id)

str(data)
rm(data1)
rm(data2)

# Set seed for reproducibility
set.seed(123)

# Split data into training and testing sets
train_sub = sample(nrow(data), 7/10 * nrow(data))
Train_data = data[train_sub,]
rownames(Train_data) = Train_data[[1]]
Train_data = Train_data[,-1]

Test_data = data[-train_sub,]
rownames(Test_data) = Test_data[[1]]
Test_data = Test_data[,-1]

# Prepare feature matrices
x_train = data.matrix(Train_data[,29:2947])
x_test = data.matrix(Test_data[,29:2947])

# Normalize data
x_train = scale(x_train, center = TRUE, scale = TRUE)
x_test = scale(x_test, center = TRUE, scale = TRUE)

# Prepare outcome variables
library(glmnet)
y_train <- Train_data$LCAD
y_test = Test_data$LCAD

# LASSO with 10-fold cross-validation
cvfit = cv.glmnet(x = x_train, y = y_train,
                  family = 'binomial',
                  alpha = 1,
                  nfolds = 10)

# Extract coefficients at lambda.1se
coef_mat <- coef(cvfit, s = "lambda.1se")
coef_mat <- as.matrix(coef_mat)

# Identify non-zero coefficient features
lasso_features <- rownames(coef_mat)[which(as.numeric(coef_mat) != 0)]
lasso_features <- setdiff(lasso_features, c("(Intercept)", "Intercept"))
lasso_features

# Define protein lists
PCAD_port_list = c(
  "acrv1", "adamtsl2", "adipoq", "angptl3",
  "apoc1", "asgr1", "bmp6", "bst2",
  "c8b", "ccl21", "cd101", "cd200",
  "cga", "ckmt1a_ckmt1b", "coq7", "cst5",
  "ctsb", "dcbld2", "dpep2", "edar",
  "fabp1", "fam20a", "fdx2", "fshb",
  "furin", "gast", "ggh", "ghrl",
  "gnas", "gpr158", "hgf", "hyou1",
  "ids", "itgav", "krt18", "lgals4",
  "lpl", "ment", "mep1a", "met",
  "mog", "mpo", "mxra8", "nrp1",
  "ntprobnp", "pam", "pcsk9", "pigr",
  "pla2g7", "pnpt1", "prrt3", "prss8",
  "pspn", "ren", "sdc4", "selenop",
  "sell", "serpina6", "sit1", "spint3",
  "tnni3"
)

LCAD_port_list = c(
  "ace2", "acrv1", "adgrg2", "adm", "ager", "agrp", "amot",
  "angptl3", "apoc1", "b3gnt7", "bcan", "ca14", "ccdc80", "cd300lg",
  "cdh2", "cdh3", "chga", "chrdl2", "clec3b", "cntn5", "cpa4",
  "ctsb", "ctsh", "ctsl", "cxcl17", "dcbld2", "dsg2", "dsg4",
  "eda2r", "egfr", "eln", "enpp2", "enpp6", "epha4", "faslg",
  "galnt2", "gast", "gdf15", "gfap", "gimap8", "havcr1", "hepacam2",
  "hs6st1", "hsd11b1", "hspb6", "igfbp3", "igfbp6", "il1rl1", "klk15",
  "klk3", "lcp1", "lect2", "lgals4", "lrrn1", "lrtm2", "lto1",
  "lypd3", "megf10", "ment", "met", "mmp12", "nefl", "ntprobnp",
  "odam", "paep", "pcsk9", "pgf", "plat", "pon3", "prap1",
  "prcp", "prrt3", "prss8", "pspn", "reck", "ren", "rnase6",
  "s100p", "serpina6", "tacstd2", "tex101", "vgf", "vwc2"
)

# Find intersection
inter_port <- intersect(PCAD_port_list, LCAD_port_list)
write.table(inter_port, file = 'inter_port.tsv', sep = '\t')

# Create Venn diagram
library(ggplot2)
library(UpSetR)
library(ggvenn)

venn_list <- list(
  PCAD_port_list = PCAD_port_list,
  LCAD_port_list = LCAD_port_list
)

p = ggvenn(
  venn_list,
  fill_color = c("#00AFBB", "#E7B800"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 5
)

ggsave("venn_plot.tiff", p, width = 10, height = 8, dpi = 300)

# Find unique proteins
PCAD_port_only = setdiff(PCAD_port_list, inter_port)
write.table(PCAD_port_only, file = 'PCAD_port_only.tsv', sep = '\t')

LCAD_port_only = setdiff(LCAD_port_list, inter_port)
write.table(LCAD_port_only, file = 'LCAD_port_only.tsv', sep = '\t')
```