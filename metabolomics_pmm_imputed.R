```r
library(data.table)

new1 <- fread("metebolomics.csv")

# Calculate missing rate per row
l1 = rowMeans(is.na(new1[, -1]))
id = which(l1 > 0.2)
new2 = new1[-id, ]

# Calculate missing rate per column
l2 = colMeans(is.na(new2[, -1]))
id2 = which(l2 > 0.2)
new3 = new2

# Perform multiple imputation
library(mice)
library(missForest)

new4 = new3[, -1]

id3 = c(1:ncol(new4))
r = cor(new4, use = "pairwise.complete.obs")

# Convert to data.frame
new4 <- as.data.frame(new4)
rm(new1)
rm(new2)

# Imputation process
while (length(id3) > 0) {
  i <- id3[1]
  batch_size <- 40  
  
  cat("Processing column", i, "(", colnames(new4)[i], ") - Batch size:", batch_size, 
      "- Remaining:", length(id3), "\n")
  
  # Select columns by correlation
  index1 <- order(abs(r[i, ]), decreasing = TRUE)
  candidate_cols <- colnames(new4)[index1]
  name1 <- candidate_cols[1:batch_size]
  
  # Ensure current column is included
  current_col <- colnames(new4)[i]
  if (!(current_col %in% name1)) {
    name1 <- c(current_col, name1[1:(batch_size-1)])
  }
  
  # Extract subset data
  data2 <- new4[, name1, drop = FALSE]
  
  # Multiple imputation
  imp1 <- mice(data2, m = 1, method = "pmm", maxit = 10, printFlag = FALSE)
  completed_data <- complete(imp1)
  
  # Update original data
  new4[, name1] <- completed_data
  
  # Update processing list
  processed_indices <- which(colnames(new4) %in% name1)
  id3 <- id3[!id3 %in% processed_indices]
}

cat("Imputation completed! Remaining columns with missing values:", length(id3), "\n")

fwrite(cbind(new3[, 1], new4), "pcad_metabolomics_imputed_pmm.csv")
```