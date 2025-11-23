library(data.table)

# ivn function
inverse_normal_transform <- function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# Step 1: Load the proteome data and do the inverse normal transform
proteome_data <- read.table("proteome_for_plink.txt", header = T, sep = "\t")
rownames(proteome_data) <- proteome_data[, 2]
proteome_data <- proteome_data[, 3:367]
ivn_proteome_data <- sapply(proteome_data, inverse_normal_transform)
ivn_proteome_data <- as.data.frame(ivn_proteome_data)
print(paste("step1 - complete case num is", sum(complete.cases(ivn_proteome_data) == TRUE)))
write.table(ivn_proteome_data, "ivn_proteome_data.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Step 2: Load the covar data and do the linear regression to get the residuals
residuals_proteome_data <- data.frame(row.names = rownames(ivn_proteome_data))
covar <- read.table("covar.txt", sep = "\t", header = T)
covar$plate <- as.factor(covar$plate)
covar$sex <- as.factor(covar$sex)
for (i in colnames(ivn_proteome_data)) {
    x <- ivn_proteome_data[, i]
    tmpdata <- cbind(x, covar)
    res <- lm(x ~ ., data = tmpdata)
    tmprlt <- data.frame(residuals(res))
    colnames(tmprlt) <- i
    residuals_proteome_data <- cbind(residuals_proteome_data, tmprlt)
}
print(paste("step2 - complete case num is", sum(complete.cases(residuals_proteome_data) == TRUE)))
write.table(residuals_proteome_data, "residuals_proteome_data.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# Step 3: Run the pQTL mapping
cmd <- paste("plink2 --bfile CAS_genotype_data --pheno residuals_proteome_data.txt --no-psam-pheno --glm allow-no-covars --threads 12 --no-input-missing-phenotype --out pqtl_mapping")
system(cmd)
