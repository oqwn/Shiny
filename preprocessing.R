library(GEOquery)
library(dplyr)
library(stringr)
library(sva)
library(Biobase)
library(caret)
library(cvTools)

# setwd("C:/Users/shrey/Desktop/DATA3888/appfiles/")
set.seed(1)

# ------------------------------------- BATCH CORRECTION -----------------------
# our primary dataset for binary model
gse26 = getGEO("GSE26578")
gse26 = gse26$GSE26578_series_matrix.txt.gz

# our secondary dataset for batch correction
gse14 = getGEO("GSE14328")
gse14 = gse14$GSE14328_series_matrix.txt.gz

# 411 common genes
length(intersect(fData(gse26)$`Gene`, fData(gse14)$`Gene Symbol`))

# prepping pData for combining
# introduce a new column called graft_status with binary outcomes for gse26
gse26_pData = pData(gse26) %>% 
  mutate(graft_status = case_when(str_detect(pData(gse26)$`disease status:ch1`, fixed("stable")) ~ "stable transplant",
                                  str_detect(pData(gse26)$`disease status:ch1`, fixed("acute")) ~ "acute rejection",
                                  str_detect(pData(gse26)$`disease status:ch1`, fixed("dysfunction")) ~ "acute rejection"))


# manually adding batch
gse26_pData$batch = rep("A", nrow(gse26_pData))
gse26_pData = gse26_pData %>% select(c(geo_accession, graft_status, batch))

# prepping fData data for combining
common_genes = intersect(fData(gse26)$`Gene`, fData(gse14)$`Gene Symbol`)
gse26_index = fData(gse26) %>% filter(`Gene` %in% common_genes)

# for fData, we take only the first uniquely occurring gene - this gives us the 411 common genes 
gse26_unique = gse26_index %>% group_by(`Gene`) %>% slice_head(n = 1)
gse26_fData = gse26_unique %>% mutate(ID = substring(ID, 1, 5))
gse26_fData = gse26_fData %>% select(c(ID, Gene, GB_ACC, GeneDes))
gse26_unique_index = gse26_fData$ID

# prepping exprs for combining
# taking aggregrate of the non-unique gene exprs
gse26_exprs = exprs(gse26) %>% as.data.frame()
gse26_exprs$full_geneID = rownames(gse26_exprs)
gse26_exprs_grouped = gse26_exprs %>% mutate(geneID = substring(full_geneID, 1, 5)) %>% select(-full_geneID)
gse26_exprs_grouped = gse26_exprs_grouped %>% group_by(geneID) %>% summarise_all(mean)
gse26_exprs_grouped = gse26_exprs_grouped %>% filter(geneID %in% gse26_unique_index) %>% as.data.frame()
rownames(gse26_exprs_grouped) = gse26_exprs_grouped$geneID
gse26_exprs_grouped = gse26_exprs_grouped %>% select(-geneID)

# for gse26, fData = gse26_fData, pData = gse26_pData, exprs = gse26_exprs_grouped

# same thing for gse14
# prepping pData for combining
gse14_pData = pData(gse14) %>%
  mutate(graft_status = case_when(str_detect(pData(gse14)$`characteristics_ch1`, fixed("stable")) ~ "stable transplant",
                                  str_detect(pData(gse14)$`characteristics_ch1`, fixed("acute")) ~ "acute rejection"))

# manually adding batch
gse14_pData$batch = rep("B", nrow(gse14_pData))
gse14_pData$geo_accession = rownames(gse14_pData)
gse14_pData = gse14_pData %>% select(c(geo_accession, graft_status, batch))

# prepping fData for combining
gse14_index = fData(gse14) %>% filter(`Gene Symbol` %in% common_genes)

# for fData, we take only the first uniquely occurring gene - this gives us the 411 common genes 
gse14_unique = gse14_index %>% group_by(`Gene Symbol`) %>% slice_head(n = 1)
gse14_fData = gse14_unique %>% select(c(ID, `Gene Symbol`, GB_ACC, `Gene Title`))
gse14_unique_index = gse14_unique$ID

# prepping exprs for combining 
# taking aggregate of the non-unique exprs
gse14_exprs = exprs(gse14) %>% as.data.frame()
gse14_exprs$geneID = rownames(gse14_exprs) 
gse14_exprs = gse14_exprs %>% filter(geneID %in% gse14_unique_index) %>% select(-geneID)

# for gse14, pData = gse14_pData, fData = gse14_fData, exprs = gse14_exprs

# combining
combined_pData = rbind(gse26_pData, gse14_pData)
combined_fData = gse26_fData %>% as.data.frame()
combined_exprs = cbind(gse26_exprs_grouped, gse14_exprs)

# building ExpressionSet object
combined_pData = AnnotatedDataFrame(data = combined_pData)
combined_gse = ExpressionSet(assayData = as.matrix(combined_exprs), phenoData = combined_pData, fData = combined_fData)
fData(combined_gse) = combined_fData

# batch correction
pheno = pData(combined_gse)
edata = exprs(combined_gse)
batch = pheno$batch
mod = model.matrix(~as.factor(graft_status), data = pheno)

# non-parametric adjustment, mean-only version
combat = sva::ComBat(dat = edata, batch = batch, mod = NULL, par.prior = FALSE, mean.only = TRUE)
corrected_gse = ExpressionSet(assayData = as.matrix(combat), phenoData = combined_pData)
fData(corrected_gse) = combined_fData

# removing the one NaN row from exprs
exprs_df = exprs(corrected_gse) %>% as.data.frame()
df_without_na = exprs_df[complete.cases(exprs_df), ]

fData_without_na = filter(fData(corrected_gse), ID != "ZH354")

corrected_gse = ExpressionSet(assayData = as.matrix(df_without_na), phenoData = combined_pData)
fData(corrected_gse) = fData_without_na


# ----------------------------------------- BINARY MODEL -----------------------------------------

folds = 5
repeats = 50

cv_binary = c()
se = c()
sp = c()
balanced_accs = c()

for (i in 1:repeats) {
  cvSets = cvFolds(nrow(pData(gse26)), folds)
  cv_each = c()
  sensitivities = c()
  specificities = c()
  balanced_acc = c()
  
  for (j in 1:folds) {
    test_id = cvSets$subsets[cvSets$which == j]
    X_train = t(exprs(corrected_gse))[-test_id, ]
    X_test = t(exprs(corrected_gse))[test_id, ]
    y_train = corrected_gse$graft_status[-test_id]
    y_test = corrected_gse$graft_status[test_id]
    
    design = model.matrix(~y_train)
    fit = lmFit(t(X_train), design)
    fit = eBayes(fit)
    top = topTable(fit, n = Inf)
    top = top %>% filter(adj.P.Val < 0.05)
    DE_genes = rownames(top)
    
    X_train = X_train[, DE_genes]
    X_test = X_test[, DE_genes]
    
    knn_fit = class::knn(train = X_train, test = X_test, cl = y_train, k = 1)
    cv_each[j] = mean(knn_fit == y_test)
    
    num_positives = sum(knn_fit == "acute rejection")
    num_negatives = sum(knn_fit == "stable transplant")
    
    sensitivities[j] = sum(knn_fit == "acute rejection" & y_test == "acute rejection")/num_positives
    specificities[j] = sum(knn_fit == "stable transplant" & y_test == "stable transplant")/num_negatives
    balanced_acc[j] = (sensitivities[j] + specificities[j])/2
  }
  cv_binary = append(cv_binary, mean(cv_each))
  se = append(se, mean(sensitivities))
  sp = append(sp, mean(specificities))
  balanced_accs = append(balanced_accs, mean(balanced_acc))
}

# store cv accuracies to make boxplots
cv_binary_metrics = cbind(cv_binary, se, sp, balanced_accs)
cv_binary_metrics %>% write.csv("cv_binary_metrics.csv")

# store the model input data
model_data = exprs(corrected_gse) %>% t() %>% as.data.frame()
knn_model = list(train_data = model_data, cl = corrected_gse$graft_status, k = 1)
saveRDS(knn_model, file = "binary_model.rds")

# ------------------------- IFTA SCORE ---------------------------------

gse25 = getGEO("GSE25902")
gse25 = gse25$GSE25902_series_matrix.txt.gz

# finding the top DE genes
ifta = gse25$`i-ifta grade:ch1` %>% as.factor()
design = model.matrix(~0 + ifta)
fit = lmFit(exprs(gse25), design)
fit = eBayes(fit)
tT = topTable(fit, number = Inf)
tT = tT %>% filter(adj.P.Val < 0.05) 
top_genes_ifta = rownames(tT[1:500, ])

# filter exprs to get top 500 genes
gse25_exprs = exprs(gse25) %>% as.data.frame()
gse25_exprs$geneID = rownames(exprs(gse25))
gse25_top500 = gse25_exprs %>% filter(geneID %in% top_genes_ifta)

# format exprs for SVM + CV
training_data = gse25_top500 %>% t() %>% as.data.frame()
training_data = training_data[1:120, ]

# convert all exprs values to numeric
training_data = apply(training_data, 2, as.numeric) %>% as.data.frame()

# normalisation
training_data = log2(training_data + 1)
training_data = apply(training_data, 2, function(x) {
  (x - mean(x))/sd(x)
})
training_data = data.frame(training_data, row.names = row.names(training_data))

training_data$ifta = as.numeric(ifta)

ctrl = trainControl(method = "repeatedcv", number = 5, repeats = 50)
model = train(ifta ~ ., data = training_data, method = "svmRadial", trControl = ctrl)
res = model$resample

model = model$finalModel

# save the model as a .RDS file
saveRDS(model, file = "ifta_model.rds")

# save MAE to a file so we can make boxplot
ifta_acc = data.frame(mae = res$MAE, rmse = res$RMSE, rsquared = res$Rsquared)
ifta_acc %>% write.csv("cv_ifta.csv")

# -------------------------------- CADI SCORE ---------------------------------------

# finding the top DE genes
cadi = gse25$`cadi score:ch1` %>% as.factor()
design = model.matrix(~0 + cadi)
fit = lmFit(exprs(gse25), design)
fit = eBayes(fit)
tT = topTable(fit, number = Inf)
tT = tT %>% filter(adj.P.Val < 0.05) 
top_genes_cadi = rownames(tT[1:500, ])

# filter exprs to get top 500 genes
gse25_exprs = exprs(gse25) %>% as.data.frame()
gse25_exprs$geneID = rownames(exprs(gse25))
gse25_top500 = gse25_exprs %>% filter(geneID %in% top_genes_cadi)

# format exprs for SVM + CV
training_data = gse25_top500 %>% t() %>% as.data.frame()
training_data = training_data[1:120, ]

# convert all exprs values to numeric
training_data = apply(training_data, 2, as.numeric) %>% as.data.frame()

# normalisation
training_data = log2(training_data + 1)
training_data = apply(training_data, 2, function(x) {
  (x - mean(x))/sd(x)
})
training_data = data.frame(training_data, row.names = row.names(training_data))

training_data$cadi = as.numeric(cadi)


ctrl = trainControl(method = "repeatedcv", number = 5, repeats = 50)
model = train(cadi ~ ., data = training_data, method = "svmRadial", trControl = ctrl)
res = model$resample

model = model$finalModel

# save the model as a .RDS file
saveRDS(model, file = "cadi_model.rds")

# save MAE to a file so we can make boxplot
cadi_acc = data.frame(mae = res$MAE, rmse = res$RMSE, rsquared = res$Rsquared)
cadi_acc %>% write.csv("cv_cadi.csv")

# ------------------------- LIST OF GENES NEEDED FOR PREDICTION -----------------------

binary_genes = fData(corrected_gse)
ifta_genes = fData(gse25) %>% filter(ID %in% top_genes_ifta)
cadi_genes = fData(gse25) %>% filter(ID %in% top_genes_cadi)

ifta_genes = ifta_genes %>% select(c(`ID`, `Gene Symbol`, `GB_ACC`, `Gene Ontology Molecular Function`)) %>% 
  rename(Gene = `Gene Symbol`, GeneDes = `Gene Ontology Molecular Function`)

cadi_genes = cadi_genes %>% select(c(ID, `Gene Symbol`, GB_ACC, `Gene Ontology Molecular Function`)) %>% 
  rename(Gene  = `Gene Symbol`, GeneDes = `Gene Ontology Molecular Function`)

binary_genes %>% write.csv("binary_genes.csv")
cadi_genes %>% write.csv("cadi_genes.csv")
ifta_genes %>% write.csv("ifta_genes.csv")

# --------------------------- PCA FOR BATCH CORRECTION -------------------------------------

before_pca = prcomp(t(exprs(combined_gse)))
before_toplot = data.frame(combined_gse$batch, pc1 = before_pca$x[, 1], pc2 = before_pca$x[, 2])
ggplot(before_toplot, aes(x = pc1, y = pc2, color = combined_gse$batch)) + geom_point(size = 3, shape = 18) + theme_minimal() +
  labs(x = "PC1", y = "PC2") + scale_color_discrete(name = "Batch") + theme(text = element_text(size = 20))

after_pca = prcomp(t(exprs(corrected_gse)))
after_toplot = data.frame(corrected_gse$batch, pc1 = after_pca$x[, 1], pc2 = after_pca$x[, 2])
ggplot(after_toplot, aes(x = pc1, y = pc2, color = corrected_gse$batch)) + geom_point(size = 3, shape = 18) + theme_minimal() +
  labs(x = "PC1", y = "PC2") + scale_color_discrete(name = "Batch") + theme(text = element_text(size = 20))

# --------------------------- ROBUSTNESS TEST -----------------------------------------------

# gse36 = getGEO("GSE36059")
# gse36 = gse36$GSE36059_series_matrix.txt.gz

# remove nephrectomies because that's not relevant
removed = which(pData(gse36)$`diagnosis (tcmr,abmr,mixed,non-rejecting):ch1` != "Nephrectomy")
filtered_pData = pData(gse36)[removed, ]
filtered_exprs = exprs(gse36) %>% as.data.frame()
filtered_exprs = filtered_exprs[, removed]

filtered_pData = AnnotatedDataFrame(data = filtered_pData)
filtered_gse = ExpressionSet(assayData = as.matrix(filtered_exprs), phenoData = filtered_pData)
fData(filtered_gse) = fData(gse36)

pData(filtered_gse) = pData(filtered_gse) %>% 
  mutate(graft_status = case_when(pData(filtered_gse)$`diagnosis (tcmr,abmr,mixed,non-rejecting):ch1` == "non-rejecting" ~ "stable transplant",
                                  pData(filtered_gse)$`diagnosis (tcmr,abmr,mixed,non-rejecting):ch1` == "ABMR" ~ "acute rejection",
                                  pData(filtered_gse)$`diagnosis (tcmr,abmr,mixed,non-rejecting):ch1` == "TCMR" ~ "acute rejection", 
                                  pData(filtered_gse)$`diagnosis (tcmr,abmr,mixed,non-rejecting):ch1` == "MIXED" ~ "acute rejection"))

# finding the top DE genes
design = model.matrix(~0 + pData(filtered_gse)$graft_status)
fit = lmFit(exprs(filtered_gse), design)
fit = eBayes(fit)
tT = topTable(fit, number = Inf)
tT = tT %>% filter(adj.P.Val < 0.05) 
top_genes_gse36 = rownames(tT[1:410, ])

common_genes_gse36 = fData(filtered_gse) %>% filter(ID %in% top_genes_gse36)

# only 7 common genes
intersect(common_genes_gse36$`Gene Symbol`, binary_genes$Gene)

robust_test = exprs(filtered_gse) %>% as.data.frame()
robust_test$ID = rownames(robust_test)
robust_test = robust_test %>% filter(ID %in% top_genes_gse36) %>% select(-ID) %>% t()
robust_output = class::knn(train = knn_model$train_data, test = robust_test, cl = knn_model$cl, k = knn_model$k)

robust_acc = mean(robust_output == filtered_gse$graft_status)

num_positives = sum(robust_output == "acute rejection")
num_negatives = sum(robust_output == "stable transplant")

sens = sum(robust_output == "acute rejection" & filtered_gse$graft_status == "acute rejection")/num_positives
spec = sum(robust_output == "stable transplant" & filtered_gse$graft_status == "stable transplant")/num_negatives
balanced_acc_robust = (sens + spec)/2

# --------------------------------------- TEMPLATE FILES FOR EACH MODEL WITH GENE SYMBOLS ----------------------
col_names = fData_without_na$Gene

# Sample list of column names
columnNames <- c(col_names)

# Create a data frame with empty rows
dataFrame <- data.frame(matrix(ncol = length(columnNames), nrow = 0))
colnames(dataFrame) <- columnNames

# Save the data frame to a CSV file
write.csv(dataFrame, file = "column_names.csv", row.names = FALSE)

data <- read.csv("column_names.csv")

# Create a blank row with the same number of columns as the original data
blank_row <- data.frame(matrix(ncol = ncol(data), nrow = 1))
colnames(blank_row) <- colnames(data)

# Append the blank row to the data
data <- rbind(data, blank_row)

# Save the updated data back to the CSV file, overwriting the original content
write.csv(data, file = "column_names.csv", row.names = FALSE)

# setwd("C:/Users/shrey/Desktop/DATA3888/projectfiles/")
# extracting first gene symbol only
ifta_genes$Gene = word(ifta_genes$Gene, 1)
ifta_genes$Gene

col_names_ifta = ifta_genes$Gene

# Sample list of column names
columnNames <- c(col_names_ifta)

# Create a data frame with empty rows
dataFrame <- data.frame(matrix(ncol = length(columnNames), nrow = 0))
colnames(dataFrame) <- columnNames

# Save the data frame to a CSV file
write.csv(dataFrame, file = "column_names_ifta.csv", row.names = FALSE)

# adding blank row -----
data <- read.csv("column_names_ifta.csv")

# Create a blank row with the same number of columns as the original data
blank_row <- data.frame(matrix(ncol = ncol(data), nrow = 1))
colnames(blank_row) <- colnames(data)

# Append the blank row to the data
data <- rbind(data, blank_row)

# Save the updated data back to the CSV file, overwriting the original content
write.csv(data, file = "column_names_ifta.csv", row.names = FALSE)

setwd("C:/Users/shrey/Desktop/DATA3888/projectfiles/")
# extracting first gene symbol only
cadi_genes$Gene = word(cadi_genes$Gene, 1)
cadi_genes$Gene

col_names_cadi = cadi_genes$Gene

# Sample list of column names
columnNames <- c(col_names_cadi)

# Create a data frame with empty rows
dataFrame <- data.frame(matrix(ncol = length(columnNames), nrow = 0))
colnames(dataFrame) <- columnNames

# Save the data frame to a CSV file
write.csv(dataFrame, file = "column_names_cadi.csv", row.names = FALSE)

# adding blank row -----
data <- read.csv("column_names_cadi.csv")

# Create a blank row with the same number of columns as the original data
blank_row <- data.frame(matrix(ncol = ncol(data), nrow = 1))
colnames(blank_row) <- colnames(data)

# Append the blank row to the data
data <- rbind(data, blank_row)

# Save the updated data back to the CSV file, overwriting the original content
write.csv(data, file = "column_names_cadi.csv", row.names = FALSE)
