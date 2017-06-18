library(limma)
library(annotate)
library(lumiHumanAll.db)
library(glmnet)
library(caret)
library(Biobase)
library(lumi)
library(biomaRt)

library(Cairo)
library(openxlsx)

library(rlist)
library(stringr)
library(magrittr)
library(tidyverse)

ClassifierWorkbook <- function(rank.table, filename) {
    score.key <- colnames(rank.table) %>% str_detect("Coefficient.") 
    score.cols <- which(score.key)

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = rank.table)
    conditionalFormatting(wb, 1, cols = score.cols, rows = 1:nrow(rank.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:ncol(rank.table), widths = "auto")
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 11)
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

source("../common_functions.R")
lumi.import <- ReadRDSgz("../baseline_lumi/save/rmcov.lumi.rda")

#Create trainControl object
svm.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, savePredictions = TRUE)
svm.fitcontrol.cv <- trainControl(method = "repeatedcv", verboseIter = TRUE, number = 5, repeats = 10, savePredictions = TRUE)
svm.fitcontrol.cv.10 <- trainControl(method = "repeatedcv", verboseIter = TRUE, number = 10, repeats = 10, savePredictions = TRUE)

set.seed(12345)
status.lasso <- glmnet(t(exprs(lumi.import)), lumi.import$Status, "multinomial", type.multinomial = "grouped", nlambda = 1000, alpha = 1)
status.lasso.patient <- coef(status.lasso)$Patient[-1,1000]
status.lasso.carrier <- coef(status.lasso)$Carrier[-1,1000]
status.lasso.control <- coef(status.lasso)$Control[-1,1000]
status.lasso.df <- tibble(Symbol = names(status.lasso.coef), 
                          Coefficient.Patient = status.lasso.patient,
                          Coefficient.Carrier = status.lasso.carrier,
                          Coefficient.Control = status.lasso.control) %>% 
                   arrange(desc(abs(Coefficient.Patient)))

lumi.import.top <- lumi.import[status.lasso.patient !=0, ]
set.seed(12345)
svm.final <- train(x = t(exprs(lumi.import.top)), y = droplevels(factor(lumi.import.top$Status)), 
                   preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol)
svm.final.cv <- train(x = t(exprs(lumi.import.top)), y = droplevels(factor(lumi.import.top$Status)), 
                   preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol.cv)
svm.final.cv.10 <- train(x = t(exprs(lumi.import.top)), y = droplevels(factor(lumi.import.top$Status)), 
                   preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol.cv.10)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', 
                  values = as.character(status.lasso.df$Symbol), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")
lasso.annot <- left_join(status.lasso.df, bm.table) %>% 
    dplyr::select(Symbol, Definition, dplyr::contains("Coefficient")) 
ClassifierWorkbook(lasso.annot, "lasso.table.xlsx")

svm.preds <- svm.final$pred
p <- ggplot(svm.preds, aes(x = obs, fill = pred)) + 
    geom_bar(color = "black") + 
    facet_wrap( ~ Resample) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Predicted Diagnosis") +
    xlab("Observed Diagnosis")

CairoPDF("svm.linear.preds", height = 25, width = 15)
print(p)
dev.off()

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort
