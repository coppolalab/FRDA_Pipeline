library(caret)
library(glmnet)
library(biomaRt)
library(lumi)
library(WGCNA)
library(Metrics)
library(boot)

library(openxlsx)
library(Cairo)

library(broom)
library(rlist)
library(stringr)
library(magrittr)
library(tidyverse)

source("../../code/common_functions.R")

RegressionWorkbook <- function(rank.table, filename) {
    score.key <- colnames(rank.table) %>% str_detect("Coefficient") 
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

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

lumi.import <- ReadRDSgz("../WGCNA_GAA/save/rmcov.lumi.rda")

#SVM
svm.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, savePredictions = TRUE)
svm.fitcontrol.cv <- trainControl(method = "repeatedcv", verboseIter = TRUE, number = 5, repeats = 10, savePredictions = TRUE)
svm.fitcontrol.cv.10 <- trainControl(method = "repeatedcv", verboseIter = TRUE, number = 10, repeats = 10, savePredictions = TRUE)

set.seed(12345)
gaa.enet <- glmnet(t(exprs(lumi.import)), as.vector(lumi.import$FDS), "gaussian", nlambda = 10000, alpha = 0.5, type.gaussian = "naive")
gaa.enet.coef <- coef(gaa.enet)[-1,ncol(gaa.enet$beta)]
gaa.enet.df <- tibble(Symbol = names(gaa.enet.coef), Coefficient = signif(gaa.enet.coef, 3)) %>% arrange(desc(abs(Coefficient)))

lumi.top <- lumi.import[gaa.enet.coef != 0,]
set.seed(12345)
svm.train <- train(x = t(exprs(lumi.top)), y = as.vector(lumi.import$FDS), 
                   preProcess = c("center", "scale"), metric = "RMSE", 
                   method = "svmPoly", trControl = svm.fitcontrol) 
svm.train.cv <- train(x = t(exprs(lumi.top)), y = as.vector(lumi.import$FDS), 
                   preProcess = c("center", "scale"), metric = "RMSE", 
                   method = "svmPoly", trControl = svm.fitcontrol.cv) 
svm.train.cv.10 <- train(x = t(exprs(lumi.top)), y = as.vector(lumi.import$FDS), 
                   preProcess = c("center", "scale"), metric = "RMSE", 
                   method = "svmPoly", trControl = svm.fitcontrol.cv.10) 
svm.preds <- svm.train$pred

p <- ggplot(svm.preds, aes(x = obs, y = pred)) + 
    geom_point() + 
    facet_wrap( ~ Resample) + 
    stat_smooth(method = "loess") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Predicted Functional Stage") +
    xlab("Observed Functional Stage")

CairoPDF("svm.poly.preds", height = 25, width = 25)
print(p)
dev.off()

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(gaa.enet.df$Symbol), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")
enet.annot <- left_join(gaa.enet.df, bm.table) %>% 
    select(Symbol, Definition, Coefficient) %>% 
    as_tibble
RegressionWorkbook(enet.annot, "./elasticnet.table.xlsx")
