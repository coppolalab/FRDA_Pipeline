#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(igraph)

#Reading and writing tables
library(readr)
library(openxlsx)

#Classification
library(MASS)
library(randomForest)
library(rpart)
library(supclust)
library(class)
library(e1071)
library(pamr)
library(glmnet)
library(sda)
library(randomGLM)
library(rknn)
library(crossval)
library(WGCNA)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(parallel)
library(Biobase)
library(lumi)

#Functional programming
library(magrittr)
library(purrr)
library(functional)

saveRDS.gz <- function(object,file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

readRDS.gz <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

enrichr.wkbk <- function(subindex, full.df, new.fullpath)
{
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    dir.create(new.fullpath, showWarnings = FALSE, recursive = TRUE)
    filename = paste(file.path(new.fullpath, subindex), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

lumi.import <- readRDS.gz("../baseline_lumi/save/rmcov.collapse.rda")
#status <- model.matrix( ~ 0 + lumi.import$Status)
#colnames(status) <- c("Carrier", "Control", "Patient")

#fdata <- fData(lumi.import)
#lumi.import <- lumi.import[!is.na(fdata$SYMBOL),]
#fdata <- fData(lumi.import)
#pdata <- pData(lumi.import)
#expr.collapse <- collapseRows(exprs(lumi.import), factor(fdata$SYMBOL), rownames(lumi.import), method = "function", methodFunction = colMeans)$datETcollapsed

lumi.pca <- lumi.import[,lumi.import$Status == "Patient" | lumi.import$Status == "Carrier"]
lumi.pco <- lumi.import[,lumi.import$Status == "Patient" | lumi.import$Status == "Control"]
lumi.cc <- lumi.import[,lumi.import$Status == "Control" | lumi.import$Status == "Carrier"]
lumi.pca$Status %<>% droplevels
lumi.pco$Status %<>% droplevels
lumi.cc$Status %<>% droplevels

#Random Forest
rf.pca <- randomForest(t(exprs(lumi.pca)), droplevels(lumi.pca$Status), ntree = 1000)
#rfcv.pca <- rfcv(t(exprs(lumi.pca)), droplevels(lumi.pca$Status), ntree = 1000)
saveRDS.gz(rf.pca, "./save/rf.pca.rda")
#saveRDS.gz(rfcv.pca, "./save/rfcv.pca.rda")
rf.pca.genes <- data.frame("nuID" = rownames(rf.pca$importance), "Importance" = as.vector(rf.pca$importance)) %>% join(fdata) %>% arrange(desc(Importance))
write.xlsx(rf.pca.genes, "rf.pca.genes.xlsx")

rf.pco <- randomForest(t(exprs(lumi.pco)), droplevels(lumi.pco$Status), ntree = 1000)
saveRDS.gz(rf.pco, "./save/rf.pco.rda")
rf.pco.genes <- data.frame("nuID" = rownames(rf.pco$importance), "Importance" = as.vector(rf.pco$importance)) %>% join(fdata) %>% arrange(desc(Importance))
write.xlsx(rf.pco.genes, "rf.pco.genes.xlsx")

rf.cc <- randomForest(t(exprs(lumi.cc)), droplevels(lumi.cc$Status), ntree = 1000)
saveRDS.gz(rf.cc, "./save/rf.cc.rda")
rf.cc.genes <- data.frame("nuID" = rownames(rf.cc$importance), "Importance" = as.vector(rf.cc$importance)) %>% join(fdata) %>% arrange(desc(Importance))
write.xlsx(rf.cc.genes, "rf.cc.genes.xlsx")

rfcv.pca <- readRDS.gz("./save/rfcv.pca.rda")

rfcv.pco <- rfcv(t(exprs(lumi.pco)), droplevels(lumi.pco$Status), ntree = 1000)
saveRDS.gz(rfcv.pco, "./save/rfcv.pco.rda")
rfcv.pco <- readRDS.gz("./save/rfcv.pco.rda")

rfcv.cc <- rfcv(t(exprs(lumi.cc)), droplevels(lumi.cc$Status), ntree = 1000)
saveRDS.gz(rfcv.cc, "./save/rfcv.cc.rda")
rfcv.cc <- readRDS.gz("./save/rfcv.cc.rda")

rf.bigmtry.pca <- randomForest(t(exprs(lumi.pca)), droplevels(lumi.pca$Status), ntree = 1000, mtry = ncol(t(exprs(lumi.pca))))
saveRDS.gz(rf.bigmtry.pca, "./save/rf.bigmtry.pca.rda")
rf.bigmtry.pca.genes <- data.frame("nuID" = rownames(rf.bigmtry.pca$importance), "Importance" = as.vector(rf.bigmtry.pca$importance)) %>% join(fdata) %>% arrange(desc(Importance))
write.xlsx(rf.bigmtry.pca.genes, "rf.bigmtry.pca.genes.xlsx")

rf.bigmtry.pco <- randomForest(t(exprs(lumi.pco)), droplevels(lumi.pco$Status), ntree = 1000, mtry = ncol(t(exprs(lumi.pco))))
saveRDS.gz(rf.bigmtry.pco, "./save/rf.bigmtry.pco.rda")
rf.bigmtry.pco.genes <- data.frame("nuID" = rownames(rf.bigmtry.pco$importance), "Importance" = as.vector(rf.bigmtry.pco$importance)) %>% join(fdata) %>% arrange(desc(Importance))
write.xlsx(rf.bigmtry.pco.genes, "rf.bigmtry.pco.genes.xlsx")

rf.bigmtry.cc <- randomForest(t(exprs(lumi.cc)), droplevels(lumi.cc$Status), ntree = 1000, mtry = ncol(t(exprs(lumi.cc))))
saveRDS.gz(rf.bigmtry.cc, "./save/rf.bigmtry.cc.rda")
rf.bigmtry.cc.genes <- data.frame("nuID" = rownames(rf.bigmtry.cc$importance), "Importance" = as.vector(rf.bigmtry.cc$importance)) %>% join(fdata) %>% arrange(desc(Importance))
write.xlsx(rf.bigmtry.cc.genes, "rf.bigmtry.cc.genes.xlsx")

#Linear Discriminant Analysis
lda.predfunction <- function(train.x, train.y, test.x, test.y, negative)
{
   sda.fit <- sda(train.x, train.y)
   ynew <- predict(sda.fit, test.x)$class
   
   sda.confusion <- confusionMatrix(test.y, ynew, negative = negative)
   sda.diagnostic <- diagnosticErrors(sda.confusion)
}

dlda.predfunction <- function(train.x, train.y, test.x, test.y, negative)
{
   sda.fit <- sda(train.x, train.y, diagonal = TRUE)
   ynew <- predict(sda.fit, test.x)$class
   
   sda.confusion <- confusionMatrix(test.y, ynew, negative = negative)
   sda.diagnostic <- diagnosticErrors(sda.confusion)
}

lda.pca.cv <- crossval(lda.predfunction, t(exprs(lumi.pca)), lumi.pca$Status, K = 3, B = 1, negative = "Carrier")
lda.pca <- sda.ranking(t(expr.pca), lumi.pca$Status)
saveRDS.gz(lda.pca, "./save/lda.pca.rda")
lda.pca.df <- data.frame(rownames(lda.pca), Ranking = lda.pca[,"score"]) %>% arrange(desc(Ranking))
write.xlsx(lda.pca.df, "lda.pca.xlsx")

lda.pco.cv <- crossval(lda.predfunction, t(exprs(lumi.pco)), lumi.pco$Status, K = 3, B = 1, negative = "Control")
lda.pco <- sda.ranking(t(expr.pco), lumi.pco$Status)
saveRDS.gz(lda.pco, "./save/lda.pco.rda")
lda.pco.df <- data.frame(rownames(lda.pco), Ranking = lda.pco[,"score"]) %>% arrange(desc(Ranking))
write.xlsx(lda.pco.df, "lda.pco.xlsx")

lda.cc.cv <- crossval(lda.predfunction, t(exprs(lumi.cc)), lumi.cc$Status, K = 3, B = 1, negative = "Control")
lda.cc <- sda.ranking(t(expr.cc), lumi.cc$Status)
saveRDS.gz(lda.cc, "./save/lda.cc.rda")
lda.cc.df <- data.frame(rownames(lda.cc), Ranking = lda.cc[,"score"]) %>% arrange(desc(Ranking))
write.xlsx(lda.cc.df, "lda.cc.xlsx")

dlda.pca.cv <- crossval(dlda.predfunction, t(exprs(lumi.pca)), lumi.pca$Status, K = 3, B = 1, negative = "Carrier")
dlda.pca <- sda.ranking(t(expr.pca), lumi.pca$Status)
saveRDS.gz(dlda.pca, "./save/dlda.pca.rda")
dlda.pca.df <- data.frame(rownames(dlda.pca), Ranking = dlda.pca[,"score"]) %>% arrange(desc(Ranking))
write.xlsx(dlda.pca.df, "dlda.pca.xlsx")

dlda.pco.cv <- crossval(dlda.predfunction, t(exprs(lumi.pco)), lumi.pco$Status, K = 3, B = 1, negative = "Control")
dlda.pco <- sda.ranking(t(expr.pco), lumi.pco$Status)
saveRDS.gz(dlda.pco, "./save/dlda.pco.rda")
dlda.pco.df <- data.frame(rownames(dlda.pco), Ranking = dlda.pco[,"score"]) %>% arrange(desc(Ranking))
write.xlsx(dlda.pco.df, "dlda.pco.xlsx")

dlda.cc.cv <- crossval(dlda.predfunction, t(exprs(lumi.cc)), lumi.cc$Status, K = 3, B = 1, negative = "Control")
dlda.cc <- sda.ranking(t(expr.cc), lumi.cc$Status)
saveRDS.gz(dlda.cc, "./save/dlda.cc.rda")
dlda.cc.df <- data.frame(rownames(dlda.cc), Ranking = dlda.cc[,"score"]) %>% arrange(desc(Ranking))
write.xlsx(dlda.cc.df, "dlda.cc.xlsx")

#Random GLM
get.thinned.accuracy <- function(threshold, rglm.trained, dataset.test)
{
    print(threshold)
    thinned.object <- thinRandomGLM(rglm.object, threshold)

    thinned.correct <- which(thinned.object$y.original == thinned.object$predictedOOB) %>% length
    thinned.accuracy <- thinned.correct / length(thinned.object$y.original)
    thinned.features <- which(thinned.object$timesSelectedByForwardRegression > 0) %>% length
    return(c(threshold, thinned.accuracy, thinned.features))
}

rglm.fold <- function(train.x, train.y, test.x, test.y, negative, nThreads)
{
    #print(dim(train.x))
    #print(train.y)
    rglm.trained <- randomGLM(train.x, train.y, nThreads = nThreads)
    rglm.predicted <- predict(rglm.trained, newdata = test.x, type = "class") 
    
    accuracy.matrix <- confusionMatrix(test.y, rglm.predicted, negative = negative)
    accuracy.return <- diagnosticErrors(accuracy.matrix)
    
    return(accuracy.return)
}

nfold <- 3
#rglm.pca <- randomGLM(t(exprs(lumi.pca)), droplevels(lumi.pca$Status), nThreads = 7)
rglm.fold.pca <- crossval(rglm.fold, t(exprs(lumi.pca)), droplevels(lumi.pca$Status), K = 3, B = 1, negative = "Carrier", nThreads = 7)
saveRDS.gz(rglm.fold.pca, "./save/rglm.fold.pca.rda")
rglm.thresholds <- seq(1, 13, 1)
rglm.pca.accuracies <- map(rglm.thresholds, get.thinned.accuracy, rglm.pca) 
rglm.pca.table <- reduce(rglm.pca.accuracies, rbind) %>% data.frame
colnames(rglm.pca.table) <- c("Index", "Accuracy", "Num.Genes")

thinned.pca <- thinRandomGLM(rglm.pca, 1)
saveRDS.gz(thinned.pca, "./save/thinned.pca.rda")
pca.rglm <- thinned.pca$timesSelectedByForwardRegression 
pca.rglm.df <- data.frame(Symbol = featureNames(lumi.pca), times.selected = as.vector(t(pca.rglm))) %>% filter(times.selected > 0) %>% arrange(desc(times.selected)) 
write.xlsx(pca.rglm.df, "pca.rglm.xlsx")

#rglm.pco <- randomGLM(t(exprs(lumi.pco)), droplevels(lumi.pco$Status), nThreads = 7)
rglm.fold.pco <- crossval(rglm.fold, t(exprs(lumi.pco)), droplevels(lumi.pco$Status), K = 3, B = 1, negative = "Control", nThreads = 7)
saveRDS.gz(rglm.fold.pco, "./save/rglm.fold.pco.rda")
rglm.pco.accuracies <- map(rglm.thresholds, get.thinned.accuracy, rglm.pco) 
rglm.pco.table <- reduce(rglm.pco.accuracies, rbind) %>% data.frame
colnames(rglm.pco.table) <- c("Index", "Accuracy", "Num.Genes")

thinned.pco <- thinRandomGLM(rglm.pco, 1)
saveRDS.gz(thinned.pco, "./save/thinned.pco.rda")
pco.rglm <- thinned.pco$timesSelectedByForwardRegression 
pco.rglm.df <- data.frame(Symbol = featureNames(lumi.pco), times.selected = as.vector(t(pco.rglm))) %>% filter(times.selected > 0) %>% arrange(desc(times.selected)) 
write.xlsx(pco.rglm.df, "pco.rglm.xlsx")

#rglm.cc <- randomGLM(t(exprs(lumi.cc)), droplevels(lumi.cc$Status), nThreads = 7)
rglm.fold.cc <- crossval(rglm.fold, t(exprs(lumi.cc)), droplevels(lumi.cc$Status), K = 3, B = 1, negative = "Control", nThreads = 7)
saveRDS.gz(rglm.fold.cc, "./save/rglm.fold.cc.rda")
rglm.cc.accuracies <- map(seq(1, 13, 1), get.thinned.accuracy, rglm.cc) 
rglm.cc.table <- reduce(rglm.cc.accuracies, rbind) %>% data.frame
colnames(rglm.cc.table) <- c("Index", "Accuracy", "Num.Genes")

thinned.cc <- thinRandomGLM(rglm.cc, 1)
saveRDS.gz(thinned.cc, "./save/thinned.cc.rda")
cc.rglm <- thinned.cc$timesSelectedByForwardRegression 
cc.rglm.df <- data.frame(Symbol = featureNames(lumi.cc), times.selected = as.vector(t(cc.rglm))) %>% filter(times.selected > 0) %>% arrange(desc(times.selected)) 
write.xlsx(cc.rglm.df, "cc.rglm.xlsx")

#Support Vector Machine
#Compare to 10 fold at some point!
list.nums.init <- vector()
recursive.multiply <- function(list.nums, current.value)
{
   new.val = round(0.8 * current.value)
   if (new.val > 2)
   {
       list.nums <- c(list.nums, new.val)
       recursive.multiply(list.nums, new.val)
   }
   else
   {
       return(list.nums)
   }
}

list.nums <- recursive.multiply(list.nums.init, nrow(lumi.import))
list.nums <- c(nrow(lumi.import), list.nums)
nums.1000 <- list.nums[list.nums < 1000]

get.minerror <- function(svm.list)
{
    min(svm.list$error)
}
extract.svmlist <- function(featsweep.object)
{
    svm.list <- featsweep.object$svm.list
    map_dbl(svm.list, get.minerror) %>% mean
}
source("../../code/msvmRFE.R")

#svm.pca <- svm(x = t(expr.pca), y = droplevels(lumi.pca$Status), type = "C-classification", kernel = "linear", cross = 5, cost = 1, cachesize = 16000)
#saveRDS.gz(svm.pca, "./save/svm.pca.rda")
#svm.pca <- svmRFE(svm.pca.input, k = 10, halve.above = 100)

nfold <- 3
svm.pca.input <- cbind(as.integer(droplevels(lumi.pca$Status)), t(exprs(lumi.pca)))
nrows.pca <- nrow(svm.pca.input)
folds.pca.index <- rep(1:nfold, len=nrows.pca)[sample(nrows.pca)]
folds.pca.rows <- lapply(1:nfold, function(x) which(folds.pca.index == x))
svm.pca.results <- lapply(folds.pca.rows, svmRFE.wrap, svm.pca.input, k=3, halve.above=3)
saveRDS.gz(svm.pca.results, "svm.pca.results.rda")
svm.pca.features <- WriteFeatures(svm.pca.results, svm.pca.input, save=F)
colnames(svm.pca.features) <- "Symbol"
#svm.pca.genes <- join(svm.pca.features, fdata)
write.xlsx(svm.pca.features, "svm.pca.features.xlsx")

#svm.pco <- svm(x = t(expr.pco), y = droplevels(lumi.pco$Status), type = "C-classification", kernel = "linear", cross = 5, cost = 1, cachesize = 16000)
#saveRDS.gz(svm.pco, "./save/svm.pco.rda")

svm.pco.input <- cbind(as.integer(droplevels(lumi.pco$Status)), t(exprs(lumi.pco)))
nrows.pco <- nrow(svm.pco.input)
folds.pco.index <- rep(1:nfold, len=nrows.pco)[sample(nrows.pco)]
folds.pco.rows <- lapply(1:nfold, function(x) which(folds.pco.index == x))
svm.pco.results <- lapply(folds.pco.rows, svmRFE.wrap, svm.pco.input, k=3, halve.above=3)
svm.pco.features <- WriteFeatures(svm.pco.results, svm.pco.input, save=F)
colnames(svm.pco.features) <- "Symbol"
#svm.pco.genes <- join(svm.pco.features, fdata)
write.xlsx(svm.pco.features, "svm.pco.features.xlsx")

#svm.cc <- svm(x = t(expr.cc), y = droplevels(lumi.cc$Status), type = "C-classification", kernel = "linear", cross = 5, cost = 1, cachesize = 16000)
#saveRDS.gz(svm.cc, "./save/svm.cc.rda")

svm.cc.input <- cbind(as.integer(droplevels(lumi.cc$Status)), t(exprs(lumi.cc)))
nrows.cc <- nrow(svm.cc.input)
folds.cc.index <- rep(1:nfold, len=nrows.cc)[sample(nrows.cc)]
folds.cc.rows <- lapply(1:nfold, function(x) which(folds.cc.index == x))
svm.cc.results <- lapply(folds.cc.rows, svmRFE.wrap, svm.cc.input, k=3, halve.above=3)
svm.cc.features <- WriteFeatures(svm.cc.results, svm.cc.input, save=F)
colnames(svm.cc.features) <- "Symbol"
#svm.cc.genes <- join(svm.cc.features, fdata)
write.xlsx(svm.cc.features, "svm.cc.features.xlsx")

featsweep.pca <- lapply(list.nums, FeatSweep.wrap, svm.pca.results, svm.pca.input)
svm.pca.accuracies <- map_dbl(featsweep.pca, extract.svmlist)
svm.pca.accuracies.df <- data.frame(Num.Genes = list.nums, Accuracy = (1 - svm.pca.accuracies)) 
svm.pca.accuracies.df$comparison <- "Patient vs Carrier"
saveRDS.gz(featsweep.pca, "./save/featsweep.pca.rda")
saveRDS.gz(svm.pca.accuracies.df, "./save/svm.pca.accuracies.df.rda")

featsweep.pco <- lapply(list.nums, FeatSweep.wrap, svm.pco.results, svm.pco.input)
svm.pco.accuracies <- map_dbl(featsweep.pco, extract.svmlist)
svm.pco.accuracies.df <- data.frame(Num.Genes = list.nums, Accuracy = (1 - svm.pco.accuracies))
svm.pco.accuracies.df$comparison <- "Patient vs Control"
saveRDS.gz(featsweep.pco, "./save/featsweep.pco.rda")

featsweep.cc <- lapply(list.nums, FeatSweep.wrap, svm.cc.results, svm.cc.input)
svm.cc.accuracies <- map_dbl(featsweep.cc, extract.svmlist)
svm.cc.accuracies.df <- data.frame(Num.Genes = list.nums, Accuracy = (1 - svm.cc.accuracies))
svm.cc.accuracies.df$comparison <- "Carrier vs Control"
saveRDS.gz(featsweep.cc, "./save/featsweep.cc.rda")

svm.rfe.accuracies <- rbind(svm.pca.accuracies.df, svm.pco.accuracies.df, svm.cc.accuracies.df)
svm.rfe.accuracies$Num.Genes <- factor(svm.rfe.accuracies$Num.Genes, unique(svm.rfe.accuracies$Num.Genes, decreasing = TRUE)) 
rfe.plot <- ggplot(svm.rfe.accuracies, aes(x = Num.Genes, y = Accuracy, col = comparison, group = comparison)) + geom_line() + geom_point()
rfe.plot <- rfe.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ ylim(0, 1) + xlab("Number of Genes")
CairoPDF("svm.rfe.accuracies", width = 10, height = 5)
print(rfe.plot)
dev.off()

source("../../code/GO/enrichr.R")

svm.cc.submit <- svm.cc.features[1:500,]
colnames(svm.cc.submit) <- c("Symbol", "X1", "X2")
svm.pca.submit <- svm.pca.features[1:500,]
colnames(svm.pca.submit) <- c("Symbol", "X1", "X2")
svm.pco.submit <- svm.pco.features[1:500,]
colnames(svm.pco.submit) <- c("Symbol", "X1", "X2")

get.stringdb(svm.cc.submit, "svm.cc", "./cc")
get.stringdb(svm.pca.submit, "svm.pca", "./pca")
get.stringdb(svm.pco.submit, "svm.pco", "./pco")

enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 
svm.cc.enrichr <- map(enrichr.terms, get.enrichrdata, svm.cc.submit, FALSE)
names(svm.cc.enrichr) <- enrichr.terms
map(names(svm.cc.enrichr), enrichr.wkbk, svm.cc.enrichr, "./enrichr/cc")

svm.pca.enrichr <- map(enrichr.terms, get.enrichrdata, svm.pca.submit, FALSE)
names(svm.pca.enrichr) <- enrichr.terms
map(names(svm.pca.enrichr), enrichr.wkbk, svm.pca.enrichr, "./enrichr/pca")

svm.pco.enrichr <- map(enrichr.terms, get.enrichrdata, svm.pco.submit, FALSE)
names(svm.pco.enrichr) <- enrichr.terms
map(names(svm.pco.enrichr), enrichr.wkbk, svm.pco.enrichr, "./enrichr/pco")

#Shrunken Centroid
get.genelists <- function(threshold, ...)
{
    sc.genes <- pamr.listgenes(..., threshold)
    colnames(sc.genes)[1] <- "Symbol"
    return(sc.genes)
}

pca.data <- list(x = exprs(lumi.pca), y = lumi.pca$Status, geneid = featureNames(lumi.pca), genenames = featureNames(lumi.pca)) 
sc.pca <- pamr.train(pca.data)
sc.pca.cv <- pamr.cv(fit = sc.pca, data = pca.data, nfold = 3)
num.genes.pca <- cbind(sc.pca.cv$threshold, sc.pca.cv$size, (1 - sc.pca.cv$error)) %>% data.frame
colnames(num.genes.pca) <- c("Threshold", "Num.Genes", "Accuracy")
pca.500 <- filter(num.genes.pca, Num.Genes <= 500 & Num.Genes > 0)
pca.500.genes <- map(pca.500$Threshold, get.genelists, fit = sc.pca, data = pca.data)
write.xlsx(pca.500.genes[[which.max(pca.500$Accuracy)]], "sc.pca.xlsx")
saveRDS.gz(sc.pca, "./save/sc.pca.rda")

pco.data <- list(x = exprs(lumi.pco), y = lumi.pco$Status, geneid = featureNames(lumi.pco), genenames = featureNames(lumi.pco)) 
sc.pco <- pamr.train(pco.data)
sc.pco.cv <- pamr.cv(fit = sc.pco, data = pco.data, nfold = 3)
num.genes.pco <- cbind(sc.pco.cv$threshold, sc.pco.cv$size, (1 - sc.pco.cv$error)) %>% data.frame
colnames(num.genes.pco) <- c("Threshold", "Num.Genes", "Accuracy")
pco.500 <- filter(num.genes.pco, Num.Genes <= 500 & Num.Genes > 0)
pco.500.genes <- map(pco.500$Threshold, get.genelists, fit = sc.pco, data = pco.data)
write.xlsx(pco.500.genes[[which.max(pco.500$Accuracy)]], "sc.pco.xlsx")
saveRDS.gz(sc.pco, "./save/sc.pco.rda")

cc.data <- list(x = exprs(lumi.cc), y = lumi.cc$Status, geneid = featureNames(lumi.cc), genenames = featureNames(lumi.cc)) 
sc.cc <- pamr.train(cc.data)
sc.cc.cv <- pamr.cv(fit = sc.cc, data = cc.data, nfold = 3)
num.genes.cc <- cbind(sc.cc.cv$threshold, sc.cc.cv$size, (1 - sc.cc.cv$error)) %>% data.frame
colnames(num.genes.cc) <- c("Threshold", "Num.Genes", "Accuracy")
cc.500 <- filter(num.genes.cc, Num.Genes <= 500 & Num.Genes > 0)
cc.500.genes <- map(cc.500$Threshold, get.genelists, fit = sc.cc, data = cc.data)
write.xlsx(cc.500.genes[[which.max(cc.500$Accuracy)]], "sc.cc.xlsx")
saveRDS.gz(sc.cc, "./save/sc.cc.rda")

#Ridge Regression
ridge.pca <- cv.glmnet(x = t(exprs(lumi.pca)), y = droplevels(lumi.pca$Status), type.measure = "auc", nfolds = 3, family = "binomial", alpha = 0)
ridge.pca.betas <- ridge.pca$glmnet.fit$beta
ridge.pca.betas <- ridge.pca.betas[,ncol(ridge.pca.betas)] %>% abs
ridge.pca.df <- data.frame(Symbol = names(ridge.pca.betas), Importance = ridge.pca.betas) %>% arrange(desc(Importance))
write.xlsx(ridge.pca.df, "./ridge.pca.xlsx")
saveRDS.gz(ridge.pca, "./save/ridge.pca.rda")

ridge.pco <- cv.glmnet(x = t(exprs(lumi.pco)), y = droplevels(lumi.pco$Status), family = "binomial", type.measure = "auc", nfolds = 3, alpha = 0)
ridge.pco.betas <- ridge.pco$glmnet.fit$beta
ridge.pco.betas <- ridge.pco.betas[,ncol(ridge.pco.betas)] %>% abs
ridge.pco.df <- data.frame(Symbol = names(ridge.pco.betas), Importance = ridge.pco.betas) %>% arrange(desc(Importance))
write.xlsx(ridge.pco.df, "./ridge.pco.xlsx")
saveRDS.gz(ridge.pco, "./save/ridge.pco.rda")

ridge.cc <- cv.glmnet(x = t(exprs(lumi.cc)), y = droplevels(lumi.cc$Status), family = "binomial", type.measure = "auc", nfolds = 3, alpha = 0)
ridge.cc.betas <- ridge.cc$glmnet.fit$beta
ridge.cc.betas <- ridge.cc.betas[,ncol(ridge.cc.betas)] %>% abs
ridge.cc.df <- data.frame(Symbol = names(ridge.cc.betas), Importance = ridge.cc.betas) %>% arrange(desc(Importance))
write.xlsx(ridge.cc.df, "./ridge.cc.xlsx")
saveRDS.gz(ridge.cc, "./save/ridge.cc.rda")

#Lasso Regression
lasso.pca <- cv.glmnet(x = t(exprs(lumi.pca)), y = droplevels(lumi.pca$Status), family = "binomial", type.measure = "auc", nfolds = 3, alpha = 1)
lasso.pca.betas <- lasso.pca$glmnet.fit$beta
lasso.pca.betas <- lasso.pca.betas[,ncol(lasso.pca.betas)] %>% abs
lasso.pca.df <- data.frame(Symbol = names(lasso.pca.betas), Importance = lasso.pca.betas) %>% arrange(desc(Importance))
write.xlsx(lasso.pca.df, "./lasso.pca.xlsx")
saveRDS.gz(lasso.pca, "./save/lasso.pca.rda")

lasso.pco <- cv.glmnet(x = t(exprs(lumi.pco)), y = droplevels(lumi.pco$Status), family = "binomial", type.measure = "auc", nfolds = 3, alpha = 1)
lasso.pco.betas <- lasso.pco$glmnet.fit$beta
lasso.pco.betas <- lasso.pco.betas[,ncol(lasso.pco.betas)] %>% abs
lasso.pco.df <- data.frame(Symbol = names(lasso.pco.betas), Importance = lasso.pco.betas) %>% arrange(desc(Importance))
write.xlsx(lasso.pco.df, "./lasso.pco.xlsx")
saveRDS.gz(lasso.pco, "./save/lasso.pco.rda")

lasso.cc <- cv.glmnet(x = t(exprs(lumi.cc)), y = droplevels(lumi.cc$Status), family = "binomial", type.measure = "auc", nfolds = 3, alpha = 1)
lasso.cc.betas <- lasso.cc$glmnet.fit$beta
lasso.cc.betas <- lasso.cc.betas[,ncol(lasso.cc.betas)] %>% abs
lasso.cc.df <- data.frame(Symbol = names(lasso.cc.betas), Importance = lasso.cc.betas) %>% arrange(desc(Importance))
write.xlsx(lasso.cc.df, "./lasso.cc.xlsx")
saveRDS.gz(lasso.cc, "./save/lasso.cc.rda")

#Elastic Net Regression
elasticnet.pca <- cv.glmnet(x = t(exprs(lumi.pca)), y = droplevels(lumi.pca$Status), family = "binomial", type.measure = "auc", nfolds = 3, alpha = 0.5)
elasticnet.pca.betas <- elasticnet.pca$glmnet.fit$beta
elasticnet.pca.betas <- elasticnet.pca.betas[,ncol(elasticnet.pca.betas)] %>% abs
elasticnet.pca.df <- data.frame(Symbol = names(elasticnet.pca.betas), Importance = elasticnet.pca.betas) %>% arrange(desc(Importance))
write.xlsx(elasticnet.pca.df, "./elasticnet.pca.xlsx")
saveRDS.gz(elasticnet.pca, "./save/elasticnet.pca.rda")

elasticnet.pco <- cv.glmnet(x = t(exprs(lumi.pco)), y = droplevels(lumi.pco$Status), family = "binomial", type.measure = "auc", nfolds = 3, alpha = 0.5)
elasticnet.pco.betas <- elasticnet.pco$glmnet.fit$beta
elasticnet.pco.betas <- elasticnet.pco.betas[,ncol(elasticnet.pco.betas)] %>% abs
elasticnet.pco.df <- data.frame(Symbol = names(elasticnet.pco.betas), Importance = elasticnet.pco.betas) %>% arrange(desc(Importance))
write.xlsx(elasticnet.pco.df, "./elasticnet.pco.xlsx")
saveRDS.gz(elasticnet.pco, "./save/elasticnet.pco.rda")

elasticnet.cc <- cv.glmnet(x = t(exprs(lumi.cc)), y = droplevels(lumi.cc$Status), family = "binomial", type.measure = "auc", nfolds = 3, alpha = 0.5)
elasticnet.cc.betas <- elasticnet.cc$glmnet.fit$beta
elasticnet.cc.betas <- elasticnet.cc.betas[,ncol(elasticnet.cc.betas)] %>% abs
elasticnet.cc.df <- data.frame(Symbol = names(elasticnet.cc.betas), Importance = elasticnet.cc.betas) %>% arrange(desc(Importance))
write.xlsx(elasticnet.cc.df, "./elasticnet.cc.xlsx")
saveRDS.gz(elasticnet.cc, "./save/elasticnet.cc.rda")

#Random KNN
rknn.pca <- rknnBeg(t(exprs(lumi.pca)), droplevels(lumi.pca$Status), pk = 0.8, r = 1000 )
rknn.pca.sig <- rknn.pca$vars[[which.max(rknn.pca$mean_accuracy)]]
rknn.pca.df <- data.frame(Symbol = names(rknn.pca.sig), Importance = rknn.pca.sig) 
write.xlsx(rknn.pca.df, "rknn.pca.df.xlsx")
saveRDS.gz(rknn.pca, "./save/rknn.pca.rda")

rknn.pco <- rknnBeg(t(exprs(lumi.pco)), droplevels(lumi.pco$Status), pk = 0.8, r = 1000 )
rknn.pco.sig <- rknn.pco$vars[[which.max(rknn.pco$mean_accuracy)]]
rknn.pco.df <- data.frame(Symbol = names(rknn.pco.sig), Importance = rknn.pco.sig) 
write.xlsx(rknn.pco.df, "rknn.pco.df.xlsx")
saveRDS.gz(rknn.pco, "./save/rknn.pco.rda")

rknn.cc <- rknnBeg(t(exprs(lumi.cc)), droplevels(lumi.cc$Status), pk = 0.8, r = 1000 )
rknn.cc.sig <- rknn.cc$vars[[which.max(rknn.cc$mean_accuracy)]]
rknn.cc.df <- data.frame(Symbol = names(rknn.cc.sig), Importance = rknn.cc.sig) 
write.xlsx(rknn.cc.df, "rknn.cc.df.xlsx")
saveRDS.gz(rknn.cc, "./save/rknn.cc.rda")

#Pooled accuracy table
rf.accuracy <- list((1 - rfcv.pca$error.cv), (1 - rfcv.pco$error.cv), (1 - rfcv.cc$error.cv)) %>% map(max) %>% unlist #needs SE
rglm.accuracies <- list(rglm.pca.table$Accuracy, rglm.pco.table$Accuracy, rglm.cc.table$Accuracy) %>% map(max) %>% unlist #needs SE
svm.accuracies <- list(svm.pca$tot.accuracy, svm.pco$tot.accuracy, svm.cc$tot.accuracy) %>% map(max) %>% unlist #Calculate SE!
sc.accuracies <- list(num.genes.pca$Accuracy, num.genes.pco$Accuracy, num.genes.cc$Accuracy) %>% map(max) %>% unlist #needs SE
ridge.accuracies <- list(ridge.pca$cvm, ridge.pco$cvm, ridge.cc$cvm) %>% map(max) %>% unlist
lasso.accuracies <- list(lasso.pca$cvm, lasso.pco$cvm, lasso.cc$cvm) %>% map(max) %>% unlist
elasticnet.accuracies <- list(elasticnet.pca$cvm, elasticnet.pco$cvm, elasticnet.cc$cvm) %>% map(max) %>% unlist
rknn.accuracies <- list(rknn.pca$mean_accuracy, rknn.pco$mean_accuracy, rknn.cc$mean_accuracy) %>% map(max) %>% unlist #needs SE
lda.accuracies <- c(lda.pca.cv$stat["acc"], lda.pco.cv$stat["acc"], lda.cc.cv$stat["acc"])
dlda.accuracies <- c(dlda.pca.cv$stat["acc"], dlda.pco.cv$stat["acc"], dlda.cc.cv$stat["acc"])

accuracies.table <- rbind(svm.accuracies / 100, sc.accuracies, ridge.accuracies, lasso.accuracies, elasticnet.accuracies, rknn.accuracies, lda.accuracies, dlda.accuracies) %>% data.frame 
colnames(accuracies.table) <- c("PCA", "PCO", "CC")
accuracies.table$mean <- rowMeans(accuracies.table)
accuracies.table$method <- c("SVM", "SC", "Ridge", "Lasso", "Elastic.Net", "RKNN", "LDA", "DLDA") %>% factor
accuracies.table %<>% arrange(desc(mean))

#Accuracy plots
pca.plot <- ggplot(accuracies.table, aes(x = factor(method), y = PCA)) + geom_bar(stat = "identity") + theme_bw()
pca.plot <- pca.plot + theme(axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pca.plot <- pca.plot + ylab("Mean Accuracy") + ylim(0, 1)
CairoPDF("pca.accuracies.pdf", height = 5, width = 7)
print(pca.plot)
dev.off()

pca.plot <- ggplot(accuracies.table, aes(x = factor(method), y = PCO)) + geom_bar(stat = "identity") + theme_bw()
pca.plot <- pca.plot + theme(axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pca.plot <- pca.plot + ylab("Mean Accuracy") + ylim(0, 1)
CairoPDF("pco.accuracies", height = 6, width = 7)
print(pca.plot)
dev.off()

pca.plot <- ggplot(accuracies.table, aes(x = factor(method), y = CC)) + geom_bar(stat = "identity") + theme_bw()
pca.plot <- pca.plot + theme(axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pca.plot <- pca.plot + ylab("Mean Accuracy") + ylim(0, 1)
CairoPDF("cc.accuracies", height = 5, width = 7)
print(pca.plot)
dev.off()

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

#GO plots
pca.molec <- read.xlsx("./enrichr/pca/GO_Molecular_Function.xlsx") %>% slice(c(1,2,3))
pca.molec$Database <- "GO Molecular Function"
pca.panther <- read.xlsx("./enrichr/pca/Panther.xlsx") %>% slice(1)
pca.panther$Database <- "Panther"
pca.reactome <- read.xlsx("./enrichr/pca/Reactome_2015.xlsx") %>% slice(2)
pca.reactome$Database <- "Reactome"
pca.enrichr <- rbind(pca.molec, pca.panther, pca.reactome)
pca.enrichr$Gene.Count <- map(pca.enrichr$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
pca.enrichr$Log.pvalue <- -(log10(pca.enrichr$P.value))

pca.enrichr$GO.Term %<>% str_replace_all("\\ \\(.*$", "") %>% tolower
pca.enrichr$Format.Name <- paste(pca.enrichr$Database, ": ", pca.enrichr$GO.Term, " (", pca.enrichr$Gene.Count, ")", sep = "")
pca.enrichr.plot <- select(pca.enrichr, Format.Name, Log.pvalue) %>% melt(id.vars = "Format.Name") 

p <- ggplot(pca.enrichr.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = pca.enrichr$Format.Name, hjust = "left", aes(y = 0.1))
p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "FALSE",  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
CairoPDF("pca.enrichr", height = 5, width = 9)
print(p)
dev.off()

pco.biol <- read.xlsx("./enrichr/pco/GO_Biological_Process.xlsx") %>% slice(c(1,2))
pco.biol$Database <- "GO Biological Process"
pco.enrichr <- pco.biol
pco.enrichr$Gene.Count <- map(pco.enrichr$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
pco.enrichr$Log.pvalue <- -(log10(pco.enrichr$P.value))

pco.enrichr$GO.Term %<>% str_replace_all("\\ \\(.*$", "") %>% tolower
pco.enrichr$Format.Name <- paste(pco.enrichr$Database, ": ", pco.enrichr$GO.Term, " (", pco.enrichr$Gene.Count, ")", sep = "")
pco.enrichr.plot <- select(pco.enrichr, Format.Name, Log.pvalue) %>% melt(id.vars = "Format.Name") 

p <- ggplot(pco.enrichr.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = pco.enrichr$Format.Name, hjust = "left", aes(y = 0.1))
p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "FALSE",  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
CairoPDF("pco.enrichr", height = 5, width = 9)
print(p)
dev.off()

cc.biol <- read.xlsx("./enrichr/cc/GO_Biological_Process.xlsx") %>% slice(2)
cc.biol$Database <- "GO Biological Process"
cc.molec <- read.xlsx("./enrichr/cc/GO_Molecular_Function.xlsx") %>% slice(1)
cc.molec$Database <- "GO Molecular Function"
cc.enrichr <- rbind(cc.biol, cc.molec)
cc.enrichr$Gene.Count <- map(cc.enrichr$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
cc.enrichr$Log.pvalue <- -(log10(cc.enrichr$P.value))

cc.enrichr$GO.Term %<>% str_replace_all("\\ \\(.*$", "") %>% tolower
cc.enrichr$Format.Name <- paste(cc.enrichr$Database, ": ", cc.enrichr$GO.Term, " (", cc.enrichr$Gene.Count, ")", sep = "")
cc.enrichr.plot <- select(cc.enrichr, Format.Name, Log.pvalue) %>% melt(id.vars = "Format.Name") 

p <- ggplot(cc.enrichr.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = cc.enrichr$Format.Name, hjust = "left", aes(y = 0.1))
p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "FALSE",  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
CairoPDF("cc.enrichr", height = 5, width = 9)
print(p)
dev.off()
