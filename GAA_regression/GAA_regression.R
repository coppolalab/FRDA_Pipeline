#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)

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

lumi.cleaned <- readRDS.gz("../WGCNA_GAA/save/export.lumi.rda")
expr.gaa <- readRDS.gz("../WGCNA_GAA/save/expr.collapse.rda")

#Random Forest
rf.gaa <- randomForest(t(exprs(lumi.cleaned)), lumi.cleaned$GAA1, ntree = 1000)
saveRDS.gz(rf.gaa, "./save/rf.gaa.rda")
rf.gaa.genes <- data.frame("Symbol" = rownames(rf.gaa$importance), "Importance" = as.vector(rf.gaa$importance)) %>% arrange(desc(Importance))
write.xlsx(rf.gaa.genes, "rf.gaa.genes.xlsx")

rfcv.gaa <- rfcv(t(exprs(lumi.cleaned)), lumi.cleaned$GAA1, step = 0.8, ntree = 1000)
saveRDS.gz(rfcv.gaa, "./save/rfcv.gaa.rda")

#rf.bigmtry.pca <- randomForest(t(exprs(lumi.pca)), droplevels(lumi.pca$Status), ntree = 1000, mtry = ncol(t(exprs(lumi.pca))))
#saveRDS.gz(rf.bigmtry.pca, "./save/rf.bigmtry.pca.rda")
#rf.bigmtry.pca.genes <- data.frame("nuID" = rownames(rf.bigmtry.pca$importance), "Importance" = as.vector(rf.bigmtry.pca$importance)) %>% join(fdata) %>% arrange(desc(Importance))
#write.xlsx(rf.bigmtry.pca.genes, "rf.bigmtry.pca.genes.xlsx")

#Random GLM
get.thinned.accuracy <- function(threshold, rglm.object)
{
    thinned.object <- thinRandomGLM(rglm.object, threshold)
    thinned.mse <- (thinned.object$y.original - thinned.object$predictedOOB)^2 %>% mean
    thinned.features <- which(thinned.object$timesSelectedByForwardRegression > 0) %>% length
    return(c(threshold, thinned.mse, thinned.features))
}

rglm.gaa <- randomGLM(t(exprs(lumi.cleaned)), lumi.cleaned$GAA1,  nThreads = 7)
saveRDS.gz(rglm.gaa, "./save/rglm.gaa.rda")
rglm.thresholds <- seq(1, 10, 1)
rglm.gaa.accuracies <- map(rglm.thresholds, get.thinned.accuracy, rglm.gaa) 
rglm.gaa.table <- reduce(rglm.gaa.accuracies, rbind) %>% data.frame
colnames(rglm.gaa.table) <- c("Index", "MSE", "Num.Genes")

rglm.gaa.table$Num.Genes <- factor(rglm.gaa.table$Num.Genes, sort(rglm.gaa.table$Num.Genes, decreasing = TRUE)) 
gaa.plot <- ggplot(rglm.gaa.table, aes(x = Num.Genes, y = MSE, group = 1)) + geom_line() + geom_point()
gaa.plot <- gaa.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, max(rglm.gaa.table$MSE) + 1000)
CairoPDF("rglm.gaa.accuracies", width = 8, height = 5)
print(gaa.plot)
dev.off()

thinned.gaa <- thinRandomGLM(rglm.gaa, 1)
saveRDS.gz(thinned.gaa, "./save/thinned.gaa.rda")
gaa.rglm <- thinned.gaa$timesSelectedByForwardRegression 
gaa.rglm.df <- data.frame(Symbol = rownames(exprs(lumi.cleaned)), times.selected = as.vector(t(gaa.rglm))) %>% filter(times.selected > 0) %>% arrange(desc(times.selected)) 
write.xlsx(gaa.rglm.df, "gaa.rglm.xlsx")

num.genes <- as.character(rglm.gaa.table$Num.Genes) %>% as.numeric
min.mse.num <- num.genes[which.min(rglm.gaa.table$MSE)]
gaa.glm.submit <- gaa.rglm.df[1:min.mse.num,]
source("../../code/GO/enrichr.R")

get.stringdb(gaa.glm.submit, "rglm.gaa")
enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 

rglm.enrichr <- map(enrichr.terms, get.enrichrdata, gaa.glm.submit, FALSE)
names(rglm.enrichr) <- enrichr.terms
map(names(rglm.enrichr), enrichr.wkbk, rglm.enrichr, "rglm")

enrichr.wkbk <- function(database, full.df, subdir)
{
    dataset <- full.df[[database]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    
    dir.create(file.path("./enrichr", subdir), recursive = TRUE)
    filepath = paste(file.path("./enrichr", subdir, database), ".xlsx", sep = "")
    saveWorkbook(wb, filepath, overwrite = TRUE) 
}

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

list.nums <- recursive.multiply(list.nums.init, nrow(exprs(lumi.cleaned)))
list.nums <- c(nrow(exprs(lumi.cleaned)), list.nums)
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
source("../../code/msvmRFE_reg.R")

svm.pca <- svm(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, type = "eps-regression", kernel = "radial", cross = 3, cost = 1, cachesize = 16000)
saveRDS.gz(svm.gaa, "./save/svm.gaa.rda")
svm.pca <- svmRFE(svm.pca.input, k = 10, halve.above = 100)

nfold = 3
svm.gaa.input <- cbind(lumi.cleaned$GAA1, t(exprs(lumi.cleaned)))
nrows.gaa = nrow(svm.gaa.input)
folds.gaa.index = rep(1:nfold, len=nrows.gaa)[sample(nrows.gaa)]
folds.gaa.rows = lapply(1:nfold, function(x) which(folds.gaa.index == x))
svm.gaa.results = lapply(folds.gaa.rows, svmRFE.wrap, svm.gaa.input, k=3, halve.above=3)
saveRDS.gz(svm.gaa.results, "svm.gaa.results.rda")
svm.gaa.features = WriteFeatures(svm.gaa.results, svm.gaa.input, save=F)
colnames(svm.gaa.features) <- "Symbol"
#svm.gaa.genes <- join(svm.gaa.features, fdata)
write.xlsx(svm.gaa.features, "svm.gaa.features.xlsx")

featsweep.gaa = map(list.nums, FeatSweep.wrap, svm.gaa.results, svm.gaa.input)
svm.gaa.accuracies <- map_dbl(featsweep.gaa, extract.svmlist)
svm.gaa.accuracies.df <- data.frame(Num.Genes = list.nums, Accuracy = svm.gaa.accuracies)
saveRDS.gz(featsweep.gaa, "./save/featsweep.gaa.rda")
saveRDS.gz(svm.gaa.accuracies.df, "./save/svm.gaa.accuracies.df.rda")

svm.gaa.accuracies.df$Num.Genes <- factor(svm.gaa.accuracies.df$Num.Genes, sort(svm.gaa.accuracies.df$Num.Genes, decreasing = TRUE)) 
#levels(svm.gaa.accuracies.df$Num.Genes) <- svm.gaa.accuracies.df$Num.Genes
gaa.plot <- ggplot(svm.gaa.accuracies.df, aes(x = Num.Genes, y = Accuracy, group = 1)) + geom_line() + geom_point()
gaa.plot <- gaa.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0, max(svm.gaa.accuracies) + 1000)
CairoPDF("svm.gaa.accuracies", width = 8, height = 5)
print(gaa.plot)
dev.off()

svm.numgenes <- as.character(svm.gaa.accuracies.df$Num.Genes) %>% as.numeric
svm.minmse <- svm.numgenes[which.min(svm.gaa.accuracies.df$Accuracy)]

svm.features.submit <- svm.gaa.features[1:svm.minmse,]
get.stringdb(svm.features.submit, "svm.gaa")

colnames(svm.features.submit) <- c("Symbol", "X1", "X2")
svm.enrichr <- map(enrichr.terms, get.enrichrdata, svm.features.submit, FALSE)
names(svm.enrichr) <- enrichr.terms
map(names(svm.enrichr), enrichr.wkbk, svm.enrichr, "svm")

#Ridge Regression
ridge.gaa <- cv.glmnet(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, type.measure = "mse", nfolds = 3, family = "gaussian", alpha = 0)
ridge.gaa.betas <- ridge.gaa$glmnet.fit$beta
ridge.gaa.betas <- ridge.gaa.betas[,ncol(ridge.gaa.betas)] %>% abs
ridge.gaa.df <- data.frame(Symbol = names(ridge.gaa.betas), Importance = ridge.gaa.betas) %>% arrange(desc(Importance))
write.xlsx(ridge.gaa.df, "./ridge.gaa.xlsx")
saveRDS.gz(ridge.gaa, "./save/ridge.gaa.rda")

#Lasso Regression
lasso.gaa <- cv.glmnet(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, type.measure = "mse", nfolds = 3, family = "gaussian", alpha = 1)
lasso.gaa.betas <- lasso.gaa$glmnet.fit$beta
lasso.gaa.betas <- lasso.gaa.betas[,ncol(lasso.gaa.betas)] %>% abs
lasso.gaa.df <- data.frame(Symbol = names(lasso.gaa.betas), Importance = lasso.gaa.betas) %>% arrange(desc(Importance))
write.xlsx(lasso.gaa.df, "./lasso.gaa.xlsx")
saveRDS.gz(lasso.gaa, "./save/lasso.gaa.rda")

#Elastic Net Regression
elasticnet.gaa <- cv.glmnet(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, type.measure = "mse", nfolds = 3, family = "gaussian", alpha = 0.5)
elasticnet.gaa.betas <- elasticnet.gaa$glmnet.fit$beta
elasticnet.gaa.betas <- elasticnet.gaa.betas[,ncol(elasticnet.gaa.betas)] %>% abs
elasticnet.gaa.df <- data.frame(Symbol = names(elasticnet.gaa.betas), Importance = elasticnet.gaa.betas) %>% arrange(desc(Importance))
write.xlsx(elasticnet.gaa.df, "./elasticnet.gaa.xlsx")
saveRDS.gz(elasticnet.gaa, "./save/elasticnet.gaa.rda")

#Random KNN
rknn.gaa <- rknnBeg(t(exprs(lumi.cleaned)), lumi.cleaned$GAA1, pk = 0.8, r = 1000)
rknn.gaa.sig <- rknn.gaa$vars[[which.max(rknn.gaa$mean_accuracy)]]
rknn.gaa.df <- data.frame(Symbol = names(rknn.gaa.sig), Importance = rknn.gaa.sig) 
write.xlsx(rknn.gaa.df, "rknn.gaa.df.xlsx")
saveRDS.gz(rknn.gaa, "./save/rknn.gaa.rda")

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

#Pooled accuracy table
rf.accuracy <- min(rfcv.gaa$error.cv) #needs SE
rglm.accuracy <- min(rglm.gaa.table$MSE) #needs SE
svm.accuracy <- min(svm.gaa.accuracies) #Calculate SE!
ridge.accuracy <- min(ridge.gaa$cvm)
lasso.accuracy <- min(lasso.gaa$cvm)
elasticnet.accuracy <- min(elasticnet.gaa$cvm)
#rknn.accuracies <- max()

accuracies.table <- data.frame(MSE = c(rglm.accuracy, svm.accuracy, ridge.accuracy, lasso.accuracy, elasticnet.accuracy), method = c("RGLM", "SVM", "Ridge", "Lasso", "Elastic.Net")) 
accuracies.table %<>% arrange(MSE)

#Accuracy plots
gaa.plot <- ggplot(accuracies.table, aes(x = method, y = MSE)) + geom_bar(stat = "identity") + theme_bw()
gaa.plot <- gaa.plot + theme(axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
gaa.plot <- gaa.plot + ylab("Mean Squared Error") 
CairoPDF("mses", height = 5, width = 7)
print(gaa.plot)
dev.off()

gaa.biol <- read.xlsx("./enrichr/svm/GO_Biological_Process.xlsx") %>% slice(c(1,2))
gaa.biol$Database <- "GO Biological Process"
gaa.molec <- read.xlsx("./enrichr/svm/GO_Molecular_Function.xlsx") %>% slice(1)
gaa.molec$Database <- "GO Molecular Function"
gaa.reactome <- read.xlsx("./enrichr/svm/Reactome_2015.xlsx") %>% slice(45)
gaa.reactome$Database <- "Reactome"
gaa.enrichr <- rbind(gaa.biol, gaa.molec, gaa.reactome)
gaa.enrichr$Gene.Count <- map(gaa.enrichr$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
gaa.enrichr$Log.pvalue <- -(log10(gaa.enrichr$P.value))

gaa.enrichr$GO.Term %<>% str_replace_all("\\ \\(.*$", "") %>% tolower
gaa.enrichr$Format.Name <- paste(gaa.enrichr$Database, ": ", gaa.enrichr$GO.Term, " (", gaa.enrichr$Gene.Count, ")", sep = "")
gaa.enrichr.plot <- select(gaa.enrichr, Format.Name, Log.pvalue) %>% melt(id.vars = "Format.Name") 

p <- ggplot(gaa.enrichr.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = gaa.enrichr$Format.Name, hjust = "left", aes(y = 0.1))
p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "FALSE",  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
CairoPDF("gaa.enrichr", height = 5, width = 9)
print(p)
dev.off()

#SDA
gaa.ranking <- sda.ranking(t(exprs(lumi.cleaned)), factor(lumi.cleaned$Two.Bin), diagonal = FALSE)
gaa.ranking.df <- data.frame(Symbol = names(gaa.ranking[,"score"]), Score = gaa.ranking[,"score"]) %>% arrange(desc(Score))
gaa.ranking.top <- slice(gaa.ranking.df, 1:500)
write.xlsx(gaa.ranking.top, "top.sda.xlsx")

threebin.ranking <- sda.ranking(t(exprs(lumi.cleaned)), factor(lumi.cleaned$Three.Bin), diagonal = FALSE)
threebin.df <- as.matrix(threebin.ranking)
dim(threebin.df) <- dim(threebin.ranking)
threebin.ranking.df <- data.frame(Symbol = names(threebin.ranking[,"score"]), Score = threebin.ranking[,"score"]) %>% arrange(desc(Score))
threebin.ranking.top <- slice(threebin.ranking.df)

sda.gaa.linear <- train(x = t(exprs(lumi.cleaned)), y = factor(lumi.cleaned$Repeat.Size), preProcess = c("center", "scale"), metric = "Kappa", method = "sda", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.gaa.poly, "./save/svm.gaa.poly.rda")

#SC
pamr.gaa.linear <- train(x = t(exprs(lumi.cleaned)), y = factor(lumi.cleaned$Repeat.Size), preProcess = c("center", "scale"), metric = "Kappa", method = "pam", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.gaa.poly, "./save/svm.gaa.poly.rda")

#PenalizedLDA
penalizedlda.gaa.linear <- train(x = t(exprs(lumi.cleaned)), y = factor(lumi.cleaned$Repeat.Size), preProcess = c("center", "scale"), metric = "Kappa", method = "PenalizedLDA", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.gaa.poly, "./save/svm.gaa.poly.rda")

#MARS
earth.fitcontrol <- trainControl(method = "cv", verboseIter = TRUE, number = 3, trim = TRUE, allowParallel = FALSE)
earth.gaa <- train(x = t(exprs(lumi.cleaned)), y = factor(lumi.cleaned$Repeat.Size), method = "bagEarthGCV", Use.beta.cache = TRUE, trControl = earth.fitcontrol) 

#ExtraTree
extraTree.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, allowParallel = FALSE, savePredictions = TRUE)
extraTree.fitcontrol <- trainControl(method = "repeatedcv", verboseIter = TRUE, number = 3, repeats = 10,  allowParallel = FALSE, savePredictions = TRUE)
extraTree.gaa <- train(x = t(exprs(lumi.cleaned)), y = factor(lumi.cleaned$Repeat.Size), method = "extraTrees", ntree = 1000, numThreads = 8, tuneGrid = expand.grid(mtry = floor(sqrt(nrow(lumi.long))), numRandomCuts = 1), trControl = extraTree.fitcontrol)
extraTree.varImp <- varImp(extraTree.gaa)

#Random Fern
rferns.gaa <- train(x = data.frame(t(exprs(lumi.cleaned))), y = factor(lumi.cleaned$Repeat.Size), method = "rFerns", trControl = extraTree.fitcontrol)

#RGLM
rglm.fitcontrol <- trainControl(method = "cv", verboseIter = TRUE, number = 3, allowParallel = FALSE, savePredictions = TRUE)
rglm.gaa <- train(x = t(exprs(lumi.cleaned)), y = factor(lumi.cleaned$Repeat.Size), metric = "Kappa", method = "randomGLM", nThreads = 8, tuneGrid = data.frame(maxInteractionOrder = 1), trControl = rglm.fitcontrol)
saveRDS.gz(rglm.gaa, "./save/rglm.gaa.rda")

#GLMnet
glmnet.gaa <- train(x = t(exprs(lumi.cleaned)), y = factor(lumi.cleaned$Repeat.Size), metric = "Kappa", method = "glmnet", trControl = rglm.fitcontrol)
saveRDS.gz(glmnet.gaa, "./save/glmnet.gaa.rda")

#Darch
darch.gaa <- darch(x = t(exprs(lumi.cleaned)), y = factor(lumi.cleaned$Two.Bin), retainData = T)

#RFE
rfe.size <- c(nrow(lumi.cleaned), recursive.multiply(current.value = nrow(lumi.cleaned)))
svm.gaa.rfe <- map(rfe.size, svm.cutoffs, gaa.ranking.df, lumi.cleaned)
svm.gaa.kappa <- map(svm.gaa.rfe, getElement, "results") %>% map(getElement, "Kappa") %>% map_dbl(max)
svm.gaa.kappasd <- map(svm.gaa.rfe, getElement, "results") %>% map(getElement, "KappaSD") %>% map_dbl(max)
svm.rfe.df <- data.frame(Num.Genes = rfe.size, Kappa = svm.gaa.kappa, Kappa.SD = svm.gaa.kappasd)
svm.rfe.df$Num.Genes %<>% factor(levels = sort(rfe.size, decreasing = TRUE))
saveRDS.gz(svm.gaa.rfe, "./save/svm.gaa.rfe.rda")

extraTree.gaa.rfe <- map(rfe.size, extraTree.cutoffs, care.df, lumi.cleaned)
saveRDS.gz(extraTree.gaa.rfe, "./save/extraTree.gaa.rfe.rda")
extraTree.gaa.rmse <- map(extraTree.gaa.rfe, getElement, "results") %>% map(getElement, "RMSE") %>% map_dbl(max)

sda.gaa.rfe <- map(rfe.size, sda.cutoffs, gaa.ranking.df, lumi.cleaned)
saveRDS.gz(sda.gaa.rfe, "./save/sda.gaa.rfe.rda")
sda.gaa.kappa <- map(sda.gaa.rfe, getElement, "results") %>% map(getElement, "Kappa") %>% map_dbl(max)
sda.gaa.kappasd <- map(sda.gaa.rfe, getElement, "results") %>% map(getElement, "KappaSD") %>% map_dbl(max)
sda.rfe.df <- data.frame(Num.Genes = rfe.size, Kappa = sda.gaa.kappa, Kappa.SD = sda.gaa.kappasd)
sda.rfe.df$Num.Genes %<>% factor(levels = sort(rfe.size, decreasing = TRUE))

kappa.plot <- ggplot(svm.rfe.df, aes(x = Num.Genes, y = Kappa, group = 1)) + ylim(c(0,1))
kappa.plot <- kappa.plot + geom_errorbar(aes(ymax = Kappa + Kappa.SD, ymin = Kappa - Kappa.SD), width = 0.25) + geom_line() + geom_point()
kappa.plot <- kappa.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Number of Genes") + ylab("Kappa") 
CairoPDF("svm.kappa", width = 10, height = 5)
print(kappa.plot)
dev.off()

sda.top.expr <- lumi.cleaned[gaa.ranking.df$Symbol[1:5],]
CairoPDF("sda.top", height = 6, width = 18)
transparentTheme(trans = 0.9)
featurePlot(x = t(exprs(sda.top.expr)), y = sda.top.expr$Two.Bin, plot = "density", scales = list(x = list(relation = "free"), y = list(relation = "free")), adjust = 1.5, pch = "|", layout = c(5,1), auto.key = list(columns = 2), lwd = 2)
dev.off()

svm.threebin.rfe <- map(rfe.size, svm.threebin, threebin.ranking.df, lumi.cleaned)
svm.threebin.kappa <- map(svm.threebin.rfe, getElement, "results") %>% map(getElement, "Kappa") %>% map_dbl(max)
svm.threebin.kappasd <- map(svm.threebin.rfe, getElement, "results") %>% map(getElement, "KappaSD") %>% map_dbl(max)
svm.threebin.df <- data.frame(Num.Genes = rfe.size, Kappa = svm.threebin.kappa, Kappa.SD = svm.threebin.kappasd)
svm.threebin.df$Num.Genes %<>% factor(levels = sort(rfe.size, decreasing = TRUE))
saveRDS.gz(svm.threebin.rfe, "./save/svm.threebin.rfe.rda")

#Needs some work to figure this out
sda.preds <- sda.gaa.rfe[[14]]$pred
p <- ggplot(sda.preds, aes(x = obs, y = pred)) + geom_point() + facet_wrap( ~ Resample)
#p <- p + ylim(c(0, 2000)) + stat_smooth()
CairoPDF("sda.linear.preds", height = 30, width = 40)
print(p)
dev.off()
