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
options(java.parameters = "-Xmx16g")
library(caret)
library(MASS) #<-- this fool keeps overriding dplyr
library(AppliedPredictiveModeling)

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
library(broom)

source("../common_functions.R")
lumi.cleaned <- readRDS.gz("../WGCNA_GAA/save/rmcov.baseline.rda")
full.cor.df <- read.xlsx("../WGCNA_GAA/gaa.cor.xlsx")
lumi.top <- lumi.cleaned[full.cor.df$Symbol[1:500],]

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

full.cor.df <- gaa.cor(lumi.cleaned, "gaa.cor.xlsx")

svm.gaa.linear <- train(x = t(exprs(lumi.cleaned)), y = factor(lumi.cleaned$Repeat.Size), preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.gaa.poly, "./save/svm.gaa.poly.rda")

lasso.gaa <- train(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, preProcess = c("center", "scale"), metric = "RMSE", method = "lasso", trControl = svm.fitcontrol) #slow
blasso.gaa <- train(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, preProcess = c("center", "scale"), metric = "RMSE", method = "blasso", trControl = svm.fitcontrol) #slow
rqlasso.gaa <- train(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, preProcess = c("center", "scale"), metric = "RMSE", method = "rqlasso", trControl = svm.fitcontrol) #slow
relaxo.gaa <- train(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, preProcess = c("center", "scale"), metric = "RMSE", method = "relaxo", trControl = svm.fitcontrol) #slow

svm.gaa.linear.fail <- train(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, preProcess = c("center", "scale"), metric = "RMSE", method = "svmLinear", trControl = svm.fitcontrol.boot)

model.design <- model.matrix(~ Repeat.Size, pData(lumi.cleaned)) #Make covariate matrix for limma
fit <- lmFit(lumi.cleaned, model.design)
fit$df.residual <- (fit$df.residual - 3)
fitb <- eBayes(fit)

top.object <- topTable(fitb, coef = 2, n = Inf)
top.object$Symbol <- rownames(top.object)

svm.continuous.rfe.radial <- map(rfe.size, svm.continuous.radial, gaa.ranking.df, lumi.cleaned)
saveRDS.gz(svm.continuous.rfe.radial, "./save/svm.continuous.rfe.radial.rda")

#Radial
svm.gaa.radial <- train(x = t(exprs(lumi.cleaned)), y = lumi.top$GAA1, metric = "RMSE", method = "svmRadial", trControl = svm.fitcontrol)
saveRDS.gz(svm.gaa.radial, "./save/svm.gaa.radial.rda")

lm.gaa <- train(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, metric = "RMSE", method = "lm", trControl = svm.fitcontrol)

pcr.gaa <- train(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, metric = "RMSE", method = "pcr", trControl = svm.fitcontrol)

dnn.gaa <- train(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, metric = "RMSE", method = "dnn", trControl = svm.fitcontrol) #terrible

#Polynomial
svm.gaa.poly.long <- train(x = t(exprs(lumi.top.long)), y = lumi.top.long$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmPoly", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.gaa.poly.long, "./save/svm.gaa.poly.long.rda")

svm.gaa.poly <- train(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmPoly", trControl = svm.fitcontrol.boot.par)
saveRDS.gz(svm.gaa.poly, "./save/svm.gaa.poly.rda")

#Random Forest Methods
rf.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, allowParallel = FALSE, savePredictions = TRUE)

#Ranger
rf.gaa <- train(x = t(exprs(lumi.top.long)), y = lumi.top.long$GAA1, metric = "Rsquared", method = "ranger", num.trees = 1000, tuneGrid = data.frame(mtry = floor(sqrt(nrow(lumi.top.long)))), trControl = rf.fitcontrol)
saveRDS.gz(rf.gaa, "./save/rf.gaa.rda")

#ExtraTrees
extraTree.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, allowParallel = FALSE, savePredictions = TRUE)
extraTree.gaa <- train(x = t(exprs(lumi.top)), y = lumi.top$GAA1, method = "extraTrees", ntree = 1000, numThreads = 8, tuneGrid = expand.grid(mtry = nrow(lumi.top), numRandomCuts = 1), trControl = extraTree.fitcontrol)
saveRDS.gz(extraTree.gaa, "./save/extraTree.gaa.rda")

extraTree.gaa.varImp <- varImp(extraTree.gaa)
extraTree.gaa.df <- extraTree.gaa.varImp$importance
extraTree.gaa.df$Symbol <- rownames(extraTree.gaa.df)
extraTree.gaa.df %<>% arrange(desc(Overall))
write.xlsx(extraTree.gaa.df,  "extratrees.gaa.xlsx")

extraTree.gaa.submit <- extraTree.gaa.df[1:500,]
enrichr.terms <- list("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "Panther_2016") 
extraTree.gaa.enrichr <- map(enrichr.terms, get.enrichrdata, extraTree.gaa.submit, FALSE)
names(extraTree.gaa.enrichr) <- enrichr.terms
map(names(extraTree.gaa.enrichr), enrichr.wkbk, extraTree.gaa.enrichr, "./enrichr/gaa")

extraTree.preds <- extraTree.gaa$pred

p <- ggplot(extraTree.preds, aes(x = obs, y = pred)) + geom_point() + facet_wrap( ~ Resample) + stat_smooth()
#p <- p + ylim(c(0, 2000)) + stat_smooth()
CairoPDF("extraTree.preds", height = 30, width = 40)
print(p)
dev.off()

gaa.gobiol.file <- "./enrichr/gaa/GO_Biological_Process_2015.xlsx"
gaa.gobiol <- read.xlsx(gaa.gobiol.file)
gaa.gobiol$Database <- "GO Biological Process"
gaa.gobiol$Num.Genes <- map(gaa.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
gaa.gobiol.filter <- filter(gaa.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(gaa.gobiol.filter, extraTree.gaa.submit$Symbol, file_path_sans_ext(gaa.gobiol.file))

gaa.gomolec.file <- "./enrichr/gaa/GO_Molecular_Function_2015.xlsx"
gaa.gomolec <- read.xlsx(gaa.gomolec.file)
gaa.gomolec$Database <- "GO Molecular Function"
gaa.gomolec$Num.Genes <- map(gaa.gomolec$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
gaa.gomolec.filter <- filter(gaa.gomolec, Num.Genes > 4) %>% filter(P.value < 0.05)

#Bagged MARS
earth.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, allowParallel = FALSE)
earth.gaa <- train(x = t(exprs(lumi.top.long)), y = lumi.top.long$GAA1, method = "bagEarthGCV", Use.beta.cache = TRUE, trControl = earth.fitcontrol)
earth.gaa.resample <- earth.gaa$resample
saveRDS.gz(earth.gaa, "./save/earth.gaa.rda")

#RGLM
rglm.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, allowParallel = FALSE)
rglm.gaa <- train(x = t(exprs(lumi.top)), y = lumi.top$GAA1, method = "randomGLM", nThreads = 8, tuneGrid = data.frame(maxInteractionOrder = 1), trControl = rglm.fitcontrol)
saveRDS.gz(rglm.gaa, "./save/rglm.gaa.rda")

#GLMNET - meh
glmnet.gaa <- train(x = t(exprs(lumi.top.long)), y = lumi.top.long$GAA1, metric = "Rsquared", method = "glmnet", trControl = rf.fitcontrol)
saveRDS.gz(glmnet.gaa, "./save/glmnet.gaa.rda")

#QRF
qrf.gaa <- train(x = t(exprs(lumi.top.long)), y = lumi.top.long$GAA1, method = "qrf", nthreads = 8, ntree = 1000, tuneGrid = data.frame(mtry = floor(sqrt(nrow(lumi.top.long)))), trControl = rf.fitcontrol) 
saveRDS.gz(qrf.gaa, "./save/qrf.gaa.rda")

#RRF - maybe...
rrf.gaa <- train(x = t(exprs(lumi.top.long)), y = lumi.top.long$GAA1, method = "RRFglobal", ntree = 1000, tuneGrid = expand.grid(mtry = floor(sqrt(nrow(lumi.top.long))), coefReg = c(0.01, 0.5, 1)), trControl = rf.fitcontrol)
saveRDS.gz(rrf.gaa, "rrf.gaa.rda")

#LARS - meh
lars.gaa <- train(x = t(exprs(lumi.top.long)), y = lumi.top.long$GAA1, method = "lars", use.Gram = FALSE, trace = TRUE, trControl = rf.fitcontrol)
saveRDS.gz(lars.gaa, "lars.gaa.rda")

#FRBS - forget it for now
gfs.thrift.gaa <- train(x = t(exprs(lumi.top)), y = lumi.top$GAA1, method = "GFS.THRIFT", trControl = rglm.fitcontrol)

svm.featureplot <- function(svm.object, lumi.object, filename)
{
    svm.gaa.model <- svm.object$finalModel
    svm.gaa.imp <- t(svm.gaa.model@coef) %*% svm.gaa.model@xmatrix
    svm.gaa.weight <- svm.gaa.imp[1,] ^ 2
    svm.gaa.df <- data.frame(Symbol = names(svm.gaa.weight), Importance = svm.gaa.weight) %>% arrange(desc(Importance))

    svm.gaa.top5 <- svm.gaa.df$Symbol[1:5]
    svm.gaa.top5.index <- match(svm.gaa.top5, featureNames(lumi.object))
    svm.gaa.top5.expr <- t(exprs(lumi.object[svm.gaa.top5.index,]))

    CairoPDF(filename, height = 6, width = 18)
        theme1 <- trellis.par.get()
        theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
        theme1$plot.symbol$pch = 16
        theme1$plot.line$col = rgb(1, 0, 0, .7)
        theme1$plot.line$lwd <- 2
        trellis.par.set(theme1)
        fp <- featurePlot(x = svm.gaa.top5.expr, y = lumi.object$GAA1, plot = "scatter", type = c("p", "smooth"), span = 0.5, layout = c(5,1))
        print(fp)
    dev.off()
    return(svm.gaa.df)
}

svm.cutoffs <- function(cutoff, top.object, lumi.object)
{
    top.genes <- top.object$Symbol[1:cutoff]
    lumi.top <- lumi.object[top.genes,]
    print(dim(lumi.top))
    svm.linear <- train(x = t(exprs(lumi.top)), y = lumi.top$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmLinear", trControl = svm.fitcontrol.boot.par) 
    return(svm.linear)
}

lumi.short <- lumi.cleaned[,lumi.cleaned$GAA1 <= 466]

gaa.densityplot(lumi.long$GAA1, "gaa_density_long")

short.cor.df <- gaa.cor(lumi.short, "gaa.cor.short.xlsx") 

#Linear
cutoffs <- c(seq(100, 500, 100), seq(1000, 3000, 500))
svm.gaa.cutoffs.long <- map(cutoffs, svm.cutoffs, long.cor.df, lumi.long)
svm.gaa.cutoffs.short <- map(cutoffs, svm.cutoffs, short.cor.df, lumi.short)
svm.gaa.cutoffs.intermediate <- map(cutoffs, svm.cutoffs, intermediate.cor.df, lumi.intermediate)
#svm.gaa.cutoffs.intermediate2 <- map(cutoffs, svm.cutoffs, intermediate2.cor.df, lumi.intermediate2)
svm.gaa.cutoffs.all <- map(cutoffs, svm.cutoffs, full.cor.df, lumi.cleaned)

svm.gaa.rsquared.long <- map(svm.gaa.cutoffs.long, getElement, "results") %>% map(getElement, "Rsquared") %>% map_dbl(max)
svm.gaa.rsquared.short <- map(svm.gaa.cutoffs.short, getElement, "results") %>% map(getElement, "Rsquared") %>% map_dbl(max)
svm.gaa.rsquared.intermediate <- map(svm.gaa.cutoffs.intermediate, getElement, "results") %>% map(getElement, "Rsquared") %>% map_dbl(max)
#svm.gaa.rsquared.intermediate2 <- map(svm.gaa.cutoffs.intermediate2, getElement, "results") %>% map(getElement, "Rsquared") %>% map_dbl(max)
svm.gaa.rsquared.all <- map(svm.gaa.cutoffs.all, getElement, "results") %>% map(getElement, "Rsquared") %>% map_dbl(max)

svm.gaa.rsquared.df <- cbind(svm.gaa.rsquared.long, svm.gaa.rsquared.short, svm.gaa.rsquared.intermediate, svm.gaa.rsquared.all) %>% data.frame
colnames(svm.gaa.rsquared.df) <- c("Long", "Short", "Intermediate", "All")
svm.gaa.rsquared.df$Num.Genes <- cutoffs
svm.rsquared.gather <- gather(svm.gaa.rsquared.df, Repeat.Length, Rsquared, -Num.Genes)

svm.gaa.rmse.long <- map(svm.gaa.cutoffs.long, getElement, "results") %>% map(getElement, "RMSE") %>% map_dbl(max)
svm.gaa.rmse.short <- map(svm.gaa.cutoffs.short, getElement, "results") %>% map(getElement, "RMSE") %>% map_dbl(max)
svm.gaa.rmse.intermediate <- map(svm.gaa.cutoffs.intermediate, getElement, "results") %>% map(getElement, "RMSE") %>% map_dbl(max)
#svm.gaa.rmse.intermediate2 <- map(svm.gaa.cutoffs.intermediate2, getElement, "results") %>% map(getElement, "RMSE") %>% map_dbl(max)
svm.gaa.rmse.all <- map(svm.gaa.cutoffs.all, getElement, "results") %>% map(getElement, "RMSE") %>% map_dbl(max)

svm.gaa.rmse.df <- cbind(svm.gaa.rmse.long, svm.gaa.rmse.short, svm.gaa.rmse.intermediate, svm.gaa.rmse.all) %>% data.frame
colnames(svm.gaa.rmse.df) <- c("Long", "Short", "Intermediate", "All")
svm.gaa.rmse.df$Num.Genes <- cutoffs
svm.rmse.gather <- gather(svm.gaa.rmse.df, Repeat.Length, RMSE, -Num.Genes)

svm.gaa.rmse.long.sd <- map(svm.gaa.cutoffs.long, getElement, "results") %>% map(getElement, "RMSESD") %>% map_dbl(max)
svm.gaa.rmse.short.sd <- map(svm.gaa.cutoffs.short, getElement, "results") %>% map(getElement, "RMSESD") %>% map_dbl(max)
svm.gaa.rmse.intermediate.sd <- map(svm.gaa.cutoffs.intermediate, getElement, "results") %>% map(getElement, "RMSESD") %>% map_dbl(max)
#svm.gaa.rmse.intermediate2.sd <- map(svm.gaa.cutoffs.intermediate2, getElement, "results") %>% map(getElement, "RMSESD") %>% map_dbl(max)
svm.gaa.rmse.all.sd <- map(svm.gaa.cutoffs.all, getElement, "results") %>% map(getElement, "RMSESD") %>% map_dbl(max)

svm.gaa.rmse.sd.df <- cbind(svm.gaa.rmse.long.sd, svm.gaa.rmse.short.sd, svm.gaa.rmse.intermediate.sd, svm.gaa.rmse.all.sd) %>% data.frame
colnames(svm.gaa.rmse.sd.df) <- c("Long", "Short", "Intermediate", "All")
svm.gaa.rmse.sd.df$Num.Genes <- cutoffs
svm.rmse.sd.gather <- gather(svm.gaa.rmse.sd.df, Repeat.Length, RMSE.SD, -Num.Genes)

svm.gaa.rsquared.long.sd <- map(svm.gaa.cutoffs.long, getElement, "results") %>% map(getElement, "RsquaredSD") %>% map_dbl(max)
svm.gaa.rsquared.short.sd <- map(svm.gaa.cutoffs.short, getElement, "results") %>% map(getElement, "RsquaredSD") %>% map_dbl(max)
svm.gaa.rsquared.intermediate.sd <- map(svm.gaa.cutoffs.intermediate, getElement, "results") %>% map(getElement, "RsquaredSD") %>% map_dbl(max)
#svm.gaa.rsquared.intermediate2.sd <- map(svm.gaa.cutoffs.intermediate2, getElement, "results") %>% map(getElement, "RsquaredSD") %>% map_dbl(max)
svm.gaa.rsquared.all.sd <- map(svm.gaa.cutoffs.all, getElement, "results") %>% map(getElement, "RsquaredSD") %>% map_dbl(max)

svm.gaa.rsquared.sd.df <- cbind(svm.gaa.rsquared.long.sd, svm.gaa.rsquared.short.sd, svm.gaa.rsquared.intermediate.sd, svm.gaa.rsquared.all.sd) %>% data.frame
colnames(svm.gaa.rsquared.sd.df) <- c("Long", "Short", "Intermediate", "All")
svm.gaa.rsquared.sd.df$Num.Genes <- cutoffs
svm.rsquared.sd.gather <- gather(svm.gaa.rsquared.sd.df, Repeat.Length, Rsquared.SD, -Num.Genes)

svm.cutoffs.df <- join(svm.rsquared.gather, svm.rsquared.sd.gather) %>% join(svm.rmse.gather) %>% join(svm.rmse.sd.gather)
svm.cutoffs.df$Num.Genes %<>% factor(levels = sort(cutoffs, decreasing = TRUE))

rsquared.plot <- ggplot(svm.cutoffs.df, aes(x = Num.Genes, y = Rsquared, group = Repeat.Length, color = Repeat.Length))  
rsquared.plot <- rsquared.plot + geom_errorbar(aes(ymax = Rsquared + Rsquared.SD, ymin = Rsquared - Rsquared.SD), width = 0.25) + geom_line() + geom_point()
rsquared.plot <- rsquared.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Number of Genes") + ylab("Rsquared")
CairoPDF("svm.rsquared", width = 10, height = 5)
print(rsquared.plot)
dev.off()

rmse.plot <- ggplot(svm.cutoffs.df, aes(x = Num.Genes, y = RMSE, group = Repeat.Length, color = Repeat.Length))  
rmse.plot <- rmse.plot + geom_errorbar(aes(ymax = RMSE + RMSE.SD, ymin = RMSE - RMSE.SD), width = 0.25) + geom_line() + geom_point()
rmse.plot <- rmse.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Number of Genes") + ylab("rmse") 
#rmse.plot <- rmse.plot + ylim(0, 500) + scale_y_log10(breaks = c(seq(10, 100, 10), seq(150, 500, 50)))
CairoPDF("svm.rmse", width = 10, height = 5)
print(rmse.plot)
dev.off()

lumi.top <- lumi.cleaned[full.cor.df$Symbol[1:1000],]
lumi.top.long <- lumi.long[long.cor.df$Symbol[1:500],]
lumi.top.short <- lumi.short[short.cor.df$Symbol[1:500],]
lumi.top.intermediate <- lumi.intermediate[intermediate.cor.df$Symbol[1:500],]

lumi.mads <- apply(exprs(lumi.cleaned), 1, mad) %>% sort(decreasing = TRUE)
lumi.top <- lumi.cleaned[names(lumi.mads[1:2500]),]

saveRDS.gz(lumi.top.long, "./save/lumi.top.long.rda")
saveRDS.gz(lumi.top.intermediate, "./save/lumi.top.intermediate.rda")
saveRDS.gz(lumi.top.short, "./save/lumi.top.short.rda")
long.cor.df <- gaa.cor(lumi.long, "gaa.cor.long.xlsx")
intermediate.cor.df <- gaa.cor(lumi.intermediate, "gaa.cor.intermediate.xlsx")
#intermediate2.cor.df <- gaa.cor(lumi.intermediate2, "gaa.cor.intermediate2.xlsx")
gaa.densityplot(lumi.short$GAA1, "gaa_density_short")
gaa.densityplot(lumi.intermediate$GAA1, "gaa_density_intermediate")
#gaa.densityplot(lumi.intermediate2$GAA1, "gaa_density_intermediate2")
lumi.intermediate <- lumi.cleaned[,lumi.cleaned$GAA1 > 466 & lumi.cleaned$GAA1 <= 800]
#lumi.intermediate2 <- lumi.cleaned[,lumi.cleaned$GAA1 > 650 & lumi.cleaned$GAA1 <= 800]
lumi.long <- lumi.cleaned[,lumi.cleaned$GAA1 > 800]
#gaa.genes <- read.xlsx("../WGCNA_GAA/all.gaa.xlsx")

svm.gaa.linear.long <- train(x = t(exprs(lumi.top.long)), y = lumi.top.long$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmLinear", trControl = svm.fitcontrol.boot)

saveRDS.gz(svm.gaa.linear.long, "./save/svm.gaa.linear.long.rda")

svm.gaa.linear.short <- train(x = t(exprs(lumi.top.short)), y = lumi.top.short$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmLinear", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.gaa.linear.short, "./save/svm.gaa.linear.short.rda")

svm.gaa.linear.intermediate <- train(x = t(exprs(lumi.top.intermediate)), y = lumi.top.intermediate$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmLinear", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.gaa.linear.intermediate, "./save/svm.gaa.linear.intermediate.rda")

#RFE - is this necessary?
#extraTree.rfeFitControl <- trainControl(method = "none", verboseIter = TRUE, allowParallel = FALSE)
#extraTree.rfeControl <- rfeControl(functions = list(summary = caretFuncs$summary, fit = caretFuncs$fit, rank = caretFuncs$rank, pred = caretFuncs$pred, selectSize = pickSizeTolerance, selectVar = caretFuncs$selectVar), rerank = FALSE, method = "repeatedcv", number = 3, repeats = 10, verbose = TRUE, allowParallel = FALSE)

#rfe.size <- c(nrow(lumi.cleaned), recursive.multiply(current.value = nrow(lumi.cleaned)))

#extraTree.rfeRun.gaa <- rfe(x = t(exprs(lumi.cleaned)), y = lumi.cleaned$GAA1, sizes = rfe.size, rfeControl = extraTree.rfeControl, method = "extraTrees", ntree = 1000, numThreads = 8, tuneGrid = expand.grid(mtry = floor(sqrt(nrow(lumi.cleaned))), numRandomCuts = 1), trControl = extraTree.rfeFitControl)
#saveRDS.gz(extraTree.rfeRun.gaa, "./save/extraTree.gaa.rfe.rda")
#extraTree.rfe.gaa.top5 <- extraTree.rfeRun.gaa$optVariables[1:5] 
#extraTree.rfe.gaa.top5.index <- match(extraTree.rfe.gaa.top5, featureNames(lumi.cleaned))
#extraTree.rfe.gaa.top5.expr <- t(exprs(lumi.cleaned[extraTree.rfe.gaa.top5.index,]))

#CairoPDF("extratrees.gaa.top", height = 6, width = 18)
#theme1 <- trellis.par.get()
#theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
#theme1$plot.symbol$pch = 16
#theme1$plot.line$col = rgb(1, 0, 0, .7)
#theme1$plot.line$lwd <- 2
#trellis.par.set(theme1)
#fp <- featurePlot(x = extraTree.rfe.gaa.top5.expr, y = lumi.cleaned$GAA1, plot = "scatter", type = c("p", "smooth"), span = 0.5, layout = c(5,1))
#print(fp)
#dev.off()

#extraTree.rfe.accuracies <- extraTree.rfeRun.gaa$results
#extraTree.rfe.accuracies$Variables <- factor(extraTree.rfe.accuracies$Variables, sort(unique(extraTree.rfe.accuracies$Variables), decreasing = TRUE)) 

#rfe.plot <- ggplot(extraTree.rfe.accuracies, aes(x = Variables, y = RMSE, group = 1))  
#rfe.plot <- rfe.plot + geom_errorbar(aes(ymax = RMSE + RMSESD, ymin = RMSE - RMSESD), width = 0.25, color = "black") + geom_line() + geom_point()
#rfe.plot <- rfe.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(100,260) + xlab("Number of Genes") + ylab("RMSE (# GAA repeats)")
#CairoPDF("extraTree.rfe.accuracies", width = 10, height = 5)
#print(rfe.plot)
#dev.off()

svm.long.df <- svm.featureplot(svm.gaa.linear.long, lumi.top.long, "svm.gaa.long.top5")
svm.intermediate.df <- svm.featureplot(svm.gaa.linear.intermediate, lumi.top.intermediate, "svm.gaa.intermediate.top5")
svm.short.df <- svm.featureplot(svm.gaa.linear.short, lumi.top.short, "svm.gaa.short.top5")

lumi.second.intersect <- intersect(featureNames(lumi.cleaned), featureNames(lumi.second)) #Get intersect of two time points
lumi.baseline.intersect <- lumi.cleaned[lumi.second.intersect,] #get baseline that contains intersect of two timepoints
lumi.two.intersect <- lumi.second[lumi.second.intersect,] #get second that contains intersect of two timepoints

lumi.bst <- lumi.baseline.intersect[,lumi.cleaned$GAA1 <= 466] #Long repeats from baseline intersect
lumi.bit <- lumi.baseline.intersect[,lumi.cleaned$GAA1 > 466 & lumi.cleaned$GAA1 <= 800] #Intermediate repeats from baseline intersect
lumi.blt <- lumi.baseline.intersect[,lumi.cleaned$GAA1 > 800] #Short repeats from baseline intersect

lumi.second.short <- lumi.two.intersect[,lumi.second$GAA1 <= 475] #Short repeats from second intersect
lumi.second.intermediate <- lumi.two.intersect[,lumi.second$GAA1 > 475 & lumi.second$GAA1 <= 800] #Intermediate repeats from second intersect
lumi.second.long <- lumi.two.intersect[,lumi.second$GAA1 > 800] #Long repeats from second intersect

short.cor.second <- gaa.cor(lumi.bst, "second.cor.short.xlsx") #Short correlation from baseline intersect
long.cor.second <- gaa.cor(lumi.blt, "second.cor.long.xlsx") #Long correlation from baseline intersect
intermediate.cor.second <- gaa.cor(lumi.bit, "second.cor.intermediate.xlsx") #Intermediate correlation from baseline intersect

lumi.top.sl <- lumi.blt[long.cor.second$Symbol[1:500],] #Subset long baseline intersect for top 300
lumi.top.si <- lumi.bit[intermediate.cor.second$Symbol[1:500],] #Subset intermediate baseline intersect for top 300
lumi.top.ss <- lumi.bst[short.cor.second$Symbol[1:500],] #Subset short baseline intersect for top 300

lumi.test.sl <- lumi.second.long[long.cor.second$Symbol[1:500],] #Subset long baseline intersect for test 300
lumi.test.si <- lumi.second.intermediate[intermediate.cor.second$Symbol[1:500],] #Subset intermediate baseline intersect for test 300
lumi.test.ss <- lumi.second.short[short.cor.second$Symbol[1:500],] #Subset short baseline intersect for test 300

svm.second.linear.long <- train(x = t(exprs(lumi.blt)), y = lumi.blt$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmLinear", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.second.linear.long, "./save/svm.second.linear.long.rda")

svm.second.linear.intermediate <- train(x = t(exprs(lumi.bit)), y = lumi.bit$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmLinear", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.second.linear.intermediate, "./save/svm.second.linear.intermediate.rda")

svm.second.linear.short <- train(x = t(exprs(lumi.bst)), y = lumi.bst$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmLinear", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.second.linear.short, "./save/svm.second.linear.short.rda")

lumi.baseline.longer <- lumi.baseline.intersect[,lumi.baseline.intersect$GAA1 > 466]
svm.second.linear <- train(x = t(exprs(lumi.baseline.intersect)), y = lumi.baseline.intersect$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmLinear", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.second.linear, "./save/svm.second.linear.rda")

svm.second.longer <- train(x = t(exprs(lumi.baseline.longer)), y = lumi.baseline.longer$GAA1, preProcess = c("center", "scale"), metric = "Rsquared", method = "svmLinear", trControl = svm.fitcontrol.boot)
saveRDS.gz(svm.second.longer, "./save/svm.second.longer.rda")

lumi.second.longer <- lumi.two.intersect[,lumi.two.intersect$GAA1 > 466]
svm.second.long.predict <- predict(svm.second.linear.long$finalModel, t(exprs(lumi.second.long)))
svm.second.intermediate.predict <- predict(svm.second.linear.intermediate$finalModel, t(exprs(lumi.second.intermediate)))
svm.second.short.predict <- predict(svm.second.linear.short$finalModel, t(exprs(lumi.second.short)))

svm.second.longer.predict <- predict(svm.second.longer$finalModel, t(exprs(lumi.second.short)))
svm.second.predict <- predict(svm.second.linear$finalModel, t(exprs(lumi.second.longer)))
