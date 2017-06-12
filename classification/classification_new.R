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

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', 
                  values = as.character(status.lasso.df$Symbol), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")
lasso.annot <- left_join(status.lasso.df, bm.table) %>% 
    dplyr::select(Symbol, Definition, dplyr::contains("Coefficient")) 
ClassifierWorkbook(lasso.annot, "lasso.table.xlsx")
STOP

#ranger.final <- train(x = data.frame(t(exprs(lumi.import.top))), y = droplevels(factor(lumi.import.top$Status)), 
                      #preProcess = c("center", "scale"), metric = "Kappa", method = "ranger", trControl = ranger.fitcontrol, 
                      #num.trees = 1000, tuneGrid = data.frame(mtry = nrow(lumi.import.top)))
#ranger.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, allowParallel = FALSE, savePredictions = TRUE)
#lumi.pco.top5 <- data.frame(Status = pco.collapse$Status, t(exprs(pco.collapse[pco.lasso.df$Symbol[1:5],]))) %>% gather(Gene, Expression, -Status)
#lumi.pco.top5$Gene %<>% factor(levels = pco.lasso.df$Symbol[1:5])
#top5.plot <- ggplot(lumi.pco.top5, aes(x = Expression, color = Status)) + geom_density()
#top5.plot <- top5.plot + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#top5.plot <- top5.plot + facet_wrap( ~ Gene, ncol = 5, scales = "free") + scale_color_manual(values = c("red", "blue"))
#CairoPDF("svm.pco.top5", width = 20, height = 5)
#print(top5.plot)
#dev.off()

#pca.lasso <- glmnet(t(exprs(pca.collapse)), pca.collapse$Status, "binomial", nlambda = 1000, alpha = 1)
#pca.lasso.coef <- coef(pca.lasso)[-1,1000]
#pca.lasso.df <- data.frame(Symbol = names(pca.lasso.coef), Coefficient.pca = pca.lasso.coef) %>% arrange(desc(abs(Coefficient.pca)))

#pca.collapse.top <- pca.collapse[pca.lasso.coef != 0,]
#svm.pca.final <- train(x = data.frame(t(exprs(pca.collapse.top))), y = droplevels(factor(pca.collapse.top$Status)), preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol)

#pca.collapse.top5 <- data.frame(Status = pca.collapse$Status, t(exprs(pca.collapse[pca.lasso.df$Symbol[1:5],]))) %>% gather(Gene, Expression, -Status)
#pca.collapse.top5$Gene %<>% factor(levels = str_replace(pca.lasso.df$Symbol[1:5], "-", "."))
#top5.plot <- ggplot(pca.collapse.top5, aes(x = Expression, color = Status)) + geom_density()
#top5.plot <- top5.plot + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#top5.plot <- top5.plot + facet_wrap( ~ Gene, ncol = 5, scales = "free") + scale_color_manual(values = c("darkgreen", "blue"))
#CairoPDF("svm.pca.top5", width = 20, height = 5)
#print(top5.plot)
#dev.off()

#cc.collapse$Status %<>% factor(levels = c("Control", "Carrier"))
#cc.lasso <- glmnet(t(exprs(cc.collapse)), cc.collapse$Status, "binomial", nlambda = 1000, alpha = 1)
#cc.lasso.coef <- coef(cc.lasso)[-1,1000] 
#cc.lasso.df <- data.frame(Symbol = names(cc.lasso.coef), Coefficient.cc = cc.lasso.coef) %>% arrange(desc(abs(Coefficient.cc)))

#cc.collapse.top <- cc.collapse[cc.lasso.coef != 0,]
#svm.cc.final <- train(x = t(exprs(cc.collapse.top)), y = droplevels(factor(cc.collapse.top$Status)), preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol)

#cc.collapse.top5 <- data.frame(Status = cc.collapse$Status, t(exprs(cc.collapse[cc.lasso.df$Symbol[1:5],]))) %>% gather(Gene, Expression, -Status)
#cc.collapse.top5$Gene %<>% factor(levels = str_replace(cc.lasso.df$Symbol[1:5], "-", "."))
#top5.plot <- ggplot(cc.collapse.top5, aes(x = Expression, color = Status)) + geom_density()
#top5.plot <- top5.plot + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#top5.plot <- top5.plot + facet_wrap( ~ Gene, ncol = 5, scales = "free") + scale_color_manual(values = c("darkgreen", "red"))
#CairoPDF("svm.cc.top5", width = 20, height = 5)
#print(top5.plot)
#dev.off()


objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

#orig.pco <- isoMDS(dist(t(exprs(pco.collapse)), method = "manhattan"))
#orig.pco.points <- data.frame(orig.pco$points)
#colnames(orig.pco.points) <- c("PC1", "PC2")
#orig.pco.df <- mutate(orig.pco.points, Status = pco.collapse$Status)

#top.pco <- isoMDS(dist(t(exprs(pco.collapse.top)), method = "manhattan"))
#top.pco.points <- data.frame(top.pco$points)
#colnames(top.pco.points) <- c("PC1", "PC2")
#top.pco.df <- mutate(top.pco.points, Status = pco.collapse.top$Status)

#PCAPlot(orig.pco.df, "Status", "Status", "mds_pco_all")
#PCAPlot(top.pco.df, "Status", "Status", "mds_pco_top")

#orig.pca <- isoMDS(dist(t(exprs(pca.collapse)), method = "manhattan"))
#orig.pca.points <- data.frame(orig.pca$points)
#colnames(orig.pca.points) <- c("PC1", "PC2")
#orig.pca.df <- mutate(orig.pca.points, Status = pca.collapse$Status)

#top.pca <- isoMDS(dist(t(exprs(pca.collapse.top)), method = "manhattan"))
#top.pca.points <- data.frame(top.pca$points)
#colnames(top.pca.points) <- c("PC1", "PC2")
#top.pca.df <- mutate(top.pca.points, Status = pca.collapse.top$Status)

#PCAPlot(orig.pca.df, "Status", "Status", "mds_pca_all")
#PCAPlot(top.pca.df, "Status", "Status", "mds_pca_top")

#orig.cc <- isoMDS(dist(t(exprs(cc.collapse)), method = "manhattan"))
#orig.cc.points <- data.frame(orig.cc$points)
#colnames(orig.cc.points) <- c("PC1", "PC2")
#orig.cc.df <- mutate(orig.cc.points, Status = cc.collapse$Status)

#top.cc <- isoMDS(dist(t(exprs(cc.collapse.top)), method = "manhattan"))
#top.cc.points <- data.frame(top.cc$points)
#colnames(top.cc.points) <- c("PC1", "PC2")
#top.cc.df <- mutate(top.cc.points, Status = cc.collapse.top$Status)

#PCAPlot(orig.cc.df, "Status", "Status", "mds_cc_all")
#PCAPlot(top.cc.df, "Status", "Status", "mds_cc_top")

#kappa.vector <- c(svm.pco.final$results$Kappa, svm.pca.final$results$Kappa, svm.cc.final$results$Kappa)
#kappasd.vector <- c(svm.pco.final$results$KappaSD, svm.pca.final$results$KappaSD, svm.cc.final$results$KappaSD)
#comparison.vector <- c("Patient vs. Control", "Patient vs. Carrier", "Carrier vs. Control") 
#comparison.vector %<>% factor(levels = comparison.vector)
#svm.df <- data.frame(Comparison = comparison.vector, Kappa = kappa.vector, KappaSD = kappasd.vector)

#rfe.plot <- ggplot(svm.df, aes(x = Comparison, y = Kappa, fill = Comparison))  
#rfe.plot <- rfe.plot + geom_col(color = "black") + geom_errorbar(aes(ymax = Kappa + KappaSD, ymin = Kappa - KappaSD), width = 0.25, color = "black")  
#rfe.plot <- rfe.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
#rfe.plot <- rfe.plot + theme(plot.background = element_blank(), legend.position = "none", panel.border = element_rect(color = "black", size = 1))
#rfe.plot <- rfe.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
#CairoPDF("svm.rfe", width = 3, height = 4)
#print(rfe.plot)
#dev.off()


##Model Testing
#pbmc.final <- readRDS.gz("../../pbmc/save/pbmc.collapse.rda")
#pbmc.pca <- pbmc.final[,pbmc.final$Status == "FRDA" | pbmc.final$Status == "Carrier"]
#pbmc.pco <- pbmc.final[,pbmc.final$Status == "FRDA" | pbmc.final$Status == "Normal"]
#pbmc.cc <- pbmc.final[,pbmc.final$Status == "Carrier" | pbmc.final$Status == "Normal"]

#pca.overlap <- intersect(featureNames(pbmc.final), featureNames(pca.collapse.top))
#pco.overlap <- intersect(featureNames(pbmc.final), featureNames(lumi.pco.top))
#cc.overlap <- intersect(featureNames(pbmc.final), featureNames(lumi.cc.top))

#lumi.pca.pbmc <- lumi.pca.top[featureNames(lumi.pca.top) %in% pca.overlap, ]
#pca.pbmc.reduce <- pbmc.pca[featureNames(pbmc.pca) %in% pca.overlap, ]
#lumi.pco.pbmc <- lumi.pco.top[featureNames(lumi.pco.top) %in% pco.overlap, ]
#pco.pbmc.reduce <- pbmc.pco[featureNames(pbmc.pco) %in% pco.overlap, ]
#lumi.cc.pbmc <- lumi.cc.top[featureNames(lumi.cc.top) %in% cc.overlap, ]
#cc.pbmc.reduce <- pbmc.cc[featureNames(pbmc.cc) %in% cc.overlap, ]

#svm.pca.pbmc <- train(x = t(exprs(lumi.pca.pbmc)), y = droplevels(factor(lumi.pca.pbmc$Status)), preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol) 
#svm.pco.pbmc <- train(x = t(exprs(lumi.pco.pbmc)), y = droplevels(factor(lumi.pco.pbmc$Status)), preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol)
#svm.cc.pbmc <- train(x = t(exprs(lumi.cc.pbmc)), y = droplevels(factor(lumi.cc.pbmc$Status)), preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol)

#saveRDS.gz(svm.pca.pbmc, "./save/svm.pca.pbmc.rda")
#saveRDS.gz(svm.pco.pbmc, "./save/svm.pco.pbmc.rda")
#saveRDS.gz(svm.cc.pbmc, "./save/svm.cc.pbmc.rda")

#pca.pbmc.preprocess <- scale(t(exprs(pca.pbmc.reduce)), center = TRUE)
#pco.pbmc.preprocess <- scale(t(exprs(pco.pbmc.reduce)), center = TRUE)
#cc.pbmc.preprocess <- scale(t(exprs(cc.pbmc.reduce)), center = TRUE)

#pbmc.pca.predict <- predict(svm.pca.pbmc, newdata = pca.pbmc.preprocess)
#pbmc.pco.predict <- predict(svm.pco.pbmc, newdata = pco.pbmc.preprocess)
#pbmc.cc.predict <- predict(svm.cc.pbmc, newdata = cc.pbmc.preprocess)

#pbmc.pca.status <- factor(pbmc.pca$Status) %>% revalue(c(FRDA = "Patient"))
#pbmc.pco.status <- factor(pbmc.pco$Status) %>% revalue(c(FRDA = "Patient", Normal = "Control")) %>% factor(levels = c("Control", "Patient"))
#pbmc.cc.status <- factor(pbmc.cc$Status) %>% revalue(c(Normal = "Control"))

#pbmc.pca.confMat <- confusionMatrix(pbmc.pca.predict, pbmc.pca.status)
#pbmc.pco.confMat <- confusionMatrix(pbmc.pco.predict, pbmc.pco.status)
#pbmc.cc.confMat <- confusionMatrix(pbmc.cc.predict, pbmc.cc.status)

##Van Houten
#vanhouten.final <- readRDS.gz("../../VanHauten/baseline/save/intensities.collapse.rda")
#vanhouten.targets <- readRDS.gz("../../VanHauten/baseline/save/targets.final.rda")

#vanhouten.overlap <- intersect(rownames(vanhouten.final), featureNames(lumi.pco.top))
#lumi.pco.vanhouten <- lumi.pco[featureNames(lumi.pco) %in% vanhouten.overlap, ]
#vanhouten.reduce <- vanhouten.final[rownames(vanhouten.final) %in% vanhouten.overlap, ]

#extraTree.pco.vanhouten <- train(x = t(exprs(lumi.pco.vanhouten)), y = droplevels(factor(lumi.pco.vanhouten$Status)), preProcess = c("center", "scale"), method = "extraTrees", ntree = 1000, numThreads = 8, tuneGrid = expand.grid(mtry = floor(sqrt(nrow(lumi.pco.vanhouten))), numRandomCuts = 1), trControl = extraTree.fitcontrol)
#svm.pco.vanhouten <- train(x = t(exprs(lumi.pco.vanhouten)), y = droplevels(factor(lumi.pco.vanhouten$Status)), preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol) #run me

#vanhouten.preprocess <- scale(t(vanhouten.reduce), center = TRUE)
#vanhouten.predict <- predict(svm.pco.vanhouten, vanhouten.preprocess)

#vanhouten.confMat <- confusionMatrix(vanhouten.predict, vanhouten.targets$Status)

#all.samples <- readRDS.gz("../baseline_lumi/save/rmcov.collapse.rda")

#rfe.size <- c(nrow(pca.collapse), RecursiveMultiply(current.value = nrow(pca.collapse)))

#pco.ranking <- sda.ranking(t(exprs(pco.collapse)), pco.collapse$Status, diagonal = FALSE)
#pco.ranking.df <- data.frame(Symbol = names(pco.ranking[,"score"]), Score.pco = signif(pco.ranking[,"score"], 3)) #%>% arrange(desc(Score))
#write.xlsx(pco.ranking.df, "pco.sda.xlsx")

#pco.ranking.list <- ReduceCAT(lumi.object = pco.collapse)
#pco.ranking.median <- map(pco.ranking.list, select, Score) %>% map(unlist) %>% reduce(cbind) %>% apply(1, median)
#pco.ranking.median.df <- data.frame(Symbol = featureNames(pco.collapse), Score = pco.ranking.median) %>% arrange(desc(Score))

#pco.enet <- glmnet(t(exprs(pco.collapse)), pco.collapse$Status, "binomial", alpha = 0.5)
#pco.enet.coef <- coef(pco.enet)[-1,100]

#svm.pco.rfe <- map(rfe.size, SVMCutoffs, pco.ranking.median.df, pco.collapse)
#saveRDS.gz(svm.pco.rfe, "./save/svm.pco.rfe.rda")

#svm.pco.kappa <- map(svm.pco.rfe, getElement, "results") %>% map(getElement, "Kappa") %>% map_dbl(max)
#svm.pco.kappasd <- map(svm.pco.rfe, getElement, "results") %>% map(getElement, "KappaSD") %>% map_dbl(max)
#svm.pco.df <- data.frame(Num.Genes = rfe.size, Kappa = svm.pco.kappa, KappaSD = svm.pco.kappasd)
#svm.pco.df$Num.Genes %<>% factor(levels = sort(rfe.size, decreasing = TRUE))
#svm.pco.df$Comparison <- "Patient vs. Control"

#rfe.plot <- ggplot(svm.pco.df, aes(x = Num.Genes, y = Kappa, group = 1))  
#rfe.plot <- rfe.plot + geom_errorbar(aes(ymax = Kappa + KappaSD, ymin = Kappa - KappaSD), width = 0.25, color = "black") + geom_line() + geom_point()
#rfe.plot <- rfe.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Number of Genes")
#CairoPDF("svm.pco.rfe.accuracies", width = 10, height = 5)
#print(rfe.plot)
#dev.off()

#pca.ranking <- sda.ranking(t(exprs(pca.collapse)), pca.collapse$Status, diagonal = FALSE)
#pca.ranking.df <- data.frame(Symbol = names(pca.ranking[,"score"]), Score.pca = signif(pca.ranking[,"score"], 3)) #%>% arrange(desc(Score))
#write.xlsx(pca.ranking.df, "pca.sda.xlsx")

#pca.ranking.list <- ReduceCAT(lumi.object = pca.collapse)
#pca.ranking.median <- map(pca.ranking.list, select, Score) %>% map(unlist) %>% reduce(cbind) %>% apply(1, median)
#pca.ranking.median.df <- data.frame(Symbol = featureNames(pca.collapse), Score = pca.ranking.median) %>% arrange(desc(Score))

#pca.enet <- glmnet(t(exprs(pca.collapse)), pca.collapse$Status, "binomial", alpha = 0.75)
#pca.enet.coef <- coef(pca.enet)[-1,100]

#pca.ridge <- glmnet(t(exprs(pca.collapse)), pca.collapse$Status, "binomial", alpha = 0)
#pca.ridge.coef <- coef(pca.ridge)[-1,100]

#svm.pca.rfe <- map(rfe.size, SVMCutoffs, pca.ranking.median.df, pca.collapse)
#saveRDS.gz(svm.pca.rfe, "./save/svm.pca.rfe.rda")
#svm.pca.kappa <- map(svm.pca.rfe, getElement, "results") %>% map(getElement, "Kappa") %>% map_dbl(max)
#svm.pca.kappasd <- map(svm.pca.rfe, getElement, "results") %>% map(getElement, "KappaSD") %>% map_dbl(max)
#svm.pca.df <- data.frame(Num.Genes = rfe.size, Kappa = svm.pca.kappa, KappaSD = svm.pca.kappasd)
#svm.pca.df$Num.Genes %<>% factor(levels = sort(rfe.size, decreasing = TRUE))
#svm.pca.df$Comparison <- "Patient vs. Carrier"

#rfe.plot <- ggplot(svm.pca.df, aes(x = Num.Genes, y = Kappa, group = 1))  
#rfe.plot <- rfe.plot + geom_errorbar(aes(ymax = Kappa + KappaSD, ymin = Kappa - KappaSD), width = 0.25, color = "black") + geom_line() + geom_point()
#rfe.plot <- rfe.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Number of Genes") + ylim(c(0,1))
#CairoPDF("svm.pca.rfe.accuracies", width = 10, height = 5)
#print(rfe.plot)
#dev.off()

#cc.ranking <- sda.ranking(t(exprs(cc.collapse)), cc.collapse$Status, diagonal = FALSE)
#cc.ranking.df <- data.frame(Symbol = names(cc.ranking[,"score"]), Score.cc = signif(cc.ranking[,"score"], 3)) #%>% arrange(desc(Score))
#write.xlsx(cc.ranking.df, "cc.sda.xlsx")

#svm.cc.rfe <- map(rfe.size, SVMCutoffs, cc.ranking.df, cc.collapse)
#saveRDS.gz(svm.cc.rfe, "./save/svm.cc.rfe.rda")
#svm.cc.kappa <- map(svm.cc.rfe, getElement, "results") %>% map(getElement, "Kappa") %>% map_dbl(max)
#svm.cc.kappasd <- map(svm.cc.rfe, getElement, "results") %>% map(getElement, "KappaSD") %>% map_dbl(max)
#svm.cc.df <- data.frame(Num.Genes = rfe.size, Kappa = svm.cc.kappa, KappaSD = svm.cc.kappasd)
#svm.cc.df$Num.Genes %<>% factor(levels = sort(rfe.size, decreasing = TRUE))
#svm.cc.df$Comparison <- "Control vs. Carrier"

#rfe.plot <- ggplot(svm.cc.df, aes(x = Num.Genes, y = Kappa, group = 1))  
#rfe.plot <- rfe.plot + geom_errorbar(aes(ymax = Kappa + KappaSD, ymin = Kappa - KappaSD), width = 0.25, color = "black") + geom_line() + geom_point()
#rfe.plot <- rfe.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Number of Genes")
#CairoPDF("svm.cc.rfe.accuracies", width = 10, height = 5)
#print(rfe.plot)
#dev.off()

#svm.all.df <- rbind(svm.pco.df, svm.pca.df, svm.cc.df)
#rfe.plot <- ggplot(svm.all.df, aes(x = Num.Genes, y = Kappa, group = Comparison, color = Comparison))  
#rfe.plot <- rfe.plot + geom_errorbar(aes(ymax = Kappa + KappaSD, ymin = Kappa - KappaSD), width = 0.25, color = "black") + geom_line() + geom_point()
#rfe.plot <- rfe.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + xlab("Number of Genes") + ylim(c(0,1))
#rfe.plot <- rfe.plot + theme(plot.background = element_blank(), legend.background = element_blank(), panel.border = element_rect(color = "black", size = 1))
#rfe.plot <- rfe.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#CairoPDF("svm.rfe.accuracies", width = 10, height = 3.5, bg = "transparent")
#print(rfe.plot)
#dev.off()

#SVMCutoffs <- function(cutoff, top.object, lumi.object) {
    #top.genes <- top.object$Symbol[1:cutoff]
    #lumi.top <- lumi.object[top.genes,]
    #print(dim(lumi.top))
    #svm.linear <- train(x = t(exprs(lumi.top)), y = droplevels(factor(lumi.top$Status)), preProcess = c("center", "scale"), metric = "Kappa", method = "svmLinear", trControl = svm.fitcontrol) #run me
    #return(svm.linear)
#}

##RFE
#RecursiveMultiply <- function(list.nums = vector(), current.value) {
   #new.val = round(0.8 * current.value)
   #if (new.val > 9) {
       #list.nums <- c(list.nums, new.val)
       #RecursiveMultiply(list.nums, new.val)
   #} else {
       #return(list.nums)
   #}
#}

#BootstrapCAT <- function(lumi.object) {
    #boot.sample <- sample.int(ncol(lumi.object), size = ncol(lumi.object), replace = TRUE) 
    #lumi.boot <- lumi.object[,boot.sample]
    #boot.ranking <- sda.ranking(t(exprs(lumi.boot)), lumi.boot$Status, diagonal = FALSE)
    #boot.ranking.df <- data.frame(Symbol = names(boot.ranking[,"score"]), Score = signif(boot.ranking[,"score"], 3)) %>% arrange(Symbol)
    #boot.ranking.df
#}

#ReduceCAT <- function(count.boot = 0, list.boots = list(), lumi.object, boot.count = 100) {
    #count.boot <- count.boot + 1
    #if (count.boot < boot.count) {
        #boot.new <- BootstrapCAT(lumi.object)
        #list.boots <- list.append(list.boots, boot.new)
        #ReduceCAT(count.boot = count.boot, list.boots = list.boots, lumi.object = lumi.object)
    #} else {
        #return(list.boots)
    #}
#}

