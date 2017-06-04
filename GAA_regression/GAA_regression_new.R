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

source("../common_functions.R")

GAADensityPlot <- function(gaa.vector, filename) {
    gaa.density <- density(gaa.vector) #compute density kernel of GAA distribution
    density.df <- data.frame(x = gaa.density$x, y = gaa.density$y) #extract density coordinates and put in data frame
    gaa.summary <- summary(gaa.vector) %>% tidy #get summary info
    gaa.quartiles <- unlist(gaa.summary[c(2,3,5)]) #extract quartiles 1 and 3 and median
    density.df$quant <- factor(findInterval(density.df$x, gaa.quartiles)) #assign each x-coordinate of density to the correct quartile
    names(gaa.quartiles) <- gaa.quartiles #name the quartiles vector the same as its contents, because ggplot is weird like that

    #Explanation of this plot: by supplying x and y, geom_line prints the lines at the correct vertical heights and x locations - used instead of vline
    #geom_ribbon is like density except it seems to support having multiple colors 
    
    p <- ggplot(density.df, aes(x, y)) + geom_line() + geom_ribbon(aes(ymin = 0, ymax = y, fill = quant)) 
    p <- p + scale_x_continuous(breaks = gaa.quartiles) + scale_fill_brewer(palette = "Set1", guide = "none")
    p <- p + xlab("Repeat Length") + ylab("Density")
    CairoPDF(filename, width = 6, height = 6)
    print(p)
    dev.off()
}

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

lumi.import <- ReadRDSgz("../WGCNA_GAA/save/export.lumi.rda")
#lumi.import <- lumi.import[,lumi.import$FDS != 6 & lumi.import$FDS != 1]

all.model <- model.matrix( ~ Age + RIN, pData(lumi.import))[,-1] %>% data.frame
all.rmcov <- removeBatchEffect(lumi.import, covariates = all.model) %>% t
#all.rmcov <- empiricalBayesLM(t(exprs(lumi.import)), removedCovariates = all.model) %>%
    #extract2("adjustedData") 

#SVM
svm.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, savePredictions = TRUE)
ranger.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, allowParallel = FALSE, savePredictions = TRUE)

set.seed(12345)
gaa.enet <- glmnet(all.rmcov, as.vector(scale(lumi.import$FDS)), "gaussian", nlambda = 10000, alpha = 0.5, type.gaussian = "naive")
gaa.enet.coef <- coef(gaa.enet)[-1,ncol(gaa.enet$beta)]
gaa.enet.df <- tibble(Symbol = names(gaa.enet.coef), Coefficient = signif(gaa.enet.coef, 3)) %>% arrange(desc(abs(Coefficient))
full.predict <- predict(gaa.enet, all.rmcov)
predict.plot <- data.frame(Original = scale(lumi.import$FDS), Predicted = full.predict[,ncol(full.predict)])
full.rmse <- rmse(full.predict[,ncol(full.predict)], scale(lumi.import$FDS))
lm(full.predict[,ncol(full.predict)] ~ scale(lumi.import$FDS)) %>% summary %>% str

p <- ggplot(predict.plot, aes(x = Original, y = Predicted)) + 
    geom_jitter() + 
    stat_smooth(method = "loess") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

CairoPDF("glm.pred.original", height = 4, width = 4)
print(p)
dev.off()

FitGLM <- function(expr.matrix, sample.index, y.vector) {
    boot.sample <- expr.matrix[sample.index,]
    print("Fitting penalized regression model...")
    gaa.enet <- glmnet(boot.sample, as.vector(scale(y.vector)), "gaussian", nlambda = 10000, alpha = 0.5, type.gaussian = "naive")
    print("..model finished")
    predict.boot <- predict(gaa.enet, boot.sample)
    boot.rmse <- rmse(predict.boot[,ncol(predict.boot)], scale(y.vector))
    boot.rsquared <- lm(predict.boot[,ncol(predict.boot)] ~ scale(y.vector)) %>% 
        summary %>% extract2("r.squared")
    c(RMSE = boot.rmse, Rsquared = boot.rsquared)
}

boot.test <- boot(all.rmcov, FitGLM, R = 25, y.vector = lumi.import$FDS)
all.rmcov.df <- data.frame(FDS.scaled = scale(lumi.import$FDS), all.rmcov)
h2o.test <- h2o.deeplearning(y = "FDS.scaled", training_frame = all.rmcov.df)

lumi.top <- all.rmcov[,gaa.enet.coef != 0]
set.seed(12345)
svm.train <- train(x = lumi.top, y = as.vector(scale(lumi.import$FDS)), 
                   preProcess = c("center", "scale"), metric = "RMSE", 
                   method = "svmPoly", trControl = svm.fitcontrol) 
glm.train <- train(x = all.rmcov, y = as.vector(scale(lumi.import$FDS)), 
                   preProcess = c("center", "scale"), metric = "RMSE", 
                   method = "glmnet", family = "gaussian", nlambda = 10000,
                   tuneGrid = data.frame(alpha = 0.5, lambda = 0.0163), trControl = svm.fitcontrol) 
#ranger.train <- train(x = all.rmcov, y = lumi.import$FDS, preProcess = c("center", "scale"), metric = "RMSE", method = "ranger", trControl = ranger.fitcontrol, num.trees = 1000, tuneGrid = data.frame(mtry = nrow(lumi.top)))
svm.preds <- svm.train$pred
glm.preds <- glm.train$pred

temp <- predict(glm.train$finalModel, all.rmcov)
full.rmse <- rmse(temp[,ncol(temp)], scale(lumi.import$FDS))
lm(temp[,ncol(temp)] ~ scale(lumi.import$FDS)) %>% summary

p <- ggplot(svm.preds, aes(x = obs, y = pred)) + 
    geom_jitter() + 
    facet_wrap( ~ Resample) + 
    stat_smooth(method = "loess") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

CairoPDF("svm.poly.preds", height = 30, width = 40)
print(p)
dev.off()

p <- ggplot(glm.preds, aes(x = obs, y = pred)) + 
    geom_jitter() + 
    facet_wrap( ~ Resample) + 
    stat_smooth(method = "loess") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

CairoPDF("glm.preds", height = 30, width = 40)
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

enet.positive <- filter(gaa.enet.df, Coefficient > 0)
enet.negative <- filter(gaa.enet.df, Coefficient < 0)

svm.top5.positive <- lumi.cleaned[enet.positive$Symbol[1:5],]
svm.top5.negative <- lumi.cleaned[enet.negative$Symbol[1:5],]

svm.top5.positive <- data.frame(GAA1 = lumi.cleaned$GAA1, t(exprs(lumi.cleaned[enet.positive$Symbol[1:5],]))) %>% gather(Gene, Expression, -GAA1)
svm.top5.positive$Gene %<>% factor(levels = str_replace(enet.positive$Symbol[1:5], "-", "."))
top5.plot <- ggplot(svm.top5.positive, aes(x = Expression, y = GAA1)) + geom_point() + stat_smooth()
top5.plot <- top5.plot + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
top5.plot <- top5.plot + facet_wrap( ~ Gene, ncol = 5, scales = "free" )
CairoPDF("svm.top5.positive", width = 20, height = 4)
print(top5.plot)
dev.off()

svm.top5.negative <- data.frame(GAA1 = lumi.cleaned$GAA1, t(exprs(lumi.cleaned[enet.negative$Symbol[1:5],]))) %>% gather(Gene, Expression, -GAA1)
svm.top5.negative$Gene %<>% factor(levels = str_replace(enet.negative$Symbol[1:5], "-", "."))
top5.plot <- ggplot(svm.top5.negative, aes(x = Expression, y = GAA1)) + geom_point() + stat_smooth()
top5.plot <- top5.plot + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
top5.plot <- top5.plot + facet_wrap( ~ Gene, ncol = 5, scales = "free" )
CairoPDF("svm.top5.negative", width = 20, height = 4)
print(top5.plot)
dev.off()

regressed.genes <- read.xlsx("../WGCNA_GAA/ebam.annot.xlsx")
regressed.sig <- filter(regressed.genes, Posterior > 0.9)

bf.genes <- read.xlsx("../WGCNA_GAA/bf.fds.xlsx") %>% as_tibble
bf.sig <- filter(bf.genes, Log.Bayes.Factor > 0.5 & Symbol %in% colnames(lumi.top)) %>% 
    arrange(desc(abs(Median)))

intersect(regressed.sig$Symbol, featureNames(lumi.top))
intersect(bf.sig$Symbol, colnames(lumi.top))

STOP
#registerDoMC(cores = 6)
#svm.continuous.rfe <- map(rfe.size, SVMContinuous, care.ranking.median.df, lumi.cleaned)
#SaveRDSgz(svm.continuous.rfe, "./save/svm.continuous.rfe.rda")

#svm.continuous.cost <- map(svm.continuous.rfe, extract2, "finalModel") %>% map(attr, "param") %>% map_dbl(extract2, "C") 
#svm.continuous.degree <- map(svm.continuous.rfe, extract2, "finalModel") %>% map(attr, "kernelf") %>% map(attr, "kpar") %>% map_dbl(extract2, "degree")
#svm.continuous.scale <- map(svm.continuous.rfe, extract2, "finalModel") %>% map(attr, "kernelf") %>% map(attr, "kpar") %>% map_dbl(extract2, "scale")

#svm.continuous.results <- map(svm.continuous.rfe, extract2, "results")
#svm.continuous.filtered <- map2(svm.continuous.results, svm.continuous.cost, FilterResults, "C") %>% 
    #map2(svm.continuous.degree, FilterResults, "degree") %>% 
    #map2(svm.continuous.scale, FilterResults, "scale")

#svm.continuous.rmse <- map_dbl(svm.continuous.filtered, extract2, "RMSE")
#svm.continuous.rmsesd <- map_dbl(svm.continuous.filtered, extract2, "RMSESD")
#svm.continuous.rsquared <- map_dbl(svm.continuous.filtered, extract2, "Rsquared")
#svm.continuous.rsquaredsd <- map_dbl(svm.continuous.filtered, extract2, "RsquaredSD")

#svm.continuous.df <- data.frame(Num.Genes = rfe.size, RMSE = svm.continuous.rmse, RMSE.SD = svm.continuous.rmsesd, Rsquared = svm.continuous.rsquared, Rsquared.SD = svm.continuous.rsquaredsd)
#svm.continuous.df$Num.Genes %<>% factor(levels = sort(rfe.size, decreasing = TRUE))

#rmse.plot <- ggplot(svm.continuous.df, aes(x = Num.Genes, y = RMSE, group = 1))  
#rmse.plot <- rmse.plot + geom_errorbar(aes(ymax = RMSE + RMSE.SD, ymin = RMSE - RMSE.SD), width = 0.25, color = "black") + geom_line() + geom_point()
#rmse.plot <- rmse.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Number of Genes") + ylim(c(0,400))
#rmse.plot <- rmse.plot + theme(plot.background = element_blank())
#CairoPDF("svm.rmse", width = 8, height = 3.5, bg = "transparent")
#print(rmse.plot)
#dev.off()

#rsquared.plot <- ggplot(svm.continuous.df, aes(x = Num.Genes, y = Rsquared, group = 1))  
#rsquared.plot <- rsquared.plot + geom_errorbar(aes(ymax = Rsquared + Rsquared.SD, ymin = Rsquared - Rsquared.SD), width = 0.25, color = "black") + geom_line() + geom_point()
#rsquared.plot <- rsquared.plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Number of Genes") + ylim(c(0,1))
#rsquared.plot <- rsquared.plot + theme(plot.background = element_blank(), panel.border = element_rect(color = "black", size = 1))
#rsquared.plot <- rsquared.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#rsquared.plot <- rsquared.plot + ylab(expression(R^2))
#CairoPDF("svm.rsquared", width = 8, height = 3.5, bg = "transparent")
#print(rsquared.plot)
#dev.off()

#care.positive <- arrange(care.ranking.median.df, desc(Score)) 
#care.negative <- arrange(care.ranking.median.df, Score)

#RecursiveMultiply <- function(current.value, list.nums = vector()) {
   #new.val = round(0.8 * current.value)
   #if (new.val > 9)
   #{
       #list.nums <- c(list.nums, new.val)
       #RecursiveMultiply(new.val, list.nums)
   #}
   #else
   #{
       #return(list.nums)
   #}
#}

#SVMContinuous <- function(cutoff, top.object, lumi.object) {
    #top.genes <- top.object$Symbol[1:cutoff]
    #lumi.top <- lumi.object[top.genes,]
    #print(dim(lumi.top))
    #svm.fitcontrol <- trainControl(method = "boot632", verboseIter = TRUE, number = 25, savePredictions = TRUE, allowParallel = FALSE)
    #svm.train <- train(x = t(exprs(lumi.top)), y = lumi.top$GAA1, preProcess = c("center", "scale"), metric = "RMSE", method = "svmPoly", trControl = svm.fitcontrol)  
    #return(svm.train)
#}

#FilterResults <- function(results.df, value, column) {
    #filtered.df <- filter_(results.df, str_c(column, "==", value))
    #filtered.df
#}

#BootstrapCAR <- function(lumi.object) {
    #boot.sample <- sample.int(ncol(lumi.object), size = ncol(lumi.object), replace = TRUE) 
    #lumi.boot <- lumi.object[,boot.sample]
    #care.ranking <- carscore(t(exprs(lumi.boot)), lumi.boot$GAA1, lambda = 0, verbose = TRUE)
    #boot.ranking.df <- data.frame(Symbol = names(care.ranking), Score = as.vector(care.ranking)) %>% arrange(Symbol)
    #boot.ranking.df
#}

#ReduceCAR <- function(count.boot = 0, list.boots = list(), lumi.object, boot.count = 1000) {
    #if (count.boot < boot.count) {
        #count.boot <- count.boot + 1
        #boot.new <- BootstrapCAR(lumi.object)
        #list.boots <- list.append(list.boots, boot.new)
        #ReduceCAR(count.boot = count.boot, list.boots = list.boots, lumi.object = lumi.object)
    #} else {
        #return(list.boots)
    #}
#}

#care.ranking <- carscore(t(exprs(lumi.cleaned)), lumi.cleaned$GAA1, lambda = 0, verbose = TRUE)
#care.df <- data.frame(Symbol = names(care.ranking), Score = signif(as.vector(care.ranking), 3)) %>% arrange(desc(abs(Score)))

#care.ranking.list <- ReduceCAR(lumi.object = lumi.cleaned)
#care.ranking.median <- map(care.ranking.list, select, Score) %>% map(unlist) %>% reduce(cbind) %>% apply(1, median)
#care.ranking.median.df <- data.frame(Symbol = featureNames(lumi.cleaned), Score = care.ranking.median) %>% arrange(desc(abs(Score)))

#svm.continuous.final <- svm.continuous.rfe[[17]]
#care.df.top <- slice(care.df, 1:as.numeric(as.character(svm.continuous.df$Num.Genes[17])))
#rfe.size <- c(nrow(lumi.cleaned), RecursiveMultiply(nrow(lumi.cleaned)))
