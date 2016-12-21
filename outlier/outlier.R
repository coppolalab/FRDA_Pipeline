#Microarray
library(lumi)
library(matrixStats)

#Functional programming
library(magrittr)
library(reshape2)
library(dplyr)
library(tidyr)

#String operations
library(stringr)

#Importing data
library(readr)
library(openxlsx)
library(ggplot2)
library(Cairo)

Threshold <- function(normalized, pdata, threshold) {
    normalized <- data.frame(Status = pdata$Status, t(normalized)) 
    all.status <- summary(pdata$Status)

    normalized.melt <- gather(normalized, Symbol, Z.score, -Status)
    normalized.above <- filter(normalized.melt, Z.score > threshold) %>% dcast(Symbol ~ Status) 
    normalized.below <- filter(normalized.melt, Z.score < -(threshold)) %>% dcast(Symbol ~ Status) 

    normalized.above$Control <- normalized.above$Control / all.status["Control"]
    normalized.above$Carrier <- normalized.above$Carrier / all.status["Carrier"]
    normalized.above$Patient <- normalized.above$Patient / all.status["Patient"]
    normalized.above$Diff.Control <- normalized.above$Patient - normalized.above$Control
    normalized.above$Diff.Carrier <- normalized.above$Carrier - normalized.above$Control
    normalized.above %<>% arrange(desc(Diff.Control), desc(Diff.Carrier))

    normalized.below$Control <- normalized.below$Control / all.status["Control"]
    normalized.below$Carrier <- normalized.below$Carrier / all.status["Carrier"]
    normalized.below$Patient <- normalized.below$Patient / all.status["Patient"]
    normalized.below$Diff.Control <- normalized.below$Patient - normalized.below$Control
    normalized.below$Diff.Carrier <- normalized.below$Carrier - normalized.below$Control
    normalized.below %<>% arrange(desc(Diff.Control), desc(Diff.Carrier))

    Workbook(normalized.above, paste("./normalized.above.", threshold, ".xlsx", sep = ""))
    Workbook(normalized.below, paste("./normalized.below.", threshold, ".xlsx", sep = ""))
    return(list(normalized.above, normalized.below))
}

Workbook <- function(dataset, filename) {
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 1, cols = 1:ncol(dataset), widths = "auto")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

OutlierPlot <- function(gene.symbol, import.lumi, direction, ecdf.fun, num.sd = 3) {
    gene.expr <- import.lumi[match(gene.symbol, featureNames(import.lumi)),]
    gene.plot <- data.frame(Sample.Name = import.lumi$Sample.Name, Status = import.lumi$Status, Expression = as.vector(exprs(gene.expr)))
    gene.plot$Status %<>% factor(levels = c("Patient", "Carrier", "Control"))
    gene.plot %<>% arrange(Status)
    gene.plot$Sample.Name %<>% factor(levels = gene.plot$Sample.Name)

    if (direction == "above")
    {
        intercept.value <- mean(gene.plot$Expression) + (sd(gene.plot$Expression) * num.sd)
    }
    else
    {
        intercept.value <- mean(gene.plot$Expression) - (sd(gene.plot$Expression) * num.sd)
    }

    quantile.value <- ecdf.fun(mean(gene.plot$Expression)) %>% signif(2)

    p <- ggplot(gene.plot, aes(Sample.Name, Expression, col = Status)) + geom_point() + theme_bw()
    p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) 
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + geom_hline(yintercept = mean(gene.plot$Expression), col = "black", size = 2)
    p <- p + geom_hline(yintercept = intercept.value, col = "darkorange", alpha = 0.5, size = 2)
    p <- p + ggtitle(str_c(gene.symbol, " (", quantile.value*100, "% quantile)"))
    CairoPDF(gene.symbol, width = 7, height = 4)
    print(p)
    dev.off()
}

source("../common_functions.R")
import.lumi <- readRDS.gz("../baseline_lumi/save/rmcov.collapse.rda")

means <- rowMeans(exprs(import.lumi))
sds <- rowSds(exprs(import.lumi))
normalized <- sweep(exprs(import.lumi), 1, means) %>% sweep(1, sds, "/")
normalized.3 <- Threshold(normalized, pData(import.lumi), 3)
normalized.4 <- Threshold(normalized, pData(import.lumi), 4)
normalized.5 <- Threshold(normalized, pData(import.lumi), 5)
saveRDS.gz(normalized, file = "./save/normalized.rda")

ecdf.data <- ecdf(means)
OutlierPlot("RPS2", import.lumi, "below", ecdf.data)
OutlierPlot("TRAPPC3", import.lumi, "below", ecdf.data)
OutlierPlot("STON1", import.lumi, "above", ecdf.data)
OutlierPlot("GUCY1B3", import.lumi, "above", ecdf.data)

