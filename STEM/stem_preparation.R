#Functional programming
library(magrittr)
library(purrr)
library(functional)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doBy)

#String operations
library(stringr)

#Reading and writing tables
library(readr)

#For plotting
library(ggplot2)

#For microarray stuff
library(Biobase)
library(matrixStats)

find.baseline <- function(data.vector)
{
    if (min(data.vector) > 1)
    {
        return(TRUE)
    }
    else
    {
        return(FALSE)
    }
}

find.skips <- function(data.vector)
{
    if (sum(data.vector) != length(data.vector))
    {
        return(TRUE)
    }
    else
    {
        return(FALSE)
    }
}

gen.longtable <- function(dataset, targets, timepoints, threshold)
{
    PIDN.long <- filter(targets, Sample.Num == as.character(timepoints)) %>% select(PIDN) %>% as.matrix 
    targets.long <- filter(targets, PIDN %in% PIDN.long) %>% filter(Sample.Num <= timepoints)
    intensities.long <- select(dataset, one_of(targets.long$Sample.Name))
    intensities.long <- rbind(Sample.Num = targets.long$Sample.Num, PIDN = as.numeric(targets.long$PIDN), intensities.long) %>% t %>% data.frame
    intensities.melt <- melt(intensities.long, id.vars = c("Sample.Num", "PIDN"))
    intensities.summed <- dcast(intensities.melt, PIDN + variable ~ Sample.Num, value.var = "value")

    colnames(intensities.summed)[-1] <- c("Gene", paste("Time.", 1:timepoints, sep = ""))
    diffs <- combn(intensities.summed[,-c(1,2)], 2, simplify = FALSE, Compose(as.matrix, rowDiffs, abs)) %>% reduce(cbind)
    colnames(diffs) <- paste("Diff.", 1:ncol(diffs), sep = "")
    threshold.col <- (diffs > threshold) %>% apply(1, reduce, `|`)
    diff.df <- data.frame(diffs, "Threshold" = threshold.col)
    intensities.summed <- data.frame(intensities.summed, diff.df)
    
    intensities.filtered <- filter(intensities.summed, Threshold == TRUE)
    annot.reduce.long <- select(annot.reduce, Probe_Id, Symbol)
    intensities.annot <- merge(intensities.filtered, annot.reduce.long, by.x = "Gene", by.y = "Probe_Id")
    intensities.annot$Symbol[intensities.annot$Symbol == ""] <- "Unknown"
    intensities.annot$Gene %<>% paste(intensities.annot$PIDN, sep = "_")

    intensities.export <- select(intensities.annot, Gene, Symbol, contains("Time"))
    colnames(intensities.export)[1:2] <- c("Spot", "Gene")
    status <- unique(targets$Status)
    filename <- paste("stem_data_", status, "_", timepoints, "_", threshold, ".tsv", sep = "")
    write.table(intensities.export, file = filename, row.names = FALSE, sep = "\t")

    createconfig.call <- paste("python ./create_config.py", status, timepoints, threshold)
    system(createconfig.call)
}

gen.stem <- function(targets.subset, dataset)
{
    gen.longtable(dataset, targets.subset, 3, 0.8)
    gen.longtable(dataset, targets.subset, 3, 0.6)
    gen.longtable(dataset, targets.subset, 4, 0.8)
    gen.longtable(dataset, targets.subset, 4, 0.6)
}

load(file = "../WGCNA/save/intensities.combat.rda")
load(file = "../WGCNA/save/targets.final.known.rda")
load(file = "../WGCNA/save/annot.reduce.rda")

PEER.factors <- read_csv("../WGCNA/factor_10.csv")
PEER.weights <- read_csv("../WGCNA/weight_10.csv")
PEER.factors %<>% select(-(X1:X2), -X5, -X7)
PEER.weights %<>% select(-(X1:X2), -X5, -X7)

PEER.residuals <- as.matrix(PEER.factors) %*% t(as.matrix(PEER.weights)) %>% t
intensities.PEER <- intensities.combat - PEER.residuals
colnames(intensities.PEER) <- colnames(intensities.combat)
rownames(intensities.PEER) <- rownames(intensities.combat)
intensities.PEER %<>% data.frame
save(intensities.PEER, file = "./save/intensities.peer.rda")

targets.final.known$Sample.Num %<>% str_replace("1r", "1")
targets.reduce <- select(targets.final.known, Sample.Name, Sample.Num)

skips <- by(targets.final.known, targets.final.known$PIDN, select, Sample.Num) %>% laply(Compose(as.numeric, sort, diff, find.skips))
baselines <- by(targets.final.known, targets.final.known$PIDN, select, Sample.Num) %>% laply(Compose(as.numeric, find.baseline))
baseline.ids <- sort(unique(targets.final.known$PIDN))[baselines]
skip.ids <- sort(unique(targets.final.known$PIDN))[skips]
all.ids <- unique(c(skip.ids, baseline.ids))

#Replace aberrant sample ranges.  Needs some improvement
targets.fix <- filter(targets.final.known, PIDN %in% all.ids)
fixed.samplenums <- by(targets.fix, targets.fix$PIDN, select, Sample.Num) %>% llply(function(x){ return(1:length(x)) }) %>% reduce(c) %>% as.character #Use a proper lambda!
targets.final.known %<>% arrange(PIDN)
targets.final.known[targets.final.known$PIDN %in% all.ids,]$Sample.Num <- fixed.samplenums
targets.final.known$Sample.Num %<>% as.numeric
save(targets.final.known, file = "./save/targets.final.known.rda")

by(targets.final.known, targets.final.known$Status, gen.stem, intensities.PEER)
system("./stem_bin.sh")

source("../GO/enrichr.R")
enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 
cluster.parallel <- makeForkCluster()

cluster.files <- list.files("./output", full.names = TRUE) %>% str_subset("genetable") 
l_ply(cluster.files, module.split)

module.split <- function(filename)
{
    dataset <- read_tsv(filename)
    folder.name <- str_split(filename, "_")[[1]][3:5] %>% paste(collapse = "_")
    colnames(dataset)[1] <- "Symbol"
    dataset.split <- split(dataset, factor(dataset$Profile))  
    module.names <- names(dataset.split)
    parLapply(cluster.parallel, module.names, submit.enrichr, dataset.split, enrichr.terms, folder.name)
}

submit.enrichr <- function(index, dataset, enrichr.terms, prefix)
{
    data.subset <- dataset[[index]]
    new.fullpath <- file.path("./enrichr", prefix, index)
    dir.create(new.fullpath, showWarnings = FALSE, recursive = TRUE)
    enrichr.data <- lapply(enrichr.terms, get.enrichrdata, data.subset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    lapply(names(enrichr.data), enrichr.wkbk, enrichr.data, new.fullpath)
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

    filename = paste(file.path(new.fullpath, subindex), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

stem.clusters <- read_tsv("./output/stem_defaults_Patient_4_0.6_genetable.tsv")
profile.15 <- filter(stem.clusters, Profile == "15")
profile.15.plot <- select(profile.15, Gene, Time.1, Time.2, Time.3, Time.4) %>% melt
p <- ggplot(profile.15.plot, aes(x = variable, y = value, group = Gene, col = Gene)) + geom_line()
p <- p + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Log2 Normalized Expression") 
CairoPDF("profile15", width = 8, height = 6)
print(p)
dev.off()

profile.22 <- filter(stem.clusters, Profile == "22")
profile.22.plot <- select(profile.22, Gene, Time.1, Time.2, Time.3, Time.4) %>% melt
p <- ggplot(profile.22.plot, aes(x = variable, y = value, group = Gene, col = Gene)) + geom_line()
p <- p + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Log2 Normalized Expression") 
CairoPDF("profile22", width = 8, height = 6)
print(p)
dev.off()

profile.49 <- filter(stem.clusters, Profile == "49")
profile.49.plot <- select(profile.49, Gene, Time.1, Time.2, Time.3, Time.4) %>% melt
p <- ggplot(profile.49.plot, aes(x = variable, y = value, group = Gene, col = Gene)) + geom_line()
p <- p + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Log2 Normalized Expression") 
CairoPDF("profile49", width = 8, height = 6)
print(p)
dev.off()

#intensities.summed <- by(intensities.long, intensities.long$Sample.Num, colMeans) %>% reduce(rbind) %>% t %>% data.frame %>% slice(-1)
#intensities.summed <- by(intensities.long, intensities.long$Sample.Num, Compose(as.matrix, as.vector)) %>% llply(`[`, -c(1:(nrow(intensities.long)/3))) %>% reduce(rbind) %>% t %>% data.frame
#rownames(intensities.summed) <- rownames(intensities.PEER)

#by.gene <- summaryBy(Gene ~ PIDN, intensities.annot, FUN = length)
#by.pidn <- summaryBy(PIDN ~ Symbol, intensities.annot, FUN = length) %>% arrange(PIDN.length)

#p <- ggplot(by.gene, aes(x = Gene.length)) + geom_histogram(stat = "bin", binwidth = 5) + theme_bw()
#p <- ggplot(by.pidn, aes(x = PIDN.length)) + geom_histogram(stat = "bin", binwidth = 1) + theme_bw()
