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
library(openxlsx)

#For plotting
library(ggplot2)

#For microarray stuff
library(Biobase)
library(matrixStats)
library(abind)

library(dtw)
library(fastICA)
library(minerva)
library(minet)
library(parallel)
library(moments)

library(Cairo)
library(WGCNA)
library(Rgraphviz)
library(heatmap.plus)
library(flashClust)
enableWGCNAThreads()

load(file = "../STEM/save/intensities.peer.rda")
load(file = "../STEM/save/targets.final.known.rda")
load(file = "../WGCNA/save/annot.reduce.rda")

gen.long <- function(dataset, targets, timepoints)
{
    PIDN.long <- filter(targets, Sample.Num == as.character(timepoints)) %>% select(PIDN) %>% as.matrix 
    targets.long <- filter(targets, PIDN %in% PIDN.long) %>% filter(Sample.Num <= timepoints)
    intensities.long <- select(dataset, one_of(targets.long$Sample.Name))
    intensities.long %<>% t %>% data.frame
    intensities.long <- data.frame(Sample.Num = targets.long$Sample.Num, PIDN = targets.long$PIDN, Status = targets.long$Status, intensities.long) 
    intensities.melt <- melt(intensities.long, id.vars = c("Sample.Num", "PIDN", "Status"))
    intensities.summed <- dcast(intensities.melt, PIDN + variable + Status ~ Sample.Num, value.var = "value")

    return(intensities.summed)
}

gen.median <- function(dataset)
{
    intensities.median <- by(dataset, factor(dataset$variable), select, -PIDN, -variable) %>% lapply(Compose(as.matrix, colMedians)) %>% reduce(rbind) %>% data.frame
    intensities.median$Probe_Id <- unique(dataset$variable)
    return(intensities.median)
}

gen.mean <- function(dataset)
{
    intensities.mean <- by(dataset, factor(dataset$variable), select, -PIDN, -variable) %>% lapply(Compose(as.matrix, colMeans)) %>% reduce(rbind) %>% data.frame
    intensities.mean$Probe_Id <- unique(dataset$variable)
    return(intensities.mean)
}

gen.dist <- function(dataset)
{
    dtws <- mapply(dtw, split(dataset[[1]], row(dataset[[1]])), split(dataset[[2]], row(dataset[[2]])))
    colnames(dtws) <- rownames(dataset)
    dists <- apply(dtws, 2, `[[`, 'normalizedDistance')
    return(dists)
}

gen.dtw <- function(dataset)
{
    combinations <- combn(dataset, 2, simplify = FALSE) %>% lapply(function(x) {lapply(x, select, -Probe_Id)} )
    dists <- lapply(combinations, gen.dist)
    return(dists)
}

extract.ica <- function(index, dataset)
{
    icm.sig <- which(abs(dataset[index]) > 3)
    dataset.sig <- dataset[icm.sig,] #%>% arrange_(index)
    dataset.sig.sorted <- dataset.sig[order(abs(dataset.sig[,index]), decreasing = TRUE),]
    return(as.numeric(rownames(dataset.sig.sorted)))
}

gen.ica <- function(dataset)
{
    factors <- names(dataset)
    ica.factors <- lapply(factors, extract.ica, dataset)
    return(ica.factors)
}

slice.ica <- function(index, dataset)
{
    return(slice(dataset, index))
}

subset.ica <- function(raw.data, ica.members)
{
   ica.genes <- lapply(ica.members, slice.ica, raw.data) 
   ica.joined <- lapply(ica.genes, merge, annot.reduce)
   ica.joined %<>% lapply(select, Probe_Id, Accession:Definition, contains("X"))
   return(list(ica.joined))
}

enrichr.ica <- function(status.index, data.raw, enrichr.terms, prefix)
{
    data.subset <- data.raw[[status.index]]
    full.path <- file.path("./enrichr", prefix, status.index)
    dir.create(full.path, showWarnings = FALSE, recursive = TRUE)
    names(data.subset) <- paste("X", 1:length(data.subset), sep = "")
    lapply(names(data.subset), enrichr.submit, data.subset, full.path, enrichr.terms, FALSE)
}

enrichr.submit <- function(index, full.df, full.path, enrichr.terms, use.weights)
{
    dataset <- full.df[[index]]
    new.fullpath <- file.path(full.path, index)
    dir.create(new.fullpath, showWarnings = FALSE, recursive = TRUE)
    enrichr.data <- lapply(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    parLapply(cluster.parallel, names(enrichr.data), enrichr.wkbk, enrichr.data, new.fullpath)
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

gen.small.workbook <- function(dataset, filename)
{
    #coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    #colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")
    dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    #conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:3, widths = "auto")
    #setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 4, widths = 45)
    setColWidths(wb, 1, cols = 5:7, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

gen.icaplots <- function(index, dataset)
{
    data.subset <- dataset[[index]]
    names(data.subset) <- paste("ICM", 1:length(data.subset))
    dataset.normalized <- lapply(data.subset, normalize.data)

    dataset.melt <- melt(dataset.normalized, id.vars = "Probe_Id")
    dataset.melt$variable %<>% str_replace("X", "")
    p <- ggplot(dataset.melt, aes(x = variable, y = value, group = Probe_Id, color = Probe_Id)) + geom_line(size = 1.5) 
    p <- p + theme_bw() + theme(legend.position = "none") + ylab("Normalized log2 Expression") + xlab("Time") + facet_wrap(~ L1, scale = "free_y")
    CairoPDF(index, width = 12, height = 8)
    print(p)
    dev.off()
}

normalize.data <- function(dataset)
{
    rownames(dataset) <- dataset$Probe_Id
    dataset %<>% select(-Probe_Id)
    dataset.normalized <- sweep(dataset, 1, dataset$X1)
    dataset.normalized$Probe_Id <- rownames(dataset.normalized) 
    return(dataset.normalized[1:10,])
}

ica.wkbk <- function(index, dataset, new.fullpath)
{
    data.subset <- dataset[[index]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = data.subset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1:3, 5:ncol(data.subset)), widths = "auto")
    setColWidths(wb, 1, cols = 4, widths = 45)

    filename = paste(file.path(new.fullpath, index), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

gen.icatables <- function(index, dataset, prefix)
{
    data.subset <- dataset[[index]]
    full.path <- file.path("./modules", prefix, index)
    dir.create(full.path, showWarnings = FALSE, recursive = TRUE)
    names(data.subset) <- paste("ICM", 1:length(data.subset), sep = "")
    print(str(data.subset))
    l_ply(names(data.subset), ica.wkbk, data.subset, full.path)
}

intensities.3 <- gen.long(data.frame(intensities.PEER), targets.final.known, 3)
intensities.4 <- gen.long(data.frame(intensities.PEER), targets.final.known, 4)
intensities.8 <- gen.long(data.frame(intensities.PEER), targets.final.known, 8)

intensities.3.median <- by(intensities.3, as.factor(intensities.3$Status), select, -Status) %>% lapply(gen.median)
intensities.4.median <- by(intensities.4, as.factor(intensities.4$Status), select, -Status) %>% lapply(gen.median)
intensities.3.mean <- by(intensities.3, as.factor(intensities.3$Status), select, -Status) %>% lapply(gen.mean)
intensities.4.mean <- by(intensities.4, as.factor(intensities.4$Status), select, -Status) %>% lapply(gen.mean)
intensities.8.median <- gen.median(select(intensities.8, -Status))

intensities.3.ica <- lapply(intensities.3.median, select, -Probe_Id) %>% lapply(fastICA, 3) %>% lapply(`[[`, "S") %>% lapply(data.frame) %>% lapply(gen.ica)
intensities.4.ica <- lapply(intensities.4.median, select, -Probe_Id) %>% lapply(fastICA, 4) %>% lapply(`[[`, "S") %>% lapply(data.frame) %>% lapply(gen.ica)
intensities.3.ica.m <- lapply(intensities.3.mean, select, -Probe_Id) %>% lapply(fastICA, 3) %>% lapply(`[[`, "S") %>% lapply(data.frame) %>% lapply(gen.ica)
intensities.4.ica.m <- lapply(intensities.4.mean, select, -Probe_Id) %>% lapply(fastICA, 4) %>% lapply(`[[`, "S") %>% lapply(data.frame) %>% lapply(gen.ica)
intensities.8.ica <- fastICA(select(intensities.8.median, -Probe_Id), 8)[["S"]] %>% data.frame %>% gen.ica

intensities.3.ica.genes <- mapply(intensities.3.median, intensities.3.ica, FUN = subset.ica)
intensities.4.ica.genes <- mapply(intensities.4.median, intensities.4.ica, FUN = subset.ica)
intensities.3.ica.genes.m <- mapply(intensities.3.mean, intensities.3.ica, FUN = subset.ica)
intensities.4.ica.genes.m <- mapply(intensities.4.mean, intensities.4.ica, FUN = subset.ica)
intensities.8.ica.genes <- subset.ica(intensities.8.median, intensities.8.ica)[[1]]

names(intensities.3.ica.genes) <- names(intensities.3.ica)
names(intensities.4.ica.genes) <- names(intensities.4.ica)
names(intensities.3.ica.genes.m) <- names(intensities.3.ica.m)
names(intensities.4.ica.genes.m) <- names(intensities.4.ica.m)

l_ply(names(intensities.4.ica.genes), gen.icaplots, intensities.4.ica.genes)

l_ply(names(intensities.3.ica.genes), gen.icatables, intensities.3.ica.genes, "ica.3.median")
l_ply(names(intensities.4.ica.genes), gen.icatables, intensities.4.ica.genes, "ica.4.median")
l_ply(names(intensities.3.ica.genes.m), gen.icatables, intensities.3.ica.genes.m, "ica.3.mean")
l_ply(names(intensities.4.ica.genes.m), gen.icatables, intensities.4.ica.genes.m, "ica.4.mean")

source("../GO/enrichr.R")
enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 
cluster.parallel <- makeForkCluster()

lapply(names(intensities.3.ica.genes), enrichr.ica, intensities.3.ica.genes, enrichr.terms, "intensities.3.median")
lapply(names(intensities.4.ica.genes), enrichr.ica, intensities.4.ica.genes, enrichr.terms, "intensities.4.median")
lapply(names(intensities.3.ica.genes.m), enrichr.ica, intensities.3.ica.genes.m, enrichr.terms, "intensities.3.mean")
lapply(names(intensities.4.ica.genes.m), enrichr.ica, intensities.4.ica.genes.m, enrichr.terms, "intensities.4.mean")

intensities.4.dists <- gen.dtw(intensities.4.mean) %>% reduce(cbind) %>% data.frame
rownames(intensities.4.dists) <- rownames(intensities.PEER)
colnames(intensities.4.dists) <- coltitles
intensities.4.dists$Probe_Id <- rownames(intensities.4.dists)
dtw.4.joined <- join(intensities.4.dists, annot.reduce)
dtw.4.joined %<>% select(Probe_Id:Definition, Carrier.Control:Control.Patient)
dtw.cc.4 <- arrange(dtw.4.joined, desc(Carrier.Control))
dtw.pca.4 <- select(dtw.4.joined, Probe_Id:Definition, Carrier.Patient, Control.Patient, Carrier.Control) %>% arrange(desc(Carrier.Patient))
dtw.pco.4 <- select(dtw.4.joined, Probe_Id:Definition, Control.Patient, Carrier.Patient, Carrier.Control) %>% arrange(desc(Control.Patient))
gen.small.workbook(dtw.cc.4, "carrier.control.4timepoints.median.xlsx")
gen.small.workbook(dtw.pca.4, "patient.carrier.4timepoints.median.xlsx")
gen.small.workbook(dtw.pco.4, "patient.control.4timepoints.median.xlsx")

intensities.4.patient <- intensities.4.median$Patient
ica.patient <- intensities.4.ica.genes$Patient %>% llply(select, Probe_Id) %>% reduce(c) %>% llply(as.character) %>% reduce(c) %>% unique
dtw.patient <- unique(c(dtw.pca.4[1:2000,]$Probe_Id, dtw.pco.4[1:2000,]$Probe_Id))
all.patient <- unique(c(ica.patient, dtw.patient))

intensities.4.carrier <- intensities.4.median$Carrier
ica.carrier <- intensities.4.ica.genes$Carrier %>% llply(select, Probe_Id) %>% reduce(c) %>% llply(as.character) %>% reduce(c) %>% unique
dtw.carrier <- unique(c(dtw.pca.4[1:2000,]$Probe_Id, dtw.cc.4[1:2000,]$Probe_Id))
all.carrier <- unique(c(ica.carrier, dtw.carrier))

intensities.4.control <- intensities.4.median$Control
ica.control <- intensities.4.ica.genes$Control %>% llply(select, Probe_Id) %>% reduce(c) %>% llply(as.character) %>% reduce(c) %>% unique
dtw.control <- unique(c(dtw.cc.4[1:2000,]$Probe_Id, dtw.pco.4[1:2000,]$Probe_Id))
all.control <- unique(c(ica.control, dtw.control))

genes.list <- list(all.patient, all.carrier, all.control)
names(genes.list) <- c("Patient", "Carrier", "Control")
mi.networks <- lapply(names(genes.list), gen.mi, intensities.4.median, genes.list)
names(mi.networks) <- names(genes.list)
saveRDS.gz(mi.networks, file = "./save/mi.networks.rda")

gen.mi <- function(status, expr, genes)
{
    expr.subset <- expr[[status]]
    genes.subset <- genes[[status]]
    mi.expr <- slice(expr.subset, match(genes.subset, Probe_Id)) 
    rownames(mi.expr) <- mi.expr$Probe_Id
    mi.expr %<>% select(-Probe_Id) %>% t

    mi.mrnet <- minet(mi.expr, method = "mrnet", estimator = "mi.shrink", disc = "equalfreq")
    return(mi.mrnet)
}

gen.mi.tables <- function(status, mi.networks)
{
    mi.connectivity <- colSums(mi.networks[[status]]) %>% sort
    mi.genes <- data.frame(Probe_Id = names(mi.connectivity), Connectivity = mi.connectivity)
    mi.genes %<>% merge(annot.reduce) %>% arrange(desc(Connectivity))
    filename <- paste("mi.", status, ".network.xlsx", sep = "")
    write.xlsx(mi.genes, filename)
}
lapply(names(mi.networks), gen.mi.tables, mi.networks)

normalize.expr <- function(dataset)
{
    PIDNS <- dataset$PIDN
    dataset %<>% select(-Probe_Id, -PIDN, -Status)
    dataset.normalized <- sweep(dataset, 1, dataset$`1`)
    dataset.normalized$PIDN <- PIDNS
    return(dataset.normalized)
}
colnames(intensities.4)[2] <- "Probe_Id"

slc22a7 <- filter(intensities.4, Probe_Id == "ILMN_1653200" & Status == "Patient") %>% normalize.expr %>% melt(id.vars = "PIDN")
p <- ggplot(slc22a7, aes(x = factor(variable), y = value, col = factor(PIDN), group = factor(PIDN))) + geom_line() + facet_wrap(~ PIDN)
p <- p + theme(axis.title.x = element_blank(), legend.position = "none") + ylab("Log2 Normalized Expression")
#p <- p + stat_summary(aes(group = NULL), fun.y = median, color = "red", geom = "point", size = 2)
CairoPDF("slc22a7", width = 21, height = 14)
print(p)
dev.off()

ISCA1L <- filter(intensities.4, Probe_Id == "ILMN_1672024" & Status == "Patient") %>% normalize.expr %>% melt(id.vars = "PIDN")
p <- ggplot(ISCA1L, aes(x = factor(variable), y = value, col = factor(PIDN), group = factor(PIDN))) + geom_line() + facet_wrap(~ PIDN)
#p <- p + stat_summary(aes(group = NULL), fun.y = median, color = "red", geom = "point", size = 2)
p <- p + theme(axis.title.x = element_blank(), legend.position = "none") + ylab("Log2 Normalized Expression")
CairoPDF("isca1l", width = 21, height = 14)
print(p)
dev.off()

atp5ep2 <- filter(intensities.4, Probe_Id == "ILMN_1756674" & Status == "Patient") %>% normalize.expr %>% melt(id.vars = "PIDN")
p <- ggplot(atp5ep2, aes(x = factor(variable), y = value, col = factor(PIDN), group = factor(PIDN))) + geom_line() + facet_wrap(~ PIDN, scale = "free_y")
#p <- p + stat_summary(aes(group = NULL), fun.y = median, color = "red", geom = "point", size = 2)
p <- p + theme(axis.title.x = element_blank(), legend.position = "none") + ylab("Log2 Normalized Expression")
CairoPDF("atp5ep2", width = 21, height = 14)
print(p)
dev.off()

timm22 <- filter(intensities.4, Probe_Id == "ILMN_1765332" & Status == "Patient") %>% normalize.expr %>% melt(id.vars = "PIDN")
p <- ggplot(timm22, aes(x = factor(variable), y = as.numeric(value), col = factor(PIDN), group = factor(PIDN))) + geom_line() + facet_wrap(~ PIDN)
#p <- p + stat_summary(aes(group = NULL), fun.y = median, color = "red", geom = "point", size = 2)
p <- p + theme(axis.title.x = element_blank(), legend.position = "none") + ylab("Log2 Normalized Expression")
CairoPDF("timm22", width = 21, height = 14)
print(p)
dev.off()

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
intensities.3.dists <- gen.dtw(intensities.3.mean) %>% reduce(cbind) %>% data.frame
rownames(intensities.3.dists) <- rownames(intensities.PEER)
colnames(intensities.3.dists) <- coltitles
intensities.3.dists$Probe_Id <- rownames(intensities.3.dists)
dtw.3.joined <- join(intensities.3.dists, annot.reduce)
dtw.3.joined %<>% select(Probe_Id:Definition, Carrier.Control:Control.Patient)
dtw.cc.3 <- arrange(dtw.3.joined, desc(Carrier.Control))
dtw.pca.3 <- select(dtw.3.joined, Probe_Id:Definition, Carrier.Patient, Control.Patient, Carrier.Control) %>% arrange(desc(Carrier.Patient))
dtw.pco.3 <- select(dtw.3.joined, Probe_Id:Definition, Control.Patient, Carrier.Patient, Carrier.Control) %>% arrange(desc(Control.Patient))
gen.small.workbook(dtw.cc.3, "carrier.control.3timepoints.median.xlsx")
gen.small.workbook(dtw.pca.3, "patient.carrier.3timepoints.median.xlsx")
gen.small.workbook(dtw.pco.3, "patient.control.3timepoints.median.xlsx")

intensities.3.dists.m <- gen.dtw(intensities.3.mean) %>% reduce(cbind) %>% data.frame
rownames(intensities.3.dists.m) <- rownames(intensities.PEER)
colnames(intensities.3.dists.m) <- coltitles
intensities.3.dists.m$Probe_Id <- rownames(intensities.3.dists.m)
dtw.3.joined.m <- join(intensities.3.dists.m, annot.reduce)
dtw.3.joined.m %<>% select(Probe_Id:Definition, Carrier.Control:Control.Patient)
dtw.cc.3.m <- arrange(dtw.3.joined.m, desc(Carrier.Control))
dtw.pca.3.m <- select(dtw.3.joined.m, Probe_Id:Definition, Carrier.Patient, Control.Patient, Carrier.Control) %>% arrange(desc(Carrier.Patient))
dtw.pco.3.m <- select(dtw.3.joined.m, Probe_Id:Definition, Control.Patient, Carrier.Patient, Carrier.Control) %>% arrange(desc(Control.Patient))
gen.small.workbook(dtw.cc.3.m, "carrier.control.3timepoints.mean.xlsx")
gen.small.workbook(dtw.pca.3.m, "patient.carrier.3timepoints.mean.xlsx")
gen.small.workbook(dtw.pco.3.m, "patient.control.3timepoints.mean.xlsx")

intensities.4.dists.m <- gen.dtw(intensities.4.mean) %>% reduce(cbind) %>% data.frame
rownames(intensities.4.dists.m) <- rownames(intensities.PEER)
colnames(intensities.4.dists.m) <- coltitles
intensities.4.dists.m$Probe_Id <- rownames(intensities.4.dists.m)
dtw.4.joined.m <- join(intensities.4.dists.m, annot.reduce)
dtw.4.joined.m %<>% select(Probe_Id:Definition, Carrier.Control:Control.Patient)
dtw.cc.4.m <- arrange(dtw.4.joined.m, desc(Carrier.Control))
dtw.pca.4.m <- select(dtw.4.joined.m, Probe_Id:Definition, Carrier.Patient, Control.Patient, Carrier.Control) %>% arrange(desc(Carrier.Patient))
dtw.pco.4.m <- select(dtw.4.joined.m, Probe_Id:Definition, Control.Patient, Carrier.Patient, Carrier.Control) %>% arrange(desc(Control.Patient))
gen.small.workbook(dtw.cc.4.m, "carrier.control.4timepoints.mean.xlsx")
gen.small.workbook(dtw.pca.4.m, "patient.carrier.4timepoints.mean.xlsx")
gen.small.workbook(dtw.pco.4.m, "patient.control.4timepoints.mean.xlsx")

#mi.network.4 <- build.mim(mi.expr.4, estimator = "mi.shrink", disc = "equalfreq")
#mi.dist.4 <- 1 - mi.network.4
#mi.scale.4 <- log2(exp(1)) * (max(mi.dist.4) - min(mi.dist.4)) * mi.dist.4
#mi.mrnet.4 <- mrnet(mi.scale.4)

#saveRDS.gz(mi.mrnet.4, file = "./save/mi.mrnet.4.rda")
#mi.network.aracne.4 <- minet(ica.expr.4, method = "aracne", estimator = "mi.shrink", disc = "equalfreq")
#mi.network.clr.4 <- minet(ica.expr.4, method = "clr", estimator = "mi.shrink", disc = "equalfreq")
#saveRDS.gz(mi.network.mrnet.4, file = "./save/mi.network.mrnet.4")
#saveRDS.gz(mi.network.aracne.4, file = "./save/mi.network.aracne.4")
#saveRDS.gz(mi.network.clr.4, file = "./save/mi.network.clr.4")

#mi.network.mine.4 <- mine(as.matrix(mi.expr.4), master = 1:ncol(mi.expr.4), alpha = 1)
#mi.network.mine.mic.4 <- mi.network.mine.4$MIC
#saveRDS.gz(mi.network.mine.4, file = "./save/mi.network.mine.4")
#saveRDS.gz(mi.network.mine.mic.4, file = "./save/mi.network.mine.mic.4")

#mine.dist.4 <- 1 - mi.network.mine.mic.4
#mine.scale.4 <- log2(exp(1)) * mine.dist.4
#mine.mrnet.4 <- mrnet(mine.scale.4)
#mine.genes <- colSums(mine.mrnet.4) #%>% sort #%>% Filter(f = Curry(`>`, e2 = 1))
##Hub Genes
#mine.genes.hub <- mine.genes[mine.genes > 5] %>% names

##SLC25A3
#mine.genes.net <- mine.mrnet.4[,mine.genes.hub] %>% data.frame
#mine.genes.net$Probe_Id <- rownames(mine.genes.net)
#test.matrix <- filter(mine.genes.net, ILMN_2332713 > 1)
#test.joined <- merge(test.matrix, annot.reduce)
#write.xlsx(test.joined, "test1.xlsx")

#saveRDS.gz(mine.mrnet.4, file = "./save/mine.mrnet.4")
#mine.clr.4 <- clr(mi.network.mine.mic.4)
#saveRDS.gz(mine.clr.4, file = "./save/mine.clr.4")
#mine.aracne.4 <- aracne(mi.network.mine.mic.4)
#saveRDS.gz(mine.aracne.4, file = "./save/mine.aracne.4")

#intensities.carrier <- intensities.4.median[["Carrier"]]
#ica.genes.carrier <- intensities.4.ica.genes[["Carrier"]] %>% llply(select, Probe_Id) %>% reduce(c) %>% llply(as.character) %>% reduce(c) %>% unique
#dtw.genes.carrier <- unique(c(dtw.pca.4[1:300,]$Probe_Id, dtw.cc.4[1:300,]$Probe_Id))

#all.genes.carrier <- unique(c(ica.genes.carrier, dtw.genes.carrier)) #%>% paste(collapse = "|")

#mi.expr.carrier <- slice(intensities.carrier, match(all.genes.carrier, Probe_Id)) #%>% join(annot.reduce)
#rownames(mi.expr.carrier) <- mi.expr.carrier$Probe_Id
##mi.expr.carrier %<>% select(-(Probe_Id:Definition)) %>% t
#mi.expr.carrier %<>% select(-Probe_Id) %>% t
#saveRDS.gz(mi.expr.carrier, file = "./save/mi.expr.4.rda")

#mi.network.mine.carrier <- mine(mi.expr.carrier, master = 1:ncol(mi.expr.carrier), alpha = 1)
#mi.network.mine.mic.carrier <- mi.network.mine.carrier$MIC
#saveRDS.gz(mi.network.mine.carrier, file = "./save/mi.network.mine.carrier")
#saveRDS.gz(mi.network.mine.mic.carrier, file = "./save/mi.network.mine.mic.carrier")

#mine.dist.carrier <- 1 - mi.network.mine.mic.carrier
#mine.scale.carrier <- log2(exp(1)) * mine.dist.carrier
#mine.mrnet.carrier <- mrnet(mine.scale.carrier)

#intensities.control <- intensities.4.median[["Control"]]
#ica.genes.control <- intensities.4.ica.genes[["Control"]] %>% llply(select, Probe_Id) %>% reduce(c) %>% llply(as.character) %>% reduce(c) %>% unique
#dtw.genes.control <- unique(c(dtw.pco.4[1:300,]$Probe_Id, dtw.cc.4[1:300,]$Probe_Id))

#all.genes.control <- unique(c(ica.genes.control, dtw.genes.control)) #%>% paste(collapse = "|")

#mi.expr.control <- slice(intensities.control, match(all.genes.control, Probe_Id)) #%>% join(annot.reduce)
#rownames(mi.expr.control) <- mi.expr.control$Probe_Id
##mi.expr.control %<>% select(-(Probe_Id:Definition)) %>% t
#mi.expr.control %<>% select(-Probe_Id) %>% t
#saveRDS.gz(mi.expr.control, file = "./save/mi.expr.4.rda")

#mi.network.mine.control <- mine(mi.expr.control, master = 1:ncol(mi.expr.control), alpha = 1)
#mi.network.mine.mic.control <- mi.network.mine.control$MIC
#saveRDS.gz(mi.network.mine.control, file = "./save/mi.network.mine.control")
#saveRDS.gz(mi.network.mine.mic.control, file = "./save/mi.network.mine.mic.control")

#mine.dist.control <- 1 - mi.network.mine.mic.control
#mine.scale.control <- log2(exp(1)) * mine.dist.control
#mine.mrnet.control <- mrnet(mine.scale.control)
#intensities.8.median$Probe_Id %<>% as.character
#ica.genes <- intensities.8.ica.genes %>% llply(select, Probe_Id) %>% reduce(c) %>% llply(as.character) %>% reduce(c) %>% unique
#dtw.genes <- unique(c(dtw.pca.4[1:300,]$Probe_Id, dtw.pco.4[1:300,]$Probe_Id))

#skewed.genes <- apply(select(intensities.8.median, -Probe_Id), 1, anscombe.test) %>% llply(`[`, "p.value") %>% reduce(c) %>% reduce(c)
#skewed.genes.index <- which(skewed.genes < 0.05)
#anscombe.genes <- slice(intensities.8.median, skewed.genes.index)$Probe_Id

#skewed.genes.j <- apply(select(intensities.8.median, -Probe_Id), 1, jarque.test) %>% llply(`[`, "p.value") %>% reduce(c) %>% reduce(c)
#skewed.genes.j.index <- which(skewed.genes.j < 0.05)
#jarque.genes <- slice(intensities.8.median, skewed.genes.j.index)$Probe_Id
#all.genes <- unique(c(ica.genes, dtw.genes, anscombe.genes, jarque.genes)) #%>% paste(collapse = "|")

#ica.expr.8 <- slice(intensities.8.median, match(all.genes, Probe_Id)) %>% join(annot.reduce)
#rownames(ica.expr) <- ica.expr$Probe_Id
#ica.expr.8 %<>% select(-(Probe_Id:Definition)) %>% t
#save(ica.expr.8, file = "./save/icr.expr.rda")

#mi.network.mrnet <- minet(ica.expr, method = "mrnet", estimator = "mi.shrink", disc = "equalfreq")
#mi.network.aracne <- minet(ica.expr, method = "aracne", estimator = "mi.shrink", disc = "equalfreq")
#mi.network.clr <- minet(ica.expr, method = "clr", estimator = "mi.shrink", disc = "equalfreq")
#saveRDS.gz(mi.network.mrnet, file = "./save/mi.network.mrnet")
#saveRDS.gz(mi.network.aracne, file = "./save/mi.network.aracne")
#saveRDS.gz(mi.network.clr, file = "./save/mi.network.clr")
#CairoPDF("test.pdf", width = 6, height = 6)
#plot(as(mi.network.mrnet, "graphNEL"))
#dev.off()

#mi.network.mine <- mine(as.matrix(ica.expr), master = 1:ncol(ica.expr), alpha = 1)
#mi.network.mine.mic <- mi.network.mine$MIC
#saveRDS.gz(mi.network.mine, file = "./save/mi.network.mine")
#saveRDS.gz(mi.network.mine.mic, file = "./save/mi.network.mine.mic")
#mine.mrnet <- mrnet(mi.network.mine.mic)
#saveRDS.gz(mine.mrnet, file = "./save/mine.mrnet")
#mine.clr <- clr(mi.network.mine.mic)
#saveRDS.gz(mine.clr, file = "./save/mine.clr")
#mine.aracne <- aracne(mi.network.mine.mic)
#saveRDS.gz(mine.aracne, file = "./save/mine.aracne")
