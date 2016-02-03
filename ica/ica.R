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
library(lumi)
library(lumiHumanAll.db)
library(annotate)

#Longitudinal Analysis
library(fastICA)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(vadr)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doBy)

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

extract.ica <- function(index, dataset)
{
    icm.sig <- which(abs(dataset[index]) > 3)
    dataset.sig <- dataset[icm.sig,] #%>% arrange_(index)
    dataset.sig.sorted <- dataset.sig[order(abs(dataset.sig[,index]), decreasing = TRUE),]
    return(as.numeric(rownames(dataset.sig.sorted)))
}

gen.ica <- function(dataset)
{
    dataset.abs <- abs(dataset) 
    ica.subset <- apply((dataset.abs > 3), 1, any)
    return(dataset.abs[ica.subset,])
}

split.ica <- function(data.vector)
{
    dataset <- data.frame(Symbol = names(data.vector), Module = data.vector)
    dataset.split <- split(dataset, dataset$Module) 
    names(dataset.split) <- paste("X", 1:length(dataset.split), sep = "")
    return(dataset.split)
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

gen.icaplots <- function(index, dataset)
{
    data.subset <- dataset[[index]]
    names(data.subset) <- paste("ICM", 1:length(data.subset))
    dataset.normalized <- lapply(data.subset, normalize.data)

    dataset.melt <- melt(dataset.normalized, id.vars = "Symbol")
    dataset.melt$variable %<>% str_replace("X", "")
    p <- ggplot(dataset.melt, aes(x = variable, y = value, group = Symbol, color = Symbol)) + geom_line(size = 1.5) 
    p <- p + theme_bw() + theme(legend.position = "none") + ylab("Normalized log2 Expression") + xlab("Time") + facet_wrap(~ L1, scale = "free_y")
    CairoPDF(index, width = 12, height = 8)
    print(p)
    dev.off()
}

normalize.data <- function(dataset)
{
    rownames(dataset) <- dataset$Symbol
    dataset %<>% select(-Symbol)
    dataset.normalized <- sweep(dataset, 1, dataset$X1)
    dataset.normalized$Symbol <- rownames(dataset.normalized) 
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

stringdb.ica <- function(status.index, data.raw, prefix)
{
    data.subset <- data.raw[[status.index]]
    full.path <- file.path("./stringdb", prefix, status.index)
    dir.create(full.path, showWarnings = FALSE, recursive = TRUE)
    names(data.subset) <- paste("X", 1:length(data.subset), sep = "")
    map(names(data.subset), stringdb.submit, data.subset, full.path)
}

stringdb.submit <- function(component, all.components, full.path)
{
    dataset <- all.components[[component]]
    get.stringdb(dataset, component, full.path, 400)
}

setnames <- function(dataset, symbols.vector)
{
    rownames(dataset) <- symbols.vector
    return(dataset)
}

seed.ICA <- function(intensities, ica.list = list(), iter.count = 0)
{
   if (iter.count < 250) 
   {
        ica.new <- map(intensities, select, -Symbol) %>% map(fastICA, 4) %>% map(`[[`, "S") %>% map(data.frame) 
        iter.count <- iter.count + 1
        ica.list <- c(ica.list, ica.new)
        seed.ICA(intensities, ica.list, iter.count)
   }
   else
   {
        return(ica.list)
   }
}

ica.all <- seed.ICA
intensities.means <- readRDS.gz("../dtw/save/intensities.means.rda")

intensities.ica <- map(intensities.means, select, -Symbol) %>% map(fastICA, 4) %>% map(`[[`, "S") %>% map(data.frame) 
intensities.ica %<>% map(setnames, intensities.means$Patient$Symbol) %>% map(gen.ica)

intensities.ica.genes <- map(intensities.ica, apply, 1, which.max) %>% map(split.ica)

#l_ply(names(intensities.ica.genes), gen.icatables, intensities.ica.genes, "ica.mean")

source("../../code/GO/enrichr.R")
enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 

map(names(intensities.ica.genes), enrichr.ica, intensities.ica.genes, enrichr.terms, "intensities.mean")
map(names(intensities.ica.genes), stringdb.ica, intensities.ica.genes, "means")
saveRDS.gz(intensities.ica.genes, "./save/intensities.ica.genes.rda")
