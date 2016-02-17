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

#Network Analysis
library(minerva)
library(minet)
library(parallel)
library(moments)
library(Rgraphviz)
library(flashClust)
library(WGCNA)
enableWGCNAThreads()

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

library(Cairo)

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

intensities.means <- readRDS.gz("../dtw/save/intensities.means.rda")
intensities.ica.genes <- readRDS.gz("../ica/save/intensities.ica.genes.rda")
dtw.cc <- readRDS.gz("../dtw/save/dtw.cc.rda")
dtw.pca <- readRDS.gz("../dtw/save/dtw.pca.rda")
dtw.pco <- readRDS.gz("../dtw/save/dtw.pco.rda")
fdata <- readRDS.gz("../dtw/save/fdata.rda")

genes <- list()
intensities.patient <- intensities.means$Patient
ica.patient <- intensities.ica.genes$Patient %>% map(select, Symbol) %>% map(Compose(unlist, as.character)) %>% reduce(c) %>% unique
dtw.patient <- unique(c(dtw.pca[1:2000,]$Symbol, dtw.pco[1:2000,]$Symbol)) %>% str_replace("\\-", "\\.")
genes$Patient <- unique(c(ica.patient, dtw.patient))
saveRDS.gz(all.patient, "./save/all.patient.rda")

intensities.carrier <- intensities.means$Carrier
ica.carrier <- intensities.ica.genes$Carrier %>% map(select, Symbol) %>% map(Compose(unlist, as.character)) %>% reduce(c) %>% unique
dtw.carrier <- unique(c(dtw.pca[1:2000,]$Symbol, dtw.cc[1:2000,]$Symbol)) %>% str_replace("\\-", "\\.")
genes$Carrier <- unique(c(ica.carrier, dtw.carrier))
saveRDS.gz(all.carrier, "./save/all.carrier.rda")

intensities.control <- intensities.means$Control
ica.control <- intensities.ica.genes$Control %>% map(select, Symbol) %>% map(Compose(unlist, as.character)) %>% reduce(c) %>% unique
dtw.control <- unique(c(dtw.cc[1:2000,]$Symbol, dtw.pco[1:2000,]$Symbol)) %>% str_replace("\\-", "\\.")
genes$Control <- unique(c(ica.control, dtw.control))
saveRDS.gz(all.control, "./save/all.control.rda")

gen.mi <- function(status, expr, genes)
{
    expr.subset <- expr[[status]]
    genes.subset <- genes[[status]]
    
    mi.expr.df <- slice(expr.subset, match(genes.subset, Probe_Id)) 
    print(dim(mi.expr.df))
    rownames(mi.expr.df) <- mi.expr.df$Probe_Id
    mi.expr <- select(mi.expr.df, -Probe_Id) %>% t

    mi.mrnet <- minet(mi.expr, method = "mrnet", estimator = "mi.shrink", disc = "equalwidth")
    TOM.PEER <- TOMsimilarity(mi.mrnet, verbose = 5)
    dissimilarity.TOM <- 1 - TOM.PEER

    CairoPDF(paste("scalefree", status, sep = "_"), height = 6, width = 6)
    scaleFreePlot(TOM.PEER) #Meh
    dev.off()

    geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")

    CairoPDF(file = paste("./genecluster", status, sep = "_"), height = 10, width = 15)
    plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
    dev.off()

    min.module.size <- 50
    #Identify modules using dynamic tree cutting with hybrid clustering
    dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, cutHeight = 0.995, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
    dynamic.colors <- labels2colors(dynamic.modules)
    saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

    CairoPDF(paste("./gene_dendrogram_and_module_colors_min50", status, sep = "_"), height = 10, width = 15)
    plotDendroAndColors(geneTree, dynamic.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    dev.off()

    #Calculate module eigengenes
    ME.list <- moduleEigengenes(mi.expr, colors = dynamic.colors, softPower = softPower, nPC = 1)
    ME.genes <- ME.list$eigengenes
    MEDiss <- 1 - cor(ME.genes)
    METree <- flashClust(as.dist(MEDiss), method = "average")
    saveRDS.gz(METree, file = "./save/me.tree.rda")

    CairoPDF(paste("./module_eigengene_clustering_min50", status, sep = "_"), height = 10, width = 15)
    plot(METree, xlab = "", sub = "", main = "")
    dev.off()

    ##Check if any modules are too similar and merge them.  Possibly not working.
    #ME.dissimilarity.threshold <- 0.20
    #merge.all <- mergeCloseModules(mi.expr, dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3) #PC analysis may be failing because of Intel MKL Lapack routine bug.  Test with openBLAS in R compiled with gcc.
    #merged.colors <- merge.all$colors
    #merged.genes <- merge.all$newMEs

    #CairoPDF(paste("module_eigengene_clustering_min50", status, sep = "_"), height = 10, width = 15)
    #plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
    #dev.off()

    #Use merged eigengenes 
    #module.colors <- merged.colors
    #saveRDS.gz(module.colors, file = "./save/module.colors.rda")
    #color.order <- c("grey", standardColors(50))
    #modules.labels <- match(module.colors, color.order)
    #saveRDS.gz(modules.labels, file = "./save/modules.labels.rda")
    #ME.genes <- merged.genes
    #saveRDS.gz(ME.genes, file = "./save/me.genes.rda")

    #all.degrees <- intramodularConnectivity(adjacency.PEER, module.colors)
    #fdata <- featureData(lumi.import)@data
    #colnames(fdata) %<>% tolower %>% capitalize
    #gene.info.join <- data.frame(nuID = featureNames(lumi.import), select(fdata, Accession, Symbol, Definition))
    #gene.info <- mutate(gene.info.join, module.colors = module.colors, mean.count = apply(expr.data.PEER, 2, mean)) %>% data.frame(all.degrees)

    #CairoPDF(paste("eigengenes", status, sep = "_"), height = 10, width = 18)
    #par(cex = 0.7)
    #plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.adjacency = 0.3, cex.preservation = 0.3, plotPreservation = "standard")
    #dev.off()

    symbol.vector <- mi.expr.df$Probe_Id %>% str_replace("\\.", "\\-")
    module.expr.out <- data.frame(Symbol = symbol.vector, module.color = dynamic.colors)
    write.xlsx(module.expr.out, paste("./module", status, "xlsx", sep = "."))
    return(module.expr.out)
}

mi.networks <- lapply(names(genes), gen.mi, intensities.means, genes)
mi.network <- gen.mi("Patient", intensities.means, genes)
#names(mi.networks) <- names(genes.list)
#saveRDS.gz(mi.networks, file = "./save/mi.networks.rda")

module.patients <- read.xlsx("./module.Patient.xlsx")
module.symbols <- split(module.patients, module.patients$module.color)

source("../../code/GO/enrichr.R")
stringdb.submit <- function(module.color, module.df)
{
    module.genes <- module.df[[module.color]]
    get.stringdb(module.genes, module.color)
}
map(names(module.symbols), stringdb.submit, module.symbols)

enrichr.submit <- function(index, full.df, enrichr.terms, use.weights)
{
    dataset <- full.df[[index]]
    dir.create(file.path("./enrichr", index), showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), enrichr.wkbk, enrichr.data, index)
}

enrichr.wkbk <- function(subindex, full.df, index)
{
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    dir.create(file.path("./enrichr", index))
    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 
trap1 <- map(names(module.symbols), enrichr.submit, module.symbols, enrichr.terms, FALSE)

normalize.expr <- function(dataset)
{
    PIDNS <- dataset$PIDN
    dataset %<>% select(-Probe_Id, -PIDN, -Status)
    dataset.normalized <- sweep(dataset, 1, dataset$`1`)
    dataset.normalized$PIDN <- PIDNS
    return(dataset.normalized)
}

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#mi.expr.df <- slice(intensities.patient, match(genes$Patient, Probe_Id)) 
#rownames(mi.expr.df) <- mi.expr.df$Probe_Id
#mi.expr <- select(mi.expr.df, -Probe_Id) %>% t

#mi.mrnet <- minet(mi.expr, method = "mrnet", estimator = "mi.shrink", disc = "equalwidth")
#TOM.PEER <- TOMsimilarity(mi.mrnet, verbose = 5)
#colnames(TOM.PEER) <- mi.expr.df$Probe_Id
#rownames(TOM.PEER) <- mi.expr.df$Probe_Id
#TOM.pruned <- TOM.PEER
#TOM.pruned[TOM.pruned < 0.001] <- 0

#scaleFreePlot(TOM.pruned)

#mi.graph <- graph_from_adjacency_matrix(TOM.PEER, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = TRUE)
#vertex_attr(mi.graph)$name <- mi.expr.df$Probe_Id
#edge.weights <- edge_attr(mi.graph, "weight")
#pruned.subnet <- delete.edges(mi.graph, which(edge.weights < .1))
#num.edges <- map(1:vcount(pruned.subnet), incident, graph = pruned.subnet) %>% map_dbl(length) 
#pruned.subnet2 <- delete.vertices(pruned.subnet, which(num.edges == 0))
#vertex.colors <- rainbow(vcount(pruned.subnet2))
#V(pruned.subnet2)$color <- vertex.colors
##edge.df <- data.frame(edge_attr(symbols.subnet))
##edge.thickness <- edge.df$combined_score / 200
#comp.membership <- components(pruned.subnet2)$membership
#component1 <- delete.vertices(pruned.subnet2, which(comp.membership != 1))

#CairoPDF("mi.network", width = 30, height = 30)
#plot.igraph(component1, vertex.size = 2, vertex.label.dist = 0.12, vertex.label.degree = -pi/2, vertex.label.font = 2, vertex.label.color = "black", edge.color = "#0000FF99")
#dev.off()
