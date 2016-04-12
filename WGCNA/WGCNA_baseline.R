#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()

#For baseline processing
library(limma)
library(sva)
library(R.utils)
library(Biobase)

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

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(igraph)
library(TeachingDemos)

#Reading and writing tables
library(readr)
library(openxlsx)
library(parallel)

match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|"))

saveRDS.gz <- function(object,file,threads=parallel::detectCores()) 
{
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

readRDS.gz <- function(file,threads=parallel::detectCores()) 
{
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

gen.pcaplot <- function(filename, dataset, facet.bool, size.height, size.width)
{
    colnames(dataset)[2] <- "Module"
    dataset$Module %<>% str_replace("ME", "") 
    p <- ggplot(dataset, aes(x = as.numeric(x), y = value, fill = Module, color = Module)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("First Principal Component")
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(dataset$x)))
    p <- p + scale_color_manual(values = sort(unique(dataset$Module)))
    if (facet.bool == TRUE)
    {
        p <- p + facet_wrap(~ Module)
        p <- p + theme(legend.position = "none")
    } 
    CairoPDF(filename, height = size.height, width = size.width)
    print(p)
    dev.off()
}

gen.heatmap <- function(dataset, ME.genes)
{
    color <- as.character(unique(dataset$module.colors))
    dataset %<>% select(-module.colors) %>% scale
    max.dataset <- max(abs(dataset))
    print(dim(dataset))
    CairoPDF(paste("./modules/", color, sep = ""), width = 21, height = 12)
    par(mar = c(3.5,3,2,3))
    par(oma = c(4,0,2,0))
    plotMat(dataset, zlim = c(-max.dataset, max.dataset), main = paste(color, " (", nrow(dataset), ")", sep = ""))

    ME.genes.plot <- select(ME.genes, Sample.ID, matches(color))
    p <- ggplot(ME.genes.plot, aes_string(x = "Sample.ID", y = color))
    p <- p + geom_bar(stat = "identity") + xlab("Eigengene Expression")#+ ylim(c(-6, 16)) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 2))  
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    print(p)
    dev.off()
}

enrichr.submit <- function(index, full.df, enrichr.terms, use.weights)
{
    dataset <- filter(full.df, module.colors == index)
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

    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

plot.eigencor <- function(module.traits.pval, status.col, status.vector)
{
    sig.col <- paste(status.col, ".p.value", sep = "")
    cor.status.labeled <- data.frame(Color = rownames(module.traits.pval), select_(data.frame(module.traits.pval), sig.col))
    filter.cond <- paste(sig.col, "< 0.05")
    status.sig <- filter_(cor.status.labeled, filter.cond)
    me.genes.status <- select(ME.genes, one_of(as.character(status.sig$Color)))
    me.genes.status$Status <- status.vector
    split.cols <- str_split(status.col, "\\.")[[1]]
    me.status.melt <- melt(me.genes.status, id.vars = "Status") %>% filter(Status == split.cols[1] | Status == split.cols[2])
    colnames(me.status.melt)[2] <- "Module"

    sum.fun <- function(data.vector){ data.frame(ymin = min(data.vector), ymax = max(data.vector), y = mean(data.vector)) }
    me.status.melt$Module %<>% as.character
    p <- ggplot(me.status.melt, aes(x = factor(Status), y = value, col = Module)) + geom_point(position = "jitter")
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
    p <- p + theme(axis.title.x = element_blank()) + stat_summary(aes(group = 1), fun.y = mean, geom = "line", col = "black", position = position_dodge(width = 0.9))
    p <- p + scale_color_manual(values = sort(unique(me.status.melt$Module)))
    p <- p + facet_wrap(~ Module, scales = "free_y")
    p <- p + theme(legend.position = "none")

    filename <- paste(split.cols[1], "_", split.cols[2], "_eigengenes_05", sep = "")
    CairoPDF(filename, height = 13, width = 20)
    print(p)
    dev.off()
}

lumi.import <- readRDS.gz(file = "../baseline_lumi/save/export.lumi.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

fdata <- fData(lumi.import)
lumi.import <- lumi.import[!is.na(fdata$SYMBOL),]
fdata <- fData(lumi.import)
expr.data.PEER <- exprs(lumi.import)
expr.collapse <- collapseRows(expr.data.PEER, factor(fdata$SYMBOL), rownames(expr.data.PEER))$datETcollapsed %>% t
saveRDS.gz(expr.collapse, "./save/expr.collapse.rda")

sft.PEER <- pickSoftThreshold(expr.collapse, powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.PEER.df <- sft.PEER$fitIndices
saveRDS.gz(sft.PEER, file = "./save/sft.peer.rda")

#Plot scale indendence and mean connectivity as functions of power
sft.PEER.df$multiplied <- sft.PEER.df$SFT.R.sq * -sign(sft.PEER.df$slope)
p <- ggplot(sft.PEER.df, aes(x = Power,  y = multiplied, label = Power)) + geom_point() + geom_text(vjust = -0.6, size = 4, col = "red")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(aes(yintercept = 0.9))
p <- p + xlab("Soft Threshold") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale Independence")
CairoPDF(file = "./scaleindependence", width = 6, height = 6)
print(p)
dev.off()

p <- ggplot(sft.PEER.df, aes(x = Power,  y = mean.k.)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Soft Threshold") + ylab("Mean Connectivity") + ggtitle("Mean Connectivity")
CairoPDF(file = "./meanconnectivity", width = 6, height = 6)
print(p)
dev.off()

softPower <- 10
adjacency.PEER <- adjacency(expr.collapse, power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
saveRDS.gz(adjacency.PEER, file = "./save/adjacency.PEER.rda")

TOM.PEER <- TOMsimilarity(adjacency.PEER, verbose = 5)
dissimilarity.TOM <- 1 - TOM.PEER
saveRDS.gz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 20

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, cutHeight = 0.995, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.collapse, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50_tree", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(expr.collapse, dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) #PC analysis may be failing because of Intel MKL Lapack routine bug.  Test with openBLAS in R compiled with gcc.
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
saveRDS.gz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
saveRDS.gz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
saveRDS.gz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 10, width = 18)
par(cex = 0.7)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.adjacency = 0.3, cex.preservation = 0.3, plotPreservation = "standard")
dev.off()

modules.out <- data.frame(Symbol = colnames(expr.collapse), Module = module.colors)
write.xlsx(modules.out, "modules_out.xlsx")

eigengene.model <- function(ME.vector, trait.vector)
{
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.aov <- aov(ME ~ Trait, trait.df) %>% TukeyHSD
    trait.tukey <- data.frame(trait.aov$Trait)
    colnames(trait.tukey) %<>% str_replace(" ", ".")
    return(select(trait.tukey, diff, p.adj))
}

eigengene.anova <- function(ME.vector, trait.vector)
{
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- lm(ME ~ Trait, trait.df) %>% anova
    return(trait.anova$'Pr(>F)'[1])
}

targets.final.known <- pData(lumi.import)
targets.final.known$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))
anova.status <- map_dbl(ME.genes, eigengene.anova, targets.final.known$Status) %>% p.adjust("fdr")

cor.status <- map(ME.genes, eigengene.model, targets.final.known$Status)
status.diff <- map(cor.status, select, diff) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.diff) <- names(cor.status)
status.pval <- map(cor.status, select, p.adj) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.pval) <- names(cor.status)

pval.adjust <- map(status.pval, p.adjust, method = "fdr") %>% reduce(cbind) %>% data.frame
rownames(pval.adjust) <- rownames(status.pval)
colnames(pval.adjust) <- paste(colnames(status.pval), ".pval")
status.pval$Module <- rownames(status.pval)

text.matrix.traits <- paste(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(pval.adjust), 1), ')', sep = '')
dim(text.matrix.traits) = dim(status.diff)

heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
gen.text.heatmap(as.matrix(status.diff), text.matrix.traits, colnames(status.diff), colnames(ME.genes), "", "module-trait relationships", heatmap.range)

all.degrees <- intramodularConnectivity(adjacency.PEER, module.colors)
colnames(fdata) %<>% tolower %>% capitalize
gene.info <- data.frame(Symbol = rownames(all.degrees), module.color = module.colors, all.degrees)

write_csv(data.frame(table(module.colors)), path = "./final_eigengenes.csv") 
gene.info$kscaled <- by(gene.info, gene.info$module.color, select, kWithin) %>% llply(function(x) { x / max (x) }) %>% reduce(c)
saveRDS.gz(gene.info, file = "./save/gene.info.rda")

gene.module.membership <- as.data.frame(bicor(expr.collapse, ME.genes, maxPOutliers = 0.05))
module.membership.pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.collapse)))
names(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")

module.membership <- cbind(select(gene.info, Symbol:module.color), gene.module.membership, module.membership.pvalue)
write_csv(module.membership, "module_membership.csv")

colnames(gene.module.membership) %<>% str_replace("MM.", "")
colnames(module.membership.pvalue) %<>% str_replace("MM.pvalue.", "")
gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

color.values <- unique(module.colors)
color.key <- paste(head(color.values, n = 1), tail(color.values, n = 1), sep = ":")
gene.module.membership.long <- gather_(data = gene.module.membership, "module.comparison", "correlation", color.values)
module.membership.pvalue.long <- gather_(data = module.membership.pvalue, "module.comparison", "p.value", color.values)
membership.join <- join(gene.module.membership.long, module.membership.pvalue.long)
eigengene.connectivity <- join(membership.join, gene.info) %>% select(Symbol, module.color:kscaled, module.comparison:p.value)
write_csv(eigengene.connectivity, "eigengene_connectivity.csv")
all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y")
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$x <- as.factor(1:nrow(smooth.df))
smooth.plot <- melt(smooth.df, id.vars = "x")

gen.pcaplot("all_principal_components", smooth.plot, FALSE, 10, 15)
gen.pcaplot("facet_principal_components", smooth.plot, TRUE, 13, 25)

sample.ids <- factor(rownames(expr.collapse), levels = rownames(expr.collapse))
colnames(ME.genes) %<>% str_replace("ME", "")
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.collapse), module.colors)
cluster <- makeForkCluster(8)
#split(expr.data.plot, expr.data.plot$module.colors) %>% parLapply(cl = cluster, gen.heatmap, ME.genes.plot)
#by(expr.data.plot, expr.data.plot$module.colors, gen.heatmap, ME.genes.plot)
green.plot <- filter(expr.data.plot, module.colors == "green")
green.table <- filter(module.membership, module.color == "green") #%>% select(Symbol, MM.green)
gen.heatmap(green.plot, ME.genes.plot)

#modules.out <- select(gene.info, Symbol, module.color)
#targets.final.known$Sample.Name %<>% str_replace(" ", "")

source("../../code/GO/enrichr.R")

enrichr.terms <- c("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "HumanCyc_2016", "NCI-Nature_2016", "Panther_2016") 
#enrichr.submit("grey", modules.out, enrichr.terms, FALSE)
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names, enrichr.submit, modules.out, enrichr.terms, FALSE)
#blue.test <- filter(modules.out, module.color == "blue")
#test <- get.enrichrdata("GO_Biological_Process", blue.test, FALSE)
#test <- enrichr.submit("blue", modules.out, enrichr.terms, FALSE)

split(modules.out, modules.out$module.color) %>% map(submit.stringdb)
green.submit <- filter(modules.out, Module == "green")
green.ppi <- get.stringdb(green.submit, "green", "./stringdb")

green.ppi.edges <- attr(E(green.ppi), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_")
green.ppi.df <- data.frame(Edge = green.ppi.edges, Weight = edge_attr(green.ppi, "combined_score"))
green.ppi.clustered <- names(V(green.ppi))[clusters(green.ppi)$membership == 1]

green.members <- filter(modules.out, Module == "green")$Symbol %>% map_chr(match.exact) %>% paste(collapse = "|")
adjacency.green <- adjacency.PEER[grepl(green.members, rownames(adjacency.PEER)), grepl(green.members, colnames(adjacency.PEER))]
green.igraph <- graph_from_adjacency_matrix(adjacency.green, mode = "undirected", weighted = TRUE, diag = FALSE)

num.edges.green <- map(1:vcount(green.ppi), incident, graph = green.ppi) %>% map_dbl(length) #Computer number of edges for each vertex
names(num.edges.green) <- names(V(green.ppi)) #Label number edges with 
vertices.green <- names(V(green.ppi)) %>% map_chr(match.exact) %>% paste(collapse = "|")
missing.edges.green <- !grepl(vertices.green, names(V(green.igraph)))
pruned.green <- delete.vertices(green.igraph, which(missing.edges.green))

num.edges.green <- num.edges.green[match(names(V(pruned.green)), names(num.edges.green))]
pruned.green <- delete.vertices(pruned.green, which(num.edges.green == 0))
pruned.green.ppi <- delete.vertices(green.ppi, which(num.edges.green == 0))

green.communities.optimal <- cluster_optimal(pruned.green.ppi, weights = edge_attr(pruned.green.ppi, "combined_score")/1000)

cluster.todf <- function(cluster.element)
{
    returnval <- data.frame(cluster.element[1], "Group" = names(cluster.element))
    names(returnval) <- c("Members", "Group")
    return(list(returnval))
}
green.communities.df <- communities(green.communities.optimal) %>% lmap(cluster.todf) %>% reduce(rbind)

pruned.green.edges <- attr(E(pruned.green), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_") 
pruned.green.df <- data.frame(Edge = pruned.green.edges, Weight = edge_attr(pruned.green, "weight")) 

green.filter <- map_chr(green.ppi.df$Edge, match.exact) %>% paste(collapse = "|") %>% grepl(pruned.green.df$Edge)
pruned.green.df[green.filter,]$Weight <- green.ppi.df$Weight / 200
pruned.green.df[!green.filter,]$Weight <- 0
pruned.green.df$Color <- "#dddddd99"
pruned.green.df[green.filter,]$Color <- "#0000FF99"

green.colors <- rainbow(length(unique(green.communities.df$Group)))
green.communities.sort <- green.communities.df[match(names(V(pruned.green)), green.communities.df$Members),]
green.colors.assign <- green.colors[green.communities.sort$Group]

V(pruned.green)$color <- green.colors.assign
edge.df <- data.frame(edge_attr(pruned.green))
#edge.thickness <- edge.df$combined_score / 20
green.nchars <- vertex_attr(pruned.green, "name") %>% map_int(nchar)
dist.green <- .14 / (max(green.nchars) - min(green.nchars))
green.dists <- ((max(green.nchars) - green.nchars)^1.3 * dist.green) + 0.26

CairoPDF("green_igraph.pdf", width = 10, height = 10)
par(mar=c(0,0,0,0) + 0.1)
plot.igraph(pruned.green, vertex.size = 6, vertex.label.dist = green.dists, vertex.label.degree = pi/2, vertex.label.font = 2, vertex.label.color = "black", margin = 0, edge.color = pruned.green.df$Color, edge.width = pruned.green.df$Weight)
dev.off()

green.ppi.colors <- rainbow(length(unique(green.communities.df$Group)))
green.ppi.communities.sort <- green.communities.df[match(names(V(pruned.green.ppi)), green.communities.df$Members),]
green.ppi.colors.assign <- green.ppi.colors[green.ppi.communities.sort$Group]
V(pruned.green.ppi)$color <- green.ppi.colors.assign
CairoPDF("green_igraph_original.pdf", width = 10, height = 10)
par(mar=c(0,0,0,0) + 0.1)
plot.igraph(pruned.green.ppi, vertex.size = 6, vertex.label.degree = pi/2, vertex.label.font = 2, vertex.label.color = "black", margin = 0, edge.color = "#dddddd99", edge.width = green.ppi.df$Weight / 200)
dev.off()

submit.stringdb <- function(module.subset)
{
    catch <- get.stringdb(module.subset, unique(module.subset$module.color), "./stringdb")
}

rownames(ME.genes) <- rownames(expr.collapse)
PEER.factors <- read_csv("../baseline_lumi/factor_8.csv") %>% select(-(X1:X6))
rownames(PEER.factors) <- rownames(expr.collapse)
colnames(PEER.factors) <- paste("X", 1:ncol(PEER.factors), sep = "")

source("../common_functions.R")

#targets.age <- select(targets.final.known, Sample.Name, Draw.Age) #%>% filter(!is.na(Draw.Age))
#targets.age$Draw.Age %<>% as.numeric
#cor.age <- gen.cor(ME.genes, targets.age)

#targets.sex <- filter(targets.final.known, Sex != "UNKNOWN")
#targets.sex.m <- model.matrix( ~ 0 + factor(targets.sex$Sex) )[,-1] %>% data.frame 
#colnames(targets.sex.m) <- c("Sex")
#targets.sex.m %<>% mutate(Sample.Name = targets.final.known$Sample.Name)
#cor.sex <- gen.cor(ME.genes, targets.sex.m)

#targets.rin <- select(targets.final.known, Sample.Name, RIN) 
#cor.rin <- gen.cor(ME.genes, targets.rin)

#targets.onset <- select(targets.final.known, Sample.Name, Onset) %>% filter(!is.na(Onset))
#targets.onset$Onset %<>% as.numeric
#cor.onset.PEER <- gen.cor(PEER.factors, targets.onset)

#PEER.factors.plot <- select(read_csv("../baseline_lumi/factor_8.csv"), -(X1:X6))
#names(PEER.factors.plot) <- paste("X", 1:ncol(PEER.factors.plot), sep = "")
#PEER.factors.df <- mutate(PEER.factors.plot, Sample.Name = targets.final.known$Sample.Name) 
#cor.PEER <- gen.cor(ME.genes, PEER.factors.df) %>% data.frame

#module.traits.all <- cbind(cor.status, cor.age, cor.sex, cor.rin) %>% data.frame
#module.traits.pval <- select(module.traits.all, contains("p.value")) %>% as.matrix %>% apply(2, p.adjust, "fdr")
#module.traits.cor <- select(module.traits.all, -contains("p.value")) %>% as.matrix

#module.PEER.all <- cbind(cor.gaa.PEER, cor.onset.PEER) %>% data.frame
#module.PEER.pval <- select(module.PEER.all, contains("p.value")) %>% as.matrix %>% apply(2, p.adjust, "fdr")
#module.PEER.cor <- select(module.PEER.all, -contains("p.value")) %>% as.matrix

#cor.PEER.cor <- select(cor.PEER, -contains("p.value")) %>% as.matrix
#cor.PEER.pval <- select(cor.PEER, contains("p.value")) %>% as.matrix

#module.trait.out <- data.frame(Module = rownames(module.traits.cor), module.traits.cor, module.traits.pval)
#module.PEER.out <- data.frame(Factor = rownames(module.PEER.cor), module.PEER.cor, module.PEER.pval)
#cor.PEER.out <- data.frame(Module = rownames(cor.PEER), cor.PEER.cor, cor.PEER.pval)

#write_csv(module.trait.out, "module_trait_cor.csv")
#write_csv(module.PEER.out, "module_PEER_cor.csv")
#write_csv(cor.PEER.out, "cor_PEER_cor.csv")

#text.matrix.PEER <- paste(signif(module.PEER.cor, 2), '\n(', signif(module.PEER.pval, 1), ')', sep = '')
#text.matrix.PEERcor <- paste(signif(cor.PEER.cor, 2), '\n(', signif(cor.PEER.pval, 1), ')', sep = '')
#dim(text.matrix.PEER) = dim(module.PEER.cor)
#dim(text.matrix.PEERcor) = dim(cor.PEER.cor)

#gen.text.heatmap(module.PEER.cor, text.matrix.PEER, colnames(module.PEER.cor), rownames(module.PEER.cor), "", "PEER factor-trait relationships")
#gen.text.heatmap(cor.PEER.cor, text.matrix.PEERcor, colnames(cor.PEER.cor), colnames(ME.genes), "", "module-PEER factor relationships")

plot.eigencor(module.traits.pval, "Control.Carrier", lumi.import$Status)
plot.eigencor(module.traits.pval, "Control.Patient", lumi.import$Status)
plot.eigencor(module.traits.pval, "Carrier.Patient", lumi.import$Status)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#Final plots
lightgreenyellow.biol <- read.xlsx("./enrichr/lightgreen/lightgreen_GO_Biological_Process.xlsx") %>% slice(c(1,2))
lightgreen.biol$Database <- "GO Biological Process"
lightgreen.kegg <- read.xlsx("./enrichr/lightgreen/lightgreen_KEGG_2015.xlsx") %>% slice(c(1,3))
lightgreen.kegg$Database <- "KEGG"
lightgreen.reactome <- read.xlsx("./enrichr/lightgreen/lightgreen_Reactome_2015.xlsx") %>% slice(1)
lightgreen.reactome$Database <- "Reactome"
lightgreen.enrichr <- rbind(lightgreen.biol, lightgreen.kegg, lightgreen.reactome)
lightgreen.enrichr$Gene.Count <- map(lightgreen.enrichr$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
lightgreen.enrichr$Log.pvalue <- -(log10(lightgreen.enrichr$P.value))

lightgreen.enrichr$GO.Term %<>% str_replace_all("\\ \\(.*$", "") %>% tolower
lightgreen.enrichr$Format.Name <- paste(lightgreen.enrichr$Database, ": ", lightgreen.enrichr$GO.Term, " (", lightgreen.enrichr$Gene.Count, ")", sep = "")
lightgreen.enrichr.plot <- select(lightgreen.enrichr, Format.Name, Log.pvalue) %>% melt(id.vars = "Format.Name") 

p <- ggplot(lightgreen.enrichr.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = lightgreen.enrichr$Format.Name, hjust = "left", aes(y = 0.1))
p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "FALSE",  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
CairoPDF("lightgreen.enrichr", height = 5, width = 8)
print(p)
dev.off()

greenyellow.molec <- read.xlsx("./enrichr/greenyellow/greenyellow_GO_Molecular_Function.xlsx") %>% slice(c(1,2))
greenyellow.molec$Database <- "GO Molecular Function"
greenyellow.kegg <- read.xlsx("./enrichr/greenyellow/greenyellow_KEGG_2015.xlsx") %>% slice(c(1,3,6))
greenyellow.kegg$Database <- "KEGG"
greenyellow.reactome <- read.xlsx("./enrichr/greenyellow/greenyellow_Reactome_2015.xlsx") %>% slice(c(1,36,38,39))
greenyellow.reactome$Database <- "Reactome"
greenyellow.enrichr <- rbind(greenyellow.molec, greenyellow.kegg, greenyellow.reactome)
greenyellow.enrichr$Gene.Count <- map(greenyellow.enrichr$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
greenyellow.enrichr$Log.pvalue <- -(log10(greenyellow.enrichr$P.value))

greenyellow.enrichr$GO.Term %<>% str_replace_all("\\ \\(.*$", "") %>% tolower
greenyellow.enrichr$Format.Name <- paste(greenyellow.enrichr$Database, ": ", greenyellow.enrichr$GO.Term, " (", greenyellow.enrichr$Gene.Count, ")", sep = "")
greenyellow.enrichr.plot <- select(greenyellow.enrichr, Format.Name, Log.pvalue) %>% melt(id.vars = "Format.Name") 

p <- ggplot(greenyellow.enrichr.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = greenyellow.enrichr$Format.Name, hjust = "left", aes(y = 0.1))
p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "FALSE",  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
CairoPDF("greenyellow.enrichr", height = 5, width = 9)
print(p)
dev.off()

lightgreen.eigen <- data.frame(Eigengene = ME.genes$lightgreen, Status = lumi.import$Status) %>% filter(Status != "Carrier")
p <- ggplot(lightgreen.eigen, aes(x = Status, y = Eigengene)) + geom_point(position = "jitter", col = "lightgreen")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
p <- p + theme(axis.title.x = element_blank()) + stat_summary(aes(group = 1), fun.y = mean, geom = "line", col = "black", position = position_dodge(width = 0.9))
q <- p + theme(legend.position = "none")

CairoPDF("lightgreen.eigengene", height = 5, width = 7)
print(p)
dev.off()

greenyellow.eigen <- data.frame(Eigengene = ME.genes$greenyellow, Status = lumi.import$Status) %>% filter(Status != "Carrier")
p <- ggplot(greenyellow.eigen, aes(x = Status, y = Eigengene)) + geom_point(position = "jitter", col = "greenyellow")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
p <- p + theme(axis.title.x = element_blank()) + stat_summary(aes(group = 1), fun.y = mean, geom = "line", col = "black", position = position_dodge(width = 0.9))
p <- p + theme(legend.position = "none")

CairoPDF("greenyellow.eigengene", height = 5, width = 7)
print(p)
dev.off()

#get.stringdb(data.frame(Symbol = str_split(lightgreen.biol$Genes[2], ",")[[1]]), "lightgreen_celldeath")
match.exact <- paste %<<<% "^" %<<% c("$", sep = "")

pcd.submit <- data.frame(Symbol = str_split(lightgreen.biol$Genes[2], ",")[[1]])
pcd.ppi <- get.stringdb(pcd.submit, "pcd.ppi")
pcd.ppi.edges <- attr(E(pcd.ppi), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_")
pcd.ppi.df <- data.frame(Edge = pcd.ppi.edges, Weight = edge_attr(pcd.ppi, "combined_score"))
pcd.ppi.clustered <- names(V(pcd.ppi))[clusters(pcd.ppi)$membership == 1]

pcd.members <- str_split(lightgreen.biol$Genes[2], ",")[[1]] %>% map_chr(match.exact) %>% paste(collapse = "|")
#adjacency.pcd <- adjacency.PEER[grepl(pcd.members, rownames(adjacency.PEER)), grepl(pcd.members, colnames(adjacency.PEER))]
pcd.igraph <- graph_from_adjacency_matrix(adjacency.pcd, mode = "undirected", weighted = TRUE, diag = FALSE)

num.edges.pcd <- map(1:vcount(pcd.ppi), incident, graph = pcd.ppi) %>% map_dbl(length) 
names(num.edges.pcd) <- names(V(pcd.ppi))
num.edges.pcd <- num.edges.pcd[match(names(V(pcd.igraph)), names(num.edges.pcd))]
pruned.pcd <- delete.vertices(pcd.igraph, which(num.edges.pcd == 0))

pruned.pcd.edges <- attr(E(pruned.pcd), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_") 
pruned.pcd.df <- data.frame(Edge = pruned.pcd.edges, Weight = edge_attr(pruned.pcd, "weight")) 

pcd.filter <- map_chr(pcd.ppi.df$Edge, match.exact) %>% paste(collapse = "|") %>% grepl(pruned.pcd.df$Edge)
pruned.pcd.df[pcd.filter,]$Weight <- pcd.ppi.df$Weight / 200
pruned.pcd.df[!pcd.filter,]$Weight <- 0
pruned.pcd.df$Color <- "#dddddd99"
pruned.pcd.df[pcd.filter,]$Color <- "#0000FF99"

pcd.colors <- rep("magenta", length(V(pruned.pcd)))
pcd.colors[match(pcd.ppi.clustered, names(V(pruned.pcd)))] <- "turquoise"
V(pruned.pcd)$color <- pcd.colors

pcd.nchars <- vertex_attr(pruned.pcd, "name") %>% map_int(nchar)
dist.pcd <- .14 / (max(pcd.nchars) - min(pcd.nchars))
pcd.dists <- ((max(pcd.nchars) - pcd.nchars)^1.3 * dist.pcd) + 0.76

CairoPDF("pcd_igraph.pdf", width = 6, height = 6)
par(mar=c(0,0,0,0) + 0.1)
plot.igraph(pruned.pcd, vertex.size = 20, vertex.label.dist = pcd.dists, vertex.label.degree = pi/2, vertex.label.font = 2, vertex.label.color = "black", edge.color = pruned.pcd.df$Color, edge.width = pruned.pcd.df$Weight)
dev.off()

oxo.submit <- data.frame(Symbol = str_split(greenyellow.molec$Genes[1], ",")[[1]])
oxo.ppi <- get.stringdb(oxo.submit, "oxo.ppi")
oxo.ppi.edges <- attr(E(oxo.ppi), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_")
oxo.ppi.df <- data.frame(Edge = oxo.ppi.edges, Weight = edge_attr(oxo.ppi, "combined_score"))
oxo.ppi.clustered <- names(V(oxo.ppi))[clusters(oxo.ppi)$membership == 1]

oxo.members <- str_split(greenyellow.molec$Genes[1], ",")[[1]] %>% map_chr(match.exact) %>% paste(collapse = "|")
#adjacency.oxo <- adjacency.PEER[grepl(oxo.members, rownames(adjacency.PEER)), grepl(oxo.members, colnames(adjacency.PEER))]
oxo.igraph <- graph_from_adjacency_matrix(adjacency.oxo, mode = "undirected", weighted = TRUE, diag = FALSE)

num.edges.oxo <- map(1:vcount(oxo.ppi), incident, graph = oxo.ppi) %>% map_dbl(length) 
names(num.edges.oxo) <- names(V(oxo.ppi))
num.edges.oxo <- num.edges.oxo[match(names(V(oxo.igraph)), names(num.edges.oxo))]
pruned.oxo <- delete.vertices(oxo.igraph, which(num.edges.oxo == 0))

pruned.oxo.edges <- attr(E(pruned.oxo), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_") 
pruned.oxo.df <- data.frame(Edge = pruned.oxo.edges, Weight = edge_attr(pruned.oxo, "weight")) 

oxo.filter <- map_chr(oxo.ppi.df$Edge, match.exact) %>% paste(collapse = "|") %>% grepl(pruned.oxo.df$Edge)
pruned.oxo.df[oxo.filter,]$Weight <- oxo.ppi.df$Weight / 200
pruned.oxo.df[!oxo.filter,]$Weight <- 0
pruned.oxo.df$Color <- "#dddddd99"
pruned.oxo.df[oxo.filter,]$Color <- "#0000FF99"

oxo.colors <- rep("magenta", length(V(pruned.oxo)))
oxo.colors[match(oxo.ppi.clustered, names(V(pruned.oxo)))] <- "turquoise"

V(pruned.oxo)$color <- oxo.colors
edge.df <- data.frame(edge_attr(pruned.oxo))
#edge.thickness <- edge.df$combined_score / 20
oxo.nchars <- vertex_attr(pruned.oxo, "name") %>% map_int(nchar)
dist.oxo <- .14 / (max(oxo.nchars) - min(oxo.nchars))
oxo.dists <- ((max(oxo.nchars) - oxo.nchars)^1.3 * dist.oxo) + 0.76

CairoPDF("oxo_igraph.pdf", width = 6, height = 6)
par(mar=c(0,0,0,0) + 0.1)
plot.igraph(pruned.oxo, vertex.size = 20, vertex.label.dist = oxo.dists, vertex.label.degree = pi/2, vertex.label.font = 2, vertex.label.color = "black", edge.color = pruned.oxo.df$Color, edge.width = pruned.oxo.df$Weight)
dev.off()


ME.lightgreen.plot <- data.frame(Sample.ID = rownames(ME.genes), select(ME.genes, lightgreen), Status = lumi.import$Status) %>% filter(Status != "Carrier") %>% arrange(Status, Sample.ID)
ME.lightgreen.plot$Sample.ID <- factor(ME.lightgreen.plot$Sample.ID, levels = ME.lightgreen.plot$Sample.ID)
p <- ggplot(ME.lightgreen.plot, aes(x = Sample.ID, y = lightgreen)) + geom_bar(stat = "identity", aes(fill = Status))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
CairoPDF("lightgreen_pco", width = 8, height = 5)
print(p)
dev.off()

ME.greenyellow.plot <- data.frame(Sample.ID = rownames(ME.genes), select(ME.genes, greenyellow), Status = lumi.import$Status) %>% filter(Status != "Carrier") %>% arrange(Status, Sample.ID)
ME.greenyellow.plot$Sample.ID <- factor(ME.greenyellow.plot$Sample.ID, levels = ME.greenyellow.plot$Sample.ID)
p <- ggplot(ME.greenyellow.plot, aes(x = Sample.ID, y = greenyellow)) + geom_bar(stat = "identity", aes(fill = Status))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
CairoPDF("greenyellow_pco", width = 8, height = 5)
print(p)
dev.off()
