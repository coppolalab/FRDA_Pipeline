library(WGCNA)
library(flashClust)
enableWGCNAThreads()
library(BayesFactor)
library(limma)
library(R.utils)
library(lumi)
library(biomaRt)

library(igraph)
library(UpSetR)
library(openxlsx)
library(Cairo)

#String operations
library(tools)
library(broom)
library(stringr)
library(magrittr)
library(tidyverse)

EigengeneANOVA <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- lm(ME ~ Trait, trait.df) %>% anova %>% tidy
    trait.anova$p.value[1]
}

EigengeneBayes <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- anovaBF(ME ~ Trait, data = trait.df) 
    trait.anova
}

EigengeneAOV <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.aov <- aov(ME ~ Trait, trait.df) %>% TukeyHSD 
    aov.contrasts <- rownames(trait.aov[[1]])
    trait.tidy <- tidy(trait.aov) %>% select(estimate, adj.p.value)
    rownames(trait.tidy) <- aov.contrasts
    trait.tidy
}

EigengeneBoxplot <- function(eigengene.df, status.vector, color) {
    color.column <- str_c("ME", color)
    gene.df <- data.frame(Status = status.vector, Expression = as.vector(eigengene.df[[color.column]]))
    gene.df$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(gene.df, aes(x = Status, y = Expression, fill = Status)) + geom_violin(scale = "width", trim = FALSE) 
    p <- p + geom_boxplot(width = 0.25, outlier.shape = NA) + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(plot.background = element_blank(), axis.title.x = element_blank()) + ggtitle(str_c(capitalize(color), " Module"))
    p <- p + theme(panel.border = element_rect(color = "black", size = 1), plot.margin = unit(c(1,1,1,1), "lines"))
    p <- p + ylab("Eigengene Value")
    CairoPDF(color, width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

EstimateViolinPlot <- function(estimate.df, color, ylimits) {
    estimate.extract <- as.matrix(estimate.df) %>% data.frame 
    estimate.groups <- select(estimate.extract, matches("^Trait")) 
    estimate.add <- sweep(estimate.groups, 1, estimate.extract$mu, "+")
    estimate.plot <- gather(estimate.add, Group, Estimate) 
    estimate.plot$Group %<>% str_replace("Trait\\.", "") %>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(estimate.plot, aes(x = Group, y = Estimate, fill = Group)) + geom_violin(scale = "width", trim = FALSE, draw_quantiles = c(0.05, 0.5, 0.95)) + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(color = "black", size = 1))
    p <- p + theme(plot.background = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) + ylim(ylimits) 
    p <- p + ylab("Eigengene Estimate") + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(str_c(capitalize(color), "Module", sep = " ")) 
    CairoPDF(str_c(color, "estimate", sep = "."), width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

ModuleWorkbook <- function(module.table, filename) {
    pval.key <- colnames(module.table) %>% str_detect("pvalue") 
    pval.cols <- which(pval.key)
    MM.key <- colnames(module.table) %>% str_detect("MM.") 
    MM.cols <- which(MM.key & !pval.key)

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = module.table)
    #sig.pvalues <- createStyle(fontColour = "red")
    #conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(module.table), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = MM.cols, rows = 1:nrow(module.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:ncol(module.table), widths = "auto")
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 11)
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

PCAPlot <- function(filename, dataset, facet.bool, size.height, size.width) {
    dataset$Module %<>% str_replace("ME", "") 
    p <- ggplot(dataset, aes(x = Gene, y = PCA1, fill = Module, color = Module)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("First Principal Component")
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
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

MEHeatmap <- function(dataset, ME.genes) {
    color <- as.character(unique(dataset$Module)) %>% str_replace("ME", "")
    colnames(ME.genes) %<>% str_replace("ME", "")
    dataset %<>% select(-Module) %>% scale
    max.dataset <- max(abs(dataset))

    CairoPDF(paste("./modules/", color, sep = ""), width = 21, height = 12)
    par(mar = c(3.5,3,2,3))
    par(oma = c(4,0,2,0))
    plotMat(dataset, zlim = c(-max.dataset, max.dataset))

    ME.genes.plot <- select(ME.genes, Sample.ID, matches(color))
    p <- ggplot(ME.genes.plot, aes_string(x = "Sample.ID", y = color))
    p <- p + geom_bar(stat = "identity") + xlab("Eigengene Expression") 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.x = element_blank())  
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    print(p)
    dev.off()
}

EnrichrWorkbook <- function(database, full.df, colname) {
    dataset <- full.df[[database]]

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 1, cols = 1, widths = 45)
    setColWidths(wb, 1, cols = c(2:ncol(dataset)), widths = "auto")
    
    dir.create(file.path("./enrichr", colname), recursive = TRUE)
    filename = str_c(file.path("./enrichr", colname, database), ".xlsx")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

EnrichrPlot <- function(enrichr.df, filename, plot.title, plot.height = 5, plot.width = 8, color = "default") {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    #enrichr.df$Log.Bayes.Factor <- log10(enrichr.df$Bayes.Factor)
    enrichr.df$Term %<>% str_replace_all("\\ \\(GO.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% str_replace_all(",.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.Bayes.Factor)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Bayes.Factor)  

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Log.Bayes.Factor)) 
    if (color == "default") {
        p <- p + geom_bar(stat = "identity") 
    } else {
        p <- p + geom_bar(stat = "identity", fill = color) 
    }
    p <- p + coord_flip() + theme_bw() + theme(legend.position = "none") 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' Bayes Factor'))) + theme(plot.background = element_blank())
    p <- p + theme(panel.border = element_rect(color = "black", size = 1))
    p <- p + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(plot.title)
    CairoPDF(str_c(filename, ".enrichr"), height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

GetKscaled <- function(gene.list, module.membership) {
    filter(module.membership, is.element(Symbol, gene.list)) %>% select(Symbol, kscaled) %>% arrange(desc(kscaled))
}

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

source("../../code/common_functions.R")
lumi.import <- ReadRDSgz(file = "../baseline_lumi/save/rmout.collapse.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.collapse <- exprs(lumi.import) %>% t
sft <- pickSoftThreshold(expr.collapse, powerVector = powers, verbose = 5, networkType = "signed", corFnc = bicor, corOptions = list(maxPOutliers = 0.05))
sft.df <- sft$fitIndices
saveRDS.gz(sft, file = "./save/sft.rda")

#Plot scale indendence and mean connectivity as functions of power
sft.df$multiplied <- sft.df$SFT.R.sq * -sign(sft.df$slope)
p <- ggplot(sft.df, aes(x = Power,  y = multiplied, label = Power)) + geom_point() + geom_text(vjust = -0.6, size = 4, col = "red")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(aes(yintercept = 0.9))
p <- p + xlab("Soft Threshold") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale Independence")
CairoPDF(file = "./scaleindependence", width = 6, height = 6)
print(p)
dev.off()

p <- ggplot(sft.df, aes(x = Power,  y = mean.k.)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Soft Threshold") + ylab("Mean Connectivity") + ggtitle("Mean Connectivity")
CairoPDF(file = "./meanconnectivity", width = 6, height = 6)
print(p)
dev.off()

softPower <- 6
adjacency.expr <- adjacency(expr.collapse, power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
gc()

TOM <- TOMsimilarity(adjacency.expr, verbose = 5)
dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
gc()
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

#Identify modules using dynamic tree cutting with hybrid clustering
min.module.size <- 20
dynamic.modules <- cutreeDynamic(dendro = geneTree, cutHeight = 0.995, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.collapse, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_tree", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(expr.collapse, dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3, corFnc = bicor, corOptions = list(maxPOutliers = 0.05)) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering", height = 8, width = 10, bg = "transparent")
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("", ""), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "All subjects", axes = FALSE, ylab = "", cex.main = 4.0)
dev.off()

CairoPDF("module_eigengene_clustering_new", height = 8, width = 10, bg = "transparent")
plot(geneTree)
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
saveRDS.gz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
saveRDS.gz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
saveRDS.gz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 8)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

#Create heatmap of eigengene significance
color.values <- unique(module.colors)
lumi.pdata <- pData(lumi.import)
lumi.pdata$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))

model.cov <- model.matrix(~ Sex + Age + RIN, data = lumi.pdata) #Create model matrix of covariates to be removed
ME.genes.rmcov <- removeBatchEffect(t(ME.genes), covariates = model.cov[,-1]) %>% 
    t %>% data.frame #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences anova.status <- map_dbl(ME.genes.rmcov, EigengeneANOVA, lumi.pdata$Status) %>% p.adjust("fdr") %>% signif(3)
bayes.status <- map(ME.genes.rmcov, EigengeneBayes, lumi.pdata$Status) 
bf.status <- map(bayes.status, extractBF) %>% map_dbl(extract2, "bf") %>% map_dbl(log10)
posterior.status <- map(bayes.status, posterior, iterations = 100000) 

aov.status <- map(ME.genes, EigengeneAOV, lumi.pdata$Status)
status.diff <- map(aov.status, select, estimate) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.diff) <- names(aov.status)
status.pval <- map(aov.status, select, adj.p.value) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.pval) <- names(aov.status)

pval.adjust <- map(status.pval, p.adjust, method = "fdr", n = length(color.values)) %>% reduce(cbind) %>% data.frame
rownames(pval.adjust) <- rownames(status.pval)
colnames(pval.adjust) <- paste(colnames(status.pval), ".pval")

text.matrix.traits <- str_c(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(pval.adjust), 1), ')')
dim(text.matrix.traits) = dim(status.diff)

heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
width.dynamic <- 3 + (1 * ncol(text.matrix.traits))

CairoPDF("module_trait_relationships", width = width.dynamic, height = 10)
par(mar = c(8, 8, 3, 3))
labeledHeatmap(Matrix = as.matrix(status.diff), xLabels = colnames(status.diff), yLabels = colnames(ME.genes), ySymbols = colnames(ME.genes), yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix.traits, setStdMargins = F, zlim = heatmap.range, main = "")
dev.off()

#Eigengene plots
EigengeneBoxplot(ME.genes, lumi.pdata$Status, "black")
EigengeneBoxplot(ME.genes, lumi.pdata$Status, "green")

EstimateViolinPlot(posterior.status$MEgreen, "green", c(-0.05,0))
EstimateViolinPlot(posterior.status$MEblack, "black", c(-0.03,0.025))

#Generate network statistics
all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Symbol = rownames(all.degrees), Module = module.colors, all.degrees, stringsAsFactors = FALSE) %>% arrange(Module)
gene.info$kscaled <- by(gene.info, gene.info$Module, select, kWithin) %>% map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)
gene.info[,3:ncol(gene.info)] <- signif(gene.info[,3:ncol(gene.info)], 3)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = gene.info$Symbol, mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")

gene.info.annot <- left_join(gene.info, bm.table) %>% 
    select(Symbol, Definition, Module:kscaled) %>% 
    filter(!duplicated(Symbol))

gene.module.membership <- data.frame(bicor(expr.collapse, ME.genes, maxPOutliers = 0.05)) %>% signif(3)
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
gene.module.membership$Symbol <- rownames(gene.module.membership)
#module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.collapse)))
#colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
#module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

module.membership <- left_join(gene.info.annot, gene.module.membership) %>% arrange(Module, desc(kscaled)) #%>% left_join(module.membership.pvalue, by = "Symbol")
ModuleWorkbook(module.membership, "module_membership.xlsx")

#Enrichr
source("../../code/GO/enrichr.R")

modules.only <- select(module.membership, Symbol, Module)
color.names <- unique(module.colors) %>% sort
enrichr.terms <- c("GO Biological Process", "GO Molecular Function", "KEGG", "Reactome") 

black.only <- filter(modules.only, Module == "black")
green.only <- filter(modules.only, Module == "green")

black.enrichr <- map(enrichr.terms, GetHyper, black.only$Symbol, modules.only$Symbol)
names(black.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, black.enrichr, "black")

green.enrichr <- map(enrichr.terms, GetHyper, green.only$Symbol, modules.only$Symbol)
names(green.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, green.enrichr, "green")

#Final Enrichr plots
black.gobiol.file <- "./enrichr/black/GO Biological Process.xlsx"
black.gobiol <- read.xlsx(black.gobiol.file) 
GetKappaCluster(file_path_sans_ext(black.gobiol.file), black.gobiol, black.only$Symbol, "Log.Bayes.Factor")
black.gobiol$Database <- "GO BP"

black.gomole.file <- "./enrichr/black/GO Molecular Function.xlsx"
black.gomole <- read.xlsx(black.gomole.file)
black.gomole$Database <- "GO MF"

black.reactome.file <- "./enrichr/black/Reactome.xlsx"
black.reactome <- read.xlsx(black.reactome.file)
GetKappaCluster(file_path_sans_ext(black.reactome.file), black.reactome, black.only$Symbol, "Log.Bayes.Factor")
black.reactome$Database <- "Reactome"

black.kegg.file <- "./enrichr/black/KEGG.xlsx"
black.kegg <- read.xlsx(black.kegg.file)
GetKappaCluster(file_path_sans_ext(black.kegg.file), black.kegg, black.only$Symbol, "Log.Bayes.Factor")
black.kegg$Database <- "KEGG"

black.gobiol.final <- slice(black.gobiol, c(1, 2, 13, 7))
black.gomole.final <- slice(black.gomole, 1)
black.kegg.final <- slice(black.kegg, c(2, 7, 11))
black.reactome.final <- slice(black.reactome, c(2, 6))

black.enrichr <- rbind(black.gobiol.final, black.gomole.final, black.kegg.final, black.reactome.final)
EnrichrPlot(black.enrichr, "black", "Black Module", plot.height = 4, plot.width = 5.5)

green.gobiol.file <- "./enrichr/green/GO Biological Process.xlsx"
green.gobiol <- read.xlsx(green.gobiol.file) 
#green.gobiol.reduce <- slice(green.gobiol, 1:200)
GetKappaCluster(file_path_sans_ext(green.gobiol.file), green.gobiol, green.only$Symbol, "Log.Bayes.Factor")
green.gobiol$Database <- "GO BP"

green.gomole.file <- "./enrichr/green/GO Molecular Function.xlsx"
green.gomole <- read.xlsx(green.gomole.file)
GetKappaCluster(file_path_sans_ext(green.gomole.file), green.gomole, green.only$Symbol, "Log.Bayes.Factor")
green.gomole$Database <- "GO MF"

green.reactome.file <- "./enrichr/green/Reactome.xlsx"
green.reactome <- read.xlsx(green.reactome.file)
GetKappaCluster(file_path_sans_ext(green.reactome.file), green.reactome, green.only$Symbol, "Log.Bayes.Factor")
green.reactome$Database <- "Reactome"

green.kegg.file <- "./enrichr/green/KEGG.xlsx"
green.kegg <- read.xlsx(green.kegg.file)
#GetKappaCluster(file_path_sans_ext(green.kegg.file), green.kegg, green.only$Symbol, "Log.Bayes.Factor")
green.kegg$Database <- "KEGG"

green.gobiol.final <- slice(green.gobiol, c(1, 7, 13))
green.gomole.final <- slice(green.gomole, 4)
green.reactome.final <- slice(green.reactome, c(2, 3))

green.enrichr <- rbind(green.gobiol.final, green.gomole.final, green.reactome.final)
EnrichrPlot(green.enrichr, "green", "green Module", plot.height = 4, plot.width = 5.5, color = "green")

#PPI for whole module
#black.ppi.all <- GetPPI(black.only$Symbol)
#black.all.graph <- graph_from_edgelist(as.matrix(black.ppi.all))
#black.incident <- map(V(black.all.graph), incident, graph = black.all.graph) %>% map_int(length)

#brown.ppi.all <- GetPPI(brown.only$Symbol)
#brown.all.graph <- graph_from_edgelist(as.matrix(brown.ppi.all))
#brown.incident <- map(V(brown.all.graph), incident, graph = brown.all.graph) %>% map_int(length)

##Get genes in GO categories
#brown.genes <- map(brown.enrichr$Genes, str_split, ",") %>% map(extract2, 1) 
#brown.symbols <- reduce(brown.genes, c) %>% unique 

#black.genes <- map(black.enrichr$Genes, str_split, ",") %>% map(extract2, 1) 
#black.symbols <- reduce(black.genes, c) %>% unique 

##PPI
#brown.ppi <-  map(brown.genes, GetPPI)
#black.ppi <-  map(black.genes, GetPPI)

#PlotPPI(adjacency.expr, brown.enrichr$Genes[[4]], brown.ppi[[4]], "brown_test", plot.width = 20, plot.height = 20, vertex.size = 15)

##Kscaled
#black.kscaled <- map(black.genes, GetKscaled, module.membership)
#black.kscaled.df <- reduce(black.kscaled, rbind) %>% distinct
#black.kscaled.df %<>% slice(match(black.symbols, black.kscaled.df$Symbol))

#brown.kscaled <- map(brown.genes, GetKscaled, module.membership)
#brown.kscaled.df <- reduce(brown.kscaled, rbind) %>% distinct
#brown.kscaled.df %<>% slice(match(brown.symbols, brown.kscaled.df$Symbol))

##GO Overlap
#black.intersect <- fromList(black.genes) %>%  data.frame 
#black.weighted <- apply(black.intersect, 2, multiply_by, black.kscaled.df$kscaled) %>% t %>% data.frame
#colnames(black.weighted) <- black.symbols
#black.weighted$Term <- str_replace_all(black.enrichr$Term, "\\ \\(GO.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% str_replace_all(",.*$", "")  
#black.weighted$Format.Name <- str_c(black.enrichr$Database, ": ", black.weighted$Term)
#black.weighted$Format.Name %<>% factor(levels = black.weighted$Format.Name)
#black.plot <- gather(black.weighted, Gene, Kscaled, -Format.Name, -Term)
#black.plot$Gene %<>% factor(levels = black.symbols)

#CairoPDF("black.heatmap", width = 8, height = 30)
#p <- ggplot(black.plot, aes(Format.Name, Gene, fill = Kscaled)) + geom_raster() + theme_bw()
#p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
#p <- p + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(panel.background = element_blank(), panel.border = element_blank())
#p <- p + theme(plot.background = element_blank())
#p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#print(p)
#dev.off()

#brown.intersect <- fromList(brown.genes) %>%  data.frame 
#brown.weighted <- apply(brown.intersect, 2, multiply_by, brown.kscaled.df$kscaled) %>% t %>% data.frame
#colnames(brown.weighted) <- brown.symbols
#brown.weighted$Term <- str_replace_all(brown.enrichr$Term, "\\ \\(GO.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% str_replace_all(",.*$", "")  
#brown.weighted$Format.Name <- str_c(brown.enrichr$Database, ": ", brown.weighted$Term)
#brown.weighted$Format.Name %<>% factor(levels = brown.weighted$Format.Name)
#brown.plot <- gather(brown.weighted, Gene, Kscaled, -Format.Name, -Term)
#brown.plot$Gene %<>% factor(levels = brown.symbols)

#CairoPDF("brown.heatmap", width = 8, height = 30)
#p <- ggplot(brown.plot, aes(Format.Name, Gene, fill = Kscaled)) + geom_raster() + theme_bw()
#p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
#p <- p + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(panel.background = element_blank(), panel.border = element_blank())
#p <- p + theme(plot.background = element_blank())
#p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#print(p)
#dev.off()
