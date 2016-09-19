#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()

#For baseline processing
library(limma)
library(R.utils)
library(lumi)
library(biomaRt)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
#library(igraph)
#library(TeachingDemos)

#Reading and writing tables
library(readr)
library(openxlsx)

#Functional programming
library(magrittr)
library(purrr)

#Data arrangement
library(dplyr)
library(tidyr)
library(broom)

#String operations
library(stringr)

EigengeneANOVA <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- lm(ME ~ Trait, trait.df) %>% anova %>% tidy
    return(trait.anova$p.value[1])
}

EigengeneAOV <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.aov <- aov(ME ~ Trait, trait.df) %>% TukeyHSD 
    aov.contrasts <- rownames(trait.aov[[1]])
    trait.tidy <- tidy(trait.aov) %>% select(estimate, adj.p.value)
    rownames(trait.tidy) <- aov.contrasts
    return(trait.tidy)
}

EigengeneBoxplot <- function(eigengene.df, status.vector, color) {
    color.column <- str_c("ME", color)
    gene.df <- data.frame(Status = status.vector, Expression = as.vector(eigengene.df[[color.column]]))
    gene.df$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(gene.df, aes(x = Status, y = Expression, fill = Status)) + geom_violin() + geom_boxplot(width = 0.1) + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    CairoPDF(color, width = 6, height = 6)
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
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(module.table), rule = "<0.05", style = sig.pvalues)
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

EnrichrSubmit <- function(index, full.df, enrichr.terms, use.weights) {
    dataset <- filter(full.df, module.colors == index)
    dir.create(file.path("./enrichr", index), showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), enrichr.wkbk, enrichr.data, index)
}

EnrichrWorkbook <- function(subindex, full.df, index) {
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

EnrichrPlot <- function(enrichr.df, filename, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
    enrichr.df$Log.pvalue <- -(log10(enrichr.df$P.value))
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- paste(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")", sep = "")
    enrichr.df %<>% arrange(Log.pvalue)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.pvalue) %>% melt(id.vars = "Format.Name") 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = enrichr.df$Format.Name, hjust = "left", aes(y = 0.1)) + coord_flip() + theme_bw() + theme(legend.position = "none") 
    p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
    CairoPDF(filename, height = plot.height, width = plot.width)
    print(p)
    dev.off()
}

FilterEnrichr <- function(enrichr.table, p.value, filename, cluster = FALSE, ebam.df = data.frame()) {
    enrichr.table$Num.Genes <- map(enrichr.table$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
    enrichr.table %<>% filter(Num.Genes > 4) %>% filter(P.value < p.value)
    if(cluster == TRUE) {
        GetKappaCluster(enrichr.table, ebam.df$Symbol, file_path_sans_ext(filename))
    }
    enrichr.table
}

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

source("../common_functions.R")
lumi.import <- readRDS.gz(file = "../baseline_lumi/save/rmcov.collapse.rda")

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

TOM <- TOMsimilarity(adjacency.expr, verbose = 5)
dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
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

CairoPDF("module_eigengene_clustering", height = 10, width = 15)
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

CairoPDF("eigengenes", height = 6, width = 8)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

#Create heatmap of eigengene significance
color.values <- unique(module.colors)
lumi.pdata <- pData(lumi.import)
lumi.pdata$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))
anova.status <- map_dbl(ME.genes, EigengeneANOVA, lumi.pdata$Status) %>% p.adjust("fdr") %>% signif(3)

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
EigengeneBoxplot(ME.genes, lumi.pdata$Status, "brown")

#Generate network statistics
all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Symbol = rownames(all.degrees), Module = module.colors, all.degrees, stringsAsFactors = FALSE) %>% arrange(Module)
gene.info$kscaled <- by(gene.info, gene.info$Module, select, kWithin) %>% map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = gene.info$Symbol, mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")

gene.info.annot <- left_join(gene.info, bm.table) %>% select(Symbol, Definition, Module:kscaled)

gene.module.membership <- data.frame(bicor(expr.collapse, ME.genes, maxPOutliers = 0.05))
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.collapse)))
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

module.membership <- left_join(gene.info.annot, gene.module.membership) %>% left_join(module.membership.pvalue, by = "Symbol")
ModuleWorkbook(module.membership, "module_membership.xlsx")

#PCA plots
all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% map(extract2, "y")
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$Gene <- as.factor(1:nrow(smooth.df))
smooth.plot <- gather(smooth.df, Module, PCA1, -Gene)

PCAPlot("all_principal_components", smooth.plot, FALSE, 10, 15)
PCAPlot("facet_principal_components", smooth.plot, TRUE, 13, 20)

#Heatmaps
sample.ids <- factor(rownames(expr.collapse), levels = rownames(expr.collapse))
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.collapse), Module = module.colors)
split(expr.data.plot, expr.data.plot$Module) %>% map(MEHeatmap, ME.genes.plot)

#Enrichr
source("../../code/GO/enrichr.R")

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "HumanCyc_2016", "NCI-Nature_2016", "Panther_2016") 
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names[-6], enrichr.submit, modules.out, enrichr.terms, FALSE)

green.only <- filter(modules.out, Module == "green")
pink.only <- filter(modules.out, Module == "pink")

#Final plots
green.gobiol.file <- "./enrichr/green/green_GO_Biological_Process_2015.xlsx"
green.gobiol <- read.xlsx(green.gobiol.file) 
green.gobiol$Num.Genes <- map(green.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
green.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(green.gobiol, green.only$Symbol, file_path_sans_ext(green.gobiol.file))
green.gobiol$Database <- "GO Biological Process"
green.gobiol.final <- slice(green.gobiol, c(1, 25, 35))

green.molec <- read.xlsx("./enrichr/green/green_GO_Molecular_Function_2015.xlsx") %>% slice(14)
green.molec$Database <- "GO Molecular Function"
green.molec$Num.Genes <- map(green.molec$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)

green.reactome <- read.xlsx("./enrichr/green/green_Reactome_2016.xlsx") %>% slice(2)
green.reactome$Database <- "Reactome"
green.reactome$Num.Genes <- map(green.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)

green.enrichr <- rbind(green.gobiol.final, green.molec, green.reactome)
gen.enrichrplot(green.enrichr, "green.enrichr")

pink.gobiol.file <- "./enrichr/pink/pink_GO_Biological_Process_2015.xlsx"
pink.gobiol <- read.xlsx(pink.gobiol.file) 
pink.gobiol$Num.Genes <- map(pink.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pink.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
get.kappa.cluster(pink.gobiol, pink.only$Symbol, file_path_sans_ext(pink.gobiol.file))
pink.gobiol$Database <- "GO Biological Process"
pink.gobiol.final <- slice(pink.gobiol, c(4, 53, 41, 6))

pink.reactome <- read.xlsx("./enrichr/pink/pink_Reactome_2016.xlsx") %>% slice(1)
pink.reactome$Database <- "Reactome"
pink.reactome$Num.Genes <- map(pink.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)

pink.enrichr <- rbind(pink.gobiol.final, pink.reactome)

gen.enrichrplot(pink.gobiol.final, "pink.enrichr")

green.gobiol.resp <- str_split(green.gobiol.final[2,]$Genes, ",")[[1]]
pink.gobiol.toll <- str_split(pink.gobiol.final[4,]$Genes, ",")[[1]]
pink.gobiol.apop <- str_split(pink.gobiol.final[2,]$Genes, ",")[[1]]

ggr.ppi <- get.ppi(green.gobiol.resp)
pgt.ppi <- get.ppi(pink.gobiol.toll)
pga.ppi <- get.ppi(pink.gobiol.apop)

plot.ppi(adjacency.expr, green.gobiol.resp, ggr.ppi, "green_igraph", TRUE)
plot.ppi(adjacency.expr, pink.gobiol.apop, pga.ppi, "pga_igraph", TRUE)
#plot.ppi(adjacency.expr, pink.gobiol.toll, pgt.unique.final, "pgt_igraph", FALSE)

test.minet <- minet(expr.collapse, method = "aracne")
