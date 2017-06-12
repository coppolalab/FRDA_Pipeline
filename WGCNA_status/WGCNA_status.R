library(WGCNA)
library(flashClust)
enableWGCNAThreads()
library(BayesFactor)
library(limma)
library(R.utils)
library(lumi)
library(biomaRt)
library(matrixStats)

library(igraph)
library(UpSetR)
library(openxlsx)
library(Cairo)

library(tools)
library(broom)
library(stringr)
library(magrittr)
library(tidyverse)

EigengeneBoxplot <- function(eigengene.df, status.vector, color) {
    color.column <- str_c("ME", color)
    gene.df <- data.frame(Status = status.vector, Expression = as.vector(eigengene.df[[color.column]]))
    gene.df$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(gene.df, aes(x = Status, y = Expression, fill = Status)) + 
         geom_violin(scale = "width", trim = FALSE) + 
         geom_boxplot(width = 0.25, outlier.shape = NA) + 
         theme_bw() + 
         theme(legend.position = "none", 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1), 
              plot.background = element_blank(), 
              plot.margin = unit(c(1,1,1,1), "lines")
              axis.title.x = element_blank()) + 
         ggtitle(str_c(capitalize(color), " Module")) +
         ylab("Eigengene Value")

    CairoPDF(color, width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

DEBayes <- function(gene.vector, trait.df) {
    trait.df <- data.frame(trait.df, Gene = gene.vector)
    set.seed(12345)
    trait.anova <- lmBF(Gene ~ Status + RIN, data = trait.df) 
    notrait.anova <- lmBF(Gene ~ RIN, data = trait.df) 
    combined <- c(trait.anova, notrait.anova)
    combined
}

EstimateViolinPlot <- function(estimate.df, color, module.size, ylimits) {
    estimate.extract <- as.matrix(estimate.df) %>% data.frame 
    estimate.groups <- select(estimate.extract, matches("^Status")) 
    estimate.add <- sweep(estimate.groups, 1, estimate.extract$mu, "+")
    estimate.plot <- gather(estimate.add, Group, Estimate) 
    estimate.plot$Group %<>% str_replace("Status\\.", "") %>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(estimate.plot, aes(x = Group, y = Estimate, fill = Group)) + 
         geom_violin(scale = "width", trim = FALSE, draw_quantiles = c(0.05, 0.5, 0.95)) + 
         theme_bw() + 
         theme(legend.position = "none", 
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(), 
             panel.border = element_rect(color = "black", size = 1),
             plot.background = element_blank(), 
             plot.title = element_text(hjust = 0.5),
             axis.title.x = element_blank(), 
             axis.ticks.x = element_blank()) + ylim(ylimits) + 
         ylab("Eigengene Estimate") + 
         ggtitle(str_c(capitalize(color), " Module (", module.size, " Transcripts)")) 

    CairoPDF(str_c(color, "estimate", sep = "."), width = 5, height = 4.5, bg = "transparent")
    print(p)
    dev.off()
}

ModuleWorkbook <- function(module.table, filename) {
    MM.key <- colnames(module.table) %>% str_detect("MM.") 

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = module.table)

    conditionalFormatting(wb, 1, cols = MM.cols, rows = 1:nrow(module.table), 
                          style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:ncol(module.table), widths = "auto")

    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
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
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% 
        map(unlist) %>% map_int(length)
    enrichr.df$Term %<>% str_replace_all("\\ \\(GO.*$", "") %>% 
        str_replace_all("\\_Homo.*$", "") %>% 
        str_replace_all(",.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
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

    p <- p + coord_flip() + 
         theme_bw() + 
         theme(legend.position = "none", 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(color = "black", size = 1),
               plot.background = element_blank(),
               plot.title = element_text(hjust = 0.5),
               axis.title.y = element_blank(), 
               axis.ticks.y = element_blank()) + 
        ylab(expression(paste(Log[10], ' Bayes Factor'))) + 
        ggtitle(plot.title)

    CairoPDF(str_c(filename, ".enrichr"), height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

GetKscaled <- function(gene.list, module.membership) {
    filter(module.membership, is.element(Symbol, gene.list)) %>% 
        select(Symbol, kscaled) %>% 
        arrange(desc(kscaled))
}

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

source("../../code/common_functions.R")
lumi.import <- ReadRDSgz(file = "../baseline_lumi/save/combat.collapse.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.collapse <- exprs(lumi.import) %>% t
sft <- pickSoftThreshold(expr.collapse, powerVector = powers, 
                         verbose = 5, networkType = "signed", 
                         corFnc = bicor, corOptions = list(maxPOutliers = 0.05))
sft.df <- sft$fitIndices
SaveRDSgz(sft, file = "./save/sft.rda")

#Plot scale indendence and mean connectivity as functions of power
sft.df$multiplied <- sft.df$SFT.R.sq * -sign(sft.df$slope)
p <- ggplot(sft.df, aes(x = Power,  y = multiplied, label = Power)) + 
     geom_point() + 
     geom_text(vjust = -0.6, size = 4, col = "red") + 
     geom_hline(aes(yintercept = 0.9)) +
     theme_bw() + 
     theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
     xlab("Soft Threshold") + 
     ylab("Scale Free Topology Model Fit, signed R^2") + 
     ggtitle("Scale Independence")

CairoPDF(file = "./scaleindependence", width = 6, height = 6)
print(p)
dev.off()

p <- ggplot(sft.df, aes(x = Power,  y = mean.k.)) + 
     geom_point() + 
     theme_bw() + 
     theme(panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           plot.title = element_text(hjust = 0.5)) + 
     xlab("Soft Threshold") + 
     ylab("Mean Connectivity") + 
     ggtitle("Mean Connectivity")

CairoPDF(file = "./meanconnectivity", width = 6, height = 6)
print(p)
dev.off()

softPower <- 6
adjacency.expr <- adjacency(expr.collapse, power = softPower, type = "signed", 
                            corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
gc()

TOM <- TOMsimilarity(adjacency.expr, verbose = 5)
dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
gc()
SaveRDSgz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", 
     main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)
dev.off()

#Identify modules using dynamic tree cutting with hybrid clustering
min.module.size <- 20
dynamic.modules <- cutreeDynamic(dendro = geneTree, cutHeight = 0.995, 
                                 method = "hybrid", distM = dissimilarity.TOM, 
                                 pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
SaveRDSgz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.collapse, colors = dynamic.colors, 
                            softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
SaveRDSgz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_tree", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(expr.collapse, dynamic.colors, 
                               cutHeight = ME.dissimilarity.threshold, verbose = 3, 
                               corFnc = bicor, corOptions = list(maxPOutliers = 0.05)) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering", height = 8, width = 10, bg = "transparent")
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("", ""), 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, main = "All subjects", axes = FALSE, 
                    ylab = "", cex.main = 4.0)
dev.off()

CairoPDF("module_eigengene_clustering_new", height = 8, width = 10, bg = "transparent")
plot(geneTree)
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
SaveRDSgz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
SaveRDSgz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
SaveRDSgz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 8)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

#Create heatmap of eigengene significance
color.values <- unique(module.colors)
lumi.pdata <- pData(lumi.import)
lumi.pdata$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))

model.cov <- model.matrix(~ Sex + Age, data = lumi.pdata)[,-1] %>% t #Create model matrix of covariates to be removed
covar.bayes <- map(data.frame(ME.genes), CovarBayes, select(lumi.pdata, Sex, Age))
posterior.covar <- map(covar.bayes, posterior, iterations = 10000) 
covar.coefs <- map(posterior.covar, colMedians) %>% reduce(rbind)
colnames(covar.coefs) <- colnames(posterior.covar[[1]])

ME.genes.rmcov <- data.frame(ME.genes - (covar.coefs[,c("Sex-MALE","Age-Age")] %*% model.cov))
SaveRDSgz(ME.genes.rmcov, "./save/ME.genes.rmcov.rda")

bayes.status <- map(ME.genes.rmcov, DEBayes, select(lumi.pdata, Status, RIN)) 
bf.status <- map(bayes.status, extractBF) %>% 
    map(extract2, "bf") %>% 
    map_dbl(reduce, divide_by) %>%
    map_dbl(log10)
posterior.status <- map(bayes.status, posterior, index = 1, iterations = 10000) 

#Generate network statistics
all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Symbol = rownames(all.degrees), 
                        Module = module.colors, all.degrees, stringsAsFactors = FALSE) %>% arrange(Module)
gene.info$kscaled <- by(gene.info, gene.info$Module, select, kWithin) %>% 
    map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)
gene.info[,3:ncol(gene.info)] <- signif(gene.info[,3:ncol(gene.info)], 3)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), 
                  filters = 'hgnc_symbol', 
                  values = gene.info$Symbol, 
                  mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")

gene.info.annot <- left_join(gene.info, bm.table) %>% 
    select(Symbol, Definition, Module:kscaled) %>% 
    filter(!duplicated(Symbol))

gene.module.membership <- data.frame(bicor(expr.collapse, ME.genes, maxPOutliers = 0.05)) %>% signif(3)
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
gene.module.membership$Symbol <- rownames(gene.module.membership)

module.membership <- left_join(gene.info.annot, gene.module.membership) %>% arrange(Module, desc(kscaled)) 
ModuleWorkbook(module.membership, "module_membership.xlsx")

#text.matrix.traits <- str_c(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(pval.adjust), 1), ')')
#dim(text.matrix.traits) = dim(status.diff)

#heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
#width.dynamic <- 3 + (1 * ncol(text.matrix.traits))

#CairoPDF("module_trait_relationships", width = width.dynamic, height = 10)
#par(mar = c(8, 8, 3, 3))
#labeledHeatmap(Matrix = as.matrix(status.diff), xLabels = colnames(status.diff), 
               #yLabels = colnames(ME.genes), ySymbols = colnames(ME.genes), 
               #yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix.traits, 
               #setStdMargins = F, zlim = heatmap.range, main = "")
#dev.off()

#Eigengene plots
modules.only <- select(module.membership, Symbol, Module)
magenta.only <- filter(modules.only, Module == "magenta")

EigengeneBoxplot(ME.genes.rmcov, lumi.pdata$Status, "magenta")
posterior.magenta.plot <- filter(data.frame(posterior.status$MEmagenta), Status.Patient > -0.01)
EstimateViolinPlot(posterior.magenta.plot, "magenta", nrow(magenta.only), c(-0.03, 0.015))

#Enrichr
source("../../code/GO/enrichr.R")

color.names <- unique(module.colors) %>% sort
enrichr.terms <- c("GO Biological Process", "GO Molecular Function", "KEGG", "Reactome") 

magenta.enrichr <- map(enrichr.terms, GetHyper, magenta.only$Symbol, modules.only$Symbol)
names(magenta.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, magenta.enrichr, "magenta")

#Final Enrichr plots
magenta.gobiol.file <- "./enrichr/magenta/GO Biological Process.xlsx"
magenta.gobiol <- read.xlsx(magenta.gobiol.file) 
magenta.gobiol$Database <- "GO BP"
GetKappaCluster(file_path_sans_ext(magenta.gobiol.file), magenta.gobiol, magenta.only$Symbol, "Log.Bayes.Factor")

magenta.gomole.file <- "./enrichr/magenta/GO Molecular Function.xlsx"
magenta.gomole <- read.xlsx(magenta.gomole.file)
magenta.gomole$Database <- "GO MF"

magenta.reactome.file <- "./enrichr/magenta/Reactome.xlsx"
magenta.reactome <- read.xlsx(magenta.reactome.file)
magenta.reactome$Database <- "Reactome"

magenta.kegg.file <- "./enrichr/magenta/KEGG.xlsx"
magenta.kegg <- read.xlsx(magenta.kegg.file)
magenta.kegg$Database <- "KEGG"

magenta.gobiol.final <- slice(magenta.gobiol, c(1,2,9,12))
magenta.gomole.final <- slice(magenta.gomole, 2)
magenta.kegg.final <- slice(magenta.kegg, c(3,6))
magenta.reactome.final <- slice(magenta.reactome, c(3,12))

magenta.enrichr <- rbind(magenta.gobiol.final, magenta.gomole.final, magenta.kegg.final, magenta.reactome.final)
EnrichrPlot(magenta.enrichr, "magenta", "Magenta Module", plot.height = 4, plot.width = 5.5)
