library(WGCNA)
library(flashClust)
enableWGCNAThreads()

library(lumi)
library(lumiHumanAll.db)
library(annotate)
library(biomaRt)
library(BayesFactor)
library(matrixStats)

library(Cairo)
library(heatmap.plus)
library(openxlsx)

library(tools)
library(R.utils)
library(broom)
library(stringr)
library(magrittr)
library(pryr)
library(tidyverse)
library(forcats)

EigengeneBayes <- function(gene.vector, trait.df) {
    trait.df <- data.frame(trait.df, Gene = gene.vector)
    set.seed(12345)
    trait.anova <- lmBF(Gene ~ FDS + Sex + RIN, data = trait.df) 
    notrait.anova <- lmBF(Gene ~ Sex + RIN, data = trait.df) 
    combined <- c(trait.anova, notrait.anova)
    combined
}

EigengeneScatterplot <- function(eigengene.df, fds.vector, module.size, color) {
    color.column <- str_c("ME", color)
    gene.df <- data.frame(FDS = fds.vector, Expression = as.vector(eigengene.df[[color.column]]))

    p <- ggplot(gene.df, aes(x = FDS, y = Expression)) + 
         geom_point() + 
         stat_smooth(method = "lm", se = TRUE) +
         theme_bw() + 
         theme(legend.position = "none", 
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.border = element_rect(color = "black", size = 1),
               axis.ticks.x = element_blank(),
               plot.margin = unit(c(1,1,1,1), "lines"),
               plot.background = element_blank(),
               plot.title = element_text(hjust = 0.5)) +
         xlab("Functional Stage") + ylab("Eigengene Value") + 
         ggtitle(str_c(capitalize(color), " Module (", module.size, " Transcripts)")) 

    CairoPDF(str_c(color, ".pdf"), width = 5, height = 5, bg = "transparent")
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

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Log.Bayes.Factor)) + 
         geom_bar(stat = "identity") + 
         coord_flip() + 
         theme_bw() + 
         theme(legend.position = "none", 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(color = "black", size = 1),
               plot.background = element_blank(),
               plot.title = element_text(hjust = 0.5),
               axis.title.y = element_blank(), 
               axis.ticks.y = element_blank()) + 
        ylab("logBF") + 
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
combat.collapse <- ReadRDSgz("../phenotype_regression/save/combat.collapse.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.data <- t(exprs(combat.collapse))
sft <- pickSoftThreshold(expr.data, powerVector = powers, verbose = 5, 
                         networkType = "signed", corFnc = bicor, 
                         corOptions = list(maxPOutliers = 0.05))
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
adjacency.expr <- adjacency(expr.data, power = softPower, type = "signed", 
                            corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
gc()

TOM <- TOMsimilarity(adjacency.expr)
dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
SaveRDSgz(geneTree, file = "./save/gene.tree.rda")
gc()

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
ME.list <- moduleEigengenes(expr.data, colors = dynamic.colors, 
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
merge.all <- mergeCloseModules(expr.data, dynamic.colors, 
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

#Determine eigengene significance
color.values <- unique(module.colors)
lumi.pdata <- pData(combat.collapse)

model.cov <- model.matrix(~ Age, data = lumi.pdata)[,-1] %>% t #Create model matrix of covariates to be removed
covar.bayes <- map(data.frame(ME.genes), CovarBayes, select(lumi.pdata, Age))
posterior.covar <- map(covar.bayes, posterior, iterations = 10000) 
covar.coefs <- map(posterior.covar, colMedians) %>% reduce(rbind)
colnames(covar.coefs) <- colnames(posterior.covar[[1]])

ME.genes.rmcov <- data.frame(ME.genes - (covar.coefs[,c("Age")] %*% model.cov))
SaveRDSgz(ME.genes.rmcov, "./save/ME.genes.rmcov.rda")

bayes.fds <- map(ME.genes.rmcov, EigengeneBayes, select(lumi.pdata, FDS, Sex, RIN)) 
bf.fds <- map(bayes.fds, extractBF) %>% 
    map(extract2, "bf") %>% 
    map_dbl(reduce, divide_by) %>%
    map_dbl(log10)
set.seed(12345)
posterior.fds <- map(bayes.fds, posterior, index = 1, iterations = 10000) 
posterior.fds.medians <- map(posterior.fds, colMedians) %>% reduce(rbind) %>% 
    set_colnames(make.names(colnames(posterior.fds[[1]]))) %>%
    as_tibble %>% mutate(Module = names(posterior.fds), Log.Bayes.Factor = bf.fds) %>%
    select(Module, Log.Bayes.Factor, mu:g_continuous) %>% mutate_if(is.numeric, signif, digits = 3) %>%
    arrange(desc(Log.Bayes.Factor))

#Generate network statistics
all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Symbol = rownames(all.degrees), Module = module.colors, 
                        all.degrees, stringsAsFactors = FALSE) %>% arrange(Module)
gene.info$kscaled <- by(gene.info, gene.info$Module, select, kWithin) %>% 
    map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)
gene.info[3:ncol(gene.info)] <- signif(gene.info[3:ncol(gene.info)], 3)
gene.info.annot <- arrange(gene.info, Module, desc(kscaled)) %>% 
    left_join(bm.table) %>% 
    select(Symbol, Definition, Module:kscaled)

gene.module.membership <- data.frame(bicor(expr.data, ME.genes, maxPOutliers = 0.05)) %>% 
    signif(3)
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
gene.module.membership$Symbol <- rownames(gene.module.membership)

module.membership <- left_join(gene.info.annot, gene.module.membership) 
ModuleWorkbook(module.membership, "module_membership.xlsx")

#text.matrix.traits <- str_c(signif(cor.gaa$Correlation, 2), '\n(', signif(cor.gaa$Adj.P.value, 1), ')')
#dim(text.matrix.traits) = dim(as.matrix(cor.gaa$Correlation))

#CairoPDF("module_trait_relationships", width = 4, height = 10)
#par(mar = c(8, 8, 3, 3))
#labeledHeatmap(Matrix = as.matrix(cor.gaa$Correlation), xLabels = "GAA1", yLabels = colnames(ME.genes), 
               #ySymbols = colnames(ME.genes), yColorLabels = TRUE, colors = greenWhiteRed(50), 
               #textMatrix = text.matrix.traits, zlim = c(-1,1), main = "")
#dev.off()

#Eigengene plots
modules.only <- select(module.membership, Module, Symbol)
enrichr.terms <- c("GO Biological Process", "GO Molecular Function", "KEGG", "Reactome") 
color.names <- unique(modules.only$Module) %>% sort

red.only <- filter(modules.only, Module == "red")
yellow.only <- filter(modules.only, Module == "yellow")
magenta.only <- filter(modules.only, Module == "magenta")

EigengeneScatterplot(ME.genes.rmcov, lumi.pdata$FDS, nrow(red.only), "red")
EigengeneScatterplot(ME.genes.rmcov, lumi.pdata$FDS, nrow(yellow.only), "yellow")
EigengeneScatterplot(ME.genes.rmcov, lumi.pdata$FDS, nrow(magenta.only), "magenta")

#Enrichr
source("../../code/GO/enrichr.R")

red.enrichr <- map(enrichr.terms, GetHyper, red.only$Symbol, modules.only$Symbol)
names(red.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, red.enrichr, "red")

yellow.enrichr <- map(enrichr.terms, GetHyper, yellow.only$Symbol, modules.only$Symbol)
names(yellow.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, yellow.enrichr, "yellow")

magenta.enrichr <- map(enrichr.terms, GetHyper, magenta.only$Symbol, modules.only$Symbol)
names(magenta.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, magenta.enrichr, "magenta")

red.gobiol.file <- "./enrichr/red/GO Biological Process.xlsx"
red.gobiol <- read.xlsx(red.gobiol.file) 
red.gobiol$Database <- "GO BP"
GetKappaCluster(file_path_sans_ext(red.gobiol.file), red.gobiol, red.only$Symbol, "Log.Bayes.Factor")

red.gomole.file <- "./enrichr/red/GO Molecular Function.xlsx"
red.gomole <- read.xlsx(red.gomole.file)
red.gomole$Database <- "GO MF"

red.reactome.file <- "./enrichr/red/Reactome.xlsx"
red.reactome <- read.xlsx(red.reactome.file)
red.reactome$Database <- "Reactome"
GetKappaCluster(file_path_sans_ext(red.reactome.file), red.reactome, red.only$Symbol, "Log.Bayes.Factor")

red.kegg.file <- "./enrichr/red/KEGG.xlsx"
red.kegg <- read.xlsx(red.kegg.file)
red.kegg$Database <- "KEGG"

red.gobiol.final <- slice(red.gobiol, c(1,7))
red.gomole.final <- slice(red.gomole, 2)
red.kegg.final <- slice(red.kegg, 1)
red.reactome.final <- slice(red.reactome, c(2,5))

red.enrichr.final <- rbind(red.gobiol.final, red.gomole.final, red.kegg.final, red.reactome.final)
EnrichrPlot(red.enrichr.final, "red", "Red Module", plot.height = 3, plot.width = 5.5, color = "black")

yellow.gobiol.file <- "./enrichr/yellow/GO Biological Process.xlsx"
yellow.gobiol <- read.xlsx(yellow.gobiol.file) 
yellow.gobiol$Database <- "GO BP"
GetKappaCluster(file_path_sans_ext(yellow.gobiol.file), yellow.gobiol, yellow.only$Symbol, "Log.Bayes.Factor")

yellow.gomole.file <- "./enrichr/yellow/GO Molecular Function.xlsx"
yellow.gomole <- read.xlsx(yellow.gomole.file)
yellow.gomole$Database <- "GO MF"
GetKappaCluster(file_path_sans_ext(yellow.gomole.file), yellow.gomole, yellow.only$Symbol, "Log.Bayes.Factor")

yellow.reactome.file <- "./enrichr/yellow/Reactome.xlsx"
yellow.reactome <- read.xlsx(yellow.reactome.file)
yellow.reactome$Database <- "Reactome"
GetKappaCluster(file_path_sans_ext(yellow.reactome.file), yellow.reactome, yellow.only$Symbol, "Log.Bayes.Factor")

yellow.kegg.file <- "./enrichr/yellow/KEGG.xlsx"
yellow.kegg <- read.xlsx(yellow.kegg.file)
yellow.kegg$Database <- "KEGG"

yellow.gobiol.final <- slice(yellow.gobiol, c(1,3,7,9,18))
yellow.gomole.final <- slice(yellow.gomole, c(6,11))
yellow.kegg.final <- slice(yellow.kegg, 7)
yellow.reactome.final <- slice(yellow.reactome, c(1,3,6))

yellow.enrichr.final <- rbind(yellow.gobiol.final, yellow.gomole.final, yellow.kegg.final, yellow.reactome.final)
EnrichrPlot(yellow.enrichr.final, "yellow", "Yellow Module", plot.height = 4, plot.width = 5, color = "black")

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

magenta.gobiol.final <- slice(magenta.gobiol, c(1,2,4,6,8))
magenta.gomole.final <- slice(magenta.gomole, 1)
magenta.reactome.final <- slice(magenta.reactome, c(1,3))
magenta.kegg.final <- slice(magenta.kegg, c(2,4))

magenta.enrichr.final <- rbind(magenta.gobiol.final, magenta.gomole.final, magenta.reactome.final, magenta.kegg.final)
EnrichrPlot(magenta.enrichr.final, "magenta", "Magenta Module", plot.height = 4, plot.width = 5, color = "black")
