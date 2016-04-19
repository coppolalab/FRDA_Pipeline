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

gen.enrichrplot <- function(enrichr.df, filename, plot.height = 5, plot.width = 8)
{
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

lumi.import <- readRDS.gz(file = "../baseline_lumi/save/rmcov.collapse.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.collapse <- exprs(lumi.import) %>% t
sft <- pickSoftThreshold(expr.collapse, powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
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
saveRDS.gz(adjacency, file = "./save/adjacency.rda")

TOM <- TOMsimilarity(adjacency.expr, verbose = 5)
dissimilarity.TOM <- 1 - TOM
saveRDS.gz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 20

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
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

CairoPDF("eigengenes", height = 6, width = 8)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
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
anova.status <- map_dbl(ME.genes, eigengene.anova, targets.final.known$Status) %>% p.adjust("fdr") %>% signif(3)

cor.status <- map(ME.genes, eigengene.model, targets.final.known$Status)
status.diff <- map(cor.status, select, diff) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.diff) <- names(cor.status)
status.pval <- map(cor.status, select, p.adj) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.pval) <- names(cor.status)

pval.adjust <- map(status.pval, p.adjust, method = "fdr") %>% reduce(cbind) %>% data.frame
rownames(pval.adjust) <- rownames(status.pval)
colnames(pval.adjust) <- paste(colnames(status.pval), ".pval")
status.pval$Module <- rownames(status.pval)

text.matrix.traits <- paste(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(status.pval), 1), ')', sep = '')
dim(text.matrix.traits) = dim(status.diff)

heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
gen.text.heatmap(as.matrix(status.diff), text.matrix.traits, colnames(status.diff), colnames(ME.genes), "", "module-trait relationships", heatmap.range)

all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
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
split(expr.data.plot, expr.data.plot$module.colors) %>% map(gen.heatmap, ME.genes.plot)
by(expr.data.plot, expr.data.plot$module.colors, gen.heatmap, ME.genes.plot)

#modules.out <- select(gene.info, Symbol, module.color)
#targets.final.known$Sample.Name %<>% str_replace(" ", "")

source("../../code/GO/enrichr.R")

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "HumanCyc_2016", "NCI-Nature_2016", "Panther_2016") 
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names[-6], enrichr.submit, modules.out, enrichr.terms, FALSE)

#skipped
split(modules.out, modules.out$module.color) %>% map(submit.stringdb)

submit.stringdb <- function(module.subset)
{
    catch <- get.stringdb(module.subset, unique(module.subset$module.color), "./stringdb")
}

rownames(ME.genes) <- rownames(expr.collapse)

source("../common_functions.R")

ME.df <- data.frame(Status = targets.final.known$Status, ME.genes)
p <- ggplot(ME.df, aes(factor(Status), magenta, col = "magenta") ) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), legend.position = "none")
p <- p + ylab("Eigengene") + ggtitle(paste("P < ", anova.status["magenta"]))
CairoPDF("magenta.eigengene", width = 4, height = 4)
print(p)
dev.off()

p <- ggplot(ME.df, aes(factor(Status), black ) ) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), legend.position = "none")
p <- p + ylab("Eigengene") + ggtitle(paste("P < ", anova.status["black"]))
CairoPDF("black.eigengene", width = 4, height = 4)
print(p)
dev.off()

plot.eigencor(module.traits.pval, "Control.Carrier", lumi.import$Status)
plot.eigencor(module.traits.pval, "Control.Patient", lumi.import$Status)
plot.eigencor(module.traits.pval, "Carrier.Patient", lumi.import$Status)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#Final plots
black.biol <- read.xlsx("./enrichr/black/black_GO_Biological_Process_2015.xlsx") %>% slice(c(3,31))
black.biol$Database <- "GO Biological Process"
black.molec <- read.xlsx("./enrichr/black/black_GO_Molecular_Function_2015.xlsx") %>% slice(c(1,2))
black.molec$Database <- "GO Molecular Function"
black.kegg <- read.xlsx("./enrichr/black/black_KEGG_2016.xlsx") %>% slice(c(5,6,13))
black.kegg$Database <- "KEGG"
black.reactome <- read.xlsx("./enrichr/black/black_Reactome_2016.xlsx") %>% slice(13)
black.reactome$Database <- "Reactome"
black.enrichr <- rbind(black.biol, black.molec, black.kegg, black.reactome)

gen.enrichrplot(black.enrichr, "black.enrichr")

magenta.biol <- read.xlsx("./enrichr/magenta/magenta_GO_Biological_Process_2015.xlsx") %>% slice(c(3,6))
magenta.biol$Database <- "GO Biological Function"
magenta.molec <- read.xlsx("./enrichr/magenta/magenta_GO_Molecular_Function_2015.xlsx") %>% slice(13)
magenta.molec$Database <- "GO Molecular Function"
magenta.reactome <- read.xlsx("./enrichr/magenta/magenta_Reactome_2016.xlsx") %>% slice(c(2,3))
magenta.reactome$Database <- "Reactome"
magenta.enrichr <- rbind(magenta.molec, magenta.biol, magenta.reactome)

gen.enrichrplot(magenta.enrichr, "magenta.enrichr")
