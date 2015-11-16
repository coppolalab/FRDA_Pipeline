#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()

#For baseline processing
library(limma)
library(sva)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(lambda.r)

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

#Reading and writing tables
library(readr)
library(openxlsx)
library(parallel)

CV <- function(data.vector)
{
    coef.var <- sd(data.vector) / mean(data.vector) * 100
    return(coef.var)
}

gen.cv <- function(dataset)
{
    cv <- apply(dataset, 1, CV) %>% as.vector
    dataset.cv <- mutate(data.frame(dataset), Coef.var = cv, Probe_Id = rownames(dataset)) %>% arrange(desc(Coef.var))
    dataset.arr <- select(dataset.cv, -Coef.var) %>% slice(1:25000) 
    rownames(dataset.arr) <- dataset.arr$Probe_Id
    dataset.arr <- select(dataset.arr, -Probe_Id) %>% t
    return(dataset.arr)
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

load(file = "../WGCNA/save/intensities.rmreps.rda")
load(file = "../WGCNA/save/targets.final.known.rda")
load(file = "../WGCNA/save/annot.reduce.rda")
source("../common_functions.R")

#repeat1 <- apply(select(targets.final.gaa, GAA1:GAA2), 1, min)
#repeat2 <- apply(select(targets.final.gaa, GAA1:GAA2), 1, max)
targets.final.known$GAA1 %<>% as.numeric
targets.final.known$GAA2 %<>% as.numeric

targets.final.gaa <- filter(targets.final.known, (!is.na(GAA1) | !is.na(GAA2)) & !is.na(Onset))
hell.conditions <- targets.final.gaa$GAA1 < 40 & !is.na(targets.final.gaa$GAA1)
#targets.final.gaa[hell.conditions,]$GAA1 <- NA
targets.final.gaa[is.na(targets.final.gaa$GAA1),]$GAA1 <- targets.final.gaa[is.na(targets.final.gaa$GAA1),]$GAA2

intensities.gaa <- select(data.frame(intensities.rmreps), one_of(targets.final.gaa$Sample.Name)) %>% as.matrix %>% normalizeBetweenArrays

targets.final.gaa$Sample.Num %<>% str_replace("1r", "1")

plot.cor <- function(expr.row, gaa.repeats)
{
    symbol <- expr.row["Symbol"]
    symbol.index <- match(expr.row["Symbol"], expr.row)
    probe.index <- match(expr.row["Probe_Id"], expr.row)
    expr.row.clean <- expr.row[-c(symbol.index, probe.index)]

    expr.plot <- cbind(GAA = as.numeric(gaa.repeats), Expression = as.numeric(expr.row.clean)) %>% data.frame
    cor.value <- cor(expr.plot$GAA, expr.plot$Expression)
    cor.pvalue <- corPvalueStudent(cor.value, nrow(expr.plot))

    p <- ggplot(expr.plot, aes(x = as.numeric(GAA), y = as.numeric(Expression))) + geom_point()   
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ylab("Log 2 normalized expression") + xlab("GAA Repeat Length") + stat_smooth(method = "lm", se = FALSE)
    p <- p + ggtitle(paste(symbol, '\n', "Correlation:", round(cor.value, 2), "(P < ", format(cor.pvalue, digits = 3), ")"))
    CairoPDF(file.path("./correlations", symbol), height = 6, width = 10)
    print(p)
    dev.off()
}

gaa.cor.spearman <- apply(intensities1.gaa, 1, cor, targets1.gaa$GAA1, method = "spearman")
gaa.cor.bicor <- apply(intensities1.gaa, 1, bicor, targets1.gaa$GAA1, nThreads = 7)

model.sex.full <- model.matrix(~ 0 + factor(targets.final.gaa$Sex))[,-1]
targets.final.gaa$Draw.Age %<>% as.numeric 
draw.age.cov <- log2(targets.final.gaa$Draw.Age)
targets.final.gaa$Diff.Baseline %<>% as.numeric
diff.baseline.cov <- targets.final.gaa$Diff.Baseline
diff.baseline.cov[diff.baseline.cov != 0] %<>% log2
targets.final.gaa$Onset %<>% as.numeric

covariates.full <- cbind(Sex = model.sex.full, Age = draw.age.cov, GAA1 = log2(targets.final.gaa$GAA1), Onset = log2(targets.final.gaa$Onset), Diff.Baseline = diff.baseline.cov, RIN = log2(targets.final.gaa$RIN))

intensities.combat <- ComBat(dat = intensities.gaa, batch = factor(targets.final.gaa$Batch), mod = covariates.full)
saveRDS.gz(intensities.combat, file = "./save/intensities.combat.rda")

#change to 6
gen.peer(10, intensities.combat, TRUE, covariates.full)
PEER.precision <- read_csv("./precision_10.csv")
PEER.precision$Factor <- rownames(PEER.precision)
data.plot <- PEER.precision[-c(1:6),]
colnames(data.plot)[1] <- "Precision"
p <- ggplot(data.plot, aes(x = as.numeric(Factor), y = Precision, group = 1)) + geom_line()
CairoPDF("./precision_reduced")
print(p)
dev.off()

PEER.weights.plot <- select(read_csv("./weight_10.csv"), -(X1:X7))
PEER.weights.sums <- colSums(abs(PEER.weights.plot)) %>% data.frame
PEER.weights.sums$Factor <- 1:nrow(PEER.weights.sums)
colnames(PEER.weights.sums)[1] <- "Weight"

p <- ggplot(PEER.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
p <- p + theme_bw() + xlab("Factor") + ylab("Weigtht")
CairoPDF("./PEER_weights", height = 4, width = 6)
print(p)
dev.off()

PEER.weights <- read_csv("./weight_10.csv") %>% select(-(X3:X5), -X7)
PEER.factors <- read_csv("./factor_10.csv") %>% select(-(X3:X5), -X7)
PEER.residuals <- as.matrix(PEER.factors) %*% t(as.matrix(PEER.weights)) %>% t
intensities.PEER <- intensities.combat - PEER.residuals
saveRDS.gz(intensities.PEER, file = "./save/intensities.peer.rda")

targets1.gaa <- filter(targets.final.gaa, Sample.Num == "1")
intensities1.gaa <- select(data.frame(intensities.PEER), one_of(targets1.gaa$Sample.Name)) %>% as.matrix %>% normalizeBetweenArrays

gaa.cor <- apply(intensities1.gaa, 1, cor, targets1.gaa$GAA1)
gaa.cor.pval <- corPvalueStudent(gaa.cor, length(targets1.gaa$GAA1))
cor.df <- cbind(gaa.cor, gaa.cor.pval) %>% data.frame
colnames(cor.df) <- c("Correlation", "P.value")
cor.df$Probe_Id <- rownames(intensities1.gaa)
cor.df %<>% arrange(P.value)
cor.joined <- join(cor.df, annot.reduce)

cor.probes <- c("ILMN_1726589", "ILMN_1659688", "ILMN_1703524", "ILMN_1765796", "ILMN_2184373", "ILMN_1745172", "ILMN_1747344") %>% paste(collapse = "|")
intensities1.cor <- data.frame(Probe_Id = rownames(intensities1.gaa), intensities1.gaa) %>% filter(grepl(cor.probes, Probe_Id))
annot.symbol <- select(annot.reduce, Probe_Id, Symbol)
intensities1.joined <- join(intensities1.cor, annot.symbol)
apply(intensities1.joined, 1, plot.cor, targets1.gaa$GAA1)
#Calculate coefficient of variation for all genes and then rank them using this
expr.data.PEER <- gen.cv(intensities.PEER)
saveRDS.gz(expr.data.PEER, file = "./save/expr.data.PEER.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

sft.PEER <- pickSoftThreshold(expr.data.PEER, powerVector = powers, verbose = 5, networkType = "signed")
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

softPower <- 7
adjacency.PEER <- adjacency(expr.data.PEER, power = softPower, type = "signed")
saveRDS.gz(adjacency.PEER, file = "./save/adjacency.PEER.rda")

TOM.PEER <- TOMsimilarity(adjacency.PEER)
dissimilarity.TOM <- 1 - TOM.PEER
saveRDS.gz(TOM.PEER, file = "./save/tom.PEER.rda")
saveRDS.gz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./1-genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 50

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, cutHeight = 0.995, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./1-gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.data.PEER, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - cor(ME.genes)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./3-module_eigengene_clustering_min50", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.05
merge.all <- mergeCloseModules(expr.data.PEER, dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3) #PC analysis may be failing because of Intel MKL Lapack routine bug.  Test with openBLAS in R compiled with gcc.
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("4-module_eigengene_clustering_min50", height = 10, width = 15)
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

all.degrees <- intramodularConnectivity(adjacency.PEER, module.colors)
expr.data.t <- mutate(data.frame(t(expr.data.PEER)), Probe_Id = colnames(expr.data.PEER)) 
gene.info.join <- join(expr.data.t, annot.reduce)
gene.info <- select(gene.info.join, Probe_Id:Definition) %>% mutate(module.colors = module.colors, mean.count = apply(expr.data.PEER, 2, mean)) %>% data.frame(all.degrees)

CairoPDF("5-eigengenes", height = 10, width = 18)
par(cex = 0.7)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.adjacency = 0.3, cex.preservation = 0.3, plotPreservation = "standard")
dev.off()

write_csv(data.frame(table(module.colors)), path = "./final_eigengenes.csv") 
gene.info$kscaled <- by(gene.info, gene.info$module.colors, select, kWithin) %>% llply(function(x) { x / max (x) }) %>% reduce(c)
saveRDS.gz(gene.info, file = "./save/gene.info.rda")

gene.module.membership <- as.data.frame(cor(expr.data.PEER, ME.genes, use = "p"))
module.membership.pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.data.PEER)))
names(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")

module.membership <- cbind(select(gene.info, Probe_Id:Definition), gene.module.membership, module.membership.pvalue)
write_csv(module.membership, "module_membership.csv")

colnames(gene.module.membership) %<>% str_replace("MM.", "")
colnames(module.membership.pvalue) %<>% str_replace("MM.pvalue.", "")
gene.module.membership$Probe_Id <- rownames(gene.module.membership)
module.membership.pvalue$Probe_Id <- rownames(module.membership.pvalue)

gene.module.membership.long <- gather(gene.module.membership, module.comparison, correlation, lightgreen:darkred)
module.membership.pvalue.long <- gather(module.membership.pvalue, module.comparison, p.value, lightgreen:darkred)
membership.join <- join(gene.module.membership.long, module.membership.pvalue.long)
eigengene.connectivity <- join(membership.join, gene.info) %>% select(Probe_Id, Accession:kscaled, module.comparison:p.value)
write_csv(eigengene.connectivity, "eigengene_connectivity.csv")

all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y") 
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$x <- as.factor(1:nrow(smooth.df))
smooth.plot <- melt(smooth.df, id.vars = "x")

gen.pcaplot("all_principal_components", smooth.plot, FALSE, 10, 15)
gen.pcaplot("facet_principal_components", smooth.plot, TRUE, 13, 25)

sample.ids <- factor(rownames(expr.data.PEER), levels = rownames(expr.data.PEER))
colnames(ME.genes) %<>% str_replace("ME", "")
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.data.PEER), module.colors)
cluster <- makeForkCluster(8)
split(expr.data.plot, expr.data.plot$module.colors) %>% parLapply(cl = cluster, gen.heatmap, ME.genes.plot)

modules.out <- select(expr.data.plot, module.colors)
modules.out$Probe_Id <- rownames(modules.out)
modules.out %<>% join(annot.reduce)
write.xlsx(modules.out, "modules_out.xlsx")
source("../GO/enrichr.R")

enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 
#enrichr.submit("blue", modules.out, enrichr.terms, FALSE)
color.names <- unique(module.colors) %>% sort
l_ply(color.names, enrichr.submit, modules.out, enrichr.terms, FALSE)

enrichr.submit <- function(index, full.df, enrichr.terms, use.weights)
{
    dataset <- filter(full.df, module.colors == index)
    dir.create(file.path("./enrichr", index), showWarnings = FALSE)
    enrichr.data <- parLapply(cluster, enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    parLapply(cluster, names(enrichr.data), enrichr.wkbk, enrichr.data, index)
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

targets.final.gaa$Sample.Name %<>% str_replace(" ", "")
rownames(ME.genes) <- rownames(expr.data.PEER)
rownames(PEER.factors) <- rownames(expr.data.PEER)
colnames(PEER.factors) <- paste("X", 1:ncol(PEER.factors), sep = "")

targets.age <- select(targets.final.gaa, Sample.Name, Draw.Age) %>% filter(!is.na(Draw.Age))
targets.age$Draw.Age %<>% as.numeric
cor.age <- gen.cor(ME.genes, targets.age)

targets.sex <- filter(targets.final.gaa, Sex != "UNKNOWN")
targets.sex.m <- model.matrix( ~ 0 + factor(targets.sex$Sex) )[,-1] %>% data.frame #%>% mutate(Sample.Name = targets.sex$Sample.Name)
colnames(targets.sex.m) <- "Sex"
targets.sex.m %<>% mutate(Sample.Name = targets.final.gaa$Sample.Name)
cor.sex <- gen.cor(ME.genes, targets.sex.m)

targets.gaa <- select(targets.final.gaa, Sample.Name, GAA1) %>% filter(!is.na(GAA1)) #%>% filter(!is.na(GAA2))
cor.gaa <- gen.cor(ME.genes, targets.gaa)

targets.onset <- select(targets.final.gaa, Sample.Name, Onset) %>% filter(!is.na(Onset))
cor.onset <- gen.cor(ME.genes, targets.onset)

targets.diff <- select(targets.final.gaa, Sample.Name, Diff.Baseline)
cor.diff <- gen.cor(ME.genes, targets.diff)

PEER.factors.plot <- select(read_csv("./factor_10.csv"), -(X1:X7))
names(PEER.factors.plot) <- paste("X", 1:ncol(PEER.factors.plot), sep = "")
PEER.factors.df <- mutate(PEER.factors.plot, Sample.Name = targets.final.gaa$Sample.Name) 
cor.PEER <- gen.cor(ME.genes, PEER.factors.df) %>% data.frame

module.traits.all <- cbind(cor.age, cor.sex, cor.gaa, cor.onset, cor.diff) %>% data.frame
module.traits.pval <- select(module.traits.all, contains("p.value")) %>% as.matrix
module.traits.cor <- select(module.traits.all, -contains("p.value")) %>% as.matrix

#module.PEER.all <- cor.gaa.PEER %>% data.frame
#module.PEER.pval <- select(module.PEER.all, contains("p.value")) %>% as.matrix
#module.PEER.cor <- select(module.PEER.all, -contains("p.value")) %>% as.matrix

cor.PEER.cor <- select(cor.PEER, -contains("p.value")) %>% as.matrix
cor.PEER.pval <- select(cor.PEER, contains("p.value")) %>% as.matrix

module.trait.out <- data.frame(Module = rownames(module.traits.cor), module.traits.cor, module.traits.pval)
#module.PEER.out <- data.frame(Factor = rownames(module.PEER.cor), module.PEER.cor, module.PEER.pval)
cor.PEER.out <- data.frame(Module = rownames(cor.PEER), cor.PEER.cor, cor.PEER.pval)

write_csv(module.trait.out, "module_trait_cor.csv")
#write_csv(module.PEER.out, "module_PEER_cor.csv")
write_csv(cor.PEER.out, "cor_PEER_cor.csv")

text.matrix.traits <- paste(signif(module.traits.cor, 2), '\n(', signif(module.traits.pval, 1), ')', sep = '')
#text.matrix.PEER <- paste(signif(module.PEER.cor, 2), '\n(', signif(module.PEER.pval, 1), ')', sep = '')
text.matrix.PEERcor <- paste(signif(cor.PEER.cor, 2), '\n(', signif(cor.PEER.pval, 1), ')', sep = '')
dim(text.matrix.traits) = dim(module.traits.cor)
#dim(text.matrix.PEER) = dim(module.PEER.cor)
dim(text.matrix.PEERcor) = dim(cor.PEER.cor)

gen.text.heatmap(module.traits.cor, text.matrix.traits, colnames(module.traits.cor), colnames(ME.genes), "", "module-trait relationships")
#gen.text.heatmap(module.PEER.cor, text.matrix.PEER, colnames(module.PEER.cor), rownames(module.PEER.cor), "", "PEER factor-trait relationships")
gen.text.heatmap(cor.PEER.cor, text.matrix.PEERcor, colnames(cor.PEER.cor), colnames(ME.genes), "", "module-PEER factor relationships")

cor.gaa.labeled <- data.frame(Color = rownames(cor.gaa), cor.gaa)
gaa.sig <- filter(data.frame(cor.gaa.labeled), GAA1.p.value < 0.05)
me.genes.gaa <- select(ME.genes, one_of(as.character(gaa.sig$Color)))
me.genes.gaa$GAA <- targets.final.gaa$GAA1
me.gaa.melt <- melt(me.genes.gaa, id.vars = "GAA")

colnames(me.gaa.melt)[2] <- "Module"
me.gaa.melt$Module %<>% as.character
p <- ggplot(me.gaa.melt, aes(x = GAA, y = value, col = Module)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(me.gaa.melt$GAA)))
p <- p + scale_color_manual(values = sort(unique(me.gaa.melt$Module)))
p <- p + facet_wrap(~ Module)
p <- p + theme(legend.position = "none")

CairoPDF("gaa_eigengenes_05", height = 13, width = 20)
print(p)
dev.off()

cor.gaa.labeled <- data.frame(Color = rownames(cor.gaa), cor.gaa)
gaa.sig <- filter(data.frame(cor.gaa.labeled), GAA1.p.value < 0.005)
me.genes.gaa <- select(ME.genes, one_of(as.character(gaa.sig$Color)))
me.genes.gaa$GAA <- targets.final.gaa$GAA1
me.gaa.melt <- melt(me.genes.gaa, id.vars = "GAA")

test.gaa <- select(me.genes.gaa, salmon, GAA)
colnames(me.gaa.melt)[2] <- "Module"
p <- ggplot(me.gaa.melt, aes(x = GAA, y = value, fill = Module, color = Module)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(me.gaa.melt$GAA)))
p <- p + scale_color_manual(values = sort(unique(me.gaa.melt$Module)))
p <- p + facet_wrap(~ Module)
p <- p + theme(legend.position = "none")

CairoPDF("gaa_eigengenes_005", height = 13, width = 20)
print(p)
dev.off()

cor.onset.labeled <- data.frame(Color = rownames(cor.onset), cor.onset)
onset.sig <- filter(data.frame(cor.onset.labeled), Onset.p.value < 0.05)
me.genes.onset <- select(ME.genes, one_of(as.character(onset.sig$Color)))
me.genes.onset$Onset <- targets.final.gaa$Onset
me.onset.melt <- melt(me.genes.onset, id.vars = "Onset")

colnames(me.onset.melt)[2] <- "Module"
me.onset.melt$Module %<>% as.character
p <- ggplot(me.onset.melt, aes(x = as.numeric(Onset), y = value, col = Module)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(me.onset.melt$Onset)))
p <- p + scale_color_manual(values = sort(unique(me.onset.melt$Module)))
p <- p + facet_wrap(~ Module)
p <- p + theme(legend.position = "none")

CairoPDF("onset_eigengenes_05", height = 13, width = 20)
print(p)
dev.off()

cor.onset.labeled <- data.frame(Color = rownames(cor.onset), cor.onset)
onset.sig <- filter(data.frame(cor.onset.labeled), Onset.p.value < 0.005)
me.genes.onset <- select(ME.genes, one_of(as.character(onset.sig$Color)))
me.genes.onset$Onset <- targets.final.gaa$Onset
me.onset.melt <- melt(me.genes.onset, id.vars = "Onset")

colnames(me.onset.melt)[2] <- "Module"
p <- ggplot(me.onset.melt, aes(x = Onset, y = value, fill = Module, color = Module)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(me.onset.melt$Onset)))
p <- p + scale_color_manual(values = sort(unique(me.onset.melt$Module)))
p <- p + facet_wrap(~ Module)
p <- p + theme(legend.position = "none")

CairoPDF("onset_eigengenes_005", height = 13, width = 20)
print(p)
dev.off()

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
