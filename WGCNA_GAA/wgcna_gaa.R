#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()

#For baseline processing
library(limma)
library(sva)
library(lumi)
library(lumiHumanAll.db)
library(annotate)
library(R.utils)

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

plot.eigencor <- function(module.traits.pval, col.name, pheno.vector)
{
    sig.col <- paste(col.name, ".p.value", sep = "")
    cor.status.labeled <- data.frame(Color = rownames(module.traits.pval), select_(data.frame(module.traits.pval), sig.col))
    filter.cond <- paste(sig.col, "< 0.05")
    colors.sig <- filter_(cor.status.labeled, filter.cond)
    me.genes.subset <- select(ME.genes, one_of(as.character(colors.sig$Color)))
    me.genes.subset[[col.name]] <- pheno.vector
    me.subset.melt <- melt(me.genes.subset, id.vars = col.name) 
    colnames(me.subset.melt)[2] <- "Module"

    me.subset.melt$Module %<>% as.character
    p <- ggplot(me.subset.melt, aes_string(x = col.name, y = "value", col = "Module")) + geom_point(position = "jitter")
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene") + stat_smooth(method = "lm", se = FALSE)
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(pheno.vector)))
    p <- p + scale_color_manual(values = sort(unique(me.subset.melt$Module)))
    p <- p + facet_wrap(~ Module, scales = "free_y") + stat_smooth(method = "lm", se = TRUE)
    p <- p + theme(legend.position = "none")

    filename <- paste(col.name, "_eigengenes_05", sep = "")
    CairoPDF(filename, height = 13, width = 20)
    print(p)
    dev.off()
}

gen.boxplot <- function(filename, lumi.object, colorscheme, maintext, ylabtext)
{
    #dataset %<>% t %>% data.frame
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Batch = lumi.object$Batch)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Batch"))
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value, fill = factor(Batch))) + geom_boxplot() + theme_bw()
    p <- p + scale_fill_manual(values = colorscheme)
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Oxygen", width = 20 , height = 8)
}

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

lumi.import <- readRDS.gz(file = "../baseline_lumi/save/lumi.baseline.rda")
remove.all <- readRDS.gz(file = "../baseline_lumi/save/remove.all.rda")
source("../common_functions.R")

lumi.rmreps <- lumi.import[,remove.all] #Remove the previously dropped samples
patient.key <- !(is.na(lumi.rmreps$GAA1)) & !(is.na(lumi.rmreps$Onset)) #Select only patients with known GAA repeat length and Onset
batch.key <- !(lumi.rmreps$Batch == 8) & !(lumi.rmreps$Batch == 14) #Remove batches 8 and 14 because they only have one array each and this messes up ComBat
rin.key <- !(lumi.rmreps$RIN < 5.0) #Remove low RIN samples.  Needs to be done at the main preprocessing level instead
remove.new <- patient.key & batch.key & rin.key
lumi.patient <- lumi.rmreps[,remove.new]

saveRDS.gz(lumi.patient, file = "./save/lumi.patient.rda")
lumi.patient.norm <- lumiN(lumi.patient, method = "rsn") #Normalize with robust spline regression
lumi.patient.qual <- lumiQ(lumi.patient.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05

qcsum <- lumi.patient.qual@QC$sampleSummary %>% t %>% data.frame
colnames(qcsum) %<>% str_replace("\\.0\\.01\\.", "")
qcsum$Sample.Name <- rownames(qcsum)
qcsum$RIN <- lumi.patient.qual$RIN
qcsum$Sample.Num <- lumi.patient.qual$Sample.Num
#plot(lumi.patient.qual, what = 'sampleRelation')

lumi.patient.cutoff <- detectionCall(lumi.patient.qual) #Get the count of probes which passed the detection threshold per sample
lumi.patient.expr <- lumi.patient.qual[which(lumi.patient.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.patient <- getSYMBOL(rownames(lumi.patient.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.patient.annot <- lumi.patient.expr[!symbols.lumi.patient,] #Drop any probe which is not annotated
saveRDS.gz(lumi.patient.annot, file = "./save/lumi.patient.annot.rda")

#Use ComBat for batch effect correction
model.sex <- model.matrix( ~ 0 + factor(lumi.patient.annot$Sex) )
colnames(model.sex) <- c("Female", "Male")
model.sex.reduce <- model.sex[,-1]

model.combat <- cbind(Male = model.sex.reduce, Age = as.numeric(lumi.patient.annot$Draw.Age), RIN = lumi.patient.annot$RIN, GAA = lumi.patient.annot$GAA1, Onset = as.numeric(lumi.patient.annot$Onset))

expr.combat <- ComBat(dat = exprs(lumi.patient.annot), batch = factor(lumi.patient.annot$Batch), mod = model.combat)
lumi.combat <- lumi.patient.annot
exprs(lumi.combat) <- expr.combat
saveRDS.gz(lumi.combat, file = "./save/lumi.combat.rda")

#Run PEER analysis and correlate to known covariates
gen.peer(8, exprs(lumi.combat), TRUE, model.combat)
PEER.weights.plot <- read_csv("./weight_8.csv") %>% select(-(X1:X6))
PEER.weights.sums <- colSums(abs(PEER.weights.plot)) %>% data.frame
PEER.weights.sums$Factor <- 1:nrow(PEER.weights.sums)
colnames(PEER.weights.sums)[1] <- "Weight"

p <- ggplot(PEER.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
p <- p + theme_bw() + xlab("Factor") + ylab("Weigtht")
CairoPDF("./PEER_weights", height = 4, width = 6)
print(p)
dev.off()

model.PEER_covariate <- read_csv("./factor_8.csv") %>% select(-(X1:X6))
model.cov <- cbind(Male = model.sex.reduce, Age = lumi.combat$Draw.Age, RIN = lumi.combat$RIN)
model.full.cov <- cbind(model.cov, model.PEER_covariate)
cleaned.expr <- removeBatchEffect(exprs(lumi.combat), covariates = model.full.cov)
lumi.cleaned <- lumi.combat
exprs(lumi.cleaned) <- cleaned.expr
saveRDS.gz(lumi.cleaned, file = "./save/lumi.cleaned.rda")

batch.colors <- c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green")
gen.boxplot("baseline_intensity_corrected.jpg", lumi.cleaned, batch.colors, "Covariate-corrected intensity", "Intensity")

fdata <- fData(lumi.cleaned)
lumi.cleaned <- lumi.cleaned[!is.na(fdata$SYMBOL),]
fdata <- fData(lumi.cleaned)
pdata <- pData(lumi.cleaned)
expr.collapse <- collapseRows(exprs(lumi.cleaned), factor(fdata$SYMBOL), rownames(lumi.cleaned), method = "function", methodFunction = colMeans)$datETcollapsed
saveRDS.gz(expr.collapse, "./save/expr.collapse.rda")

#Plot GSTM3 vs GAA
gstm.cor <- filter(cor.df, Symbol == "GSTM3")
gstm.df <- data.frame(Expression = expr.collapse["GSTM3",], GAA = pdata$GAA1)
p <- ggplot(gstm.df, aes(x = GAA, y = Expression)) + geom_point() + geom_smooth(method = lm) + xlab("GAA1 (# of GAA repeats)") + ylab("VST normalized Expression")
CairoPDF("gstm_gaa", width = 7, height = 5)
print(p)
dev.off()

#Correlations
gaa.cor <- apply(expr.collapse, 1, cor, lumi.cleaned$GAA1)
gaa.cor.pval <- corPvalueStudent(gaa.cor, length(lumi.cleaned$GAA1)) %>% p.adjust("fdr")
cor.df <- cbind(gaa.cor, gaa.cor.pval) %>% data.frame
colnames(cor.df) <- c("Correlation", "P.value")
#cor.df$nuId <- featureNames(lumi.cleaned)
cor.df$Symbol <- rownames(expr.collapse)
#cor.df$Definition <- featureData(lumi.cleaned)@data$DEFINITION
cor.df %<>% arrange(P.value)
cor.df.sg <- filter(cor.df, P.value < 0.05)
write.xlsx(cor.df.sg, "significant_gaa.xlsx")

#cor.probes <- c("ILMN_1726589", "ILMN_1659688", "ILMN_1703524", "ILMN_1765796", "ILMN_2184373", "ILMN_1745172", "ILMN_1747344") %>% paste(collapse = "|")
#intensities1.cor <- data.frame(Probe_Id = rownames(intensities1.gaa), intensities1.gaa) %>% filter(grepl(cor.probes, Probe_Id))
#annot.symbol <- select(annot.reduce, Probe_Id, Symbol)
#intensities1.joined <- join(intensities1.cor, annot.symbol)
#apply(intensities1.joined, 1, plot.cor, targets1.gaa$GAA1)

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.data.PEER <- t(expr.collapse)
saveRDS.gz(expr.data.PEER, file = "./save/expr.data.PEER.rda")
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

softPower <- 10
adjacency.PEER <- adjacency(expr.data.PEER, power = softPower, type = "signed")
saveRDS.gz(adjacency.PEER, file = "./save/adjacency.PEER.rda")

TOM.PEER <- TOMsimilarity(adjacency.PEER)
dissimilarity.TOM <- 1 - TOM.PEER
saveRDS.gz(TOM.PEER, file = "./save/tom.PEER.rda")
saveRDS.gz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 50

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, cutHeight = 0.995, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.data.PEER, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - cor(ME.genes)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.05
merge.all <- mergeCloseModules(expr.data.PEER, dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3) #PC analysis may be failing because of Intel MKL Lapack routine bug.  Test with openBLAS in R compiled with gcc.
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

all.degrees <- intramodularConnectivity(adjacency.PEER, module.colors)
colnames(fdata) %<>% tolower %>% capitalize
gene.info <- data.frame(Symbol = rownames(all.degrees), module.color = module.colors, all.degrees)

write_csv(data.frame(table(module.colors)), path = "./final_eigengenes.csv") 
gene.info$kscaled <- by(gene.info, gene.info$module.color, select, kWithin) %>% llply(function(x) { x / max (x) }) %>% reduce(c)
saveRDS.gz(gene.info, file = "./save/gene.info.rda")

gene.module.membership <- as.data.frame(cor(expr.data.PEER, ME.genes, use = "p"))
module.membership.pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.data.PEER)))
names(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")

module.membership <- cbind(select(gene.info, Symbol:module.color), gene.module.membership, module.membership.pvalue)
write_csv(module.membership, "module_membership.csv")

colnames(gene.module.membership) %<>% str_replace("MM.", "")
colnames(module.membership.pvalue) %<>% str_replace("MM.pvalue.", "")
gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

color.values <- unique(module.colors)
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

sample.ids <- factor(rownames(expr.data.PEER), levels = rownames(expr.data.PEER))
colnames(ME.genes) %<>% str_replace("ME", "")
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.data.PEER), module.colors)
cluster <- makeForkCluster(8)
split(expr.data.plot, expr.data.plot$module.colors) %>% parLapply(cl = cluster, gen.heatmap, ME.genes.plot)

modules.out <- select(gene.info, Symbol, module.color)
write.xlsx(modules.out, "modules_out.xlsx")

source("../../code/GO/enrichr.R")

enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 
#enrichr.submit("blue", modules.out, enrichr.terms, FALSE)
color.names <- unique(module.colors) %>% sort
grey.location <- str_detect(color.names, "grey")
color.names.reduce <- color.names[!grey.location]
l_ply(color.names.reduce, enrichr.submit, modules.out, enrichr.terms, FALSE)

modules.out.reduce <- filter(modules.out, module.color != "grey")
modules.out.reduce$module.color %<>% droplevels

split(modules.out.reduce, modules.out.reduce$module.color) %>% map(submit.stringdb)

submit.stringdb <- function(module.subset)
{
    get.stringdb(module.subset, unique(module.subset$module.color), "./stringdb")
}

targets.final.gaa <- pData(lumi.cleaned)
#targets.final.gaa$Sample.Name %<>% str_replace(" ", "")
rownames(ME.genes) <- rownames(expr.data.PEER)
PEER.factors <- read_csv("./factor_8.csv") %>% select(-(X1:X6)) 
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

targets.gaa <- select(targets.final.gaa, Sample.Name, GAA1) 
cor.gaa <- gen.cor(ME.genes, targets.gaa)

targets.final.gaa$Onset %<>% as.numeric
targets.onset <- select(targets.final.gaa, Sample.Name, Onset) 
cor.onset <- gen.cor(ME.genes, targets.onset)

targets.rin <- select(targets.final.gaa, Sample.Name, RIN)
cor.rin <- gen.cor(ME.genes, targets.rin)

PEER.factors.df <- mutate(PEER.factors, Sample.Name = targets.final.gaa$Sample.Name) 
cor.PEER <- gen.cor(ME.genes, PEER.factors.df) %>% data.frame

module.traits.all <- cbind(cor.age, cor.sex, cor.gaa, cor.onset, cor.rin) %>% data.frame
module.traits.pval <- select(module.traits.all, contains("p.value")) %>% as.matrix %>% apply(2, p.adjust, "fdr")
module.traits.cor <- select(module.traits.all, -contains("p.value")) %>% as.matrix

cor.PEER.cor <- select(cor.PEER, -contains("p.value")) %>% as.matrix
cor.PEER.pval <- select(cor.PEER, contains("p.value")) %>% as.matrix %>% apply(2, p.adjust, "fdr")

module.trait.out <- data.frame(Module = rownames(module.traits.cor), module.traits.cor, module.traits.pval)
cor.PEER.out <- data.frame(Module = rownames(cor.PEER), cor.PEER.cor, cor.PEER.pval)

write_csv(module.trait.out, "module_trait_cor.csv")
write_csv(cor.PEER.out, "cor_PEER_cor.csv")

text.matrix.traits <- paste(signif(module.traits.cor, 2), '\n(', signif(module.traits.pval, 1), ')', sep = '')
text.matrix.PEERcor <- paste(signif(cor.PEER.cor, 2), '\n(', signif(cor.PEER.pval, 1), ')', sep = '')
dim(text.matrix.traits) = dim(module.traits.cor)
dim(text.matrix.PEERcor) = dim(cor.PEER.cor)

gen.text.heatmap(module.traits.cor, text.matrix.traits, colnames(module.traits.cor), colnames(ME.genes), "", "module-trait relationships")
gen.text.heatmap(cor.PEER.cor, text.matrix.PEERcor, colnames(cor.PEER.cor), colnames(ME.genes), "", "module-PEER factor relationships")

plot.eigencor(module.traits.pval, "GAA1", lumi.cleaned$GAA1)
plot.eigencor(module.traits.pval, "Onset", lumi.cleaned$Onset)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort


