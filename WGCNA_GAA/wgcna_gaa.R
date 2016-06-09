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
library(irr)
library(tools)

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
library(heatmap.plus)

#Reading and writing tables
library(readr)
library(openxlsx)
library(parallel)

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

plot.eigencor <- function(module.color, col.name, ME.genes, pheno.vector)
{
    me.genes.subset <- select(ME.genes, one_of(as.character(module.color)))
    me.genes.subset[[col.name]] <- pheno.vector
    me.subset.melt <- melt(me.genes.subset, id.vars = col.name) 
    colnames(me.subset.melt)[2] <- "Module"

    me.subset.melt$Module %<>% as.character
    p <- ggplot(me.subset.melt, aes_string(x = col.name, y = "value", col = "Module")) + geom_point(position = "jitter")
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene") + stat_smooth(method = "lm", se = FALSE)
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(pheno.vector)))
    p <- p + scale_color_manual(values = sort(unique(me.subset.melt$Module)))
    p <- p + stat_smooth(method = "lm", se = TRUE)
    p <- p + theme(legend.position = "none")

    filename <- paste(col.name, module.color, "eigengenes_05", sep = "_")
    CairoPDF(filename, height = 5, width = 6)
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
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    map(names(enrichr.data), enrichr.wkbk, enrichr.data, index)
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

lumi.import <- readRDS.gz(file = "../baseline_lumi/save/lumi.baseline.rda")
remove.all <- readRDS.gz(file = "../baseline_lumi/save/remove.all.rda")

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

lumi.patient.cutoff <- detectionCall(lumi.patient.qual) #Get the count of probes which passed the detection threshold per sample
lumi.patient.expr <- lumi.patient.qual[which(lumi.patient.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.patient <- getSYMBOL(rownames(lumi.patient.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.patient.annot <- lumi.patient.expr[!symbols.lumi.patient,] #Drop any probe which is not annotated
saveRDS.gz(lumi.patient.annot, file = "./save/lumi.patient.annot.rda")

#Use ComBat for batch effect correction
model.sex <- model.matrix( ~ 0 + factor(lumi.patient.annot$Sex) )
colnames(model.sex) <- c("Female", "Male")
model.sex.reduce <- model.sex[,-1]

model.combat <- cbind(Male = model.sex.reduce, Age = as.numeric(lumi.patient.annot$Draw.Age), RIN = lumi.patient.annot$RIN, GAA = lumi.patient.annot$GAA1) %>% data.frame
age.sex <- lm(model.combat$Age ~ model.combat$Male) %>% anova
age.gaa <- lm(model.combat$Age ~ model.combat$GAA) %>% anova
sex.gaa <- lm(model.combat$Male ~ model.combat$GAA) %>% anova

expr.combat <- ComBat(dat = exprs(lumi.patient.annot), batch = factor(lumi.patient.annot$Site), mod = model.combat)
lumi.combat <- lumi.patient.annot
exprs(lumi.combat) <- expr.combat

expr.combat <- ComBat(dat = exprs(lumi.combat), batch = factor(lumi.combat$Batch), mod = model.combat)
exprs(lumi.combat) <- expr.combat
saveRDS.gz(lumi.combat, file = "./save/lumi.combat.rda")

model.cov <- cbind(Male = model.sex.reduce, Age = lumi.combat$Draw.Age, RIN = lumi.combat$RIN)
cleaned.expr <- removeBatchEffect(exprs(lumi.combat), covariates = model.cov)
lumi.cleaned <- lumi.combat
exprs(lumi.cleaned) <- cleaned.expr
saveRDS.gz(lumi.cleaned, file = "./save/lumi.cleaned.rda")

#batch.colors <- c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green")
#gen.boxplot("baseline_intensity_corrected.jpg", lumi.cleaned, batch.colors, "Covariate-corrected intensity", "Intensity")

pdata <- pData(lumi.combat)
expr.collapse <- collapseRows(exprs(lumi.cleaned), getSYMBOL(featureNames(lumi.cleaned), 'lumiHumanAll.db'), rownames(lumi.cleaned))$datETcollapsed
export.lumi <- ExpressionSet(assayData = expr.collapse, phenoData = phenoData(lumi.combat))
saveRDS.gz(expr.collapse, "./save/expr.collapse.rda")
saveRDS.gz(export.lumi, "./save/export.lumi.rda")

DGKD.cor <- filter(cor.df, Symbol == "DGKD")
DGKD.df <- data.frame(Expression = expr.collapse["DGKD",], GAA = pdata$GAA1)
p <- ggplot(DGKD.df, aes(x = GAA, y = Expression)) + geom_point() + geom_smooth(method = lm) 
p <- p + xlab("GAA1 (# of GAA repeats)") + ylab("VST normalized Expression") 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", size = 2))
CairoPDF("DGKD_gaa", width = 7, height = 5)
print(p)
dev.off()

CDC42EP1.cor <- filter(cor.df, Symbol == "CDC42EP1")
CDC42EP1.df <- data.frame(Expression = expr.collapse["CDC42EP1",], GAA = pdata$GAA1)
p <- ggplot(CDC42EP1.df, aes(x = GAA, y = Expression)) + geom_point() + geom_smooth(method = lm) 
p <- p + xlab("GAA1 (# of GAA repeats)") + ylab("VST normalized Expression") 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", size = 2))
CairoPDF("CDC42EP1_gaa", width = 7, height = 5)
print(p)
dev.off()

RANBP2.cor <- filter(cor.df, Symbol == "RANBP2")
RANBP2.df <- data.frame(Expression = expr.collapse["RANBP2",], GAA = pdata$GAA1)
p <- ggplot(RANBP2.df, aes(x = GAA, y = Expression)) + geom_point() + geom_smooth(method = lm) 
p <- p + xlab("GAA1 (# of GAA repeats)") + ylab("VST normalized Expression") 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", size = 2))
CairoPDF("RANBP2_gaa", width = 7, height = 5)
print(p)
dev.off()

UBE2D3.cor <- filter(cor.df, Symbol == "UBE2D3")
UBE2D3.df <- data.frame(Expression = expr.collapse["UBE2D3",], GAA = pdata$GAA1)
p <- ggplot(UBE2D3.df, aes(x = GAA, y = Expression)) + geom_point() + geom_smooth(method = lm) 
p <- p + xlab("GAA1 (# of GAA repeats)") + ylab("VST normalized Expression") 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", size = 2)) 
CairoPDF("UBE2D3_gaa", width = 7, height = 5)
print(p)
dev.off()

#Correlations
gaa.cor <- apply(expr.collapse, 1, cor, lumi.cleaned$GAA1)
gaa.cor.pval <- corPvalueStudent(gaa.cor, length(lumi.cleaned$GAA1)) 
gaa.cor.adjpval <- p.adjust(gaa.cor.pval, "fdr")
cor.df <- cbind(gaa.cor, gaa.cor.pval, gaa.cor.adjpval) %>% data.frame
colnames(cor.df) <- c("Correlation", "P.value", "Adj.P.value")
cor.df$Symbol <- rownames(expr.collapse)
cor.df %<>% arrange(P.value)
cor.df.sg <- filter(cor.df, Adj.P.value < 0.05)
write.xlsx(cor.df.sg, "significant_gaa.xlsx")
write.xlsx(cor.df, "all.gaa.xlsx")

cor.enrichr <- map(enrichr.terms, get.enrichrdata, cor.df.sg, FALSE)
names(cor.enrichr) <- enrichr.terms
map(names(cor.enrichr), enrichr.wkbk, cor.enrichr, "gaa_cor")

get.updown <- function(filter.vector, enrichr.df)
{
    grep.vector <- str_replace_all(filter.vector, ",", "|")
    enrichr.filter <- filter(enrichr.df, grepl(grep.vector, Symbol))
    num.up <- which(enrichr.filter$Correlation > 0) %>% length
    num.down <- which(enrichr.filter$Correlation < 0) %>% length
    return(c("Up" = num.up, "Down" = num.down))
}

gaa.gobiol <- read.xlsx("./enrichr/gaa_cor/gaa_cor_GO_Biological_Process.xlsx") %>% select(GO.Term, P.value, Genes) %>% slice(c(1, 8, 9, 18))
gaa.gobiol$Database <- "GO Biological Process"
gaa.gomole <- read.xlsx("./enrichr/gaa_cor/gaa_cor_GO_Molecular_Function.xlsx") %>% select(GO.Term, P.value, Genes) %>% slice(c(1, 8))
gaa.gomole$Database <- "GO Molecular Process"
gaa.reactome <- read.xlsx("./enrichr/gaa_cor/gaa_cor_Reactome_2015.xlsx") %>% select(GO.Term, P.value, Genes) %>% slice(c(23, 25, 34))
gaa.reactome$Database <- "Reactome"
gaa.enrichr <- rbind(gaa.gobiol, gaa.gomole, gaa.reactome)

gaa.enrichr$Gene.Count <- map(gaa.enrichr$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
gaa.enrichr$Log.pvalue <- -(log10(gaa.enrichr$P.value))

gaa.updown <- map(gaa.enrichr$Genes, get.updown, cor.df.sg) %>% reduce(rbind)
#colnames(gaa.updown) <- c("Up", "Down")
gaa.enrichr <- cbind(gaa.enrichr, gaa.updown)
gaa.enrichr$Log.Up <- gaa.enrichr$Log.pvalue * gaa.enrichr$Up / gaa.enrichr$Gene.Count
gaa.enrichr$Log.Down <- gaa.enrichr$Log.pvalue * gaa.enrichr$Down / gaa.enrichr$Gene.Count
gaa.enrichr$GO.Term %<>% str_replace_all("\\ \\(.*$", "") %>% tolower
gaa.enrichr$Format.Name <- paste(gaa.enrichr$Database, ": ", gaa.enrichr$GO.Term, " (", gaa.enrichr$Gene.Count, ")", sep = "")
gaa.enrichr.plot <- select(gaa.enrichr, Format.Name, Log.Up, Log.Down) %>% melt(id.vars = "Format.Name") 

p <- ggplot(gaa.enrichr.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = c(gaa.enrichr$Format.Name, rep("", length(gaa.enrichr$Format.Name))), hjust = "left", aes(y = 0.1))
p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "FALSE",  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
CairoPDF("gaa.enrichr", height = 5, width = 8)
print(p)
dev.off()

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.data <- t(expr.collapse)
saveRDS.gz(expr.data, file = "./save/expr.data.rda")
sft <- pickSoftThreshold(expr.data, powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.df <- sft$fitIndices
saveRDS.gz(sft, file = "./save/sft.peer.rda")

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
adjacency.expr <- adjacency(expr.data, power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
saveRDS.gz(adjacency.expr, file = "./save/adjacency.expr.rda")

TOM <- TOMsimilarity(adjacency.expr)
dissimilarity.TOM <- 1 - TOM
saveRDS.gz(TOM, file = "./save/tom.rda")
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
ME.list <- moduleEigengenes(expr.data, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.20
merge.all <- mergeCloseModules(expr.data, dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) #PC analysis may be failing because of Intel MKL Lapack routine bug.  Test with openBLAS in R compiled with gcc.
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

CairoPDF("eigengenes", height = 6, width = 10)
par(cex = 0.7)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.adjacency = 0.3, cex.preservation = 0.3, plotPreservation = "standard")
dev.off()

all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
colnames(fdata) %<>% tolower %>% capitalize
gene.info <- data.frame(Symbol = rownames(all.degrees), module.color = module.colors, all.degrees)

write_csv(data.frame(table(module.colors)), path = "./final_eigengenes.csv") 
gene.info$kscaled <- by(gene.info, gene.info$module.color, select, kWithin) %>% llply(function(x) { x / max (x) }) %>% reduce(c)
saveRDS.gz(gene.info, file = "./save/gene.info.rda")

gene.module.membership <- as.data.frame(bicor(expr.data, ME.genes, maxPOutliers = 0.05))
module.membership.pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.data)))
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

sample.ids <- factor(rownames(expr.data), levels = rownames(expr.data))
colnames(ME.genes) %<>% str_replace("ME", "")
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.data), module.colors)
cluster <- makeForkCluster(8)
split(expr.data.plot, expr.data.plot$module.colors) %>% map(gen.heatmap, ME.genes.plot)

modules.out <- select(gene.info, Symbol, module.color)
write.xlsx(modules.out, "modules_out.xlsx")

source("../../code/GO/enrichr.R")
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "HumanCyc_2016", "NCI-Nature_2016", "Panther_2016")
enrichr.submit("pink", modules.out, enrichr.terms, FALSE)
color.names <- unique(module.colors) %>% sort
l_ply(color.names, enrichr.submit, modules.out, enrichr.terms, FALSE)

#modules.out.reduce <- filter(modules.out, module.color != "grey")
#modules.out.reduce$module.color %<>% droplevels

#split(modules.out.reduce, modules.out.reduce$module.color) %>% map(submit.stringdb)

#submit.stringdb <- function(module.subset)
#{
    #get.stringdb(module.subset, unique(module.subset$module.color), "./stringdb")
#}

connectivity.reduce <- filter(eigengene.connectivity, module.color != "grey" & module.comparison != "grey" & module.color == module.comparison)

targets.final.gaa <- pData(lumi.cleaned)
rownames(ME.genes) <- rownames(expr.data)

source("../common_functions.R")
targets.gaa <- select(targets.final.gaa, Sample.Name, GAA1) 
cor.gaa <- gen.cor(ME.genes, targets.gaa)

module.traits.all <- data.frame(cor.gaa)
module.traits.pval <- select(module.traits.all, contains("p.value")) %>% as.matrix %>% apply(2, p.adjust, "fdr")
module.traits.cor <- select(module.traits.all, -contains("p.value")) %>% as.matrix

module.trait.out <- data.frame(Module = rownames(module.traits.cor), module.traits.cor, module.traits.pval)

write_csv(module.trait.out, "module_trait_cor.csv")

text.matrix.traits <- paste(signif(module.traits.cor, 2), '\n(', signif(module.traits.pval, 1), ')', sep = '')
dim(text.matrix.traits) = dim(module.traits.cor)

gen.text.heatmap(module.traits.cor, text.matrix.traits, colnames(module.traits.cor), colnames(ME.genes), "", "module-trait relationships")

ME.genes.plot <- ME.genes
colnames(ME.genes.plot) %<>% str_replace_all("ME", "")
plot.eigencor("yellow", "GAA1", ME.genes.plot, lumi.cleaned$GAA1)
plot.eigencor("brown", "GAA1", ME.genes.plot, lumi.cleaned$GAA1)
plot.eigencor("greenyellow", "GAA1", ME.genes.plot, lumi.cleaned$GAA1)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

get.kappa <- function(term.current, all.terms)
{
    map(all.terms, cbind, term.current) %>% map(kappa2) %>% map_dbl(getElement, "value")
}

get.kappa.cluster <- function(enrichr.output, gene.names, filename)
{
    num.genes <- length(gene.names)
    enrichr.list <- map(enrichr.output$Genes, str_split, ",") %>% map(getElement, 1) 
    enrichr.match <- map(enrichr.list, is.element, el = toupper(gene.names)) %>% reduce(rbind) %>% t
    rownames(enrichr.match) <- toupper(gene.names)
    colnames(enrichr.match) <- enrichr.output$Term
    enrichr.match.df <- data.frame(enrichr.match)

    enrichr.kappa <- map(enrichr.match.df, get.kappa, enrichr.match.df) %>% reduce(rbind)
    enrichr.kappa[enrichr.kappa < 0] <- 0

    rownames(enrichr.kappa) <- colnames(enrichr.kappa) <- enrichr.output$Term

    CairoPDF(str_c(filename, "heatmap", sep = "."), width = 30, height = 30)
    heatmap.plus(enrichr.kappa, col = heat.colors(40), symm = TRUE, margins = c(20,20))
    dev.off()

    kappa.dist <- dist(enrichr.kappa, method = "manhattan")
    kappa.clust <- hclust(kappa.dist, method = "average")

    CairoPDF(str_c(filename, "clust", sep = "."), height = 30, width = 30)
    plot(kappa.clust)
    dev.off()

    kappa.modules <- cutreeDynamic(kappa.clust, minClusterSize = 2, method = "tree")
    kappa.modules.TOM <- cutreeDynamic(kappa.clust, distM = TOMdist(enrichr.kappa), minClusterSize = 2, method = "hybrid")
    kappa.modules.df <- data.frame(Term = rownames(enrichr.kappa), Module = kappa.modules, Module.TOM = kappa.modules.TOM)

    enrichr.output$Module <- kappa.modules
    enrichr.output$Module.TOM <- kappa.modules.TOM
    enrichr.output %<>% select(Index:Combined.Score, Module:Module.TOM, Genes)
    
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = enrichr.output, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(enrichr.output)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    
    saveWorkbook(wb, str_c(filename, "table.xlsx", sep = "."), overwrite = TRUE) 
}

yellow.only <- filter(modules.out, module.color == "yellow")
yellow.gobiol.file <- "./enrichr/yellow/yellow_GO_Biological_Process_2015.xlsx"
yellow.gobiol <- read.xlsx(yellow.gobiol.file) 
yellow.gobiol$Num.Genes <- map(yellow.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
yellow.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
yellow.gobiol$Database <- "GO Biological Process"
get.kappa.cluster(yellow.gobiol, yellow.only$Symbol, file_path_sans_ext(yellow.gobiol.file))
yellow.gobiol.final <- slice(yellow.gobiol, c(52, 63))

yellow.gomole.file <- "./enrichr/yellow/yellow_GO_Molecular_Function_2015.xlsx"
yellow.gomole <- read.xlsx(yellow.gomole.file) 
yellow.gomole$Num.Genes <- map(yellow.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
yellow.gomole %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
yellow.gomole$Database <- "GO Molecular Function"
get.kappa.cluster(yellow.gomole, yellow.only$Symbol, file_path_sans_ext(yellow.gomole.file))
yellow.gomole.final <- slice(yellow.gomole, c(12))

yellow.reactome.file <- "./enrichr/yellow/yellow_Reactome_2016.xlsx"
yellow.reactome <- read.xlsx(yellow.reactome.file) 
yellow.reactome$Num.Genes <- map(yellow.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
yellow.reactome %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
yellow.reactome$Database <- "Reactome"
get.kappa.cluster(yellow.reactome, yellow.only$Symbol, file_path_sans_ext(yellow.reactome.file))
yellow.reactome.final <- slice(yellow.reactome, 29)

yellow.kegg <- read.xlsx("./enrichr/yellow/yellow_KEGG_2016.xlsx") %>% slice(2)
yellow.kegg$Num.Genes <- map(yellow.kegg$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
yellow.kegg$Database <- "KEGG"

yellow.enrichr <- rbind(yellow.gobiol.final, yellow.gomole.final, yellow.reactome.final, yellow.kegg)
gen.enrichrplot(yellow.enrichr, "yellow.enrichr")

#get kegg

brown.only <- filter(modules.out, module.color == "brown")
brown.gobiol.file <- "./enrichr/brown/brown_GO_Biological_Process_2015.xlsx"
brown.gobiol <- read.xlsx(brown.gobiol.file) 
brown.gobiol$Num.Genes <- map(brown.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
brown.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
brown.gobiol$Database <- "GO Biological Process"
get.kappa.cluster(brown.gobiol, brown.only$Symbol, file_path_sans_ext(brown.gobiol.file))
brown.gobiol.final <- slice(brown.gobiol, c(31, 43))

brown.gomole.file <- "./enrichr/brown/brown_GO_Molecular_Function_2015.xlsx"
brown.gomole <- read.xlsx(brown.gomole.file) 
brown.gomole$Num.Genes <- map(brown.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
brown.gomole %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
brown.gomole$Database <- "GO Molecular Function"
get.kappa.cluster(brown.gomole, brown.only$Symbol, file_path_sans_ext(brown.gomole.file))
brown.gomole.final <- slice(brown.gomole, c(18, 14))

brown.reactome.file <- "./enrichr/brown/brown_Reactome_2016.xlsx"
brown.reactome <- read.xlsx(brown.reactome.file) 
brown.reactome$Num.Genes <- map(brown.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
brown.reactome %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
brown.reactome$Database <- "Reactome"
get.kappa.cluster(brown.reactome, brown.only$Symbol, file_path_sans_ext(brown.reactome.file))
brown.reactome.final <- slice(brown.reactome, c(11))

#get kegg

brown.enrichr <- rbind(brown.gobiol.final, brown.gomole.final, brown.reactome.final)
gen.enrichrplot(brown.enrichr, "brown.enrichr")

greenyellow.only <- filter(modules.out, module.color == "greenyellow")
greenyellow.gobiol.file <- "./enrichr/greenyellow/greenyellow_GO_Biological_Process_2015.xlsx"
greenyellow.gobiol <- read.xlsx(greenyellow.gobiol.file) 
greenyellow.gobiol$Num.Genes <- map(greenyellow.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
greenyellow.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
greenyellow.gobiol$Database <- "GO Biological Process"
get.kappa.cluster(greenyellow.gobiol, greenyellow.only$Symbol, file_path_sans_ext(greenyellow.gobiol.file))
greenyellow.gobiol.final <- slice(greenyellow.gobiol, c(5, 34, 29, 39))

greenyellow.gomole.file <- "./enrichr/greenyellow/greenyellow_GO_Molecular_Function_2015.xlsx"
greenyellow.gomole <- read.xlsx(greenyellow.gomole.file) 
greenyellow.gomole$Num.Genes <- map(greenyellow.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
greenyellow.gomole %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
greenyellow.gomole$Database <- "GO Molecular Function"
get.kappa.cluster(greenyellow.gomole, greenyellow.only$Symbol, file_path_sans_ext(greenyellow.gomole.file))
greenyellow.molec.final <- slice(greenyellow.gomole, c(1, 3))

greenyellow.reactome.file <- "./enrichr/greenyellow/greenyellow_Reactome_2016.xlsx"
greenyellow.reactome <- read.xlsx(greenyellow.reactome.file) 
greenyellow.reactome$Num.Genes <- map(greenyellow.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
greenyellow.reactome %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
greenyellow.reactome$Database <- "Reactome"
get.kappa.cluster(greenyellow.reactome, greenyellow.only$Symbol, file_path_sans_ext(greenyellow.reactome.file))

greenyellow.final <- rbind(greenyellow.gobiol.final, greenyellow.molec.final)
gen.enrichrplot(greenyellow.final, "greenyellow.enrichr")

#PPI
biogrid.ppi <- read_tsv("~/Downloads/BIOGRID-ORGANISM-3.4.136.tab2/BIOGRID-ORGANISM-Homo_sapiens-3.4.136.tab2.txt") %>% data.frame
biogrid.ppi.reduce <- select(biogrid.ppi, contains("Official"))

inweb.ppi <- read_tsv("../../Dementia Project/mapt/InWeb3_HC_NonRed.txt", col_names = FALSE)
colnames(inweb.ppi) <- c("Interactor.A", "Interactor.B")
iw.hugo <- read_tsv("../../Dementia Project/mapt/IWtoHugo.txt", col_names = FALSE)
colnames(iw.hugo) <- c("IW.ID", "HUGO")

yellow.kegg.als <- str_split(yellow.kegg$Genes, ",")[[1]]
yellow.reactome.oxid <- str_split(yellow.reactome$Genes, ",")[[1]]
yellow.gobiol.tetra <- str_split(yellow.gobiol.final[2,]$Genes, ",")[[1]]

brown.gobiol.acet <- str_split(brown.gobiol.final[1,]$Genes, ",")[[1]]
brown.reactome.tca <- str_split(brown.reactome.final$Genes, ",")[[1]]
brown.gomole.iron <- str_split(brown.gomole.final[1,]$Genes, ",")[[1]]

yka.ppi <- get.ppi(yellow.kegg.als) #Maybe don't prune
yro.ppi <- get.ppi(yellow.reactome.oxid) 
ygt.ppi <- get.ppi(yellow.gobiol.tetra) #Don't prune

bga.ppi <- get.ppi(brown.gobiol.acet) #Don't prune
brt.ppi <- get.ppi(brown.reactome.tca) #Hueg
bgi.ppi <- get.ppi(brown.gomole.iron) #Don't prune?

adjacency.expr <- readRDS.gz("./save/adjacency.expr.rda")

get.ppi <- function(gene.list)
{
    ppi.id <- filter(iw.hugo, is.element(HUGO, gene.list))$IW.ID
    ppi.inweb.interactors <- filter(inweb.ppi, is.element(Interactor.A, ppi.id)) %>% filter(is.element(Interactor.B, ppi.id))
    ppi.inweb.symbols <- merge(ppi.inweb.interactors, iw.hugo, by.x = "Interactor.A", by.y = "IW.ID") %>% merge(iw.hugo, by.x = "Interactor.B", by.y = "IW.ID")
    ppi.inweb.final <- select(ppi.inweb.symbols, contains("HUGO"))
    colnames(ppi.inweb.final) <- c("Symbol.A", "Symbol.B")

    ppi.biogrid <- filter(biogrid.ppi.reduce, Official.Symbol.Interactor.A %in% gene.list) %>% filter(Official.Symbol.Interactor.B %in% gene.list)
    colnames(ppi.biogrid) <- c("Symbol.A", "Symbol.B")

    ppi.combined <- rbind(ppi.inweb.final, ppi.biogrid)
    ppi.combined$unique <- str_c(ppi.combined$Symbol.A, ppi.combined$Symbol.B, sep = ".")
    ppi.unique <- filter(ppi.combined, !duplicated(unique)) %>% select(-unique)
    ppi.self <- apply(ppi.unique, 1, reduce, identical)
    ppi.unique.final <- filter(ppi.unique, !ppi.self)
    return(ppi.unique.final)
}

plot.ppi(adjacency.expr, yellow.kegg.als, yka.ppi, "yka_igraph")
yro.igraph <- plot.ppi(adjacency.expr, yellow.reactome.oxid, yro.ppi, "yro_igraph", prune = TRUE, clust.keep = 2)
yro.incident <- map(1:vcount(yro.igraph), incident, graph = yro.igraph) %>% map_dbl(length)
names(yro.incident) <- V(yro.igraph)$name
plot.ppi(adjacency.expr, yellow.gobiol.tetra, ygt.ppi, "ygt_igraph")

plot.ppi(adjacency.expr, brown.gobiol.acet, bga.ppi, "bga_igraph")
brt.igraph <- plot.ppi(adjacency.expr, brown.reactome.tca, brt.ppi, "brt_igraph", TRUE)
brt.incident <- map(1:vcount(brt.igraph), incident, graph = brt.igraph) %>% map_dbl(length)
names(brt.incident) <- V(brt.igraph)$name
plot.ppi(adjacency.expr, brown.gomole.iron, bgi.ppi, "bgi_igraph")

#PPI Plot
plot.ppi <- function(adjacency.expr, gene.list, ppi.edges, filename, prune = FALSE, clust.keep = 1, plot.width = 7, plot.height = 7)
{
    expr.adjacency <- adjacency.expr[gene.list,gene.list]
    ppi.adjacency <- matrix(0, ncol = ncol(expr.adjacency), nrow = nrow(expr.adjacency), dimnames = list(rownames(expr.adjacency), colnames(expr.adjacency)))
    ppi.adjacency[as.matrix(ppi.edges)] <- 1
    final.adjacency <- expr.adjacency * ppi.adjacency
    final.mins <- apply(final.adjacency, 2, min)
    final.maxs <- apply(final.adjacency, 2, max)
    final.scaled <- scale(final.adjacency, center = final.mins, scale = final.maxs - final.mins)
    final.scaled[is.nan(final.scaled)] <- 0

    final.igraph <- graph_from_adjacency_matrix(expr.adjacency, mode = "undirected", weighted = TRUE, diag = FALSE)
    ppi.igraph <- graph_from_adjacency_matrix(final.scaled, mode = "undirected", weighted = TRUE, diag = FALSE)

    if (prune == TRUE)
    {
        #num.edges <- map(1:vcount(ppi.igraph), incident, graph = ppi.igraph) %>% map_dbl(length)

        final.igraph <- delete.vertices(final.igraph, which(clusters(ppi.igraph)$membership != clust.keep))
        ppi.igraph <- delete.vertices(ppi.igraph, which(clusters(ppi.igraph)$membership != clust.keep))
    }

    #ppi.colors <- rainbow(length(unique(clusters(ppi.igraph)$membership)))
    #final.colors <- ppi.colors[clusters(ppi.igraph)$membership]

    final.edges <- attr(E(final.igraph), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_")
    final.df <- data.frame(Edges = final.edges, Weight = edge_attr(final.igraph, "weight"))
    if (sum(final.scaled) > 0)
    {
        ppi.edges <- attr(E(ppi.igraph), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_")
        ppi.df <- data.frame(Edges = ppi.edges, Weight = edge_attr(ppi.igraph, "weight"))

        final.filter <- is.element(final.df$Edges, ppi.df$Edges)
        final.df[final.filter,]$Weight <- ppi.df$Weight 
        final.df[!final.filter,]$Weight <- 0
        final.df$Color <- "#dddddd99"
        final.df[final.filter,]$Color <- "#0000FF99"
    }
    else
    {
        final.df <- list(Color = "#0000FF99", Weight = 0)
    }

    CairoPDF(filename, width = plot.width, height = plot.height)
    par(mar=c(0,0,0,0) + 0.5)
    plot.igraph(final.igraph, layout = layout_nicely(final.igraph), vertex.size = 35, vertex.label.degree = pi/2, vertex.label.font = 2, vertex.label.color = "black", edge.color = final.df$Color, edge.width = 5*final.df$Weight)
    dev.off()

    return(ppi.igraph)
}

