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

load(file = "../baseline/save/intensities.rda")
load(file = "../baseline/save/targets.final.rda")
load(file = "../baseline/save/annot.rda")

#targets.final$Sample.Name %<>% str_replace(" ", "")
#colnames(intensities) %<>% str_replace(" ", "")
annot.reduce <- select(annot, Probe_Id, Accession, Symbol, Definition)
saveRDS.gz(annot.reduce, file = "./save/annot.reduce.rda")

source("../common_functions.R")
#Preprocessing of intenstities
targets.final.known <- filter(targets.final, Status != "Unknown" & !is.na(Draw.Age)) %>% droplevels
intensities.known <- select(intensities, one_of(targets.final.known$Sample.Name))
intensities.norm <- log2(intensities.known) %>% as.matrix %>% normalizeBetweenArrays
IAC.all <- gen.IACcluster("./IAC", intensities.norm, "All intensities normalized")
sd.all <- gen.sdplot("./sd", IAC.all, "All intensities normalized")

remove.names <- c("CHOP_122_1_Pat", "CHOP_63_2_Car", "FA_6438_2_Con")
remove.names.key <- paste(remove.names, collapse = "|")
intensities.rmout <- select(data.frame(intensities.norm), -(one_of(remove.names))) %>% as.matrix %>% normalizeBetweenArrays 
IAC.rmout <- gen.IACcluster("./IAC_rmout", intensities.rmout, "All intensities normalized")
sd.rmout <- gen.sdplot("./sd_rmout", IAC.rmout, "All intensities normalized")

reps <- sd.rmout[str_detect(names(sd.rmout), "1r")] %>% abs
orig.key <- str_replace_all(names(reps), "1r", "1")
orig <- sd.rmout[orig.key] %>% abs

reps.final <- data.frame(orig, reps)
max.key <- apply(reps.final, 1, which.max)
reps.names <- data.frame(orig.key, names(reps))
remove.names.reps <- reps.names[cbind(seq_along(max.key), max.key)]

remove.key.reps <- paste(remove.names.reps, collapse = "|")
intensities.rmreps <- select(data.frame(intensities.rmout), -one_of(remove.names.reps)) %>% as.matrix %>% normalizeBetweenArrays

remove.all <- paste(remove.names.key, remove.key.reps, sep = "|")
targets.final.known %<>% filter(!grepl(remove.all, Sample.Name)) %>% arrange(PIDN)

gen.days<- function(data.vector)
{
    return(as.numeric(data.vector - data.vector[1]))
}
diff.received <- by(targets.final.known, targets.final.known$PIDN, select, Date.Drawn) %>% llply(gen.days) %>% reduce(c)
diff.received[diff.received < 0] <- 0
targets.final.known$Diff.Baseline <- diff.received

order.key <- match(colnames(intensities.rmreps), targets.final.known$Sample.Name)
targets.final.known <- targets.final.known[order.key,]

targets.final.known %<>% filter(!(is.na(DOB)) & Sex != "UNKNOWN")
#targets.final.gaa <- filter(targets.final.known, !is.na(GAA1) & !is.na(GAA2))
saveRDS.gz(targets.final.known, file = "./save/targets.final.known.rda")
#saveRDS.gz(targets.final.gaa, file = "./save/targets.final.gaa.rda")

model.status.full <- model.matrix( ~ 0 + factor(targets.final.known$Status) )[,-2]
colnames(model.status.full) <- c("Carrier", "Patient")
model.sex.full <- model.matrix(~ 0 + factor(targets.final.known$Sex))[,-1]
targets.final.known$Draw.Age %<>% as.numeric 
draw.age.cov <- log2(targets.final.known$Draw.Age)
targets.final.known$Diff.Baseline %<>% as.numeric
#targets.final.known[targets.final.known$Diff.Baseline != 0,]$Diff.Baseline %<>% log2
diff.baseline.cov <- targets.final.known$Diff.Baseline
diff.baseline.cov[diff.baseline.cov != 0] %<>% log2
targets.final.known$RIN %<>% as.numeric
rin.cov <- log2(targets.final.known$RIN)
covariates.full <- cbind(Patient = model.status.full, Sex = model.sex.full, Age = draw.age.cov, Diff.Baseline = diff.baseline.cov, RIN = rin.cov)

intensities.rmreps <- select(data.frame(intensities.rmreps), one_of(targets.final.known$Sample.Name)) %>% as.matrix %>% normalizeBetweenArrays
saveRDS.gz(intensities.rmreps, file = "./save/intensities.rmreps.rda")
#intensities.gaa <- select(data.frame(intensities.rmreps), one_of(targets.final.gaa$Sample.Name)) %>% as.matrix %>% normalizeBetweenArrays
#saveRDS.gz(intensities.gaa, file = "./save/intensities.gaa.rda")
intensities.combat <- ComBat(dat = intensities.rmreps, batch = factor(targets.final.known$Batch), mod = covariates.full)
saveRDS.gz(intensities.combat, file = "./save/intensities.combat.rda")

gen.peer(10, intensities.combat, TRUE, covariates.full)
PEER.precision <- read_csv("./precision_10.csv")
PEER.precision$Factor <- rownames(PEER.precision)
data.plot <- PEER.precision[-c(1:5, 7),]
colnames(data.plot)[1] <- "Precision"
p <- ggplot(data.plot, aes(x = as.numeric(Factor), y = Precision, group = 1)) + geom_line()
CairoPDF("precision_reduced")
print(p)
dev.off()

PEER.weights.plot <- select(read_csv("./weight_10.csv"), -(X1:X7))
PEER.weights.sums <- colSums(abs(PEER.weights.plot)) %>% data.frame
PEER.weights.sums$Factor <- 1:nrow(PEER.weights.sums)
colnames(PEER.weights.sums)[1] <- "Weight"

p <- ggplot(PEER.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
p <- p + theme_bw() + xlab("Factor") + ylab("Weight")
CairoPDF("./PEER_weights", height = 4, width = 6)
print(p)
dev.off()

PEER.weights <- select(read_csv("./weight_10.csv"), -(X1:X5), -X7)
PEER.factors <- select(read_csv("./factor_10.csv"), -(X1:X5), -X7)
PEER.residuals <- as.matrix(PEER.factors) %*% t(as.matrix(PEER.weights)) %>% t
intensities.PEER <- intensities.combat - PEER.residuals
saveRDS.gz(intensities.PEER, file = "./save/intensities.peer.rda")

#Calculate coefficient of variation for all genes and then rank them using this
#expr.data.combat <- gen.cv(intensities.combat)
#saveRDS.gz(expr.data.combat, file = "./save/expr.data.combat.rda")
expr.data.PEER <- gen.cv(intensities.PEER)
saveRDS.gz(expr.data.PEER, file = "./save/expr.data.PEER.rda")
expr.data.out <- t(expr.data.PEER)
expr.data.out <- data.frame(Probe_Id = rownames(expr.data.out), expr.data.out)
write.csv(expr.data.out, "./expr.data.out.csv", row.names=FALSE)

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))
#sft.combat <- pickSoftThreshold(expr.data.combat, powerVector = powers, verbose = 5, networkType = "signed")
#sft.combat.df <- sft.combat$fitIndices
#saveRDS.gz(sft.combat, file = "./save/sft.combat.rda")

sft.PEER <- pickSoftThreshold(expr.data.PEER, powerVector = powers, verbose = 5, networkType = "signed")
sft.PEER.df <- sft.PEER$fitIndices
saveRDS.gz(sft.PEER, file = "./save/sft.peer.rda")

#sft.combat.df <- sft.combat$fitIndices
#sft.combat.df <- sft.PEER.df
#sft.combat.df <- sft.combat.df.bicor
#sft.combat.df <- sft.PEER.df.bicor

#Plot scale indendence and mean connectivity as functions of power
sft.PEER.df$multiplied <- sft.PEER.df$SFT.R.sq * -sign(sft.PEER.df$slope)
sft.combat.df$multiplied <- sft.combat.df$SFT.R.sq * -sign(sft.combat.df$slope)
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

#adjacency.PEER.bicor <- adjacency(expr.data.PEER, power = 12, type = "signed", corFnc = "bicor")
#saveRDS.gz(adjacency.PEER.bicor, file = "./save/adjacency.PEER.bicor.rda")

#adjacency.PEER.splineReg <- adjacency.splineReg(expr.data.PEER)

TOM.PEER <- TOMsimilarity(adjacency.PEER, verbose = 5)
dissimilarity.TOM <- 1 - TOM.PEER
saveRDS.gz(TOM.PEER, file = "./save/tom.PEER.rda")
saveRDS.gz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

#TOM.PEER.bicor <- TOMsimilarity(adjacency.PEER.bicor)
#dissimilarity.TOM.bicor <- 1 - TOM.PEER.bicor
#saveRDS.gz(TOM.PEER.bicor, file = "./save/tom.PEER.bicor.rda")
#saveRDS.gz(dissimilarity.TOM.bicor, file = "./save/dissimilarity.TOM.bicor.rda")

#geneTree.bicor = flashClust(as.dist(dissimilarity.TOM.bicor), method = "average")
#saveRDS.gz(geneTree.bicor, file = "./save/gene.tree.bicor.rda")

CairoPDF(file = "./1-genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

#CairoPDF(file = "./1-genecluster.bicor", height = 10, width = 15)
#plot(geneTree.bicor, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
#dev.off()

min.module.size <- 50

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, cutHeight = 0.995, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./1-gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#dynamic.modules.bicor <- cutreeDynamic(dendro = geneTree.bicor, method = "hybrid", distM = dissimilarity.TOM.bicor, cutHeight = 0.995, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
#dynamic.colors.bicor <- labels2colors(dynamic.modules.bicor)
#saveRDS.gz(dynamic.colors.bicor, file = "./save/dynamic.colors.bicolor.rda")

#CairoPDF(file = "./1-gene_dendrogram_and_module_colors_min100_bicor", height = 10, width = 15)
#plotDendroAndColors(geneTree.bicor, dynamic.colors.bicor, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene Dendrogram and Module Colors")
#dev.off()

#softPower.bicor <- 12
#softPower <- softPower.bicor
#dynamic.colors <- dynamic.colors.bicor
#adjacency.PEER <- adjacency.PEER.bicor
#geneTree <- geneTree.bicor

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

gene.module.membership.long <- gather(gene.module.membership, module.comparison, correlation, black:steelblue)
module.membership.pvalue.long <- gather(module.membership.pvalue, module.comparison, p.value, black:steelblue)
membership.join <- join(gene.module.membership.long, module.membership.pvalue.long)
eigengene.connectivity <- join(membership.join, gene.info) %>% select(Probe_Id, Accession:kscaled, module.comparison:p.value)
write_csv(eigengene.connectivity, "eigengene_connectivity.csv")

all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y") smooth.df <- data.frame(all.smooth)
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
#by(expr.data.plot, expr.data.plot$module.colors, gen.heatmap, ME.genes.plot)

modules.out <- select(expr.data.plot, module.colors)
modules.out$Probe_Id <- rownames(modules.out)
modules.out %<>% join(annot.reduce)
write.xlsx(modules.out, "modules_out.xlsx")
targets.final.known$Sample.Name %<>% str_replace(" ", "")

source("../GO/enrichr.R")

enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther", "Chromosome_Location") 
#enrichr.submit("darkmagenta", modules.out, enrichr.terms, FALSE)
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

rownames(ME.genes) <- rownames(expr.data.PEER)
rownames(PEER.factors) <- rownames(expr.data.PEER)
colnames(PEER.factors) <- paste("X", 1:ncol(PEER.factors), sep = "")

targets.status <- model.matrix( ~ 0 + factor(targets.final.known$Status) )[,-2] %>% data.frame
#targets.status <- targets.status[,c(2,1,3)]
colnames(targets.status) <- c("Carrier", "Patient")
targets.status %<>% mutate(Sample.Name = targets.final.known$Sample.Name)
cor.status <- gen.cor(ME.genes, targets.status)

targets.age <- select(targets.final.known, Sample.Name, Draw.Age) %>% filter(!is.na(Draw.Age))
targets.age$Draw.Age %<>% as.numeric
cor.age <- gen.cor(ME.genes, targets.age)
#cor.age.PEER <- gen.cor(PEER.factors, targets.age)

targets.sex <- filter(targets.final.known, Sex != "UNKNOWN")
targets.sex.m <- model.matrix( ~ 0 + factor(targets.sex$Sex) )[,-1] %>% data.frame #%>% mutate(Sample.Name = targets.sex$Sample.Name)
colnames(targets.sex.m) <- c("Sex")
targets.sex.m %<>% mutate(Sample.Name = targets.final.known$Sample.Name)
cor.sex <- gen.cor(ME.genes, targets.sex.m)
#cor.sex.PEER <- gen.cor(PEER.factors, targets.sex.m)

targets.gaa <- select(targets.final.known, Sample.Name, GAA1 ) %>% filter(!is.na(GAA1))# %>% filter(!is.na(GAA2))
#cor.gaa <- gen.cor(ME.genes, targets.gaa)
cor.gaa.PEER <- gen.cor(PEER.factors, targets.gaa)

targets.onset <- select(targets.final.known, Sample.Name, Onset) %>% filter(!is.na(Onset))
cor.onset.PEER <- gen.cor(PEER.factors, targets.onset)

targets.diff <- select(targets.final.known, Sample.Name, Diff.Baseline)
cor.diff <- gen.cor(ME.genes, targets.diff)

PEER.factors.plot <- select(read_csv("./factor_10.csv"), -(X1:X7))
names(PEER.factors.plot) <- paste("X", 1:ncol(PEER.factors.plot), sep = "")
PEER.factors.df <- mutate(PEER.factors.plot, Sample.Name = targets.final.known$Sample.Name) 
cor.PEER <- gen.cor(ME.genes, PEER.factors.df) %>% data.frame

module.traits.all <- cbind(cor.status, cor.age, cor.sex, cor.diff) %>% data.frame
module.traits.pval <- select(module.traits.all, contains("p.value")) %>% as.matrix
module.traits.cor <- select(module.traits.all, -contains("p.value")) %>% as.matrix

module.PEER.all <- cbind(cor.gaa.PEER, cor.onset.PEER) %>% data.frame
module.PEER.pval <- select(module.PEER.all, contains("p.value")) %>% as.matrix
module.PEER.cor <- select(module.PEER.all, -contains("p.value")) %>% as.matrix

cor.PEER.cor <- select(cor.PEER, -contains("p.value")) %>% as.matrix
cor.PEER.pval <- select(cor.PEER, contains("p.value")) %>% as.matrix

module.trait.out <- data.frame(Module = rownames(module.traits.cor), module.traits.cor, module.traits.pval)
module.PEER.out <- data.frame(Factor = rownames(module.PEER.cor), module.PEER.cor, module.PEER.pval)
cor.PEER.out <- data.frame(Module = rownames(cor.PEER), cor.PEER.cor, cor.PEER.pval)

write_csv(module.trait.out, "module_trait_cor.csv")
write_csv(module.PEER.out, "module_PEER_cor.csv")
write_csv(cor.PEER.out, "cor_PEER_cor.csv")

text.matrix.traits <- paste(signif(module.traits.cor, 2), '\n(', signif(module.traits.pval, 1), ')', sep = '')
text.matrix.PEER <- paste(signif(module.PEER.cor, 2), '\n(', signif(module.PEER.pval, 1), ')', sep = '')
text.matrix.PEERcor <- paste(signif(cor.PEER.cor, 2), '\n(', signif(cor.PEER.pval, 1), ')', sep = '')
dim(text.matrix.traits) = dim(module.traits.cor)
dim(text.matrix.PEER) = dim(module.PEER.cor)
dim(text.matrix.PEERcor) = dim(cor.PEER.cor)

gen.text.heatmap(module.traits.cor, text.matrix.traits, colnames(module.traits.cor), colnames(ME.genes), "", "module-trait relationships")
gen.text.heatmap(module.PEER.cor, text.matrix.PEER, colnames(module.PEER.cor), rownames(module.PEER.cor), "", "PEER factor-trait relationships")
gen.text.heatmap(cor.PEER.cor, text.matrix.PEERcor, colnames(cor.PEER.cor), colnames(ME.genes), "", "module-PEER factor relationships")

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
