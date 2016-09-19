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
library(biomaRt)
library(samr)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(heatmap.plus)

#Reading and writing tables
library(readr)
library(openxlsx)
library(parallel)

#Functional programming
library(magrittr)
library(purrr)

#Data arrangement
library(dplyr)
library(tidyr)
library(broom)

#Boxplot
Boxplot <- function(filename, lumi.object, maintext, ylabtext) {
    #dataset %<>% t %>% data.frame
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Batch = lumi.object$Batch)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Batch"))
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value, fill = factor(Batch))) + geom_boxplot() + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Oxygen", width = 20 , height = 8)
}

#Histogram
Histogram <- function(filename, lumi.object) {
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Name = sampleNames(lumi.object), Status = lumi.object$Status)
    dataset.m <- gather(dataset.addvars, nuID, Expression, -Sample.Name, -Status)

    p <- ggplot(dataset.m, aes(Expression, group = Sample.Name, col = factor(Status))) + geom_density() + theme_bw()
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("Histogram of VST Expression") + ylab("Density") + xlab("VST Expression") 
    CairoPDF(filename, height = 5, width = 9)
    print(p)
    dev.off()
}

#MDS function 
MDSPlot <- function(filename, dataset, targetset, colorscheme = "none", variablename) {
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$Slide.ID, factor(targetset[[variablename]]))
    colnames(target.data) <- c("Slide.ID", variablename)
    colnames(dataset.plot) <- c("Slide.ID", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if (colorscheme != "none") {
        p <- p + scale_color_manual(values = colorscheme) 
    }
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename)
    CairoPDF(file = filename, height = 7, width = 7)
    print(p)
    dev.off()
}

#Samplewise connectivity plot
ConnectivityPlot <- function(filename, dataset, maintitle) {
    norm.adj <- (0.5 + 0.5 * bicor(exprs(dataset)))
    colnames(norm.adj) <- dataset$Slide.ID
    rownames(norm.adj) <- dataset$Slide.ID
    net.summary <- fundamentalNetworkConcepts(norm.adj)
    net.connectivity <- net.summary$Connectivity
    connectivity.zscore <- (net.connectivity - mean(net.connectivity)) / sqrt(var(net.connectivity))

    connectivity.plot <- data.frame(Slide.ID = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))
    p <- ggplot(connectivity.plot, aes(x = Sample.Num, y = Z.score, label = Slide.ID) )
    p <- p + geom_text(size = 4, colour = "red")
    p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
    CairoPDF(filename, width = 10, height = 10)
    print(p)
    dev.off()
    return(connectivity.zscore)
}

CorrelationWorkbook <- function(de.table, filename) { 
    pval.cols <- colnames(de.table) %>% str_detect("P.value") %>% which
    cor.cols <- colnames(de.table) %>% str_detect("Correlation") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = de.table)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(de.table), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = cor.cols, rows = 1:nrow(de.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:6, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 11)
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

GeneScatterplot <- function(lumi.object, gene.symbol) {
    gene.expr <- as.vector(exprs(lumi.object[gene.symbol,]))
    gene.df <- data.frame(GAA1 = lumi.object$GAA1, Expression = gene.expr)

    p <- ggplot(gene.df, aes(x = GAA1, y = Expression)) + geom_point() + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
    p <- p + theme(plot.margin = unit(c(1,1,1,1), "lines")) + ggtitle(gene.symbol)
    p <- p + stat_smooth(method = "lm", se = TRUE) + scale_x_continuous(as.numeric(unique(lumi.object$GAA1)))
    CairoPDF(gene.symbol, width = 6, height = 6)
    print(p)
    dev.off()
}

EigengeneCorrelation <- function(dataset, trait.df) {
    module.trait.cor <- bicor(dataset, trait.df, maxPOutliers = 0.05)
    module.trait.cor.pval <- corPvalueStudent(module.trait.cor, length(dataset))
    return(c(module.trait.cor, module.trait.cor.pval))
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

EigengeneScatterplot <- function(eigengene.df, gaa.vector, color) {
    color.column <- str_c("ME", color)
    gene.df <- data.frame(GAA1 = gaa.vector, Expression = as.vector(eigengene.df[[color.column]]))

    p <- ggplot(gene.df, aes_string(x = "GAA1", y = "Expression")) + geom_point(col = color) + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank(), axis.title.x = element_blank()) 
    p <- p + stat_smooth(method = "lm", se = TRUE) + scale_x_continuous(as.numeric(unique(gaa.vector)))
    CairoPDF(color, width = 6, height = 6)
    print(p)
    dev.off()
}

#Enrichr 
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

source("../common_functions.R")

lumi.import <- readRDS.gz(file = "../baseline_lumi/save/lumi.baseline.rda")
lumi.patient <- lumi.import[,lumi.import$Status == "Patient" & !is.na(lumi.import$GAA1) & lumi.import$Age < 35 & lumi.import$Age > 11]
batch.colors <- data.frame(Batch = factor(1:19), Color = c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green", "blue4", "red4", "green4")) #Assign batch colors
lumi.patient$Batch.Color <- left_join(select(pData(lumi.patient), Batch), batch.colors) %>% select(Color) #add variable for batch color

#gen.boxplot("patient_intensity_raw.jpg", lumi.vst, "VST Transformed Intensity", "Intensity")
Histogram("patient_histogram", lumi.patient)
mds.patient <- exprs(lumi.patient) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("patient_mds_vst", mds.patient, pData(lumi.patient), as.character(batch.colors$Color), "Batch") #label PCs by vst

lumi.patient.norm <- lumiN(lumi.patient, method = "rsn") #Normalize with robust spline regression
lumi.patient.cutoff <- detectionCall(lumi.patient.norm) #Get the count of probes which passed the detection threshold per sample
lumi.patient.expr <- lumi.patient.norm[which(lumi.patient.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.patient <- getSYMBOL(rownames(lumi.patient.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.patient.annot <- lumi.patient.expr[!symbols.lumi.patient,] #Drop any probe which is not annotated
SaveRDSgz(lumi.patient.annot, file = "./save/lumi.patient.annot.rda")

#Regenerate plots
#gen.boxplot("baseline_intensity_norm.jpg", lumi.expr.annot, batch.colors, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
Histogram("histogram_norm", lumi.patient.annot)
mds.norm <- exprs(lumi.patient.annot) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_batch_norm", mds.norm, pData(lumi.patient.annot), as.character(batch.colors$Color), "Batch") #label PCs by batch

#Use ComBat for batch effect correction
lumi.patient.annot$Sex %<>% droplevels
model.combat <- model.matrix(~ GAA1 + Sex + Age + RIN, pData(lumi.patient.annot)) %>% data.frame

expr.combat <- ComBat(dat = exprs(lumi.patient.annot), batch = factor(lumi.patient.annot$Site))
lumi.combat <- lumi.patient.annot
exprs(lumi.combat) <- expr.combat

expr.combat <- ComBat(dat = exprs(lumi.combat), batch = factor(lumi.combat$Batch))
exprs(lumi.combat) <- expr.combat
SaveRDSgz(lumi.combat, file = "./save/lumi.combat.rda")

#Regenerate plots
#gen.boxplot("patient_intensity_combat.jpg", lumi.combat, "VST Transformed Intensity (Batch Corrected)", "Intensity")
Histogram("patient_histogram_combat", lumi.combat)
mds.patient.combat <- exprs(lumi.combat) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_batch_combat", mds.patient.combat, pData(lumi.combat), as.character(batch.colors$Color), "Batch") #label PCs by batch

#Remove outliers
combat.connectivity <- ConnectivityPlot("patient_connectivity_combat", lumi.combat, "")
connectivity.outlier <- combat.connectivity[combat.connectivity < -2]
outlier.names <- names(connectivity.outlier) %>% is.element(el = sampleNames(lumi.combat))

lumi.rmout <- lumi.combat[,!outlier.names]

#Regenerate plots
#gen.boxplot("patient_intensity_rmout.jpg", lumi.rmout, "VST Transformed Intensity (Outliers Removed)", "Intensity")
Histogram("patient_histogram_rmout", lumi.rmout)
mds.patient.rmout <- exprs(lumi.rmout) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_batch_rmout", mds.patient.rmout, pData(lumi.rmout), as.character(batch.colors$Color), "Batch") #label PCs by batch

#Check for collinearity
age.gaa <- lm(as.numeric(Age) ~ GAA1, pData(lumi.rmout)) %>% anova %>% tidy
sex.gaa <- lm(as.integer(Sex) ~ GAA1, pData(lumi.rmout)) %>% anova %>% tidy
site.gaa <- lm(as.integer(factor(Site)) ~ GAA1, pData(lumi.rmout)) %>% anova %>% tidy
batch.gaa <- lm(as.integer(Batch) ~ GAA1, pData(lumi.rmout)) %>% anova %>% tidy
RIN.gaa <- lm(RIN ~ GAA1, pData(lumi.rmout)) %>% anova %>% tidy

p <- ggplot(pData(lumi.rmout), aes(x = GAA1, y = Age)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("Age (p < ", signif(age.gaa$p.value[1], 3), "*)", sep = ""))
CairoPDF("./age_scatterplot.pdf", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = GAA1, y = RIN)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("RIN (p < ", signif(RIN.gaa$p.value[1], 3), ")", sep = ""))
CairoPDF("./rin_boxplot.pdf", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = factor(Sex), y = GAA1, fill = factor(Sex))) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Sex (p < ", signif(sex.gaa$p.value[1], 3), "*)", sep = "")) 
CairoPDF("./sex_barplot.pdf", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = factor(Batch), y = GAA1, fill = factor(Batch))) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.position = "none") 
p <- p + ggtitle(paste("Batch (p < ", signif(batch.gaa$p.value[1], 3), "*)", sep = "")) 
CairoPDF("./batch_barplot.pdf", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = factor(Site), y = GAA1, fill = factor(Site))) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.position = "none") 
p <- p + ggtitle(paste("Site (p < ", signif(site.gaa$p.value[1], 3), "*)", sep = "")) 
CairoPDF("./site_barplot.pdf", height = 6, width = 6)
plot(p)
dev.off()

model.rmout <- model.matrix(~ Sex + Age + RIN, pData(lumi.rmout)) %>% data.frame

#Remove effects of covariates
cleaned.expr <- removeBatchEffect(exprs(lumi.rmout), covariates = model.rmout[,2:4])
lumi.cleaned <- lumi.rmout
exprs(lumi.cleaned) <- cleaned.expr
SaveRDSgz(lumi.cleaned, file = "./save/lumi.cleaned.rda")

#Collapse by gene symbol
expr.collapse <- collapseRows(exprs(lumi.cleaned), getSYMBOL(featureNames(lumi.cleaned), 'lumiHumanAll.db'), rownames(lumi.cleaned))$datETcollapsed
export.lumi <- ExpressionSet(assayData = expr.collapse, phenoData = phenoData(lumi.cleaned))
SaveRDSgz(export.lumi, "./save/export.lumi.rda")

#Calculate raw gene correlations with GAA1
cor.gaa <- apply(exprs(export.lumi), 1, cor, export.lumi$GAA1)
cor.gaa.pval <- corPvalueStudent(cor.gaa, length(export.lumi$GAA1)) 
cor.gaa.adjpval <- p.adjust(cor.gaa.pval, "fdr")

cor.gaa.df <- cbind(cor.gaa, cor.gaa.pval, cor.gaa.adjpval) %>% data.frame
colnames(cor.gaa.df) <- c("Correlation", "P.value", "Adj.P.value")
cor.gaa.df$Symbol <- featureNames(export.lumi)
cor.gaa.df %<>% arrange(P.value)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(featureNames(export.lumi)), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")
cor.gaa.annot <- left_join(cor.gaa.df, bm.table) %>% select(Symbol, Definition, Correlation:Adj.P.value)

CorrelationWorkbook(cor.gaa.df, "gaa.cor.xlsx")

samr.gaa <- SAM(x = exprs(export.lumi), y = export.lumi$GAA1, nperms = 1000, resp.type = "Quantitative", center.arrays = TRUE, genenames = featureNames(export.lumi), testStatistic = "standard", random.seed = 12345, logged2 = TRUE, fdr.output = 0.05)
write.xlsx(samr.gaa$siggenes.table$genes.up, "samr.gaa.up.xlsx")
write.xlsx(samr.gaa$siggenes.table$genes.lo, "samr.gaa.down.xlsx")
SaveRDSgz(samr.gaa, "./samr.gaa.rda")

#Top GAA1 correlated genes
#Negative
GeneScatterplot(export.lumi, "DGKD")
GeneScatterplot(export.lumi, "DBNDD1")

#Positive
GeneScatterplot(export.lumi, "EIF1")
GeneScatterplot(export.lumi, "TMOD1")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.data <- t(expr.collapse)
sft <- pickSoftThreshold(expr.data, powerVector = powers, verbose = 5, networkType = "signed", corFnc = bicor, corOptions = list(maxPOutliers = 0.05))
sft.df <- sft$fitIndices
SaveRDSgz(sft, file = "./save/sft.rda")

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

softPower <- 7
adjacency.expr <- adjacency(expr.data, power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")

TOM <- TOMsimilarity(adjacency.expr)
dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
SaveRDSgz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

#Identify modules using dynamic tree cutting with hybrid clustering
min.module.size <- 20
dynamic.modules <- cutreeDynamic(dendro = geneTree, cutHeight = 0.995, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
SaveRDSgz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.data, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
SaveRDSgz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.20
merge.all <- mergeCloseModules(expr.data, dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3, corFnc = bicor, corOptions = list(maxPOutliers = 0.05)) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
SaveRDSgz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
SaveRDSgz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
SaveRDSgz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 10)
par(cex = 0.7)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.adjacency = 0.3, cex.preservation = 0.3, plotPreservation = "standard")
dev.off()

#Determine eigengene significance
color.values <- unique(module.colors)
lumi.pdata <- pData(export.lumi)

cor.gaa <- map(ME.genes, EigengeneCorrelation, lumi.pdata$GAA1) %>% reduce(rbind) %>% data.frame
colnames(cor.gaa) <- c("Correlation", "P.value")
cor.gaa$Adj.P.value <- p.adjust(cor.gaa$P.value, method = "fdr")
cor.gaa$Module <- colnames(ME.genes)

text.matrix.traits <- str_c(signif(cor.gaa$Correlation, 2), '\n(', signif(cor.gaa$Adj.P.value, 1), ')')
dim(text.matrix.traits) = dim(as.matrix(cor.gaa$Correlation))

CairoPDF("module_trait_relationships", width = 4, height = 10)
par(mar = c(8, 8, 3, 3))
labeledHeatmap(Matrix = as.matrix(cor.gaa$Correlation), xLabels = "GAA1", yLabels = colnames(ME.genes), ySymbols = colnames(ME.genes), yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix.traits, zlim = c(-1,1), main = "")
dev.off()

#Eigengene plots
EigengeneScatterplot(ME.genes, lumi.pdata$GAA1, "brown")

#Generate network statistics
all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Symbol = rownames(all.degrees), Module = module.colors, all.degrees, stringsAsFactors = FALSE) %>% arrange(Module)
gene.info$kscaled <- by(gene.info, gene.info$Module, select, kWithin) %>% map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)
gene.info.annot <- left_join(gene.info, bm.table) %>% select(Symbol, Definition, Module:kscaled)

gene.module.membership <- data.frame(bicor(expr.data, ME.genes, maxPOutliers = 0.05))
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.data)))
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

module.membership <- left_join(gene.info.annot, gene.module.membership) %>% left_join(module.membership.pvalue, by = "Symbol")
ModuleWorkbook(module.membership, "module_membership.xlsx")

#PCA Plots
all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% map(extract2, "y")
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$Gene <- as.factor(1:nrow(smooth.df))
smooth.plot <- gather(smooth.df, Module, PCA1, -Gene)

PCAPlot("all_principal_components", smooth.plot, FALSE, 10, 15)
PCAPlot("facet_principal_components", smooth.plot, TRUE, 13, 20)

#Heatmaps
sample.ids <- factor(rownames(expr.data), levels = rownames(expr.data))
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.data), Module = module.colors)
split(expr.data.plot, expr.data.plot$Module) %>% map(MEHeatmap, ME.genes.plot)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#Enrichr
source("../../code/GO/enrichr.R")
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "HumanCyc_2016", "NCI-Nature_2016", "Panther_2016")
enrichr.submit("pink", modules.out, enrichr.terms, FALSE)
color.names <- unique(module.colors) %>% sort
l_ply(color.names, enrichr.submit, modules.out, enrichr.terms, FALSE)

black.only <- filter(modules.out, module.color == "black")
black.gobiol.file <- "./enrichr/black/black_GO_Biological_Process_2015.xlsx"
black.gobiol <- read.xlsx(black.gobiol.file) 
black.gobiol$Num.Genes <- map(black.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
black.gobiol$Database <- "GO Biological Process"
get.kappa.cluster(black.gobiol, black.only$Symbol, file_path_sans_ext(black.gobiol.file))
black.gobiol.final <- slice(black.gobiol, c(8, 23, 6, 14, 13, 24))

black.gomole.file <- "./enrichr/black/black_GO_Molecular_Function_2015.xlsx"
black.gomole <- read.xlsx(black.gomole.file) 
black.gomole$Num.Genes <- map(black.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.gomole %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
black.gomole$Database <- "GO Molecular Function"
get.kappa.cluster(black.gomole, black.only$Symbol, file_path_sans_ext(black.gomole.file))
black.gomole.final <- slice(black.gomole, c(10))

black.reactome.file <- "./enrichr/black/black_Reactome_2016.xlsx"
black.reactome <- read.xlsx(black.reactome.file) 
black.reactome$Num.Genes <- map(black.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.reactome %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
black.reactome$Database <- "Reactome"
get.kappa.cluster(black.reactome, black.only$Symbol, file_path_sans_ext(black.reactome.file))
black.reactome.final <- slice(black.reactome, c(28, 39))

black.kegg <- read.xlsx("./enrichr/black/black_KEGG_2016.xlsx") %>% slice(39)
black.kegg$Num.Genes <- map(black.kegg$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.kegg$Database <- "KEGG"

black.enrichr <- rbind(black.gobiol.final, black.gomole.final, black.reactome.final, black.kegg)
gen.enrichrplot(black.enrichr, "black.enrichr")

#get kegg

blue.only <- filter(modules.out, module.color == "blue")
blue.gobiol.file <- "./enrichr/blue/blue_GO_Biological_Process_2015.xlsx"
blue.gobiol <- read.xlsx(blue.gobiol.file) 
blue.gobiol$Num.Genes <- map(blue.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
blue.gobiol$Database <- "GO Biological Process"
blue.gobiol.filter <- filter(blue.gobiol, Num.Genes > 4) %>% filter(P.value < 0.01)
get.kappa.cluster(blue.gobiol.filter, blue.only$Symbol, file_path_sans_ext(blue.gobiol.file))
blue.gobiol.final <- slice(blue.gobiol, c(19, 7, 6))

blue.gomole.file <- "./enrichr/blue/blue_GO_Molecular_Function_2015.xlsx"
blue.gomole <- read.xlsx(blue.gomole.file) 
blue.gomole$Num.Genes <- map(blue.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
blue.gomole$Database <- "GO Molecular Function"
blue.gomole.filter <- filter(blue.gomole, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(blue.gomole.filter, blue.only$Symbol, file_path_sans_ext(blue.gomole.file))
blue.gomole.final <- slice(blue.gomole, c(5, 13, 3, 15))

blue.reactome.file <- "./enrichr/blue/blue_Reactome_2016.xlsx"
blue.reactome <- read.xlsx(blue.reactome.file) 
blue.reactome$Num.Genes <- map(blue.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
blue.reactome$Database <- "Reactome"
blue.reactome.filter <- filter(blue.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(blue.reactome.filter, blue.only$Symbol, file_path_sans_ext(blue.reactome.file))
blue.reactome.final <- slice(blue.reactome, c(7, 1, 24, 38))

#get kegg

blue.enrichr <- rbind(blue.gobiol.final, blue.gomole.final, blue.reactome.final)
gen.enrichrplot(blue.enrichr, "blue.enrichr")

pink.only <- filter(modules.out, module.color == "pink")
pink.gobiol.file <- "./enrichr/pink/pink_GO_Biological_Process_2015.xlsx"
pink.gobiol <- read.xlsx(pink.gobiol.file) 
pink.gobiol$Num.Genes <- map(pink.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pink.gobiol$Database <- "GO Biological Process"
pink.gobiol.final <- filter(pink.gobiol, Num.Genes > 4) %>% filter(P.value < 0.01)
get.kappa.cluster(pink.gobiol.final, pink.only$Symbol, file_path_sans_ext(pink.gobiol.file))
pink.gobiol.final <- slice(pink.gobiol, c(1, 72, 116))

pink.gomole.file <- "./enrichr/pink/pink_GO_Molecular_Function_2015.xlsx"
pink.gomole <- read.xlsx(pink.gomole.file) 
pink.gomole$Num.Genes <- map(pink.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pink.gomole$Database <- "GO Molecular Function"
pink.gomole.filter <- filter(pink.gomole, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(pink.gomole.filter, pink.only$Symbol, file_path_sans_ext(pink.gomole.file))
pink.molec.final <- slice(pink.gomole, c(1))

#pink.reactome.file <- "./enrichr/pink/pink_Reactome_2016.xlsx"
#pink.reactome <- read.xlsx(pink.reactome.file) 
#pink.reactome$Num.Genes <- map(pink.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
#pink.reactome$Database <- "Reactome"
#pink.reactome.final <- filter(pink.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
#get.kappa.cluster(pink.reactome.final, pink.only$Symbol, file_path_sans_ext(pink.reactome.file))

pink.kegg.file <- "./enrichr/pink/pink_KEGG_2016.xlsx"
pink.kegg <- read.xlsx(pink.kegg.file) 
pink.kegg$Num.Genes <- map(pink.kegg$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pink.kegg$Database <- "kegg"
pink.kegg.final <- slice(pink.kegg, c(5, 10))

pink.final <- rbind(pink.gobiol.final, pink.molec.final, pink.kegg.final)
gen.enrichrplot(pink.final, "pink.enrichr")

#PPI - fix me!
biogrid.ppi <- read_tsv("../../code/BIOGRID-ORGANISM-3.4.136.tab2/BIOGRID-ORGANISM-Homo_sapiens-3.4.136.tab2.txt") %>% data.frame
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

