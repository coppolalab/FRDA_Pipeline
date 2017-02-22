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
library(BayesFactor)
library(siggenes)
library(BayesianFirstAid)

#String operations
library(stringr)
library(tools)
library(R.utils)

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
library(plyr)
library(tidyr)
library(broom)
library(UpSetR)

#Boxplot
Boxplot <- function(filename, lumi.object, maintext, ylabtext) {
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
    #cor.cols <- colnames(de.table) %>% str_detect("Correlation") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = de.table)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(de.table), rule = "<0.05", style = sig.pvalues)
    #conditionalFormatting(wb, 1, cols = cor.cols, rows = 1:nrow(de.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:ncol(de.table), widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 11)
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

BayesWorkbook <- function(de.table, filename) { 
    pval.cols <- colnames(de.table) %>% str_detect("Bayes.Factor") %>% which
    cor.cols <- colnames(de.table) %>% str_detect("Median") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = de.table)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(de.table), rule = ">3", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = cor.cols, rows = 1:nrow(de.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:ncol(de.table), widths = 15)
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
    p <- p + stat_smooth(method = "lm", se = TRUE) 
    p <- p + theme(axis.ticks.x = element_blank()) + xlab("GAA1 Expansion Size (bp)") + ylab("VST Normalized Expression")
    p <- p + theme(plot.margin = unit(c(1,1,1,1), "lines"), plot.background = element_blank()) + ggtitle(gene.symbol)
    p <- p + ggtitle(gene.symbol) + theme(panel.border = element_rect(color = "black", size = 1))
    CairoPDF(str_c(gene.symbol, ".pdf"), width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

GetBF <- function(gene.vector, model.design) {
    final.df <- mutate(model.design, Gene = scale(gene.vector))
    bf.all <- regressionBF(Gene ~ GAA1, final.df)
    bf.all
}

GetLM <- function(gene.vector, model.design) {
    final.df <- mutate(model.design, Gene = scale(gene.vector))
    lm.all <- lm(Gene ~ GAA1, final.df)
    lm.all
}

EigengeneCorrelation <- function(dataset, trait.df) {
    module.trait.cor <- bicor(dataset, trait.df, maxPOutliers = 0.05)
    module.trait.cor.pval <- corPvalueStudent(module.trait.cor, length(dataset))
    c(module.trait.cor, module.trait.cor.pval)
}

GeneCorrelation <- function(dataset, trait.df) {
    module.trait.cor <- bicor(dataset, trait.df, maxPOutliers = 0.05)
    module.trait.cor
}

EigengeneBayes <- function(dataset, trait.df) {
    module.trait.cor <- bayes.cor.test(dataset, trait.df)
    #module.trait.cor.pval <- corPvalueStudent(module.trait.cor, length(dataset))
    #c(module.trait.cor, module.trait.cor.pval)
    return(module.trait.cor)
}

EigengeneBF <- function(dataset, trait.df) {
    lm.df <- data.frame(Eigengene = scale(dataset), Trait = trait.df)
    module.trait.cor <- regressionBF(Eigengene ~ Trait, lm.df)
    return(module.trait.cor)
}

ModuleWorkbook <- function(module.table, filename) {
    pval.key <- colnames(module.table) %>% str_detect("pvalue") 
    pval.cols <- which(pval.key)
    MM.key <- colnames(module.table) %>% str_detect("MM.") 
    MM.cols <- which(MM.key)

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

EigengeneScatterplot <- function(eigengene.df, gaa.vector, color) {
    color.column <- str_c("ME", color)
    gene.df <- data.frame(GAA1 = gaa.vector, Expression = as.vector(eigengene.df[[color.column]]))

    p <- ggplot(gene.df, aes_string(x = "GAA1", y = "Expression")) + geom_point(col = color) + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("GAA1 Expansion Size (bp)") + ylab("Eigengene Value")
    p <- p + theme(axis.ticks.x = element_blank(), plot.margin = unit(c(1,1,1,1), "lines")) 
    p <- p + stat_smooth(method = "lm", se = TRUE) 
    p <- p + ggtitle(str_c(capitalize(color), " Module")) + theme(plot.background = element_blank(), panel.border = element_rect(color = "black", size = 1))
    CairoPDF(str_c(color, ".pdf"), width = 4, height = 4, bg = "transparent")
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

EnrichrPlot <- function(enrichr.df, filename, plot.height = 5, plot.width = 8, maintitle = "auto") {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Log.Bayes.Factor <- log10(enrichr.df$Bayes.Factor)
    enrichr.df$Term %<>% str_replace_all("\\ \\(GO.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% str_replace_all(",.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.Bayes.Factor)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Bayes.Factor)  

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Log.Bayes.Factor)) + geom_bar(stat = "identity", fill = "brown") + coord_flip() + theme_bw() + theme(legend.position = "none") 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' Bayes Factor')))
    p <- p + theme(panel.border = element_rect(color = "black", size = 1), plot.background = element_blank(), plot.title = element_text(hjust = 0.5))  
    if (maintitle == "auto"){
        p <- p + ggtitle(str_c(capitalize(filename), " Module"))
    } else {
        p <- p + ggtitle(maintitle)
    }
    CairoPDF(str_c(filename, ".enrichr"), height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

GetKscaled <- function(gene.list, module.membership) {
    filter(module.membership, is.element(Symbol, gene.list)) %>% select(Symbol, kscaled) %>% arrange(desc(kscaled))
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

#Check for collinearity - fix tomorrow
age.gaa <- regressionBF(Age ~ GAA1, data = pData(lumi.rmout)) 
sex.gaa <- anovaBF(Sex ~ GAA1, pData(lumi.rmout))
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
cleaned.expr <- removeBatchEffect(exprs(lumi.rmout), covariates = model.rmout[,1:4])
lumi.cleaned <- lumi.rmout
exprs(lumi.cleaned) <- cleaned.expr
SaveRDSgz(lumi.cleaned, file = "./save/lumi.cleaned.rda")

#Collapse by gene symbol
expr.collapse <- collapseRows(exprs(lumi.cleaned), getSYMBOL(featureNames(lumi.cleaned), 'lumiHumanAll.db'), rownames(lumi.cleaned))$datETcollapsed
export.lumi <- ExpressionSet(assayData = expr.collapse, phenoData = phenoData(lumi.cleaned))
SaveRDSgz(export.lumi, "./save/export.lumi.rda")

model.design <- model.matrix(~ scale(GAA1), pData(export.lumi)) #Make covariate matrix for limma
colnames(model.design)[2] <- "GAA1"
SaveRDSgz(model.design, file = "./save/model.design.rda")

##Bicor
#bicor.gaa <- apply(expr.collapse, 1, GeneCorrelation, export.lumi$GAA1)
#bicor.df <- data.frame(Symbol = featureNames(export.lumi), Bicor = bicor.gaa) %>% arrange(desc(Bicor))
#write.xlsx(bicor.df, "bicor.xlsx")
#skew.gaa <- apply(expr.collapse, 1, skewness)
#skew.df <- data.frame(Symbol = names(skew.gaa), Skewness = skew.gaa)
#skew.quantile <- quantile(skew.gaa, c(0.025, 0.975))
#skew.filter <- filter(skew.df, Skewness > skew.quantile[1] & Skewness < skew.quantile[2]) 
#SaveRDSgz(skew.filter, "./save/skew.filter.rda")

##BF test
bf.gaa <- apply(expr.collapse, 1, GetBF, data.frame(model.design))
SaveRDSgz(bf.gaa, "./save/bf.gaa.rda")
bf.extract <- map(bf.gaa, extractBF) %>% map_dbl(extract2, "bf")
bf.df <- data.frame(Symbol = featureNames(export.lumi), Bayes.Factor = bf.extract) %>% arrange(desc(Bayes.Factor))
bf.sig <- filter(bf.df, Bayes.Factor > 3)$Symbol
bf.posterior <- map(bf.gaa, posterior, iterations = 10000)
SaveRDSgz(bf.posterior, "./save/bf.posterior.rda")
bf.quantile <- map(bf.posterior, magrittr::extract, TRUE, 1) %>% map(quantile, c(0.025, 0.975)) %>% reduce(rbind) %>% set_colnames(c("CI_2.5", "CI_97.5")) %>% data.frame
bf.quantile$Median <- map(bf.posterior, magrittr::extract, TRUE, 1) %>% map_dbl(median)
bf.quantile$Symbol <- featureNames(export.lumi)
bf.posterior.df <- left_join(bf.df, bf.quantile)
write.xlsx(bf.posterior.df, "bf.gaa.xlsx")

#lumi.filter <- export.lumi[skew.filter$Symbol,]
##Limma fit
#fit <- lmFit(scale(exprs(lumi.filter)), model.design)
#fit$df.residual <- (fit$df.residual - 5)
#fitb <- eBayes(fit)
#fitb$t[,2] %<>% abs #Something weird going on here?

#top.gaa <- topTable(fitb, coef = 2, n = Inf)
#top.gaa$Symbol <- rownames(top.gaa)
#ebam.gaa <- limma2ebam(fitb, coef = 2, moderate = TRUE, delta = 0.9)  
#ebam2excel(ebam.gaa, 0.9, "ebam.gaa.csv")

#ebam.gaa.df <- data.frame(Z.score = signif(ebam.gaa@z, 3), Posterior = signif(ebam.gaa@posterior, 3), Coefficient = signif(fitb$coef[,2], 3))
#ebam.gaa.df$Symbol <- rownames(ebam.gaa.df)
#ebam.gaa.df %<>% arrange(desc(Posterior))
#write.xlsx(ebam.gaa.df, "ebam.gaa.xlsx")
#SaveRDSgz(ebam.gaa.df, "ebam.gaa.df.rda")
#ebam.annot <- left_join(ebam.gaa.df, bm.table) %>% select(Symbol, Definition, Coefficient, Posterior, Z.score)

#Annotate top table
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(featureNames(export.lumi)), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")
bf.gaa.annot <- left_join(bf.posterior.df, bm.table) %>% select(Symbol, Definition, Bayes.Factor, Median, CI_2.5, CI_97.5)
BayesWorkbook(bf.gaa.annot, "bf.gaa.xlsx")

bf.gaa.plot <- mutate(bf.posterior.df, Log.Bayes.Factor = log10(bf.posterior.df$Bayes.Factor))

BayesPlot <- function(siggene, filename, threshold, plot.title, posterior.column = "Posterior", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Posterior Probability") {
    siggene$Significant <- siggene[[posterior.column]] > threshold 
    siggene$Significant %<>% factor(levels = c(TRUE, FALSE)) %>% revalue(c("TRUE" = "Yes", "FALSE" = "No"))
    p <- ggplot(siggene, aes_string(x = log.column, y = posterior.column)) + geom_point(aes(color = Significant))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_blank())
    p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"), plot.title = element_text(hjust = 0.5))
    p <- p + xlab(xlabel) + ylab(ylabel) + ggtitle(plot.title) + scale_color_manual(values = c("orange", "darkgreen"), name = "Significant", labels = c("Yes", "No"))
    CairoPDF(filename, width = 5, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

BayesPlot(bf.gaa.plot, "gaa_volcano", 0.5, "GAA1", "Log.Bayes.Factor", "Median", "Regression Coefficient", "Log Bayes Factor")

#Enrichr
gaa.gobiol.file <- "./enrichr/gaa/GO Biological Process.xlsx"
gaa.gobiol <- read.xlsx(gaa.gobiol.file) 
GetKappaCluster(file_path_sans_ext(gaa.gobiol.file), gaa.gobiol, modules.only$Symbol, "Bayes.Factor")

gaa.gobiol.final <- slice(gaa.gobiol, c(1, 3, 13, 17))
EnrichrPlot(gaa.gobiol.final, "gaa", plot.height = 3, plot.width = 6, "GAA1")

#Top GAA1 correlated genes
#Positive
GeneScatterplot(export.lumi, "EIF1")
GeneScatterplot(export.lumi, "TMOD1")

#Negative
GeneScatterplot(export.lumi, "LIN7B")
GeneScatterplot(export.lumi, "MKL2")

#Biggest logFC (figure out!)
GeneScatterplot(export.lumi, "CA1")
GeneScatterplot(export.lumi, "ALAS2")

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

CairoPDF(file = "./gene_dendrogram_and_module_colors", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.data, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
SaveRDSgz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.20
merge.all <- mergeCloseModules(expr.data, dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3, corFnc = bicor, corOptions = list(maxPOutliers = 0.05)) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering", height = 8, width = 10, bg = "transparent")
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("", ""), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Patients Only", cex.main = 4.0,  axes = FALSE, ylab = "")
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

bayes.gaa <- map(ME.genes, EigengeneBF, scale(lumi.pdata$GAA1))
bf.eigengene <- map(bayes.gaa, extractBF) %>% map_dbl(extract2, "bf")
posterior.eigengene <- map(bayes.gaa, posterior, iterations = 10000)
quantile.eigengene <- map(posterior.eigengene, magrittr::extract, TRUE, 1) %>% map(quantile, c(0.025, 0.975)) %>% reduce(rbind) %>% set_colnames(c("CI_2.5", "CI_97.5")) %>% data.frame
quantile.eigengene$Median <- map(posterior.eigengene, magrittr::extract, TRUE, 1) %>% map_dbl(median)
quantile.eigengene$Color <- colnames(ME.genes) 
posterior.eigengene.df <- data.frame(Bayes.Factor = bf.eigengene, quantile.eigengene)

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
gene.info[3:ncol(gene.info)] <- signif(gene.info[3:ncol(gene.info)], 3)
gene.info.annot <- arrange(gene.info, Module, desc(kscaled)) %>% left_join(bm.table) %>% select(Symbol, Definition, Module:kscaled)

gene.module.membership <- data.frame(bicor(expr.data, ME.genes, maxPOutliers = 0.05)) %>% signif(3)
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.data)))
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

module.membership <- left_join(gene.info.annot, gene.module.membership) #%>% left_join(module.membership.pvalue, by = "Symbol")
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
modules.only <- select(module.membership, Module, Symbol)
enrichr.terms <- c("GO Biological Process", "GO Molecular Function", "GO Cellular Component", "KEGG", "Reactome", "GTEx Up", "GTEx Down") 
color.names <- unique(modules.only$Module) %>% sort

gaa.only <- filter(bf.annot, Bayes.Factor > 3)$Symbol
gaa.enrichr <- map(enrichr.terms, GetHyper, gaa.only, bf.annot$Symbol)
names(gaa.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, gaa.enrichr, "gaa")

gaa.gobiol.file <- "./enrichr/gaa/GO Biological Process.xlsx"
gaa.gobiol <- read.xlsx(gaa.gobiol.file) 
gaa.gobiol$Database <- "GO Biological Process"
gaa.gobiol.final <- slice(gaa.gobiol, c(1, 3, 6, 13, 15))
EnrichrPlot(gaa.gobiol.final, "gaa.gobiol")

brown.only <- filter(modules.only, Module == "brown")
brown.enrichr <- map(enrichr.terms, GetHyper, brown.only$Symbol, modules.only$Symbol)
names(brown.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, brown.enrichr, "brown")

brown.gobiol.file <- "./enrichr/brown/GO Biological Process.xlsx"
brown.gobiol <- read.xlsx(brown.gobiol.file) 
GetKappaCluster(file_path_sans_ext(brown.gobiol.file), brown.gobiol, modules.only$Symbol)
brown.gobiol$Database <- "GO BP"

brown.gomole.file <- "./enrichr/brown/GO Molecular Function.xlsx"
brown.gomole <- read.xlsx(brown.gomole.file)
GetKappaCluster(file_path_sans_ext(brown.gomole.file), brown.gomole, modules.only$Symbol)
brown.gomole$Database <- "GO MF"

brown.reactome.file <- "./enrichr/brown/Reactome.xlsx"
brown.reactome <- read.xlsx(brown.reactome.file)
GetKappaCluster(file_path_sans_ext(brown.reactome.file), brown.reactome, modules.only$Symbol)
brown.reactome$Database <- "Reactome"

brown.kegg.file <- "./enrichr/brown/KEGG.xlsx"
brown.kegg <- read.xlsx(brown.kegg.file)
brown.kegg$Database <- "KEGG"

brown.gobiol.final <- slice(brown.gobiol, c(10, 6, 5, 11, 93, 109))
brown.gomole.final <- slice(brown.gomole, 15)
brown.reactome.final <- slice(brown.reactome, c(4, 25))

brown.enrichr.final <- rbind(brown.gobiol.final, brown.gomole.final, brown.reactome.final)
EnrichrPlot(brown.enrichr, "brown", plot.height = 4, plot.width = 6)

#PPI for whole module
all.ppi <- GetPPI(brown.only$Symbol)
all.graph <- graph_from_edgelist(as.matrix(all.ppi))
brown.incident <- map(V(all.graph), incident, graph = all.graph) %>% map_int(length)

#Get genes in GO categories
brown.genes <- map(brown.enrichr$Genes, str_split, ",") %>% map(extract2, 1) 
brown.symbols <- reduce(brown.genes, c) %>% unique 

#Get PPI in GO genes
brown.ppi <-  map(brown.genes, GetPPI)
PlotPPI(adjacency.expr, brown.enrichr$Genes[[6]], brown.ppi[[6]], "reactome_insulin", plot.width = 20, plot.height = 20, vertex.size = 15)
PlotPPI(adjacency.expr, brown.enrichr$Genes[[7]], brown.ppi[[7]], "reactome_axon", plot.width = 20, plot.height = 20, vertex.size = 15)

#Get hub genes in GO genes
brown.kscaled <- map(brown.genes, GetKscaled, module.membership)
brown.kscaled.df <- reduce(brown.kscaled, rbind) %>% distinct
brown.kscaled.df %<>% slice(match(brown.symbols, brown.kscaled.df$Symbol))

#Plot shared genes in GO categories
brown.intersect <- fromList(brown.genes) %>%  data.frame 
brown.weighted <- apply(brown.intersect, 2, multiply_by, brown.kscaled.df$kscaled) %>% t %>% data.frame
colnames(brown.weighted) <- brown.symbols
brown.weighted$Term <- str_replace_all(brown.enrichr$Term, "\\ \\(GO.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% str_replace_all(",.*$", "")  
brown.weighted$Format.Name <- str_c(brown.enrichr$Database, ": ", brown.weighted$Term)
brown.weighted$Format.Name %<>% factor(levels = brown.weighted$Format.Name)
brown.plot <- gather(brown.weighted, Gene, Kscaled, -Format.Name, -Term)
brown.plot$Gene %<>% factor(levels = brown.symbols)

CairoPDF("enrichr.heatmap", width = 8, height = 30)
p <- ggplot(brown.plot, aes(Format.Name, Gene, fill = Kscaled)) + geom_raster() + theme_bw()
p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(panel.background = element_blank(), panel.border = element_blank())
p <- p + theme(plot.background = element_blank())
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
print(p)
dev.off()

elasticnet.table <- read.xlsx("../GAA_regression/elasticnet.table.xlsx") %>% select(-Definition)
elasticnet.gaa <- filter(elasticnet.table, Coefficient > 0 | Coefficient < 0)
gaa.combined <- left_join(elasticnet.gaa, bf.posterior.df)
gaa.combined.filter <- filter(gaa.combined, Bayes.Factor > 3)
