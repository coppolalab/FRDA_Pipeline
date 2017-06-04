#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()

library(limma)
library(lumi)
library(lumiHumanAll.db)
library(annotate)
library(biomaRt)
library(BayesFactor)
library(sva)

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
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(de.table), rule = ">0.5", style = sig.pvalues)
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
    CairoPDF(str_c(gene.symbol, ".pdf"), width = 6, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

GetGAABF <- function(gene.vector, model.design) {
    set.seed(12345)
    final.df <- mutate(model.design, Gene = gene.vector)
    bf.all <- regressionBF(Gene ~ GAA1.scaled, final.df)
    bf.all
}

GetDurationBF <- function(gene.vector, model.design) {
    set.seed(12345)
    final.df <- mutate(model.design, Gene = gene.vector)
    bf.all <- regressionBF(Gene ~ Duration.scaled, final.df)
    bf.all
}

GetFDSBF <- function(gene.vector, model.design) {
    set.seed(12345)
    final.df <- mutate(model.design, Gene = gene.vector)
    bf.all <- regressionBF(Gene ~ FDS.scaled, final.df)
    bf.all
}

EigengeneBF <- function(dataset, trait.df) {
    set.seed(12345)
    lm.df <- data.frame(Eigengene = dataset, Trait = trait.df)
    module.trait.cor <- regressionBF(Eigengene ~ Trait, lm.df)
    return(module.trait.cor)
}

#GetLM <- function(gene.vector, model.design) {
    #final.df <- mutate(model.design, Gene = gene.vector)
    #lm.all <- lmBF(Gene ~ GAA1 + Age + Sex + RIN , final.df)
    #lm.nogaa <- lmBF(Gene ~ Age + Sex + RIN, final.df)
    #c(lm.all, lm.nogaa)
#}

#GeneCorrelation <- function(dataset, trait.df) {
    #module.trait.cor <- bicor(dataset, trait.df, maxPOutliers = 0.05)
    #module.trait.cor
#}

#EigengeneCorrelation <- function(dataset, trait.df) {
    #module.trait.cor <- bicor(dataset, trait.df, maxPOutliers = 0.05)
    #module.trait.cor.pval <- corPvalueStudent(module.trait.cor, length(dataset))
    #c(module.trait.cor, module.trait.cor.pval)
#}

#EigengeneBayes <- function(dataset, trait.df) {
    #module.trait.cor <- bayes.cor.test(dataset, trait.df)
    ##module.trait.cor.pval <- corPvalueStudent(module.trait.cor, length(dataset))
    ##c(module.trait.cor, module.trait.cor.pval)
    #return(module.trait.cor)
#}

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

    p <- ggplot(gene.df, aes_string(x = "GAA1", y = "Expression")) + geom_jitter(col = color) + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("FDS Score") + ylab("Eigengene Value")
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

EnrichrPlot <- function(enrichr.df, filename, plot.height = 5, plot.width = 8, maintitle = "auto", color = "brown") {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% 
        map(unlist) %>% 
        map_int(length)
    #enrichr.df$Log.Bayes.Factor <- log10(enrichr.df$Bayes.Factor)
    enrichr.df$Term %<>% str_replace_all("\\ \\(GO.*$", "") %>% 
        str_replace_all("\\_Homo.*$", "") %>% 
        str_replace_all(", [a-z].*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.Bayes.Factor)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Bayes.Factor)  

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Log.Bayes.Factor)) + 
        geom_bar(stat = "identity", fill = color) + 
        coord_flip() + 
        theme_bw() + 
        theme(legend.position = "none", 
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1),
              plot.background = element_blank(), 
              plot.title = element_text(hjust = 0.5)) +
        ylab(expression(paste(Log[10], ' Bayes Factor')))
    if (maintitle == "auto"){
        p <- p + ggtitle(str_c(capitalize(filename), " Module"))
    } else {
        p <- p + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(maintitle)
    }
    CairoPDF(str_c(filename, ".enrichr"), height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

BayesPlot <- function(siggene, filename, threshold, plot.title, posterior.column = "Posterior", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Posterior Probability") {
    siggene$Significant <- siggene[[posterior.column]] > threshold 
    p <- ggplot(siggene, aes_string(x = log.column, y = posterior.column)) + 
        geom_point(aes(color = Significant)) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              legend.position = "none", 
              plot.background = element_blank(), 
              panel.border = element_rect(size = 1, color = "black"), 
              plot.title = element_text(hjust = 0.5)) + 
        xlab(xlabel) + ylab(ylabel) + 
        ggtitle(plot.title) + 
        scale_color_manual(values = c("orange", "darkgreen"), 
                           name = "Significant", 
                           labels = c("Yes", "No"))
    CairoPDF(filename, width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

GetKscaled <- function(gene.list, module.membership) {
    filter(module.membership, is.element(Symbol, gene.list)) %>% select(Symbol, kscaled) %>% arrange(desc(kscaled))
}

CollinearityScatterplot <- function(pdata.df, x.variable, y.variable, plot.prefix, plot.width = 6, plot.height = 6) {
    scatterplot <- ggplot(pdata.df, aes_string(x = x.variable, y = y.variable)) + 
        geom_point() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank()) +
        theme(panel.grid.minor = element_blank()) + 
        theme(plot.background = element_blank()) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(str_c(x.variable, " vs. ", y.variable))

    filename <- str_c(plot.prefix, x.variable, "vs.", y.variable, "scatterplot.pdf", sep = "_")
    CairoPDF(filename, height = plot.width, width = plot.height, bg = "transparent")
    plot(scatterplot)
    dev.off()
}

CollinearityBoxplot <- function(pdata.df, x.variable, y.variable, plot.prefix, plot.width = 6, plot.height = 6) {
    scatterplot <- ggplot(pdata.df, aes_string(x = x.variable, y = y.variable)) + 
        geom_boxplot() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank()) +
        theme(panel.grid.minor = element_blank()) + 
        theme(plot.background = element_blank()) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(str_c(x.variable, " vs. ", y.variable))

    filename <- str_c(plot.prefix, x.variable, "vs.", y.variable, "boxplot.pdf", sep = "_")
    CairoPDF(filename, height = plot.width, width = plot.height, bg = "transparent")
    plot(scatterplot)
    dev.off()
}

DensityPlot <- function(pdata.df, variable.name) {
    p <- ggplot(pdata.df, aes_string(variable.name)) + 
        geom_density() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        xlab(variable.name)
     
    CairoPDF(str_c(variable.name, "_density"), height = 5, width = 9, bg = "transparent")
    print(p)
    dev.off()
}

HistogramPlot <- function(pdata.df, variable.name) {
    p <- ggplot(pdata.df, aes_string(variable.name)) + 
        stat_count() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        xlab(variable.name)
     
    CairoPDF(str_c(variable.name, "_density"), height = 5, width = 9, bg = "transparent")
    print(p)
    dev.off()
}

GetCollinearity <- function(pdata.df) {
    age.gaa <- regressionBF(GAA1 ~ Age, data = pdata.df) %>% extractBF %>% extract2("bf")
    sex.gaa <- anovaBF(GAA1 ~ Sex, pdata.df) %>% extractBF %>% extract2("bf")
    site.gaa <- anovaBF(GAA1 ~ Site, pdata.df) %>% extractBF %>% extract2("bf")
    batch.gaa <- anovaBF(GAA1 ~ Batch, pdata.df) %>% extractBF %>% extract2("bf")
    RIN.gaa <- regressionBF(GAA1 ~ RIN, pdata.df) %>% extractBF %>% extract2("bf")
    FDS.gaa <- regressionBF(GAA1 ~ FDS, pdata.df) %>% extractBF %>% extract2("bf")
    duration.gaa <- regressionBF(GAA1 ~ Duration, pdata.df) %>% extractBF %>% extract2("bf")

    age.fds <- regressionBF(FDS ~ Age, data = pdata.df) %>% extractBF %>% extract2("bf")
    duration.fds <- regressionBF(FDS ~ Duration, data = pdata.df) %>% extractBF %>% extract2("bf")

    duration.age <- regressionBF(Duration ~ Age, data = pdata.df) %>% extractBF %>% extract2("bf")

    comparison.vector <- c("Age x GAA1", "Sex x GAA1", "Site x GAA1", "Batch x GAA1", 
                           "RIN x GAA1", "FDS x GAA1", "Duration x GAA1", 
                           "Age x FDS", "Duration x FDS", "Duration x Age")

    BF.vector <- c(age.gaa, sex.gaa, site.gaa, batch.gaa, RIN.gaa, FDS.gaa, duration.gaa, 
                   age.fds, duration.fds, duration.age) %>% log10 %>% signif(3)
    bf.tibble <- tibble(Comparison = comparison.vector, BF = BF.vector)
    bf.tibble
}

ExtractBF <- function(bf.list, gene.names) {
    bf.extract <- map(bf.list, extractBF) %>% 
        map(extract2, "bf") %>% 
        map_dbl(reduce, divide_by) %>% 
        map_dbl(log10) %>% 
        signif(3)
    bf.df <- tibble(Symbol = gene.names, Log.Bayes.Factor = bf.extract) %>% arrange(desc(Log.Bayes.Factor))
    bf.df
}

source("../../code/common_functions.R")

lumi.import <- ReadRDSgz(file = "../baseline_lumi/save/lumi.baseline.rda")
patient.pheno <- pData(lumi.import) %>% as_tibble %>% filter(Status == "Patient")
write.xlsx(patient.pheno, "patient_pheno.xlsx")
lumi.patient <- lumi.import[,lumi.import$Status == "Patient" & 
                            !is.na(lumi.import$GAA1) & 
                            !is.na(lumi.import$GAA2) &
                            !is.na(lumi.import$Duration) &
                            lumi.import$Duration >= 0 &
                            !is.na(lumi.import$FDS)]  
batch.colors <- data.frame(Batch = factor(1:19), Color = c("black","navy","blue","red","orange","cyan","tan","purple",
                                                           "lightcyan","lightyellow","darkseagreen","brown","salmon",
                                                           "gold4","pink","green", "blue4", "red4", "green4")) #Assign batch colors
lumi.patient$Batch.Color <- left_join(select(pData(lumi.patient), Batch), batch.colors) %>% 
    select(Color) %>% unlist #add variable for batch color

#gen.boxplot("patient_intensity_raw.jpg", lumi.vst, "VST Transformed Intensity", "Intensity")
Histogram("patient_histogram", lumi.patient)
mds.patient <- exprs(lumi.patient) %>% t %>% 
    dist(method = "manhattan") %>% 
    cmdscale(eig = TRUE) #get first two principle components
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
mds.norm <- exprs(lumi.patient.annot) %>% t %>% 
    dist(method = "manhattan") %>% 
    cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_batch_norm", mds.norm, pData(lumi.patient.annot), as.character(batch.colors$Color), "Batch") #label PCs by batch

#Use ComBat for batch effect correction
lumi.patient.annot$Sex %<>% droplevels
model.combat <- model.matrix(~ GAA1 + Sex + Age + RIN, pData(lumi.patient.annot)) %>% data.frame

expr.combat <- ComBat(dat = exprs(lumi.patient.annot), batch = factor(lumi.patient.annot$Site))
lumi.combat <- lumi.patient.annot
exprs(lumi.combat) <- expr.combat

lumi.combat <- lumi.combat[,lumi.combat$Batch != "8" & lumi.combat$Batch != "9"]
expr.combat <- ComBat(dat = exprs(lumi.combat), batch = factor(lumi.combat$Batch))
exprs(lumi.combat) <- expr.combat
SaveRDSgz(lumi.combat, file = "./save/lumi.combat.rda")

#Regenerate plots
#gen.boxplot("patient_intensity_combat.jpg", lumi.combat, "VST Transformed Intensity (Batch Corrected)", "Intensity")
Histogram("patient_histogram_combat", lumi.combat)
mds.patient.combat <- exprs(lumi.combat) %>% t %>% 
    dist(method = "manhattan") %>% 
    cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_batch_combat", mds.patient.combat, pData(lumi.combat), as.character(batch.colors$Color), "Batch") #label PCs by batch

#Remove outliers
combat.connectivity <- ConnectivityPlot("patient_connectivity_combat", lumi.combat, "")
connectivity.outlier <- combat.connectivity[combat.connectivity < -2]
outlier.names <- names(connectivity.outlier) %>% is.element(el = sampleNames(lumi.combat))

lumi.rmout <- lumi.combat[,!outlier.names & lumi.combat$GAA1 < 2000]
lumi.rmout$Site %<>% factor
lumi.rmout$Batch %<>% factor

#Regenerate plots
#gen.boxplot("patient_intensity_rmout.jpg", lumi.rmout, "VST Transformed Intensity (Outliers Removed)", "Intensity")
Histogram("patient_histogram_rmout", lumi.rmout)
mds.patient.rmout <- exprs(lumi.rmout) %>% t %>% 
    dist(method = "manhattan") %>% 
    cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_batch_rmout", mds.patient.rmout, pData(lumi.rmout), as.character(batch.colors$Color), "Batch") #label PCs by batch

#Check for collinearity 
#Subset
#tertiles of FDS - 1-2.5, 3-4.5, 5-6
#tertiles of Duration - 0-7, 8-14, 15+

#DensityPlot(pData(lumi.rmout), "Duration")

#rm(short.duration)
#short.duration <- filter(pData(lumi.rmout), Duration < 8)
#rm(medium.duration)
#medium.duration <- filter(pData(lumi.rmout), Duration < 15 & Duration >= 8)
#rm(long.duration)
#long.duration <- filter(pData(lumi.rmout), Duration >= 15)

#DensityPlot(short.duration, "Duration")
#DensityPlot(short.duration, "GAA1")
#DensityPlot(short.duration, "Age")

#rm(short.duration.bf)
#short.duration.bf <- GetCollinearity(short.duration)

#rm(medium.duration.bf)
#medium.duration.bf <- GetCollinearity(medium.duration)

#rm(long.duration.bf)
#long.duration.bf <- GetCollinearity(long.duration)

#rm(all.bf)
#all.bf <- GetCollinearity(pData(lumi.rmout))

#CollinearityScatterplot(pData(lumi.rmout), "GAA1", "Age", "all_patients")
#CollinearityScatterplot(short.duration, "GAA1", "Age", "short_duration")
#CollinearityScatterplot(medium.duration, "GAA1", "Age", "medium_duration")
#CollinearityScatterplot(long.duration, "GAA1", "Age", "long_duration")

#CollinearityBoxplot(short.duration, "Batch", "GAA1", "short_duration")
#CollinearityBoxplot(medium.duration, "Batch", "GAA1", "medium_duration")
#CollinearityBoxplot(long.duration, "Batch", "GAA1", "long_duration")

#HistogramPlot(fds.known, "FDS")
#FDS
#fds.known <- filter(pData(lumi.rmout), !is.na(FDS))
#low.fds <- filter(fds.known, FDS <= 4)
##medium.fds <- filter(fds.known, FDS < 5 & FDS >= 3)
#high.fds <- filter(fds.known, FDS > 4)

#low.fds.bf <- GetCollinearity(low.fds)
##medium.fds.bf <- GetCollinearity(medium.fds)
#high.fds.bf <- GetCollinearity(high.fds)
#all.fds.bf <- GetCollinearity(fds.known)

#Collapse by gene symbol
expr.collapse <- collapseRows(exprs(lumi.rmout), getSYMBOL(featureNames(lumi.rmout), 'lumiHumanAll.db'), 
                              rownames(lumi.rmout))$datETcollapsed
export.lumi <- ExpressionSet(assayData = expr.collapse, phenoData = phenoData(lumi.rmout))
export.lumi$GAA1.scaled <- scale(export.lumi$GAA1) %>% as.vector
export.lumi$FDS.scaled <- scale(export.lumi$FDS) %>% as.vector
export.lumi$Duration.scaled <- scale(export.lumi$Duration) %>% as.vector
SaveRDSgz(export.lumi, "./save/export.lumi.rda")

#Annotate top table
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(featureNames(export.lumi)), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")

all.model <- model.matrix( ~ Age + RIN, pData(export.lumi))[,-1] %>% data.frame
all.rmcov <- removeBatchEffect(export.lumi, covariates = all.model) %>% t
all.rmcov.scaled <- scale(all.rmcov)
#all.rmcov <- empiricalBayesLM(t(exprs(export.lumi)), removedCovariates = all.model) %>%
    #extract2("adjustedData")

bf.fds.all <- apply(all.rmcov.scaled, 2, GetFDSBF, model.design = pData(export.lumi))
SaveRDSgz(bf.fds.all, file = "./save/bf.fds.all.rda")
bf.fds.df <- ExtractBF(bf.fds.all, featureNames(export.lumi))
bf.fds.posterior <- map(bf.fds.all, posterior, iterations = 10000)
SaveRDSgz(bf.fds.posterior, "./save/bf.fds.posterior.rda")
bf.fds.posterior.df <- map(bf.fds.posterior, magrittr::extract, TRUE, 1) %>% 
    map(quantile, c(0.025, 0.5, 0.975)) %>% 
    reduce(rbind) %>% as.tibble %>%
    set_colnames(c("CI_2.5", "Median", "CI_97.5")) %>%
    mutate(Symbol = featureNames(export.lumi))
bf.fds.final.df <- left_join(bf.fds.posterior.df, bf.fds.df) %>% 
    arrange(desc(Log.Bayes.Factor)) %>%
    select(Symbol, Log.Bayes.Factor, CI_2.5, Median, CI_97.5) #%>% left_join(lumi.fds.cor.df)
bf.fds.sig <- filter(bf.fds.final.df, Log.Bayes.Factor > 0.5) %>% arrange(desc(abs(Median)))
SaveRDSgz(bf.fds.posterior.df, "./save/bf.fds.posterior.df.rda")

bf.fds.annot <- left_join(bf.fds.final.df, bm.table) %>% 
    select(Symbol, Definition, Log.Bayes.Factor, Median, CI_2.5, CI_97.5)
BayesWorkbook(bf.fds.annot, "bf.fds.xlsx")

bf.duration.all <- apply(all.rmcov.scaled, 2, GetDurationBF, model.design = pData(export.lumi))
SaveRDSgz(bf.duration.all, file = "./save/bf.duration.all.rda")
bf.duration.df <- ExtractBF(bf.duration.all, featureNames(export.lumi))
bf.duration.posterior <- map(bf.duration.all, posterior, iterations = 10000)
SaveRDSgz(bf.duration.posterior, "./save/bf.duration.posterior.rda")
bf.duration.posterior.df <- map(bf.duration.posterior, magrittr::extract, TRUE, 1) %>% 
    map(quantile, c(0.025, 0.5, 0.975)) %>% 
    reduce(rbind) %>% as.tibble %>%
    set_colnames(c("CI_2.5", "Median", "CI_97.5")) %>%
    mutate(Symbol = featureNames(export.lumi))
bf.duration.final.df <- left_join(bf.duration.posterior.df, bf.duration.df) %>% 
    arrange(desc(Log.Bayes.Factor)) %>%
    select(Symbol, Log.Bayes.Factor, CI_2.5, Median, CI_97.5) #%>% left_join(lumi.duration.cor.df)
bf.duration.sig <- filter(bf.duration.final.df, Log.Bayes.Factor > 0.5) %>% arrange(desc(abs(Median)))
SaveRDSgz(bf.duration.posterior.df, "./save/bf.duration.posterior.df.rda")

bf.duration.annot <- left_join(bf.duration.final.df, bm.table) %>% select(Symbol, Definition, Log.Bayes.Factor, Median, CI_2.5, CI_97.5)
BayesWorkbook(bf.duration.annot, "bf.duration.xlsx")

bf.gaa.all <- apply(all.rmcov.scaled, 2, GetGAABF, model.design = pData(export.lumi))
SaveRDSgz(bf.gaa.all, file = "./save/bf.gaa.all.rda")
bf.gaa.df <- ExtractBF(bf.gaa.all, featureNames(export.lumi))
bf.gaa.posterior <- map(bf.gaa.all, posterior, iterations = 10000)
SaveRDSgz(bf.gaa.posterior, "./save/bf.gaa.posterior.rda")
bf.gaa.posterior.df <- map(bf.gaa.posterior, magrittr::extract, TRUE, 1) %>% 
    map(quantile, c(0.025, 0.5, 0.975)) %>% 
    reduce(rbind) %>% as.tibble %>%
    set_colnames(c("CI_2.5", "Median", "CI_97.5")) %>%
    mutate(Symbol = featureNames(export.lumi))
bf.gaa.final.df <- left_join(bf.gaa.posterior.df, bf.gaa.df) %>% 
    arrange(desc(Log.Bayes.Factor)) %>%
    select(Symbol, Log.Bayes.Factor, CI_2.5, Median, CI_97.5) #%>% left_join(lumi.gaa.cor.df)
bf.gaa.sig <- filter(bf.gaa.final.df, Log.Bayes.Factor > 0.5) %>% arrange(desc(abs(Median)))
SaveRDSgz(bf.gaa.posterior.df, "./save/bf.gaa.posterior.df.rda")

bf.gaa.annot <- left_join(bf.gaa.final.df, bm.table) %>% select(Symbol, Definition, Log.Bayes.Factor, Median, CI_2.5, CI_97.5)
BayesWorkbook(bf.gaa.annot, "bf.gaa.xlsx")
BayesPlot(bf.fds.final.df, "fds_uplot", 0.5, "FDS", "Log.Bayes.Factor", "Median", "Regression Coefficient", "Log Bayes Factor")

enrichr.terms <- c("GO Biological Process", "GO Molecular Function", "KEGG", "Reactome")
bf.fds.enrichr <- map(enrichr.terms, GetHyper, bf.fds.sig$Symbol, bf.fds.final.df$Symbol)
names(bf.fds.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, bf.fds.enrichr, "fds")

#Enrichr
fds.gobiol.file <- "./enrichr/fds/GO Biological Process.xlsx"
fds.gobiol <- read.xlsx(fds.gobiol.file) 
GetKappaCluster(file_path_sans_ext(fds.gobiol.file), fds.gobiol, modules.only$Symbol, "Log.Bayes.Factor")
fds.gobiol$Database <- "GO BP"
fds.gobiol.final <- slice(fds.gobiol, c(2, 23, 31, 58))

fds.gomole.file <- "./enrichr/fds/GO Molecular Function.xlsx"
fds.gomole <- read.xlsx(fds.gomole.file) 
GetKappaCluster(file_path_sans_ext(fds.gomole.file), fds.gomole, modules.only$Symbol, "Log.Bayes.Factor")
fds.gomole$Database <- "GO MF"
fds.gomole.final <- slice(fds.gomole, c(8, 12, 13))

fds.reactome.file <- "./enrichr/fds/Reactome.xlsx"
fds.reactome <- read.xlsx(fds.reactome.file) 
GetKappaCluster(file_path_sans_ext(fds.reactome.file), fds.reactome, modules.only$Symbol, "Log.Bayes.Factor")
fds.reactome$Database <- "Reactome"
fds.reactome.final <- slice(fds.reactome, c(2, 4, 13))

fds.enrichr.final <- rbind(fds.gobiol.final, fds.gomole.final, fds.reactome.final)
EnrichrPlot(fds.enrichr.final, "fds", plot.height = 4, plot.width = 6, maintitle = "FDS", "blue")
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

expr.data <- t(exprs(export.lumi))
sft <- pickSoftThreshold(expr.data, powerVector = powers, verbose = 5, 
                         networkType = "signed", corFnc = bicor, 
                         corOptions = list(maxPOutliers = 0.05))
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

softPower <- 6
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

cor.gaa <- map(ME.genes, EigengeneCorrelation, lumi.pdata$GAA1) %>% 
    reduce(rbind) %>% data.frame
colnames(cor.gaa) <- c("Correlation", "P.value")
cor.gaa$Adj.P.value <- p.adjust(cor.gaa$P.value, method = "fdr")
cor.gaa$Module <- colnames(ME.genes)

ME.genes.model <- model.matrix( ~ Age + RIN, pData(export.lumi))[,-1]
ME.genes.rmcov <- removeBatchEffect(t(ME.genes), covariates = ME.genes.model) %>% t
ME.genes.rmcov.scaled <- scale(ME.genes.rmcov)
bayes.gaa <- apply(ME.genes.rmcov.scaled, 2, EigengeneBF, lumi.pdata$FDS.scaled)
bf.eigengene <- map(bayes.gaa, extractBF) %>% map_dbl(extract2, "bf") %>% map_dbl(log10)
posterior.eigengene <- map(bayes.gaa, posterior, iterations = 10000)
quantile.eigengene <- map(posterior.eigengene, magrittr::extract, TRUE, 1) %>% 
    map(quantile, c(0.025, 0.5, 0.975)) %>% 
    reduce(rbind) %>% 
    set_colnames(c("CI_2.5", "Median", "CI_97.5")) %>% 
    data.frame
quantile.eigengene$Color <- colnames(ME.genes) 
posterior.eigengene.df <- data.frame(Log.Bayes.Factor = bf.eigengene, quantile.eigengene) %>%
    mutate_if(is.numeric, signif, digits = 3)

text.matrix.traits <- str_c(signif(cor.gaa$Correlation, 2), '\n(', signif(cor.gaa$Adj.P.value, 1), ')')
dim(text.matrix.traits) = dim(as.matrix(cor.gaa$Correlation))

CairoPDF("module_trait_relationships", width = 4, height = 10)
par(mar = c(8, 8, 3, 3))
labeledHeatmap(Matrix = as.matrix(cor.gaa$Correlation), xLabels = "GAA1", yLabels = colnames(ME.genes), ySymbols = colnames(ME.genes), yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix.traits, zlim = c(-1,1), main = "")
dev.off()

#Eigengene plots
EigengeneScatterplot(ME.genes, lumi.pdata$FDS, "green")
EigengeneScatterplot(ME.genes, lumi.pdata$FDS, "blue")
EigengeneScatterplot(ME.genes, lumi.pdata$FDS, "pink")
EigengeneScatterplot(ME.genes, lumi.pdata$FDS, "purple")

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
#module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.data)))
#colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
#module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

module.membership <- left_join(gene.info.annot, gene.module.membership) #%>% left_join(module.membership.pvalue, by = "Symbol")
ModuleWorkbook(module.membership, "module_membership.xlsx")

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#Enrichr
source("../../code/GO/enrichr.R")
modules.only <- select(module.membership, Module, Symbol)
enrichr.terms <- c("GO Biological Process", "GO Molecular Function", "KEGG", "Reactome") 
color.names <- unique(modules.only$Module) %>% sort

blue.only <- filter(modules.only, Module == "blue")
green.only <- filter(modules.only, Module == "green")
pink.only <- filter(modules.only, Module == "pink")
purple.only <- filter(modules.only, Module == "purple")

blue.enrichr <- map(enrichr.terms, GetHyper, blue.only$Symbol, modules.only$Symbol)
names(blue.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, blue.enrichr, "blue")

green.enrichr <- map(enrichr.terms, GetHyper, green.only$Symbol, module.only$Symbol)
names(green.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, green.enrichr, "green")

pink.enrichr <- map(enrichr.terms, GetHyper, pink.only$Symbol, modules.only$Symbol)
names(pink.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, pink.enrichr, "pink")

purple.enrichr <- map(enrichr.terms, GetHyper, purple.only$Symbol, modules.only$Symbol)
names(purple.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, purple.enrichr, "purple")

blue.gobiol.file <- "./enrichr/blue/GO Biological Process.xlsx"
blue.gobiol <- read.xlsx(blue.gobiol.file) 
GetKappaCluster(file_path_sans_ext(blue.gobiol.file), blue.gobiol, blue.only$Symbol, "Log.Bayes.Factor")
blue.gobiol$Database <- "GO BP"

blue.gomole.file <- "./enrichr/blue/GO Molecular Function.xlsx"
blue.gomole <- read.xlsx(blue.gomole.file)
GetKappaCluster(file_path_sans_ext(blue.gomole.file), blue.gomole, blue.only$Symbol, "Log.Bayes.Factor")
blue.gomole$Database <- "GO MF"

blue.reactome.file <- "./enrichr/blue/Reactome.xlsx"
blue.reactome <- read.xlsx(blue.reactome.file)
GetKappaCluster(file_path_sans_ext(blue.reactome.file), blue.reactome, blue.only$Symbol, "Log.Bayes.Factor")
blue.reactome$Database <- "Reactome"

blue.kegg.file <- "./enrichr/blue/KEGG.xlsx"
blue.kegg <- read.xlsx(blue.kegg.file)
blue.kegg$Database <- "KEGG"

blue.gobiol.final <- slice(blue.gobiol, c(4, 5, 27, 16))
blue.gomole.final <- slice(blue.gomole, 2)
blue.kegg.final <- slice(blue.kegg, 3)
blue.reactome.final <- slice(blue.reactome, c(21))

blue.enrichr.final <- rbind(blue.gobiol.final, blue.gomole.final, blue.kegg.final, blue.reactome.final)
EnrichrPlot(blue.enrichr.final, "blue", plot.height = 4, plot.width = 6)

green.gobiol.file <- "./enrichr/green/GO Biological Process.xlsx"
green.gobiol <- read.xlsx(green.gobiol.file) 
GetKappaCluster(file_path_sans_ext(green.gobiol.file), green.gobiol, green.only$Symbol, "Log.Bayes.Factor")
green.gobiol$Database <- "GO BP"

green.gomole.file <- "./enrichr/green/GO Molecular Function.xlsx"
green.gomole <- read.xlsx(green.gomole.file)
GetKappaCluster(file_path_sans_ext(green.gomole.file), green.gomole, green.only$Symbol, "Log.Bayes.Factor")
green.gomole$Database <- "GO MF"

green.reactome.file <- "./enrichr/green/Reactome.xlsx"
green.reactome <- read.xlsx(green.reactome.file)
GetKappaCluster(file_path_sans_ext(green.reactome.file), green.reactome, green.only$Symbol, "Log.Bayes.Factor")
green.reactome$Database <- "Reactome"

green.kegg.file <- "./enrichr/green/KEGG.xlsx"
green.kegg <- read.xlsx(green.kegg.file)
green.kegg$Database <- "KEGG"

green.gobiol.final <- slice(green.gobiol, c(15, 1, 16, 3, 53))
green.gomole.final <- slice(green.gomole, c(1, 5, 16))
green.kegg.final <- slice(green.kegg, 23)
green.reactome.final <- slice(green.reactome, c(2, 10, 20))

green.enrichr.final <- rbind(green.gobiol.final, green.gomole.final, green.kegg.final, green.reactome.final)
EnrichrPlot(green.enrichr.final, "green", plot.height = 4, plot.width = 6)

pink.gobiol.file <- "./enrichr/pink/GO Biological Process.xlsx"
pink.gobiol <- read.xlsx(pink.gobiol.file) 
GetKappaCluster(file_path_sans_ext(pink.gobiol.file), pink.gobiol, pink.only$Symbol, "Log.Bayes.Factor")
pink.gobiol$Database <- "GO BP"

pink.gomole.file <- "./enrichr/pink/GO Molecular Function.xlsx"
pink.gomole <- read.xlsx(pink.gomole.file)
GetKappaCluster(file_path_sans_ext(pink.gomole.file), pink.gomole, pink.only$Symbol, "Log.Bayes.Factor")
pink.gomole$Database <- "GO MF"

pink.reactome.file <- "./enrichr/pink/Reactome.xlsx"
pink.reactome <- read.xlsx(pink.reactome.file)
GetKappaCluster(file_path_sans_ext(pink.reactome.file), pink.reactome, pink.only$Symbol, "Log.Bayes.Factor")
pink.reactome$Database <- "Reactome"

pink.kegg.file <- "./enrichr/pink/KEGG.xlsx"
pink.kegg <- read.xlsx(pink.kegg.file)
GetKappaCluster(file_path_sans_ext(pink.kegg.file), pink.kegg, pink.only$Symbol, "Log.Bayes.Factor")
pink.kegg$Database <- "KEGG"

pink.gobiol.final <- slice(pink.gobiol, c(1, 6, 8, 11, 13, 17))
pink.gomole.final <- slice(pink.gomole, 2)
pink.reactome.final <- slice(pink.reactome, c(2, 11))
pink.kegg.final <- slice(pink.kegg, c(12, 5, 2, 13))

pink.enrichr.final <- rbind(pink.gobiol.final, pink.gomole.final, pink.reactome.final, pink.kegg.final)
EnrichrPlot(pink.enrichr.final, "pink", plot.height = 4, plot.width = 6)

purple.gobiol.file <- "./enrichr/purple/GO Biological Process.xlsx"
purple.gobiol <- read.xlsx(purple.gobiol.file) 
GetKappaCluster(file_path_sans_ext(purple.gobiol.file), purple.gobiol, purple.only$Symbol, "Log.Bayes.Factor")
purple.gobiol$Database <- "GO BP"

purple.gomole.file <- "./enrichr/purple/GO Molecular Function.xlsx"
purple.gomole <- read.xlsx(purple.gomole.file)
GetKappaCluster(file_path_sans_ext(purple.gomole.file), purple.gomole, purple.only$Symbol, "Log.Bayes.Factor")
purple.gomole$Database <- "GO MF"

purple.reactome.file <- "./enrichr/purple/Reactome.xlsx"
purple.reactome <- read.xlsx(purple.reactome.file)
GetKappaCluster(file_path_sans_ext(purple.reactome.file), purple.reactome, purple.only$Symbol, "Log.Bayes.Factor")
purple.reactome$Database <- "Reactome"

purple.kegg.file <- "./enrichr/purple/KEGG.xlsx"
purple.kegg <- read.xlsx(purple.kegg.file)
purple.kegg$Database <- "KEGG"

purple.gobiol.final <- slice(purple.gobiol, c(5, 9, 24))
purple.reactome.final <- slice(purple.reactome, c(2, 11, 6))

purple.enrichr.final <- rbind(purple.gobiol.final, purple.reactome.final)
EnrichrPlot(purple.enrichr.final, "purple", plot.height = 4, plot.width = 6)

elasticnet.table <- read.xlsx("../GAA_regression/elasticnet.table.xlsx") %>% select(-Definition)
elasticnet.gaa <- filter(elasticnet.table, Coefficient > 0 | Coefficient < 0)
gaa.combined <- left_join(elasticnet.gaa, bf.posterior.df)
gaa.combined.filter <- filter(gaa.combined, Log.Bayes.Factor > 0.5) %>% arrange(desc(abs(Median)))

#PPI for whole module
#all.ppi <- GetPPI(brown.only$Symbol)
#all.graph <- graph_from_edgelist(as.matrix(all.ppi))
#brown.incident <- map(V(all.graph), incident, graph = all.graph) %>% map_int(length)

##Get genes in GO categories
#brown.genes <- map(brown.enrichr$Genes, str_split, ",") %>% 
    #map(extract2, 1) 
#brown.symbols <- reduce(brown.genes, c) %>% unique 

##Get PPI in GO genes
#brown.ppi <-  map(brown.genes, GetPPI)
#PlotPPI(adjacency.expr, brown.enrichr$Genes[[6]], brown.ppi[[6]], "reactome_insulin", plot.width = 20, plot.height = 20, vertex.size = 15)
#PlotPPI(adjacency.expr, brown.enrichr$Genes[[7]], brown.ppi[[7]], "reactome_axon", plot.width = 20, plot.height = 20, vertex.size = 15)

##Get hub genes in GO genes
#brown.kscaled <- map(brown.genes, GetKscaled, module.membership)
#brown.kscaled.df <- reduce(brown.kscaled, rbind) %>% distinct
#brown.kscaled.df %<>% slice(match(brown.symbols, brown.kscaled.df$Symbol))

##Plot shared genes in GO categories
#brown.intersect <- fromList(brown.genes) %>%  data.frame 
#brown.weighted <- apply(brown.intersect, 2, multiply_by, brown.kscaled.df$kscaled) %>% t %>% data.frame
#colnames(brown.weighted) <- brown.symbols
#brown.weighted$Term <- str_replace_all(brown.enrichr$Term, "\\ \\(GO.*$", "") %>% 
    #str_replace_all("\\_Homo.*$", "") %>% 
    #str_replace_all(",.*$", "")  
#brown.weighted$Format.Name <- str_c(brown.enrichr$Database, ": ", brown.weighted$Term)
#brown.weighted$Format.Name %<>% factor(levels = brown.weighted$Format.Name)
#brown.plot <- gather(brown.weighted, Gene, Kscaled, -Format.Name, -Term)
#brown.plot$Gene %<>% factor(levels = brown.symbols)

#CairoPDF("enrichr.heatmap", width = 8, height = 30)
#p <- ggplot(brown.plot, aes(Format.Name, Gene, fill = Kscaled)) + 
    #geom_raster() + 
    #theme_bw() + 
    #theme(axis.title.x = element_blank(), 
          #axis.title.y = element_blank(),
          #axis.ticks.x = element_blank(), 
          #axis.ticks.y = element_blank(), 
          #axis.text.x = element_text(angle = 45, hjust = 1), 
          #panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(),
          #panel.background = element_blank(), 
          #panel.border = element_blank(),
          #plot.background = element_blank(),
          #legend.position = "none")
#print(p)
#dev.off()


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
#lumi.short.duration <- export.lumi[,export.lumi$Duration < 12]
#lumi.short.model <- model.matrix( ~ Age + Sex + RIN, pData(lumi.short.duration))[,-1]
#short.duration.rmcov <- empiricalBayesLM(t(exprs(lumi.short.duration)), removedCovariates = lumi.short.model) %>%
    #extract2("adjustedData")

#short.duration.rmcov.scaled <- scale(short.duration.rmcov)
#bf.gaa.short <- apply(short.duration.rmcov.scaled, 2, GetBF, model.design = pData(lumi.short.duration))
#SaveRDSgz(bf.gaa.short, "./save/bf.gaa.short.rda")

#lumi.short.cor <- bicor(short.duration.rmcov, lumi.short.duration$GAA1) %>% as.vector
#lumi.short.cor.df <- tibble(Symbol = featureNames(lumi.short.duration), 
                            #Correlation = lumi.short.cor, 
                            #P.Value = corPvalueStudent(lumi.short.cor, ncol(lumi.short.duration))) %>%
                     #arrange(desc(abs(Correlation)))
#lumi.short.cor.df$Adj.P.Value <- p.adjust(lumi.short.cor.df$P.Value, method = "fdr")

##tcl.robust.df <- mutate(pData(lumi.medium.duration), Gene = medium.duration.robust[,"TRAF2"])
#tcl.rmcov <- medium.duration.rmcov[,"TCL6"]
#tcl.full <- GetBF(medium.duration.rmcov.scaled[,"TCL6"], pData(lumi.medium.duration)) 

#lumi.medium.duration <- export.lumi[,export.lumi$Duration >= 9 & export.lumi$Duration < 15]
#medium.duration.robust <- exprs(lumi.medium.duration) %>% t %>% 
    #sweep(2, colMedians(medium.duration.expr)) %>% 
    #sweep(2, colMads(medium.duration.expr), divide_by)

#lumi.medium.model <- model.matrix( ~ Age + Sex + RIN, pData(lumi.medium.duration))[,-1]
#medium.duration.rmcov <- empiricalBayesLM(t(exprs(lumi.medium.duration)), removedCovariates = lumi.medium.model) %>%
    #extract2("adjustedData")
#medium.duration.rmcov.scaled <- scale(medium.duration.rmcov)

#lumi.medium.cor <- bicor(medium.duration.rmcov, lumi.medium.duration$GAA1, maxPOutliers = 0.05) %>% as.vector
#lumi.medium.cor.df <- tibble(Symbol = featureNames(lumi.medium.duration), 
                            #Correlation = lumi.medium.cor, 
                            #P.Value = corPvalueStudent(lumi.medium.cor, ncol(lumi.medium.duration))) %>%
                     #arrange(desc(abs(Correlation)))
#lumi.medium.cor.df$Adj.P.Value <- p.adjust(lumi.medium.cor.df$P.Value, method = "fdr")

#bf.gaa.medium <- apply(medium.duration.rmcov.scaled, 2, GetBF, model.design = pData(lumi.medium.duration))
#SaveRDSgz(bf.gaa.medium, "./save/bf.gaa.medium.rda")

#lumi.long.duration <- export.lumi[,export.lumi$Duration >= 12]
#lumi.long.model <- model.matrix( ~ Age + Sex + RIN, pData(lumi.long.duration))[,-1]
#long.duration.rmcov <- empiricalBayesLM(t(exprs(lumi.long.duration)), removedCovariates = lumi.long.model) %>%
    #extract2("adjustedData")
#long.duration.rmcov.scaled <- scale(long.duration.rmcov)

#bf.gaa.long <- apply(long.duration.rmcov.scaled, 2, GetBF, model.design = pData(lumi.long.duration))
#SaveRDSgz(bf.gaa.long, "./save/bf.gaa.long.rda")

#bf.short.df <- ExtractBF(bf.gaa.short, featureNames(export.lumi))
#bf.medium.df <- ExtractBF(bf.gaa.medium, featureNames(export.lumi))
#bf.long.df <- ExtractBF(bf.gaa.long, featureNames(export.lumi))

#bf.short.posterior <- map(bf.gaa.short, posterior, iterations = 10000)
#SaveRDSgz(bf.short.posterior, "./save/bf.short.posterior.rda")
#bf.short.posterior.df <- map(bf.short.posterior, magrittr::extract, TRUE, 1) %>% 
    #map(quantile, c(0.025, 0.5, 0.975)) %>% 
    #reduce(rbind) %>% as.tibble %>%
    #set_colnames(c("CI_2.5", "Median", "CI_97.5")) %>%
    #mutate(Symbol = featureNames(export.lumi))
#bf.short.final.df <- left_join(bf.short.posterior.df, bf.short.df) %>% 
    #arrange(desc(Log.Bayes.Factor)) %>%
    #select(Symbol, Log.Bayes.Factor, CI_2.5, Median, CI_97.5) #%>% left_join(lumi.short.cor.df)
#bf.short.sig <- filter(bf.short.final.df, Log.Bayes.Factor > 0.5) %>% arrange(desc(abs(Median)))
#SaveRDSgz(bf.short.posterior.df, "./save/bf.short.posterior.df.rda")

#bf.medium.posterior <- map(bf.gaa.medium, posterior, iterations = 10000)
#SaveRDSgz(bf.medium.posterior, "./save/bf.medium.posterior.rda")
#bf.medium.posterior.df <- map(bf.medium.posterior, magrittr::extract, TRUE, 1) %>% 
    #map(quantile, c(0.025, 0.5, 0.975)) %>% 
    #reduce(rbind) %>% as.tibble %>%
    #set_colnames(c("CI_2.5", "Median", "CI_97.5")) %>%
    #mutate(Symbol = featureNames(export.lumi))
#bf.medium.final.df <- left_join(bf.medium.posterior.df, bf.medium.df) %>% 
    #arrange(desc(Log.Bayes.Factor)) %>%
    #select(Symbol, Log.Bayes.Factor, CI_2.5, Median, CI_97.5) #%>% left_join(lumi.medium.cor.df)
#bf.medium.sig <- filter(bf.medium.final.df, Log.Bayes.Factor > 0.5) %>% arrange(desc(abs(Median)))
#SaveRDSgz(bf.medium.posterior.df, "./save/bf.medium.posterior.df.rda")

#bf.long.posterior <- map(bf.gaa.long, posterior, iterations = 10000)
#SaveRDSgz(bf.long.posterior, "./save/bf.long.posterior.rda")
#bf.long.posterior.df <- map(bf.long.posterior, magrittr::extract, TRUE, 1) %>% 
    #map(quantile, c(0.025, 0.5, 0.975)) %>% 
    #reduce(rbind) %>% as.tibble %>%
    #set_colnames(c("CI_2.5", "Median", "CI_97.5")) %>%
    #mutate(Symbol = featureNames(export.lumi))
#bf.long.final.df <- left_join(bf.long.posterior.df, bf.long.df) %>% 
    #arrange(desc(Log.Bayes.Factor)) %>%
    #select(Symbol, Log.Bayes.Factor, CI_2.5, Median, CI_97.5) #%>% left_join(lumi.long.cor.df)
#bf.long.sig <- filter(bf.long.final.df, Log.Bayes.Factor > 0.5) %>% arrange(desc(abs(Median)))
#SaveRDSgz(bf.long.posterior.df, "./save/bf.long.posterior.df.rda")

#GeneScatterplot(lumi.medium.duration, "ALDH1B1")
#GeneScatterplot(lumi.medium.duration, "TNFRSF19")

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
