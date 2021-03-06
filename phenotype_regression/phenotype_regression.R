library(limma)
library(lumi)
library(lumiHumanAll.db)
library(annotate)
library(biomaRt)
library(BayesFactor)
library(sva)
library(matrixStats)
library(WGCNA)

library(Cairo)
library(openxlsx)

library(tools)
library(R.utils)
library(broom)
library(stringr)
library(magrittr)
library(pryr)
library(forcats)
library(tidyverse)

#Boxplot function
BoxPlot <- function(filename, lumi.object, colorscheme, maintext, ylabtext) {
    expr.df <- exprs(lumi.object) %>% t %>% as_tibble
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Batch = lumi.object$Batch)
    dataset.gather <- gather(dataset.addvars, -Sample.Status, -Batch)

    p <- ggplot(dataset.gather, aes(x = Sample.Status, y = value, fill = factor(Batch))) + 
        geom_boxplot() + 
        theme_bw() + 
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5)) + 
        ggtitle(maintext) + ylab(ylabtext) + xlab("Sample")  + 
        scale_fill_manual(values = colorscheme)
    ggsave(filename = filename, plot = p, family = "Noto Sans", width = 20 , height = 8)
}

#Histogram
Histogram <- function(filename, lumi.object) {
    expr.df <- exprs(lumi.object) %>% t %>% as_tibble
    dataset.addvars <- mutate(expr.df, Sample.Name = sampleNames(lumi.object), Status = lumi.object$Status)
    dataset.gather <- gather(dataset.addvars, nuID, Expression, -Sample.Name, -Status)

    p <- ggplot(dataset.gather, aes(Expression, group = Sample.Name, col = factor(Status))) + 
        geom_density() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) + 
        ggtitle("Histogram of VST Expression") + 
        ylab("Density") + xlab("VST Expression") 

    CairoPDF(filename, height = 5, width = 9)
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

    connectivity.plot <- tibble(Slide.ID = names(connectivity.zscore), 
                                Z.score = connectivity.zscore, 
                                Sample.Num = 1:length(connectivity.zscore))

    p <- ggplot(connectivity.plot, aes(x = Sample.Num, y = Z.score, label = Slide.ID) ) + 
        geom_text(size = 4, colour = "red") + 
        geom_hline(aes(yintercept = -2)) + 
        geom_hline(yintercept = -3) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) + 
        xlab("Sample Number") + ylab("Z score") + 
        ggtitle(maintitle)

    CairoPDF(filename, width = 10, height = 10)
    print(p)
    dev.off()
    connectivity.zscore
}

#MDS function 
MDSPlot <- function(filename, dataset, targetset, colorscheme = "none", variablename) {
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$Slide.ID, factor(targetset[[variablename]]))
    colnames(target.data) <- c("Slide.ID", variablename)
    colnames(dataset.plot) <- c("Slide.ID", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)

    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + 
        geom_point() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) + 
        xlab("Component 1") + ylab("Component 2") + 
        ggtitle(variablename)

    if (colorscheme != "none") {
        p <- p + scale_color_manual(values = colorscheme) 
    }

    CairoPDF(file = filename, height = 7, width = 7)
    print(p)
    dev.off()
}

CollinearityScatterplot <- function(pdata.df, x.variable, y.variable, plot.width = 6, plot.height = 6) {
    scatterplot <- ggplot(pdata.df, aes_string(x = x.variable, y = y.variable)) + 
        geom_point() + 
        geom_smooth(method = "lm") +
        theme_bw() + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        ggtitle(str_c(x.variable, " vs. ", y.variable))

    filename <- str_c(x.variable, "vs.", y.variable, "scatterplot.pdf", sep = "_")
    CairoPDF(filename, height = plot.width, width = plot.height, bg = "transparent")
    plot(scatterplot)
    dev.off()
}

CollinearityViolinPlot <- function(pdata.df, x.variable, y.variable, plot.width = 5.5, plot.height = 5) {
    violinplot <- ggplot(pdata.df, aes_string(x = x.variable, y = y.variable, fill = x.variable)) + 
        geom_violin(scale = "width", trim = FALSE, draw_quantiles = c(0.05, 0.5, 0.95)) +
        theme_bw() + 
        theme(axis.title.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        ggtitle(str_c(x.variable, " vs. ", y.variable))

    filename <- str_c(x.variable, "vs.", y.variable, "boxplot.pdf", sep = "_")
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    plot(violinplot)
    dev.off()
}

CollinearityBoxplot <- function(pdata.df, x.variable, y.variable, plot.width = 5.5, plot.height = 5) {
    p <- ggplot(pdata.df, aes_string(x = x.variable, y = y.variable, fill = x.variable)) + 
        geom_boxplot() +
        theme_bw() + 
        theme(axis.title.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.position = "none") +
        ggtitle(str_c(x.variable, " vs. ", y.variable))

    filename <- str_c(x.variable, "vs.", y.variable, "boxplot.pdf", sep = "_")
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    plot(p)
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
        xlab(variable.name) + ylab("Number of Subjects")
     
    CairoPDF(str_c(variable.name, "_histogram"), height = 5, width = 9, bg = "transparent")
    print(p)
    dev.off()
}

GetCollinearity <- function(pdata.df) {
    age.FDS <- regressionBF(FDS ~ Age, data = pdata.df) %>% extractBF %>% extract2("bf")
    sex.FDS <- anovaBF(FDS ~ Sex, pdata.df) %>% extractBF %>% extract2("bf")
    site.FDS <- anovaBF(FDS ~ Site, pdata.df) %>% extractBF %>% extract2("bf")
    batch.FDS <- anovaBF(FDS ~ Batch, pdata.df) %>% extractBF %>% extract2("bf")
    RIN.FDS <- regressionBF(FDS ~ RIN, pdata.df) %>% extractBF %>% extract2("bf")
    duration.FDS <- regressionBF(FDS ~ Duration, pdata.df) %>% extractBF %>% extract2("bf")

    comparison.vector <- c("Age x FDS", "Sex x FDS", "Site x FDS", "Batch x FDS", 
                           "RIN x FDS", "Duration x FDS")

    BF.vector <- c(age.FDS, sex.FDS, site.FDS, batch.FDS, RIN.FDS, duration.FDS) %>% log10 %>% signif(3)
    bf.tibble <- tibble(Comparison = comparison.vector, BF = BF.vector)
    bf.tibble
}

CovarBayes <- function(gene.name, expr.matrix, trait.df) {
    gene.vector <- expr.matrix[gene.name,]
    trait.df <- data.frame(trait.df, Gene = gene.vector)
    set.seed(12345)
    trait.anova <- regressionBF(Gene ~ Age, data = trait.df) 
    set.seed(12345)
    trait.posterior <- posterior(trait.anova, iterations = 10000) 
    saveRDS(trait.posterior, str_c("./save/posterior/", gene.name), compress = FALSE)
}

FDSBayes <- function(gene.name, expr.matrix, trait.df) {
    gene.vector <- expr.matrix[gene.name,]
    trait.df <- data.frame(trait.df, Gene = gene.vector)
    set.seed(12345)
    trait.anova <- lmBF(Gene ~ FDS + Sex + RIN, data = trait.df) 
    notrait.anova <- lmBF(Gene ~ Sex + RIN, data = trait.df) 
    combined <- c(trait.anova, notrait.anova)
    combined
}

GAABayes <- function(gene.name, expr.matrix, trait.df) {
    gene.vector <- expr.matrix[gene.name,]
    trait.df <- data.frame(trait.df, Gene = gene.vector)
    set.seed(12345)
    trait.anova <- lmBF(Gene ~ GAA1 + Sex + RIN, data = trait.df) 
    notrait.anova <- lmBF(Gene ~ Sex + RIN, data = trait.df) 
    combined <- c(trait.anova, notrait.anova)
    combined
}

DurationBayes <- function(gene.name, expr.matrix, trait.df) {
    gene.vector <- expr.matrix[gene.name,]
    trait.df <- data.frame(trait.df, Gene = gene.vector)
    set.seed(12345)
    trait.anova <- lmBF(Gene ~ Duration + Sex + RIN, data = trait.df) 
    notrait.anova <- lmBF(Gene ~ Sex + RIN, data = trait.df) 
    combined <- c(trait.anova, notrait.anova)
    combined
}

GetPosteriorProbability <- function(model.pheno, variable, n.iterations) {
    variable.column <- str_c(variable, "-", variable)
    set.seed(12345)
    posterior.model <- posterior(model.pheno, index = 1, iterations = n.iterations)
    posterior.variable <- posterior.model[,variable.column]
    median.posterior <- median(posterior.variable)

    if (sign(median.posterior) > 0) {
        pp.variable <- length(which(posterior.variable > 0)) / n.iterations
    } else {
        pp.variable <- length(which(posterior.variable < 0)) / n.iterations
    }
    c(coef.variable = median.posterior, pp.variable = pp.variable)
}

BayesPlot <- function(gene.df, filename, threshold, plot.title, posterior.column = "Posterior", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Posterior Probability") {
    gene.df$Significant <- gene.df[[posterior.column]] > threshold 
    sig.df <- filter(gene.df, Significant)
    plot.title <- str_c(plot.title, " (", nrow(sig.df), "/", nrow(gene.df), " Transcripts)")

    p <- ggplot(gene.df, aes_string(x = log.column, y = posterior.column)) + 
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
        scale_color_manual(values = c("gray", "blue"))

    CairoPDF(filename, width = 5, height = 5, bg = "transparent")
    print(p)
    dev.off()
}

BayesWorkbook <- function(de.table, filename) { 
    bf.cols <- colnames(de.table) %>% str_detect("logBF") %>% which
    pp.cols <- colnames(de.table) %>% str_detect("pp") %>% which
    cor.cols <- colnames(de.table) %>% str_detect("coef") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = de.table)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = bf.cols, rows = 1:nrow(de.table), 
                          rule = ">0.5", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = pp.cols, rows = 1:nrow(de.table), 
                          rule = ">0.95", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = cor.cols, rows = 1:nrow(de.table), 
                          style = c("#63BE7B", "white", "red"), type = "colourScale")

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
    gene.df <- data.frame(FDS = lumi.object$FDS, Expression = gene.expr)

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
         xlab("FDS") + ylab("VST Normalized Expression") + 
         ggtitle(gene.symbol) 

    CairoPDF(str_c(gene.symbol, ".pdf"), width = 6, height = 6, bg = "transparent")
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

EnrichrPlot <- function(enrichr.df, enrichr.expr, prefix, plot.title, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% 
        map(unlist) %>% map_int(length)
    enrichr.updown <- map(enrichr.df$Genes, UpDown, enrichr.expr) %>% reduce(rbind)
    colnames(enrichr.updown) <- c("Down", "Up")

    enrichr.df <- cbind(enrichr.df, enrichr.updown)
    enrichr.df$Log.Up <- enrichr.df$Log.Bayes.Factor * enrichr.df$Up / enrichr.df$Gene.Count
    enrichr.df$Log.Down <- enrichr.df$Log.Bayes.Factor * enrichr.df$Down / enrichr.df$Gene.Count
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% 
        str_replace_all("\\_Homo.*$", "") #%>% tolower #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.Bayes.Factor)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Up, Log.Down) %>% 
        gather(Direction, Length, -Format.Name) 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Length, fill = Direction))  + 
        geom_bar(stat = "identity", color = "black")  + 
        scale_fill_discrete(name = "Direction", labels = c("Up", "Down")) + 
        coord_flip() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),  
              plot.background = element_blank(), 
              plot.title = element_text(hjust = 0.5),
              legend.background = element_blank(), 
              axis.title.y = element_blank(), 
              axis.ticks.y = element_blank(), 
              axis.text.y = element_text(size = 12),
              axis.line.x = element_line()) +
        ylab("logBF") + ggtitle(plot.title)

    CairoPDF(str_c(prefix, ".enrichr"), height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

UpDown <- function(filter.vector, enrichr.df) {
    split.vector <- str_split(filter.vector, ",")[[1]]
    enrichr.filter <- filter(enrichr.df, Symbol %in% split.vector)
    enrichr.vector <- c("Up" = length(which(sign(enrichr.filter$Median) == 1)), 
                        "Down" = length(which(sign(enrichr.filter$Median) == -1)))
    enrichr.vector
}

source("../../code/common_functions.R")

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

lumi.import <- ReadRDSgz(file = "../differential_expression/save/lumi.baseline.rda")
patient.pheno <- pData(lumi.import) %>% as_tibble %>% filter(Status == "Patient")
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

#Remove outliers
combat.connectivity <- ConnectivityPlot("patient_connectivity_combat", lumi.patient.annot, "")
connectivity.outlier <- combat.connectivity[combat.connectivity < -2]
outlier.names <- names(connectivity.outlier) %>% is.element(el = sampleNames(lumi.patient.annot))

lumi.rmout <- lumi.patient.annot[,!outlier.names & lumi.patient.annot$GAA1 < 2000]
lumi.rmout$Site %<>% factor
lumi.rmout$Batch %<>% factor
lumi.rmout$Sex %<>% droplevels

#Regenerate plots
#gen.boxplot("patient_intensity_rmout.jpg", lumi.rmout, "VST Transformed Intensity (Outliers Removed)", "Intensity")
Histogram("patient_histogram_rmout", lumi.rmout)
mds.patient.rmout <- exprs(lumi.rmout) %>% t %>% 
    dist(method = "manhattan") %>% 
    cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_batch_rmout", mds.patient.rmout, pData(lumi.rmout), as.character(batch.colors$Color), "Batch") #label PCs by batch

#Use ComBat for batch effect correction
expr.combat <- ComBat(dat = exprs(lumi.rmout), batch = factor(lumi.rmout$Site))
lumi.combat <- lumi.rmout
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

#Collapse by gene symbol
combat.collapse.expr <- collapseRows(exprs(lumi.combat), getSYMBOL(featureNames(lumi.combat), 'lumiHumanAll.db'), 
                              rownames(lumi.combat))$datETcollapsed
combat.collapse <- ExpressionSet(assayData = combat.collapse.expr, phenoData = phenoData(lumi.combat))
SaveRDSgz(combat.collapse, "./save/combat.collapse.rda")


#Check for collinearity 
collinearity <- GetCollinearity(pData(combat.collapse))
DensityPlot(pData(combat.collapse), "Duration")
DensityPlot(pData(combat.collapse), "GAA1")
HistogramPlot(pData(combat.collapse), "FDS")

CollinearityScatterplot(pData(combat.collapse), "GAA1", "Age")
CollinearityScatterplot(pData(combat.collapse), "Duration", "Age")
CollinearityScatterplot(pData(combat.collapse), "FDS", "Age")
CollinearityScatterplot(pData(combat.collapse[,combat.collapse$RIN > 5]), "FDS", "RIN")

CollinearityViolinPlot(pData(combat.collapse), "Sex", "FDS")
CollinearityViolinPlot(pData(combat.collapse), "Site", "FDS")
CollinearityBoxplot(pData(combat.collapse), "Batch", "FDS")

#Remove effects of covariates
catch <- mclapply(featureNames(combat.collapse), CovarBayes, exprs(combat.collapse), select(pData(combat.collapse), Age), mc.cores = 8)
posterior.covar <- list.files("./save/posterior", full.names = TRUE) %>% map(readRDS)
SaveRDSgz(posterior.covar, "./save/posterior.covar.rda")
covar.coefs <- map(posterior.covar, colMedians) %>% reduce(rbind)
colnames(covar.coefs) <- colnames(posterior.covar[[1]])
SaveRDSgz(covar.coefs, "./save/covar.coefs.rda")

covar.remove <- model.matrix(~ Age, pData(combat.collapse))[,-1] %>% t
rmcov.collapse.expr <- exprs(combat.collapse) - (covar.coefs[,c("Age")] %*% covar.remove)
rmcov.lumi <- combat.collapse #Make a copy of lumi object
exprs(rmcov.lumi) <- rmcov.collapse.expr #Transfer cleaned expression values into new lumi object
SaveRDSgz(rmcov.lumi, "./save/rmcov.lumi.rda")

#Get FDS regression
model.fds <- mclapply(featureNames(rmcov.lumi), FDSBayes, exprs(rmcov.lumi), 
                  select(pData(rmcov.lumi), FDS, Sex, RIN), mc.cores = 8)
names(model.fds) <- featureNames(rmcov.lumi)
SaveRDSgz(model.fds, "./save/model.fds.rda")
bf.fds <- map(model.fds, extractBF) %>%
    map(extract2, "bf") %>% 
    map_dbl(reduce, divide_by) %>% log10
bf.fds.df <- tibble(Symbol = featureNames(rmcov.lumi), logBF = bf.fds)

fds.pp <- mclapply(model.fds, GetPosteriorProbability, "FDS", 10000, mc.cores = 8) %>%
    reduce(rbind) %>% as_tibble
colnames(fds.pp) <- c("coef.FS", "pp.FS")
fds.pp$Symbol <- featureNames(rmcov.lumi)

#Annotate top table
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(featureNames(rmcov.lumi)), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")
bm.table %<>% filter(!duplicated(Symbol))

bf.fds.annot <- left_join(bf.fds.df, fds.pp) %>% left_join(bm.table) %>%
    select(Symbol, Definition, logBF, coef.FS, pp.FS) %>%
    mutate_at(vars(logBF, coef.FS), signif, digits = 3) %>%
    arrange(desc(logBF))
SaveRDSgz(bf.fds.annot, "./save/bf.fds.annot.rda")
bf.fds.sig <- filter(bf.fds.annot, logBF > 0.5 & pp.FS > 0.95)
bf.fds.fdr <- mean(1 - bf.fds.sig$pp.FS)
BayesWorkbook(bf.fds.annot, "bf.fds.xlsx")
BayesPlot(bf.fds.annot, "fds_uplot", 0.5, "Functional Stage", "logBF", "coef.FS", 
          "Regression Coefficient", "logBF")

#Get GAA regression
catch <- mclapply(featureNames(rmcov.lumi), GAABayes, exprs(rmcov.lumi), 
                  select(pData(rmcov.lumi), GAA1, Sex, RIN), mc.cores = 8)
model.gaa <- list.files("./save/model", full.names = TRUE) %>% map(readRDS) 
SaveRDSgz(model.gaa, "./save/model.gaa.rda")
names(model.gaa) <- featureNames(rmcov.lumi)
bf.gaa <- map(model.gaa, extractBF) %>%
    map(extract2, "bf") %>% 
    map_dbl(reduce, divide_by) %>% log10
bf.gaa.df <- tibble(Symbol = featureNames(rmcov.lumi), logBF = bf.gaa)

gaa.pp <- mclapply(model.gaa, GetPosteriorProbability, "GAA1", 10000, mc.cores = 8) %>%
    reduce(rbind) %>% as_tibble
colnames(gaa.pp) <- c("coef.GAA", "pp.GAA")
gaa.pp$Symbol <- featureNames(rmcov.lumi)

bf.gaa.annot <- left_join(bf.gaa.df, gaa.pp) %>% left_join(bm.table) %>%
    select(Symbol, Definition, logBF, coef.GAA, pp.GAA) %>%
    arrange(desc(logBF))
bf.gaa.sig <- filter(bf.gaa.annot, logBF > 0.5 & pp.GAA > 0.95)
bf.gaa.fdr <- mean(1 - bf.gaa.sig$pp.GAA)
SaveRDSgz(bf.gaa.annot, "./save/bf.gaa.annot.rda")

BayesWorkbook(bf.gaa.annot, "bf.gaa.xlsx")
BayesPlot(bf.gaa.final.df, "gaa_uplot", 0.5, "GAA1", "Log.Bayes.Factor", "Median", 
          "Regression Coefficient", "Log Bayes Factor")

#Get duration regression
catch <- mclapply(featureNames(rmcov.lumi), DurationBayes, exprs(rmcov.lumi), select(pData(rmcov.lumi), Duration, Sex, RIN), mc.cores = 8)
model.duration <- list.files("./save/model", full.names = TRUE) %>% map(readRDS) 
SaveRDSgz(model.duration, "./save/model.duration.rda")
bf.duration <- map(model.duration, extractBF) %>%
    map(extract2, "bf") %>% 
    map_dbl(reduce, divide_by) %>% log10
bf.duration.df <- tibble(Symbol = featureNames(rmcov.lumi), logBF = bf.duration)

duration.pp <- mclapply(model.duration, GetPosteriorProbability, "Duration", 10000, mc.cores = 8) %>%
    reduce(rbind) %>% as_tibble
colnames(duration.pp) <- c("coef.Duration", "pp.Duration")
duration.pp$Symbol <- featureNames(rmcov.lumi)

bf.duration.annot <- left_join(bf.duration.df, duration.pp) %>% left_join(bm.table) %>%
    select(Symbol, Definition, logBF, coef.Duration, pp.Duration) %>%
    arrange(desc(logBF))
bf.duration.sig <- filter(bf.duration.annot, logBF > 0.5 & pp.Duration > 0.95)
bf.duration.fdr <- mean(1 - bf.duration.sig$pp.Duration)
SaveRDSgz(bf.duration.annot, "./save/bf.duration.annot.rda")

BayesWorkbook(bf.duration.annot, "bf.duration.xlsx")
BayesPlot(bf.duration.final.df, "duration_uplot", 0.5, "Duration", "Log.Bayes.Factor", "Median", "Regression Coefficient", "Log Bayes Factor")

module.membership <- read.xlsx("../WGCNA_phenotypes/module_membership.xlsx")
module.red <- filter(module.membership, Module == "red") 
fds.red <- filter(bf.fds.annot, Log.Bayes.Factor > 0.5) %>%
    inner_join(module.red) %>% arrange(desc(kscaled))
module.yellow <- filter(module.membership, Module == "yellow") 
fds.yellow <- filter(bf.fds.annot, Log.Bayes.Factor > 0.5) %>%
    inner_join(module.yellow) %>% arrange(desc(kscaled))
module.magenta <- filter(module.membership, Module == "magenta") 
fds.magenta <- filter(bf.fds.annot, Log.Bayes.Factor > 0.5) %>%
    inner_join(module.magenta) %>% arrange(desc(kscaled))

elasticnet.table <- read.xlsx("../GAA_regression/elasticnet.table.xlsx") %>% select(-Definition)
elasticnet.fds <- filter(elasticnet.table, Coefficient != 0)
fds.combined <- left_join(elasticnet.fds, bf.fds.final.df)
fds.combined.filter <- filter(fds.combined, Log.Bayes.Factor > 0.5) %>% arrange(desc(abs(Median)))

#Enrichr
source("../../code/GO/enrichr.R")

enrichr.terms <- c("GO Biological Process", "GO Molecular Function", "KEGG", "Reactome")
bf.fds.enrichr <- map(enrichr.terms, GetHyper, bf.fds.sig$Symbol, bf.fds.final.df$Symbol)
names(bf.fds.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, bf.fds.enrichr, "fds")

fds.gobiol.file <- "./enrichr/fds/GO Biological Process.xlsx"
fds.gobiol <- read.xlsx(fds.gobiol.file) 
fds.gobiol$Database <- "GO BP"
GetKappaCluster(file_path_sans_ext(fds.gobiol.file), fds.gobiol, bf.fds.sig$Symbol, "Log.Bayes.Factor")

fds.gomole.file <- "./enrichr/fds/GO Molecular Function.xlsx"
fds.gomole <- read.xlsx(fds.gomole.file) 
fds.gomole$Database <- "GO MF"
GetKappaCluster(file_path_sans_ext(fds.gomole.file), fds.gomole, bf.fds.sig$Symbol, "Log.Bayes.Factor")

fds.reactome.file <- "./enrichr/fds/Reactome.xlsx"
fds.reactome <- read.xlsx(fds.reactome.file) 
fds.reactome$Database <- "Reactome"
GetKappaCluster(file_path_sans_ext(fds.reactome.file), fds.reactome, bf.fds.sig$Symbol, "Log.Bayes.Factor")

fds.kegg.file <- "./enrichr/fds/KEGG.xlsx"
fds.kegg <- read.xlsx(fds.kegg.file) 
fds.kegg$Database <- "KEGG"

fds.gobiol.final <- slice(fds.gobiol, c(1,4,6,22))
fds.gomole.final <- slice(fds.gomole, c(7,20))
fds.reactome.final <- slice(fds.reactome, c(1,4))
fds.kegg.final <- slice(fds.kegg, 3)

fds.enrichr.final <- rbind(fds.gobiol.final, fds.gomole.final, fds.reactome.final, fds.kegg.final)
EnrichrPlot(fds.enrichr.final, bf.fds.final.df,  "fds", plot.height = 4, plot.width = 8, plot.title = "Functional Stage")

