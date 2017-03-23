#For DE
library(lumi)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(limma)
library(annotate)
library(biomaRt)
library(siggenes)
library(sva)
library(broom)
library(WGCNA)
library(tools)

#String operations
library(stringr)

#Reading and writing tables
library(readr)
library(openxlsx)

#Plotting
library(ggplot2)
library(Cairo)
library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)

#Data arrangement
library(plyr)
library(dplyr)
library(tidyr)

#Functional programming
library(magrittr)
library(purrr)
library(rlist)

#Boxplot function
BoxPlot <- function(filename, lumi.object, colorscheme, maintext, ylabtext) {
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Batch = lumi.object$Batch)
    dataset.m <- gather(dataset.addvars, id = c("Sample.Status", "Batch"))
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value, fill = factor(Batch))) + geom_boxplot() + theme_bw()
    p <- p + scale_fill_manual(values = colorscheme)
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Noto Sans", width = 20 , height = 8)
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
    connectivity.zscore
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

#Make Excel spreadsheet
DEWorkbook <- function(de.table, filename) { 
    zscore.cols <- colnames(de.table) %>% str_detect("Z.score") %>% which
    posterior.cols <- colnames(de.table) %>% str_detect("Posterior") %>% which
    logFC.cols <- colnames(de.table) %>% str_detect("logFC") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = de.table)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = posterior.cols, rows = 1:nrow(de.table), rule = ">0.9", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = zscore.cols, rows = 1:nrow(de.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    conditionalFormatting(wb, 1, cols = logFC.cols, rows = 1:nrow(de.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:12, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 11)
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

BayesPlot <- function(siggene, filename, threshold, plot.title, posterior.column = "Posterior", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Posterior Probability") {
    siggene$Significant <- siggene[[posterior.column]] > threshold 
    siggene$Significant %<>% factor(levels = c(TRUE, FALSE)) %>% revalue(c("TRUE" = "Yes", "FALSE" = "No"))
    p <- ggplot(siggene, aes_string(x = log.column, y = posterior.column)) + geom_point(aes(color = Significant))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
    p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + xlab(xlabel) + ylab(ylabel) + ggtitle(plot.title) + scale_color_manual(values = c("orange", "darkgreen"), name = "Significant", labels = c("Yes", "No"))
    CairoPDF(filename, width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

GeneBoxplot <- function(lumi.object, gene.symbol) {
    gene.expr <- as.vector(exprs(lumi.object[gene.symbol,]))
    gene.df <- data.frame(Status = lumi.object$Status, Expression = gene.expr)
    gene.df$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(gene.df, aes(x = Status, y = Expression, fill = Status)) + geom_violin() + geom_boxplot(width = 0.1, outlier.shape = NA) + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle(gene.symbol) + theme(plot.background = element_blank(), panel.border = element_rect(color = "black", size = 1))
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + ylab("VST Normalized Expression")
    CairoPDF(str_c(gene.symbol, ".pdf"), width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

#GO functions
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

EnrichrPlot <- function(enrichr.df, enrichr.expr, filename, plot.title, plot.height = 5, plot.width = 8, opposite = FALSE) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Log.Bayes.Factor <- log10(enrichr.df$Bayes.Factor)
    enrichr.updown <- map(enrichr.df$Genes, UpDown, enrichr.expr, opposite) %>% reduce(rbind)
    colnames(enrichr.updown) <- c("Down", "Up")
    enrichr.df <- cbind(enrichr.df, enrichr.updown)
    enrichr.df$Log.Up <- enrichr.df$Log.Bayes.Factor * enrichr.df$Up / enrichr.df$Gene.Count
    enrichr.df$Log.Down <- enrichr.df$Log.Bayes.Factor * enrichr.df$Down / enrichr.df$Gene.Count
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.Bayes.Factor)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Up, Log.Down) %>% gather(Direction, Length, -Format.Name) 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Length, fill = Direction)) 
    p <- p + geom_bar(stat = "identity", size = 1) 
    #p <- p + geom_text(label = c(as.character(enrichr.df$Format.Name), rep("", nrow(enrichr.df))), hjust = "left", aes(y = 0.1)) 
    p <- p + scale_fill_discrete(name = "Direction", labels = c("Up", "Down")) 
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' Bayes Factor')))
    p <- p + theme(plot.background = element_blank(), legend.background = element_blank(), axis.text.y = element_text(size = 12), panel.border = element_rect(color = "black", size = 1))
    p <- p + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(plot.title)
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

UpDown <- function(filter.vector, enrichr.df, opposite = FALSE) {
    split.vector <- str_split(filter.vector, ",")[[1]]
    enrichr.filter <- filter(enrichr.df, is.element(Symbol, split.vector))
    enrichr.vector <- c("Up" = length(which(sign(enrichr.filter$Z.score) == 1)), "Down" = length(which(sign(enrichr.filter$Z.score) == -1)))
    enrichr.vector
}

#Code for getting size of objects in memory
objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

source("../common_functions.R") #Load shared functions
targets.final <- ReadRDSgz("../phenotypedata/targets.final.rda") #Load phenotype data

#Drop probes not found in all batches - this step is NOT normally necessary, was needed due to some inconsistencies in the 
#settings used with GenomeStudio which cannot be resolved easily
filenames <- list.files("../raw_data/", pattern = "*.txt|*.tsv", full.names = TRUE) #Get list of batch files
filenames.csv <- list.files("../raw_data/", pattern = "*.csv", full.names = TRUE)
intensities.list <- map(filenames, read_tsv) #Read files in to a list
intensities.csv <- read_csv(filenames.csv)
intensities.list %<>% list.append(intensities.csv)
intensities.mat <- reduce(intensities.list, merge) #Merge all of the batches into one table
write.table(intensities.mat, "../raw_data/all_batches.tsv", sep = "\t", row.names = FALSE) #Save table to disk

#Read in raw data and remove unrelated samples
lumi.raw <- lumiR("../raw_data/all_batches.tsv", lib.mapping = "lumiHumanIDMapping", checkDupId = TRUE, convertNuID = TRUE, QC = FALSE) #Read in the joined batches
frda.key <- match(targets.final$Slide.ID, colnames(lumi.raw)) #Match sample names in targets table to samples in lumi object
lumi.raw <- lumi.raw[,frda.key] #Subset lumi object to get rid of unrelated samples

rownames(targets.final) <- targets.final$Slide.ID #Set the row names to be the same as the sampleID (may be unnecessary)
pData(lumi.raw) <- targets.final #Add phenotype data
SaveRDSgz(lumi.raw, "./save/lumi.raw.rda")

#remove some duplicate arrays
duplicate.slides <- c("6303230097_F", "7649540101_A", "6303248029_G", "9534190037_D", "3998573091_J", "3998573091_L", "3998582035_G", "3998573091_B", "3998573091_C", "200654510056_C", "3998573091_K", "6116733118_B", "200654510052_A", "200661380016_J", "3998573126_J", "3998573117_C", "5406958033_C", "3998582041_D") %>% paste(collapse = "|") 
duplicate.key <- !grepl(duplicate.slides, lumi.raw$Slide.ID) #Find which samples are not duplicates
knownstatus.key <- !grepl("Unknown|UNKNOWN", lumi.raw$Status) #Find which samples don't have missing FA status
notsuspect.key <- !grepl("6165|6172|6174", lumi.raw$PIDN) #Find which samples don't have PIDNs with suspect FA status
knownage.key <- !is.na(lumi.raw$Age) #Find which samples have a valid date of birth
knownrin.key <- !is.na(lumi.raw$RIN) #Find which samples have a valid RIN
knownsex.key <- !is.na(lumi.raw$Sex) #Find which samples have a valid Sex
combined.key <- duplicate.key & knownstatus.key & notsuspect.key & knownage.key & knownrin.key & knownsex.key #Find which samples pass all of the above filters
lumi.known <- lumi.raw[,combined.key] #Subset the lumi object according the values in combined.key
lumi.known$Status %<>% droplevels #remove UNKNOWN from possible values of Status
lumi.known$Sex %<>% droplevels #remove UNKNOWN from possible values of Sex
lumi.known$Sample.Num %<>% str_replace("r", "") %>% as.numeric

lumi.vst <- lumiT(lumi.known) #Perform variance stabilized transformation
SaveRDSgz(lumi.vst, file = "./save/lumi.vst.rda")

batch.colors <- data.frame(Batch = factor(1:19), Color = c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green", "blue4", "red4", "green4")) #Assign batch colors
lumi.vst$Batch.Color <- left_join(select(pData(lumi.vst), Batch), batch.colors) %>% select(Color) #add variable for batch color
#gen.boxplot("vst_intensity_notnorm.jpg", lumi.vst,  batch.colors, "VST Transformed Intensity not normalized", "Intensity") #Make boxplot of transformed intensities

min.samplenum <- group_by(pData(lumi.known), PIDN) %>% summarise(min(Sample.Num)) %>% data.frame
colnames(min.samplenum)[2] <- "Sample.Num"
pdata.baseline <- left_join(min.samplenum, pData(lumi.known))

lumi.baseline <- lumi.vst[,pdata.baseline$Slide.ID]
SaveRDSgz(lumi.baseline, "./save/lumi.baseline.rda")
lumi.norm <- lumiN(lumi.baseline, method = "rsn") #Normalize with robust spline regression
lumi.cutoff <- detectionCall(lumi.norm) #Get the count of probes which passed the detection threshold per sample
lumi.expr <- lumi.norm[which(lumi.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi <- getSYMBOL(rownames(lumi.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.expr.annot <- lumi.expr[!symbols.lumi,] #Drop any probe which is not annotated
SaveRDSgz(lumi.expr.annot, file = "./save/lumi.expr.annot.rda")

#Regenerate plots
#gen.boxplot("baseline_intensity_norm.jpg", lumi.expr.annot, batch.colors, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
Histogram("histogram_norm", lumi.expr.annot)
mds.norm <- exprs(lumi.expr.annot) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_status_norm", mds.norm, pData(lumi.expr.annot), "none", "Status") #label PCs by status
MDSPlot("mds_batch_norm", mds.norm, pData(lumi.expr.annot), as.character(batch.colors$Color), "Batch") #label PCs by batch

#Remove batch effect
model.combat <- model.matrix(~ Status + Sex + Age + RIN, pData(lumi.expr.annot))

expr.combat <- ComBat(dat = exprs(lumi.expr.annot), batch = factor(lumi.expr.annot$Site), mod = model.combat) #Run ComBat
lumi.combat <- lumi.expr.annot #Create a new lumi object as a copy of lumi.expr.annot
exprs(lumi.combat) <- expr.combat #Update the expression values of the new lumi object to include the new, corrected intensities

expr.combat <- ComBat(dat = exprs(lumi.combat), batch = factor(lumi.combat$Batch), mod = model.combat) #Run ComBat
exprs(lumi.combat) <- expr.combat #Update the expression values of the new lumi object to include the new, corrected intensities

SaveRDSgz(lumi.combat, file = "./save/lumi.combat.rda")

#Regenerate plots
#gen.boxplot("combat_intensity_norm.jpg", lumi.expr.annot, batch.colors, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
Histogram("histogram_combat", lumi.expr.annot)
mds.combat <- exprs(lumi.expr.annot) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_status_combat", mds.combat, pData(lumi.expr.annot), "none", "Status") #label PCs by status
MDSPlot("mds_batch_combat", mds.combat, pData(lumi.expr.annot), as.character(batch.colors$Color), "Batch") #label PCs by batch

#Remove outlier 
combat.connectivity <- ConnectivityPlot("connectivity", lumi.combat, "")
connectivity.outlier <- names(combat.connectivity[combat.connectivity < -2])
remove.key <- !(sampleNames(lumi.combat) %in% connectivity.outlier)

lumi.rmout <- lumi.combat[,remove.key]
SaveRDSgz(lumi.rmout, "./save/lumi.rmout.rda")

#Check for collinearity
age.all <- lm(Age ~ Status, pData(lumi.rmout)) %>% anova %>% tidy
age.aov <- aov(Age ~ Status, pData(lumi.rmout)) %>% TukeyHSD
sex.all <- lm(as.integer(Sex) ~ Status, pData(lumi.rmout)) %>% anova %>% tidy
sex.aov <- aov(as.integer(Sex) ~ Status, pData(lumi.rmout)) %>% TukeyHSD
batch.all <- lm(as.integer(Batch) ~ Status, pData(lumi.rmout)) %>% anova %>% tidy
RIN.all <- lm(RIN ~ Status, pData(lumi.rmout)) %>% anova %>% tidy
site.all <- lm(as.integer(factor(Site)) ~ Status, pData(lumi.rmout)) %>% anova %>% tidy

p <- ggplot(pData(lumi.rmout), aes(x = Status, y = Age, fill = Status)) + geom_violin() + geom_boxplot(width = 0.1) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("Age (p < ", signif(age.all$p.value[1], 3), "*)", sep = ""))
CairoPDF("./age_boxplot.pdf", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = Status, y = RIN)) + geom_boxplot(width = 0.25) + geom_violin()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("RIN (p < ", round(RIN.all$p.value[1], 3), ")", sep = ""))
CairoPDF("./rin_boxplot.pdf", height = 3, width = 3)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = Status, fill = factor(Sex))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Sex (p < ", round(sex.all$p.value[1], 3), "*)", sep = "")) 
CairoPDF("./sex_barplot.pdf", height = 3, width = 4)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = Status, fill = factor(Batch))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.position = "none") 
p <- p + ggtitle(paste("Batch (p < ", round(batch.all$p.value[1], 3), "*)", sep = "")) 
CairoPDF("./batch_barplot.pdf", height = 3, width = 3)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = Status, fill = factor(Site))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.position = "none") 
p <- p + ggtitle(paste("Batch (p < ", round(batch.all$p.value[1], 3), "*)", sep = "")) 
CairoPDF("./site_barplot.pdf", height = 3, width = 3)
plot(p)
dev.off()

#Remove effects of covariates
model.cov <- model.matrix(~ Sex + Age + RIN, data = pData(lumi.rmout)) #Create model matrix of covariates to be removed

rmcov.expr <- removeBatchEffect(exprs(lumi.rmout), covariates = model.cov[,-1]) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
rmcov.lumi <- lumi.rmout #Make a copy of lumi object
exprs(rmcov.lumi) <- rmcov.expr #Transfer cleaned expression values into new lumi object
SaveRDSgz(rmcov.lumi, file = "./save/rmcov.lumi.rda")

#Collapse the data by symbol
rmcov.collapse.expr <- collapseRows(exprs(rmcov.lumi), getSYMBOL(featureNames(rmcov.lumi), 'lumiHumanAll.db'), rownames(rmcov.lumi)) #collapseRows by symbol
rmcov.collapse <- ExpressionSet(assayData = rmcov.collapse.expr$datETcollapsed, phenoData = phenoData(rmcov.lumi))
SaveRDSgz(rmcov.collapse, file = "./save/rmcov.collapse.rda")

model.design <- model.matrix(~ 0 + Status, pData(rmcov.collapse)) #Make covariate matrix for limma
colnames(model.design) %<>% str_replace("Status", "")
SaveRDSgz(model.design, file = "./save/model.design.rda")

#Limma fit
fit <- lmFit(rmcov.collapse, model.design)
fit$df.residual <- (fit$df.residual - 5)

contrasts.anova <- makeContrasts(cc = Carrier - Control, pco = Patient - Control, pca = Patient - Carrier, levels = model.design)
fit.anova <- contrasts.fit(fit, contrasts.anova)
fitb <- eBayes(fit.anova)

#Create DE table
#Make top tables for each coefficient
toptable.cc <- topTable(fitb, coef = 1, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.cc)[2:ncol(toptable.cc)] %<>% str_c(".cc")
toptable.cc$Symbol <- rownames(toptable.cc)

toptable.pco <- topTable(fitb, coef = 2, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.pco)[2:ncol(toptable.pco)] %<>% str_c(".pco")
toptable.pco$Symbol <- rownames(toptable.pco)

toptable.pca <- topTable(fitb, coef = 3, n = Inf) %>% select(AveExpr, logFC, t:B)
colnames(toptable.pca)[2:ncol(toptable.pca)] %<>% str_c(".pca")
toptable.pca$Symbol <- rownames(toptable.pca)

SaveRDSgz(toptable.pco, "./save/toptable.pco.rda")
SaveRDSgz(toptable.pca, "./save/toptable.pca.rda")
SaveRDSgz(toptable.cc, "./save/toptable.cc.rda")

#Breakdown into separate comparisons
lumi.pca <- lumi.rmout[,grepl("Patient|Carrier", lumi.rmout$Status)]
lumi.pco <- lumi.rmout[,grepl("Patient|Control", lumi.rmout$Status)]
lumi.cc <- lumi.rmout[,grepl("Carrier|Control", lumi.rmout$Status)]

lumi.pca$Status %<>% droplevels
lumi.pco$Status %<>% droplevels
lumi.cc$Status %<>% droplevels

model.pca <- model.matrix(~ Sex + Age + RIN, data = pData(lumi.pca)) #Create model matrix of covariates to be removed
rmcov.pca.expr <- removeBatchEffect(exprs(lumi.pca), covariates = model.pca[,-1]) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
rmcov.pca <- lumi.pca #Make a copy of lumi object
exprs(rmcov.pca) <- rmcov.pca.expr #Transfer cleaned expression values into new lumi object

#Collapse the data by symbol
pca.collapse.expr <- collapseRows(exprs(rmcov.pca), getSYMBOL(featureNames(rmcov.pca), 'lumiHumanAll.db'), rownames(rmcov.pca)) #collapseRows by symbol
pca.collapse <- ExpressionSet(assayData = pca.collapse.expr$datETcollapsed, phenoData = phenoData(rmcov.pca))
SaveRDSgz(pca.collapse, file = "./save/pca.collapse.rda")

model.pco <- model.matrix(~ Sex + Age + RIN, data = pData(lumi.pco)) #Create model matrix of covariates to be removed
rmcov.pco.expr <- removeBatchEffect(exprs(lumi.pco), covariates = model.pco[,-1]) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
rmcov.pco <- lumi.pco #Make a copy of lumi object
exprs(rmcov.pco) <- rmcov.pco.expr #Transfer cleaned expression values into new lumi object

#Collapse the data by symbol
pco.collapse.expr <- collapseRows(exprs(rmcov.pco), getSYMBOL(featureNames(rmcov.pco), 'lumiHumanAll.db'), rownames(rmcov.pco)) #collapseRows by symbol
pco.collapse <- ExpressionSet(assayData = pco.collapse.expr$datETcollapsed, phenoData = phenoData(rmcov.pco))
SaveRDSgz(pco.collapse, file = "./save/pco.collapse.rda")

model.cc <- model.matrix(~ Sex + Age + RIN, data = pData(lumi.cc)) #Create model matrix of covariates to be removed
rmcov.cc.expr <- removeBatchEffect(exprs(lumi.cc), covariates = model.cc[,-1]) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
rmcov.cc <- lumi.cc #Make a copy of lumi object
exprs(rmcov.cc) <- rmcov.cc.expr #Transfer cleaned expression values into new lumi object

#Collapse the data by symbol
cc.collapse.expr <- collapseRows(exprs(rmcov.cc), getSYMBOL(featureNames(rmcov.cc), 'lumiHumanAll.db'), rownames(rmcov.cc)) #collapseRows by symbol
cc.collapse <- ExpressionSet(assayData = cc.collapse.expr$datETcollapsed, phenoData = phenoData(rmcov.cc))
SaveRDSgz(cc.collapse, file = "./save/cc.collapse.rda")

#EBAM test
a0.pco <- find.a0(pco.collapse, as.integer(pco.collapse$Status), B = 1000, rand = 12345)
ebam.pco <- ebam(a0.pco)
ebam.pco.df <- data.frame(Z.score = ebam.pco@z, Posterior = ebam.pco@posterior)
ebam.pco.df$Symbol <- rownames(ebam.pco.df)
pco.ebam <- ebam.pco.df
ebam.pco.df %<>% arrange(desc(Posterior))
SaveRDSgz(ebam.pco.df, "./save/ebam.pco.df.rda")
write.xlsx(ebam.pco.df, "ebam.pco.xlsx")

a0.pca <- find.a0(pca.collapse, as.integer(pca.collapse$Status), B = 1000, rand = 12345)
ebam.pca <- ebam(a0.pca)
ebam.pca.df <- data.frame(Z.score = ebam.pca@z, Posterior = ebam.pca@posterior)
ebam.pca.df$Symbol <- rownames(ebam.pca.df)
pca.ebam <- ebam.pca.df
ebam.pca.df %<>% arrange(desc(Posterior))
SaveRDSgz(ebam.pca.df, "./save/ebam.pca.df.rda")
write.xlsx(ebam.pca.df, "ebam.pca.xlsx")

cc.collapse$Status %<>% factor(levels = c("Control", "Carrier"))
a0.cc <- find.a0(cc.collapse, as.integer(cc.collapse$Status), B = 1000, rand = 12345)
ebam.cc <- ebam(a0.cc)
ebam.cc.df <- data.frame(Z.score = ebam.cc@z, Posterior = ebam.cc@posterior)
ebam.cc.df$Symbol <- rownames(ebam.cc.df)
ebam.cc.df %<>% arrange(desc(Posterior))
SaveRDSgz(ebam.cc.df, "./save/ebam.cc.df.rda")

#Retrieve annotation information from Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(toptable.all$Symbol), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")

colnames(ebam.pca.df)[1:2] <- str_c("pca.", colnames(ebam.pca.df)[1:2])
pca.join <- left_join(ebam.pca.df, toptable.pca) %>% select(pca.Z.score:AveExpr)

colnames(ebam.pco.df)[1:2] <- str_c("pco.", colnames(ebam.pco.df)[1:2])
pco.join <- left_join(ebam.pco.df, toptable.pco) %>% select(pco.Z.score:AveExpr)

colnames(ebam.cc.df)[1:2] <- str_c("cc.", colnames(ebam.cc.df)[1:2])
cc.join <- left_join(ebam.cc.df, toptable.cc) %>% select(cc.Z.score:AveExpr)

ebam.join <- left_join(pca.join, pco.join) %>% left_join(cc.join)

table.annot <- left_join(ebam.join, bm.table) %>% 
    select(Symbol, Definition, AveExpr, dplyr::contains("logFC"), dplyr::contains("Posterior"), dplyr::contains("Z.score")) %>%
    arrange(desc(Posterior.pco))
table.annot[,3:ncol(table.annot)] <- signif(table.annot[,3:ncol(table.annot)], 3)

DEWorkbook(table.annot, "table_annotated.xlsx")

BayesPlot(table.annot, "ebam_pco", 0.9, "Patient vs. Control", posterior.column = "Posterior.pco", log.column = "logFC.pco")
BayesPlot(table.annot, "ebam_pca", 0.9, "Patient vs. Carrier", posterior.column = "Posterior.pca", log.column = "logFC.pca")

#Anova heatmaps
#Calculate ratios for use in tables
pco.controls <- exprs(pco.collapse[,pco.collapse$Status == "Control"])
pco.controls.means <- rowMeans(pco.controls)
pco.patients <- exprs(pco.collapse[,pco.collapse$Status == "Patient"])

pca.carriers <- exprs(pca.collapse[,pca.collapse$Status == "Carrier"])
pca.carriers.means <- rowMeans(pca.carriers)
pca.patients <- exprs(pca.collapse[,pca.collapse$Status == "Patient"])

coef.pco <- pco.patients - pco.controls.means
coef.pca <- pca.patients - pca.carriers.means

SaveRDSgz(coef.pca, "coef.pca.rda")
SaveRDSgz(coef.pco, "coef.pco.rda")

#Adjust with ebam
pco.plot <- coef.pco[ebam.pco.df$pco.Posterior > 0.9,]
pca.plot <- coef.pca[ebam.pca.df$pca.Posterior > 0.9,]

CairoPDF("anova_heatmap_patient_vs_control", width = 18, height = 18)
heatmap.object <- heatmap.2(as.matrix(pco.plot), col = rev(redgreen(48)), breaks = c(seq(-2, -1.25, 0.25), seq(-1, 1, 0.05), seq(1.25, 2, 0.25)), trace = "none", labCol = "", labRow = "", keysize = 0.9)
dev.off()

CairoPDF("anova_heatmap_patient_vs_carrier", width = 18, height = 18)
heatmap.object <- heatmap.2(as.matrix(pca.plot), col = rev(redgreen(48)), breaks = c(seq(-2, -1.25, 0.25), seq(-1, 1, 0.05), seq(1.25, 2, 0.25)), trace = "none", labCol = "", labRow = "", keysize = 0.9)
dev.off()

pca.venn.df <- arrange(ebam.pca.df, Symbol)
pca.venn.df$Significant <- pca.venn.df$pca.Posterior > 0.9 
pca.venn.df$Decide <- as.integer(pca.venn.df$Significant) * sign(pca.venn.df$pca.Z.score) 
pco.venn.df <- arrange(ebam.pco.df, Symbol)
pco.venn.df$Significant <- pco.venn.df$pco.Posterior > 0.9 
pco.venn.df$Decide <- as.integer(pco.venn.df$Significant) * sign(pco.venn.df$pco.Z.score) 

venn.matrix <- cbind("Patient vs. Control" = pco.venn.df$Decide, "Patient vs. Carrier" = pca.venn.df$Decide)
rownames(venn.matrix) <- featureNames(pca.collapse)

#Venn Diagram - adjust with ebam
CairoPDF("venn_diagram", width = 10, height = 10)
vennDiagram(venn.matrix, include = c("up", "down"), counts.col = c("red", "green"), show.include = "FALSE")
dev.off()

#Patient vs. Control Posterior
GeneBoxplot(rmcov.collapse, "ABCA1")
GeneBoxplot(rmcov.collapse, "MARC1")
GeneBoxplot(rmcov.collapse, "CD27")
GeneBoxplot(rmcov.collapse, "SERPINE2")

#Patient vs. Carrier Posterior
GeneBoxplot(rmcov.collapse, "NOV")
GeneBoxplot(rmcov.collapse, "CORO1C")
GeneBoxplot(rmcov.collapse, "LBH")
GeneBoxplot(rmcov.collapse, "RALGDS")

#Shared Posterior
GeneBoxplot(rmcov.collapse, "MARC1")
GeneBoxplot(rmcov.collapse, "NUP214")

#Shared Log Fold
GeneBoxplot(rmcov.collapse, "MMP9")
GeneBoxplot(rmcov.collapse, "ANPEP")

#Miscellaneous
GeneBoxplot(rmcov.collapse, "ZNF746")

#Submit genes to Enrichr
source('../../code/GO/enrichr.R')
enrichr.terms <- c("GO Biological Process", "GO Molecular Function", "GO Cellular Component", "KEGG", "Reactome", "GTEx Up", "GTEx Down") 

pca.sig <- filter(ebam., Posterior > 0.9)$Symbol
pca.enrichr <- map(enrichr.terms, GetHyper, pca.sig, pca.ebam$Symbol)
names(pca.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, pca.enrichr, "pca")

pca.gobiol.file <- "./enrichr/pca/GO Biological Process.xlsx" 
pca.gobiol <- read.xlsx(pca.gobiol.file) 
GetKappaCluster(file_path_sans_ext(pca.gobiol.file), pca.enrichr[["GO Biological Process"]], pca.ebam$Symbol)
pca.gobiol$Database <- "GO BP"

pca.gomole.file <- "./enrichr/pca/GO Molecular Function.xlsx"
pca.gomole <- read.xlsx(pca.gomole.file) 
pca.gomole$Database <- "GO MF"

pca.reactome.file <- "./enrichr/pca/Reactome.xlsx"
pca.reactome <- read.xlsx(pca.reactome.file) 
GetKappaCluster(file_path_sans_ext(pca.reactome.file), pca.enrichr[["Reactome"]], pca.ebam$Symbol)
pca.reactome$Database <- "Reactome"

pca.kegg.file <- "./enrichr/pca/KEGG.xlsx"
pca.kegg <- read.xlsx(pca.kegg.file) 
pca.kegg$Database <- "KEGG"

pca.gobiol.final <- slice(pca.gobiol, c(104, 8, 46, 58, 82))
pca.gomole.final <- slice(pca.gomole, 1)
pca.kegg.final <- slice(pca.kegg, 9)
pca.reactome.final <- slice(pca.reactome, 7)

pca.enrichr <- rbind(pca.gobiol.final, pca.gomole.final, pca.kegg.final, pca.reactome.final)
EnrichrPlot(pca.enrichr, pca.ebam, "pca.enrichr", "Patient vs. Carrier", plot.height = 4, plot.width = 8.5)

#Patient vs. Control
pco.sig <- filter(ebam.pco.df, Posterior.pco > 0.9)$Symbol
pco.enrichr <- map(enrichr.terms, GetHyper, pco.sig, pco.ebam$Symbol)
names(pco.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, pco.enrichr, "pco")

pco.gobiol.file <- "./enrichr/pco/GO Biological Process.xlsx"
pco.gobiol <- read.xlsx(pco.gobiol.file) 
GetKappaCluster(file_path_sans_ext(pco.gobiol.file), pco.enrichr[["GO Biological Process"]], pco.ebam$Symbol)
pco.gobiol$Database <- "GO BP"

pco.gomole.file <- "./enrichr/pco/GO Molecular Function.xlsx"
pco.gomole <- read.xlsx(pco.gomole.file) 
pco.gomole$Database <- "GO MF"

pco.reactome.file <- "./enrichr/pco/Reactome.xlsx"
pco.reactome <- read.xlsx(pco.reactome.file) 
pco.reactome$Database <- "Reactome"

pco.gobiol.final <- slice(pco.gobiol, c(1, 2, 3, 48))
pco.gomole.final <- slice(pco.gomole, 2)
pco.reactome.final <- slice(pco.reactome, 2)

pco.enrichr <- rbind(pco.gobiol.final, pco.gomole.final, pco.reactome.final)
EnrichrPlot(pco.enrichr, pco.ebam, "pco.enrichr", "Patient vs. Control", plot.height = 3.25, plot.width = 7.5)  

#Pheno table
pdata.final <- pData(rmcov.collapse)[,-ncol(pData(rmcov.collapse))] 
pdata.age <- group_by(pdata.final, Status) %>% summarise(Age = mean(Age)) 
status.all <- as.matrix(table(pdata.final$Status))
status.male <- filter(pdata.final, Sex == "MALE") %>% group_by(Status) %>% summarise(Male = n())
status.female <- filter(pdata.final, Sex == "FEMALE") %>% group_by(Status) %>% summarise(Female = n())
table.csv <- data.frame(Group = c("Carrier", "Control", "Patient", "Total"), Total = c(status.all, sum(status.all)), Male = c(status.male$Male, sum(status.male$Male)), Female = c(status.female$Female, sum(status.female$Female)), Age = c(round(pdata.age$Age, 2), round(mean(pdata.age$Age), 2)))

write_csv(table.csv, "../../FRDA poster/subjects.csv")

lasso.table <- read.xlsx("../classification/lasso.table.xlsx") %>% select(-Definition)
lasso.pca <- filter(lasso.table, Coefficient.pca > 0)
lasso.pco <- filter(lasso.table, Coefficient.pco > 0)
pca.combined <- left_join(lasso.pca, table.annot)
pca.combined.final <- filter(pca.combined, Posterior.pco > 0.9) 
pco.combined <- left_join(lasso.pco, table.annot)
pco.combined.final <- filter(pco.combined, Posterior.pco > 0.9) 

bf.gaa.table <- read.xlsx("../WGCNA_GAA/bf.gaa.xlsx") 
bf.gaa.sig <- filter(bf.gaa.table, Bayes.Factor > 3)
gaa.combined <- left_join(bf.gaa.sig, table.annot)
gaa.combined.filter <- filter(gaa.combined, Posterior.pco > 0.9)

dlk1.interactors <- read_tsv("../dlk1_interactors.txt") %>% data.frame %>% select(Official.Symbol.Interactor.A, Official.Symbol.Interactor.B)
dlk1.unique <- c(dlk1.interactors$Official.Symbol.Interactor.A, dlk1.interactors$Official.Symbol.Interactor.B) %>% unique

de.dlk1 <- filter(table.annot, Symbol %in% dlk1.unique)  
gaa.dlk1 <- filter(bf.gaa.table, Symbol %in% dlk1.unique)  

