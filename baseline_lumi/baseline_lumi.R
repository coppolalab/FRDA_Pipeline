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
library(dplyr)
library(tidyr)

#Functional programming
library(magrittr)
library(purrr)
library(rlist)

#Boxplot function - must fix
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
    pval.cols <- colnames(de.table) %>% str_detect("P.Value") %>% which
    adj.pval.cols <- colnames(de.table) %>% str_detect("adj.P.Val") %>% which
    coef.cols <- colnames(de.table) %>% str_detect("logFC") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = de.table)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(de.table), rule = "<0.005", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = adj.pval.cols, rows = 1:nrow(de.table), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(de.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:12, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 11)
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Plot number of genes at each threshold
DecidePlot <- function(decide.plot, filename, width.plot = 6, height.plot = 7) {
    decide.ggplot <- ggplot() + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = Contrast, y = Num.Genes), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    decide.ggplot <- decide.ggplot + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Contrast, y = Num.Genes, ymax = max(Num.Genes) + 110, label = Num.Genes), hjust = -0.3, position = position_dodge(width = 1))
    decide.ggplot <- decide.ggplot + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = Contrast, y = Num.Genes), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    decide.ggplot <- decide.ggplot + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Contrast, y = Num.Genes, ymax = min(Num.Genes) - 110, label = abs(Num.Genes)), hjust = 1.3, position = position_dodge(width = 1))
    decide.ggplot <- decide.ggplot + facet_grid(Test + Num ~ .) 
    decide.ggplot <- decide.ggplot + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    decide.ggplot <- decide.ggplot + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0))
    CairoPDF(filename, width = width.plot, height = height.plot)
    print(decide.ggplot)
    dev.off()
}

#Volcano plot
VolcanoPlot <- function(top.table, filename, pval.column = "P.Value", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Log.Pvalue")
{
    top.table$Log.Pvalue <- -log10(top.table[[pval.column]])
    p <- ggplot(top.table, aes_string(x = log.column, y = "Log.Pvalue")) + geom_point(aes(color = Significant))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(legend.position = "none")
    p <- p + xlab(xlabel) + ylab(ylabel)
    CairoPNG(str_c(filename, ".png"), width = 300, height = 300)
    print(p)
    dev.off()
}

GeneBoxplot <- function(lumi.object, gene.symbol) {
    gene.expr <- as.vector(exprs(lumi.object[gene.symbol,]))
    gene.df <- data.frame(Status = lumi.object$Status, Expression = gene.expr)
    gene.df$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(gene.df, aes(x = Status, y = Expression, fill = Status)) + geom_violin() + geom_boxplot(width = 0.1) + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    CairoPDF(str_c(gene.symbol, ".pdf"), width = 6, height = 6)
    print(p)
    dev.off()
}

#GO functions
EnrichrSubmit <- function(dataset, enrichr.terms, colname) {
    dir.create(file.path("./enrichr", colname), showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, GetEnrichrData, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]

    names(enrichr.data) <- enrichr.names

    map(names(enrichr.data), EnrichrWorkbook, enrichr.data, colname)
    enrichr.data
}

EnrichrWorkbook <- function(database, full.df, colname) {
    dataset <- full.df[[database]]

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    
    dir.create(file.path("./enrichr", colname), recursive = TRUE)
    filename = str_c(file.path("./enrichr", colname, database), ".xlsx")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

EnrichrPlot <- function(enrichr.df, enrichr.expr, filename, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Log.pvalue <- -(log10(enrichr.df$Adj.P.value))
    enrichr.updown <- map(enrichr.df$Genes, UpDown, enrichr.expr) %>% reduce(rbind)
    colnames(enrichr.updown) <- c("Up", "Down")
    enrichr.df <- cbind(enrichr.df, enrichr.updown)
    enrichr.df$Log.Up <- enrichr.df$Log.pvalue * enrichr.df$Up / enrichr.df$Gene.Count
    enrichr.df$Log.Down <- enrichr.df$Log.pvalue * enrichr.df$Down / enrichr.df$Gene.Count
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- paste(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")", sep = "")
    enrichr.df %<>% arrange(Log.pvalue)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Up, Log.Down) %>% gather(Direction, Length, -Format.Name) 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Length, fill = Direction)) + geom_bar(stat = "identity") + geom_text(label = c(as.character(enrichr.df$Format.Name), rep("", nrow(enrichr.df))), hjust = "left", aes(y = 0.1)) + scale_fill_discrete(name = "Direction", labels = c("Up", "Down"))
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
    CairoPDF(filename, height = plot.height, width = plot.width)
    print(p)
    dev.off()
}

UpDown <- function(filter.vector, enrichr.df) {
    split.vector <- str_split(filter.vector, ",")[[1]]
    enrichr.filter <- filter(enrichr.df, is.element(Symbol, split.vector))
    enrichr.vector <- c("Up" = length(which(sign(enrichr.filter$z.value) == 1)), "Down" = length(which(sign(enrichr.filter$z.value) == -1)))
    enrichr.vector
}

FilterEnrichr <- function(enrichr.table, p.value, filename, cluster = FALSE, ebam.df = data.frame()) {
    enrichr.table$Num.Genes <- map(enrichr.table$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
    enrichr.table %<>% filter(Num.Genes > 4) %>% filter(P.value < p.value)
    if(cluster == TRUE) {
        GetKappaCluster(enrichr.table, ebam.df$Symbol, file_path_sans_ext(filename))
    }
    enrichr.table
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
toptable.cc <- topTable(fitb, coef = 1, n = Inf) 
colnames(toptable.cc) %<>% str_c(".cc")
toptable.cc$Symbol <- rownames(toptable.cc)

toptable.pco <- topTable(fitb, coef = 2, n = Inf) 
colnames(toptable.pco) %<>% str_c(".pco")
toptable.pco$Symbol <- rownames(toptable.pco)

toptable.pca <- topTable(fitb, coef = 3, n = Inf) 
colnames(toptable.pca) %<>% str_c(".pca")
toptable.pca$Symbol <- rownames(toptable.pca)

SaveRDSgz(toptable.pco, "./save/toptable.pco.rda")
SaveRDSgz(toptable.pca, "./save/toptable.pca.rda")
SaveRDSgz(toptable.cc, "./save/toptable.cc.rda")

toptable.all <- left_join(toptable.cc, toptable.pco) %>% left_join(toptable.pca)

#Retrieve annotation information from Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(toptable.all$Symbol), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")

toptable.annot <- left_join(toptable.all, bm.table) %>% 
    select(Symbol, Definition, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val"), dplyr::matches("t\\."), dplyr::matches("B\\."), dplyr::contains("AveExpr")) %>%
    arrange(P.Value.pco)

DEWorkbook(toptable.annot, "toptable_annotated.xlsx")

#Siggenes
siggene.sam.cc <- limma2sam(fitb, coef = 1) 
sam2excel(siggene.sam.cc, 0.7, "siggene.sam.cc.csv")
siggene.sam.pco <- limma2sam(fitb, coef = 2) 
sam2excel(siggene.sam.pco, 0.6, "siggene.sam.pco.csv")
siggene.sam.pca <- limma2sam(fitb, coef = 3) 
sam2excel(siggene.sam.pca, 1.1, "siggene.sam.pca.csv")

siggene.ebam.cc <- limma2ebam(fitb, coef = 1) 
ebam2excel(siggene.ebam.cc, 0.85, "siggene.ebam.cc.csv")
siggene.ebam.pco <- limma2ebam(fitb, coef = 2, delta = 0.85) 
ebam2excel(siggene.ebam.pco, 0.85, "siggene.ebam.pco.csv")
siggene.ebam.pca <- limma2ebam(fitb, coef = 3, delta = 0.85) 
ebam2excel(siggene.ebam.pca, 0.85, "siggene.ebam.pca.csv")

ebam.pco.df <- read_csv("./siggene.ebam.pco.csv", skip = 15)
ebam.pca.df <- read_csv("./siggene.ebam.pca.csv", skip = 15)
colnames(ebam.pco.df)[5] <- "Symbol"
colnames(ebam.pca.df)[5] <- "Symbol"

#Generate statisical cutoff plot
decide <- list(c("fdr", 0.05), c("fdr", 0.10), c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot <- map_df(decide, Decide, fitb, FALSE) %>% gather(Contrast, Num.Genes, -Test, -Num, -Direction) #Compute significance cutoffs

DecidePlot(decide.plot, "threshold_selection")
decide.final <- decideTests(fitb, adjust.method = "none", p.value = 0.001) 
decide.final.df <- data.frame(decide.final)
decide.final.df$Symbol <- rownames(decide.final.df)

#Make volcano plots
toptable.pco$Significant <- factor(toptable.pco$Symbol %in% ebam.pco.df$Symbol)
toptable.pca$Significant <- factor(toptable.pca$Symbol %in% ebam.pca.df$Symbol)

VolcanoPlot(toptable.pco, "volcano_pco", pval.column = "P.Value.pco", log.column = "logFC.pco")
VolcanoPlot(toptable.pca, "volcano_pca", pval.column = "P.Value.pca", log.column = "logFC.pca")
#VolcanoPlot(toptable.cc, "volcano_cc", cutoff = 0.001, cutoff.column = "P.Value.cc", log.column = "logFC.cc")

#Anova heatmaps
#Calculate ratios for use in tables
all.controls <- exprs(rmcov.collapse[,rmcov.collapse$Status == "Control"])
all.controls.means <- rowMeans(all.controls)
all.carriers <- exprs(rmcov.collapse[,rmcov.collapse$Status == "Carrier"])
all.carriers.means <- rowMeans(all.carriers)
all.patients <- exprs(rmcov.collapse[,rmcov.collapse$Status == "Patient"])

coef.cc <- all.carriers - all.controls.means
coef.pco <- all.patients - all.controls.means
coef.pca <- all.patients - all.carriers.means

SaveRDSgz(coef.cc, "coef.cc.rda")
SaveRDSgz(coef.pca, "coef.pca.rda")
SaveRDSgz(coef.pco, "coef.pco.rda")

decide.pco <- filter(decide.final.df, pco != 0)$Symbol
decide.pca <- filter(decide.final.df, pca != 0)$Symbol
decide.cc <- filter(decide.final.df, cc != 0)$Symbol

#Adjust with ebam
pco.plot <- coef.pco[decide.pco,]
pca.plot <- coef.pca[decide.pca,]
cc.plot <- coef.cc[decide.cc,]

CairoPDF("anova_heatmap_patient_vs_control", width = 18, height = 18)
heatmap.object <- heatmap.2(as.matrix(pco.plot), col = rev(redgreen(48)), breaks = c(seq(-2, -1.25, 0.25), seq(-1, 1, 0.05), seq(1.25, 2, 0.25)), trace = "none", labCol = "", labRow = "", keysize = 0.9)
dev.off()

CairoPDF("anova_heatmap_patient_vs_carrier", width = 18, height = 18)
heatmap.object <- heatmap.2(as.matrix(pca.plot), col = rev(redgreen(48)), breaks = c(seq(-2, -1.25, 0.25), seq(-1, 1, 0.05), seq(1.25, 2, 0.25)), trace = "none", labCol = "", labRow = "", keysize = 0.9)
dev.off()

CairoPDF("anova_heatmap_carrier_vs_control", width = 18, height = 18)
heatmap.object <- heatmap.2(as.matrix(cc.plot), col = rev(redgreen(48)), breaks = c(seq(-2, -1.25, 0.25), seq(-1, 1, 0.05), seq(1.25, 2, 0.25)), trace = "none", labCol = "", labRow = "", keysize = 0.9)
dev.off()

#P-value density plot
fit.pvals <- data.frame(fitb$p.value)
colnames(fit.pvals) <- c("Carrier_vs._Control", "Patient_vs._Control", "Patient_vs._Carrier")
fit.pvals.plot <- gather(fit.pvals, Contrasts, P.Value)
fit.pvals.plot$Contrasts %<>% str_replace_all("_", " ")
p <- ggplot(fit.pvals.plot, aes(x = P.Value, fill = Contrasts)) + geom_density() + facet_grid(. ~ Contrasts)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title.y = element_blank(), legend.position = "none") 
CairoPDF("pvalue_density", height = 4, width = 12)
print(p)
dev.off()

#Venn Diagram - adjust with ebam
CairoPDF("venn_diagram", width = 6, height = 6)
vennDiagram(decide.final, include = c("up", "down"), counts.col = c("red", "green"), show.include = "FALSE")
dev.off()

#Gene boxplots
GeneBoxplot(rmcov.collapse, "SERPINE2")

#Patient vs. Control P-value Up
GeneBoxplot(rmcov.collapse, "NLRP12")
GeneBoxplot(rmcov.collapse, "NUP214")

#Patient vs. Control P-value Down
GeneBoxplot(rmcov.collapse, "CD27")
GeneBoxplot(rmcov.collapse, "SKAP1")

#Patient vs. Control Fold Change Up
GeneBoxplot(rmcov.collapse, "MMP9")
GeneBoxplot(rmcov.collapse, "ABCA1")

#Patient vs. Control Fold Change Down
GeneBoxplot(rmcov.collapse, "CD8A")

#Submit genes to Enrichr
source('../../code/GO/enrichr.R')
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016") 

EnrichrSubmit(ebam.pco.df, enrichr.terms, "pco")
EnrichrSubmit(ebam.pca.df, enrichr.terms, "pca")

pca.gobiol.file <- "./enrichr/pca/GO_Biological_Process_2015.xlsx" 
pca.gobiol <- read.xlsx(pca.gobiol.file) 
pca.gobiol.filter <- FilterEnrichr(pca.gobiol, 0.01, pca.gobiol.file, TRUE, ebam.pca.df)
pca.gobiol$Database <- "GO Biological Process"

pca.gomole.file <- "./enrichr/pca/GO_Molecular_Function_2015.xlsx"
pca.gomole <- read.xlsx(pca.gomole.file) 
pca.gomole.filter <- FilterEnrichr(pca.gomole, 0.05, pca.gobiol.file, TRUE, ebam.pca.df)
pca.gomole$Database <- "GO Molecular Function"

pca.reactome.file <- "./enrichr/pca/Reactome_2016.xlsx"
pca.reactome <- read.xlsx(pca.reactome.file) 
pca.reactome.filter <- FilterEnrichr(pca.reactome, 0.05, pca.reactome.file, TRUE, ebam.pca.df)
pca.reactome$Database <- "Reactome"

pca.kegg.file <- "./enrichr/pca/KEGG_2016.xlsx"
pca.kegg <- read.xlsx(pca.kegg.file) 
pca.kegg.filter <- FilterEnrichr(pca.kegg, 0.05, pca.kegg.file, TRUE, ebam.pca.df)
pca.kegg$Database <- "KEGG"

pca.gobiol.final <- slice(pca.gobiol, c(1, 24, 62))
pca.gomole.final <- slice(pca.gomole, c(1, 3))
pca.reactome.final <- slice(pca.reactome, c(4, 6))

pca.enrichr <- rbind(pca.gobiol.final, pca.gomole.final, pca.reactome.final)
pca.enrichr$Adj.P.value <- p.adjust(pca.enrichr$P.value, method = "fdr")
EnrichrPlot(pca.enrichr, ebam.pca.df, "pca.enrichr")

#Patient vs. Control
pco.gobiol.file <- "./enrichr/pco/GO_Biological_Process_2015.xlsx"
pco.gobiol <- read.xlsx(pco.gobiol.file) 
pco.gobiol.filter <- FilterEnrichr(pco.gobiol, 0.05, pco.gobiol.file, TRUE, ebam.pco.df)
pco.gobiol$Database <- "GO Biological Process"

pco.gomole.file <- "./enrichr/pco/GO_Molecular_Function_2015.xlsx"
pco.gomole <- read.xlsx(pco.gomole.file) 
pco.gomole.filter <- FilterEnrichr(pco.gomole, 0.05, pco.gomole.file)
pco.gomole$Database <- "GO Molecular Function"

pco.gobiol.final <- slice(pco.gobiol, c(18, 7, 42))
pco.gomole.final <- slice(pco.gomole, 9)

pco.enrichr <- rbind(pco.gobiol.final, pco.gomole.final)
pco.enrichr$Adj.P.value <- p.adjust(pco.enrichr$P.value, method = "fdr")
EnrichrPlot(pco.enrichr, ebam.pco.df, "pco.enrichr") 

