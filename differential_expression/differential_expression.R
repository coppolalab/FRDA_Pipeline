library(lumi)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(annotate)
library(biomaRt)
library(sva)
library(WGCNA)
library(BayesFactor)
library(matrixStats)
library(limma)

library(openxlsx)
library(Cairo)
library(RColorBrewer)
library(tools)
library(R.utils)

library(pryr)
library(stringr)
library(broom)
library(plyr)
library(magrittr)
library(rlist)
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

CovarBayes <- function(gene.name, expr.matrix, trait.df) {
    gene.vector <- expr.matrix[gene.name,]
    trait.df <- data.frame(trait.df, Gene = gene.vector)
    set.seed(12345)
    trait.anova <- lmBF(Gene ~ Sex + Age, data = trait.df) 
    set.seed(12345)
    trait.posterior <- posterior(trait.anova, iterations = 10000) 
    saveRDS(trait.posterior, str_c("./save/posterior/", gene.name), compress = FALSE)
}

DEBayes <- function(gene.name, expr.matrix, trait.df) {
    gene.vector <- expr.matrix[gene.name,]
    trait.df <- data.frame(trait.df, Gene = gene.vector)
    set.seed(12345)
    trait.anova <- lmBF(Gene ~ Status + RIN, data = trait.df) 
    notrait.anova <- lmBF(Gene ~ RIN, data = trait.df) 
    combined <- c(trait.anova, notrait.anova)
    saveRDS(combined, str_c("./save/model/", gene.name), compress = FALSE)
    set.seed(12345)
    trait.posterior <- posterior(trait.anova, iterations = 10000) 
    saveRDS(trait.posterior, str_c("./save/posterior/", gene.name), compress = FALSE)
}

#Make Excel spreadsheet
DEWorkbook <- function(de.table, filename) { 
    bf.cols <- colnames(de.table) %>% str_detect("Log.Bayes.Factor") %>% which
    logFC.cols <- colnames(de.table) %>% str_detect("logFC") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = de.table)

    sig.bf <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = bf.cols, 
                          rows = 1:nrow(de.table), rule = ">0.5", style = sig.bf)
    conditionalFormatting(wb, 1, cols = logFC.cols, 
                          rows = 1:nrow(de.table), style = c("#63BE7B", "white", "red"), 
                          type = "colourScale")

    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:ncol(de.table), widths = "auto")

    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

BayesPlot <- function(bayes.df, filename, plot.title, diff.column, log.column = "logFC", xlabel = "Log Fold Change") {
    sig.df <- filter_(bayes.df, str_c("Log.Bayes.Factor > 0.5 & ", diff.column, " == TRUE"))
    notsig.df <- filter_(bayes.df, str_c("Log.Bayes.Factor < 0.5 | ", diff.column, " == FALSE"))
    plot.title <- str_c(plot.title, " (", nrow(sig.df), "/", nrow(bayes.df), " Transcripts)")

    p <- ggplot() +
         geom_point(aes_string(x = log.column, y = "Log.Bayes.Factor"), color = "gray", notsig.df) + 
         geom_point(aes_string(x = log.column, y = "Log.Bayes.Factor"), color = "mediumblue", sig.df) + 
         theme_bw() + 
         theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               panel.border = element_rect(size = 1, color = "black"),
               legend.position = "none",
               plot.background = element_blank(),
               plot.title = element_text(hjust = 0.5)) + 
         xlab(xlabel) + ylab(expression(paste(Log[10], " Bayes Factor"))) + 
         ggtitle(plot.title) 

    CairoPDF(filename, width = 6, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

EstimateViolinPlot <- function(estimate.df, gene, ylimits = NA) {
    estimate.extract <- as.matrix(estimate.df) %>% data.frame 
    estimate.groups <- select(estimate.extract, matches("^Status")) 
    estimate.add <- sweep(estimate.groups, 1, estimate.extract$mu, "+")
    estimate.plot <- gather(estimate.add, Group, Estimate) 
    estimate.plot$Group %<>% str_replace("Status\\.", "") %>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(estimate.plot, aes(x = Group, y = Estimate, fill = Group)) + 
        geom_violin(scale = "width", trim = FALSE, draw_quantiles = c(0.05, 0.5, 0.95)) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1), 
              legend.position = "none",
              plot.background = element_blank(), 
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank(), 
              axis.ticks.x = element_blank()) +
        ylab("Expression Estimate") + ggtitle(capitalize(gene))
    if (!is.na(ylimits)) {
        p <- p + ylim(ylimits) 
    }
    CairoPDF(str_c(gene, "estimate", sep = "."), width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

GeneBoxplot <- function(lumi.object, gene.symbol) {
    gene.expr <- as.vector(exprs(lumi.object[gene.symbol,]))
    gene.df <- tibble(Status = lumi.object$Status, Expression = gene.expr)
    gene.df$Status %<>% factor(levels = c("Control", "Carrier", "Patient"))

    p <- ggplot(gene.df, aes(x = Status, y = Expression, fill = Status)) + 
        geom_violin() + 
        geom_boxplot(width = 0.1, outlier.shape = NA) + 
        theme_bw() + 
        theme(legend.position = "none", 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", size = 1),
              plot.background = element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_blank(), 
              axis.ticks.x = element_blank()) + 
        ggtitle(gene.symbol) + ylab("VST Normalized Expression")

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

EnrichrPlot <- function(enrichr.df, enrichr.expr, prefix, plot.title, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% 
        map(unlist) %>% map_int(length)
    enrichr.updown <- map(enrichr.df$Genes, UpDown, enrichr.expr, prefix) %>% reduce(rbind)
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
        geom_bar(stat = "identity", size = 1)  + 
        scale_fill_discrete(name = "Direction", labels = c("Up", "Down")) + 
        coord_flip() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", size = 1),  
              plot.background = element_blank(), 
              plot.title = element_text(hjust = 0.5),
              legend.background = element_blank(), 
              axis.title.y = element_blank(), 
              axis.ticks.y = element_blank(), 
              axis.text.y = element_text(size = 12)) +
        ylab("logBF") + ggtitle(plot.title)

    CairoPDF(str_c(prefix, ".enrichr"), height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

UpDown <- function(filter.vector, enrichr.df, prefix) {
    split.vector <- str_split(filter.vector, ",")[[1]]
    enrichr.filter <- filter(enrichr.df, Symbol %in% split.vector)
    log.column <- str_c("logFC.", prefix)
    enrichr.vector <- c("Up" = length(which(sign(enrichr.filter[[log.column]]) == 1)), 
                        "Down" = length(which(sign(enrichr.filter[[log.column]]) == -1)))
    enrichr.vector
}

#Code for getting size of objects in memory
objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

source("../../code/common_functions.R") #Load shared functions
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
duplicate.slides <- c("6303230097_F", "7649540101_A", "6303248029_G", "9534190037_D", "3998573091_J", "3998573091_L", 
                      "3998582035_G", "3998573091_B", "3998573091_C", "200654510056_C", "3998573091_K", "6116733118_B", 
                      "200654510052_A", "200661380016_J", "3998573126_J", "3998573117_C", "5406958033_C", "3998582041_D") %>% paste(collapse = "|") 
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

#batch.colors <- data.frame(Batch = factor(1:19), Color = c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green", "blue4", "red4", "green4")) #Assign batch colors
#lumi.vst$Batch.Color <- left_join(select(pData(lumi.vst), Batch), batch.colors) %>% select(Color) #add variable for batch color
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

#Remove outlier 
expr.annot.connectivity <- ConnectivityPlot("connectivity", lumi.expr.annot, "")
connectivity.outlier <- names(expr.annot.connectivity[expr.annot.connectivity < -2])
remove.key <- !(sampleNames(lumi.expr.annot) %in% connectivity.outlier)

lumi.rmout <- lumi.expr.annot[,remove.key]
SaveRDSgz(lumi.rmout, "./save/lumi.rmout.rda")
lumi.rmout$Site %<>% factor

#Remove batch effect
expr.combat <- ComBat(dat = exprs(lumi.rmout), batch = factor(lumi.rmout$Site)) #Run ComBat
lumi.combat <- lumi.rmout #Create a new lumi object as a copy of lumi.expr.annot
exprs(lumi.combat) <- expr.combat #Update the expression values of the new lumi object to include the new, corrected intensities

expr.combat <- ComBat(dat = exprs(lumi.combat), batch = factor(lumi.combat$Batch)) #Run ComBat
exprs(lumi.combat) <- expr.combat #Update the expression values of the new lumi object to include the new, corrected intensities
SaveRDSgz(lumi.combat, file = "./save/lumi.combat.rda")

#Regenerate plots
#gen.boxplot("combat_intensity_norm.jpg", lumi.expr.annot, batch.colors, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
Histogram("histogram_combat", lumi.combat)
mds.combat <- exprs(lumi.combat) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_status_combat", mds.combat, pData(lumi.expr.annot), "none", "Status") #label PCs by status
MDSPlot("mds_batch_combat", mds.combat, pData(lumi.expr.annot), as.character(batch.colors$Color), "Batch") #label PCs by batch

combat.collapse.expr <- collapseRows(exprs(lumi.combat), getSYMBOL(featureNames(lumi.combat), 'lumiHumanAll.db'), rownames(lumi.combat)) #collapseRows by symbol
combat.collapse <- ExpressionSet(assayData = combat.collapse.expr$datETcollapsed, phenoData = phenoData(lumi.combat))
SaveRDSgz(combat.collapse, file = "./save/combat.collapse.rda")

#Check for collinearity
age.all <- anovaBF(Age ~ Status, pData(combat.collapse)) %>% extractBF
RIN.all <- anovaBF(RIN ~ Status, pData(combat.collapse)) %>% extractBF

sex.all <- contingencyTableBF(Sex ~ Status, pData(combat.collapse)) %>% extractBF
batch.all <- contingencyTableBF(Batch ~ Status, pData(combat.collapse)) %>% extractBF
site.all <- contingencyTableBF(Site ~ Status, pData(combat.collapse)) %>% extractBF

p <- ggplot(pData(combat.collapse), aes(x = Status, y = Age, fill = Status)) + geom_violin() + geom_boxplot(width = 0.1) 
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
catch <- mclapply(featureNames(combat.collapse), CovarBayes, exprs(combat.collapse), select(pData(combat.collapse), Sex, Age), mc.cores = 8)
posterior.covar <- list.files("./save/posterior", full.names = TRUE) %>% map(readRDS)
SaveRDSgz(posterior.covar, "./save/posterior.covar.rda")
covar.coefs <- map(posterior.covar, colMedians) %>% reduce(rbind)
colnames(covar.coefs) <- colnames(posterior.covar[[1]])
SaveRDSgz(covar.coefs, "./save/covar.coefs.rda")

covar.remove <- model.matrix(~ Sex + Age, pData(combat.collapse))[,-1] %>% t
rmcov.collapse.expr <- exprs(combat.collapse) - (covar.coefs[,c("Sex-MALE","Age-Age")] %*% covar.remove)
rmcov.lumi <- combat.collapse #Make a copy of lumi object
exprs(rmcov.lumi) <- rmcov.collapse.expr #Transfer cleaned expression values into new lumi object
SaveRDSgz(rmcov.lumi, "./save/rmcov.lumi.rda")

#Differential Expression
catch <- mclapply(featureNames(rmcov.lumi), DEBayes, exprs(rmcov.lumi), select(pData(rmcov.lumi), Status, RIN), mc.cores = 8)
model.de <- list.files("./save/model", full.names = TRUE) %>% map(readRDS) 
SaveRDSgz(model.de, "./save/model.de.rda")
de.bf <- map(model.de, extractBF) %>% 
    map(extract2, "bf") %>% 
    map_dbl(reduce, divide_by) %>% 
    map_dbl(log10) %>% 
    signif(3)
bf.df <- tibble(Symbol = featureNames(rmcov.lumi), Log.Bayes.Factor = de.bf) 
posterior.status <- list.files("./save/posterior", full.names = TRUE) %>% map(readRDS)
SaveRDSgz(posterior.status, "./save/posterior.status.rda")

posterior.mu.df <- map(posterior.status, magrittr::extract, TRUE, 1) %>% 
    map(quantile, c(0.025, 0.5, 0.975)) %>% 
    reduce(rbind) %>% as.tibble %>%
    set_colnames(c("mu_CI_2.5", "mu_Median", "mu_CI_97.5")) 
posterior.carrier.df <- map(posterior.status, magrittr::extract, TRUE, 2) %>% 
    map(quantile, c(0.025, 0.5, 0.975)) %>% 
    reduce(rbind) %>% as.tibble %>%
    set_colnames(c("Carrier_CI_2.5", "Carrier_Median", "Carrier_CI_97.5")) 
posterior.control.df <- map(posterior.status, magrittr::extract, TRUE, 3) %>% 
    map(quantile, c(0.025, 0.5, 0.975)) %>% 
    reduce(rbind) %>% as.tibble %>%
    set_colnames(c("Control_CI_2.5", "Control_Median", "Control_CI_97.5"))
posterior.patient.df <- map(posterior.status, magrittr::extract, TRUE, 4) %>% 
    map(quantile, c(0.025, 0.5, 0.975)) %>% 
    reduce(rbind) %>% as.tibble %>%
    set_colnames(c("Patient_CI_2.5", "Patient_Median", "Patient_CI_97.5")) 

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), 
                 filters = 'hgnc_symbol', 
                 values = as.character(featureNames(rmcov.lumi)), 
                 mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")
bm.table %<>% filter(!duplicated(Symbol))

posterior.final.df <- cbind(bf.df, posterior.patient.df, posterior.carrier.df, posterior.control.df, posterior.mu.df) %>% 
    as_tibble %>% 
    mutate(logFC.pco = Patient_Median - Control_Median, 
           logFC.pca = Patient_Median - Carrier_Median, 
           logFC.cc = Carrier_Median - Control_Median, 
           pco.diff = Patient_CI_2.5 > Control_CI_97.5 | Patient_CI_97.5 < Control_CI_2.5,
           pca.diff = Patient_CI_2.5 > Carrier_CI_97.5 | Patient_CI_97.5 < Carrier_CI_2.5,
           cc.diff = Carrier_CI_2.5 > Control_CI_97.5 | Carrier_CI_97.5 < Control_CI_2.5) %>%
    mutate_if(is.numeric, signif, digits = 3) %>%
    left_join(bm.table) %>%
    select(Symbol, Definition, Log.Bayes.Factor, contains("logFC"), contains("diff"), 
           contains("mu"), contains("Patient"), contains("Control"), contains("Carrier")) 
bf.df.sig <- filter(posterior.final.df, Log.Bayes.Factor > 0.5)
SaveRDSgz(posterior.final.df, "./save/posterior.final.df.rda")
DEWorkbook(arrange(posterior.final.df, desc(Log.Bayes.Factor)), "table_annotated.xlsx")

BayesPlot(posterior.final.df, "pco.volcano", "Patient vs. Control", "pco.diff", "logFC.pco")
BayesPlot(posterior.final.df, "pca.volcano", "Patient vs. Carrier", "pca.diff", "logFC.pca")
BayesPlot(posterior.final.df, "cc.volcano", "Carrier vs. Control", "cc.diff", "logFC.cc")

pco.sig <- rep(0, nrow(posterior.final.df))
pco.sig[posterior.final.df$Log.Bayes.Factor > 0.5 & posterior.final.df$Patient_CI_2.5 > posterior.final.df$Control_CI_97.5] <- 1
pco.sig[posterior.final.df$Log.Bayes.Factor > 0.5 & posterior.final.df$Patient_CI_97.5 < posterior.final.df$Control_CI_2.5] <- -1
pca.sig <- rep(0, nrow(posterior.final.df))
pca.sig[posterior.final.df$Log.Bayes.Factor > 0.5 & posterior.final.df$Patient_CI_2.5 > posterior.final.df$Carrier_CI_97.5] <- 1
pca.sig[posterior.final.df$Log.Bayes.Factor > 0.5 & posterior.final.df$Patient_CI_97.5 < posterior.final.df$Carrier_CI_2.5] <- -1
cc.sig <- rep(0, nrow(posterior.final.df))
cc.sig[posterior.final.df$Log.Bayes.Factor > 0.5 & posterior.final.df$Carrier_CI_2.5 > posterior.final.df$Control_CI_97.5] <- 1
cc.sig[posterior.final.df$Log.Bayes.Factor > 0.5 & posterior.final.df$Carrier_CI_97.5 < posterior.final.df$Control_CI_2.5] <- -1

venn.matrix <- cbind("Patient vs. Control" = pco.sig, "Patient vs. Carrier" = pca.sig, "Control vs. Carrier" = cc.sig)
rownames(venn.matrix) <- featureNames(rmcov.lumi)

#Venn Diagram - adjust with ebam
CairoPDF("venn_diagram", width = 10, height = 10)
vennDiagram(venn.matrix, include = c("up", "down"), counts.col = c("red", "green"), show.include = "FALSE")
dev.off()

#Load FDS regression
bf.fds.table <- read.xlsx("../WGCNA_GAA/bf.fds.xlsx")
bf.fds.sig <- filter(bf.fds.table, Log.Bayes.Factor > 0.5)
colnames(bf.fds.sig)[3] <- "Log.Bayes.Factor.fds"
shared.sig <- inner_join(bf.df.sig, bf.fds.sig) %>% 
    arrange(desc(abs(logFC.pco)))
write.xlsx(shared.sig, "shared.sig.xlsx")
write.xlsx(arrange(shared.sig, desc(abs(Median))), "shared.sig.fs.xlsx")

#Load WGCNA
module.membership <- read.xlsx("../WGCNA/module_membership.xlsx")
module.magenta <- filter(module.membership, Module == "magenta") 
pco.magenta <- filter(posterior.final.df, Log.Bayes.Factor > 0.5 & pco.diff == TRUE) %>%
    inner_join(module.magenta) %>% arrange(desc(kscaled))

#Load Classification
lasso.table <- read.xlsx("../classification/lasso.table.xlsx") %>% select(-Definition)
lasso.keep <- filter(lasso.table, Coefficient.Patient > 0)
lasso.combined <- left_join(lasso.keep, posterior.final.df) %>% 
    filter(Log.Bayes.Factor > 0.5) %>% arrange(desc(abs(logFC.pco)))

#Submit genes to Enrichr
source('../../code/GO/enrichr.R')
enrichr.terms <- c("GO Biological Process", "GO Molecular Function", "KEGG", "Reactome") 

pca.sig.df <- filter(posterior.final.df, Log.Bayes.Factor > 0.5 & pca.diff == TRUE)
pca.enrichr <- map(enrichr.terms, GetHyper, pca.sig.df$Symbol, posterior.final.df$Symbol)
names(pca.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, pca.enrichr, "pca")

pca.gobiol.file <- "./enrichr/pca/GO Biological Process.xlsx" 
pca.gobiol <- read.xlsx(pca.gobiol.file) 
pca.gobiol$Database <- "GO BP"
GetKappaCluster(file_path_sans_ext(pca.gobiol.file), pca.enrichr[["GO Biological Process"]], pca.sig.df$Symbol, "Log.Bayes.Factor")

pca.gomole.file <- "./enrichr/pca/GO Molecular Function.xlsx"
pca.gomole <- read.xlsx(pca.gomole.file) 
pca.gomole$Database <- "GO MF"

pca.reactome.file <- "./enrichr/pca/Reactome.xlsx"
pca.reactome <- read.xlsx(pca.reactome.file) 
pca.reactome$Database <- "Reactome"
GetKappaCluster(file_path_sans_ext(pca.reactome.file), pca.enrichr[["Reactome"]], pca.sig.df$Symbol, "Log.Bayes.Factor")

pca.kegg.file <- "./enrichr/pca/KEGG.xlsx"
pca.kegg <- read.xlsx(pca.kegg.file) 
pca.kegg$Database <- "KEGG"

pca.gobiol.final <- slice(pca.gobiol, c(1,4,6,21))
pca.gomole.final <- slice(pca.gomole, 4)
pca.kegg.final <- slice(pca.kegg, c(1,5))
pca.reactome.final <- slice(pca.reactome, c(1,8))

pca.enrichr.df <- rbind(pca.gobiol.final, pca.gomole.final, pca.kegg.final, pca.reactome.final)
EnrichrPlot(pca.enrichr.df, posterior.final.df, "pca", "Patient vs. Carrier", plot.height = 4, plot.width = 7.5)

#Patient vs. Control
pco.sig.df <- filter(posterior.final.df, Log.Bayes.Factor > 0.5 & pco.diff == TRUE)
pco.enrichr <- map(enrichr.terms, GetHyper, pco.sig.df$Symbol, posterior.final.df$Symbol)
names(pco.enrichr) <- enrichr.terms
map(enrichr.terms, EnrichrWorkbook, pco.enrichr, "pco")

pco.gobiol.file <- "./enrichr/pco/GO Biological Process.xlsx"
pco.gobiol <- read.xlsx(pco.gobiol.file) 
pco.gobiol$Database <- "GO BP"
GetKappaCluster(file_path_sans_ext(pco.gobiol.file), pco.enrichr[["GO Biological Process"]], pco.sig.df$Symbol, "Log.Bayes.Factor")

pco.gomole.file <- "./enrichr/pco/GO Molecular Function.xlsx"
pco.gomole <- read.xlsx(pco.gomole.file) 
pco.gomole$Database <- "GO MF"

pco.reactome.file <- "./enrichr/pco/Reactome.xlsx"
pco.reactome <- read.xlsx(pco.reactome.file) 
pco.reactome$Database <- "Reactome"
GetKappaCluster(file_path_sans_ext(pco.reactome.file), pco.enrichr[["Reactome"]], pco.sig.df$Symbol, "Log.Bayes.Factor")

pco.kegg.file <- "./enrichr/pco/KEGG.xlsx"
pco.kegg <- read.xlsx(pco.kegg.file) 
pco.kegg$Database <- "KEGG"

pco.gobiol.final <- slice(pco.gobiol, c(1, 2, 8, 15))
pco.gomole.final <- slice(pco.gomole, 4)
pco.kegg.final <- slice(pco.kegg, 1)
pco.reactome.final <- slice(pco.reactome, 1)

pco.enrichr <- rbind(pco.gobiol.final, pco.gomole.final, pco.kegg.final, pco.reactome.final)
EnrichrPlot(pco.enrichr, posterior.final.df, "pco", "Patient vs. Control", plot.height = 3.25, plot.width = 7)  

