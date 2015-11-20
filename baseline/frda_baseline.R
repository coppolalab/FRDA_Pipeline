#Functional programming
library(magrittr)
library(purrr)
library(functional)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)

#Reading and writing tables
library(readr)
library(openxlsx)

#For DE analysis
library(Biobase)
library(marray)
library(limma)
library(MASS)
library(matrixStats)
library(lumiHumanIDMapping)
library(lumi)

#For batch correction and PEER
library(sva)
library(peer)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(parallel)

#Boxplot function
gen.boxplot <- function(filename, dataset, targetset, colorscheme, maintext, ylabtext)
{
    dataset %<>% t %>% data.frame
    dataset.addvars <- mutate(dataset, Sample.Status = rownames(dataset), Batch = targetset$Batch)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Batch"))
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value, fill = factor(Batch))) + geom_boxplot() + theme_bw()
    p <- p + scale_fill_manual(values = colorscheme)
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Oxygen", width = 20 , height = 8)
}

#Heatmap bar function
gen.heatmapbars <- function(batch.colorscheme, diagnosis.colorscheme, targetset)
{
    batch.heatmap <- data.frame("Batch" = seq(1:length(batch.colorscheme)), "Batch.Color" = batch.colorscheme)
    diagnosis.heatmap <- data.frame("Status" = levels(targetset$Status), "Diagnosis.Color" = diagnosis.colors)
    colorscheme <- data.frame("Batch" = targetset$Batch, "Status" = targetset$Status) %>% join(batch.heatmap) %>% join(diagnosis.heatmap)
    colorscheme <- as.matrix(subset(colorscheme, select = c("Batch.Color", "Diagnosis.Color")))
    return(colorscheme)
}

#Heatmap function
gen.heatmap <- function(filename, dataset, targetset, heatmap.bars, maintitle)
{
    intensities1.cor <- cor(dataset)
    CairoPDF(filename, family = "Oxygen-Regular", width = 10, height = 10)
    heatmap.plus(intensities1.cor, col = heat.colors(40), ColSideColors = heatmap.bars, scale = "none", cexCol = 0.17, cexRow = 0.17, main = maintitle)
    dev.off()
}

#IAC detection of outliers vix 
gen.IACcluster <- function(filename, dataset, maintitle)
{
    IAC = cor(dataset, use = "p")
    cluster1 = hclust(as.dist(1 - IAC))
    CairoPDF(filename, family = "Oxygen-Sans", width = 13, height = 10)
    plot(cluster1, main = paste(maintitle, " (no = ", dim(IAC)[2], ")"))
    dev.off()
    return(IAC)
}

#Create plot of standard deviations of all interarray correlations.  
gen.sdplot <- function(filename, dataset, maintitle)
{
    meanIAC <- apply(dataset, 2, mean)
    sdCorr <- sd(meanIAC)
    numbersd <- (meanIAC - mean(meanIAC)) / sdCorr
    numbersd.plot <- data.frame("Sample.Status" = names(numbersd), "Sample.Num" = seq(1:length(numbersd)), "Z.score" = numbersd)

    p <- ggplot(numbersd.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Status) )
    p <- p + geom_text(size = 4, colour = "red")
    p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
    CairoPDF(filename, family = "Oxygen-Regular", width = 10, height = 10)
    print(p)
    dev.off()
    return(numbersd)
}

#Median absolute deviation standardization function
standardize <- function(dataset)
{
    rowmed <- apply(dataset, 1, median)
    rowmad <- apply(dataset, 1, mad)
    rv <- sweep(dataset, 1, rowmed)
    rv <- sweep(rv, 1, rowmad, "/")
    return(rv)
}

#Create heatmap of top genes
gen.topgenes <- function(filename, dataset, heatmap.bars, maintitle, rowmads, num.genes)
{
    top.genes.names <- rowmads[1:num.genes]
    top.genes.intensities <- dataset[top.genes.names,]
    top.genes.dist <- dist(t(standardize(top.genes.intensities)))
    top.genes.clust <- hclust(top.genes.dist)
    top.genes.matrix <- as.matrix(top.genes.dist)

    CairoPDF(filename, family = "Oxygen-Regular", width = 10, height = 10)
    heatmap.plus(top.genes.matrix, col = rev(heat.colors(75)), distfun = function (x) as.dist(x), main = maintitle, scale = "none", ColSideColors = heatmap.bars, cexRow = 0.08, cexCol = 0.08)
    dev.off()
    return(top.genes.dist)
}

#MDS function - may need work
gen.pca <- function(filename, dataset, targetset, colorscheme, variablename)
{
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$Sample.Name, factor(targetset[,variablename]))
    colnames(target.data) <- c("Sample.Status", variablename)
    colnames(dataset.plot) <- c("Sample.Status", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + scale_fill_manual(values = colorscheme) + xlim(-80,80)
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename)
    CairoPDF(file = filename, height = 7, width = 7)
    print(p)
    dev.off()
}

#Run statistical cutoff tests
gen.decide <- function(test, fit.object, write.results)
{
    results <- decideTests(fit.object, adjust.method = test[1], p = as.numeric(test[2])) #Run test at specified cutoff
    if(write.results == TRUE)
    {
        write.fit(file = paste("./fit_", test[1], ".tsv", sep = ""), fit.object, adjust = test[1], results = results)
    }
    num.genes <- length(which(apply(results, 1, function (x) any(x, na.rm = T))))  #Make this better
    mysum <- summary(results)[-2,] #Eliminate the row for no change in expression
    mysum[2,] <- -(mysum[2,])
    colnames(mysum) <- c("Carrier_vs._Control", "Patient_vs._Control", "Patient_vs._Carrier")
    mysum <- data.frame("Test" = paste(test[1], " p<", test[2], sep = ""), "Num" = paste(num.genes, "Genes", sep = " "), "Direction" = c("positive", "negative"), mysum)
    return(mysum)
}

#Plot statistical cutoff tests
gen.decideplot <- function(filename, decide.plot)
{
    decide.plot$variable <- str_replace_all(decide.plot$variable, "_", " ")
    p <- ggplot()
    p <- p + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = variable, y = value), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = variable, y = value, ymax = max(value) + 110, label = value), hjust = -0.3, position = position_dodge(width = 1))
    p <- p + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = variable, y = value), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = variable, y = value, ymax = min(value) - 110, label = abs(value)), hjust = 1.3, position = position_dodge(width = 1))
    if (length(unique(decide.plot$Test)) > 1)
    {
        p <- p + facet_grid(Num + Test ~ .) 
        #p <- p + ggtitle("Threshold Selection")
    }
    else
    {
        p <- p + ggtitle(paste(decide.plot$Test, "\n", decide.plot$Num))
    }
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0)) + ylab("Differentially Expressed Genes")
    CairoPDF(filename, family = "Oxygen", width = 6, height = 7)
    print(p)
    dev.off()
}

gen.pval.hist <- function(filename, fit.pvals)
{
    colnames(fit.pvals) <- c("Carrier_vs._Control", "Patient_vs._Control", "Patient_vs._Carrier")
    fit.pvals.plot <- melt(fit.pvals)
    fit.pvals.plot$Contrasts <- str_replace_all(fit.pvals.plot$Contrasts, "_", " ")
    p <- ggplot(fit.pvals.plot, aes(x = value)) + geom_histogram(binwidth = 1/80) + facet_grid(. ~ Contrasts)
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("P-value distribution across contrasts") + theme(axis.title.y = element_blank()) + xlab("p-value")
    CairoPDF(filename, height = 7, width = 21)
    print(p)
    dev.off()
}

gen.venndiagram <- function(filename, results)
{
    sumV <- colSums(summary(results)[-2,])
    v <- paste(c("Carrier vs. Control", "Patient vs. Control", "Patient vs. Carrier"), " (", sumV, ")", sep = "")
    CairoPDF(filename, width = 6, height = 6)
    vennDiagram(results, names = v, main = "", include = c("up", "down"), counts.col=c(2,3), cex = 0.8)
    dev.off()
}

#Peer analysis
gen.peer <- function(num.factors, intensities, use.covariates, covariates)
{
    model = PEER() #Instantiate PEER object
    PEER_setNk(model, num.factors) #Specify the number of factors to find
    PEER_setPhenoMean(model, as.matrix(t(intensities))) #Provide mean values for each subject
    PEER_setAdd_mean(model, TRUE) #Enable use of mean in linear model
    if (use.covariates == TRUE)
    {
        PEER_setCovariates(model, as.matrix(covariates)) #Add additional covariates
    }
    PEER_setNmax_iterations(model, 1000) #Specify the maximum number of iterations
    PEER_update(model) #Calculate factors
    residuals.PEER = t(PEER_getResiduals(model))
    rownames(residuals.PEER) = rownames(intensities)
    colnames(residuals.PEER) = colnames(intensities)

    write.csv(data.frame(residuals.PEER), file = paste("residuals_", num.factors, sep = "", ".csv"), row.names = FALSE)
    write.csv(data.frame(PEER_getX(model)), file = paste("factor_", num.factors, sep = "", ".csv"), row.names = FALSE)
    write.csv(data.frame(PEER_getW(model)), file = paste("weight_", num.factors, sep = "", ".csv"), row.names = FALSE)
    write.csv(data.frame(PEER_getAlpha(model)), file = paste("precision_", num.factors, sep = "", ".csv"), row.names = FALSE)

    if (require(Cairo))
    {
        CairoPDF(file = paste("model", num.factors, ".pdf", sep = ""), width = 10, height = 10)
        PEER_plotModel(model)
        dev.off()

        CairoPDF(file = paste("precision_", num.factors, ".pdf", sep = ""), width = 10, height = 10)
        plot(PEER_getAlpha(model), col = "red", lwd = 4, main = paste("precision", num.factors, "factor", sep = " "))
        dev.off()
    }
    else
    {
        pdf(file = paste("model", num.factors, ".pdf", sep = ""), width = 10, height = 10)
        PEER_plotModel(model)
        dev.off()

        pdf(file = paste("precision_", num.factors, ".pdf", sep = ""), width = 10, height = 10)
        plot(PEER_getAlpha(model), col = "red", lwd = 4, main = paste("precision", num.factors, "factor", sep = " "))
        dev.off()
    }    
}

#Calculate ratios.  Really needs work!
gen.ratios <- function(dataset, targetset)
{
    all.samples <- data.frame(dataset)
    all.controls <- all.samples[,targetset$Status == "Control"]
    all.controls.means <- rowMeans(all.controls)
    all.carriers <- all.samples[,targetset$Status == "Carrier"]
    all.carriers.means <- rowMeans(all.carriers)
    all.patients <- all.samples[,targetset$Status == "Patient"]

    all.coefficients.car.cont <- all.carriers - all.controls.means
    all.coefficients.pat.cont <- all.patients - all.controls.means
    all.coefficients.pat.car <- all.patients - all.carriers.means
    colnames(all.coefficients.pat.car) <- paste(colnames(all.coefficients.pat.car), "_vs_carr", sep = "")

    all.coefficients <- data.frame("Probe" = rownames(all.coefficients.car.cont), all.coefficients.car.cont, all.coefficients.pat.cont, all.coefficients.pat.car)
    all.samples <- data.frame("Probe" = rownames(all.samples), all.samples)
    colnames(all.samples)[2:length(all.samples)] <- paste(colnames(all.samples[2:length(all.samples)]), "expr", sep = ".") 
    ratio.exp <- merge(all.coefficients, all.samples)
    return(ratio.exp)
}

#Generate fit object
gen.fit <- function(dataset, model.design, correlation.vector, block.vector)
{
    fit <- lmFit(dataset, design = model.design, correlation = correlation.vector, block = block.vector)
    contrasts.anova <- makeContrasts(time1.carrier_vs_time1.control = Carrier - Control, time1.patient_vs_time1.control = Patient - Control, time1.patient_vs_time1.carrier = Patient - Carrier, levels = model.design)
    fit2.anova <- contrasts.fit(fit, contrasts.anova)
    fitb <- eBayes(fit2.anova)
}

gen.workbook <- function(dataset, filename)
{
    pval.cols <- colnames(dataset) %>% str_detect("p.value.") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")
    dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:3, widths = "auto")
    #setColWidths(wb, 1, cols = 2, widths = 15)
    setColWidths(wb, 1, cols = 4, widths = 45)
    setColWidths(wb, 1, cols = 5:7, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Create genelists
gen.tables <- function(dataset, results, annot, ratio.exp, suffix)
{
    treat.de <- data.frame("Probe_Id" = rownames(results), dataset)
    anovalist <- apply(results, 1, any, na.rm = T) %>% which
    treat.de.anova <- treat.de[anovalist,]
    fitsel.ratio <- merge(treat.de.anova, annot) %>% merge(ratio.exp, by.x = "Probe_Id", by.y = "Probe")
    fitsel.ratio2 <- select(fitsel.ratio, Probe_Id, Accession, Symbol, Definition, contains("Coef."), contains("p.value."), F, F.p.value, contains("Res."), contains("t.time1."), A, Species:GI, Protein_Product:Cytoband, Ontology_Component:Synonyms, matches("CHOP|NG")) %>% arrange(desc(F))
    fitsel.return <- fitsel.ratio2

    fitsel.ratio.all <- merge(treat.de, annot) %>% merge(ratio.exp, by.x = "Probe_Id", by.y = "Probe")
    fitsel.ratio.all <- select(fitsel.ratio.all, Probe_Id, Accession, Symbol, Definition, contains("Coef."), contains("p.value."), F, F.p.value, contains("Res."), contains("t.time1."), A, Species:GI, Protein_Product:Cytoband, Ontology_Component:Synonyms, matches("CHOP|NG")) %>% arrange(desc(F))
    coef.cols <- colnames(fitsel.ratio2) %>% str_detect("Coef.") %>% which
    colnames(fitsel.ratio2)[coef.cols] <- c("Coef.Carrier vs. Control", "Coef.Patient vs. Control", "Coef.Patient vs. Carrier")
    gen.workbook(fitsel.ratio2, paste("./significant_geneList_", suffix, "_time1.xlsx", sep = ""))

    #Carrier - Control
    cc.abs <- as.numeric(fitsel.ratio2$"Coef.Carrier vs. Control") %>% abs
    fitsel.cc <- select(fitsel.ratio2, Accession, Symbol, Definition, matches("Coef.Carrier vs. Control"), matches("Coef.Patient vs. Control"), matches("Coef.Patient vs. Carrier")) %>% mutate(cc.abs) %>% arrange(desc(cc.abs)) %>% select(-cc.abs)
    gen.small.workbook(fitsel.cc, paste("./significant_geneList_", suffix, "_time1_cc.xlsx", sep = ""))

    #Patient - Control
    pco.abs <- as.numeric(fitsel.ratio2$"Coef.Patient vs. Control") %>% abs
    fitsel.pco <- select(fitsel.ratio2, Accession, Symbol, Definition, matches("Coef.Patient vs. Control"), matches("Coef.Carrier vs. Control"), matches("Coef.Patient vs. Carrier"), contains("p.value.")) %>% mutate(pco.abs) %>% arrange(desc(pco.abs)) %>% select(-pco.abs)
    gen.small.workbook(fitsel.pco, paste("./significant_geneList_", suffix, "_time1_pco.xlsx", sep = ""))

    #Patient - Carrier
    pca.abs <- as.numeric(fitsel.ratio2$"Coef.Patient vs. Carrier") %>% abs
    fitsel.pca <- select(fitsel.ratio2, Accession, Symbol, Definition, matches("Coef.Patient vs. Carrier"), matches("Coef.Patient vs. Control"), matches("Coef.Carrier vs. Control")) %>% mutate(pca.abs) %>% arrange(desc(pca.abs)) %>% select(-pca.abs)
    gen.small.workbook(fitsel.pca, paste("./significant_geneList_", suffix, "_time1_pca.xlsx", sep = ""))

    write_csv(fitsel.ratio.all, path = paste("./complete_genelist_time1_", suffix, ".csv", sep = ""))
    return(fitsel.return)
}

gen.small.workbook <- function(dataset, filename)
{
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")
    dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 2, widths = "auto")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 3, widths = 45)
    setColWidths(wb, 1, cols = 4:6, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}
    
#Make anova objects
gen.anova <- function(dataset, suffix)
{
    plot.cc <- filter(dataset, Res.time1.carrier_vs_time1.control != 0) %>% select(matches("Car$"))
    plot.pc <- filter(dataset, Res.time1.patient_vs_time1.control != 0) %>% select(matches("Pat$"))
    plot.pcar <- filter(dataset, Res.time1.patient_vs_time1.carrier != 0) %>% select(matches("carr$"))
    gen.anova.heatmap(paste("./9_anova_heatmap_carrier_vs_control", suffix, sep = "_"), plot.cc, "Carriers vs. Controls")
    gen.anova.heatmap(paste("./9_anova_heatmap_patient_vs_control", suffix, sep = "_"), plot.pc, "Patients vs. Controls")
    gen.anova.heatmap(paste("./9_anova_heatmap_patient_vs_carrier", suffix, sep = "_"), plot.pcar, "Patients vs. Carriers")
}

#Generate anova heatmaps
gen.anova.heatmap <- function(filename, dataset, maintitle)
{ CairoPDF(filename, width = 10, height = 10)
    heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 0.3, labRow = "", keysize = 0.9)
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

#Outlier functions
split.outliers <- function(direction, full.data)
{
    dataset <- full.data[[direction]]$Probe_Id
    data.index <- match(dataset, rownames(intensities1.PEER))
    expr.filter <- data.frame(Probe_Id = as.character(rownames(intensities1.PEER)), intensities1.PEER)
    expr.outlier <- slice(data.frame(expr.filter), data.index)
    expr.names <- join(expr.outlier, annot.reduce)$Symbol
    expr.names[expr.names == ""] <- as.character(expr.outlier[expr.names == "",]$Probe_Id)
    print(expr.names)

    outlier.mean <- rowMeans(expr.outlier[,-1])
    outlier.sd <- rowSds(expr.outlier[,-1]) * 3
    outlier.sdline <- outlier.mean - outlier.sd
    outlier.summary <- data.frame(Probe_Id = expr.names, Mean = outlier.mean, SdLine = outlier.sdline)

    expr.outlier$Probe_Id <- expr.names
    expr.melt <- melt(expr.outlier, id.vars = "Probe_Id")

    expr.melt$Status <- str_split(expr.melt$variable, "_") %>% lapply(tail, n = 1) %>% reduce(c)
    expr.melt$Status %<>% str_replace("Pat", "Patient") %>% str_replace("Car", "Carrier") %>% str_replace("Con", "Control")
    expr.melt$Status %<>% factor(levels = c("Patient", "Carrier", "Control"))
    expr.melt %<>% arrange(Probe_Id, Status)

    expr.melt.list <- split(expr.melt, expr.melt$Probe_Id)
    dir.create("./outliers")
    l_ply(names(expr.melt.list), plot.outlier, expr.melt.list, direction)
}

plot.outlier <- function(index, expr.melt, direction)
{
    expr.plot <- expr.melt[[index]]
    mean.value <- mean(expr.plot$value)
    if (direction == "below")
    {
        sdline.value <- mean.value - (3 * sd(expr.plot$value))
    }
    else
    {
        sdline.value <- mean.value + (3 * sd(expr.plot$value))
    }
    expr.plot$variable %<>% factor(levels = expr.plot$variable)
    p <- ggplot(expr.plot, aes(x = variable, y = value, col = Status)) + geom_point() + geom_hline(yintercept = mean.value, color = "green", size = 1.5) + geom_hline(yintercept = sdline.value, color = "orange", size = 1.5) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("Log 2 normalized expression")
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())#+ scale_x_continuous(as.numeric(unique(dataset$x)))
    CairoPDF(file.path("./outliers", paste(index, direction, sep = "_")), height = 6, width = 10)
    print(p)
    dev.off()
}

gen.fisher.stat <- function(dataset, all.status)
{
    test.status <- all.status[colnames(dataset)]
    pvals <- alply(dataset, 1, rbind, test.status) %>% llply(fisher.test) %>% llply(`[`, "p.value") %>% reduce(c) %>% reduce(c)
    return(pvals)
}

gen.fisher <- function(normalized, all.status)
{
    rownames(normalized) <- normalized$Probe_Id 
    normalized %<>% select(-Probe_Id)
    
    normalized.pvals <- combn(normalized, 2, simplify = FALSE, gen.fisher.stat, all.status) %>% reduce(cbind)
    normalized.names <- combn(normalized, 2, simplify = FALSE, colnames) %>% lapply(paste, collapse = "_vs._") %>% reduce(c)
    colnames(normalized.pvals) <- normalized.names
    return(normalized.pvals)
}

gen.totals <- function(status, normalized, threshold)
{
    normalized$Status <- status
    normalized.melt <- melt(normalized, id.vars = "Status", variable.name = "Probe_Id", value.name = "Z.score")
    normalized.above.melt <- normalized.melt
    normalized.above.melt$Z.score <- normalized.above.melt$Z.score > threshold
    normalized.above <- dcast(normalized.above.melt, Probe_Id ~ Status, sum) 

    normalized.below.melt <- normalized.melt
    normalized.below.melt$Z.score <- normalized.below.melt$Z.score < -(threshold)
    normalized.below <- dcast(normalized.below.melt, Probe_Id ~ Status, sum) 

    return(list(normalized.above, normalized.below))
}

gen.permut <- function(normalized, threshold)
{
    status.randomized <- lapply(1:1000, function(i) {sample(normalized$Status, length(normalized$Status))})
    randomized.values <- mclapply(status.randomized, gen.totals, normalized, threshold, mc.cores = 6) 
}

gen.threshold <- function(normalized, threshold, annot.reduce)
{
    normalized <- data.frame(Status = targets1.rmreps$Status, t(normalized)) 

    normalized.melt <- melt(normalized, id.vars = "Status", variable.name = "Probe_Id", value.name = "Z.score")
    all.status <- summary(targets1.rmreps$Status)
    normalized.above <- filter(normalized.melt, Z.score > threshold) %>% dcast(Probe_Id ~ Status) 
    normalized.below <- filter(normalized.melt, Z.score < -(threshold)) %>% dcast(Probe_Id ~ Status) 
    #normalized.above.fisher <- gen.fisher(normalized.above, all.status)
    #normalized.below.fisher <- gen.fisher(normalized.below, all.status) 
    #normalized.above <- data.frame(normalized.above, normalized.above.fisher)
    #normalized.below <- data.frame(normalized.below, normalized.below.fisher)

    normalized.above %<>% join(annot.reduce) %>% select(Probe_Id, Symbol, Carrier:Patient)
    normalized.below %<>% join(annot.reduce) %>% select(Probe_Id, Symbol, Carrier:Patient)
    minimal.wkbk(normalized.above, paste("./normalized.above.", threshold, ".xlsx", sep = ""))
    minimal.wkbk(normalized.below, paste("./normalized.below.", threshold, ".xlsx", sep = ""))
    return(list(normalized.above, normalized.below))
}


minimal.wkbk <- function(dataset, filename)
{
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = 1:ncol(dataset), widths = "auto")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#GO functions
enrichr.submit <- function(colname, dataset, enrichr.terms, subdir)
{
    comparison.up <- paste(colname, '==', '1')
    dataset.up <- filter_(dataset, comparison.up)
    dir.create(file.path("./enrichr", subdir), showWarnings = FALSE)
    colname.formatted <- str_replace_all(colname, "time1\\.", "") %>% str_replace("Res\\.", "")
    write.xlsx(dataset.up, file.path("./enrichr", subdir, paste(colname.formatted, "up.xlsx", sep = "_")))
    comparison.down <- paste(colname, '==', '-1')
    dataset.down <- filter_(dataset, comparison.down)
    write.xlsx(dataset.down, file.path("./enrichr", subdir, paste(colname.formatted, "down.xlsx", sep = "_")))

    up.data <- lapply(enrichr.terms, get.enrichrdata, dataset.up, FALSE)
    down.data <- lapply(enrichr.terms, get.enrichrdata, dataset.down, FALSE)
    up.names <- enrichr.terms[!is.na(up.data)]
    up.data <- up.data[!is.na(up.data)]
    down.names <- enrichr.terms[!is.na(down.data)]
    down.data <- down.data[!is.na(down.data)]

    names(up.data) <- up.names
    names(down.data) <- down.names

    
    lapply(names(up.data), enrichr.wkbk, up.data, colname, subdir, "up")
    lapply(names(down.data), enrichr.wkbk, down.data, colname, subdir, "down")
}

enrichr.wkbk <- function(database, full.df, colname, subdir, direction)
{
    dataset <- full.df[[database]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    colname.formatted <- str_replace_all(colname, "time1\\.", "") %>% str_replace("Res\\.", "")
    
    dir.create(file.path("./enrichr", subdir, colname.formatted, direction), recursive = TRUE)
    filename = paste(file.path("./enrichr", subdir, colname.formatted, direction, database), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

load("../phenotypedata/targets.final.rda")

#Note: Soon to be deprecated in favor of lumi. 
#Read in old intensities and detection p-values and make columns names consistent targets table
intensities.old <- read_csv("../old raw/chop_rawData_639samples.csv") #Load raw .csv
intensities.old %<>% select(-1) #Remove extra column
colnames(intensities.old)[-1] %<>% str_replace("X", "") #Remove "X" from from front of Slide ids
detectionscore.old <- read_csv("../old raw/chop_detcsoPval_639samples.csv")
detectionscore.old %<>% select(-1)
colnames(detectionscore.old)[-1] %<>% str_replace("X", "")

#Read in new raw file, extract intensity and p-value columns, and rename
data.new <- read_csv("../old raw/2014-278 Sample probe profile.csv") 
#write.table(data.new, file = "../raw_data/2014-278.tsv", sep = "\t", row.names = FALSE) #Convert back to tab separated table for lumi.  

data.new %<>% select(-SYMBOL, -(SEARCH_KEY:ACCESSION))
colnames(data.new)[1] <- "Probe"
intensities.new <- select(data.new, Probe, contains("AVG_Signal"))
colnames(intensities.new)[-1] %<>% str_replace(".AVG_Signal", "")
detectionscore.new <- select(data.new, Probe, contains("Detection Pval"))
colnames(detectionscore.new)[-1] %<>% str_replace(".Detection Pval", "")

#Read even newer raw file
data.newer <- read_csv("../old raw/2015-9185 Sample probe profile.csv")
#write.table(data.newer, file = "../raw_data/2015-9185.tsv", sep = "\t", row.names = FALSE) #Convert back to tab separated table for lumi. 

data.newer %<>% select(-SYMBOL, -(SEARCH_KEY:ACCESSION))
colnames(data.newer)[1] <- "Probe"
intensities.newer <- select(data.newer, Probe, contains("AVG_Signal"))
colnames(intensities.newer)[-1] %<>% str_replace(".AVG_Signal", "")

detectionscore.newer <- select(data.newer, Probe, contains("Detection Pval"))
colnames(detectionscore.newer)[-1] %<>% str_replace(".Detection Pval", "")

#Drop duplicate arrays
bad.slides <- c("3998573091_J", "3998573091_L", "3998582035_G") %>% paste(collapse = "|")
intensities.newer %<>% select(-matches(bad.slides))
detectionscore.newer %<>% select(-matches(bad.slides))
targets.final %<>% filter(!grepl(bad.slides, Slide.ID))

#Convert another profile.  Note that this profile contains samples from other experiments!
raw.2011.018 <- read_csv("../raw_data/2011-018 sample probe profile 2015.csv")
write.table(raw.2011.018, file = "../raw_data/2011-018.tsv", sep = "\t", row.names = FALSE)

#Merge old and new intensities into one data frame.  This currently results in probes being dropped because they are missing in some files and not others.  Needs attention from Joe!
intensities <- merge(intensities.old, intensities.new) %>% merge(intensities.newer)
rownames(intensities) <- intensities$Probe
intensities %<>% select(-Probe)
detectionscore.pvalue <- merge(detectionscore.old, detectionscore.new) %>% merge(detectionscore.newer)
rownames(detectionscore.pvalue) <- detectionscore.pvalue$Probe
detectionscore.pvalue %<>% select(-Probe)
#End of deprecated code - other parts of the code will change as well and will be noted

#New analyses should start here
#Note: Probe profiles in .csv format need to be converted to .tsv files because lumi produces errors when reading .csvs
filenames <- list.files("../raw_data/", pattern = "*.txt|*.tsv", full.names = TRUE)
intensities.list <- lapply(filenames, read_tsv)
missingbeads <- function(intensities)
{
    print(str_subset(colnames(intensities), "Signal") %>% length)
    print(str_subset(colnames(intensities), "NBEADS") %>% length)
}
lapply(intensities.list, missingbeads)
intensities.mat <- reduce(intensities.list, merge)
write.table(intensities.mat, "../raw_data/all_batches.tsv", sep = "\t", row.names = FALSE)

lumi.raw <- lumiR("../raw_data/all_batches.tsv", lib.mapping = "lumiHumanIDMapping", checkDupId = TRUE, convertNuID = TRUE)
lumi.vst <- lumiT(lumi.raw)
lumi.norm <- lumiN(lumi.vst, method = "rsn")
lumi.qual <- lumiQ(lumi.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.cutoff <- detectionCall(lumi.qual)
test.expr <- exprs(lumi.qual)[which(test.cutoff > 0),]
symbols.test <- getSYMBOL(rownames(test.expr), 'lumiHumanAll.db')

#Dynamically assign column names by matching slide Ids to sample names.
colnames.new.index <- match(colnames(intensities), targets.final$Slide.ID)
colnames.new <- targets.final[colnames.new.index,]$Sample.Name
colnames(intensities) <- colnames.new
colnames(detectionscore.pvalue) <- colnames.new

#Remove suspect patients as identified by Jen's spreadsheet.
suspect.ids <- c("6165", "6172", "6174")
suspect.samplenames <- filter(targets.final, PIDN %in% suspect.ids)$Sample.Name
targets.final %<>% filter(!PIDN %in% suspect.ids)
intensities %<>% select(-one_of(suspect.samplenames))
detectionscore.pvalue %<>% select(-one_of(suspect.samplenames))
saveRDS.gz(intensities, file = "./save/intensities.rda")
saveRDS.gz(detectionscore.pvalue, file = "./save/detectionscore.pvalue.rda")
saveRDS.gz(targets.final, file = "./save/targets.final.rda")

#Remove unknowns
targets.final.known <- filter(targets.final, Status != "Unknown") %>% filter(!is.na(Draw.Age)) %>% droplevels
intensities.known <- select(intensities, one_of(targets.final.known$Sample.Name)) 
detectionscore.pvalue.known <- select(detectionscore.pvalue, -contains("Unk"))
saveRDS.gz(targets.final.known, file = "./save/targets.final.known.rda")
saveRDS.gz(intensities.known, file = "./save/intensities.known.rda")
saveRDS.gz(detectionscore.pvalue.known, file = "./save/detectionscore.pvalue.known.rda")

#Load annotated array targets data from Illumina
annot <- read_tsv("../Annotation//HumanHT-12_V4_0_R2_15002873_B.txt") 
saveRDS.gz(annot, file = "./save/annot.rda")

#Extract baseline measurements only.  There are multiple time points, based upon the sample number 
targets1 <- filter(targets.final.known, Sample.Num == "1" | Sample.Num == "1r") %>% droplevels
intensities1 <- select(intensities.known, one_of(targets1$Sample.Name))
detectionscore.pvalue1 <- select(detectionscore.pvalue.known, one_of(targets1$Sample.Name)) 
saveRDS.gz(targets1, file = "./save/targets1.rda")
saveRDS.gz(intensities1, file = "./save/intensities1.rda")
saveRDS.gz(detectionscore.pvalue1, file = "./save/detectionscore.pvalue1.rda")

#To be deprecated with lumi
#Normalization at this step will be unnecessary because it will already be normalized by lumi
#Log2 transform intensities, subtract detection score p-values from 1
intensities1.log2 <- as.matrix(log2(intensities1)) 
#End of deprecated code
detectionscore1 <- 1 - detectionscore.pvalue1
saveRDS.gz(detectionscore1, file = "./save/detectionscore1.rda")

#To be deprecated with lumi
#These figures will need to be plotted differently when lumi is used
#Make boxplots of detectionscore and raw intensities
batch.colors <- c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green", "blue4", "red4")
gen.boxplot("1_base_signal.jpg", intensities1.log2, targets1, batch.colors, "Signal intensity not normalized", "Intensity")
gen.boxplot("1_base_detectionscore.jpg", detectionscore1, targets1, batch.colors, "Detection score p-value", "p-value")

#Make heatmap of all gene intensities
diagnosis.colors <- c("magenta", "green", "darkcyan")
heatmap.bars1 <- gen.heatmapbars(batch.colors, diagnosis.colors, targets1) 
gen.heatmap("2_base_heatmap", intensities1.log2, targets1, heatmap.bars1, "Clustering Based on Inter-Array Pearson Coefficient, Not normalized")

#IAC detection of outliers
IAC <- gen.IACcluster("2_base_IAC.pdf", intensities1.log2, "All samples")
gen.sdplot("2_base_IAC_sd", IAC, "label-1 samples")

#Normalize intensities between arrays with quantile normalization
intensities1.norm <- normalizeBetweenArrays(as.matrix(intensities1.log2), method = "quantile")
saveRDS.gz(intensities1.norm, file = "./save/intensities1.norm.rda")
#End of deprecated code

#Regenerate plots
gen.boxplot("3_base_intensity_norm.jpg", intensities1.norm, targets1, batch.colors, "Quantile normalized signal intensity", "Intensity")
gen.heatmap("3_base_heatmap_norm", intensities1.norm, targets1, heatmap.bars1, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm <- gen.IACcluster("3_base_IAC_norm", intensities1.norm, "All samples")
gen.sdplot("3_base_IAC_sd_norm", IAC.norm, "All samples")

#Top 500 and 1000 genes
intensities1.mads <- apply(intensities1.norm, 1, mad)
intensities1.ordered <- order(intensities1.mads, decreasing = TRUE)

top500.dist <- gen.topgenes("4_base_heatmap_mostVariable_500", intensities1.norm, heatmap.bars1, "Clustering Based on the Top 500 Most Variable Genes", intensities1.ordered, 500)
top1000.dist <- gen.topgenes("4_base_heatmap_mostVariable_1000", intensities1.norm, heatmap.bars1, "Clustering Based on the Top 1000 Most Variable Genes", intensities1.ordered, 1000)

#Principal components analysis
cm1 <- cmdscale(top1000.dist, eig = TRUE)
gen.pca("4_base_MDS_Status", cm1, targets1, diagnosis.colors, "Status")
gen.pca("4_base_MDS_Batch", cm1, targets1, batch.colors, "Batch")

#Remove outliers
remove.samples <- "CHOP_122_1_Pat"
intensities1.rm1 <- select(data.frame(intensities1.norm), -one_of(remove.samples))
targets1.rm1 <- filter(targets1, !grepl(remove.samples, Sample.Name))
heatmap.bars1.rm1 <- slice(data.frame(heatmap.bars1), -match(remove.samples, targets1$Sample.Name)) %>% as.matrix
saveRDS.gz(intensities1.rm1, file = "./save/intensities1.rm1.rda")
saveRDS.gz(targets1.rm1, file = "./save/targets.rm1.rda")
saveRDS.gz(heatmap.bars1.rm1, file = "./save/heatmap.bars1.rm1.rda")
gen.heatmap("5_base_heatmap_rm1_notNorm", intensities1.rm1, targets1.rm1, heatmap.bars1.rm1, "Clustering Based on Inter-Array Pearson Coefficient, not normalized")

#Renormalize and regenerate plots
intensities1.norm.rm1 <- normalizeBetweenArrays(as.matrix(intensities1.rm1), method = "quantile") #This will be changed with lumi!
gen.boxplot("5_base_intensity_norm_rm1.jpg", intensities1.norm.rm1, targets1.rm1, batch.colors, "Quantile normalized signal intensity, exclude 1 sample", "Intensity") 
gen.heatmap("5_base_heatmap_norm_rm1", intensities1.norm.rm1, targets1.rm1, heatmap.bars1.rm1, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.rm1 <- gen.IACcluster("5_base_IAC_norm_rm1", intensities1.norm.rm1, "Outlier removed")
sd.raw.norm.rm1 <- gen.sdplot("5_base_IAC_sd_norm_rm1", IAC.norm.rm1, "Outlier removed")
saveRDS.gz(intensities1.norm.rm1, file = "./save/intensities1.norm.rm1")
saveRDS.gz(IAC.norm.rm1, file = "./save/IAC.norm.rm1.rda")
saveRDS.gz(sd.raw.norm.rm1, file = "./save/sd.raw.norm.rm1.rda")

intensities1.mads.rm1 <- apply(intensities1.norm.rm1, 1, mad)
intensities1.ordered.rm1 <- order(intensities1.mads.rm1, decreasing = TRUE)
top500.dist.norm.rm1 <- gen.topgenes("5_base_heatmap_mostVariable_500_norm_rm1", intensities1.norm.rm1, heatmap.bars1.rm1, "Clustering Based on the Top 500 Most Variable Genes", intensities1.ordered.rm1, 500)
top1000.dist.norm.rm1 <- gen.topgenes("5_base_heatmap_mostVariable_1000_norm_rm1", intensities1.norm.rm1, heatmap.bars1.rm1, "Clustering Based on the Top 1000 Most Variable Genes", intensities1.ordered.rm1, 1000)

cm1.rm1 <- cmdscale(top1000.dist.norm.rm1, eig = TRUE)
gen.pca("5_base_MDS_Status_norm_rm1", cm1.rm1, targets1.rm1, diagnosis.colors, "Status")
gen.pca("5_base_MDS_Batch_norm_rm1", cm1.rm1, targets1.rm1, batch.colors, "Batch")

#Remove replicates
reps <- sd.raw.norm.rm1[str_detect(names(sd.raw.norm.rm1), "1r")] %>% abs
orig.key <- str_replace_all(names(reps), "1r", "1")
orig <- sd.raw.norm.rm1[orig.key] %>% abs

reps.final <- data.frame(orig, reps)
max.key <- apply(reps.final, 1, which.max)
reps.names <- data.frame(orig.key, names(reps))
remove.names <- reps.names[cbind(seq_along(max.key), max.key)]

remove.key <- paste(remove.names, collapse = "|")
targets1.rmreps <- filter(targets1.rm1, !grepl(remove.key, Sample.Name))
intensities1.rmreps <- select(data.frame(intensities1.norm.rm1), -one_of(remove.names))
heatmap.bars1.rmreps <- slice(data.frame(heatmap.bars1.rm1), -match(remove.names, targets1.rm1$Sample.Name)) %>% as.matrix
saveRDS.gz(targets1.rmreps, file = "./save/targets1.rmreps.rda")

#Renormalize and regenerate plots
intensities1.rmreps <- normalizeBetweenArrays(as.matrix(intensities1.rmreps), method = "quantile") #This will be changed with Lumi
gen.boxplot("6_base_intensity_rmreps.jpg", intensities1.rmreps, targets1.rmreps, batch.colors, "Quantile normalized signal intensity, exclude replicates", "Intensity")
gen.heatmap("6_base_heatmap_rmreps", intensities1.rmreps, targets1.rmreps, heatmap.bars1.rmreps, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.rmreps <- gen.IACcluster("6_base_IAC_rmreps", intensities1.rmreps, "Replicates removed")
sd.raw.rmreps <- gen.sdplot("6_base_IAC_sd_rmreps", IAC.rmreps, "Replicates removed")
saveRDS.gz(intensities1.rmreps, file = "./save/intensities1.rmreps.rda")
saveRDS.gz(IAC.rmreps, file = "./save/IAC.rmreps.rda")
saveRDS.gz(sd.raw.rmreps, file = "./save/sd.raw.rda")

intensities1.mads.rmreps <- apply(intensities1.rmreps, 1, mad)
intensities1.ordered.rmreps <- order(intensities1.mads.rmreps, decreasing = TRUE)
top500.dist.rmreps <- gen.topgenes("6_base_heatmap_mostVariable_500_rmreps", intensities1.rmreps, heatmap.bars1.rmreps, "Clustering Based on the Top 500 Most Variable Genes", intensities1.ordered.rmreps, 500)
top1000.dist.rmreps <- gen.topgenes("6_base_heatmap_mostVariable_1000_rmreps", intensities1.rmreps, heatmap.bars1.rmreps, "Clustering Based on the Top 1000 Most Variable Genes", intensities1.ordered.rmreps, 1000)

cm1.rmreps <- cmdscale(top1000.dist.rmreps, eig = TRUE)
gen.pca("6_base_MDS_Status_rmreps", cm1.rmreps, targets1.rmreps, diagnosis.colors, "Status")
gen.pca("6_base_MDS_Batch_rmreps", cm1.rmreps, targets1.rmreps, batch.colors, "Batch")

#Everything from here should remain unchanged
#Use ComBat for batch effect correction
model.status <- model.matrix( ~ 0 + factor(targets1.rmreps$Status) )
colnames(model.status) <- c("Carrier", "Control", "Patient")
model.status.reduce <- model.status[,-2]

model.sex <- model.matrix( ~ 0 + factor(targets1.rmreps$Sex) )
colnames(model.sex) <- c("Female", "Male")
model.sex.reduce <- model.sex[,-1]

model.combat <- cbind(model.status.reduce, Male = model.sex.reduce, Age = log2(as.numeric(targets1.rmreps$Draw.Age)), RIN = log2(targets1.rmreps$RIN))

intensities1.combat <- ComBat(dat = intensities1.rmreps, batch = factor(targets1.rmreps$Batch), mod = model.combat)
saveRDS.gz(intensities1.combat, file = "./save/intensities1.combat.rda")

#Regenerate plots
gen.boxplot("7_base_intensity_combat.jpg", intensities1.combat, targets1.rmreps, batch.colors, "Quantile normalized signal intensity, batch effect corrected", "Intensity")
gen.heatmap("7_base_heatmap_combat", intensities1.combat, targets1.rmreps, heatmap.bars1.rmreps, "")
IAC.combat <- gen.IACcluster("7_base_IAC_combat", intensities1.combat, "Outlier removed")
gen.sdplot("7_base_IAC_sd_combat.pdf", IAC.combat, "Outlier removed")

intensities1.mads.combat <- apply(intensities1.combat, 1, mad)
intensities1.ordered.combat <- order(intensities1.mads.combat, decreasing = TRUE)
top500.dist <- gen.topgenes("7_base_heatmap_mostVariable_500_combat", intensities1.combat, heatmap.bars1.rmreps, "Clustering Based on the Top 500 Most Variable Genes", intensities1.ordered.combat, 500)
top1000.dist <- gen.topgenes("7_base_heatmap_mostVariable_1000_combat", intensities1.combat, heatmap.bars1.rmreps, "Clustering Based on the Top 1000 Most Variable Genes", intensities1.ordered.combat, 1000)

cm1.combat <- cmdscale(top1000.dist, eig = TRUE)
gen.pca("7_base_MDS_Status_combat", cm1.combat, targets1.rmreps, diagnosis.colors, "Status")
gen.pca("7_base_MDS_Batch_combat", cm1.combat, targets1.rmreps, batch.colors, "Batch")

#Run PEER analysis and correlate to known covariates
gen.peer(8, intensities1.combat, TRUE, model.combat)
model.PEER_covariate <- read_csv("./factor_8.csv") %>% select(-(X1:X6))
rownames(model.PEER_covariate) <- colnames(intensities1.combat)
colnames(model.PEER_covariate) <- paste("X", 1:ncol(model.PEER_covariate), sep = "")

source("../common_functions.R")

targets1.gaa <- select(targets1.rmreps, Sample.Name, GAA1) %>% filter(!is.na(GAA1))
cor.gaa <- gen.cor(model.PEER_covariate, targets1.gaa)

targets1.onset <- select(targets1.rmreps, Sample.Name, Onset) %>% filter(!is.na(Onset)) 
cor.onset <- gen.cor(model.PEER_covariate, targets1.onset)

PEER.traits.all <- cbind(cor.gaa, cor.onset) %>% data.frame
PEER.traits.pval <- select(PEER.traits.all, contains("p.value")) %>% as.matrix
PEER.traits.cor <- select(PEER.traits.all, -contains("p.value")) %>% as.matrix

text.matrix.PEER <- paste(signif(PEER.traits.cor, 2), '\n(', signif(PEER.traits.pval, 1), ')', sep = '')
dim(text.matrix.PEER) <- dim(PEER.traits.cor)
gen.text.heatmap(PEER.traits.cor, text.matrix.PEER, colnames(PEER.traits.cor), rownames(PEER.traits.cor), "", "PEER factor-trait relationships")

PEER.trait.out <- data.frame(Factor = rownames(PEER.traits.cor), PEER.traits.cor, PEER.traits.pval)
write_csv(PEER.trait.out, "PEER_trait_cor.csv")

PEER.weights <- read_csv("./weight_8.csv") %>% select(-(X1:X6))
PEER.weights.sums <- colSums(abs(PEER.weights)) %>% data.frame
PEER.weights.sums$Factor <- 1:nrow(PEER.weights.sums)
colnames(PEER.weights.sums)[1] <- "Weight"

p <- ggplot(PEER.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
p <- p + theme_bw() + xlab("Factor") + ylab("Weight")
CairoPDF("./PEER_weights", height = 4, width = 6)
print(p)
dev.off()

#Calculate ratios for use in tables
ratio.exp <- gen.ratios(intensities1.combat, targets1.rmreps)
saveRDS.gz(ratio.exp, file = "./save/ratio.exp.rda")

#Linear model fitting
model.cov <- cbind(model.status, Male = model.sex.reduce, Age = log2(as.numeric(targets1.rmreps$Draw.Age)), RIN = log2(targets1.rmreps$RIN))
model.full <- cbind(model.cov, model.PEER_covariate)

block.correlation <- duplicateCorrelation(intensities1.combat, design = select(model.full, Carrier:Patient), block = targets1.rmreps$Family)
fit.object <- gen.fit(intensities1.combat, model.full, targets1.rmreps$Family)
saveRDS.gz(fit.object, file = "./save/fit.object.rda")

#Generate statisical cutoff
decide <- list(c("fdr", 0.05), c("fdr", 0.1), c("none", 0.001), c("none", 0.005), c("none", 0.01))
decide.plot <- ldply(decide, gen.decide, fit.object, FALSE) %>% melt(id.vars = c("Test", "Num", "Direction"))
gen.decideplot("./7_threshold_selection", decide.plot)

decide.final <- gen.decide(c("none", 0.005), fit.object, TRUE) %>% melt(id.vars = c("Test", "Num", "Direction"))
gen.decideplot("./7_selected_threshold", decide.final)
results <- decideTests(fit.object, adjust.method = "none", p = 0.005)

decide.final.fdr <- gen.decide(c("fdr", 0.1), fit.object, TRUE) %>% melt(id.vars = c("Test", "Num", "Direction"))
gen.decideplot("./7_selected_threshold_fdr", decide.final.fdr)
results.fdr <- decideTests(fit.object, adjust.method = "fdr", p = 0.1)

#Make tables
de.object <- read_tsv("./fit_none.tsv")
de.object.fdr <- read_tsv("./fit_fdr.tsv")
fit.selection <- gen.tables(de.object, results, annot, ratio.exp, "pLess005")
fit.selection.fdr <- gen.tables(de.object.fdr, results.fdr, annot, ratio.exp, "fdrLess01")
saveRDS.gz(fit.selection, file = "./save/fit.object.rda")
saveRDS.gz(fit.selection.fdr, file = "./save/fit.object.fdr.rda")

#Export tables for comparison to Van Hauten
ourgenes <- select(fit.selection.fdr, Accession:Definition, Coef.time1.patient_vs_time1.control, p.value.time1.patient_vs_time1.control, p.value.adj.time1.patient_vs_time1.control, Res.time1.patient_vs_time1.control) #%>% filter(Res.time1.patient_vs_time1.control != 0)
write_csv(ourgenes, "ourgenes.csv")

#P-value histogram
gen.pval.hist("./8_hist_pvalue", fit.object$p.value)

#Venn diagram
gen.venndiagram("./8_venn", results)
gen.venndiagram("./8_venn_fdr", results.fdr)

#Anova heatmaps
gen.anova(fit.selection, "none")
gen.anova(fit.selection.fdr, "fdr")

#Enrichr
source("../GO/enrichr.R")
enrichr.nofdr <- select(fit.selection, Probe_Id, Symbol, Res.time1.carrier_vs_time1.control, Res.time1.patient_vs_time1.control, Res.time1.patient_vs_time1.carrier)
enrichr.fdr <- select(fit.selection.fdr, Probe_Id, Symbol, Res.time1.carrier_vs_time1.control, Res.time1.patient_vs_time1.control, Res.time1.patient_vs_time1.carrier)

comparison.cols <- names(enrichr.nofdr[-(1:2)])
enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 

lapply(comparison.cols, enrichr.submit, enrichr.fdr, enrichr.terms, "fdr")
lapply(comparison.cols, enrichr.submit, enrichr.nofdr, enrichr.terms, "nofdr")

#Move to a separate file
#Outlier analysis
#This needs to be changed to use limma instead of linear algebra parlor tricks!
PEER.factors <- read_csv("./factor_8.csv") %>% select(-X1, -X2, -X6)
PEER.weights.outlier <- read_csv("./weight_8.csv") %>% select(-X1, -X2, -X6)
PEER.residuals <- as.matrix(PEER.factors) %*% t(as.matrix(PEER.weights.outlier)) %>% t
intensities1.PEER <- intensities1.combat - PEER.residuals
colnames(intensities1.PEER) <- colnames(intensities1.combat)
rownames(intensities1.PEER) <- rownames(intensities1.combat)
saveRDS.gz(intensities1.PEER, file = "./save/intensities1.peer.rda")

annot.reduce <- select(annot, Symbol, Probe_Id)
means <- rowMeans(intensities1.PEER)
sds <- rowSds(intensities1.PEER)
normalized <- sweep(intensities1.PEER, 1, means) %>% sweep(1, sds, "/")
normalized.abs <- abs(normalized)
normalized.3 <- gen.threshold(normalized, 3, annot.reduce)
normalized.4 <- gen.threshold(normalized, 4, annot.reduce)
normalized.5 <- gen.threshold(normalized, 5, annot.reduce)
saveRDS.gz(normalized, file = "./save/normalized.rda")

normalized.3.sorted <- lapply(normalized.3, arrange, Control, desc(Patient), desc(Carrier)) %>% lapply(slice, c(1:4))
names(normalized.3.sorted) <- c("above", "below")
lapply(names(normalized.3.sorted), split.outliers, normalized.3.sorted)

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

#PCA Analysis of DE genes
fdr.pca <- filter(fit.selection.fdr, Res.time1.patient_vs_time1.carrier != 0) 
fdr.pca.coef <- select(fdr.pca, contains("_vs_carr")) %>% dist %>% cmdscale(eig = TRUE)
fdr.pca.res <- select(fdr.pca, Res.time1.patient_vs_time1.carrier) %>% as.matrix %>% as.vector
de.pca(fdr.pca.coef, fdr.pca.res, "fdr_patient_vs_carrier")

fdr.pco <- filter(fit.selection.fdr, Res.time1.patient_vs_time1.control != 0) 
fdr.pco.coef <- select(fdr.pco, contains("_Pat")) %>% dist %>% cmdscale(eig = TRUE)
fdr.pco.res <- select(fdr.pco, Res.time1.patient_vs_time1.control) %>% as.matrix %>% as.vector
de.pca(fdr.pco.coef, fdr.pco.res, "fdr_patient_vs_control")

fdr.cc <- filter(fit.selection.fdr, Res.time1.carrier_vs_time1.control != 0) 
fdr.cc.coef <- select(fdr.cc, contains("_Car")) %>% dist %>% cmdscale(eig = TRUE)
fdr.cc.res <- select(fdr.cc, Res.time1.carrier_vs_time1.control) %>% as.matrix %>% as.vector
de.pca(fdr.cc.coef, fdr.cc.res, "fdr_carrier_vs_control")

nofdr.pca <- filter(fit.selection, Res.time1.patient_vs_time1.carrier != 0) 
nofdr.pca.coef <- select(nofdr.pca, contains("_vs_carr")) %>% dist %>% cmdscale(eig = TRUE)
nofdr.pca.res <- select(nofdr.pca, Res.time1.patient_vs_time1.carrier) %>% as.matrix %>% as.vector
de.pca(nofdr.pca.coef, nofdr.pca.res, "nofdr_patient_vs_carrier")

nofdr.pco <- filter(fit.selection, Res.time1.patient_vs_time1.control != 0) 
nofdr.pco.coef <- select(nofdr.pco, contains("_Pat")) %>% dist %>% cmdscale(eig = TRUE)
nofdr.pco.res <- select(nofdr.pco, Res.time1.patient_vs_time1.control) %>% as.matrix %>% as.vector
de.pca(nofdr.pco.coef, nofdr.pco.res, "nofdr_patient_vs_control")

nofdr.cc <- filter(fit.selection, Res.time1.carrier_vs_time1.control != 0) 
nofdr.cc.coef <- select(nofdr.cc, contains("_Car")) %>% dist %>% cmdscale(eig = TRUE)
nofdr.cc.res <- select(nofdr.cc, Res.time1.carrier_vs_time1.control) %>% as.matrix %>% as.vector
de.pca(nofdr.cc.coef, nofdr.cc.res, "nofdr_carrier_vs_control")

de.pca <- function(object.coef, vector.results, filename)
{
    pca.points <- data.frame(object.coef$points)
    colnames(pca.points) <- c("PCA1", "PCA2")
    vector.results[vector.results == -1] <- "Down" 
    vector.results[vector.results == "1"] <- "Up" 
    pca.points %<>% mutate(Result = vector.results)
    p <- ggplot(pca.points, aes(x = PCA1, y = PCA2, col = factor(Result))) + geom_point()
    p <- p + theme_bw() + scale_color_discrete(name = "Direction")
    CairoPDF(filename, width = 8, height = 7)
    print(p)
    dev.off()
}


