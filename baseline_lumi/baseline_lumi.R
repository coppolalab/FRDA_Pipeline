#Some utils
library(R.utils)

#For DE analysis
library(Biobase)
library(marray)
library(limma)
library(MASS)
library(matrixStats)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(annotate)
library(lumi)
library(WGCNA) #for fastcor

#For batch correction and PEER
library(sva)
library(peer)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(parallel)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(vadr)

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

match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|"))

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

#Boxplot function
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
gen.heatmap <- function(filename, lumi.object, maintitle)
{
    intensities1.cor <- corFast(exprs(lumi.object))
    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(intensities1.cor, col = heat.colors(40), ColSideColors = cbind(lumi.object$Batch.Color, lumi.object$Diagnosis.Color), scale = "none", cexCol = 0.17, cexRow = 0.17, main = maintitle)
    dev.off()
}

#IAC detection of outliers vix 
gen.IACcluster <- function(filename, dataset, maintitle)
{
    IAC = corFast(dataset, use = "p")
    cluster1 = flashClust::hclust(as.dist(1 - IAC))
    CairoPDF(filename, width = 13, height = 10)
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

    p <- ggplot(numbersd.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Status) )
    p <- p + geom_text(size = 4, colour = "red")
    p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
    CairoPDF(filename, width = 10, height = 10)
    print(p)
    dev.off()
    return(numbersd)
}

gen.connectivityplot <- function(filename, dataset, maintitle)
{
    norm.adj <- (0.5 + 0.5 * bicor(exprs(dataset)))
    colnames(norm.adj) <- dataset$Sample.Name
    rownames(norm.adj) <- dataset$Sample.Name
    net.summary <- fundamentalNetworkConcepts(norm.adj)
    net.connectivity <- net.summary$Connectivity
    connectivity.zscore <- (net.connectivity - mean(net.connectivity)) / sqrt(var(net.connectivity))

    connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))
    p <- ggplot(connectivity.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Name) )
    p <- p + geom_text(size = 4, colour = "red")
    p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
    CairoPDF(filename, width = 10, height = 10)
    print(p)
    dev.off()
    return(connectivity.zscore)
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
gen.topgenes <- function(filename, lumi.object, maintitle, rowmads, num.genes)
{
    top.genes.names <- rowmads[1:num.genes]
    dataset <- exprs(lumi.object)
    top.genes.intensities <- dataset[top.genes.names,]
    top.genes.dist <- dist(t(standardize(top.genes.intensities)))
    top.genes.clust <- flashClust::hclust(top.genes.dist)
    top.genes.matrix <- as.matrix(top.genes.dist)

    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(top.genes.matrix, col = rev(heat.colors(75)), distfun = function (x) as.dist(x), main = maintitle, scale = "none", ColSideColors = cbind(lumi.object$Batch.Color, lumi.object$Diagnosis.Color), cexRow = 0.08, cexCol = 0.08)
    dev.off()
    return(top.genes.dist)
}

#MDS function - may need work
gen.pca <- function(filename, dataset, targetset, colorscheme, variablename)
{
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$Sample.Name, factor(targetset[[variablename]]))
    colnames(target.data) <- c("Sample.Status", variablename)
    colnames(dataset.plot) <- c("Sample.Status", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + scale_fill_manual(values = colorscheme) 
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
gen.decideplot <- function(filename, decide.plot, width.plot = 6, height.plot = 7)
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
    #else
    #{
        #p <- p + ggtitle(paste(decide.plot$Test, "\n", decide.plot$Num))
    #}
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0))# + ylab("Differentially Expressed Genes")
    CairoPDF(filename, width = width.plot, height = height.plot)
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

#Calculate ratios.  
gen.ratios <- function(expr.final, pdata)
{
    all.controls <- expr.final[,pdata$Status == "Control"]
    all.controls.means <- rowMeans(all.controls)
    all.carriers <- expr.final[,pdata$Status == "Carrier"]
    all.carriers.means <- rowMeans(all.carriers)
    all.patients <- expr.final[,pdata$Status == "Patient"]

    all.coefficients.car.cont <- all.carriers - all.controls.means
    all.coefficients.pat.cont <- all.patients - all.controls.means
    all.coefficients.pat.car <- all.patients - all.carriers.means
    colnames(all.coefficients.pat.car) <- paste(colnames(all.coefficients.pat.car), "_vs_carr", sep = "")

    all.coefficients <- data.frame("Symbol" = rownames(all.coefficients.car.cont), all.coefficients.car.cont, all.coefficients.pat.cont, all.coefficients.pat.car)
    all.samples <- data.frame("Symbol" = rownames(expr.final), expr.final)
    colnames(all.samples)[2:length(all.samples)] <- paste(colnames(all.samples[2:length(all.samples)]), "expr", sep = ".") 
    ratio.exp <- merge(all.coefficients, all.samples)
    return(ratio.exp)
}

gen.group.ratios <- function(expr.final, pdata) 
{
    controls <- expr.final[,pdata$Status == "Control"]
    controls.means <- rowMeans(controls)
    carriers <- expr.final[,pdata$Status == "Carrier"]
    carriers.means <- rowMeans(carriers)
    patients <- expr.final[,pdata$Status == "Patient"]
    patients.means <- rowMeans(patients)

    cc.logratio <- carriers.means - controls.means
    pco.logratio <- patients.means - controls.means
    pca.logratio <- patients.means - carriers.means
    logratio.df <- data.frame(Carrier.Control.log2FC = cc.logratio, Patient.Control.log2FC = pco.logratio, Patient.Carrier.log2FC = pca.logratio)
    return(logratio.df)
}

#Generate fit object
gen.fit.block <- function(dataset, model.design, correlation.vector, block.vector)
{
    fit <- lmFit(dataset, design = model.design, correlation = correlation.vector, block = block.vector)
    contrasts.anova <- makeContrasts(time1.carrier_vs_time1.control = Carrier - Control, time1.patient_vs_time1.control = Patient - Control, time1.patient_vs_time1.carrier = Patient - Carrier, levels = model.design)
    fit2.anova <- contrasts.fit(fit, contrasts.anova)
    fitb <- eBayes(fit2.anova)
}

gen.fit <- function(dataset, model.design)
{
    fit <- lmFit(dataset, design = model.design)
    contrasts.anova <- makeContrasts(time1.carrier_vs_time1.control = Carrier - Control, time1.patient_vs_time1.control = Patient - Control, time1.patient_vs_time1.carrier = Patient - Carrier, levels = model.design)
    fit2.anova <- contrasts.fit(fit, contrasts.anova)
    fitb <- eBayes(fit2.anova)
}

gen.workbook <- function(dataset, filename)
{
    pval.cols <- colnames(dataset) %>% str_detect("p.value.") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")
    #dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

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
gen.tables <- function(dataset, lumi.object, ratio.exp, suffix)
{
    treat.de <- data.frame("Symbol" = rownames(dataset), dataset)
    #colnames(treat.de)[str_detect(colnames(treat.de), "Genes")] %<>% str_replace("Genes\\.", "") %>% tolower %>% capitalize
    

    fitsel.ratio.all <- merge(treat.de, ratio.exp)
    #fitsel.ratio.all$RefSeq <- lookUp(as.character(fitsel.ratio.all$nuID), "lumiHumanAll.db", "REFSEQ") %>% llply(paste, collapse = ",") %>% reduce(c)
    #fitsel.ratio.all$UniGene <- lookUp(as.character(fitsel.ratio.all$nuID), "lumiHumanAll.db", "UNIGENE") %>% llply(paste, collapse = ",") %>% reduce(c)
    #fitsel.ratio.all$EntrezID <- lookUp(as.character(fitsel.ratio.all$nuID), "lumiHumanAll.db", "ENTREZID") %>% llply(paste, collapse = ",") %>% reduce(c)
    fitsel.return.all <- select(fitsel.ratio.all, Symbol, contains("Coef."), contains("p.value."), F, F.p.value, contains("Res."), contains("t.time1."), A, matches("CHOP|NG")) %>% arrange(desc(F))

    anovalist <- apply(select(fitsel.return.all, contains("Res")), 1, any, na.rm = T) %>% which
    fitsel.return <- fitsel.return.all[anovalist,]
    coef.cols <- colnames(fitsel.return) %>% str_detect("Coef.") %>% which
    colnames(fitsel.return)[coef.cols] <- c("Coef.Carrier vs. Control", "Coef.Patient vs. Control", "Coef.Patient vs. Carrier")
    gen.workbook(fitsel.return, paste("./significant_geneList_", suffix, "_time1.xlsx", sep = ""))

    #Carrier - Control
    cc.abs <- as.numeric(fitsel.return$"Coef.Carrier vs. Control") %>% abs
    fitsel.cc <- select(fitsel.return, Symbol, matches("Coef.Carrier vs. Control"), matches("Coef.Patient vs. Control"), matches("Coef.Patient vs. Carrier")) %>% mutate(cc.abs) %>% arrange(desc(cc.abs)) %>% select(-cc.abs)
    gen.small.workbook(fitsel.cc, paste("./significant_geneList_", suffix, "_time1_cc.xlsx", sep = ""))

    #Patient - Control
    pco.abs <- as.numeric(fitsel.return$"Coef.Patient vs. Control") %>% abs
    fitsel.pco <- select(fitsel.return, Symbol, matches("Coef.Patient vs. Control"), matches("Coef.Carrier vs. Control"), matches("Coef.Patient vs. Carrier"), contains("p.value.")) %>% mutate(pco.abs) %>% arrange(desc(pco.abs)) %>% select(-pco.abs)
    gen.small.workbook(fitsel.pco, paste("./significant_geneList_", suffix, "_time1_pco.xlsx", sep = ""))

    #Patient - Carrier
    pca.abs <- as.numeric(fitsel.return$"Coef.Patient vs. Carrier") %>% abs
    fitsel.pca <- select(fitsel.return, Symbol, matches("Coef.Patient vs. Carrier"), matches("Coef.Patient vs. Control"), matches("Coef.Carrier vs. Control")) %>% mutate(pca.abs) %>% arrange(desc(pca.abs)) %>% select(-pca.abs)
    gen.small.workbook(fitsel.pca, paste("./significant_geneList_", suffix, "_time1_pca.xlsx", sep = ""))

    write.csv(fitsel.return.all, paste("./complete_genelist_time1_", suffix, ".csv", sep = ""), row.names = FALSE)
    return(fitsel.return)
}

gen.small.workbook <- function(dataset, filename)
{
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")
    #dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

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
    if (dim(plot.cc)[1] > 0)
    { gen.anova.heatmap(paste("./anova_heatmap_carrier_vs_control", suffix, sep = "_"), plot.cc, "Carriers vs. Controls") }
    gen.anova.heatmap(paste("./anova_heatmap_patient_vs_control", suffix, sep = "_"), plot.pc, "Patients vs. Controls")
    gen.anova.heatmap(paste("./anova_heatmap_patient_vs_carrier", suffix, sep = "_"), plot.pcar, "Patients vs. Carriers")
    return(list(plot.cc, plot.pc, plot.pcar))
}

#Generate anova heatmaps
gen.anova.heatmap <- function(filename, dataset, maintitle)
{ 
    CairoPDF(filename, width = 28, height = 10)
    heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 0.3, labRow = "", labCol = "", keysize = 0.9)
    dev.off()
}

#GO functions
enrichr.submit <- function(colname, dataset, enrichr.terms, subdir)
{
    filter.cond <- paste(paste(colname, '==', '1'), paste(colname, "==", '-1'), sep = "|")
    dataset.submit <- filter_(dataset, filter.cond)
    colname.formatted <- str_replace_all(colname, "time1\\.", "") %>% str_replace("Res\\.", "")
    write.xlsx(dataset.submit, file.path("./enrichr", subdir, paste(colname.formatted, "xlsx", sep = ".")))

    dir.create(file.path("./enrichr", subdir), showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset.submit, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]

    names(enrichr.data) <- enrichr.names

    map(names(enrichr.data), enrichr.wkbk, enrichr.data, colname, subdir, "enrichr")
    return(enrichr.data)
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

#PCA Analysis
subset.pca <- function(res.col, fit.selection, export.lumi, suffix)
{
    genes.key <- which(fit.selection[[res.col]] != 0)
    if (length(genes.key) > 0)
    {
        genes.nuID <- fit.selection$nuID[genes.key]
        exprs.key <- match(genes.nuID, featureNames(export.lumi))
        lumi.subset <- export.lumi[exprs.key,]
        colname.formatted <- str_replace_all(res.col, "time1\\.", "") %>% str_replace("Res\\.", "")

        #Only include the two groups being compared
        split.name <- str_split(colname.formatted, "_") %>% unlist    
        status.first <- capitalize(split.name[1])    
        status.second <- capitalize(split.name[3])    
        status.key <- export.lumi$Status == status.first | export.lumi$Status == status.second
        lumi.subset.status <- lumi.subset[,status.key]
        exprs.subset <- exprs(lumi.subset.status)
            
        exprs.dist <- standardize(exprs.subset) %>% t %>% dist
        mds.expr <- cmdscale(exprs.dist, eig = TRUE)

        filename <- paste("subgroup_pca", colname.formatted, suffix, sep = "_")
        #CairoPDF(filename, height = 20, width = 20)
        #scatterplot3d(mds.expr$points, color = lumi.subset.status$Diagnosis.Color, pch = 1, main = "3 Dimensional MDS", type = "h", lty.hplot = 2, grid = TRUE, box = TRUE)
        #dev.off()
        gen.pca(filename, mds.expr, phenoData(lumi.subset.status), diagnosis.colors, "Status")
        return(mds.expr)
    }
    else
    {
        return(NA)     
    }
}

gen.enrichrplot <- function(enrichr.df, enrichr.expr, filename, plot.height = 5, plot.width = 8)
{
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
    enrichr.df$Log.pvalue <- -(log10(enrichr.df$P.value))
    enrichr.updown <- map(enrichr.df$Genes, get.updown, enrichr.expr) %>% reduce(rbind)
    colnames(enrichr.updown) <- c("Up", "Down")
    enrichr.df <- cbind(enrichr.df, enrichr.updown)
    enrichr.df$Log.Up <- enrichr.df$Log.pvalue * enrichr.df$Up / enrichr.df$Gene.Count
    enrichr.df$Log.Down <- enrichr.df$Log.pvalue * enrichr.df$Down / enrichr.df$Gene.Count
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- paste(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")", sep = "")
    enrichr.df %<>% arrange(Log.pvalue)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Up, Log.Down) %>% melt(id.vars = "Format.Name") 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = c(as.character(enrichr.df$Format.Name), rep("", nrow(enrichr.df))), hjust = "left", aes(y = 0.1))
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "FALSE",  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
    CairoPDF(filename, height = plot.height, width = plot.width)
    print(p)
    dev.off()
}

get.updown <- function(filter.vector, enrichr.df)
{
    split.vector <- str_split(filter.vector, ",")[[1]]
    grep.vector <- match.exact(split.vector)
    enrichr.filter <- filter(enrichr.df, grepl(grep.vector, Symbol))
    enrichr.vector <- factor(enrichr.filter[,2]) %>% summary
    return(enrichr.vector)
}

load("../phenotypedata/targets.final.rda")
colnames(targets.final)[4] <- "sampleID"

#remove some duplicate arrays
bad.slides <- c("3998573091_J", "3998573091_L", "3998582035_G") %>% paste(collapse = "|") 
targets.final %<>% filter(!grepl(bad.slides, sampleID))

#Note: Probe profiles in .csv format need to be converted to .tsv files because lumi produces errors when reading .csvs
#Drop probes not found in all batches - this step is NOT normally necessary, was needed due to some inconsistencies in the 
#settings used with GenomeStudio which cannot be resolved easily
filenames <- list.files("../raw_data/", pattern = "*.txt|*.tsv", full.names = TRUE) #Get list of batch files
intensities.list <- map(filenames, read_tsv) #Read files in to a list
intensities.mat <- reduce(intensities.list, merge) #Merge all of the batches into one table
write.table(intensities.mat, "../raw_data/all_batches.tsv", sep = "\t", row.names = FALSE) #Save table to disk

#Read in raw data and remove unrelated samples
lumi.raw <- lumiR("../raw_data/all_batches.tsv", lib.mapping = "lumiHumanIDMapping", checkDupId = TRUE, convertNuID = TRUE, QC = FALSE) #Read in the joined batches
frda.key <- match(targets.final$sampleID, colnames(lumi.raw)) #Match sample names in targets table to samples in lumi object
lumi.raw <- lumi.raw[,frda.key] #Subset lumi object to get rid of unrelated samples

rownames(targets.final) <- targets.final$sampleID #Set the row names to be the same as the sampleID (may be unnecessary)
pData(lumi.raw) <- targets.final #Add phenotype data
sampleNames(lumi.raw) <- lumi.raw$Sample.Name #Add sample names

knownstatus.key <- !grepl("Unknown|UNKNOWN", lumi.raw$Status) #Find which samples don't have missing FA status
notsuspect.key <- !grepl("6165|6172|6174", lumi.raw$PIDN) #Find which samples don't have PIDNs with suspect FA status
knownage.key <- !is.na(lumi.raw$Draw.Age) #Find which samples have a valid date of birth
combined.key <- knownstatus.key & notsuspect.key & knownage.key #Find which samples pass all three of the above filters
lumi.known <- lumi.raw[,combined.key] #Subset the lumi object according the values in combined.key
lumi.known$Status %<>% droplevels #remove UNKNOWN from possible values of Status

lumi.vst <- lumiT(lumi.known) #Perform variance stabilized transformation
saveRDS.gz(lumi.known, file = "./save/lumi.known.rda")

baseline.key <- grepl("^1$|1r", lumi.vst$Sample.Num) #Find which samples are baseline
lumi.baseline <- lumi.vst[,baseline.key] #Subset lumi object for baseline samples

diagnosis.colors <- c("magenta", "green", "darkcyan")
heatmap.bars <- gen.heatmapbars(batch.colors, diagnosis.colors, phenoData(lumi.baseline)) %>% data.frame #create colors for heatmap bar
lumi.baseline$Batch.Color <- as.character(heatmap.bars$Batch.Color) #add variable for batch color
lumi.baseline$Diagnosis.Color <- as.character(heatmap.bars$Diagnosis.Color) #add variable for diagnosis color
saveRDS.gz(lumi.baseline, file = "./save/lumi.baseline.rda")

batch.colors <- c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green", "blue4", "red4") #Assign batch colors
gen.boxplot("baseline_intensity_notnorm.jpg", lumi.baseline,  batch.colors, "VST Transformed Intensity not normalized", "Intensity") #Make boxplot of transformed intensities

lumi.norm <- lumiN(lumi.baseline, method = "rsn") #Normalize with robust spline regression
lumi.qual <- lumiQ(lumi.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.cutoff <- detectionCall(lumi.qual) #Get the count of probes which passed the detection threshold per sample
lumi.expr <- lumi.qual[which(lumi.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi <- getSYMBOL(rownames(lumi.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.expr.annot <- lumi.expr[!symbols.lumi,] #Drop any probe which is not annotated
saveRDS.gz(lumi.expr.annot, file = "./save/lumi.expr.annot.rda")

#The quality control summary can be used to get an idea of what outliers are in the data
qcsum <- lumi.expr.annot@QC$sampleSummary %>% t %>% data.frame #Extract the quality control summary object
colnames(qcsum) %<>% str_replace("\\.0\\.01\\.", "")
qcsum$Sample.Name <- rownames(qcsum) #Add sample name
qcsum$RIN <- lumi.expr.annot$RIN #Add RIN
qcsum$Sample.Num <- lumi.expr.annot$Sample.Num #Add sample number

#Regenerate plots
gen.boxplot("baseline_intensity_norm.jpg", lumi.expr.annot, batch.colors, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
gen.heatmap("baseline_heatmap_norm", lumi.expr.annot, "Clustering Based on Inter-Array Pearson Coefficient, RSN normalized") #Make heatmap of sample clustering

#These QC plots may be somewhat redundant if you have already looked at the QC object
IAC.norm <- gen.IACcluster("baseline_IAC_norm", exprs(lumi.expr.annot), "All samples") #Make clustering dendrogram of samples
gen.sdplot("baseline_IAC_sd", IAC.norm, "All samples") #Make plot of standard deviations
saveRDS.gz(IAC.norm, file = "./save/IAC.norm.rda")
saveRDS.gz(sd.raw.norm, file = "./save/sd.raw.norm.rda")

#Top 500 and 1000 genes
lumi.mads <- apply(exprs(lumi.expr.annot), 1, mad) #Get median absolute deviations of expression values
lumi.ordered <- order(lumi.mads, decreasing = TRUE) #Use this to rank the genes by the variability of their expression

top500.dist <- gen.topgenes("baseline_heatmap_500", lumi.expr.annot, "Clustering Based on the Top 500 Most Variable Genes", lumi.ordered, 500) #heatmap of top 500 most variable genes
top1000.dist <- gen.topgenes("baseline_heatmap_1000", lumi.expr.annot, "Clustering Based on the Top 1000 Most Variable Genes", lumi.ordered, 1000) #heatmap of top 1000 most variable genes

#Principal components analysis
cm1 <- cmdscale(top1000.dist, eig = TRUE) #get first two principle components
gen.pca("baseline_mds_status", cm1, phenoData(lumi.expr.annot), diagnosis.colors, "Status") #label PCs by status
gen.pca("baseline_mds_batch", cm1, phenoData(lumi.expr.annot), batch.colors, "Batch") #label PCs by batch

gen.connectivityplot("baseline_connectivity", lumi.expr.annot, "")

#Remove outlier
#Potential suspects: CHOP_188_1_Car, CHOP_175_1_Pat
remove.sample <- !grepl("CHOP_122_1_Pat|CHOP_14_1_Pat|CHOP_540_1_Car|CHOP_124_1_Pat|CHOP_96_1_Car", sampleNames(lumi.baseline))
lumi.rmout <- lumi.baseline[,remove.sample]# %>% lumiN(method = "rsn")

#####NOTE: In talking with the Geschwind lab about this, they seem to think that renormalizing after remove such a small number of samples is unnecessary.  I have not rigorously examined this, but I think it is worth pointing out that this may be redundant
lumi.rmout.norm <- lumiN(lumi.rmout, method = "rsn") #Normalize with robust spline regression
lumi.rmout.qual <- lumiQ(lumi.rmout.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.rmout.cutoff <- detectionCall(lumi.rmout.qual) #Get the count of probes which passed the detection threshold per sample
lumi.rmout.expr <- lumi.rmout.qual[which(lumi.rmout.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.rmout <- getSYMBOL(rownames(lumi.rmout.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.rmout.annot <- lumi.rmout.expr[!symbols.lumi.rmout,] #Drop any probe which is not annotated
saveRDS.gz(lumi.rmout.annot, file = "lumi.rmout.annot.rda")

gen.heatmap("rm1_heatmap", lumi.rmout.annot, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.rm1 <- gen.IACcluster("rm1_IAC", exprs(lumi.rmout.annot), "Outlier removed")
sd.raw.norm.rm1 <- gen.sdplot("rm1_sd", IAC.norm.rm1, "Outlier removed")
saveRDS.gz(IAC.norm.rm1, file = "./save/IAC.norm.rm1.rda")
saveRDS.gz(sd.raw.norm.rm1, file = "./save/sd.raw.norm.rm1.rda")

lumi.rmout.mads <- apply(exprs(lumi.rmout.annot), 1, mad)
lumi.rmout.ordered <- order(lumi.rmout.mads, decreasing = TRUE)
top500.dist.norm.rm1 <- gen.topgenes("rm1_heatmap_500", lumi.rmout.annot, "Clustering Based on the Top 500 Most Variable Genes", lumi.rmout.ordered, 500)
top1000.dist.norm.rm1 <- gen.topgenes("rm1_heatmap_1000", lumi.rmout.annot, "Clustering Based on the Top 1000 Most Variable Genes", lumi.rmout.ordered, 1000)

cm1.rm1 <- cmdscale(top1000.dist.norm.rm1, eig = TRUE)
gen.pca("rm1_mds_diagnosis", cm1.rm1, phenoData(lumi.rmout.annot), diagnosis.colors, "Status")
gen.pca("rm1_mds_batch", cm1.rm1, phenoData(lumi.rmout.annot), batch.colors, "Batch")

qcsum.rmout <- lumi.rmout.annot@QC$sampleSummary %>% t %>% data.frame
colnames(qcsum.rmout) %<>% str_replace("\\.0\\.01\\.", "")
qcsum.rmout$Sample.Name <- rownames(qcsum.rmout)
qcsum.rmout$RIN <- lumi.rmout.annot$RIN
qcsum.rmout$Sample.Num <- lumi.rmout.annot$Sample.Num
qcsum.rmout$PIDN <- lumi.rmout.annot$PIDN
arrange(qcsum.rmout, distance.to.sample.mean)
PIDNs <- filter(qcsum.rmout, Sample.Num == "1r")$PIDN %>% paste(collapse = "|")

rmout.connectivity <- gen.connectivityplot("rmout_connectivity", lumi.rmout.annot, "")
connectivity.outlier <- rmout.connectivity[abs(rmout.connectivity) > 2] 
connectivity.sample <- names(connectivity.outlier) %>% as.character %>% match.exact %>% str_detect(sampleNames(lumi.baseline))

lumi.rmout.all <- lumi.rmout.annot[,!connectivity.sample]
rmout.connectivity <- gen.connectivityplot("rmout_all_connectivity", lumi.rmout.all, "")

#Remove replicates
#####NOTE: this is only necessary in the very unusual situation in which you have technical replicates (i.e. the same RNA was run on two microarrays)
##### Other wise you can just skip to the next section
reps <- sd.raw.norm.rm1[str_detect(names(sd.raw.norm.rm1), "1r")] %>% abs
orig.key <- str_replace_all(names(reps), "1r", "1")
orig <- sd.raw.norm.rm1[orig.key] %>% abs

reps.final <- data.frame(orig, reps)
max.key <- apply(reps.final, 1, which.max)
reps.names <- data.frame(orig.key, names(reps))
remove.names <- reps.names[cbind(seq_along(max.key), max.key)]

remove.key <- paste(remove.names, collapse = "|")
reps.samples <- !grepl(remove.key, sampleNames(lumi.baseline))
remove.all <-  remove.sample & reps.samples
saveRDS.gz(remove.all, "./save/remove.all.rda")

lumi.rmreps <- lumi.baseline[,remove.all]
saveRDS.gz(lumi.rmreps, file = "./save/lumi.rmreps.rda")
lumi.rmreps.norm <- lumiN(lumi.rmreps, method = "rsn") #Normalize with robust spline regression
#lumi.rmreps.qual <- lumiQ(lumi.rmreps.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.rmreps.cutoff <- detectionCall(lumi.rmreps.norm) #Get the count of probes which passed the detection threshold per sample
lumi.rmreps.expr <- lumi.rmreps.norm[which(lumi.rmreps.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.rmreps <- getSYMBOL(rownames(lumi.rmreps.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.rmreps.annot <- lumi.rmreps.expr[!symbols.lumi.rmreps,] #Drop any probe which is not annotated
saveRDS.gz(lumi.rmreps.annot, file = "./save/lumi.rmreps.annot.rda")

gen.heatmap("rmreps_heatmap", lumi.rmreps.annot, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.rmreps <- gen.IACcluster("rmreps_IAC", exprs(lumi.rmreps.annot), "Outlier removed")
sd.raw.norm.rmreps <- gen.sdplot("rmreps_sd", IAC.norm.rmreps, "Outlier removed")
saveRDS.gz(IAC.norm.rmreps, file = "./save/IAC.norm.rmreps.rda")
saveRDS.gz(sd.raw.norm.rmreps, file = "./save/sd.raw.norm.rmreps.rda")

lumi.rmreps.mads <- apply(exprs(lumi.rmreps.annot), 1, mad)
lumi.rmreps.ordered <- order(lumi.rmreps.mads, decreasing = TRUE)
top500.dist.norm.rmreps <- gen.topgenes("rmreps_heatmap_500", lumi.rmreps.annot, "Clustering Based on the Top 500 Most Variable Genes", lumi.rmreps.ordered, 500)
top1000.dist.norm.rmreps <- gen.topgenes("rmreps_heatmap_1000", lumi.rmreps.annot, "Clustering Based on the Top 1000 Most Variable Genes", lumi.rmreps.ordered, 1000)

cm1.rmreps <- cmdscale(top1000.dist.norm.rmreps, eig = TRUE)
gen.pca("rmreps_mds_diagnosis", cm1.rmreps, phenoData(lumi.rmreps.annot), diagnosis.colors, "Status")
gen.pca("rmreps_mds_batch", cm1.rmreps, phenoData(lumi.rmreps.annot), batch.colors, "Batch")

#Use ComBat for batch effect correction
#### NOTE: before running this adjustment, please check that the following things are true:
# 1) Other covariates are not correlated
# 2) Batch number is not confounded with another categorical variable (i.e. diagnosis, treatment, sex, etc)
# 3) Each batch has more than one array in it
model.status <- model.matrix( ~ 0 + factor(lumi.rmreps.annot$Status) ) %>% data.frame #Make dummy variables out of Status
colnames(model.status) <- c("Carrier", "Control", "Patient")
model.status.reduce <- model.status[,-2] #Remove control, because one needs n-1 dummy variables for n possible values of categorical variables

model.sex <- model.matrix( ~ 0 + factor(lumi.rmreps.annot$Sex) ) #Make dummy variables out of Sex
colnames(model.sex) <- c("Female", "Male")
model.sex.reduce <- model.sex[,-1] #Remove one column (like done with Status)

model.combat <- data.frame(model.status.reduce, Male = model.sex.reduce, Age = as.numeric(lumi.rmreps.annot$Draw.Age), RIN = lumi.rmreps.annot$RIN)  #Assemble covariate matrix with both categorical and continuous covariates

expr.combat <- ComBat(dat = exprs(lumi.rmreps.annot), batch = factor(lumi.rmreps.annot$Batch), mod = model.combat) #Run ComBat
lumi.combat <- lumi.rmreps.annot #Create a new lumi object as a copy of lumi.rmreps.annot
exprs(lumi.combat) <- expr.combat #Update the expression values of the new lumi object to include the new, corrected intensities
saveRDS.gz(lumi.combat, file = "./save/lumi.combat.rda")

#Regenerate plots
gen.heatmap("combat_heatmap", lumi.combat, "Clustering Based on Inter-Array Pearson Coefficient, Quantile normalized")
IAC.norm.combat <- gen.IACcluster("combat_IAC", exprs(lumi.combat), "Outlier removed")
sd.raw.norm.combat <- gen.sdplot("combat_sd", IAC.norm.combat, "Outlier removed")
saveRDS.gz(IAC.norm.combat, file = "./save/IAC.norm.combat.rda")
saveRDS.gz(sd.raw.norm.combat, file = "./save/sd.raw.norm.combat.rda")

lumi.combat.mads <- apply(exprs(lumi.combat), 1, mad)
lumi.combat.ordered <- order(lumi.combat.mads, decreasing = TRUE)
top500.dist.norm.combat <- gen.topgenes("combat_heatmap_500", lumi.combat, "Clustering Based on the Top 500 Most Variable Genes", lumi.combat.ordered, 500)
top1000.dist.norm.combat <- gen.topgenes("combat_heatmap_1000", lumi.combat, "Clustering Based on the Top 1000 Most Variable Genes", lumi.combat.ordered, 1000)

cm1.combat <- cmdscale(top1000.dist.norm.combat, eig = TRUE)
gen.pca("combat_mds_diagnosis", cm1.combat, phenoData(lumi.combat), diagnosis.colors, "Status")
gen.pca("combat_mds_batch", cm1.combat, phenoData(lumi.combat), batch.colors, "Batch")

#Run PEER analysis and correlate to known covariates
#### NOTE: This may be replaced by HCP in the near future
gen.peer(8, exprs(lumi.combat), TRUE, model.combat)
model.PEER_covariate <- read_csv("./factor_8.csv") %>% select(-(X1:X6))
rownames(model.PEER_covariate) <- colnames(lumi.combat)
colnames(model.PEER_covariate) <- paste("X", 1:ncol(model.PEER_covariate), sep = "")

#Plot PEER weights
PEER.weights <- read_csv("./weight_8.csv") %>% select(-(X1:X6)) #Get weights of PEER factors.  The first 6 columns are dropped because they contain the weights for the other covariates
PEER.weights.sums <- colSums(abs(PEER.weights)) %>% data.frame #Get sums of weights for each factor
PEER.weights.sums$Factor <- 1:nrow(PEER.weights.sums) #Add column to label by factor
colnames(PEER.weights.sums)[1] <- "Weight"

#ggplot of factor versus weight
p <- ggplot(PEER.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
p <- p + theme_bw() + xlab("Factor") + ylab("Weight")
CairoPDF("./PEER_weights", height = 4, width = 6)
print(p)
dev.off()

source("../../code/hidden_covariate_linear.R")

model.hcp <- standardize.hcp(model.combat)
expr.hcp <- exprs(lumi.combat) %>% t %>% standardize.hcp

#hcp.covariates <- hidden_covariate_linear(model.hcp, expr.hcp, k = 10, lambda = 1)
#hcp.factors <- hcp.covariates$Z
#colnames(hcp.factors) <- paste("X", 1:ncol(hcp.factors), sep = "")
#hcp.weights <- hcp.covariates$B %>% t

#hcp.model <- data.frame(Status = factor(lumi.combat$Status), hcp.factors)

#hcp.pvals <- map(colnames(hcp.factors), hcp.lm, hcp.model)
#hcp.lm <- function(hcp.factor, hcp.model)
#{
    #lm.formula <- paste(hcp.factor, "~", "Status") %>% as.formula
    #factor.lm <- lm(lm.formula, hcp.model) %>% anova
    #return(factor.lm$'Pr(>F)')
#}

#Plot PEER weights
hcp.weights.sums <- colSums(abs(hcp.weights)) %>% data.frame #Get sums of weights for each factor
hcp.weights.sums$Factor <- 1:nrow(hcp.weights.sums) #Add column to label by factor
colnames(hcp.weights.sums)[1] <- "Weight"

#ggplot of factor versus weight
p <- ggplot(hcp.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
p <- p + theme_bw() + xlab("Factor") + ylab("Weight")
CairoPDF("./hcp_weights", height = 4, width = 6)
print(p)
dev.off()

source("../common_functions.R")
targets1.gaa <- select(pData(lumi.combat), Sample.Name, GAA1) %>% filter(!is.na(GAA1)) #Get GAA1 value for those patients who have it
cor.gaa <- gen.cor(model.PEER_covariate, targets1.gaa) #Correlate to PEER factors

targets1.onset <- select(pData(lumi.combat), Sample.Name, Onset) %>% filter(!is.na(Onset)) #Get age of onset for those patients who have it
cor.onset <- gen.cor(model.PEER_covariate, targets1.onset) #Correlate to PEER factors

PEER.traits.all <- cbind(cor.gaa, cor.onset) %>% data.frame #Make table of correlation values
PEER.traits.pval <- select(PEER.traits.all, contains("p.value")) %>% as.matrix #Extract p values
PEER.traits.cor <- select(PEER.traits.all, -contains("p.value")) %>% as.matrix #Extract r-squared values

text.matrix.PEER <- paste(signif(PEER.traits.cor, 2), '\n(', signif(PEER.traits.pval, 1), ')', sep = '') #make text matrix
dim(text.matrix.PEER) <- dim(PEER.traits.cor)
gen.text.heatmap(PEER.traits.cor, text.matrix.PEER, colnames(PEER.traits.cor), rownames(PEER.traits.cor), "", "PEER factor-trait relationships") #Make labeled heatmap of correlations of PEER factors to traits

PEER.trait.out <- data.frame(Factor = rownames(PEER.traits.cor), PEER.traits.cor, PEER.traits.pval) #Write correlations to a table
write_csv(PEER.trait.out, "PEER_trait_cor.csv")

#Removing effects of covariates + PEER factors  !!DO NOT USE FOR LINEAR MODELING WITH CONTRASTS!!
model.cov <- cbind(Male = model.sex.reduce, Age = as.numeric(lumi.combat$Draw.Age), RIN = lumi.combat$RIN) #Create model matrix of covariates to be removed
model.full.cov <- cbind(model.cov, model.PEER_covariate) #Add PEER factors
rmcov.expr <- removeBatchEffect(exprs(lumi.combat), covariates = model.cov) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
rmcov.lumi <- lumi.combat #Make a copy of lumi object
exprs(rmcov.lumi) <- rmcov.expr #Transfer cleaned expression values into new lumi object
saveRDS.gz(rmcov.lumi, file = "./save/rmcov.lumi.rda")
gen.boxplot("baseline_intensity_corrected.jpg", rmcov.lumi, batch.colors, "Covariate-corrected intensity", "Intensity") #Replot intensities to make sure they look okay

export.expr <- removeBatchEffect(exprs(lumi.combat), covariates = model.full.cov) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
export.lumi <- lumi.combat #Make a copy of lumi object
exprs(export.lumi) <- export.expr #Transfer cleaned expression values into new lumi object
saveRDS.gz(export.lumi, file = "./save/export.lumi.rda")

#Collapse the data by symbol
fdata <- fData(rmcov.lumi) #Get featureData from lumi object
rmcov.lumi <- rmcov.lumi[!is.na(fdata$SYMBOL),] #Remove the probes with missing symbols (this isn't supposed to be necessary!)
fdata <- fData(rmcov.lumi) #Extract featureData again
pdata <- pData(rmcov.lumi) #Extract phenoData
expr.collapse <- collapseRows(exprs(rmcov.lumi), factor(fdata$SYMBOL), rownames(rmcov.lumi)) #collapseRows by symbol
saveRDS.gz(expr.collapse, file = "./save/expr.collapse.rda")

model.design <- cbind(model.status) #Make covariate matrix for limma
model.full <- cbind(model.design, model.PEER_covariate) #Add PEER factors
saveRDS.gz(model.full, file = "./save/model.full.rda")

model.collinear <- data.frame(Status = factor(lumi.combat$Status), model.full.cov)

age.all <- lm(Age ~ Status, model.collinear) %>% anova
age.aov <- aov(Age ~ Status, model.collinear) %>% TukeyHSD
sex.all <- lm(Male ~ Status, model.collinear) %>% anova
sex.aov <- aov(Male ~ Status, model.collinear) %>% TukeyHSD

#Block correlation is being ignored for now
block.correlation <- readRDS.gz("./save/block.correlation.rda") #Run on hoffman
#fit.object.block <- gen.fit.block(expr.collapse$datETcollapsed, model.full, block.correlation$consensus, lumi.combat$Family)

#Fit limma linear model 
fit.object <- gen.fit(expr.collapse$datETcollapsed, model.full)
saveRDS.gz(fit.object, file = "./save/fit.object.rda")

#Calculate ratios for use in tables
ratio.exp <- gen.ratios(expr.collapse$datETcollapsed, pdata)
saveRDS.gz(ratio.exp, file = "./save/ratio.exp.rda")
log.ratios <- gen.group.ratios(expr.collapse$datETcollapsed, pdata)
saveRDS.gz(log.ratios, "./save/log.ratios.rda")

#Generate statisical cutoff
decide <- list(c("fdr", 0.05), c("fdr", 0.1), c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot <- ldply(decide, gen.decide, fit.object, FALSE) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoffs
gen.decideplot("./threshold_selection", decide.plot) #Plot different significance cutoffs

decide.final <- gen.decide(c("none", 0.005), fit.object, TRUE) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute signifance cutoff p < 0.001, no FDR adjustment
gen.decideplot("./selected_threshold", decide.final) #Plot cutoff

decide.final.fdr <- gen.decide(c("fdr", 0.1), fit.object, TRUE) %>% melt(id.vars = c("Test", "Num", "Direction")) #Compute significance cutoff p < 0.1, FDR adjusted
gen.decideplot("./selected_threshold_fdr", decide.final.fdr, 3, 4) #plot cutoff

#Make tables
de.object <- read_tsv("./fit_none.tsv") #Read in unadjusted fit object
de.object.fdr <- read_tsv("./fit_fdr.tsv") #Read in adjusted fit object
rownames(de.object) <- rownames(expr.collapse$datETcollapsed)
rownames(de.object.fdr) <- rownames(expr.collapse$datETcollapsed)
fit.selection <- gen.tables(de.object, lumi.combat, ratio.exp, "pLess001") #create differential expression table for unadjusted fit
fit.selection.fdr <- gen.tables(de.object.fdr, lumi.combat, ratio.exp, "fdrLess01") #create differential expression table for adjusted fit
saveRDS.gz(fit.selection, file = "./save/fit.selection.rda")
saveRDS.gz(fit.selection.fdr, file = "./save/fit.selection.fdr.rda")

#P-value histogram
gen.pval.hist("./hist_pvalue", fit.object$p.value) 

#Venn diagram
gen.venndiagram("./venn", results) 
gen.venndiagram("./venn_fdr", results.fdr)

#Anova heatmaps
clust.none <- gen.anova(fit.selection, "none")
clust.fdr <- gen.anova(fit.selection.fdr, "fdr")
clust.pca <- dist(clust.fdr[[2]]) %>% (flashClust::hclust)

#Code for getting size of objects in memory
objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

#Submit genes to Enrichr
source('../../code/GO/enrichr.R')
enrichr.nofdr <- select(fit.selection, Symbol, Res.time1.carrier_vs_time1.control, Res.time1.patient_vs_time1.control, Res.time1.patient_vs_time1.carrier)
enrichr.fdr <- select(fit.selection.fdr, Symbol, Res.time1.carrier_vs_time1.control, Res.time1.patient_vs_time1.control, Res.time1.patient_vs_time1.carrier)

enrichr.terms <- c("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "HumanCyc_2016", "NCI-Nature_2016", "Panther_2016") 

trap1 <- enrichr.submit("Res.time1.patient_vs_time1.control", enrichr.fdr, enrichr.terms, "fdr")
comparison.cols <- names(enrichr.nofdr[-1])
trap2 <- map(comparison.cols, enrichr.submit, enrichr.nofdr, enrichr.terms, "nofdr")

##PCA Analysis of differentially expressed genes
#res.cols <- str_subset(colnames(fit.selection), "Res")
#mds.none <- lapply(res.cols, subset.pca, fit.selection, export.lumi,  "none")
#mds.fdr <- lapply(res.cols, subset.pca, fit.selection.fdr, export.lumi, "fdr")


#### NOTE: These are GO plots for this specific data.  The tables used and rows chosen were hand selected and must be updated for any new data
#GO plots
pca.gobiol <- read.xlsx("./enrichr/nofdr/patient_vs_carrier/enrichr/GO_Biological_Process.xlsx") %>% select(Term, P.value, Genes) %>% slice(c(1, 6))
pca.gobiol$Database <- "GO Biological Process"
pca.gomole <- read.xlsx("./enrichr/nofdr/patient_vs_carrier/enrichr/GO_Molecular_Function.xlsx") %>% select(Term, P.value, Genes) %>% slice(2)
pca.gomole$Database <- "GO Molecular Process"
pca.reactome <- read.xlsx("./enrichr/nofdr/patient_vs_carrier/enrichr/Reactome_2016.xlsx") %>% select(Term, P.value, Genes) %>% slice(1)
pca.reactome$Database <- "Reactome"
pca.enrichr <- rbind(pca.gobiol, pca.gomole, pca.reactome)

nofdr.pca <- select(enrichr.nofdr, Symbol, Res.time1.patient_vs_time1.carrier) 
gen.enrichrplot(pca.enrichr, nofdr.pca, "pca.enrichr")

pco.gobiol <- read.xlsx("./enrichr/nofdr/patient_vs_control/enrichr/GO_Biological_Process.xlsx") %>% select(Term, P.value, Genes) %>% slice(c(1,2,12))
pco.gobiol$Database <- "GO Biological Process"
pco.gomole <- read.xlsx("./enrichr/nofdr/patient_vs_control/enrichr/GO_Molecular_Function.xlsx") %>% select(Term, P.value, Genes) %>% slice(1)
pco.gomole$Database <- "GO Molecular Process"
pco.gomole <- read.xlsx("./enrichr/nofdr/patient_vs_control/enrichr/Reactome_2016.xlsx") %>% select(Term, P.value, Genes) %>% slice(1)
pco.gomole$Database <- "Reactome"
pco.enrichr <- rbind(pco.gobiol, pco.gomole)

nofdr.pco <- select(enrichr.fdr, Symbol, Res.time1.patient_vs_time1.control) 
gen.enrichrplot(pco.enrichr, nofdr.pco, "pco.enrichr") 

cc.gobiol <- read.xlsx("./enrichr/nofdr/carrier_vs_control/enrichr/GO_Biological_Process.xlsx") %>% select(Term, P.value, Genes) %>% slice(c(1,2))
cc.gobiol$Database <- "GO Biological Process"
cc.kegg <- read.xlsx("./enrichr/nofdr/carrier_vs_control/enrichr/KEGG_2016.xlsx") %>% select(Term, P.value, Genes) %>% slice(c(1))
cc.kegg$Database <- "KEGG"
cc.reactome <- read.xlsx("./enrichr/nofdr/carrier_vs_control/enrichr/Reactome_2016.xlsx") %>% select(Term, P.value, Genes) %>% slice(c(1))
cc.reactome$Database <- "Reactome"
cc.enrichr <- rbind(cc.gobiol, cc.kegg, cc.reactome)

nofdr.cc <- select(enrichr.fdr, Symbol, Res.time1.carrier_vs_time1.control) 
gen.enrichrplot(cc.enrichr, nofdr.cc, "cc.enrichr") 

