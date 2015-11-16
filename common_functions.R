library(peer)
library(WGCNA)
enableWGCNAThreads()

gen.cor <- function(dataset, trait.df)
{
    dataset %<>% mutate(Sample.Name = rownames(dataset))
    sample.key <- paste(trait.df$Sample.Name, collapse = "|")
    dataset.reduce <- filter(dataset, grepl(sample.key, Sample.Name))
    module.trait.cor <- cor(select(dataset.reduce, -Sample.Name), select(trait.df, -Sample.Name), use = "p")
    module.trait.cor.pval <- corPvalueStudent(module.trait.cor, nrow(dataset.reduce))
    colnames(module.trait.cor.pval) %<>% paste(".p.value", sep = "") 
    return(cbind(module.trait.cor, module.trait.cor.pval))
}

gen.text.heatmap <- function(cor.dataset, text.matrix, x.names, y.names, maintitle, filename)
{
    width.dynamic <- 3 + (1 * ncol(text.matrix))
    CairoPDF(filename, width = width.dynamic, height = 10)
    par(mar = c(8, 8, 3, 3))
    labeledHeatmap(Matrix = cor.dataset, xLabels = x.names, yLabels = y.names, ySymbols = y.names, yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix, setStdMargins = F, cex.text = 0.5, zlim = c(-1,1), main = maintitle)
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

gen.peer <- function(num.factors, intensities, use.covariates, covariates)
{
    model = PEER()
    PEER_setNk(model,num.factors)
    PEER_setPhenoMean(model, as.matrix(t(intensities)))
    PEER_setAdd_mean(model, TRUE)
    if (use.covariates == TRUE)
    {
        PEER_setCovariates(model, as.matrix(covariates))
    }
    PEER_setNmax_iterations(model, 1000)
    PEER_update(model)
    residuals.PEER = t(PEER_getResiduals(model))
    rownames(residuals.PEER) = rownames(intensities)
    colnames(residuals.PEER) = colnames(intensities)

    write_csv(data.frame(residuals.PEER), path = paste("residuals_", num.factors, sep = "", ".csv"))
    write_csv(data.frame(PEER_getX(model)), path = paste("factor_", num.factors, sep = "", ".csv"))
    write_csv(data.frame(PEER_getW(model)), path = paste("weight_", num.factors, sep = "", ".csv"))
    write_csv(data.frame(PEER_getAlpha(model)), path = paste("precision_", num.factors, sep = "", ".csv"))

    CairoPDF(file = paste("model", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    PEER_plotModel(model)
    dev.off()

    CairoPDF(file = paste("precision_", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    plot(PEER_getAlpha(model), col = "red", lwd = 4, main = paste("precision", num.factors, "factor", sep = " "))
    dev.off()
}

