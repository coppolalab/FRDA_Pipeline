library(irr)
library(tools)
library(heatmap.plus)
library(WGCNA)
library(reshape2)
enableWGCNAThreads()

gen.cor <- function(dataset, trait.df)
{
    dataset %<>% mutate(Sample.Name = rownames(dataset))
    module.trait.cor <- bicor(select(dataset, -Sample.Name), select(trait.df, -Sample.Name), maxPOutliers = 0.05)
    module.trait.cor.pval <- corPvalueStudent(module.trait.cor, nrow(dataset))
    colnames(module.trait.cor.pval) %<>% paste(".p.value", sep = "") 
    return(cbind(module.trait.cor, module.trait.cor.pval))
}

gen.text.heatmap <- function(cor.dataset, text.matrix, x.names, y.names, maintitle, filename, zlim.plot = c(-1,1))
{
    width.dynamic <- 3 + (1 * ncol(text.matrix))
    CairoPDF(filename, width = width.dynamic, height = 10)
    par(mar = c(8, 8, 3, 3))
    labeledHeatmap(Matrix = cor.dataset, xLabels = x.names, yLabels = y.names, ySymbols = y.names, yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix, setStdMargins = F, cex.text = 0.5, zlim = zlim.plot, main = maintitle)
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
{ meanIAC <- apply(dataset, 2, mean)
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

get.kappa <- function(term.current, all.terms)
{
    map(all.terms, cbind, term.current) %>% map(kappa2) %>% map_dbl(getElement, "value")
}

get.kappa.cluster <- function(enrichr.output, gene.names, filename)
{
    num.genes <- length(gene.names)
    enrichr.list <- map(enrichr.output$Genes, str_split, ",") %>% map(getElement, 1) 
    enrichr.match <- map(enrichr.list, is.element, el = toupper(gene.names)) %>% reduce(rbind) %>% t
    rownames(enrichr.match) <- toupper(gene.names)
    colnames(enrichr.match) <- enrichr.output$Term
    enrichr.match.df <- data.frame(enrichr.match)

    enrichr.kappa <- map(enrichr.match.df, get.kappa, enrichr.match.df) %>% reduce(rbind)
    enrichr.kappa[enrichr.kappa < 0] <- 0

    rownames(enrichr.kappa) <- colnames(enrichr.kappa) <- enrichr.output$Term

    CairoPDF(str_c(filename, "heatmap", sep = "."), width = 30, height = 30)
    heatmap.plus(enrichr.kappa, col = heat.colors(40), symm = TRUE, margins = c(20,20))
    dev.off()

    kappa.dist <- dist(enrichr.kappa, method = "manhattan")
    kappa.clust <- hclust(kappa.dist, method = "average")

    CairoPDF(str_c(filename, "clust", sep = "."), height = 30, width = 30)
    plot(kappa.clust)
    dev.off()

    kappa.modules <- cutreeDynamic(kappa.clust, minClusterSize = 2, method = "tree")
    kappa.modules.TOM <- cutreeDynamic(kappa.clust, distM = TOMdist(enrichr.kappa), minClusterSize = 2, method = "hybrid")
    kappa.modules.df <- data.frame(Term = rownames(enrichr.kappa), Module = kappa.modules, Module.TOM = kappa.modules.TOM)

    enrichr.output$Module <- kappa.modules
    enrichr.output$Module.TOM <- kappa.modules.TOM
    enrichr.output %<>% select(Index:Combined.Score, Module:Module.TOM, Genes)
    
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = enrichr.output, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(enrichr.output)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    
    saveWorkbook(wb, str_c(filename, "table.xlsx", sep = "."), overwrite = TRUE) 
}

gen.enrichrplot <- function(enrichr.df, filename, plot.height = 5, plot.width = 8)
{
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
    enrichr.df$Log.pvalue <- -(log10(enrichr.df$P.value))
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
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
