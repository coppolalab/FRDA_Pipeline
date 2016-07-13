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
library(R.utils)
library(irr)
library(tools)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(vadr)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(broom)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(igraph)
library(TeachingDemos)
library(heatmap.plus)

#Reading and writing tables
library(readr)
library(openxlsx)
library(parallel)

gen.pcaplot <- function(filename, dataset, facet.bool, size.height, size.width)
{
    colnames(dataset)[2] <- "Module"
    dataset$Module %<>% str_replace("ME", "") 
    p <- ggplot(dataset, aes(x = as.numeric(x), y = value, fill = Module, color = Module)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("First Principal Component")
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(dataset$x)))
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

gen.heatmap <- function(dataset, ME.genes)
{
    color <- as.character(unique(dataset$module.colors))
    dataset %<>% select(-module.colors) %>% scale
    max.dataset <- max(abs(dataset))
    print(dim(dataset))
    CairoPDF(paste("./modules/", color, sep = ""), width = 21, height = 12)
    par(mar = c(3.5,3,2,3))
    par(oma = c(4,0,2,0))
    plotMat(dataset, zlim = c(-max.dataset, max.dataset), main = paste(color, " (", nrow(dataset), ")", sep = ""))

    ME.genes.plot <- select(ME.genes, Sample.ID, matches(str_c("ME", color)))
    p <- ggplot(ME.genes.plot, aes_string(x = "Sample.ID", y = str_c("ME", color)))
    p <- p + geom_bar(stat = "identity") + xlab("Eigengene Expression")#+ ylim(c(-6, 16)) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 2))  
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    print(p)
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

plot.cor <- function(expr.row, gaa.repeats)
{
    symbol <- expr.row["Symbol"]
    symbol.index <- match(expr.row["Symbol"], expr.row)
    probe.index <- match(expr.row["Probe_Id"], expr.row)
    expr.row.clean <- expr.row[-c(symbol.index, probe.index)]

    expr.plot <- cbind(GAA = as.numeric(gaa.repeats), Expression = as.numeric(expr.row.clean)) %>% data.frame
    cor.value <- cor(expr.plot$GAA, expr.plot$Expression)
    cor.pvalue <- corPvalueStudent(cor.value, nrow(expr.plot))

    p <- ggplot(expr.plot, aes(x = as.numeric(GAA), y = as.numeric(Expression))) + geom_point()   
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ylab("Log 2 normalized expression") + xlab("GAA Repeat Length") + stat_smooth(method = "lm", se = FALSE)
    p <- p + ggtitle(paste(symbol, '\n', "Correlation:", round(cor.value, 2), "(P < ", format(cor.pvalue, digits = 3), ")"))
    CairoPDF(file.path("./correlations", symbol), height = 6, width = 10)
    print(p)
    dev.off()
}

plot.eigencor <- function(module.color, col.name, ME.genes, pheno.vector)
{
    me.genes.subset <- select(ME.genes, one_of(as.character(module.color)))
    me.genes.subset[[col.name]] <- pheno.vector
    #me.subset.melt <- melt(me.genes.subset, id.vars = col.name) 
    me.genes.subset$Module <- str_replace(module.color, "ME", "")

    #me.subset.melt$Module %<>% as.character
    p <- ggplot(me.genes.subset, aes_string(x = col.name, y = module.color, col = "Module")) + geom_point(position = "jitter")
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene") + stat_smooth(method = "lm", se = FALSE)
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(pheno.vector)))
    p <- p + scale_color_manual(values = sort(unique(me.genes.subset$Module)))
    p <- p + stat_smooth(method = "lm", se = TRUE)
    p <- p + theme(legend.position = "none")

    filename <- paste(col.name, module.color, "eigengenes_05", sep = "_")
    CairoPDF(filename, height = 5, width = 6)
    print(p)
    dev.off()
}

gen.boxplot <- function(filename, lumi.object, maintext, ylabtext)
{
    #dataset %<>% t %>% data.frame
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Batch = lumi.object$Batch)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Batch"))
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value, fill = factor(Batch))) + geom_boxplot() + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Oxygen", width = 20 , height = 8)
}

enrichr.submit <- function(index, full.df, enrichr.terms, use.weights)
{
    dataset <- filter(full.df, module.colors == index)
    dir.create(file.path("./enrichr", index), showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    map(names(enrichr.data), enrichr.wkbk, enrichr.data, index)
}

enrichr.wkbk <- function(subindex, full.df, index)
{
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

gen.enrichrplot <- function(enrichr.df, filename, plot.height = 5, plot.width = 8)
{
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
    enrichr.df$Log.pvalue <- -(log10(enrichr.df$P.value))
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
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

gen.histogram <- function(filename, lumi.object)
{
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Status = lumi.object$Status)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Status"))

    p <- ggplot(dataset.m, aes(value, group = Sample.Status, col = factor(Status))) + geom_density() + theme_bw()
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("Histogram of VST Expression") + ylab("Density") + xlab("VST Expression") 
    CairoPDF(filename, height = 5, width = 9)
    print(p)
    dev.off()
}

gen.pca <- function(filename, dataset, targetset,variablename)
{
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$Sample.Name, factor(targetset[[variablename]]))
    colnames(target.data) <- c("Sample.Status", variablename)
    colnames(dataset.plot) <- c("Sample.Status", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename)
    CairoPDF(file = filename, height = 7, width = 7)
    print(p)
    dev.off()
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

gaa.cor <- function(lumi.object, filename)
{
    cor.gaa <- apply(exprs(lumi.object), 1, bicor, lumi.object$GAA1)
    cor.gaa.pval <- corPvalueStudent(cor.gaa, length(lumi.object$GAA1)) 
    cor.gaa.adjpval <- p.adjust(cor.gaa.pval, "fdr")

    cor.gaa.df <- cbind(cor.gaa, cor.gaa.pval, cor.gaa.adjpval) %>% data.frame
    colnames(cor.gaa.df) <- c("Correlation", "P.value", "Adj.P.value")
    cor.gaa.df$Symbol <- featureNames(lumi.object)
    cor.gaa.df %<>% arrange(P.value)
    write.xlsx(cor.gaa.df, filename)
    return(cor.gaa.df)
}

get.ppi <- function(gene.list)
{
    ppi.id <- filter(iw.hugo, is.element(HUGO, gene.list))$IW.ID
    ppi.inweb.interactors <- filter(inweb.ppi, is.element(Interactor.A, ppi.id)) %>% filter(is.element(Interactor.B, ppi.id))
    ppi.inweb.symbols <- merge(ppi.inweb.interactors, iw.hugo, by.x = "Interactor.A", by.y = "IW.ID") %>% merge(iw.hugo, by.x = "Interactor.B", by.y = "IW.ID")
    ppi.inweb.final <- select(ppi.inweb.symbols, contains("HUGO"))
    colnames(ppi.inweb.final) <- c("Symbol.A", "Symbol.B")

    ppi.biogrid <- filter(biogrid.ppi.reduce, Official.Symbol.Interactor.A %in% gene.list) %>% filter(Official.Symbol.Interactor.B %in% gene.list)
    colnames(ppi.biogrid) <- c("Symbol.A", "Symbol.B")

    ppi.combined <- rbind(ppi.inweb.final, ppi.biogrid)
    ppi.combined$unique <- str_c(ppi.combined$Symbol.A, ppi.combined$Symbol.B, sep = ".")
    ppi.unique <- filter(ppi.combined, !duplicated(unique)) %>% select(-unique)
    ppi.self <- apply(ppi.unique, 1, reduce, identical)
    ppi.unique.final <- filter(ppi.unique, !ppi.self)
    return(ppi.unique.final)
}

#PPI Plot
plot.ppi <- function(adjacency.expr, gene.list, ppi.edges, filename, prune = FALSE, clust.keep = 1, plot.width = 7, plot.height = 7)
{
    expr.adjacency <- adjacency.expr[gene.list,gene.list]
    ppi.adjacency <- matrix(0, ncol = ncol(expr.adjacency), nrow = nrow(expr.adjacency), dimnames = list(rownames(expr.adjacency), colnames(expr.adjacency)))
    ppi.adjacency[as.matrix(ppi.edges)] <- 1
    final.adjacency <- expr.adjacency * ppi.adjacency
    final.mins <- apply(final.adjacency, 2, min)
    final.maxs <- apply(final.adjacency, 2, max)
    final.scaled <- scale(final.adjacency, center = final.mins, scale = final.maxs - final.mins)
    final.scaled[is.nan(final.scaled)] <- 0

    final.igraph <- graph_from_adjacency_matrix(expr.adjacency, mode = "undirected", weighted = TRUE, diag = FALSE)
    ppi.igraph <- graph_from_adjacency_matrix(final.scaled, mode = "undirected", weighted = TRUE, diag = FALSE)

    if (prune == TRUE)
    {
        #num.edges <- map(1:vcount(ppi.igraph), incident, graph = ppi.igraph) %>% map_dbl(length)

        final.igraph <- delete.vertices(final.igraph, which(clusters(ppi.igraph)$membership != clust.keep))
        ppi.igraph <- delete.vertices(ppi.igraph, which(clusters(ppi.igraph)$membership != clust.keep))
    }

    #ppi.colors <- rainbow(length(unique(clusters(ppi.igraph)$membership)))
    #final.colors <- ppi.colors[clusters(ppi.igraph)$membership]

    final.edges <- attr(E(final.igraph), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_")
    final.df <- data.frame(Edges = final.edges, Weight = edge_attr(final.igraph, "weight"))
    if (sum(final.scaled) > 0)
    {
        ppi.edges <- attr(E(ppi.igraph), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_")
        ppi.df <- data.frame(Edges = ppi.edges, Weight = edge_attr(ppi.igraph, "weight"))

        final.filter <- is.element(final.df$Edges, ppi.df$Edges)
        final.df[final.filter,]$Weight <- ppi.df$Weight 
        final.df[!final.filter,]$Weight <- 0
        final.df$Color <- "#dddddd99"
        final.df[final.filter,]$Color <- "#0000FF99"
    }
    else
    {
        final.df <- list(Color = "#0000FF99", Weight = 0)
    }

    CairoPDF(filename, width = plot.width, height = plot.height)
    par(mar=c(0,0,0,0) + 0.5)
    plot.igraph(final.igraph, layout = layout_nicely(final.igraph), vertex.size = 35, vertex.label.degree = pi/2, vertex.label.font = 2, vertex.label.color = "black", edge.color = final.df$Color, edge.width = 5*final.df$Weight)
    dev.off()

    return(ppi.igraph)
}

lumi.import <- readRDS.gz(file = "../baseline_lumi/save/lumi.known.rda")

lumi.baseline <- lumi.import[,lumi.import$Status == "Patient" & lumi.import$Sample.Num == 1 & !is.na(lumi.import$GAA1)]
lumi.vst <- lumiT(lumi.baseline)
lumi.batch <- lumi.vst[,lumi.vst$Batch != 8 & lumi.vst$Batch != 14]
lumi.batch$Batch %<>% droplevels
saveRDS.gz(lumi.batch, file = "./save/lumi.batch.rda")

gen.boxplot("baseline_intensity_raw.jpg", lumi.batch, "VST Transformed Intensity", "Intensity")
gen.histogram("baseline_histogram", lumi.batch)
mds.baseline <- exprs(lumi.batch) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
gen.pca("baseline_mds_batch", mds.baseline, pData(lumi.batch), "Batch") #label PCs by batch

lumi.patient.norm <- lumiN(lumi.batch, method = "rsn") #Normalize with robust spline regression
lumi.patient.cutoff <- detectionCall(lumi.patient.norm) #Get the count of probes which passed the detection threshold per sample
lumi.patient.expr <- lumi.patient.norm[which(lumi.patient.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.patient <- getSYMBOL(rownames(lumi.patient.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.patient.annot <- lumi.patient.expr[!symbols.lumi.patient,] #Drop any probe which is not annotated
saveRDS.gz(lumi.patient.annot, file = "./save/lumi.patient.annot.rda")

#Use ComBat for batch effect correction
lumi.patient.annot$Sex %<>% droplevels
#lumi.age.stratify <- lumi.patient.annot[,lumi.patient.annot$Draw.Age >= 6689 & lumi.patient.annot$Draw.Age <= 13720]
#lumi.patient.annot <- lumi.age.stratify
model.combat <- model.matrix(~ GAA1 + Sex + Draw.Age + RIN, pData(lumi.patient.annot)) %>% data.frame

age.gaa <- lm(as.numeric(Draw.Age) ~ GAA1, pData(lumi.patient.annot)) %>% anova %>% tidy
sex.gaa <- lm(as.integer(Sex) ~ GAA1, pData(lumi.patient.annot)) %>% anova %>% tidy
site.gaa <- lm(as.integer(factor(Site)) ~ GAA1, pData(lumi.patient.annot)) %>% anova %>% tidy

expr.combat <- ComBat(dat = exprs(lumi.patient.annot), batch = factor(lumi.patient.annot$Site), mod = model.combat)
lumi.combat <- lumi.patient.annot
exprs(lumi.combat) <- expr.combat

expr.combat <- ComBat(dat = exprs(lumi.combat), batch = factor(lumi.combat$Batch), mod = model.combat)
exprs(lumi.combat) <- expr.combat
saveRDS.gz(lumi.combat, file = "./save/lumi.combat.rda")

gen.boxplot("baseline_intensity_combat.jpg", lumi.combat, "VST Transformed Intensity (Batch Corrected)", "Intensity")
gen.histogram("baseline_histogram_combat", lumi.combat)
mds.baseline.combat <- exprs(lumi.combat) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
gen.pca("baseline_mds_batch_combat", mds.baseline.combat, phenoData(lumi.combat), "Batch") #label PCs by batch

combat.connectivity <- gen.connectivityplot("baseline_connectivity_combat", lumi.combat, "")
connectivity.sorted <- sort(abs(combat.connectivity), decreasing = TRUE)

qc.df <- data.frame(Connectivity = connectivity.sorted, RIN = pData(lumi.combat)[names(connectivity.sorted),]$RIN)

connectivity.outlier <- combat.connectivity[abs(combat.connectivity) > 2]
outlier.names <- names(connectivity.outlier) %>% is.element(el = sampleNames(lumi.combat))

lumi.rmout <- lumi.combat[,!outlier.names]
rmout.connectivity <- gen.connectivityplot("rmout_connectivity", lumi.rmout, "")

gen.boxplot("baseline_intensity_rmout.jpg", lumi.rmout, "VST Transformed Intensity (Outliers Removed)", "Intensity")
gen.histogram("baseline_histogram_rmout", lumi.rmout)
mds.baseline.rmout <- exprs(lumi.rmout) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
gen.pca("baseline_mds_batch_rmout", mds.baseline.rmout, phenoData(lumi.rmout), "Batch") #label PCs by batch

model.rmout <- model.matrix(~ Sex + Draw.Age + RIN, pData(lumi.rmout)) %>% data.frame
cleaned.expr <- removeBatchEffect(exprs(lumi.rmout), covariates = model.rmout[,2:4])
lumi.cleaned <- lumi.rmout
exprs(lumi.cleaned) <- cleaned.expr
saveRDS.gz(lumi.cleaned, file = "./save/lumi.cleaned.rda")

#batch.colors <- c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green")
#gen.boxplot("baseline_intensity_corrected.jpg", lumi.cleaned, batch.colors, "Covariate-corrected intensity", "Intensity")

expr.collapse <- collapseRows(exprs(lumi.cleaned), getSYMBOL(featureNames(lumi.cleaned), 'lumiHumanAll.db'), rownames(lumi.cleaned))$datETcollapsed
export.lumi <- ExpressionSet(assayData = expr.collapse, phenoData = phenoData(lumi.cleaned))
saveRDS.gz(expr.collapse, "./save/expr.collapse.rda")
saveRDS.gz(export.lumi, "./save/export.lumi.rda")

lumi.short <- export.lumi[,export.lumi$GAA1 <= 466]
lumi.intermediate <- export.lumi[,export.lumi$GAA1 <= 800]
lumi.long <- export.lumi[,export.lumi$GAA1 > 800]

full.cor.df <- gaa.cor(export.lumi, "gaa.cor.xlsx")
short.cor.df <- gaa.cor(lumi.short, "gaa.cor.short.xlsx") 
long.cor.df <- gaa.cor(lumi.long, "gaa.cor.long.xlsx")
intermediate.cor.df <- gaa.cor(lumi.intermediate, "gaa.cor.intermediate.xlsx")

#Correlations
#cor.enrichr <- map(enrichr.terms, get.enrichrdata, cor.df.sg, FALSE)
#names(cor.enrichr) <- enrichr.terms
#map(names(cor.enrichr), enrichr.wkbk, cor.enrichr, "gaa_cor")

#get.updown <- function(filter.vector, enrichr.df)
#{
    #grep.vector <- str_replace_all(filter.vector, ",", "|")
    #enrichr.filter <- filter(enrichr.df, grepl(grep.vector, Symbol))
    #num.up <- which(enrichr.filter$Correlation > 0) %>% length
    #num.down <- which(enrichr.filter$Correlation < 0) %>% length
    #return(c("Up" = num.up, "Down" = num.down))
#}

#gaa.gobiol <- read.xlsx("./enrichr/gaa_cor/gaa_cor_GO_Biological_Process.xlsx") %>% select(GO.Term, P.value, Genes) %>% slice(c(1, 8, 9, 18))
#gaa.gobiol$Database <- "GO Biological Process"
#gaa.gomole <- read.xlsx("./enrichr/gaa_cor/gaa_cor_GO_Molecular_Function.xlsx") %>% select(GO.Term, P.value, Genes) %>% slice(c(1, 8))
#gaa.gomole$Database <- "GO Molecular Process"
#gaa.reactome <- read.xlsx("./enrichr/gaa_cor/gaa_cor_Reactome_2015.xlsx") %>% select(GO.Term, P.value, Genes) %>% slice(c(23, 25, 34))
#gaa.reactome$Database <- "Reactome"
#gaa.enrichr <- rbind(gaa.gobiol, gaa.gomole, gaa.reactome)

#gaa.enrichr$Gene.Count <- map(gaa.enrichr$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
#gaa.enrichr$Log.pvalue <- -(log10(gaa.enrichr$P.value))

#gaa.updown <- map(gaa.enrichr$Genes, get.updown, cor.df.sg) %>% reduce(rbind)
##colnames(gaa.updown) <- c("Up", "Down")
#gaa.enrichr <- cbind(gaa.enrichr, gaa.updown)
#gaa.enrichr$Log.Up <- gaa.enrichr$Log.pvalue * gaa.enrichr$Up / gaa.enrichr$Gene.Count
#gaa.enrichr$Log.Down <- gaa.enrichr$Log.pvalue * gaa.enrichr$Down / gaa.enrichr$Gene.Count
#gaa.enrichr$GO.Term %<>% str_replace_all("\\ \\(.*$", "") %>% tolower
#gaa.enrichr$Format.Name <- paste(gaa.enrichr$Database, ": ", gaa.enrichr$GO.Term, " (", gaa.enrichr$Gene.Count, ")", sep = "")
#gaa.enrichr.plot <- select(gaa.enrichr, Format.Name, Log.Up, Log.Down) %>% melt(id.vars = "Format.Name") 

#p <- ggplot(gaa.enrichr.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = c(gaa.enrichr$Format.Name, rep("", length(gaa.enrichr$Format.Name))), hjust = "left", aes(y = 0.1))
#p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "FALSE",  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
#CairoPDF("gaa.enrichr", height = 5, width = 8)
#print(p)
#dev.off()

lumi.top <- export.lumi[,export.lumi$GAA1 <= 800]

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.data <- t(exprs(export.lumi))
saveRDS.gz(expr.data, file = "./save/expr.data.rda")
sft <- pickSoftThreshold(expr.data, powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.df <- sft$fitIndices
saveRDS.gz(sft, file = "./save/sft.peer.rda")

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
saveRDS.gz(adjacency.expr, file = "./save/adjacency.expr.rda")

TOM <- TOMsimilarity(adjacency.expr)
dissimilarity.TOM <- 1 - TOM
saveRDS.gz(TOM, file = "./save/tom.rda")
saveRDS.gz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 20

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.data, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.20
merge.all <- mergeCloseModules(expr.data, dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) #PC analysis may be failing because of Intel MKL Lapack routine bug.  Test with openBLAS in R compiled with gcc.
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
saveRDS.gz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
saveRDS.gz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
saveRDS.gz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 10)
par(cex = 0.7)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.adjacency = 0.3, cex.preservation = 0.3, plotPreservation = "standard")
dev.off()

all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Symbol = rownames(all.degrees), module.color = module.colors, all.degrees)

write_csv(data.frame(table(module.colors)), path = "./final_eigengenes.csv") 
gene.info$kscaled <- by(gene.info, gene.info$module.color, select, kWithin) %>% llply(function(x) { x / max (x) }) %>% reduce(c)
saveRDS.gz(gene.info, file = "./save/gene.info.rda")

gene.module.membership <- as.data.frame(bicor(expr.data, ME.genes, maxPOutliers = 0.05))
module.membership.pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.data)))
names(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")

module.membership <- cbind(select(gene.info, Symbol:module.color), gene.module.membership, module.membership.pvalue)
write_csv(module.membership, "module_membership.csv")

colnames(gene.module.membership) %<>% str_replace("MM.", "")
colnames(module.membership.pvalue) %<>% str_replace("MM.pvalue.", "")
gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

color.values <- unique(module.colors)
gene.module.membership.long <- gather_(data = gene.module.membership, "module.comparison", "correlation", color.values)
module.membership.pvalue.long <- gather_(data = module.membership.pvalue, "module.comparison", "p.value", color.values)
membership.join <- join(gene.module.membership.long, module.membership.pvalue.long)
eigengene.connectivity <- join(membership.join, gene.info) %>% select(Symbol, module.color:kscaled, module.comparison:p.value)
write_csv(eigengene.connectivity, "eigengene_connectivity.csv")

all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y") 
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$x <- as.factor(1:nrow(smooth.df))
smooth.plot <- melt(smooth.df, id.vars = "x")

gen.pcaplot("all_principal_components", smooth.plot, FALSE, 10, 15)
gen.pcaplot("facet_principal_components", smooth.plot, TRUE, 13, 25)

sample.ids <- factor(rownames(expr.data), levels = rownames(expr.data))
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.data), module.colors)
cluster <- makeForkCluster(8)
split(expr.data.plot, expr.data.plot$module.colors) %>% map(gen.heatmap, ME.genes.plot) 
modules.out <- select(gene.info, Symbol, module.color)
write.xlsx(modules.out, "modules_out.xlsx")

source("../../code/GO/enrichr.R")
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "HumanCyc_2016", "NCI-Nature_2016", "Panther_2016")
enrichr.submit("pink", modules.out, enrichr.terms, FALSE)
color.names <- unique(module.colors) %>% sort
l_ply(color.names, enrichr.submit, modules.out, enrichr.terms, FALSE)

#modules.out.reduce <- filter(modules.out, module.color != "grey")
#modules.out.reduce$module.color %<>% droplevels

#split(modules.out.reduce, modules.out.reduce$module.color) %>% map(submit.stringdb)

#submit.stringdb <- function(module.subset)
#{
    #get.stringdb(module.subset, unique(module.subset$module.color), "./stringdb")
#}

connectivity.reduce <- filter(eigengene.connectivity, module.color != "grey" & module.comparison != "grey" & module.color == module.comparison)

targets.final.gaa <- pData(export.lumi)
rownames(ME.genes) <- rownames(expr.data)

eigengene.model <- function(ME.vector, trait.vector)
{
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.aov <- aov(ME ~ Trait, trait.df) %>% TukeyHSD
    trait.tukey <- data.frame(trait.aov$Trait)
    colnames(trait.tukey) %<>% str_replace(" ", ".")
    return(select(trait.tukey, diff, p.adj))
}

eigengene.anova <- function(ME.vector, trait.vector)
{
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- lm(ME ~ Trait, trait.df) %>% anova
    return(trait.anova$'Pr(>F)'[1])
}

targets.final.gaa$Repeat.Size <- factor(targets.final.gaa$GAA1 > 640)
targets.final.gaa$Repeat.Size %<>% revalue(c("TRUE" = "Long", "FALSE" = "Short"))
cor.anova <- map(ME.genes, eigengene.anova, targets.final.gaa$Repeat.Size)
cor.status <- map(ME.genes, eigengene.model, targets.final.gaa$Repeat.Size)
status.diff <- map(cor.status, select, diff) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.diff) <- names(cor.status)
status.pval <- map(cor.status, select, p.adj) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.pval) <- names(cor.status)

pval.adjust <- map(status.pval, p.adjust, method = "fdr") %>% reduce(cbind) %>% data.frame
rownames(pval.adjust) <- rownames(status.pval)
colnames(pval.adjust) <- paste(colnames(status.pval), ".pval")

text.matrix.traits <- paste(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(status.pval), 1), ')', sep = '')
dim(text.matrix.traits) = dim(status.diff)

status.pval$Module <- rownames(status.pval)
heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
gen.text.heatmap(as.matrix(status.diff), text.matrix.traits, colnames(status.diff), colnames(ME.genes), "", "module-trait relationships", heatmap.range)

source("../common_functions.R")
targets.gaa <- select(targets.final.gaa, Sample.Name, GAA1) 
cor.gaa <- gen.cor(ME.genes, targets.gaa)

module.traits.all <- data.frame(cor.gaa)
module.traits.pval <- select(module.traits.all, contains("p.value")) %>% as.matrix %>% apply(2, p.adjust, "fdr")
module.traits.cor <- select(module.traits.all, -contains("p.value")) %>% as.matrix

module.trait.out <- data.frame(Module = rownames(module.traits.cor), module.traits.cor, module.traits.pval)

write_csv(module.trait.out, "module_trait_cor.csv")

text.matrix.traits <- paste(signif(module.traits.cor, 2), '\n(', signif(module.traits.pval, 1), ')', sep = '')
dim(text.matrix.traits) = dim(module.traits.cor)

gen.text.heatmap(module.traits.cor, text.matrix.traits, colnames(module.traits.cor), colnames(ME.genes), "", "module-trait relationships")

#ME.genes.plot <- ME.genes
#colnames(ME.genes.plot) %<>% str_replace_all("ME", "")
plot.eigencor("MEblack", "GAA1", ME.genes.plot, lumi.top$GAA1)
plot.eigencor("MEpink", "GAA1", ME.genes.plot, lumi.top$GAA1)
plot.eigencor("MEblue", "GAA1", ME.genes.plot, lumi.top$GAA1)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

black.only <- filter(modules.out, module.color == "black")
black.gobiol.file <- "./enrichr/black/black_GO_Biological_Process_2015.xlsx"
black.gobiol <- read.xlsx(black.gobiol.file) 
black.gobiol$Num.Genes <- map(black.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
black.gobiol$Database <- "GO Biological Process"
get.kappa.cluster(black.gobiol, black.only$Symbol, file_path_sans_ext(black.gobiol.file))
black.gobiol.final <- slice(black.gobiol, c(8, 23, 6, 14, 13, 24))

black.gomole.file <- "./enrichr/black/black_GO_Molecular_Function_2015.xlsx"
black.gomole <- read.xlsx(black.gomole.file) 
black.gomole$Num.Genes <- map(black.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.gomole %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
black.gomole$Database <- "GO Molecular Function"
get.kappa.cluster(black.gomole, black.only$Symbol, file_path_sans_ext(black.gomole.file))
black.gomole.final <- slice(black.gomole, c(10))

black.reactome.file <- "./enrichr/black/black_Reactome_2016.xlsx"
black.reactome <- read.xlsx(black.reactome.file) 
black.reactome$Num.Genes <- map(black.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.reactome %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
black.reactome$Database <- "Reactome"
get.kappa.cluster(black.reactome, black.only$Symbol, file_path_sans_ext(black.reactome.file))
black.reactome.final <- slice(black.reactome, c(28, 39))

black.kegg <- read.xlsx("./enrichr/black/black_KEGG_2016.xlsx") %>% slice(39)
black.kegg$Num.Genes <- map(black.kegg$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.kegg$Database <- "KEGG"

black.enrichr <- rbind(black.gobiol.final, black.gomole.final, black.reactome.final, black.kegg)
gen.enrichrplot(black.enrichr, "black.enrichr")

#get kegg

blue.only <- filter(modules.out, module.color == "blue")
blue.gobiol.file <- "./enrichr/blue/blue_GO_Biological_Process_2015.xlsx"
blue.gobiol <- read.xlsx(blue.gobiol.file) 
blue.gobiol$Num.Genes <- map(blue.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
blue.gobiol$Database <- "GO Biological Process"
blue.gobiol.filter <- filter(blue.gobiol, Num.Genes > 4) %>% filter(P.value < 0.01)
get.kappa.cluster(blue.gobiol.filter, blue.only$Symbol, file_path_sans_ext(blue.gobiol.file))
blue.gobiol.final <- slice(blue.gobiol, c(19, 7, 6))

blue.gomole.file <- "./enrichr/blue/blue_GO_Molecular_Function_2015.xlsx"
blue.gomole <- read.xlsx(blue.gomole.file) 
blue.gomole$Num.Genes <- map(blue.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
blue.gomole$Database <- "GO Molecular Function"
blue.gomole.filter <- filter(blue.gomole, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(blue.gomole.filter, blue.only$Symbol, file_path_sans_ext(blue.gomole.file))
blue.gomole.final <- slice(blue.gomole, c(5, 13, 3, 15))

blue.reactome.file <- "./enrichr/blue/blue_Reactome_2016.xlsx"
blue.reactome <- read.xlsx(blue.reactome.file) 
blue.reactome$Num.Genes <- map(blue.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
blue.reactome$Database <- "Reactome"
blue.reactome.filter <- filter(blue.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(blue.reactome.filter, blue.only$Symbol, file_path_sans_ext(blue.reactome.file))
blue.reactome.final <- slice(blue.reactome, c(7, 1, 24, 38))

#get kegg

blue.enrichr <- rbind(blue.gobiol.final, blue.gomole.final, blue.reactome.final)
gen.enrichrplot(blue.enrichr, "blue.enrichr")

pink.only <- filter(modules.out, module.color == "pink")
pink.gobiol.file <- "./enrichr/pink/pink_GO_Biological_Process_2015.xlsx"
pink.gobiol <- read.xlsx(pink.gobiol.file) 
pink.gobiol$Num.Genes <- map(pink.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pink.gobiol$Database <- "GO Biological Process"
pink.gobiol.final <- filter(pink.gobiol, Num.Genes > 4) %>% filter(P.value < 0.01)
get.kappa.cluster(pink.gobiol.final, pink.only$Symbol, file_path_sans_ext(pink.gobiol.file))
pink.gobiol.final <- slice(pink.gobiol, c(1, 72, 116))

pink.gomole.file <- "./enrichr/pink/pink_GO_Molecular_Function_2015.xlsx"
pink.gomole <- read.xlsx(pink.gomole.file) 
pink.gomole$Num.Genes <- map(pink.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pink.gomole$Database <- "GO Molecular Function"
pink.gomole.filter <- filter(pink.gomole, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(pink.gomole.filter, pink.only$Symbol, file_path_sans_ext(pink.gomole.file))
pink.molec.final <- slice(pink.gomole, c(1))

#pink.reactome.file <- "./enrichr/pink/pink_Reactome_2016.xlsx"
#pink.reactome <- read.xlsx(pink.reactome.file) 
#pink.reactome$Num.Genes <- map(pink.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
#pink.reactome$Database <- "Reactome"
#pink.reactome.final <- filter(pink.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
#get.kappa.cluster(pink.reactome.final, pink.only$Symbol, file_path_sans_ext(pink.reactome.file))

pink.kegg.file <- "./enrichr/pink/pink_KEGG_2016.xlsx"
pink.kegg <- read.xlsx(pink.kegg.file) 
pink.kegg$Num.Genes <- map(pink.kegg$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pink.kegg$Database <- "kegg"
pink.kegg.final <- slice(pink.kegg, c(5, 10))

pink.final <- rbind(pink.gobiol.final, pink.molec.final, pink.kegg.final)
gen.enrichrplot(pink.final, "pink.enrichr")

#PPI - fix me!
biogrid.ppi <- read_tsv("../../code/BIOGRID-ORGANISM-3.4.136.tab2/BIOGRID-ORGANISM-Homo_sapiens-3.4.136.tab2.txt") %>% data.frame
biogrid.ppi.reduce <- select(biogrid.ppi, contains("Official"))

inweb.ppi <- read_tsv("../../Dementia Project/mapt/InWeb3_HC_NonRed.txt", col_names = FALSE)
colnames(inweb.ppi) <- c("Interactor.A", "Interactor.B")
iw.hugo <- read_tsv("../../Dementia Project/mapt/IWtoHugo.txt", col_names = FALSE)
colnames(iw.hugo) <- c("IW.ID", "HUGO")

yellow.kegg.als <- str_split(yellow.kegg$Genes, ",")[[1]]
yellow.reactome.oxid <- str_split(yellow.reactome$Genes, ",")[[1]]
yellow.gobiol.tetra <- str_split(yellow.gobiol.final[2,]$Genes, ",")[[1]]

brown.gobiol.acet <- str_split(brown.gobiol.final[1,]$Genes, ",")[[1]]
brown.reactome.tca <- str_split(brown.reactome.final$Genes, ",")[[1]]
brown.gomole.iron <- str_split(brown.gomole.final[1,]$Genes, ",")[[1]]

yka.ppi <- get.ppi(yellow.kegg.als) #Maybe don't prune
yro.ppi <- get.ppi(yellow.reactome.oxid) 
ygt.ppi <- get.ppi(yellow.gobiol.tetra) #Don't prune

bga.ppi <- get.ppi(brown.gobiol.acet) #Don't prune
brt.ppi <- get.ppi(brown.reactome.tca) #Hueg
bgi.ppi <- get.ppi(brown.gomole.iron) #Don't prune?

adjacency.expr <- readRDS.gz("./save/adjacency.expr.rda")

plot.ppi(adjacency.expr, yellow.kegg.als, yka.ppi, "yka_igraph")
yro.igraph <- plot.ppi(adjacency.expr, yellow.reactome.oxid, yro.ppi, "yro_igraph", prune = TRUE, clust.keep = 2)
yro.incident <- map(1:vcount(yro.igraph), incident, graph = yro.igraph) %>% map_dbl(length)
names(yro.incident) <- V(yro.igraph)$name
plot.ppi(adjacency.expr, yellow.gobiol.tetra, ygt.ppi, "ygt_igraph")

plot.ppi(adjacency.expr, brown.gobiol.acet, bga.ppi, "bga_igraph")
brt.igraph <- plot.ppi(adjacency.expr, brown.reactome.tca, brt.ppi, "brt_igraph", TRUE)
brt.incident <- map(1:vcount(brt.igraph), incident, graph = brt.igraph) %>% map_dbl(length)
names(brt.incident) <- V(brt.igraph)$name
plot.ppi(adjacency.expr, brown.gomole.iron, bgi.ppi, "bgi_igraph")

