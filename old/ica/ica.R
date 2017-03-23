#String operations
library(stringr)
library(tools)

#Reading and writing tables
library(readr)
library(openxlsx)

#For plotting
library(ggplot2)
library(Cairo)
library(heatmap.plus)

#For microarray stuff
library(Biobase)
library(matrixStats)
library(abind)
library(lumi)
library(lumiHumanAll.db)
library(annotate)
library(WGCNA)
library(sva)
library(limma)
library(irr)

#Longitudinal Analysis
library(fastICA)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(vadr)
library(rlist)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)

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

extract.ica <- function(index, dataset)
{
    icm.sig <- which(abs(dataset[index]) > 3)
    dataset.sig <- dataset[icm.sig,] #%>% arrange_(index)
    dataset.sig.sorted <- dataset.sig[order(abs(dataset.sig[,index]), decreasing = TRUE),]
    return(as.numeric(rownames(dataset.sig.sorted)))
}

gen.ica <- function(dataset)
{
    rownames(dataset) <- dataset$Symbol
    dataset.abs <- abs(select(dataset, -Symbol)) 
    ica.subset <- apply((dataset.abs > 3), 1, any)
    return(dataset.abs[ica.subset,])
}

split.ica <- function(data.vector)
{
    dataset <- data.frame(Symbol = names(data.vector), Module = data.vector)
    dataset.split <- split(dataset, dataset$Module) 
    names(dataset.split) <- paste("X", 1:length(dataset.split), sep = "")
    return(dataset.split)
}

enrichr.ica <- function(status.index, data.raw, enrichr.terms, prefix)
{
    data.subset <- data.raw[[status.index]]
    full.path <- file.path("./enrichr", prefix, status.index)
    dir.create(full.path, showWarnings = FALSE, recursive = TRUE)
    names(data.subset) <- paste("X", 1:length(data.subset), sep = "")
    lapply(names(data.subset), enrichr.submit, data.subset, full.path, enrichr.terms, FALSE)
}

enrichr.submit <- function(index, full.df, full.path, enrichr.terms, use.weights)
{
    dataset <- full.df[[index]]
    new.fullpath <- file.path(full.path, index)
    dir.create(new.fullpath, showWarnings = FALSE, recursive = TRUE)
    enrichr.data <- lapply(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    lapply(names(enrichr.data), enrichr.wkbk, enrichr.data, new.fullpath)
}

enrichr.wkbk <- function(subindex, full.df, new.fullpath)
{
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    filename = paste(file.path(new.fullpath, subindex), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

gen.icaplots <- function(index, dataset)
{
    data.subset <- dataset[[index]]
    names(data.subset) <- paste("ICM", 1:length(data.subset))
    dataset.normalized <- lapply(data.subset, normalize.data)

    dataset.melt <- melt(dataset.normalized, id.vars = "Symbol")
    dataset.melt$variable %<>% str_replace("X", "")
    p <- ggplot(dataset.melt, aes(x = variable, y = value, group = Symbol, color = Symbol)) + geom_line(size = 1.5) 
    p <- p + theme_bw() + theme(legend.position = "none") + ylab("Normalized log2 Expression") + xlab("Time") + facet_wrap(~ L1, scale = "free_y")
    CairoPDF(index, width = 12, height = 8)
    print(p)
    dev.off()
}

normalize.data <- function(dataset)
{
    rownames(dataset) <- dataset$Symbol
    dataset %<>% select(-Symbol)
    dataset.normalized <- sweep(dataset, 1, dataset$X1)
    dataset.normalized$Symbol <- rownames(dataset.normalized) 
    return(dataset.normalized[1:10,])
}

ica.wkbk <- function(index, dataset, new.fullpath)
{
    data.subset <- dataset[[index]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = data.subset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1:3, 5:ncol(data.subset)), widths = "auto")
    setColWidths(wb, 1, cols = 4, widths = 45)

    filename = paste(file.path(new.fullpath, index), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

gen.icatables <- function(index, dataset, prefix)
{
    data.subset <- dataset[[index]]
    full.path <- file.path("./modules", prefix, index)
    dir.create(full.path, showWarnings = FALSE, recursive = TRUE)
    names(data.subset) <- paste("ICM", 1:length(data.subset), sep = "")
    print(str(data.subset))
    l_ply(names(data.subset), ica.wkbk, data.subset, full.path)
}

stringdb.ica <- function(status.index, data.raw, prefix)
{
    data.subset <- data.raw[[status.index]]
    full.path <- file.path("./stringdb", prefix, status.index)
    dir.create(full.path, showWarnings = FALSE, recursive = TRUE)
    names(data.subset) <- paste("X", 1:length(data.subset), sep = "")
    map(names(data.subset), stringdb.submit, data.subset, full.path)
}

stringdb.submit <- function(component, all.components, full.path)
{
    dataset <- all.components[[component]]
    get.stringdb(dataset, component, full.path, 400)
}

setnames <- function(dataset)
{
    rownames(dataset) <- dataset$Symbol
    return(dataset)
}

seed.ICA <- function(intensities, ica.list = list(), iter.count = 0)
{
   if (iter.count < 250) 
   {
        ica.new <- map(intensities, fastICA, 4, "deflation") %>% map(`[[`, "S") %>% map(data.frame, Symbol = rownames(intensities[[1]]))
        iter.count <- iter.count + 1
        ica.list <- c(ica.list, ica.new)
        print(iter.count)
        seed.ICA(intensities, ica.list, iter.count)
   }
   else
   {
        return(ica.list)
   }
}

collapse.ica <- function(subset.key, ica.list)
{
    ica.subset <- ica.list[which(str_detect(names(ica.list), subset.key))]
    ica.melt <- melt(ica.subset)
    ica.cast <- dcast(ica.melt, Symbol ~ variable, Compose(abs, median)) #%>% select(-Symbol)
    return(ica.cast)
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

lumi.import <- readRDS.gz("../baseline_lumi/save/rmcov.collapse.rda")
lumi.four <- readRDS.gz('../dtw/save/lumi.four.rda')
lumi.patient <- lumi.four[,lumi.four$Status == "Patient"]
lumi.norm <- lumiN(lumi.patient, method = "rsn")
lumi.qual <- lumiQ(lumi.norm, detectionTh = 0.01)

lumi.cutoff <- detectionCall(lumi.qual)
lumi.expr <- lumi.qual[which(lumi.cutoff > 0),]
symbols.lumi <- getSYMBOL(rownames(lumi.expr), 'lumiHumanAll.db') %>% is.na
lumi.expr.annot <- lumi.expr[!symbols.lumi,] #Drop any probe which is not annotated
saveRDS.gz(lumi.expr.annot, file = "./save/lumi.expr.annot.rda")

qcsum <- lumi.expr.annot@QC$sampleSummary %>% t %>% data.frame
colnames(qcsum) %<>% str_replace("\\.0\\.01\\.", "")
qcsum$Sample.Name <- rownames(qcsum)
qcsum$RIN <- lumi.expr.annot$RIN
qcsum$Sample.Num <- lumi.expr.annot$Sample.Num
qcsum %<>% arrange(distance.to.sample.mean)

qcsum.remove <- filter(qcsum, distance.to.sample.mean > 70)$Sample.Name #%>% paste(collapse = "|")

connectivity.zscore <- gen.connectivityplot("connectivity", lumi.expr.annot, "")
connectivity.remove <- connectivity.zscore[abs(connectivity.zscore) > 2] %>% names

combined.remove <- c(qcsum.remove, connectivity.remove) %>% unique %>% paste(collapse = "|")
remove.indices <- grepl(combined.remove, sampleNames(lumi.expr.annot))

lumi.rmout <- lumi.patient[,!remove.indices]
lumi.rmout.norm <- lumiN(lumi.rmout, method = "rsn") #Normalize with robust spline regression
lumi.rmout.qual <- lumiQ(lumi.rmout.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.rmout.cutoff <- detectionCall(lumi.rmout.qual) #Get the count of probes which passed the detection threshold per sample
lumi.rmout.expr <- lumi.rmout.qual[which(lumi.rmout.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.rmout <- getSYMBOL(rownames(lumi.rmout.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.rmout.annot <- lumi.rmout.expr[!symbols.lumi.rmout,] #Drop any probe which is not annotated
saveRDS.gz(lumi.rmout.annot, file = "./save/lumi.rmout.annot.rda")

qcsum.rmout <- lumi.rmout.annot@QC$sampleSummary %>% t %>% data.frame
colnames(qcsum.rmout) %<>% str_replace("\\.0\\.01\\.", "")
qcsum.rmout$Sample.Name <- rownames(qcsum.rmout)
qcsum.rmout$RIN <- lumi.rmout.annot$RIN
qcsum.rmout$Sample.Num <- lumi.rmout.annot$Sample.Num
qcsum.rmout$PIDN <- lumi.rmout.annot$PIDN
arrange(qcsum.rmout, distance.to.sample.mean)

PIDNs <- filter(qcsum.rmout, Sample.Num == "1r")$PIDN 
orig <- filter(qcsum.rmout, PIDN %in% PIDNs & Sample.Num == "1")$distance.to.sample.mean
orig.names <- filter(qcsum.rmout, PIDN %in% PIDNs & Sample.Num == "1")$Sample.Name
reps <- filter(qcsum.rmout, PIDN %in% PIDNs & Sample.Num == "1r")$distance.to.sample.mean
reps.names <- filter(qcsum.rmout, PIDN %in% PIDNs & Sample.Num == "1r")$Sample.Name
reps.final <- data.frame(orig, reps)
max.key <- apply(reps.final, 1, which.max)
reps.names <- data.frame(orig.names, reps.names)

remove.names <- reps.names[cbind(seq_along(max.key), max.key)]
remove.key <- paste(remove.names, collapse = "|")
reps.samples <- grepl(remove.key, sampleNames(lumi.patient))
remove.all <- remove.indices | reps.samples
saveRDS.gz(remove.all, "./save/remove.all.rda")

lumi.rmreps <- lumi.patient[,!remove.all]
saveRDS.gz(lumi.rmreps, "./save/lumi.rmreps.rda")

PIDN.long <- filter(pData(lumi.rmreps), Sample.Num == "4")$PIDN 
lumi.rmreps$Sample.Num[lumi.rmreps$Sample.Num == "1r"] <- 1
lumi.rmreps$Sample.Num %<>% as.character %>% as.numeric

targets.rmreps <- filter(pData(lumi.rmreps), PIDN %in% PIDN.long) 
targets.nums <- by(targets.rmreps, targets.rmreps$PIDN, nrow) %>% as.list %>% melt %>% data.frame 
missing.PIDNs <- filter(targets.nums, value < 4)$L1 
targets.final <- filter(targets.rmreps, !grepl(match.exact(missing.PIDNs), PIDN))
rmreps.index <- sampleNames(lumi.rmreps) %in% targets.final$Sample.Name 
lumi.rmreps <- lumi.rmreps[,rmreps.index]
onebatch.pidn <- filter(pData(lumi.rmreps), Batch == 5)$PIDN
lumi.rmreps <- lumi.rmreps[,lumi.rmreps$PIDN != onebatch.pidn]

lumi.rmreps.norm <- lumiN(lumi.rmreps, method = "rsn") #Normalize with robust spline regression
lumi.rmreps.qual <- lumiQ(lumi.rmreps.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.rmreps.cutoff <- detectionCall(lumi.rmreps.qual) #Get the count of probes which passed the detection threshold per sample
lumi.rmreps.expr <- lumi.rmreps.qual[which(lumi.rmreps.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.rmreps <- getSYMBOL(rownames(lumi.rmreps.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.rmreps.annot <- lumi.rmreps.expr[!symbols.lumi.rmreps,] #Drop any probe which is not annotated
saveRDS.gz(lumi.rmreps.annot, file = "./save/lumi.rmreps.annot.rda")

model.sex <- model.matrix( ~ 0 + factor(lumi.rmreps.annot$Sex) )
colnames(model.sex) <- c("Female", "Male")
model.sex.reduce <- model.sex[,-1]

model.combat <- cbind(Male = model.sex.reduce, Age = as.numeric(lumi.rmreps.annot$Draw.Age), RIN = lumi.rmreps.annot$RIN)

expr.combat <- ComBat(dat = exprs(lumi.rmreps.annot), batch = factor(lumi.rmreps.annot$Batch), mod = model.combat)
lumi.combat <- lumi.rmreps.annot
exprs(lumi.combat) <- expr.combat
saveRDS.gz(lumi.combat, "./save/lumi.combat")

model.cov <- cbind(Male = model.sex.reduce, Age = as.numeric(lumi.combat$Draw.Age), RIN = lumi.combat$RIN)
#model.full.cov <- cbind(model.cov, model.PEER_covariate)
export.expr <- removeBatchEffect(exprs(lumi.combat), covariates = model.cov)
export.lumi <- lumi.combat
exprs(export.lumi) <- export.expr
saveRDS.gz(export.lumi, file = "./save/export.lumi.rda")

lumi.collapse.expr <- collapseRows(exprs(export.lumi), getSYMBOL(featureNames(export.lumi), 'lumiHumanAll.db'), rownames(export.lumi))
lumi.collapse <- ExpressionSet(assayData = lumi.collapse.expr$datETcollapsed, phenoData = phenoData(export.lumi))
saveRDS.gz(lumi.collapse, './save/lumi.collapse.rda')

lumi.ubc <- lumi.collapse["UBC",]
lumi.ubc.df <- cbind(t(exprs(lumi.ubc[,lumi.ubc$Sample.Num == 1])), t(exprs(lumi.ubc[,lumi.ubc$Sample.Num == 2])), t(exprs(lumi.ubc[,lumi.ubc$Sample.Num == 3])), t(exprs(lumi.ubc[,lumi.ubc$Sample.Num == 4])))

lumi.ubb <- lumi.collapse["UBB",]
lumi.ubb.df <- cbind(t(exprs(lumi.ubb[,lumi.ubb$Sample.Num == 1])), t(exprs(lumi.ubb[,lumi.ubb$Sample.Num == 2])), t(exprs(lumi.ubb[,lumi.ubb$Sample.Num == 3])), t(exprs(lumi.ubb[,lumi.ubb$Sample.Num == 4])))

expr.means <- collapseRows(t(exprs(lumi.collapse)), lumi.collapse$Sample.Num, lumi.collapse$Sample.Name, method = "function", methodFunction = colMeans)$datETcollapsed %>% t
saveRDS.gz(expr.means, "./save/expr.means.rda")

gen.permuts <- function(expr.means, expr.list = list(), count.value = 1) {
    print(str_c("permutation ", count.value))
    if (count.value < 1000) {
        permut.means <- split(expr.means, row(expr.means)) %>% map(sample) %>% reduce(rbind)
        rownames(permut.means) <- rownames(expr.means)
        expr.list[[count.value]] <- permut.means
        count.value %<>% add(1)
        gen.permuts(expr.means, expr.list, count.value)
    }
    else {
        permut.means <- split(expr.means, row(expr.means)) %>% map(sample) %>% reduce(rbind)
        rownames(permut.means) <- rownames(expr.means)
        expr.list[[count.value]] <- permut.means
        count.value %<>% add(1)
        return(expr.list)
    }
}

get.pvalue <- function(test.value, dist.vector)
{
    num.above <- which(dist.vector > test.value) %>% length
    return(num.above / length(dist.vector))
}

ica.all <- seed.ICA(list(Patient = expr.means))
saveRDS.gz(ica.all, "./save/ica.all.rda")
#ica.patient <- ica.all[str_detect(names(ica.all), "Patient")]
#ica.melt <- melt(ica.patient)
#ica.snca <- filter(ica.melt, Symbol == "MMP9")
#snca.means <- by(ica.snca, factor(ica.snca$variable), select, value) %>% map(Compose(abs, median))

set.seed(12345)
ica.test1 <- fastICA(expr.means, 4, "deflation") 

permut.means <- gen.permuts(expr.means)
saveRDS.gz(permut.means, "./save/permut.means.rda")
permut.ica <- map(permut.means, fastICA, 4, "deflation")
saveRDS.gz(permut.ica, "./save/permut.ica.rda")
ica.score <- list.map(permut.ica, as.vector(S)) %>% do.call(what = rbind) %>% t
orig.score.sig  <- map2_dbl(as.vector(ica.test1$S), split(ica.score, row(ica.score)), get.pvalue) 
dim(orig.score.sig) <- c(nrow(expr.means), 4)
colnames(orig.score.sig) <- c(str_c("X", 1:4))
rownames(orig.score.sig) <- rownames(expr.means)
orig.score.df <- data.frame(orig.score.sig)
orig.score.df$Symbol <- rownames(orig.score.df)

ica.final.df <- data.frame(ica.test$S)
ica.final.df$Symbol <- rownames(expr.means)
ica.abs <- select(ica.final.df, -Symbol) %>% abs
ica.filter <- apply((ica.abs > 3), 1, any)
ica.final.filter <- ica.final.df[ica.filter,]
orig.score.filter <- orig.score.df[ica.filter,]

ica.assigned <- apply(select(ica.final.filter, -Symbol), 1, which.max)
ica.final.filter$Component <- ica.assigned

ica.X2 <- filter(ica.final.filter, Component == 2)
ica.X2.pval <- orig.score.filter[ica.X2$Symbol,]
ica.X2.final <- filter(ica.X2.pval, X2 < 0.005)
write.xlsx(ica.X2.final, "ica.X2.xlsx")

ica.X3 <- filter(ica.final.filter, Component == 3)
ica.X3.pval <- orig.score.filter[ica.X3$Symbol,]
ica.X3.final <- filter(ica.X3.pval, X3 < 0.005)
write.xlsx(ica.X3.final, "ica.X3.xlsx")

ica.X4 <- filter(ica.final.filter, Component == 4)
ica.X4.pval <- orig.score.filter[ica.X4$Symbol,]
ica.X4.final <- filter(ica.X4.pval, X4 < 0.005)
write.xlsx(ica.X4.final, "ica.X4.xlsx")

ica.genes <- c(ica.X2.final$Symbol, ica.X3.final$Symbol, ica.X4.final$Symbol)
saveRDS.gz(ica.genes, "./save/ica.genes.rda")

#names(ica.all) <- "Patient"
ica.split <- map(unique(names(ica.all)), collapse.ica, ica.all)
names(ica.split) <- unique(names(ica.all))

#intensities.ica <- map(intensities.means, select, -Symbol) %>% map(fastICA, 4) %>% map(`[[`, "S") %>% map(data.frame) 
intensities.ica <- map(ica.split, gen.ica)

intensities.ica.genes <- map(intensities.ica, apply, 1, which.max) %>% map(split.ica)
test.name <- intensities.ica.genes$Patient %>% map(select, Symbol) %>% map(Compose(unlist, as.character)) %>% reduce(c) 

l_ply(names(intensities.ica.genes), gen.icatables, intensities.ica.genes, "ica.mean")

source("../../code/GO/enrichr.R")
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "Humancyc_2016", "NCI-Nature_2016", "Panther_2016") 
map(names(intensities.ica.genes), enrichr.ica, intensities.ica.genes, enrichr.terms, "intensities.mean")
ica.ppinets <- map(names(intensities.ica.genes), stringdb.ica, intensities.ica.genes, "means")
ica.ppinets <- ica.ppinets[[1]]
names(ica.ppinets) <- names(intensities.ica.genes[[1]])

symbols.subnet <- ica.ppinets$X3
edge.threshold <- 400
edge.weights <- edge_attr(symbols.subnet, "combined_score")
pruned.subnet <- delete.edges(symbols.subnet, which(edge.weights < edge.threshold))
num.edges <- map(1:vcount(pruned.subnet), incident, graph = pruned.subnet) %>% map_dbl(length) 
pruned.subnet2 <- delete.vertices(pruned.subnet, which(num.edges == 0))
communities.optimal <- cluster_optimal(pruned.subnet2, weights = edge_attr(pruned.subnet2, "combined_score")/1000)
vertex.colors <- rainbow(vcount(pruned.subnet2))
V(pruned.subnet2)$color <- vertex.colors
edge.df <- data.frame(edge_attr(symbols.subnet))
edge.thickness <- edge.df$combined_score / 200

filepath <- file.path(prefix, plot.name)

CairoPDF(filepath, width = 30, height = 30)
plot.igraph(pruned.subnet2, vertex.size = 2, vertex.label.dist = 0.12, vertex.label.degree = -pi/2, vertex.label.font = 2, vertex.label.color = "black", edge.width = edge.thickness, edge.color = "#0000FF99")
dev.off()
return(symbols.subnet)

saveRDS.gz(intensities.ica.genes, "./save/intensities.ica.genes.rda")

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

#X1 GO
X1.gobiol.file <- "./enrichr/intensities.mean/Patient/X1/GO_Biological_Process_2015.xlsx"
X1.gobiol <- read.xlsx(X1.gobiol.file)
X1.gobiol$Database <- "GO Biological Process"
X1.gobiol$Num.Genes <- map(X1.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X1.gobiol.filter <- filter(X1.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(X1.gobiol.filter, intensities.ica.genes$Patient$X1$Symbol, file_path_sans_ext(X1.gobiol.file))
X1.gobiol.final <- slice(X1.gobiol, c(13, 18, 62, 16))

X1.gomolec.file <- "./enrichr/intensities.mean/Patient/X1/GO_Molecular_Function_2015.xlsx"
X1.gomolec <- read.xlsx(X1.gomolec.file)
X1.gomolec$Database <- "GO Molecular Function"
X1.gomolec$Num.Genes <- map(X1.gomolec$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X1.gomolec.filter <- filter(X1.gomolec, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(X1.gomolec.filter, intensities.ica.genes$Patient$X1$Symbol, file_path_sans_ext(X1.gomolec.file))
X1.gomolec.final <- slice(X1.gomolec, 6)

X1.reactome.file <- "./enrichr/intensities.mean/Patient/X1/Reactome_2016.xlsx"
X1.reactome <- read.xlsx(X1.reactome.file)
X1.reactome$Database <- "Reactome"
X1.reactome$Num.Genes <- map(X1.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X1.reactome.filter <- filter(X1.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(X1.reactome.filter, intensities.ica.genes$Patient$X1$Symbol, file_path_sans_ext(X1.reactome.file))
X1.reactome.final <- slice(X1.reactome, 1)

X1.kegg <- read.xlsx("./enrichr/intensities.mean/Patient/X1/KEGG_2016.xlsx") %>% slice(5)
X1.kegg$Num.Genes <- map(X1.kegg$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X1.kegg$Database <- "KEGG"

X1.enrichr <- rbind(X1.gobiol.final, X1.kegg, X1.reactome.final)
gen.enrichrplot(X1.enrichr, "X1.enrichr")

#X2 - probably skip
X2.gobiol.file <- "./enrichr/intensities.mean/Patient/X2/GO_Biological_Process_2015.xlsx"
X2.gobiol <- read.xlsx(X2.gobiol.file)
X2.gobiol$Database <- "GO Biological Process"
X2.gobiol$Num.Genes <- map(X2.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X2.gobiol.filter <- filter(X2.gobiol, Num.Genes > 4) %>% filter(P.value < 0.01)
get.kappa.cluster(X2.gobiol.filter, intensities.ica.genes$Patient$X2$Symbol, file_path_sans_ext(X2.gobiol.file))

#X3 GO

X3.gobiol.file <- "./enrichr/intensities.mean/Patient/X3/GO_Biological_Process_2015.xlsx"
X3.gobiol <- read.xlsx(X3.gobiol.file)
X3.gobiol$Database <- "GO Biological Process"
X3.gobiol$Num.Genes <- map(X3.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X3.gobiol.filter <- filter(X3.gobiol, Num.Genes > 4) %>% filter(P.value < 0.001)
get.kappa.cluster(X3.gobiol.filter, intensities.ica.genes$Patient$X3$Symbol, file_path_sans_ext(X3.gobiol.file))
X3.gobiol.final <- slice(X3.gobiol, c(39, 49, 61))

X3.gomolec.file <- "./enrichr/intensities.mean/Patient/X3/GO_Molecular_Function_2015.xlsx"
X3.gomolec <- read.xlsx(X3.gomolec.file)
X3.gomolec$Database <- "GO Molecular Function"
X3.gomolec$Num.Genes <- map(X3.gomolec$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X3.gomolec.filter <- filter(X3.gomolec, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(X3.gomolec.filter, intensities.ica.genes$Patient$X3$Symbol, file_path_sans_ext(X3.gomolec.file))
X3.gomolec.final <- slice(X3.gomolec, c(8,19))

X3.reactome.file <- "./enrichr/intensities.mean/Patient/X3/Reactome_2016.xlsx"
X3.reactome <- read.xlsx(X3.reactome.file)
X3.reactome$Database <- "Reactome"
X3.reactome$Num.Genes <- map(X3.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X3.reactome.filter <- filter(X3.reactome, Num.Genes > 4) %>% filter(P.value < 0.01)
get.kappa.cluster(X3.reactome.filter, intensities.ica.genes$Patient$X3$Symbol, file_path_sans_ext(X3.reactome.file))
X3.reactome.final <- slice(X3.reactome, c(36, 38))
X3.enrichr <- rbind(X3.gobiol.final, X3.gomolec.final, X3.reactome.final)
gen.enrichrplot(X3.enrichr, "X3.enrichr")

#X4
X4.gobiol.file <- "./enrichr/intensities.mean/Patient/X4/GO_Biological_Process_2015.xlsx"
X4.gobiol <- read.xlsx(X4.gobiol.file)
X4.gobiol$Database <- "GO Biological Process"
X4.gobiol$Num.Genes <- map(X4.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X4.gobiol.filter <- filter(X4.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(X4.gobiol.filter, intensities.ica.genes$Patient$X4$Symbol, file_path_sans_ext(X4.gobiol.file))
X4.gobiol.final <- slice(X4.gobiol, c(1, 14))

X4.reactome.file <- "./enrichr/intensities.mean/Patient/X4/Reactome_2016.xlsx"
X4.reactome <- read.xlsx(X4.reactome.file)
X4.reactome$Database <- "Reactome"
X4.reactome$Num.Genes <- map(X4.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
X4.reactome.filter <- filter(X4.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(X4.reactome.filter, intensities.ica.genes$Patient$X4$Symbol, file_path_sans_ext(X4.reactome.file))
X4.reactome.final <- slice(X4.reactome, c(2,3))

X4.enrichr <- rbind(X4.gobiol.final, X4.reactome.final)
gen.enrichrplot(X4.enrichr, "X4.enrichr")
