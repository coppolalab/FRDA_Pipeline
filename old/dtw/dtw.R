#String operations
library(stringr)
library(tools)

#Reading and writing tables
library(readr)
library(openxlsx)

#For plotting
library(ggplot2)

#For microarray stuff
library(Biobase)
library(matrixStats)
library(abind)
library(lumi)
library(lumiHumanAll.db)
library(annotate)
library(sva)
library(limma)

#Longitudinal analysis
library(dtw)
library(irr)

#Plotting
library(Cairo)
library(WGCNA)
library(heatmap.plus)
library(flashClust)
enableWGCNAThreads()

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(vadr)

#Data arrangement
library(data.table)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(rlist)
library(parallel)

match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|"))

gen.median <- function(dataset)
{
    intensities.median <- by(dataset, factor(dataset$variable), select, -PIDN, -variable, -Status) %>% lapply(Compose(as.matrix, colMedians)) %>% reduce(rbind) %>% data.frame
    intensities.median$Symbol <- unique(dataset$variable)
    return(intensities.median)
}

gen.mean <- function(expr.orig, pdata)
{
    id.filter <- rownames(expr.orig) %>% match.exact
    timepoints.vector <- filter(pdata, grepl(id.filter, Sample.Name))$Sample.Num
    expr.means <- collapseRows(expr.orig, timepoints.vector, rownames(expr.orig), method = "function", methodFunction = colMeans)$datETcollapsed
    return(t(expr.means))
}

gen.dist <- function(dataset)
{
    dtws <- map2(split(dataset[[1]], row(dataset[[1]])), split(dataset[[2]], row(dataset[[2]])), dtw)
    names(dtws) <- rownames(dataset[[1]])
    dists <- list.map(dtws, normalizedDistance)
    return(dists)
}

gen.dtw <- function(dataset)
{
    combinations <- combn(dataset, 2, simplify = FALSE)
    dists <- map(combinations, gen.dist)
    return(dists)
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

gen.boxplot <- function(filename, lumi.object, colorscheme, maintext, ylabtext)
{
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

gen.small.workbook <- function(dataset, filename)
{
    #coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    #colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")
    #dataset$DEFINITION %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    #conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:3, widths = "auto")
    #setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 4, widths = 45)
    setColWidths(wb, 1, cols = 5:7, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
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

    dir.create(new.fullpath, showWarnings = FALSE, recursive = TRUE)
    filename = paste(file.path(new.fullpath, subindex), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
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

gen.permuts <- function(status.dt, pdata.reduce, status.list = list(), count.value = 1) {
    print(str_c("permutation ", count.value))
    if (count.value < 1000) {
        copy.dt <- status.dt
        copy.dt$Status %<>% sample
        pdata.permuted <- join(pdata.reduce, copy.dt)
        status.list[[count.value]] <- pdata.permuted
        count.value %<>% add(1)
        gen.permuts(status.dt, pdata.reduce, status.list, count.value)
    }
    else {
        copy.dt <- status.dt
        copy.dt$Status %<>% sample
        pdata.permuted <- join(pdata.reduce, copy.dt)
        status.list[[count.value]] <- pdata.permuted
        return(status.list)
    }
}

lumi.known <- readRDS.gz("../baseline_lumi/save/lumi.known.rda")
lumi.vst <- lumiT(lumi.known)

long.key <- grepl("^1$|1r|2|3|4", lumi.vst$Sample.Num)
lumi.long <- lumi.vst[,long.key]
saveRDS.gz(lumi.long, file = "./save/lumi.long.rda")

PIDN.long <- filter(pData(lumi.long), Sample.Num == "4")$PIDN 

targets.long <- filter(pData(lumi.long), PIDN %in% PIDN.long) 
targets.nums <- by(targets.long, targets.long$PIDN, nrow) %>% as.list %>% melt %>% data.frame 

lumi.four <- lumi.long[,lumi.long$PIDN %in% PIDN.long]

lumi.norm <- lumiN(lumi.four, method = "rsn")
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
connectivity.remove <- connectivity.zscore[abs(connectivity.zscore) > 3] %>% names

combined.remove <- c(qcsum.remove, connectivity.remove) %>% unique %>% paste(collapse = "|")
remove.indices <- grepl(combined.remove, sampleNames(lumi.expr.annot))

lumi.rmout <- lumi.four[,!remove.indices]
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
reps.samples <- grepl(remove.key, sampleNames(lumi.four))
remove.all <- remove.indices | reps.samples
saveRDS.gz(remove.all, "./save/remove.all.rda")

lumi.rmreps <- lumi.four[,!remove.all]

PIDN.long <- filter(pData(lumi.rmreps), Sample.Num == "4")$PIDN 
lumi.rmreps$Sample.Num[lumi.rmreps$Sample.Num == "1r"] <- 1
lumi.rmreps$Sample.Num %<>% as.character %>% as.numeric

targets.rmreps <- filter(pData(lumi.rmreps), PIDN %in% PIDN.long) 
targets.nums <- by(targets.rmreps, targets.rmreps$PIDN, nrow) %>% as.list %>% melt %>% data.frame 
missing.PIDNs <- filter(targets.nums, value < 4)$L1 
targets.final <- filter(targets.rmreps, !grepl(match.exact(missing.PIDNs), PIDN))
rmreps.index <- sampleNames(lumi.rmreps) %in% targets.final$Sample.Name 
onebatch.pidn <- filter(pData(lumi.rmreps), Batch == 5)$PIDN
lumi.rmreps <- lumi.rmreps[,rmreps.index]
lumi.rmreps <- lumi.rmreps[,lumi.rmreps$PIDN != onebatch.pidn]

lumi.rmreps.norm <- lumiN(lumi.rmreps, method = "rsn") #Normalize with robust spline regression
lumi.rmreps.qual <- lumiQ(lumi.rmreps.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.rmreps.cutoff <- detectionCall(lumi.rmreps.qual) #Get the count of probes which passed the detection threshold per sample
lumi.rmreps.expr <- lumi.rmreps.qual[which(lumi.rmreps.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.rmreps <- getSYMBOL(rownames(lumi.rmreps.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.rmreps.annot <- lumi.rmreps.expr[!symbols.lumi.rmreps,] #Drop any probe which is not annotated
saveRDS.gz(lumi.rmreps.annot, file = "./save/lumi.rmreps.annot.rda")

model.status <- model.matrix( ~ 0 + factor(lumi.rmreps.annot$Status) )
colnames(model.status) <- c("Carrier", "Control", "Patient")
model.status.reduce <- model.status[,-2]

model.sex <- model.matrix( ~ 0 + factor(lumi.rmreps.annot$Sex) )
colnames(model.sex) <- c("Female", "Male")
model.sex.reduce <- model.sex[,-1]

model.combat <- cbind(model.status.reduce, Male = model.sex.reduce, Age = as.numeric(lumi.rmreps.annot$Draw.Age), RIN = lumi.rmreps.annot$RIN)

expr.combat <- ComBat(dat = exprs(lumi.rmreps.annot), batch = factor(lumi.rmreps.annot$Batch), mod = model.combat)
lumi.combat <- lumi.rmreps.annot
exprs(lumi.combat) <- expr.combat
saveRDS.gz(lumi.combat, "./save/lumi.combat")

source("../common_functions.R")


#Removing effects of covariates + PEER factors  !!DO NOT USE FOR LINEAR MODELING WITH CONTRASTS!!
model.cov <- cbind(Male = model.sex.reduce, Age = as.numeric(lumi.combat$Draw.Age), RIN = lumi.combat$RIN)
#model.full.cov <- cbind(model.cov, model.PEER_covariate)
export.expr <- removeBatchEffect(exprs(lumi.combat), covariates = model.cov)
export.lumi <- lumi.combat
exprs(export.lumi) <- export.expr
saveRDS.gz(export.lumi, file = "./save/export.lumi.rda")

batch.colors <- c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green", "blue4", "red4")
gen.boxplot("baseline_intensity_corrected.jpg", export.lumi, batch.colors, "Covariate-corrected intensity", "Intensity")

lumi.collapse.expr <- collapseRows(exprs(lumi.final), getSYMBOL(featureNames(lumi.final), 'lumiHumanAll.db'), rownames(lumi.final), method = "function", methodFunction = colMeans)
lumi.collapse <- ExpressionSet(assayData = lumi.collapse.expr$datETcollapsed, phenoData = phenoData(lumi.final))

saveRDS.gz(lumi.final, "./save/lumi.final.rda")
saveRDS.gz(lumi.collapse, "./save/lumi.collapse.rda")

#intensities.long <- data.frame(Sample.Num = lumi.final$Sample.Num, PIDN = lumi.final$PIDN, Status = lumi.final$Status, lumi.exprs.collapse) 
#intensities.melt <- melt(intensities.long, id.vars = c("Sample.Num", "PIDN", "Status"))
#intensities.summed <- dcast(intensities.melt, PIDN + variable + Status ~ Sample.Num, value.var = "value")
#saveRDS.gz(intensities.summed, "./save/intensities.summed.rda")

source("../../code/GO/enrichr.R")

intensities.medians <- by(intensities.summed, factor(intensities.summed$Status), gen.median)
saveRDS.gz(intensities.medians, "./save/intensities.medians.rda")
intensities.dists <- gen.dtw(intensities.medians) #%>% reduce(cbind) %>% data.frame
coltitles <- c("Carrier.vs.Control", "Patient.vs.Carrier", "Patient.vs.Control")
rownames(intensities.dists) <- colnames(lumi.exprs.collapse)
colnames(intensities.dists) <- coltitles
intensities.dists$Symbol <- rownames(intensities.dists)

dtw.cc <- arrange(intensities.dists, desc(Carrier.vs.Control))
dtw.pca <- select(intensities.dists, Symbol, Patient.vs.Carrier, Patient.vs.Control, Carrier.vs.Control) %>% arrange(desc(Patient.vs.Carrier))
dtw.pco <- select(intensities.dists, Symbol, Patient.vs.Control, Patient.vs.Carrier, Carrier.vs.Control) %>% arrange(desc(Patient.vs.Control))
gen.small.workbook(dtw.cc, "carrier.control.median.xlsx")
gen.small.workbook(dtw.pca, "patient.carrier.median.xlsx")
gen.small.workbook(dtw.pco, "patient.control.median.xlsx")

dtw.cc.submit <- dtw.cc[1:300,]
dtw.pca.submit <- dtw.pca[1:300,]
dtw.pco.submit <- dtw.pco[1:300,]

get.stringdb(dtw.cc.submit, "cc.median", "cc")
get.stringdb(dtw.pca.submit, "pca.median", "pca")
get.stringdb(dtw.pco.submit, "pco.median", "pco")

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "Humancyc_2016", "NCI-Nature_2016", "Panther_2016") 

dtw.cc.enrichr <- map(enrichr.terms, get.enrichrdata, dtw.cc.submit, FALSE)
names(dtw.cc.enrichr) <- enrichr.terms
map(names(dtw.cc.enrichr), enrichr.wkbk, dtw.cc.enrichr, "./enrichr/cc.median")

dtw.pca.enrichr <- map(enrichr.terms, get.enrichrdata, dtw.pca.submit, FALSE)
names(dtw.pca.enrichr) <- enrichr.terms
map(names(dtw.pca.enrichr), enrichr.wkbk, dtw.pca.enrichr, "./enrichr/pca.median")

dtw.pco.enrichr <- map(enrichr.terms, get.enrichrdata, dtw.pco.submit, FALSE)
names(dtw.pco.enrichr) <- enrichr.terms
map(names(dtw.pco.enrichr), enrichr.wkbk, dtw.pco.enrichr, "./enrichr/pco.median")

pdata.reduce <- select(pData(lumi.collapse), Sample.Name, Sample.Num, Status, PIDN) %>% data.table
pdata.pidn <- select(pdata.reduce, PIDN, Status) %>% unique
pidn.permuts <- gen.permuts(pdata.pidn, select(pdata.reduce, -Status))

gen.mean.wrapper <- function(pdata.dt, lumi.collapse) {
    intensities.means <- by(t(exprs(lumi.collapse)), factor(pdata.dt$Status), gen.mean, pdata.dt) 
    return(intensities.means) 
}

intensities.means <- by(t(exprs(lumi.collapse)), factor(lumi.collapse$Status), gen.mean, pData(lumi.collapse))
means.permut <- map(pidn.permuts, gen.mean.wrapper, lumi.collapse)
saveRDS.gz(means.permut, "./save/means.permut.rda")
saveRDS.gz(intensities.means, "./save/intensities.means.rda")
intensities.dists.m <- gen.dtw(intensities.means) %>% reduce(cbind) %>% data.frame

#test <- gen.dist(list(intensities.means$Carrier, intensities.means$Control))
#rownames(intensities.dists.m) <- colnames(lumi.exprs.collapse)
colnames(intensities.dists.m) <- coltitles
intensities.dists.m$Symbol <- rownames(intensities.dists.m)

dtw.cc.m <- arrange(intensities.dists.m, desc(Carrier.vs.Control))
dtw.pca.m <- select(intensities.dists.m, Symbol, Patient.vs.Carrier, Patient.vs.Control, Carrier.vs.Control) %>% arrange(desc(Patient.vs.Carrier))
dtw.pco.m <- select(intensities.dists.m, Symbol, Patient.vs.Control, Patient.vs.Carrier, Carrier.vs.Control) %>% arrange(desc(Patient.vs.Control))

system.time(dist.permuts <- mclapply(means.permut[1:10], gen.dtw, mc.cores = 8))

system.time(dist.permuts <- map(means.permut[1:10], gen.dtw, mc.cores = 8))
saveRDS.gz(dist.permuts, "./save/dist.permuts.rda")

dist.permuts <- readRDS.gz("./dist.permuts.rda")
dist.permuts.pco <- mclapply(dist.permuts, set_names, coltitles, mc.cores = 4) %>% mclapply(extract2, "Patient.vs.Control", mc.cores = 4) %>% mclapply(unlist, mc.cores = 4) %>% reduce(rbind) %>% t
dist.permuts.pca <- mclapply(dist.permuts, set_names, coltitles, mc.cores = 4) %>% mclapply(extract2, "Patient.vs.Carrier", mc.cores = 4) %>% mclapply(unlist, mc.cores = 4) %>% reduce(rbind) %>% t
dist.permuts.cc <- mclapply(dist.permuts, set_names, coltitles, mc.cores = 4) %>% mclapply(extract2, "Carrier.vs.Control", mc.cores = 4) %>% mclapply(unlist, mc.cores = 4) %>% reduce(rbind) %>% t

saveRDS.gz(dist.permuts.pco, "./save/dist.permuts.pco.rda")
saveRDS.gz(dist.permuts.pca, "./save/dist.permuts.pca.rda")
saveRDS.gz(dist.permuts.cc, "./save/dist.permuts.cc.rda")

dtw.sorted <- arrange(dtw.pca.m, Symbol)
pvalue.pca <- map2_dbl(dtw.sorted$Patient.vs.Carrier, split(dist.permuts.pca, 1:nrow(dist.permuts.pca)), get.pvalue)
pvalue.pco <- map2_dbl(dtw.sorted$Patient.vs.Control, split(dist.permuts.pco, 1:nrow(dist.permuts.pco)), get.pvalue)
pvalue.cc <- map2_dbl(dtw.sorted$Carrier.vs.Control, split(dist.permuts.cc, 1:nrow(dist.permuts.cc)), get.pvalue)

get.pvalue <- function(test.value, dist.vector)
{
    num.above <- which(dist.vector > test.value) %>% length
    return(num.above / length(dist.vector))
}

permut.pca <- readRDS.gz("./save/permut.df.pca.rda")
permut.pco <- readRDS.gz("./save/permut.df.pco.rda")
permut.cc <- readRDS.gz("./save/permut.df.cc.rda")

sig.pca <- slice(dtw.sorted, which(pvalue.pca < 0.05))
sig.pca$pvalue <- pvalue.pca[which(pvalue.pca < 0.05)]
sig.pca %<>% arrange(pvalue, desc(Patient.vs.Carrier))
write.xlsx(sig.pca, "./sig.pca.xlsx")

sig.pco <- slice(dtw.sorted, which(pvalue.pco < 0.05))
sig.pco$pvalue <- pvalue.pco[which(pvalue.pco < 0.05)]
sig.pco %<>% arrange(pvalue, desc(Patient.vs.Control))
write.xlsx(sig.pco, "./sig.pco.xlsx")

sig.cc <- slice(dtw.sorted, which(pvalue.cc < 0.05))
sig.cc$pvalue <- pvalue.cc[which(pvalue.cc < 0.05)]
sig.cc %<>% arrange(pvalue, desc(Carrier.vs.Control))
write.xlsx(sig.cc, "./sig.cc.xlsx")

saveRDS.gz(dtw.cc.m, "./save/dtw.cc.rda")
saveRDS.gz(dtw.pca.m, "./save/dtw.pca.rda")
saveRDS.gz(dtw.pco.m, "./save/dtw.pco.rda")
gen.small.workbook(dtw.cc.m, "carrier.control.mean.xlsx")
gen.small.workbook(dtw.pca.m, "patient.carrier.mean.xlsx")
gen.small.workbook(dtw.pco.m, "patient.control.mean.xlsx")

dtw.cc.m.submit <- dtw.cc.m[1:300,]
dtw.pca.m.submit <- dtw.pca.m[1:300,]
dtw.pco.m.submit <- dtw.pco.m[1:300,]
saveRDS.gz(dtw.cc.m.submit, "./save/dtw.cc.submit.rda")
saveRDS.gz(dtw.pca.m.submit, "./save/dtw.pca.submit.rda")
saveRDS.gz(dtw.pco.m.submit, "./save/dtw.pco.submit.rda")

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

symbols.only <- intensities.means$Patient$Symbol

#gen.permuts <- function(symbols.vector, symbols.list = list(), count.value = 1)
#{
   #if (count.value < 1000)
   #{
        #symbols.list[[count.value]] <- sample(symbols.vector)
        #count.value = count.value + 1
        #gen.permuts(symbols.vector, symbols.list, count.value)
   #}
   #else
   #{
        #symbols.list[[count.value]] <- sample(symbols.vector)
        #return(symbols.list)
   #}
#}

#permut.dtw <- function(symbols.new, intensities.first, intensities.second)
#{
    #intensities.first$Symbol <- symbols.new
    #intensities.second$Symbol <- symbols.new
    #intensities.first %<>% arrange(Symbol) %>% select(-Symbol)
    #intensities.second %<>% arrange(Symbol) %>% select(-Symbol)
    #dtw.permut <- map2(split(intensities.first, 1:nrow(intensities.first)), split(intensities.second, 1:nrow(intensities.second)), dtw) %>% map_dbl(`[[`, "normalizedDistance")
    #return(dtw.permut)
#}

#symbols.permut <- gen.permuts(symbols.only)
#ptm <- proc.time()
#test.permut <- map(symbols.permut[1:20], permut.dtw, intensities.means$Patient, intensities.means$Carrier)
#final.time <- proc.time() - ptm
#permut.df <- reduce(test.permut, cbind) %>% data.frame
#rownames(permut.df) <- symbols.only

#Permutations run on Hoffman because they are too slow to calculate otherwise


get.stringdb(sig.cc, "cc.mean", "cc", 900)
get.stringdb(sig.pca, "pca.mean", "pca", 900)
get.stringdb(sig.pco, "pco.mean", "pco", 900)

dtw.cc.m.enrichr <- map(enrichr.terms, get.enrichrdata, dtw.cc.m.submit, FALSE)
names(dtw.cc.m.enrichr) <- enrichr.terms
map(names(dtw.cc.m.enrichr), enrichr.wkbk, dtw.cc.m.enrichr, "./enrichr/cc.mean")

dtw.pca.m.enrichr <- map(enrichr.terms, get.enrichrdata, dtw.pca.m.submit, FALSE)
names(dtw.pca.m.enrichr) <- enrichr.terms
map(names(dtw.pca.m.enrichr), enrichr.wkbk, dtw.pca.m.enrichr, "./enrichr/pca.mean")

dtw.pco.m.enrichr <- map(enrichr.terms, get.enrichrdata, dtw.pco.m.submit, FALSE)
names(dtw.pco.m.enrichr) <- enrichr.terms
map(names(dtw.pco.m.enrichr), enrichr.wkbk, dtw.pco.m.enrichr, "./enrichr/pco.mean")

gen.enrichrplot <- function(enrichr.df, filename, plot.height = 5, plot.width = 8)
{
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
    enrichr.df$Log.pvalue <- -(log10(enrichr.df$P.value))
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*\\)", "") %>% str_replace_all("\\_.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
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

pca.gobiol.file <- "./enrichr/pca.mean/GO_Biological_Process_2015.xlsx" 
pca.gobiol <- read.xlsx(pca.gobiol.file)
pca.gobiol$Num.Genes <- map(pca.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pca.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
pca.gobiol$Database <- "GO Biological Process"
get.kappa.cluster(pca.gobiol, intensities.dists$Symbol, file_path_sans_ext(pca.gobiol.file))
pca.gobiol.final <- slice(pca.gobiol, c(1, 22, 73, 82, 35))

pca.gomole.file <- "./enrichr/pca.mean/GO_Molecular_Function_2015.xlsx" 
pca.gomole <- read.xlsx(pca.gomole.file)
pca.gomole$Num.Genes <- map(pca.gomole$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pca.gomole$Database <- "GO Molecular Function"
pca.gomole.final <- slice(pca.gomole, c(1, 13))

pca.reactome.file <- "./enrichr/pca.mean/Reactome_2016.xlsx" 
pca.reactome <- read.xlsx(pca.reactome.file)
pca.reactome$Num.Genes <- map(pca.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pca.reactome %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
pca.reactome$Database <- "Reactome"
get.kappa.cluster(pca.reactome, intensities.dists$Symbol, file_path_sans_ext(pca.reactome.file))
pca.reactome.final <- slice(pca.reactome, c(24, 48))

pca.kegg.file <- "./enrichr/pca.mean/KEGG_2016.xlsx" 
pca.kegg <- read.xlsx(pca.kegg.file)
pca.kegg$Num.Genes <- map(pca.kegg$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pca.kegg %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
pca.kegg$Database <- "KEGG"
pca.kegg.final <- slice(pca.kegg, 19)

pca.enrichr <- rbind(pca.gobiol.final, pca.gomole.final, pca.reactome.final)
gen.enrichrplot(pca.enrichr, "pca.enrichr")

pco.gobiol.file <- "./enrichr/pco.mean/GO_Biological_Process_2015.xlsx" 
pco.gobiol <- read.xlsx(pco.gobiol.file)
pco.gobiol$Database <- "GO Biological Process"
pco.gobiol$Num.Genes <- map(pco.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pco.gobiol.filter <- filter(pco.gobiol, Num.Genes > 4) %>% filter(P.value < 0.001)
get.kappa.cluster(pco.gobiol.filter, intensities.dists$Symbol, file_path_sans_ext(pco.gobiol.file))
pco.gobiol.final <- slice(pco.gobiol.filter, c(4, 21, 92))
pco.gobiol.final %<>% rbind(pco.gobiol[119,])

pco.gomolec.file <- "./enrichr/pco.mean/GO_Molecular_Function_2015.xlsx" 
pco.gomolec <- read.xlsx(pco.gomolec.file)
pco.gomolec$Database <- "GO Molecular Function"
pco.gomolec$Num.Genes <- map(pco.gomolec$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pco.gomolec.final <- slice(pco.molec, 11)

pco.reactome.file <- "./enrichr/pco.mean/Reactome_2016.xlsx" 
pco.reactome <- read.xlsx(pco.reactome.file)
pco.reactome$Database <- "Reactome"
pco.reactome$Num.Genes <- map(pco.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pco.reactome.filter <- filter(pco.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(pco.reactome.filter, intensities.dists$Symbol, file_path_sans_ext(pco.reactome.file))
pco.reactome.final <- slice(pco.reactome.filter, c(1, 17, 47))

pco.kegg <- read.xlsx("./enrichr/pco.mean/KEGG_2016.xlsx") %>% slice(9)
pco.kegg$Num.Genes <- map(pco.kegg$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pco.kegg$Database <- "KEGG"

pco.enrichr <- rbind(pco.gobiol.final, pco.gomolec.final, pco.reactome.final, pco.kegg)
gen.enrichrplot(pco.enrichr, "pco.enrichr")

cc.gobiol.file <- "./enrichr/cc.mean/GO_Biological_Process_2015.xlsx" 
cc.gobiol <- read.xlsx(cc.gobiol.file)
cc.gobiol$Database <- "GO Biological Process"
cc.gobiol$Num.Genes <- map(cc.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
cc.gobiol.filter <- filter(cc.gobiol, Num.Genes > 4) %>% filter(P.value < 0.005)
get.kappa.cluster(cc.gobiol.filter, intensities.dists$Symbol, file_path_sans_ext(cc.gobiol.file))
cc.gobiol.final <- slice(cc.gobiol, c(60, 29, 116, 20, 35))

cc.gomolec.file <- "./enrichr/cc.mean/GO_Molecular_Function_2015.xlsx" 
cc.gomolec <- read.xlsx(cc.gomolec.file)
cc.gomolec$Database <- "GO Molecular Function"
cc.gomolec$Num.Genes <- map(cc.gomolec$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
cc.gomolec.final <- slice(cc.gomolec, 10)

cc.reactome.file <- "./enrichr/cc.mean/Reactome_2016.xlsx" 
cc.reactome <- read.xlsx(cc.reactome.file)
cc.reactome$Database <- "Reactome"
cc.reactome$Num.Genes <- map(cc.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
cc.reactome.filter <- filter(cc.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(cc.reactome.filter, intensities.dists$Symbol, file_path_sans_ext(cc.reactome.file))
cc.reactome.final <- slice(cc.reactome.filter, 1)

cc.kegg <- read.xlsx("./enrichr/cc.mean/KEGG_2016.xlsx") %>% slice(26)
cc.kegg$Database <- "KEGG"
cc.kegg$Num.Genes <- map(cc.kegg$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)

cc.enrichr <- rbind(cc.gobiol.final, cc.gomolec.final, cc.reactome.final)
gen.enrichrplot(cc.enrichr, "cc.enrichr")
