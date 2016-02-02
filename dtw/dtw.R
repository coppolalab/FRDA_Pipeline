#String operations
library(stringr)

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
library(peer)
library(limma)

#Longitudinal analysis
library(dtw)

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
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doBy)

gen.median <- function(dataset)
{
    intensities.median <- by(dataset, factor(dataset$variable), select, -PIDN, -variable, -Status) %>% lapply(Compose(as.matrix, colMedians)) %>% reduce(rbind) %>% data.frame
    intensities.median$Symbol <- unique(dataset$variable)
    return(intensities.median)
}

gen.mean <- function(dataset)
{
    intensities.mean <- by(dataset, factor(dataset$variable), select, -PIDN, -variable, -Status) %>% lapply(Compose(as.matrix, colMeans)) %>% reduce(rbind) %>% data.frame
    intensities.mean$Symbol <- unique(dataset$variable)
    return(intensities.mean)
}

gen.dist <- function(dataset)
{
    dtws <- mapply(dtw, split(dataset[[1]], row(dataset[[1]])), split(dataset[[2]], row(dataset[[2]])))
    colnames(dtws) <- rownames(dataset)
    dists <- apply(dtws, 2, `[[`, 'normalizedDistance')
    return(dists)
}

gen.dtw <- function(dataset)
{
    combinations <- combn(dataset, 2, simplify = FALSE) %>% lapply(function(x) lapply(x, select, -Symbol) )
    dists <- lapply(combinations, gen.dist)
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

lumi.known <- readRDS.gz("../baseline_lumi/save/lumi.known.rda")
lumi.vst <- lumiT(lumi.known)

long.key <- grepl("^1$|1r|2|3|4", lumi.vst$Sample.Num)
lumi.long <- lumi.vst[,long.key]
saveRDS.gz(lumi.long, file = "./save/lumi.long.rda")

lumi.norm <- lumiN(lumi.long, method = "rsn")
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

qcsum.remove <- filter(qcsum, distance.to.sample.mean > 70)$Sample.Name %>% paste(collapse = "|")
remove.indices <- grepl(qcsum.remove, sampleNames(lumi.expr.annot))
lumi.rmout <- lumi.long[,!remove.indices]
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
reps.samples <- grepl(remove.key, sampleNames(lumi.long))
remove.all <- remove.indices | reps.samples
saveRDS.gz(remove.all, "./save/remove.all.rda")

lumi.rmreps <- lumi.long[,!remove.all]
lumi.rmreps.norm <- lumiN(lumi.rmreps, method = "rsn") #Normalize with robust spline regression
#lumi.rmreps.qual <- lumiQ(lumi.rmreps.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.rmreps.cutoff <- detectionCall(lumi.rmreps.norm) #Get the count of probes which passed the detection threshold per sample
lumi.rmreps.expr <- lumi.rmreps.norm[which(lumi.rmreps.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
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

source("../common_functions.R")

gen.peer(8, exprs(lumi.combat), TRUE, model.combat)
model.PEER_covariate <- read_csv("./factor_8.csv") %>% select(-(X1:X6))
rownames(model.PEER_covariate) <- colnames(lumi.combat)
colnames(model.PEER_covariate) <- paste("X", 1:ncol(model.PEER_covariate), sep = "")

targets1.gaa <- select(pData(lumi.combat), Sample.Name, GAA1) %>% filter(!is.na(GAA1))
cor.gaa <- gen.cor(model.PEER_covariate, targets1.gaa)

targets1.onset <- select(pData(lumi.combat), Sample.Name, Onset) %>% filter(!is.na(Onset)) 
cor.onset <- gen.cor(model.PEER_covariate, targets1.onset)

PEER.traits.all <- cbind(cor.gaa, cor.onset) %>% data.frame
PEER.traits.pval <- select(PEER.traits.all, contains("p.value")) %>% as.matrix
PEER.traits.cor <- select(PEER.traits.all, -contains("p.value")) %>% as.matrix

text.matrix.PEER <- paste(signif(PEER.traits.cor, 2), '\n(', signif(PEER.traits.pval, 1), ')', sep = '')
dim(text.matrix.PEER) <- dim(PEER.traits.cor)
gen.text.heatmap(PEER.traits.cor, text.matrix.PEER, colnames(PEER.traits.cor), rownames(PEER.traits.cor), "", "PEER factor-trait relationships")

PEER.trait.out <- data.frame(Factor = rownames(PEER.traits.cor), PEER.traits.cor, PEER.traits.pval)
write_csv(PEER.trait.out, "PEER_trait_cor.csv")

PEER.weights.sums <- colSums(abs(model.PEER_covariate)) %>% data.frame
PEER.weights.sums$Factor <- 1:nrow(PEER.weights.sums)
colnames(PEER.weights.sums)[1] <- "Weight"

p <- ggplot(PEER.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
p <- p + theme_bw() + xlab("Factor") + ylab("Weight")
CairoPDF("./PEER_weights", height = 4, width = 6)
print(p)
dev.off()

#Removing effects of covariates + PEER factors  !!DO NOT USE FOR LINEAR MODELING WITH CONTRASTS!!
model.cov <- cbind(Male = model.sex.reduce, Age = as.numeric(lumi.combat$Draw.Age), RIN = lumi.combat$RIN)
model.full.cov <- cbind(model.cov, model.PEER_covariate)
export.expr <- removeBatchEffect(exprs(lumi.combat), covariates = model.full.cov, design = model.status)
export.lumi <- lumi.combat
exprs(export.lumi) <- export.expr
saveRDS.gz(export.lumi, file = "./save/export.lumi.rda")

batch.colors <- c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green", "blue4", "red4")
gen.boxplot("baseline_intensity_corrected.jpg", export.lumi, batch.colors, "Covariate-corrected intensity", "Intensity")

PIDN.long <- filter(pData(export.lumi), Sample.Num == "4")$PIDN 
export.lumi$Sample.Num[export.lumi$Sample.Num == "1r"] <- 1
export.lumi$Sample.Num %<>% as.character %>% as.numeric

targets.long <- filter(pData(export.lumi), PIDN %in% PIDN.long) 
targets.nums <- by(targets.long, targets.long$PIDN, nrow) %>% as.list %>% melt %>% data.frame 
missing.PIDNs <- filter(targets.nums, value < 4)$L1 
targets.final <- filter(targets.long, !grepl(paste("^", missing.PIDNs, "$", sep = ""), PIDN))
long.index <- sampleNames(export.lumi) %in% targets.final$Sample.Name 

lumi.final <- export.lumi[,long.index]
fdata.first <- fData(lumi.final)
fdata.first$nuID <- rownames(fdata.first)

lumi.final <- lumi.final[!is.na(fdata.first$SYMBOL),]
fdata <- fData(lumi.final)
lumi.exprs <- exprs(lumi.final)
lumi.collapse <- collapseRows(lumi.exprs, factor(fdata$SYMBOL), rownames(lumi.exprs), method = "function", methodFunction = colMeans)
lumi.exprs.collapse <- lumi.collapse$datETcollapsed %>% t

saveRDS.gz(lumi.final, "./save/lumi.final.rda")
saveRDS.gz(lumi.exprs.collapse, "./save/lumi.exprs.collapse.rda")
saveRDS.gz(fdata, "./save/fdata.rda")

intensities.long <- data.frame(Sample.Num = lumi.final$Sample.Num, PIDN = lumi.final$PIDN, Status = lumi.final$Status, lumi.exprs.collapse) 
intensities.melt <- melt(intensities.long, id.vars = c("Sample.Num", "PIDN", "Status"))
intensities.summed <- dcast(intensities.melt, PIDN + variable + Status ~ Sample.Num, value.var = "value")
saveRDS.gz(intensities.summed, "./save/intensities.summed.rda")

source("../../code/GO/enrichr.R")

intensities.medians <- by(intensities.summed, factor(intensities.summed$Status), gen.median)
saveRDS.gz(intensities.medians, "./save/intensities.medians.rda")
intensities.dists <- gen.dtw(intensities.medians) %>% reduce(cbind) %>% data.frame
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

dtw.cc.submit <- dtw.cc[1:500,]
dtw.pca.submit <- dtw.pca[1:500,]
dtw.pco.submit <- dtw.pco[1:500,]

get.stringdb(dtw.cc.submit, "cc.median", "cc")
get.stringdb(dtw.pca.submit, "pca.median", "pca")
get.stringdb(dtw.pco.submit, "pco.median", "pco")

enrichr.terms <- list("GO_Biological_Process", "GO_Molecular_Function", "KEGG_2015", "WikiPathways_2015", "Reactome_2015", "BioCarta_2015", "PPI_Hub_Proteins", "HumanCyc", "NCI-Nature", "Panther") 
dtw.cc.enrichr <- map(enrichr.terms, get.enrichrdata, dtw.cc.submit, FALSE)
names(dtw.cc.enrichr) <- enrichr.terms
map(names(dtw.cc.enrichr), enrichr.wkbk, dtw.cc.enrichr, "./enrichr/cc.median")

dtw.pca.enrichr <- map(enrichr.terms, get.enrichrdata, dtw.pca.submit, FALSE)
names(dtw.pca.enrichr) <- enrichr.terms
map(names(dtw.pca.enrichr), enrichr.wkbk, dtw.pca.enrichr, "./enrichr/pca.median")

dtw.pco.enrichr <- map(enrichr.terms, get.enrichrdata, dtw.pco.submit, FALSE)
names(dtw.pco.enrichr) <- enrichr.terms
map(names(dtw.pco.enrichr), enrichr.wkbk, dtw.pco.enrichr, "./enrichr/pco.median")

intensities.means <- by(intensities.summed, factor(intensities.summed$Status), gen.mean)
saveRDS.gz(intensities.means, "./save/intensities.means.rda")
intensities.dists.m <- gen.dtw(intensities.means) %>% reduce(cbind) %>% data.frame
rownames(intensities.dists.m) <- colnames(lumi.exprs.collapse)
colnames(intensities.dists.m) <- coltitles
intensities.dists.m$Symbol <- rownames(intensities.dists.m)

dtw.cc.m <- arrange(intensities.dists.m, desc(Carrier.vs.Control))
dtw.pca.m <- select(intensities.dists.m, Symbol, Patient.vs.Carrier, Patient.vs.Control, Carrier.vs.Control) %>% arrange(desc(Patient.vs.Carrier))
dtw.pco.m <- select(intensities.dists.m, Symbol, Patient.vs.Control, Patient.vs.Carrier, Carrier.vs.Control) %>% arrange(desc(Patient.vs.Control))

saveRDS.gz(dtw.cc.m, "./save/dtw.cc.rda")
saveRDS.gz(dtw.pca.m, "./save/dtw.pca.rda")
saveRDS.gz(dtw.pco.m, "./save/dtw.pco.rda")
gen.small.workbook(dtw.cc.m, "carrier.control.mean.xlsx")
gen.small.workbook(dtw.pca.m, "patient.carrier.mean.xlsx")
gen.small.workbook(dtw.pco.m, "patient.control.mean.xlsx")

dtw.cc.m.submit <- dtw.cc.m[1:500,]
dtw.pca.m.submit <- dtw.pca.m[1:500,]
dtw.pco.m.submit <- dtw.pco.m[1:500,]


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

get.pvalue <- function(test.value, dist.vector)
{
    num.above <- which(dist.vector > test.value) %>% length
    return(num.above / length(dist.vector))
}

permut.pca <- readRDS.gz("./save/permut.df.pca.rda")
permut.pco <- readRDS.gz("./save/permut.df.pco.rda")
permut.cc <- readRDS.gz("./save/permut.df.cc.rda")

dtw.sorted <- arrange(dtw.pca, Symbol)
pvalue.pca <- map2_dbl(dtw.sorted$Patient.vs.Carrier, split(permut.pca, 1:nrow(permut.pca)), get.pvalue)
pvalue.pco <- map2_dbl(dtw.sorted$Patient.vs.Control, split(permut.pco, 1:nrow(permut.pco)), get.pvalue)
pvalue.cc <- map2_dbl(dtw.sorted$Carrier.vs.Control, split(permut.cc, 1:nrow(permut.cc)), get.pvalue)

sig.pca <- slice(dtw.sorted, which(pvalue.pca < 0.05))
sig.pca$pvalue <- pvalue.pca[which(pvalue.pca < 0.05)]
sig.pca %<>% arrange(pvalue, desc(Patient.vs.Carrier))
write.xlsx()

sig.pco <- slice(dtw.sorted, which(pvalue.pco < 0.05))
sig.pco$pvalue <- pvalue.pco[which(pvalue.pco < 0.05)]
sig.pco %<>% arrange(pvalue, desc(Patient.vs.Control))

sig.cc <- slice(dtw.sorted, which(pvalue.cc < 0.05))
sig.cc$pvalue <- pvalue.cc[which(pvalue.cc < 0.05)]
sig.cc %<>% arrange(pvalue, desc(Carrier.vs.Control))

get.stringdb(sig.cc, "cc.mean", "cc", 900)
get.stringdb(sig.pca, "pca.mean", "pca", 900)
get.stringdb(sig.pco, "pco.mean", "pco", 900)

dtw.cc.m.enrichr <- map(enrichr.terms, get.enrichrdata, sig.cc, FALSE)
names(dtw.cc.m.enrichr) <- enrichr.terms
map(names(dtw.cc.m.enrichr), enrichr.wkbk, dtw.cc.m.enrichr, "./enrichr/cc.mean")

dtw.pca.m.enrichr <- map(enrichr.terms, get.enrichrdata, sig.pca, FALSE)
names(dtw.pca.m.enrichr) <- enrichr.terms
map(names(dtw.pca.m.enrichr), enrichr.wkbk, dtw.pca.m.enrichr, "./enrichr/pca.mean")

dtw.pco.m.enrichr <- map(enrichr.terms, get.enrichrdata, sig.cc, FALSE)
names(dtw.pco.m.enrichr) <- enrichr.terms
map(names(dtw.pco.m.enrichr), enrichr.wkbk, dtw.pco.m.enrichr, "./enrichr/pco.mean")
