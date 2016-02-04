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
library(timecourse)

#Plotting
library(Cairo)
library(WGCNA)
library(heatmap.plus)
library(flashClust)
enableWGCNAThreads()

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(doBy)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(vadr)

#Notes
#Select only patients
#Drop timepoints greater than 5
#collapse rows
#factor patient IDs into replicates vector
#extract time point

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

lumi.known <- readRDS.gz("../baseline_lumi/save/lumi.known.rda")
lumi.vst <- lumiT(lumi.known)

long.key <- grepl("^1$|1r|2|3|4|5", lumi.vst$Sample.Num)
patient.key <- grepl("Patient", lumi.vst$Status)
combined.key <- long.key & patient.key
lumi.long <- lumi.vst[,combined.key]
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
#Must decide if renormalization in necessary!

lumi.rmreps.norm <- lumiN(lumi.rmreps, method = "rsn") #Normalize with robust spline regression
#lumi.rmreps.qual <- lumiQ(lumi.rmreps.norm, detectionTh = 0.01) #The detection threshold can be adjusted here.  It is probably inadvisable to use anything larger than p < 0.05
lumi.rmreps.cutoff <- detectionCall(lumi.rmreps.norm) #Get the count of probes which passed the detection threshold per sample
lumi.rmreps.expr <- lumi.rmreps.norm[which(lumi.rmreps.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi.rmreps <- getSYMBOL(rownames(lumi.rmreps.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.rmreps.annot <- lumi.rmreps.expr[!symbols.lumi.rmreps,] #Drop any probe which is not annotated
saveRDS.gz(lumi.rmreps.annot, file = "./save/lumi.rmreps.annot.rda")

lumi.rmreps.annot$Sample.Num %<>% str_replace("1r", "1") %>% as.numeric
#lumi.baseline <- lumi.rmreps.annot[,lumi.rmreps.annot$Sample.Num == "1"]
#lumi.pointtwo <- lumi.rmreps.annot[,lumi.rmreps.annot$Sample.Num == "2"]
#lumi.pointthree <- lumi.rmreps.annot[,lumi.rmreps.annot$Sample.Num == "3"]
#lumi.pointfour <- lumi.rmreps.annot[,lumi.rmreps.annot$Sample.Num == "4"]
#lumi.pointfive <- lumi.rmreps.annot[,lumi.rmreps.annot$Sample.Num == "5"]

#Needed to deal with nightmare scenario in which one or two batches per time point have only one array, meaning they can't be used with batch correction.
match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|")) #Composed function to wrap batch numbers in '^' and '$' for exact match
#For each time point, find the batch numbers with only one array
#baseline.onebatch <- vadr::chain(factor(lumi.baseline$Batch), summary, (. == 1), which, names, match.exact) 
#pointtwo.onebatch <- vadr::chain(factor(lumi.pointtwo$Batch), summary, (. == 1), which, names, match.exact)
#pointthree.onebatch <- vadr::chain(factor(lumi.pointthree$Batch), summary, (. == 1), which, names, match.exact)
#pointfour.onebatch <- vadr::chain(factor(lumi.pointfour$Batch), summary, (. == 1), which, names, match.exact)

#Then find the PIDNs corresponding to those arrays
#baseline.drop <- filter(pData(lumi.baseline), grepl(baseline.onebatch, Batch))$PIDN %>% match.exact
#pointtwo.drop <- filter(pData(lumi.pointtwo), grepl(pointtwo.onebatch, Batch))$PIDN %>% match.exact
#pointthree.drop <- filter(pData(lumi.pointthree), grepl(pointthree.onebatch, Batch))$PIDN %>% match.exact
#all.drop <- paste(baseline.drop, pointtwo.drop, pointthree.drop, sep = "|")

#lumi.baseline <- lumi.baseline[,!grepl(all.drop, lumi.baseline$PIDN)]
#lumi.pointtwo <- lumi.pointtwo[,!grepl(all.drop, lumi.pointtwo$PIDN)]
#lumi.pointthree <- lumi.pointthree[,!grepl(all.drop, lumi.pointthree$PIDN)]

source("../common_functions.R")
#correct.covar <- function(lumi.object, timepoint, peer.factors)
#{
    model.sex <- model.matrix( ~ 0 + factor(lumi.rmreps.annot$Sex) )
    colnames(model.sex) <- c("Female", "Male")
    model.sex.reduce <- model.sex[,-1]

    model.combat <- cbind(Male = model.sex.reduce, Age = as.numeric(lumi.rmreps.annot$Draw.Age), RIN = lumi.rmreps.annot$RIN)

    expr.combat <- ComBat(dat = exprs(lumi.rmreps.annot), batch = factor(lumi.rmreps.annot$Batch), mod = model.combat)
    lumi.combat <- lumi.rmreps.annot
    exprs(lumi.combat) <- expr.combat
    saveRDS.gz(lumi.combat, "./save/lumi.combat.rda")
    #batch.colors <- c("black","navy","blue","red","orange","cyan","tan","purple","lightcyan","lightyellow","darkseagreen","brown","salmon","gold4","pink","green")
    #gen.boxplot("baseline_intensity_corrected.jpg", lumi.combat, batch.colors, "Covariate-corrected intensity", "Intensity")

    gen.peer(8, exprs(lumi.combat), TRUE, model.combat)
    model.PEER_covariate <- read_csv("factor_8.csv") %>% select(-(X1:X4))
    rownames(model.PEER_covariate) <- colnames(lumi.combat)
    colnames(model.PEER_covariate) <- paste("X", 1:ncol(model.PEER_covariate), sep = "")

    #targets1.gaa <- select(pData(lumi.combat), Sample.Name, GAA1) %>% filter(!is.na(GAA1))
    #cor.gaa <- gen.cor(model.PEER_covariate, targets1.gaa)

    #targets1.onset <- select(pData(lumi.combat), Sample.Name, Onset) %>% filter(!is.na(Onset)) 
    #cor.onset <- gen.cor(model.PEER_covariate, targets1.onset)

    #PEER.traits.all <- cbind(cor.gaa, cor.onset) %>% data.frame
    #PEER.traits.pval <- select(PEER.traits.all, contains("p.value")) %>% as.matrix
    #PEER.traits.cor <- select(PEER.traits.all, -contains("p.value")) %>% as.matrix

    #text.matrix.PEER <- paste(signif(PEER.traits.cor, 2), '\n(', signif(PEER.traits.pval, 1), ')', sep = '')
    #dim(text.matrix.PEER) <- dim(PEER.traits.cor)
    #gen.text.heatmap(PEER.traits.cor, text.matrix.PEER, colnames(PEER.traits.cor), rownames(PEER.traits.cor), "", "PEER factor-trait relationships")

    #PEER.trait.out <- data.frame(Factor = rownames(PEER.traits.cor), PEER.traits.cor, PEER.traits.pval)
    #write_csv(PEER.trait.out, "PEER_trait_cor.csv")

    PEER.weights <- read_csv("./weight_8.csv") %>% select(-(X1:X4))
    PEER.weights.sums <- colSums(abs(PEER.weights)) %>% data.frame
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
    #saveRDS.gz(export.lumi, file = "./save/export.lumi.rda")
    return(export.lumi)
#}

lumi.baseline.correct <- correct.covar(lumi.baseline, 1, 8)
lumi.baseline.correct <- correct.covar(lumi.baseline, 1, 8)
lumi.baseline.correct <- correct.covar(lumi.baseline, 1, 8)
lumi.baseline.correct <- correct.covar(lumi.baseline, 1, 8)

PIDN.long <- filter(pData(export.lumi), Sample.Num == "5")$PIDN 
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
