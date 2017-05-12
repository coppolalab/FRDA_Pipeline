library(lumi)
library(WGCNA) #for fastcor
library(BayesFactor)
library(siggenes)

library(Cairo)
library(scales)
library(openxlsx)

library(plyr) #only for revalue
library(magrittr)
library(stringr)
library(tidyverse)

source("../common_functions.R")

GetHyper <- function(colname2, colname1, dataset2, dataset1, symbolname, threshold = 0.9) {
    colname1.format <- str_c(colname1, ".Sign.Posterior")
    colname2.format <- str_c(colname2, ".Sign.Posterior") 
    subset1 <- select_(dataset1, symbolname, colname1.format)
    subset2 <- select_(dataset2, symbolname, colname2.format)
    subset1.sig.up <- filter_(subset1, str_c(colname1.format, " > ", threshold))
    subset2.sig.up <- filter_(subset2, str_c(colname2.format, " > ", threshold))
    subset1.sig.down <- filter_(subset1, str_c(colname1.format, " < ", -threshold))
    subset2.sig.down <- filter_(subset2, str_c(colname2.format, " < ", -threshold))

    #Up
    up.intersect <- intersect(subset1.sig.up$Symbol, subset2.sig.up$Symbol)
    up.up <- length(up.intersect)

    up.all <- length(subset1.sig.up$Symbol) - up.up
    all.up <- intersect(subset2.sig.up$Symbol, subset1$Symbol) %>% length

    all.all <- nrow(subset1) - up.all - all.up + up.up
    up.column1 <- c(up.up, up.all)
    up.column2 <- c(all.up, all.all)
    up.matrix <- cbind(up.column1, up.column2)
    up.bf <- contingencyTableBF(up.matrix, sampleType = "hypergeom") %>% extractBF %>% extract2("bf")

    if(length(up.intersect) > 0){
        up.final <- c(signif(up.bf, 3), str_c(up.intersect, collapse = ","))
    } else {
        up.final <- c(signif(up.bf, 3), "None")
    }

    #Down
    down.intersect <- intersect(subset1.sig.down$Symbol, subset2.sig.down$Symbol)
    down.down <- length(down.intersect)

    down.all <- length(subset1.sig.down$Symbol) - down.down
    all.down <- intersect(subset2.sig.down$Symbol, subset1$Symbol) %>% length

    all.all <- nrow(subset1) - down.all - all.down + down.down
    down.column1 <- c(down.down, down.all)
    down.column2 <- c(all.down, all.all)
    down.matrix <- cbind(down.column1, down.column2)
    down.bf <- contingencyTableBF(down.matrix, sampleType = "hypergeom") %>% extractBF %>% extract2("bf")

    if(length(down.intersect) > 0){
        down.final <- c(signif(down.bf, 3), str_c(down.intersect, collapse = ","))
    } else {
        down.final <- c(signif(down.bf, 3), "None")
    }

    return.df <- rbind(up.final, down.final) %>% data.frame
    colnames(return.df) <- c("Bayes.Factor", "Genes")
    return.df$Direction <- c("Up", "Down")
    return.df$Bayes.Factor %<>% as.character %>% as.numeric
    
    return.df
}

MapHyper <- function(colname1, colname.list, dataset1, dataset2, symbolname, threshold = 0.9){
    return.hyper <- map(colname.list, GetHyper, colname1, dataset2, dataset1, symbolname, threshold)
    names(return.hyper) <- colname.list
    return.hyper
}

GetBF <- function(hyper.list, direction) {
    dir.string <- str_c("Direction == '", direction, "'")
    hyper.bf <- map(hyper.list, filter_, dir.string) %>% map_dbl(extract2, "Bayes.Factor") 
    hyper.bf
}

GetCounts <- function(hyper.list, direction) {
    dir.string <- str_c("Direction == '", direction, "'")
    hyper.genes <- map(hyper.list, filter_, dir.string) %>% map(extract2, "Genes") %>% map(str_split, ",") %>% map(unlist) 
    hyper.empty <- map(hyper.genes, str_detect, "None") %>% map_lgl(reduce, or)
    hyper.genecount <- map_int(hyper.genes, length) 
    hyper.genecount[hyper.empty] <- 0
    hyper.genecount
}

#Read in HomoloGene table
mouse.homology <- read_tsv("./HOM_MouseHumanSequence.rpt.txt") %>% data.frame %>% select(HomoloGene.ID, NCBI.Taxon.ID, Symbol)
mouse.only <- filter(mouse.homology, NCBI.Taxon.ID == 10090)
human.only <- filter(mouse.homology, NCBI.Taxon.ID == 9606)

#Read in differential expression data
human.groups <- c("cc", "pca", "pco")

ebam.pco.df <- ReadRDSgz("../baseline_lumi/save/ebam.pco.df.rda") 
colnames(ebam.pco.df)[1:2] <- str_c("pco.", colnames(ebam.pco.df)[1:2])
ebam.pco.df$pco.Sign.Posterior <- sign(ebam.pco.df$pco.Z.score) * ebam.pco.df$pco.Posterior

ebam.cc.df <- ReadRDSgz("../baseline_lumi/save/ebam.cc.df.rda")
colnames(ebam.cc.df)[1:2] <- str_c("cc.", colnames(ebam.cc.df)[1:2])
ebam.cc.df$cc.Sign.Posterior <- sign(ebam.cc.df$cc.Z.score) * ebam.cc.df$cc.Posterior

ebam.pca.df <- ReadRDSgz("../baseline_lumi/save/ebam.pca.df.rda")
colnames(ebam.pca.df)[1:2] <- str_c("pca.", colnames(ebam.pca.df)[1:2])
ebam.pca.df$pca.Sign.Posterior <- sign(ebam.pca.df$pca.Z.score) * ebam.pca.df$pca.Posterior

all.human.genes <- left_join(ebam.pca.df, ebam.pco.df) %>% left_join(ebam.cc.df)
all.human.reduce <- select(all.human.genes, Symbol, dplyr::contains("Sign.Posterior"))

#Read in GAA expansion data
ebam.gaa.df <- ReadRDSgz("../WGCNA_GAA/ebam.gaa.df.rda") %>% select(Z.score, Posterior)
ebam.gaa.df$Sign.Posterior <- sign(ebam.gaa.df$Z.score) * ebam.gaa.df$Posterior
colnames(ebam.gaa.df) <- str_c("gaa.", colnames(ebam.gaa.df))
ebam.gaa.df$Symbol <- rownames(ebam.gaa.df)

#Vijay study
vijay.heart.groups <- c("doxnd", "tgwt", "rescue")
vijay.heart.doxnd.df <- ReadRDSgz("../../Vijay_mouse/baseline/ebam.doxnd.df.rda") %>% select(Z.score, Posterior)
vijay.heart.doxnd.df$Sign.Posterior <- sign(vijay.heart.doxnd.df$Z.score) * vijay.heart.doxnd.df$Posterior
colnames(vijay.heart.doxnd.df) <- str_c("doxnd.", colnames(vijay.heart.doxnd.df))

vijay.heart.tgwt.df <- ReadRDSgz("../../Vijay_mouse/baseline/ebam.tgwt.df.rda") %>% select(Z.score, Posterior)
vijay.heart.tgwt.df$Sign.Posterior <- sign(vijay.heart.tgwt.df$Z.score) * vijay.heart.tgwt.df$Posterior
colnames(vijay.heart.tgwt.df) <- str_c("tgwt.", colnames(vijay.heart.tgwt.df))

vijay.heart.rescue.df <- ReadRDSgz("../../Vijay_mouse/baseline/ebam.rescue.df.rda") %>% select(Z.score, Posterior)
vijay.heart.rescue.df$Sign.Posterior <- -sign(vijay.heart.rescue.df$Z.score) * vijay.heart.rescue.df$Posterior
colnames(vijay.heart.rescue.df) <- str_c("rescue.", colnames(vijay.heart.rescue.df))

vijay.heart.all <- cbind(vijay.heart.doxnd.df, vijay.heart.tgwt.df, vijay.heart.rescue.df)
vijay.heart.all$Symbol <- rownames(vijay.heart.all) %>% toupper
vijay.heart.reduce <- select(vijay.heart.all, Symbol, dplyr::contains("Sign.Posterior"))

vijay.hearth.hyper <- map(human.groups, MapHyper, vijay.heart.groups, all.human.reduce, vijay.heart.reduce, "Symbol", 0.8)
names(vijay.hearth.hyper) <- str_c("human.", human.groups)
vijay.heartg.hyper <- map(vijay.heart.groups, GetHyper, "gaa", vijay.heart.reduce, ebam.gaa.df, "Symbol", 0.8)
names(vijay.heartg.hyper) <- str_c("human.", human.groups)

vijay.heart.up.bf <- map(vijay.hearth.hyper, GetBF, "Up") %>% reduce(rbind) 
colnames(vijay.heart.up.bf) <- str_c("vijay.heart.", colnames(vijay.heart.up.bf))
vijay.heart.down.bf <- map(vijay.hearth.hyper, GetBF, "Down") %>% reduce(rbind)
colnames(vijay.heart.down.bf) <- str_c("vijay.heart.", colnames(vijay.heart.down.bf))
vijay.heart.up.genecount <- map(vijay.hearth.hyper, GetCounts, "Up") %>% reduce(rbind)
colnames(vijay.heart.up.genecount) <- str_c("vijay.heart.", colnames(vijay.heart.up.genecount))
vijay.heart.down.genecount <- map(vijay.hearth.hyper, GetCounts, "Down") %>% reduce(rbind)
colnames(vijay.heart.down.genecount) <- str_c("vijay.heart.", colnames(vijay.heart.down.genecount))

#DRG
vijay.drg.groups <- c("doxnd", "tgwt", "rescue")
vijay.drg.doxnd.df <- ReadRDSgz("../../Vijay_mouse/baseline_drg/ebam.doxnd.df.rda") %>% select(Z.score, Posterior)
vijay.drg.doxnd.df$Sign.Posterior <- sign(vijay.drg.doxnd.df$Z.score) * vijay.drg.doxnd.df$Posterior
colnames(vijay.drg.doxnd.df) <- str_c("doxnd.", colnames(vijay.drg.doxnd.df))

vijay.drg.tgwt.df <- ReadRDSgz("../../Vijay_mouse/baseline_drg/ebam.tgwt.df.rda") %>% select(Z.score, Posterior)
vijay.drg.tgwt.df$Sign.Posterior <- sign(vijay.drg.tgwt.df$Z.score) * vijay.drg.tgwt.df$Posterior
colnames(vijay.drg.tgwt.df) <- str_c("tgwt.", colnames(vijay.drg.tgwt.df))

vijay.drg.rescue.df <- ReadRDSgz("../../Vijay_mouse/baseline_drg/ebam.rescue.df.rda") %>% select(Z.score, Posterior)
vijay.drg.rescue.df$Sign.Posterior <- -sign(vijay.drg.rescue.df$Z.score) * vijay.drg.rescue.df$Posterior
colnames(vijay.drg.rescue.df) <- str_c("rescue.", colnames(vijay.drg.rescue.df))

vijay.drg.all <- cbind(vijay.drg.doxnd.df, vijay.drg.tgwt.df, vijay.drg.rescue.df)
vijay.drg.all$Symbol <- rownames(vijay.drg.all) %>% toupper
vijay.drg.reduce <- select(vijay.drg.all, Symbol, dplyr::contains("Sign.Posterior"))

vijay.drgh.hyper <- map(human.groups, MapHyper, vijay.drg.groups, all.human.reduce, vijay.drg.reduce, "Symbol", 0.8)
names(vijay.drgh.hyper) <- str_c("human.", human.groups)
vijay.drgg.hyper <- map(vijay.drg.groups, GetHyper, "gaa", vijay.drg.reduce, ebam.gaa.df, "Symbol", 0.8)
names(vijay.drgg.hyper) <- str_c("human.", human.groups)

vijay.drg.up.bf <- map(vijay.drgh.hyper, GetBF, "Up") %>% reduce(rbind) 
colnames(vijay.drg.up.bf) <- str_c("vijay.drg.", colnames(vijay.drg.up.bf))
vijay.drg.down.bf <- map(vijay.drgh.hyper, GetBF, "Down") %>% reduce(rbind)
colnames(vijay.drg.down.bf) <- str_c("vijay.drg.", colnames(vijay.drg.down.bf))
vijay.drg.up.genecount <- map(vijay.drgh.hyper, GetCounts, "Up") %>% reduce(rbind)
colnames(vijay.drg.up.genecount) <- str_c("vijay.drg.", colnames(vijay.drg.up.genecount))
vijay.drg.down.genecount <- map(vijay.drgh.hyper, GetCounts, "Down") %>% reduce(rbind)
colnames(vijay.drg.down.genecount) <- str_c("vijay.drg.", colnames(vijay.drg.down.genecount))

#Cerebellum
vijay.cerebellum.groups <- c("doxnd", "tgwt", "rescue")
vijay.cerebellum.doxnd.df <- ReadRDSgz("../../Vijay_mouse/baseline_cerebellum/ebam.doxnd.df.rda") %>% select(Z.score, Posterior)
vijay.cerebellum.doxnd.df$Sign.Posterior <- sign(vijay.cerebellum.doxnd.df$Z.score) * vijay.cerebellum.doxnd.df$Posterior
colnames(vijay.cerebellum.doxnd.df) <- str_c("doxnd.", colnames(vijay.cerebellum.doxnd.df))

vijay.cerebellum.tgwt.df <- ReadRDSgz("../../Vijay_mouse/baseline_cerebellum/ebam.tgwt.df.rda") %>% select(Z.score, Posterior)
vijay.cerebellum.tgwt.df$Sign.Posterior <- sign(vijay.cerebellum.tgwt.df$Z.score) * vijay.cerebellum.tgwt.df$Posterior
colnames(vijay.cerebellum.tgwt.df) <- str_c("tgwt.", colnames(vijay.cerebellum.tgwt.df))

vijay.cerebellum.rescue.df <- ReadRDSgz("../../Vijay_mouse/baseline_cerebellum/ebam.rescue.df.rda") %>% select(Z.score, Posterior)
vijay.cerebellum.rescue.df$Sign.Posterior <- -sign(vijay.cerebellum.rescue.df$Z.score) * vijay.cerebellum.rescue.df$Posterior
colnames(vijay.cerebellum.rescue.df) <- str_c("rescue.", colnames(vijay.cerebellum.rescue.df))

vijay.cerebellum.all <- cbind(vijay.cerebellum.doxnd.df, vijay.cerebellum.tgwt.df, vijay.cerebellum.rescue.df)
vijay.cerebellum.all$Symbol <- rownames(vijay.cerebellum.all) %>% toupper
vijay.cerebellum.reduce <- select(vijay.cerebellum.all, Symbol, dplyr::contains("Sign.Posterior"))

vijay.cerebellumh.hyper <- map(human.groups, MapHyper, vijay.cerebellum.groups, all.human.reduce, vijay.cerebellum.reduce, "Symbol", 0.8)
names(vijay.cerebellumh.hyper) <- str_c("human.", human.groups)
vijay.cerebellumg.hyper <- map(vijay.cerebellum.groups, GetHyper, "gaa", vijay.cerebellum.reduce, ebam.gaa.df, "Symbol", 0.8)
names(vijay.cerebellumg.hyper) <- str_c("human.", human.groups)

vijay.cerebellum.up.bf <- map(vijay.cerebellumh.hyper, GetBF, "Up") %>% reduce(rbind) 
colnames(vijay.cerebellum.up.bf) <- str_c("vijay.cerebellum.", colnames(vijay.cerebellum.up.bf))
vijay.cerebellum.down.bf <- map(vijay.cerebellumh.hyper, GetBF, "Down") %>% reduce(rbind)
colnames(vijay.cerebellum.down.bf) <- str_c("vijay.cerebellum.", colnames(vijay.cerebellum.down.bf))
vijay.cerebellum.up.genecount <- map(vijay.cerebellumh.hyper, GetCounts, "Up") %>% reduce(rbind)
colnames(vijay.cerebellum.up.genecount) <- str_c("vijay.cerebellum.", colnames(vijay.cerebellum.up.genecount))
vijay.cerebellum.down.genecount <- map(vijay.cerebellumh.hyper, GetCounts, "Down") %>% reduce(rbind)
colnames(vijay.cerebellum.down.genecount) <- str_c("vijay.cerebellum.", colnames(vijay.cerebellum.down.genecount))

#PBMC study
pbmc.groups <- c("cc", "pco")
pbmc.cc.df <- ReadRDSgz("../../pbmc/ebam.cc.df.rda") %>% select(Z.score, Posterior)
pbmc.cc.df$Sign.Posterior <- sign(pbmc.cc.df$Z.score) * pbmc.cc.df$Posterior
colnames(pbmc.cc.df) <- str_c("cc.", colnames(pbmc.cc.df))

pbmc.pco.df <- ReadRDSgz("../../pbmc/ebam.pco.df.rda") %>% select(Z.score, Posterior)
pbmc.pco.df$Sign.Posterior <- sign(pbmc.pco.df$Z.score) * pbmc.pco.df$Posterior
colnames(pbmc.pco.df) <- str_c("pco.", colnames(pbmc.pco.df))

pbmc.all <- cbind(pbmc.cc.df, pbmc.pco.df)
pbmc.all$Symbol <- rownames(pbmc.all)
pbmc.reduce <- select(pbmc.all, Symbol, dplyr::contains("Sign.Posterior"))

pbmch.hyper <- map(human.groups, MapHyper, pbmc.groups, all.human.reduce, pbmc.reduce, "Symbol", 0.8)
names(pbmch.hyper) <- str_c("human.", human.groups)
pbmcg.hyper <- map(pbmc.groups, GetHyper, "gaa", pbmc.reduce, ebam.gaa.df, "Symbol")
names(pbmcg.hyper) <- str_c("pbmc.", pbmc.groups)

pbmc.up.bf <- map(pbmch.hyper, GetBF, "Up") %>% reduce(rbind) 
colnames(pbmc.up.bf) <- str_c("pbmc.", colnames(pbmc.up.bf))
pbmc.down.bf <- map(pbmch.hyper, GetBF, "Down") %>% reduce(rbind)
colnames(pbmc.down.bf) <- str_c("pbmc.", colnames(pbmc.down.bf))
pbmc.up.genecount <- map(pbmch.hyper, GetCounts, "Up") %>% reduce(rbind)
colnames(pbmc.up.genecount) <- str_c("pbmc.", colnames(pbmc.up.genecount))
pbmc.down.genecount <- map(pbmch.hyper, GetCounts, "Down") %>% reduce(rbind)
colnames(pbmc.down.genecount) <- str_c("pbmc.", colnames(pbmc.down.genecount))

#IPSC
#Read in IPSC data
ipsc.groups <- "Healthy.vs.FRDA"
ipsc <- ReadRDSgz("../../ipsc/ebam.ipsc.df.rda") %>% select(Z.score, Posterior) 
ipsc$Sign.Posterior <- sign(ipsc$Z.score) * ipsc$Posterior
colnames(ipsc) <- str_c("ipsc", colnames(ipsc), sep = ".")
ipsc$Symbol <- rownames(ipsc)

#Compute RRHO between DE and IPSC data
ipsch.hyper <- map(human.groups, GetHyper, colname2 = "ipsc", ipsc, all.human.reduce, "Symbol", 0.80)
ipscg.hyper <- GetHyper("ipsc", "gaa", ipsc, ebam.gaa.df, "Symbol", 0.80)

ipsc.up.bf <- GetBF(ipsch.hyper, "Up")
ipsc.down.bf <- GetBF(ipsch.hyper, "Down")
ipsc.up.genecount <- GetCounts(ipsch.hyper, "Up")
ipsc.down.genecount <- GetCounts(ipsch.hyper, "Down")

#Van Houten
#Read in Van Houten data
vanhouten.groups <- "pco"
vanhouten <- ReadRDSgz("../../VanHauten/baseline/ebam.vanhouten.df.rda") %>% select(Z.score, Posterior)
vanhouten$Sign.Posterior <- sign(vanhouten$Z.score) * vanhouten$Posterior
colnames(vanhouten) <- str_c("vanhouten", colnames(vanhouten), sep = ".")
vanhouten$Symbol <- rownames(vanhouten)

vanhoutenh.hyper <- map(human.groups, GetHyper, colname2 = "vanhouten", vanhouten, all.human.reduce, "Symbol", 0.8)
vanhouteng.hyper <- GetHyper("vanhouten", "gaa", vanhouten, ebam.gaa.df, "Symbol", 0.8)

vanhouten.up.bf <- GetBF(vanhoutenh.hyper, "Up")
vanhouten.down.bf <- GetBF(vanhoutenh.hyper, "Down")
vanhouten.up.genecount <- GetCounts(vanhoutenh.hyper, "Up")
vanhouten.down.genecount <- GetCounts(vanhoutenh.hyper, "Down")

#Combine all results
bf.up <- data.frame(vijay.heart.up.bf, vijay.drg.up.bf, vijay.cerebellum.up.bf, pbmc.up.bf, FRDA.vs.Healthy.vanhouten = vanhouten.up.bf) %>% log10 %>% signif(3) %>% slice(-1)
bf.down <- data.frame(vijay.heart.down.bf, vijay.drg.down.bf, vijay.cerebellum.down.bf, pbmc.down.bf, FRDA.vs.Healthy.vanhouten = vanhouten.down.bf) %>% log10 %>% signif(3) %>% slice(-1)
count.up <- data.frame(vijay.heart.up.genecount, vijay.drg.up.genecount, vijay.cerebellum.up.genecount, pbmc.up.genecount, FRDA.vs.Healthy.vanhouten = vanhouten.up.genecount) %>% slice(-1)
count.down <- data.frame(vijay.heart.down.genecount, vijay.drg.down.genecount, vijay.cerebellum.down.genecount, pbmc.down.genecount, FRDA.vs.Healthy.vanhouten = vanhouten.down.genecount) %>% slice(-1)

#reduced
bf.up <- data.frame(pbmc.up.bf, FRDA.vs.Healthy.vanhouten = vanhouten.up.bf) %>% log10 %>% signif(3) %>% slice(-1)
bf.down <- data.frame(pbmc.down.bf, FRDA.vs.Healthy.vanhouten = vanhouten.down.bf) %>% log10 %>% signif(3) %>% slice(-1)
count.up <- data.frame(pbmc.up.genecount, FRDA.vs.Healthy.vanhouten = vanhouten.up.genecount) %>% slice(-1)
count.down <- data.frame(pbmc.down.genecount, FRDA.vs.Healthy.vanhouten = vanhouten.down.genecount) %>% slice(-1)

up.format <- str_c(as.matrix(count.up), '\n(', as.matrix(bf.up), ')') 
dim(up.format) <- dim(bf.up)
colnames(up.format) <- colnames(bf.up)
up.format.df <- data.frame(up.format, stringsAsFactors = FALSE) 

down.format <- str_c(as.matrix(count.down), '\n(', as.matrix(bf.down), ')')
dim(down.format) <- dim(bf.down)
colnames(down.format) <- colnames(bf.down)
down.format.df <- data.frame(down.format, stringsAsFactors = FALSE) 

up.format.df$Human.Comparison <- c("Patient vs. Carrier (Upregulated)", "Patient vs. Control (Upregulated)")
#up.format.df$Human.Comparison %<>% factor(levels = up.format.df$Human.Comparison)
up.plot <- gather(up.format.df, Data.Comparison, Format.BF, -Human.Comparison)
up.plot$Data.Comparison %<>% factor %>% revalue(c(vijay.heart.doxnd = "RNAi mouse heart (DOX vs. ND)", vijay.heart.tgwt = "RNAi mouse heart (Tg vs WT)", vijay.heart.rescue = "RNAi mouse heart (Rescue)", vijay.drg.doxnd = "RNAi mouse DRG (DOX vs. ND)", vijay.drg.tgwt = "RNAi mouse DRG (Tg vs WT)", vijay.drg.rescue = "RNAi mouse DRG (Rescue)", vijay.cerebellum.doxnd = "RNAi mouse cerebellum (DOX vs. ND)", vijay.cerebellum.tgwt = "RNAi mouse cerebellum (Tg vs WT)", vijay.cerebellum.rescue = "RNAi mouse cerebellum (Rescue)", pbmc.pco = "GSE30933 (Patient vs. Control)", pbmc.cc = "GSE30933 (Carrier vs. Control)", FRDA.vs.Healthy.vanhouten = "GSE11204 (FRDA vs. Control)"))
up.plot$Data.Comparison %<>% factor(levels = unique(up.plot$Data.Comparison))
up.plot$Bayes.Factor <- as.matrix(bf.up) %>% as.numeric

down.format.df$Human.Comparison <- c("Patient vs. Carrier (Downregulated)", "Patient vs. Control (Downregulated)")
#down.format.df$Human.Comparison %<>% factor(levels = down.format.df$Human.Comparison)
down.plot <- gather(down.format.df, Data.Comparison, Format.BF, -Human.Comparison)
down.plot$Data.Comparison %<>% factor %>% revalue(c(vijay.heart.doxnd = "RNAi mouse heart (DOX vs. ND)", vijay.heart.tgwt = "RNAi mouse heart (Tg vs WT)", vijay.heart.rescue = "RNAi mouse heart (Rescue)", vijay.drg.doxnd = "RNAi mouse DRG (DOX vs. ND)", vijay.drg.tgwt = "RNAi mouse DRG (Tg vs WT)", vijay.drg.rescue = "RNAi mouse DRG (Rescue)", vijay.cerebellum.doxnd = "RNAi mouse cerebellum (DOX vs. ND)", vijay.cerebellum.tgwt = "RNAi mouse cerebellum (Tg vs WT)", vijay.cerebellum.rescue = "RNAi mouse cerebellum (Rescue)", pbmc.pco = "GSE30933 (Patient vs. Control)", pbmc.cc = "GSE30933 (Carrier vs. Control)", FRDA.vs.Healthy.vanhouten = "GSE11204 (FRDA vs. Control)"))
down.plot$Data.Comparison %<>% factor(levels = unique(down.plot$Data.Comparison))
down.plot$Bayes.Factor <- as.matrix(bf.down) %>% as.numeric

combined.plot <- rbind(up.plot, down.plot)
combined.plot$Human.Comparison %<>% factor(levels = c(down.format.df$Human.Comparison, up.format.df$Human.Comparison))

#CairoPDF("up.heatmap", width = 10, height = 5, bg = "transparent")
#p <- ggplot(rbind(up.plot, down.plot), aes(Data.Comparison, Human.Comparison, label = Format.BF, fill = Bayes.Factor)) + geom_raster() + geom_text() + theme_bw()
#p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
#p <- p + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p <- p + theme(panel.background = element_blank(), panel.border = element_blank())
#p <- p + theme(plot.background = element_blank(), plot.title = element_text(hjust = 0.5))
#p <- p + theme(axis.text.x = element_blank(), legend.position = "none")
#p <- p + scale_fill_gradient(low = "white", high = muted('red'))
##p <- p + ggtitle("Up")
#print(p)
#dev.off()

CairoPDF("combined.heatmap", width = 6.0, height = 3.5, bg = "transparent")
p <- ggplot(combined.plot, aes(Data.Comparison, Human.Comparison, label = Format.BF, fill = Bayes.Factor)) + geom_raster() + geom_text() + theme_bw()
p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(panel.background = element_blank(), panel.border = element_blank())
p <- p + theme(plot.background = element_blank(), plot.title = element_text(hjust = 0.5))
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + scale_fill_gradient(low =  "white", high = muted('red'), guide = guide_legend(title = expression(atop(Log[10], 'Bayes Factor'))))
print(p)
dev.off()

up.plot.reduce <- filter(up.plot, !grepl("Chandran mouse cerebellum|Chandran mouse DRG", Data.Comparison))
CairoPDF("up.heatmap.reduce", width = 5.5, height = 2.0, bg = "transparent")
p <- ggplot(up.plot.reduce, aes(Data.Comparison, Human.Comparison, label = Format.BF, fill = Bayes.Factor)) + geom_raster() + geom_text() + theme_bw()
p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 
p <- p + theme(plot.margin = unit(c(1,1,2,1), "lines")) + theme(plot.title = element_text(hjust = 0.5))
p <- p + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(panel.background = element_blank(), panel.border = element_blank())
p <- p + theme(plot.background = element_blank(), legend.background = element_blank())
p <- p + theme(axis.text.x = element_blank(), legend.title = element_text(hjust = 0) )
p <- p + scale_fill_gradient2(low = muted('blue'), mid = "white", high = muted('red'), guide = guide_legend(title = expression(atop(Log[10], 'Bayes Factor'))))
p <- p + ggtitle("Upregulated")
print(p)
dev.off()

down.plot.reduce <- filter(down.plot, !grepl("Chandran mouse cerebellum|Chandran mouse DRG", Data.Comparison))
CairoPDF("down.heatmap.reduce", width = 5.5, height = 3.5, bg = "transparent")
p <- ggplot(down.plot.reduce, aes(Data.Comparison, Human.Comparison, label = Format.BF, fill = Bayes.Factor)) + geom_raster() + geom_text() + theme_bw()
p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + theme(plot.margin = unit(c(1,1,1,1), "lines"))
p <- p + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) + theme(plot.title = element_text(hjust = 0.5))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(panel.background = element_blank(), panel.border = element_blank())
p <- p + theme(plot.background = element_blank(), legend.background = element_blank())
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + scale_fill_gradient2(low = muted('blue'), mid = "white", high = muted('red'), guide = guide_legend(title = expression(atop(Log[10], 'Bayes Factor'))))
p <- p + ggtitle("Downregulated")
print(p)
dev.off()
