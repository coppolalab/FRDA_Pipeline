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

GetHyper <- function(colname2, colname1, dataset2, dataset1, symbolname) {
    log.column1 <- str_c("logFC.", colname1)
    log.column2 <- str_c("logFC.", colname2)
    sig.column1 <- str_c("sig.", colname1)
    sig.column2 <- str_c("sig.", colname2)

    subset1 <- select_(dataset1, symbolname, "Log.Bayes.Factor", log.column1, sig.column1)
    subset2 <- select_(dataset2, symbolname, "Log.Bayes.Factor", log.column2, sig.column2)

    subset1.sig.up <- filter_(subset1, str_c(sig.column1, " == TRUE & sign(", 
                                             log.column1, ") == 1"))
    subset2.sig.up <- filter_(subset2, str_c("Log.Bayes.Factor > 0.5 & ", 
                                             sig.column2, " == TRUE & sign(", 
                                             log.column2, ") == 1"))
    subset1.sig.down <- filter_(subset1, str_c("Log.Bayes.Factor > 0.5 & ", 
                                             sig.column1, " == TRUE & sign(", 
                                             log.column1, ") == -1"))
    subset2.sig.down <- filter_(subset2, str_c("Log.Bayes.Factor > 0.5 & ", 
                                             sig.column2, " == TRUE & sign(", 
                                             log.column2, ") == -1"))

    #Up
    up.intersect <- intersect(subset1.sig.up[[symbolname]], subset2.sig.up[[symbolname]])
    up.up <- length(up.intersect)

    up.all <- length(subset1.sig.up[[symbolname]]) - up.up
    all.up <- intersect(subset2.sig.up[[symbolname]], subset1[[symbolname]]) %>% length

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
    down.intersect <- intersect(subset1.sig.down[[symbolname]], subset2.sig.down[[symbolname]])
    down.down <- length(down.intersect)

    down.all <- length(subset1.sig.down[[symbolname]]) - down.down
    all.down <- intersect(subset2.sig.down[[symbolname]], subset1[[symbolname]]) %>% length

    all.all <- nrow(subset1) - down.all - all.down + down.down
    down.column1 <- c(down.down, down.all)
    down.column2 <- c(all.down, all.all)
    down.matrix <- cbind(down.column1, down.column2)
    down.bf <- contingencyTableBF(down.matrix, sampleType = "hypergeom") %>% 
        extractBF %>% extract2("bf")

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

MapHyper <- function(colname1, colname.list, dataset1, dataset2, symbolname){
    return.hyper <- map(colname.list, GetHyper, colname1, dataset2, dataset1, symbolname)
    names(return.hyper) <- colname.list
    return.hyper
}

GetBF <- function(hyper.list, direction) {
    dir.string <- str_c("Direction == '", direction, "'")
    if (is.null(nrow(hyper.list))) {
        hyper.bf <- map(hyper.list, filter_, dir.string) %>% 
            map_dbl(extract2, "Bayes.Factor") 
        return(hyper.bf)
    } else {
        hyper.bf <- filter_(hyper.list, dir.string) %>% 
            extract2("Bayes.Factor")
    }
}

GetCounts <- function(hyper.list, direction) {
    dir.string <- str_c("Direction == '", direction, "'")
    if (is.null(nrow(hyper.list))) {
        hyper.genes <- map(hyper.list, filter_, dir.string) %>% 
            map(extract2, "Genes") %>% 
            map(str_split, ",") %>% map(unlist) 
        hyper.empty <- map(hyper.genes, str_detect, "None") %>% 
            map_lgl(reduce, or)
        hyper.genecount <- map_int(hyper.genes, length) 
        hyper.genecount[hyper.empty] <- 0
        return(hyper.genecount)
    } else {
        hyper.genes <- filter_(hyper.list, dir.string) %>% 
            extract2("Genes") %>% str_split(",") %>% unlist
        if (hyper.genes[1] == "None") {
            return(0)
        } else {
            return(length(hyper.genes))
        }

    }
}

#Read in HomoloGene table
mouse.homology <- read_tsv("./HOM_MouseHumanSequence.rpt.txt") %>% 
    data.frame %>% select(HomoloGene.ID, NCBI.Taxon.ID, Symbol)
mouse.only <- filter(mouse.homology, NCBI.Taxon.ID == 10090) %>% select(-NCBI.Taxon.ID)
human.only <- filter(mouse.homology, NCBI.Taxon.ID == 9606) %>% select(-NCBI.Taxon.ID)

#Read in differential expression data
human.groups <- c("pca", "pco")
all.human <- ReadRDSgz("../differential_expression/save/posterior.final.df.rda")
all.human.reduce <- select(all.human, Symbol, Log.Bayes.Factor, matches("logFC"), matches("pp")) %>%
    mutate(sig.pco = (Log.Bayes.Factor > 0.5 & pp.pco > 0.95),
           sig.pca = (Log.Bayes.Factor > 0.5 & pp.pca > 0.95),
           sig.cc = (Log.Bayes.Factor > 0.5 & pp.cc > 0.95))

#Read in GAA expansion data
fds.df <- ReadRDSgz("../phenotype_regression/save/bf.fds.annot.rda") 
colnames(fds.df)[4] <- "logFC.FDS"
fds.df$sig.FDS <- fds.df$Log.Bayes.Factor > 0.5 & fds.df$pp.FS > 0.95

#Vijay study
vijay.heart <- ReadRDSgz("../../Vijay_mouse/baseline/save/heart.df.final") 
vijay.heart.groups <- c("t5", "t4", "t3")
vijay.heart.reduce <- select(vijay.heart, Symbol, Log.Bayes.Factor, matches("t3|t4|t5"))
vijay.heart.reduce$sig.t3 <- vijay.heart.reduce$Log.Bayes.Factor > 0.5 &
    vijay.heart.reduce$pp.t3 > 0.95 &
    (sign(vijay.heart.reduce$logFC.doxnd.t3) == sign(vijay.heart.reduce$logFC.tgwt.t3))
vijay.heart.reduce$sig.t4 <- vijay.heart.reduce$Log.Bayes.Factor > 0.5 &
    vijay.heart.reduce$pp.t4 > 0.95 &
    (sign(vijay.heart.reduce$logFC.doxnd.t4) == sign(vijay.heart.reduce$logFC.tgwt.t4))
vijay.heart.reduce$sig.t5 <- vijay.heart.reduce$Log.Bayes.Factor > 0.5 &
    vijay.heart.reduce$pp.t5 > 0.95 &
    (sign(vijay.heart.reduce$logFC.doxnd.t5) == sign(vijay.heart.reduce$logFC.tgwt.t5))

vijay.heart.hom <- filter(data.frame(vijay.heart.reduce), Symbol %in% mouse.only$Symbol) %>% left_join(mouse.only)
human.hom <- filter(all.human.reduce, Symbol %in% human.only$Symbol) %>% left_join(human.only)
fds.hom <- filter(fds.df, Symbol %in% human.only$Symbol) %>% left_join(human.only)

vijay.shared.genes <- intersect(vijay.heart.hom$HomoloGene.ID, human.hom$HomoloGene.ID)
human.hom.shared <- filter(human.hom, HomoloGene.ID %in% vijay.shared.genes) %>% 
    filter(!duplicated(HomoloGene.ID)) 
vijay.heart.hom.shared <- filter(data.frame(vijay.heart.hom), HomoloGene.ID %in% vijay.shared.genes) %>%
    filter(!duplicated(HomoloGene.ID))
vijay.hearth.hyper <- map(human.groups, MapHyper, vijay.heart.groups, human.hom.shared, vijay.heart.hom.shared, "HomoloGene.ID")
names(vijay.hearth.hyper) <- str_c("human.", human.groups)

vijay.shared.genes.fds <- intersect(vijay.heart.hom$HomoloGene.ID, fds.hom$HomoloGene.ID)
fds.hom.shared <- filter(fds.hom, HomoloGene.ID %in% vijay.shared.genes.fds) %>% 
    filter(!duplicated(HomoloGene.ID))
vijay.heart.hom.shared.fds <- filter(data.frame(vijay.heart.hom), HomoloGene.ID %in% vijay.shared.genes.fds) %>% 
    filter(!duplicated(HomoloGene.ID))
vijay.heartg.hyper <- map(vijay.heart.groups, GetHyper, "FDS", vijay.heart.hom.shared.fds, fds.hom.shared, "HomoloGene.ID")
names(vijay.heartg.hyper) <- vijay.heart.groups

vijay.heart.up.bf <- map(vijay.hearth.hyper, GetBF, "Up") %>% reduce(rbind) 
vijay.heart.down.bf <- map(vijay.hearth.hyper, GetBF, "Down") %>% reduce(rbind)
vijay.heart.up.genecount <- map(vijay.hearth.hyper, GetCounts, "Up") %>% reduce(rbind)
vijay.heart.down.genecount <- map(vijay.hearth.hyper, GetCounts, "Down") %>% reduce(rbind)

vijay.heart.up.bf.fds <- map(vijay.heartg.hyper, GetBF, "Up") %>% reduce(rbind) 
rownames(vijay.heart.up.bf.fds) <- vijay.heart.groups
vijay.heart.down.bf.fds <- map(vijay.heartg.hyper, GetBF, "Down") %>% reduce(rbind) 
rownames(vijay.heart.down.bf.fds) <- vijay.heart.groups
vijay.heart.up.genecount.fds <- map(vijay.heartg.hyper, GetCounts, "Up") %>% reduce(rbind)
rownames(vijay.heart.up.genecount.fds) <- vijay.heart.groups
vijay.heart.down.genecount.fds <- map(vijay.heartg.hyper, GetCounts, "Down") %>% reduce(rbind)
rownames(vijay.heart.down.genecount.fds) <- vijay.heart.groups

#DRG
vijay.drg <- ReadRDSgz("../../Vijay_mouse/baseline/save/drg.df.final") 
vijay.drg.groups <- c("t5", "t4", "t3")
vijay.drg.reduce <- select(vijay.drg, Symbol, Log.Bayes.Factor, matches("t3|t4|t5"))
vijay.drg.reduce$sig.t3 <- vijay.drg.reduce$Log.Bayes.Factor > 0.5 &
    vijay.drg.reduce$pp.t3 > 0.95 &
    (sign(vijay.drg.reduce$logFC.doxnd.t3) == sign(vijay.drg.reduce$logFC.tgwt.t3))
vijay.drg.reduce$sig.t4 <- vijay.drg.reduce$Log.Bayes.Factor > 0.5 &
    vijay.drg.reduce$pp.t4 > 0.95 &
    (sign(vijay.drg.reduce$logFC.doxnd.t4) == sign(vijay.drg.reduce$logFC.tgwt.t4))
vijay.drg.reduce$sig.t5 <- vijay.drg.reduce$Log.Bayes.Factor > 0.5 &
    vijay.drg.reduce$pp.t5 > 0.95 &
    (sign(vijay.drg.reduce$logFC.doxnd.t5) == sign(vijay.drg.reduce$logFC.tgwt.t5))

vijay.drg.hom <- filter(data.frame(vijay.drg.reduce), Symbol %in% mouse.only$Symbol) %>% left_join(mouse.only)
vijay.shared.genes.d <- intersect(vijay.drg.hom$HomoloGene.ID, human.hom$HomoloGene.ID)
human.hom.shared.d <- filter(human.hom, HomoloGene.ID %in% vijay.shared.genes.d) %>% 
    filter(!duplicated(HomoloGene.ID)) 
vijay.drg.hom.shared <- filter(data.frame(vijay.drg.hom), HomoloGene.ID %in% vijay.shared.genes.d) %>%
    filter(!duplicated(HomoloGene.ID))
vijay.drgh.hyper <- map(human.groups, MapHyper, vijay.drg.groups, human.hom.shared.d, vijay.drg.hom.shared, "HomoloGene.ID")
names(vijay.drgh.hyper) <- str_c("human.", human.groups)

vijay.shared.genes.fds.d <- intersect(vijay.drg.hom$HomoloGene.ID, fds.hom$HomoloGene.ID)
fds.hom.shared.d <- filter(fds.hom, HomoloGene.ID %in% vijay.shared.genes.fds.d) %>% 
    filter(!duplicated(HomoloGene.ID))
vijay.drg.hom.shared.fds <- filter(data.frame(vijay.drg.hom), HomoloGene.ID %in% vijay.shared.genes.fds.d) %>% 
    filter(!duplicated(HomoloGene.ID))
vijay.drgg.hyper <- map(vijay.drg.groups, GetHyper, "FDS", vijay.drg.hom.shared.fds, fds.hom.shared.d, "HomoloGene.ID")
names(vijay.drgg.hyper) <- vijay.drg.groups

vijay.drg.up.bf <- map(vijay.drgh.hyper, GetBF, "Up") %>% reduce(rbind) 
colnames(vijay.drg.up.bf) <- str_c("drg.", vijay.drg.groups)
vijay.drg.down.bf <- map(vijay.drgh.hyper, GetBF, "Down") %>% reduce(rbind)
colnames(vijay.drg.down.bf) <- str_c("drg.", vijay.drg.groups)
vijay.drg.up.genecount <- map(vijay.drgh.hyper, GetCounts, "Up") %>% reduce(rbind)
colnames(vijay.drg.up.genecount) <- str_c("drg.", vijay.drg.groups)
vijay.drg.down.genecount <- map(vijay.drgh.hyper, GetCounts, "Down") %>% reduce(rbind)
colnames(vijay.drg.down.genecount) <- str_c("drg.", vijay.drg.groups)

vijay.drg.up.bf.fds <- map(vijay.drgg.hyper, GetBF, "Up") %>% reduce(rbind) 
rownames(vijay.drg.up.bf.fds) <- str_c("drg.", vijay.drg.groups)
vijay.drg.down.bf.fds <- map(vijay.drgg.hyper, GetBF, "Down") %>% reduce(rbind) 
rownames(vijay.drg.down.bf.fds) <- str_c("drg.", vijay.drg.groups)
vijay.drg.up.genecount.fds <- map(vijay.drgg.hyper, GetCounts, "Up") %>% reduce(rbind)
rownames(vijay.drg.up.genecount.fds) <- str_c("drg.", vijay.drg.groups)
vijay.drg.down.genecount.fds <- map(vijay.drgg.hyper, GetCounts, "Down") %>% reduce(rbind)
rownames(vijay.drg.down.genecount.fds) <- str_c("drg.", vijay.drg.groups)

#Cerebellum
vijay.cerebellum <- ReadRDSgz("../../Vijay_mouse/baseline/save/cerebellum.df.final") 
vijay.cerebellum.groups <- c("t5", "t4", "t3")
vijay.cerebellum.reduce <- select(vijay.cerebellum, Symbol, Log.Bayes.Factor, matches("t3|t4|t5"))
vijay.cerebellum.reduce$sig.t3 <- vijay.cerebellum.reduce$Log.Bayes.Factor > 0.5 &
    vijay.cerebellum.reduce$pp.t3 > 0.95 &
    (sign(vijay.cerebellum.reduce$logFC.doxnd.t3) == sign(vijay.cerebellum.reduce$logFC.tgwt.t3))
vijay.cerebellum.reduce$sig.t4 <- vijay.cerebellum.reduce$Log.Bayes.Factor > 0.5 &
    vijay.cerebellum.reduce$pp.t4 > 0.95 &
    (sign(vijay.cerebellum.reduce$logFC.doxnd.t4) == sign(vijay.cerebellum.reduce$logFC.tgwt.t4))
vijay.cerebellum.reduce$sig.t5 <- vijay.cerebellum.reduce$Log.Bayes.Factor > 0.5 &
    vijay.cerebellum.reduce$pp.t5 > 0.95 &
    (sign(vijay.cerebellum.reduce$logFC.doxnd.t5) == sign(vijay.cerebellum.reduce$logFC.tgwt.t5))

vijay.cerebellum.hom <- filter(data.frame(vijay.cerebellum.reduce), Symbol %in% mouse.only$Symbol) %>% left_join(mouse.only)
vijay.shared.genes.c <- intersect(vijay.cerebellum.hom$HomoloGene.ID, human.hom$HomoloGene.ID)
human.hom.shared.c <- filter(human.hom, HomoloGene.ID %in% vijay.shared.genes.c) %>% 
    filter(!duplicated(HomoloGene.ID)) 
vijay.cerebellum.hom.shared <- filter(data.frame(vijay.cerebellum.hom), HomoloGene.ID %in% vijay.shared.genes.c) %>%
    filter(!duplicated(HomoloGene.ID))
vijay.cerebellumh.hyper <- map(human.groups, MapHyper, vijay.cerebellum.groups, human.hom.shared.c, vijay.cerebellum.hom.shared, "HomoloGene.ID")
names(vijay.cerebellumh.hyper) <- str_c("human.", human.groups)

vijay.shared.genes.fds.c <- intersect(vijay.cerebellum.hom$HomoloGene.ID, fds.hom$HomoloGene.ID)
fds.hom.shared.c <- filter(fds.hom, HomoloGene.ID %in% vijay.shared.genes.fds.c) %>% 
    filter(!duplicated(HomoloGene.ID))
vijay.cerebellum.hom.shared.fds <- filter(data.frame(vijay.cerebellum.hom), HomoloGene.ID %in% vijay.shared.genes.fds.c) %>% 
    filter(!duplicated(HomoloGene.ID))
vijay.cerebellumg.hyper <- map(vijay.cerebellum.groups, GetHyper, "FDS", vijay.cerebellum.hom.shared.fds, fds.hom.shared.c, "HomoloGene.ID")
names(vijay.cerebellumg.hyper) <- vijay.cerebellum.groups

vijay.cerebellum.up.bf <- map(vijay.cerebellumh.hyper, GetBF, "Up") %>% reduce(rbind) 
colnames(vijay.cerebellum.up.bf) <- str_c("cerebellum.", vijay.cerebellum.groups)
vijay.cerebellum.down.bf <- map(vijay.cerebellumh.hyper, GetBF, "Down") %>% reduce(rbind)
colnames(vijay.cerebellum.down.bf) <- str_c("cerebellum.", vijay.cerebellum.groups)
vijay.cerebellum.up.genecount <- map(vijay.cerebellumh.hyper, GetCounts, "Up") %>% reduce(rbind)
colnames(vijay.cerebellum.up.genecount) <- str_c("cerebellum.", vijay.cerebellum.groups)
vijay.cerebellum.down.genecount <- map(vijay.cerebellumh.hyper, GetCounts, "Down") %>% reduce(rbind)
colnames(vijay.cerebellum.down.genecount) <- str_c("cerebellum.", vijay.cerebellum.groups)

vijay.cerebellum.up.bf.fds <- map(vijay.cerebellumg.hyper, GetBF, "Up") %>% reduce(rbind) 
rownames(vijay.cerebellum.up.bf.fds) <- str_c("cerebellum.", vijay.cerebellum.groups)
vijay.cerebellum.down.bf.fds <- map(vijay.cerebellumg.hyper, GetBF, "Down") %>% reduce(rbind) 
rownames(vijay.cerebellum.down.bf.fds) <- str_c("cerebellum.", vijay.cerebellum.groups)
vijay.cerebellum.up.genecount.fds <- map(vijay.cerebellumg.hyper, GetCounts, "Up") %>% reduce(rbind)
rownames(vijay.cerebellum.up.genecount.fds) <- str_c("cerebellum.", vijay.cerebellum.groups)
vijay.cerebellum.down.genecount.fds <- map(vijay.cerebellumg.hyper, GetCounts, "Down") %>% reduce(rbind)
rownames(vijay.cerebellum.down.genecount.fds) <- str_c("cerebellum.", vijay.cerebellum.groups)

#PBMC study
pbmc.groups <- c("cc", "pco")
pbmc.all <- ReadRDSgz("../../FRDA_overlap/pbmc/save/pbmc.df.final") %>%
    mutate(sig.pco = (Log.Bayes.Factor > 0.5 & pp.pco > 0.95),
           sig.pca = (Log.Bayes.Factor > 0.5 & pp.pca > 0.95),
           sig.cc = (Log.Bayes.Factor > 0.5 & pp.cc > 0.95))

pbmc.shared.genes <- intersect(pbmc.all$Symbol, all.human.reduce$Symbol)
human.pbmc.shared <- filter(all.human.reduce, Symbol %in% pbmc.shared.genes) 
pbmc.human.shared <- filter(pbmc.all, Symbol %in% pbmc.shared.genes) 
pbmch.hyper <- map(human.groups, MapHyper, pbmc.groups, human.pbmc.shared, pbmc.human.shared, "Symbol")
names(pbmch.hyper) <- str_c("human.", human.groups)

pbmc.shared.genes.fds <- intersect(pbmc.all$Symbol, fds.df$Symbol)
fds.pbmc.shared <- filter(fds.df, Symbol %in% pbmc.shared.genes.fds) 
pbmc.fds.shared <- filter(pbmc.all, Symbol %in% pbmc.shared.genes.fds) 
pbmcg.hyper <- map(pbmc.groups, GetHyper, "FDS", pbmc.fds.shared, fds.pbmc.shared, "Symbol")
names(pbmcg.hyper) <- str_c("pbmc.", pbmc.groups)

pbmc.up.bf <- map(pbmch.hyper, GetBF, "Up") %>% reduce(rbind) 
colnames(pbmc.up.bf) <- str_c("pbmc.", colnames(pbmc.up.bf))
pbmc.down.bf <- map(pbmch.hyper, GetBF, "Down") %>% reduce(rbind)
colnames(pbmc.down.bf) <- str_c("pbmc.", colnames(pbmc.down.bf))
pbmc.up.genecount <- map(pbmch.hyper, GetCounts, "Up") %>% reduce(rbind)
colnames(pbmc.up.genecount) <- str_c("pbmc.", colnames(pbmc.up.genecount))
pbmc.down.genecount <- map(pbmch.hyper, GetCounts, "Down") %>% reduce(rbind)
colnames(pbmc.down.genecount) <- str_c("pbmc.", colnames(pbmc.down.genecount))

pbmc.up.bf.fds <- map(pbmcg.hyper, GetBF, "Up") %>% reduce(rbind) 
rownames(pbmc.up.bf.fds) <- pbmc.groups
pbmc.down.bf.fds <- map(pbmcg.hyper, GetBF, "Down") %>% reduce(rbind)
rownames(pbmc.down.bf.fds) <- pbmc.groups
pbmc.up.genecount.fds <- map(pbmcg.hyper, GetCounts, "Up") %>% reduce(rbind)
rownames(pbmc.up.genecount.fds) <- pbmc.groups 
pbmc.down.genecount.fds <- map(pbmcg.hyper, GetCounts, "Down") %>% reduce(rbind)
rownames(pbmc.down.genecount.fds) <- pbmc.groups

#IPSC
#Read in IPSC data
#ipsc.groups <- "Healthy.vs.FRDA"
#ipsc <- ReadRDSgz("../../ipsc/ebam.ipsc.df.rda") %>% select(Z.score, Posterior) 
#ipsc$Sign.Posterior <- sign(ipsc$Z.score) * ipsc$Posterior
#colnames(ipsc) <- str_c("ipsc", colnames(ipsc), sep = ".")
#ipsc$Symbol <- rownames(ipsc)

##Compute RRHO between DE and IPSC data
#ipsch.hyper <- map(human.groups, GetHyper, colname2 = "ipsc", ipsc, all.human.reduce, "Symbol", 0.80)
#ipscg.hyper <- GetHyper("ipsc", "gaa", ipsc, ebam.gaa.df, "Symbol", 0.80)

#ipsc.up.bf <- GetBF(ipsch.hyper, "Up")
#ipsc.down.bf <- GetBF(ipsch.hyper, "Down")
#ipsc.up.genecount <- GetCounts(ipsch.hyper, "Up")
#ipsc.down.genecount <- GetCounts(ipsch.hyper, "Down")

#Van Houten
vanhouten.groups <- "pco"
vanhouten <- ReadRDSgz("../../VanHauten/baseline/save/pco.df.final")
vanhouten$sig.pco <- vanhouten$Log.Bayes.Factor > 0.5 & vanhouten$pp.pco > 0.95

vanhouten.shared.genes <- intersect(vanhouten$Symbol, all.human.reduce$Symbol)
human.vanhouten.shared <- filter(all.human.reduce, Symbol %in% vanhouten.shared.genes) 
vanhouten.human.shared <- filter(vanhouten, Symbol %in% vanhouten.shared.genes) 
vanhoutenh.hyper <- map(human.groups, GetHyper, colname2 = "pco", vanhouten.human.shared, human.vanhouten.shared, "Symbol")

vanhouten.shared.genes.fds <- intersect(vanhouten$Symbol, fds.df$Symbol)
fds.vanhouten.shared <- filter(fds.df, Symbol %in% vanhouten.shared.genes.fds) 
vanhouten.fds.shared <- filter(vanhouten, Symbol %in% vanhouten.shared.genes.fds) 
vanhouteng.hyper <- GetHyper("pco", "FDS", vanhouten.fds.shared, fds.vanhouten.shared, "Symbol")

vanhouten.up.bf <- GetBF(vanhoutenh.hyper, "Up")
vanhouten.down.bf <- GetBF(vanhoutenh.hyper, "Down")
vanhouten.up.genecount <- GetCounts(vanhoutenh.hyper, "Up")
vanhouten.down.genecount <- GetCounts(vanhoutenh.hyper, "Down")

vanhouten.up.bf.fds <- GetBF(vanhouteng.hyper, "Up")
vanhouten.down.bf.fds <- GetBF(vanhouteng.hyper, "Down")
vanhouten.up.genecount.fds <- GetCounts(vanhouteng.hyper, "Up")
vanhouten.down.genecount.fds <- GetCounts(vanhouteng.hyper, "Down")

#Combine all results
bf.up <- data.frame(vijay.heart.up.bf, pbmc.up.bf, FRDA.vs.Healthy.vanhouten = vanhouten.up.bf) %>% 
    log10 %>% signif(3) 
bf.down <- data.frame(vijay.heart.down.bf, pbmc.down.bf, FRDA.vs.Healthy.vanhouten = vanhouten.down.bf) %>% 
    log10 %>% signif(3)
bf.up.fds <- c(t(vijay.heart.up.bf.fds), pbmc.up.bf.fds, vanhouten.up.bf.fds) %>% 
    log10 %>% signif(3)
bf.down.fds <- c(t(vijay.heart.down.bf.fds), pbmc.down.bf.fds, vanhouten.down.bf.fds) %>% 
    log10 %>% signif(3)

count.up <- data.frame(vijay.heart.up.genecount, pbmc.up.genecount, FRDA.vs.Healthy.vanhouten = vanhouten.up.genecount) 
count.down <- data.frame(vijay.heart.down.genecount, pbmc.down.genecount, FRDA.vs.Healthy.vanhouten = vanhouten.down.genecount) 
count.up.fds <- c(vijay.heart.up.genecount.fds, pbmc.up.genecount.fds, vanhouten.up.genecount.fds)
count.down.fds <- c(vijay.heart.down.genecount.fds, pbmc.down.genecount.fds, vanhouten.down.genecount.fds)

bf.supp.up <- data.frame(vijay.cerebellum.up.bf, vijay.drg.up.bf) %>% 
    log10 %>% signif(3) 
bf.supp.down <- data.frame(vijay.cerebellum.down.bf, vijay.drg.down.bf) %>% 
    log10 %>% signif(3)
bf.supp.up.fds <- c(t(vijay.cerebellum.up.bf.fds), t(vijay.drg.up.bf.fds)) %>% 
    log10 %>% signif(3)
bf.supp.down.fds <- c(t(vijay.cerebellum.down.bf.fds), t(vijay.drg.down.bf.fds)) %>% 
    log10 %>% signif(3)

count.supp.up <- data.frame(vijay.cerebellum.up.genecount, vijay.drg.up.genecount) 
count.supp.down <- data.frame(vijay.cerebellum.down.genecount, vijay.drg.down.genecount) 
count.supp.up.fds <- c(vijay.cerebellum.up.genecount.fds, vijay.drg.up.genecount.fds)
count.supp.down.fds <- c(vijay.cerebellum.down.genecount.fds, vijay.drg.down.genecount.fds)

#reduced
#bf.up <- data.frame(pbmc.up.bf, FRDA.vs.Healthy.vanhouten = vanhouten.up.bf) %>% log10 %>% signif(3) %>% slice(-1)
#bf.down <- data.frame(pbmc.down.bf, FRDA.vs.Healthy.vanhouten = vanhouten.down.bf) %>% log10 %>% signif(3) %>% slice(-1)
#count.up <- data.frame(pbmc.up.genecount, FRDA.vs.Healthy.vanhouten = vanhouten.up.genecount) %>% slice(-1)
#count.down <- data.frame(pbmc.down.genecount, FRDA.vs.Healthy.vanhouten = vanhouten.down.genecount) %>% slice(-1)

up.format <- str_c(as.matrix(rbind(count.up, count.up.fds)), '\n(', 
                   as.matrix(rbind(bf.up, bf.up.fds)), ')')
dim(up.format) <- dim(rbind(bf.up, bf.up.fds))
colnames(up.format) <- colnames(bf.up)
up.format.df <- data.frame(up.format, stringsAsFactors = FALSE) 

down.format <- str_c(as.matrix(rbind(count.down, count.down.fds)), '\n(', 
                   as.matrix(rbind(bf.down, bf.down.fds)), ')')
dim(down.format) <- dim(rbind(bf.down, bf.down.fds))
colnames(down.format) <- colnames(bf.down)
down.format.df <- data.frame(down.format, stringsAsFactors = FALSE) 

up.supp.format <- str_c(as.matrix(rbind(count.supp.up, count.supp.up.fds)), '\n(', 
                   as.matrix(rbind(bf.supp.up, bf.supp.up.fds)), ')')
dim(up.supp.format) <- dim(rbind(bf.supp.up, bf.supp.up.fds))
colnames(up.supp.format) <- colnames(bf.supp.up)
up.supp.format.df <- data.frame(up.supp.format, stringsAsFactors = FALSE) 

down.supp.format <- str_c(as.matrix(rbind(count.supp.down, count.supp.down.fds)), '\n(', 
                   as.matrix(rbind(bf.supp.down, bf.supp.down.fds)), ')')
dim(down.supp.format) <- dim(rbind(bf.supp.down, bf.supp.down.fds))
colnames(down.supp.format) <- colnames(bf.supp.down)
down.supp.format.df <- data.frame(down.supp.format, stringsAsFactors = FALSE) 

up.format.df$Human.Comparison <- c("Patient vs. Carrier (Upregulated)", "Patient vs. Control (Upregulated)", "Function Score (Positive)")
up.format.df$Human.Comparison %<>% factor(levels = up.format.df$Human.Comparison)
up.plot <- gather(up.format.df, Data.Comparison, Format.BF, -Human.Comparison)
up.plot$Data.Comparison %<>% factor %>% 
    revalue(c(t5 = "RNAi mouse heart (T5)", 
              t4 = "RNAi mouse heart (T4)", 
              t3 = "RNAi mouse heart (T3)", 
              pbmc.pco = "GSE30933 (Patient vs. Control)", 
              pbmc.cc = "GSE30933 (Carrier vs. Control)", 
              FRDA.vs.Healthy.vanhouten = "GSE11204 (FRDA vs. Control)"))
up.plot$Data.Comparison %<>% factor(levels = unique(up.plot$Data.Comparison))
up.plot$Log.Bayes.Factor <- as.matrix(rbind(bf.up, bf.up.fds)) %>% as.numeric

down.format.df$Human.Comparison <- c("Patient vs. Carrier (Downregulated)", "Patient vs. Control (Downregulated)", "Function Score (Negative)")
down.format.df$Human.Comparison %<>% factor(levels = down.format.df$Human.Comparison)
down.plot <- gather(down.format.df, Data.Comparison, Format.BF, -Human.Comparison)
down.plot$Data.Comparison %<>% factor %>% 
    revalue(c(t5 = "RNAi mouse heart (T5)", 
              t4 = "RNAi mouse heart (T4)", 
              t3 = "RNAi mouse heart (T3)", 
              pbmc.pco = "GSE30933 (Patient vs. Control)", 
              pbmc.cc = "GSE30933 (Carrier vs. Control)", 
              FRDA.vs.Healthy.vanhouten = "GSE11204 (FRDA vs. Control)"))
down.plot$Data.Comparison %<>% factor(levels = unique(down.plot$Data.Comparison))
down.plot$Log.Bayes.Factor <- as.matrix(rbind(bf.down, bf.down.fds)) %>% as.numeric

combined.plot <- rbind(up.plot, down.plot)
combined.plot$Human.Comparison %<>% factor(levels = c(as.character(down.format.df$Human.Comparison), as.character(up.format.df$Human.Comparison)))

CairoPDF("up.heatmap", width = 9.0, height = 2.5, bg = "transparent")
p <- ggplot(up.plot, aes(Data.Comparison, Human.Comparison, label = Format.BF, fill = Log.Bayes.Factor)) + 
    geom_raster() + 
    geom_text() + 
    theme_bw() + 
    theme(axis.title.x = element_blank(), 
           axis.title.y = element_blank(),
           axis.ticks.x = element_blank(), 
           axis.ticks.y = element_blank(),
           axis.text.x = element_blank(),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(), 
           panel.border = element_blank(),
           plot.background = element_blank(), 
           plot.title = element_text(hjust = 0.5)) + 
    scale_fill_gradient(limits = c(0, max(up.plot$Log.Bayes.Factor)), 
                        low =  "white", high = muted('red'), guide = guide_legend(title = "logBF")) +
    ggtitle("Upregulated")
print(p)
dev.off()

CairoPDF("down.heatmap", width = 8.5, height = 4.0, bg = "transparent")
p <- ggplot(down.plot, aes(Data.Comparison, Human.Comparison, label = Format.BF, fill = Log.Bayes.Factor)) + 
    geom_raster() + 
    geom_text() + 
    theme_bw() + 
    theme(axis.title.x = element_blank(), 
           axis.title.y = element_blank(),
           axis.ticks.x = element_blank(), 
           axis.ticks.y = element_blank(),
           legend.position = "none",
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(), 
           panel.border = element_blank(),
           plot.background = element_blank(), 
           plot.title = element_text(hjust = 0.5),
           axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradient(limits = c(0, max(up.plot$Log.Bayes.Factor)), 
                        low =  "white", high = muted('red'), guide = guide_legend(title = "logBF")) +
    ggtitle("Downregulated")
print(p)
dev.off()

up.supp.format.df$Human.Comparison <- c("Patient vs. Carrier (Upregulated)", "Patient vs. Control (Upregulated)", "Function Score (Positive)")
up.supp.format.df$Human.Comparison %<>% factor(levels = up.supp.format.df$Human.Comparison)
up.supp.plot <- gather(up.supp.format.df, Data.Comparison, Format.BF, -Human.Comparison)
up.supp.plot$Data.Comparison %<>% factor %>% 
    revalue(c(cerebellum.t5 = "RNAi mouse cerebellum (T5)", 
              cerebellum.t4 = "RNAi mouse cerebellum (T4)", 
              cerebellum.t3 = "RNAi mouse cerebellum (T3)", 
              drg.t5 = "RNAi mouse drg (T5)", 
              drg.t4 = "RNAi mouse drg (T4)", 
              drg.t3 = "RNAi mouse drg (T3)"))
up.supp.plot$Data.Comparison %<>% factor(levels = unique(up.supp.plot$Data.Comparison))
up.supp.plot$Log.Bayes.Factor <- as.matrix(rbind(bf.supp.up, bf.supp.up.fds)) %>% as.numeric

down.supp.format.df$Human.Comparison <- c("Patient vs. Carrier (Downregulated)", "Patient vs. Control (Downregulated)", "Function Score (Negative)")
down.supp.format.df$Human.Comparison %<>% factor(levels = down.supp.format.df$Human.Comparison)
down.supp.plot <- gather(down.supp.format.df, Data.Comparison, Format.BF, -Human.Comparison)
down.supp.plot$Data.Comparison %<>% factor %>% 
    revalue(c(cerebellum.t5 = "RNAi mouse cerebellum (T5)", 
              cerebellum.t4 = "RNAi mouse cerebellum (T4)", 
              cerebellum.t3 = "RNAi mouse cerebellum (T3)", 
              drg.t5 = "RNAi mouse drg (T5)", 
              drg.t4 = "RNAi mouse drg (T4)", 
              drg.t3 = "RNAi mouse drg (T3)"))
down.supp.plot$Data.Comparison %<>% factor(levels = unique(down.supp.plot$Data.Comparison))
down.supp.plot$Log.Bayes.Factor <- as.matrix(rbind(bf.supp.down, bf.supp.down.fds)) %>% as.numeric

#combined.supp.plot <- rbind(up.supp.plot, down.supp.plot)
#combined.supp.plot$Human.Comparison %<>% factor(levels = c(as.character(down.supp.format.df$Human.Comparison), as.character(up.supp.format.df$Human.Comparison)))

#CairoPDF("combined.supp.heatmap", width = 10.0, height = 6, bg = "transparent")
#p <- ggplot(combined.supp.plot, aes(Data.Comparison, Human.Comparison, label = Format.BF, fill = Log.Bayes.Factor)) + 
    #geom_raster() + 
    #geom_text() + 
    #theme_bw() + 
    #theme(axis.title.x = element_blank(), 
          #axis.title.y = element_blank(),
          #axis.ticks.x = element_blank(), 
          #axis.ticks.y = element_blank(),
          #panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(),
          #panel.background = element_blank(), 
          #panel.border = element_blank(),
          #plot.background = element_blank(), 
          #plot.title = element_text(hjust = 0.5),
          #axis.text.x = element_text(angle = 45, hjust = 1)) + 
    #scale_fill_gradient(low =  "white", high = muted('red'), guide = guide_legend(title = "logBF"))
#print(p)
#dev.off()

CairoPDF("up.supp.heatmap", width = 9.0, height = 2.5, bg = "transparent")
p <- ggplot(up.supp.plot, aes(Data.Comparison, Human.Comparison, label = Format.BF, fill = Log.Bayes.Factor)) + 
    geom_raster() + 
    geom_text() + 
    theme_bw() + 
    theme(axis.title.x = element_blank(), 
           axis.title.y = element_blank(),
           axis.ticks.x = element_blank(), 
           axis.ticks.y = element_blank(),
           axis.text.x = element_blank(),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(), 
           panel.border = element_blank(),
           plot.background = element_blank(), 
           plot.title = element_text(hjust = 0.5)) + 
    scale_fill_gradient(limits = c(0, max(up.plot$Log.Bayes.Factor)), 
                        low =  "white", high = muted('red'), guide = guide_legend(title = "logBF")) +
    ggtitle("Upregulated")
print(p)
dev.off()

CairoPDF("down.supp.heatmap", width = 8.5, height = 4.0, bg = "transparent")
p <- ggplot(down.supp.plot, aes(Data.Comparison, Human.Comparison, label = Format.BF, fill = Log.Bayes.Factor)) + 
    geom_raster() + 
    geom_text() + 
    theme_bw() + 
    theme(axis.title.x = element_blank(), 
           axis.title.y = element_blank(),
           axis.ticks.x = element_blank(), 
           axis.ticks.y = element_blank(),
           legend.position = "none",
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(), 
           panel.border = element_blank(),
           plot.background = element_blank(), 
           plot.title = element_text(hjust = 0.5),
           axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradient(limits = c(0, max(up.plot$Log.Bayes.Factor)), 
                        low =  "white", high = muted('red'), guide = guide_legend(title = "logBF")) +
    ggtitle("Downregulated")
print(p)
dev.off()
