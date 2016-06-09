#Some utils
library(R.utils)

#For DE analysis
library(Biobase)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(annotate)
library(lumi)
library(WGCNA) #for fastcor
library(RRHO)

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

separate.genes <- function(group.name, dataset, operation, cutoff, symbol.vector)
{
    colname.coef <- paste("Coef", group.name, sep = ".")
    pvalue.coef <- paste("p.value", group.name, sep = ".")
    rows.de <- filter_(dataset, paste(colname.coef, operation, "0")) %>% arrange_(pvalue.coef) %>% slice(1:cutoff)
    genes.de <- rows.de[[symbol.vector]]

    return(genes.de)
}

genes.updown <- function(group.names, dataset, cutoff, symbol.vector = "HomoloGene.ID")
{
    
    genes.up <- map(group.names, separate.genes, dataset, ">=", cutoff, symbol.vector)
    genes.down <- map(group.names, separate.genes, dataset, "<", cutoff, symbol.vector)
    names(genes.up) <- group.names
    names(genes.down) <- group.names

    return(list(up = genes.up, down = genes.down))
}

compare.lists <- function(left.dataset, right.dataset)
{
    map(right.dataset, intersect, left.dataset) %>% map_int(length)
}

rrho.lists <- function(left.dataset, right.dataset)
{
    map(right.dataset, RRHO, left.dataset, alternative = "enrichment")
}

get.oddsratio <- function(overlap, cutoff, overlap.vector)
{
    or.value <- (overlap / (cutoff - overlap)) / ((cutoff - overlap) / (length(vijay.heart.overlap) - cutoff))
    return(or.value)
}

get.rrho <- function(colname1, colname2, dataset1, dataset2, symbolname, stepsize = 100)
{
    print(colname1)
    print(colname2)
    subset1 <- select_(dataset1, symbolname, str_c(colname1, ".Log.Pvalue"))
    subset2 <- select_(dataset2, symbolname, str_c(colname2, ".Log.Pvalue"))
    rrho.out <- RRHO(subset1, subset2, alternative = "enrichment", stepsize = stepsize, BY = TRUE)
    return(rrho.out)
}

dtw.rrho <- function(colname1, colname2, dataset1, dataset2, symbolname, stepsize = 50)
{
    print(colname1)
    print(colname2)
    subset1 <- select_(dataset1, symbolname, colname1)
    subset2 <- select_(dataset2, symbolname, colname2)
    rrho.out <- RRHO(subset1, subset2, alternative = "enrichment", stepsize = stepsize, BY = TRUE)
    return(rrho.out)
}

mouse.homology <- read_tsv("./HOM_MouseHumanSequence.rpt.txt") %>% data.frame %>% select(HomoloGene.ID, NCBI.Taxon.ID, Symbol)
mouse.only <- filter(mouse.homology, NCBI.Taxon.ID == 10090)
human.only <- filter(mouse.homology, NCBI.Taxon.ID == 9606)

human.groups <- c("cc", "pca", "pco")
pco.human.genes <- readRDS.gz("../baseline_lumi/save/top.object.pco.rda")
colnames(pco.human.genes) <- str_c("pco", colnames(pco.human.genes), sep = ".")
pco.human.genes$pco.Log.Pvalue <- (-log10(pco.human.genes$pco.P.Value) * sign(pco.human.genes$pco.logFC))

pca.human.genes <- readRDS.gz("../baseline_lumi/save/top.object.pca.rda")
colnames(pca.human.genes) <- str_c("pca", colnames(pca.human.genes), sep = ".")
pca.human.genes$pca.Log.Pvalue <- (-log10(pca.human.genes$pca.P.Value) * sign(pca.human.genes$pca.logFC))

cc.human.genes <- readRDS.gz("../baseline_lumi/save/top.object.cc.rda")
colnames(cc.human.genes) <- str_c("cc", colnames(cc.human.genes), sep = ".")
cc.human.genes$cc.Log.Pvalue <- (-log10(cc.human.genes$cc.P.Value) * sign(cc.human.genes$cc.logFC))

all.human.genes <- cbind(pca.human.genes, pco.human.genes, cc.human.genes)
all.human.genes$Symbol <- rownames(all.human.genes)

all.human.reduce <- select(all.human.genes, Symbol, contains("Log.Pvalue"))
human.filtered <- filter(all.human.reduce, Symbol %in% human.only$Symbol)
human.homology <- join(human.filtered, human.only)

gaa.genes <- read.xlsx("../../FRDA project/WGCNA_GAA/all.gaa.xlsx")
gaa.filtered <- filter(gaa.genes, Symbol %in% human.only$Symbol)
gaa.homology <- join(gaa.filtered, human.only)
gaa.homology$gaa.Log.Pvalue <- -log10(gaa.homology$P.value) * sign(gaa.homology$Correlation)

vijay.groups <- c("doxnd", "tgwt", "rescue")
doxnd.heart.genes <- readRDS.gz("../../Vijay_mouse/baseline/save/top.object.doxnd.rda")
colnames(doxnd.heart.genes) <- str_c("doxnd", colnames(doxnd.heart.genes), sep = ".")
doxnd.heart.genes$doxnd.Log.Pvalue <- (-log10(doxnd.heart.genes$doxnd.P.Value) * sign(doxnd.heart.genes$doxnd.logFC))

tgwt.heart.genes <- readRDS.gz("../../Vijay_mouse/baseline/save/top.object.tgwt.rda")
colnames(tgwt.heart.genes) <- str_c("tgwt", colnames(tgwt.heart.genes), sep = ".")
tgwt.heart.genes$tgwt.Log.Pvalue <- (-log10(tgwt.heart.genes$tgwt.P.Value) * sign(tgwt.heart.genes$tgwt.logFC))

rescue.heart.genes <- readRDS.gz("../../Vijay_mouse/baseline/save/top.object.rescue.rda")
colnames(rescue.heart.genes) <- str_c("rescue", colnames(rescue.heart.genes), sep = ".")
rescue.heart.genes$rescue.Log.Pvalue <- (-log10(rescue.heart.genes$rescue.P.Value) * sign(rescue.heart.genes$rescue.logFC))

vijay.heart <- cbind(doxnd.heart.genes, tgwt.heart.genes, rescue.heart.genes)
vijay.heart$Symbol <- rownames(vijay.heart)

vijay.heart.reduce <- select(vijay.heart, Symbol, contains("Log.Pvalue"))
vijay.heart.filtered <- filter(vijay.heart.reduce, Symbol %in% mouse.only$Symbol)
vijay.heart.homology <- join(vijay.heart.filtered, mouse.only)

vijay.heart.overlap <- intersect(human.homology$HomoloGene.ID, vijay.heart.homology$HomoloGene.ID)
human.vijay.heart <- filter(human.homology, HomoloGene.ID %in% vijay.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))
vijay.heart.human <- filter(vijay.heart.homology, HomoloGene.ID %in% vijay.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))

hvh.rrho <- map(vijay.groups, mkchain( map(human.groups, get.rrho, .,  human.vijay.heart, vijay.heart.human, "HomoloGene.ID"))) %>% flatten
hvh.rrho.logpval <- map(hvh.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
hvh.rrho.pval <- exp(-hvh.rrho.logpval)
dim(hvh.rrho.pval) <- c(3,3)
rownames(hvh.rrho.pval) <- human.groups
colnames(hvh.rrho.pval) <- vijay.groups

colnames(hvh.rrho.pval) %<>% str_c("heart", sep = ".")

vijay.heart.overlap.gaa <- intersect(gaa.homology$HomoloGene.ID, vijay.heart.homology$HomoloGene.ID)
gaa.vijay.heart <- filter(gaa.homology, HomoloGene.ID %in% vijay.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))
vijay.heart.gaa <- filter(vijay.heart.homology, HomoloGene.ID %in% vijay.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))

hvg.rrho <- map(vijay.groups, get.rrho, "gaa", vijay.heart.gaa, gaa.vijay.heart, "HomoloGene.ID") 
hvg.rrho.logpval <- map(hvg.rrho, getElement, "hypermat.by") %>% map_dbl(Compose(na.omit, max)) 
hvg.rrho.pval <- exp(-hvg.rrho.logpval)

doxnd.drg.genes <- readRDS.gz("../../Vijay_mouse/baseline_drg/save/top.object.doxnd.rda")
colnames(doxnd.drg.genes) <- str_c("doxnd", colnames(doxnd.drg.genes), sep = ".")
doxnd.drg.genes$doxnd.Log.Pvalue <- (-log10(doxnd.drg.genes$doxnd.P.Value) * sign(doxnd.drg.genes$doxnd.logFC))

tgwt.drg.genes <- readRDS.gz("../../Vijay_mouse/baseline_drg/save/top.object.tgwt.rda")
colnames(tgwt.drg.genes) <- str_c("tgwt", colnames(tgwt.drg.genes), sep = ".")
tgwt.drg.genes$tgwt.Log.Pvalue <- (-log10(tgwt.drg.genes$tgwt.P.Value) * sign(tgwt.drg.genes$tgwt.logFC))

rescue.drg.genes <- readRDS.gz("../../Vijay_mouse/baseline_drg/save/top.object.rescue.rda")
colnames(rescue.drg.genes) <- str_c("rescue", colnames(rescue.drg.genes), sep = ".")
rescue.drg.genes$rescue.Log.Pvalue <- (-log10(rescue.drg.genes$rescue.P.Value) * sign(rescue.drg.genes$rescue.logFC))

vijay.drg <- cbind(doxnd.drg.genes, tgwt.drg.genes, rescue.drg.genes)
vijay.drg$Symbol <- rownames(vijay.drg)

vijay.drg.reduce <- select(vijay.drg, Symbol, contains("Log.Pvalue"))
vijay.drg.filtered <- filter(vijay.drg.reduce, Symbol %in% mouse.only$Symbol)
vijay.drg.homology <- join(vijay.drg.filtered, mouse.only)

vijay.drg.overlap <- intersect(human.homology$HomoloGene.ID, vijay.drg.homology$HomoloGene.ID)
human.vijay.drg <- filter(human.homology, HomoloGene.ID %in% vijay.drg.overlap) %>% filter(!duplicated(HomoloGene.ID))
vijay.drg.human <- filter(vijay.drg.homology, HomoloGene.ID %in% vijay.drg.overlap) %>% filter(!duplicated(HomoloGene.ID))

hvd.rrho <- map(vijay.groups, mkchain( map(human.groups, get.rrho, .,  human.vijay.drg, vijay.drg.human, "HomoloGene.ID"))) %>% flatten
hvd.rrho.logpval <- map(hvd.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
hvd.rrho.pval <- exp(-hvd.rrho.logpval)
dim(hvd.rrho.pval) <- c(3,3)
rownames(hvd.rrho.pval) <- human.groups
colnames(hvd.rrho.pval) <- vijay.groups

colnames(hvd.rrho.pval) %<>% str_c("drg", sep = ".")

vijay.drg.overlap.gaa <- intersect(gaa.homology$HomoloGene.ID, vijay.drg.homology$HomoloGene.ID)
gaa.vijay.drg <- filter(gaa.homology, HomoloGene.ID %in% vijay.drg.overlap) %>% filter(!duplicated(HomoloGene.ID))
vijay.drg.gaa <- filter(vijay.drg.homology, HomoloGene.ID %in% vijay.drg.overlap) %>% filter(!duplicated(HomoloGene.ID))

gvd.rrho <- map(vijay.groups, get.rrho, "gaa", vijay.drg.gaa, gaa.vijay.drg, "HomoloGene.ID") 
gvd.rrho.logpval <- map(gvd.rrho, getElement, "hypermat.by") %>% map_dbl(Compose(na.omit, max)) 
gvd.rrho.pval <- exp(-gvd.rrho.logpval)

doxnd.cerebellum.genes <- readRDS.gz("../../Vijay_mouse/baseline_cerebellum/save/top.object.doxnd.rda")
colnames(doxnd.cerebellum.genes) <- str_c("doxnd", colnames(doxnd.cerebellum.genes), sep = ".")
doxnd.cerebellum.genes$doxnd.Log.Pvalue <- (-log10(doxnd.cerebellum.genes$doxnd.P.Value) * sign(doxnd.cerebellum.genes$doxnd.logFC))

tgwt.cerebellum.genes <- readRDS.gz("../../Vijay_mouse/baseline_cerebellum/save/top.object.tgwt.rda")
colnames(tgwt.cerebellum.genes) <- str_c("tgwt", colnames(tgwt.cerebellum.genes), sep = ".")
tgwt.cerebellum.genes$tgwt.Log.Pvalue <- (-log10(tgwt.cerebellum.genes$tgwt.P.Value) * sign(tgwt.cerebellum.genes$tgwt.logFC))

rescue.cerebellum.genes <- readRDS.gz("../../Vijay_mouse/baseline_cerebellum/save/top.object.rescue.rda")
colnames(rescue.cerebellum.genes) <- str_c("rescue", colnames(rescue.cerebellum.genes), sep = ".")
rescue.cerebellum.genes$rescue.Log.Pvalue <- (-log10(rescue.cerebellum.genes$rescue.P.Value) * sign(rescue.cerebellum.genes$rescue.logFC))

vijay.cerebellum <- cbind(doxnd.cerebellum.genes, tgwt.cerebellum.genes, rescue.cerebellum.genes)
vijay.cerebellum$Symbol <- rownames(vijay.cerebellum)

vijay.cerebellum.reduce <- select(vijay.cerebellum, Symbol, contains("Log.Pvalue"))
vijay.cerebellum.filtered <- filter(vijay.cerebellum.reduce, Symbol %in% mouse.only$Symbol)
vijay.cerebellum.homology <- join(vijay.cerebellum.filtered, mouse.only)

vijay.cerebellum.overlap <- intersect(human.homology$HomoloGene.ID, vijay.cerebellum.homology$HomoloGene.ID)
human.vijay.cerebellum <- filter(human.homology, HomoloGene.ID %in% vijay.cerebellum.overlap) %>% filter(!duplicated(HomoloGene.ID))
vijay.cerebellum.human <- filter(vijay.cerebellum.homology, HomoloGene.ID %in% vijay.cerebellum.overlap) %>% filter(!duplicated(HomoloGene.ID))

hvc.rrho <- map(vijay.groups, mkchain( map(human.groups, get.rrho, .,  human.vijay.cerebellum, vijay.cerebellum.human, "HomoloGene.ID"))) %>% flatten
hvc.rrho.logpval <- map(hvc.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
hvc.rrho.pval <- exp(-hvc.rrho.logpval)
dim(hvc.rrho.pval) <- c(3,3)
rownames(hvc.rrho.pval) <- human.groups
colnames(hvc.rrho.pval) <- vijay.groups

colnames(hvc.rrho.pval) %<>% str_c("cerebellum", sep = ".")

vijay.cerebellum.overlap.gaa <- intersect(gaa.homology$HomoloGene.ID, vijay.cerebellum.homology$HomoloGene.ID)
gaa.vijay.cerebellum <- filter(gaa.homology, HomoloGene.ID %in% vijay.cerebellum.overlap) %>% filter(!duplicated(HomoloGene.ID))
vijay.cerebellum.gaa <- filter(vijay.cerebellum.homology, HomoloGene.ID %in% vijay.cerebellum.overlap) %>% filter(!duplicated(HomoloGene.ID))

gvc.rrho <- map(vijay.groups, get.rrho, "gaa", vijay.cerebellum.gaa, gaa.vijay.cerebellum, "HomoloGene.ID") 
gvc.rrho.logpval <- map(gvc.rrho, getElement, "hypermat.by") %>% map_dbl(Compose(na.omit, max)) 
gvc.rrho.pval <- exp(-gvc.rrho.logpval)

kiko.groups <- "KIKO_vs_WT"
kiko.heart <- readRDS.gz("../../kiko/save/top.object.heart.rda")
colnames(kiko.heart) <- str_c("KIKO_vs_WT", colnames(kiko.heart), sep = ".")
kiko.heart$KIKO_vs_WT.Log.Pvalue <- (-log10(kiko.heart$KIKO_vs_WT.P.Value) * sign(kiko.heart$KIKO_vs_WT.logFC))
kiko.heart$Symbol <- rownames(kiko.heart)

kiko.heart.reduce <- select(kiko.heart, Symbol, KIKO_vs_WT.Log.Pvalue)
kiko.heart.filtered <- filter(kiko.heart.reduce, Symbol %in% mouse.only$Symbol)
kiko.heart.homology <- join(kiko.heart.filtered, mouse.only)

kiko.heart.overlap <- intersect(human.homology$HomoloGene.ID, kiko.heart.homology$HomoloGene.ID)
human.kiko.heart <- filter(human.homology, HomoloGene.ID %in% kiko.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.heart.human <- filter(kiko.heart.homology, HomoloGene.ID %in% kiko.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))

kikoh.rrho <- map(human.groups, get.rrho, "KIKO_vs_WT", human.kiko.heart, kiko.heart.human, "HomoloGene.ID")
kikoh.rrho.logpval <- map(kikoh.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
kikoh.rrho.pval <- exp(-kikoh.rrho.logpval)

kiko.heart.overlap.gaa <- intersect(gaa.homology$HomoloGene.ID, kiko.heart.homology$HomoloGene.ID)
gaa.kiko.heart <- filter(gaa.homology, HomoloGene.ID %in% kiko.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.heart.gaa <- filter(kiko.heart.homology, HomoloGene.ID %in% kiko.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))

kikohg.rrho <- get.rrho("KIKO_vs_WT", "gaa", kiko.heart.gaa, gaa.kiko.heart, "HomoloGene.ID")
kikohg.rrho.logpval <- na.omit(kikohg.rrho$hypermat.by) %>% max
kikohg.rrho.pval <- exp(-kikohg.rrho.logpval)

kiko.liver <- readRDS.gz("../../kiko/save/top.object.liver.rda")
colnames(kiko.liver) <- str_c("KIKO_vs_WT", colnames(kiko.liver), sep = ".")
kiko.liver$KIKO_vs_WT.Log.Pvalue <- (-log10(kiko.liver$KIKO_vs_WT.P.Value) * sign(kiko.liver$KIKO_vs_WT.logFC))
kiko.liver$Symbol <- rownames(kiko.liver)

kiko.liver.reduce <- select(kiko.liver, Symbol, KIKO_vs_WT.Log.Pvalue)
kiko.liver.filtered <- filter(kiko.liver.reduce, Symbol %in% mouse.only$Symbol)
kiko.liver.homology <- join(kiko.liver.filtered, mouse.only)

kiko.liver.overlap <- intersect(human.homology$HomoloGene.ID, kiko.liver.homology$HomoloGene.ID)
human.kiko.liver <- filter(human.homology, HomoloGene.ID %in% kiko.liver.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.liver.human <- filter(kiko.liver.homology, HomoloGene.ID %in% kiko.liver.overlap) %>% filter(!duplicated(HomoloGene.ID))

kikol.rrho <- map(human.groups, get.rrho, "KIKO_vs_WT", human.kiko.liver, kiko.liver.human, "HomoloGene.ID")
kikol.rrho.logpval <- map(kikol.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
kikol.rrho.pval <- exp(-kikol.rrho.logpval)

kiko.liver.overlap.gaa <- intersect(gaa.homology$HomoloGene.ID, kiko.liver.homology$HomoloGene.ID)
gaa.kiko.liver <- filter(gaa.homology, HomoloGene.ID %in% kiko.liver.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.liver.gaa <- filter(kiko.liver.homology, HomoloGene.ID %in% kiko.liver.overlap) %>% filter(!duplicated(HomoloGene.ID))

kikolg.rrho <- get.rrho("KIKO_vs_WT", "gaa", kiko.liver.gaa, gaa.kiko.liver, "HomoloGene.ID")
kikolg.rrho.logpval <- na.omit(kikolg.rrho$hypermat.by) %>% max
kikolg.rrho.pval <- exp(-kikolg.rrho.logpval)

kiko.muscle <- readRDS.gz("../../kiko/save/top.object.muscle.rda")
colnames(kiko.muscle) <- str_c("KIKO_vs_WT", colnames(kiko.muscle), sep = ".")
kiko.muscle$KIKO_vs_WT.Log.Pvalue <- (-log10(kiko.muscle$KIKO_vs_WT.P.Value) * sign(kiko.muscle$KIKO_vs_WT.logFC))
kiko.muscle$Symbol <- rownames(kiko.muscle)

kiko.muscle.reduce <- select(kiko.muscle, Symbol, KIKO_vs_WT.Log.Pvalue)
kiko.muscle.filtered <- filter(kiko.muscle.reduce, Symbol %in% mouse.only$Symbol)
kiko.muscle.homology <- join(kiko.muscle.filtered, mouse.only)

kiko.muscle.overlap <- intersect(human.homology$HomoloGene.ID, kiko.muscle.homology$HomoloGene.ID)
human.kiko.muscle <- filter(human.homology, HomoloGene.ID %in% kiko.muscle.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.muscle.human <- filter(kiko.muscle.homology, HomoloGene.ID %in% kiko.muscle.overlap) %>% filter(!duplicated(HomoloGene.ID))

kikom.rrho <- map(human.groups, get.rrho, "KIKO_vs_WT", human.kiko.muscle, kiko.muscle.human, "HomoloGene.ID")
kikom.rrho.logpval <- map(kikom.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
kikom.rrho.pval <- exp(-kikom.rrho.logpval)

kiko.muscle.overlap.gaa <- intersect(gaa.homology$HomoloGene.ID, kiko.muscle.homology$HomoloGene.ID)
gaa.kiko.muscle <- filter(gaa.homology, HomoloGene.ID %in% kiko.muscle.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.muscle.gaa <- filter(kiko.muscle.homology, HomoloGene.ID %in% kiko.muscle.overlap) %>% filter(!duplicated(HomoloGene.ID))

kikomg.rrho <- get.rrho("KIKO_vs_WT", "gaa", kiko.muscle.gaa, gaa.kiko.muscle, "HomoloGene.ID")
kikomg.rrho.logpval <- na.omit(kikomg.rrho$hypermat.by) %>% max
kikomg.rrho.pval <- exp(-kikomg.rrho.logpval)

kiko.rrho.pval <- data.frame(KIKO_vs_WT.heart = kikoh.rrho.pval, KIKO_vs_WT.muscle = kikom.rrho.pval, KIKO_vs_WT.liver = kikol.rrho.pval)
kikog.rrho.pval <- c(kikohg.rrho.pval, kikomg.rrho.pval, kikolg.rrho.pval)

#kiki!

pbmc.groups <- str_c(c("pca", "pco", "cc"), "pbmc", sep = ".")
pco.pbmc.genes <- readRDS.gz("../../pbmc/save/top.object.pco.rda")
colnames(pco.pbmc.genes) <- str_c("pco.pbmc", colnames(pco.pbmc.genes), sep = ".")
pco.pbmc.genes$pco.pbmc.Log.Pvalue <- (-log10(pco.pbmc.genes$pco.pbmc.P.Value) * sign(pco.pbmc.genes$pco.pbmc.logFC))

pca.pbmc.genes <- readRDS.gz("../../pbmc/save/top.object.pca.rda")
colnames(pca.pbmc.genes) <- str_c("pca.pbmc", colnames(pca.pbmc.genes), sep = ".")
pca.pbmc.genes$pca.pbmc.Log.Pvalue <- (-log10(pca.pbmc.genes$pca.pbmc.P.Value) * sign(pca.pbmc.genes$pca.pbmc.logFC))

cc.pbmc.genes <- readRDS.gz("../../pbmc/save/top.object.cc.rda")
colnames(cc.pbmc.genes) <- str_c("cc.pbmc", colnames(cc.pbmc.genes), sep = ".")
cc.pbmc.genes$cc.pbmc.Log.Pvalue <- (-log10(cc.pbmc.genes$cc.pbmc.P.Value) * sign(cc.pbmc.genes$cc.pbmc.logFC))

all.pbmc.genes <- cbind(pca.pbmc.genes, pco.pbmc.genes, cc.pbmc.genes)
all.pbmc.genes$Symbol <- rownames(all.pbmc.genes)

pbmc.reduce <- select(all.pbmc.genes, Symbol, contains("Log.Pvalue"))

pbmc.overlap <- intersect(all.human.reduce$Symbol, pbmc.reduce$Symbol)
human.pbmc <- filter(all.human.reduce, Symbol %in% pbmc.overlap)
pbmc.human <- filter(pbmc.reduce, Symbol %in% pbmc.overlap)

pbmch.rrho <- map(pbmc.groups, mkchain( map(human.groups, get.rrho, .,  human.pbmc, pbmc.human, "Symbol", 100))) %>% flatten
pbmch.rrho.logpval <- map(pbmch.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
pbmch.rrho.pval <- exp(-pbmch.rrho.logpval)
dim(pbmch.rrho.pval) <- c(3,3)
rownames(pbmch.rrho.pval) <- human.groups
colnames(pbmch.rrho.pval) <- pbmc.groups

pbmcg.overlap <- intersect(gaa.genes$Symbol, pbmc.reduce$Symbol)
gaa.pbmc <- filter(gaa.genes, Symbol %in% pbmcg.overlap)
gaa.pbmc$gaa.Log.Pvalue <- -(log10(gaa.pbmc$P.value)) * sign(gaa.pbmc$Correlation)
pbmc.gaa <- filter(pbmc.reduce, Symbol %in% pbmcg.overlap)

pbmcg.rrho <- map(pbmc.groups, get.rrho, "gaa", pbmc.gaa, gaa.pbmc, "Symbol", 100)
pbmcg.rrho.logpval <- map(pbmcg.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
pbmcg.rrho.pval <- exp(-pbmcg.rrho.logpval)

#IPSC
ipsc.groups <- "Healthy.vs.FRDA"
ipsc <- readRDS.gz("../../ipsc/save/top.object.frdah")
ipsc.reduce <- select(ipsc, logFC, P.Value)
ipsc.reduce$Log.Pvalue <- -(log10(ipsc.reduce$P.Value)) * sign(ipsc.reduce$logFC)
colnames(ipsc.reduce) <- str_c("ipsc", colnames(ipsc.reduce), sep = ".")
ipsc.reduce$Symbol <- rownames(ipsc.reduce)

ipsc.overlap <- intersect(all.human.reduce$Symbol, ipsc.reduce$Symbol)
human.ipsc <- filter(all.human.reduce, Symbol %in% ipsc.overlap)
ipsc.human <- filter(ipsc.reduce, Symbol %in% ipsc.overlap)

ipsch.rrho <- map(human.groups, get.rrho, "ipsc", human.ipsc, ipsc.human, "Symbol", 100)
ipsch.rrho.logpval <- map(ipsch.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
ipsch.rrho.pval <- exp(-ipsch.rrho.logpval)

ipsc.overlap.gaa <- intersect(gaa.genes$Symbol, ipsc.reduce$Symbol)
gaa.ipsc <- filter(gaa.genes, Symbol %in% ipsc.overlap)
gaa.ipsc$gaa.Log.Pvalue <- -(log10(gaa.ipsc$P.value)) * sign(gaa.ipsc$Correlation)
ipsc.gaa <- filter(ipsc.reduce, Symbol %in% ipsc.overlap)

ipscg.rrho <- get.rrho("gaa", "ipsc", gaa.ipsc, ipsc.gaa, "Symbol", 100)
ipscg.rrho.logpval <- na.omit(ipscg.rrho$hypermat.by) %>% max
ipscg.rrho.pval <- exp(-ipscg.rrho.logpval)

#Van Houten
vanhouten.groups <- "pco"
vanhouten <- readRDS.gz("../../VanHauten/baseline/save/top.object.pco.rda")
vanhouten.reduce <- select(vanhouten, logFC, P.Value)
vanhouten.reduce$Log.Pvalue <- -(log10(vanhouten.reduce$P.Value)) * sign(vanhouten.reduce$logFC)
colnames(vanhouten.reduce) <- str_c("vanhouten", colnames(vanhouten.reduce), sep = ".")
vanhouten.reduce$Symbol <- vanhouten$Symbol
vanhouten.reduce %<>% filter(!duplicated(Symbol))

vanhouten.overlap <- intersect(all.human.reduce$Symbol, vanhouten.reduce$Symbol)
human.vanhouten <- filter(all.human.reduce, Symbol %in% vanhouten.overlap)
vanhouten.human <- filter(vanhouten.reduce, Symbol %in% vanhouten.overlap)

vanhoutenh.rrho <- map(human.groups, get.rrho, "vanhouten", human.vanhouten, vanhouten.human, "Symbol")
vanhoutenh.rrho.logpval <- map(vanhoutenh.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
vanhoutenh.rrho.pval <- exp(-vanhoutenh.rrho.logpval)

vanhouten.overlap.gaa <- intersect(gaa.genes$Symbol, vanhouten.reduce$Symbol)
gaa.vanhouten <- filter(gaa.genes, Symbol %in% vanhouten.overlap)
colnames(gaa.vanhouten)[1:3] <- str_c("gaa", colnames(gaa.vanhouten)[1:3], sep = ".")
gaa.vanhouten$gaa.Log.Pvalue <- -(log10(gaa.vanhouten$gaa.P.value)) * sign(gaa.vanhouten$gaa.Correlation)
vanhouten.gaa <- filter(vanhouten.reduce, Symbol %in% vanhouten.overlap)

vanhouteng.rrho <- get.rrho("gaa", "vanhouten", gaa.vanhouten, vanhouten.gaa, "Symbol")
vanhouteng.rrho.logpval <- na.omit(vanhouteng.rrho$hypermat.by) %>% max
vanhouteng.rrho.pval <- exp(-vanhouteng.rrho.logpval)

pval.rrho <- data.frame(hvh.rrho.pval, hvd.rrho.pval, hvc.rrho.pval, kiko.rrho.pval, pbmch.rrho.pval, FRDA.vs.Healthy.ipsc = ipsch.rrho.pval, FRDA.vs.Healthy.vanhouten = vanhoutenh.rrho.pval)

gaa.pval.rrho <- c(hvg.rrho.pval, gvd.rrho.pval, gvc.rrho.pval, kikog.rrho.pval, pbmch.rrho.pval, ipscg.rrho.pval, vanhouteng.rrho.pval)
pval.rrho.all <- rbind(pval.rrho, gaa = gaa.pval.rrho)

pval.adjust.rrho <- map(pval.rrho.all, p.adjust, method = "fdr", n = (nrow(pval.rrho.all) * ncol(pval.rrho.all))) %>% reduce(cbind)
colnames(pval.adjust.rrho) <- colnames(pval.rrho.all)
rownames(pval.adjust.rrho) <- rownames(pval.rrho.all)

CairoPDF("overlap.heatmap", width = 20, height = 6)
par(mar = c(10, 8, 3, 3))
labeledHeatmap(Matrix = -log10(pval.adjust.rrho), xLabels = colnames(pval.adjust.rrho), invertColors = TRUE, yLabels = rownames(pval.adjust.rrho), ySymbols = y.names, textMatrix = signif(pval.adjust.rrho, 3), setStdMargins = F, cex.text = 1.3, cex.lab = 1.3, zlim = c(0,12.5))
dev.off()

#DTW
dtwh <- readRDS.gz("../dtw/save/dtw.pca.rda")
dtwh.filtered <- filter(dtwh, Symbol %in% human.only$Symbol)
dtwh.homology <- join(dtwh.filtered, human.only)

dtwm <- readRDS.gz("../../Vijay_mouse/dtw/save/dtw.doxnd.rda")
dtwm.filtered <- filter(dtwm, Symbol %in% mouse.only$Symbol)
dtwm.homology <- join(dtwm.filtered, mouse.only)

dtw.overlap <- intersect(dtwh.homology$HomoloGene.ID, dtwm.homology$HomoloGene.ID)
dtwh.vijay <- filter(dtwh.homology, HomoloGene.ID %in% dtw.overlap) %>% filter(!duplicated(HomoloGene.ID))
dtwm.human <- filter(dtwm.homology, HomoloGene.ID %in% dtw.overlap) %>% filter(!duplicated(HomoloGene.ID))

dtwh.groups <- colnames(dtwh.vijay)[2:4]
dtwm.groups <- colnames(dtwm.human)[1:3]

dtw.rrho <- map(dtwm.groups, mkchain( map(dtwh.groups, dtw.rrho, ., dtwh.vijay, dtwm.human, "HomoloGene.ID"))) %>% flatten
dtw.rrho.logpval <- map(dtw.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
dtw.rrho.pval <- exp(-dtw.rrho.logpval)
dim(dtw.rrho.pval) <- c(3,3)
rownames(dtw.rrho.pval) <- dtwh.groups
colnames(dtw.rrho.pval) <- dtwm.groups

#test.subset <- select(dtwh.vijay, HomoloGene.ID, Patient.vs.Control)
#test.subset <- select(dtwm.human, HomoloGene.ID, Tg.DOX.vs.Tg.ND)
#test.subset2 <- select(dtwm.human, HomoloGene.ID, Tg.ND.vs.WT.DOX)

#rrho.test <- RRHO(test.subset, test.subset2, alternative = "enrichment", stepsize = 50, BY = TRUE)
