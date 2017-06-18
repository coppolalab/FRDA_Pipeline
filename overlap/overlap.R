#For DE analysis
library(lumi)
library(WGCNA) #for fastcor
#library(RRHO)
library(BayesFactor)

library(Cairo)
library(openxlsx)

library(magrittr)
library(stringr)
library(tidyverse)

source("../common_functions.R")

GetRRHO <- function(colname1, colname2, dataset1, dataset2, symbolname, stepsize = 50) {
    subset1 <- select_(dataset1, symbolname, str_c(colname1, ".Log.Pvalue"))
    subset2 <- select_(dataset2, symbolname, str_c(colname2, ".Log.Pvalue"))
    rrho.dir <- str_c("./rrho/", colname1, "_", colname2)
    dir.create(rrho.dir, showWarnings = FALSE)
    rrho.out <- RRHO(subset1, subset2, alternative = "enrichment", stepsize = stepsize, plots = TRUE, labels = c(colname1, colname2), outputdir = rrho.dir, BY = TRUE)
    return(rrho.out)
}

#Read in HomoloGene table
mouse.homology <- read_tsv("./HOM_MouseHumanSequence.rpt.txt") %>% data.frame %>% select(HomoloGene.ID, NCBI.Taxon.ID, Symbol)
mouse.only <- filter(mouse.homology, NCBI.Taxon.ID == 10090)
human.only <- filter(mouse.homology, NCBI.Taxon.ID == 9606)

#Read in differential expression data
human.groups <- c("pca", "pco")
#pco.human.genes <- ReadRDSgz("../baseline_lumi/save/toptable.pco.rda")
#colnames(pco.human.genes) <- str_c("pco", colnames(pco.human.genes), sep = ".")
#pco.human.genes$pco.Log.Pvalue <- (-log10(pco.human.genes$pco.P.Value) * sign(pco.human.genes$pco.logFC))

#pca.human.genes <- ReadRDSgz("../baseline_lumi/save/toptable.pca.rda")
#colnames(pca.human.genes) <- str_c("pca", colnames(pca.human.genes), sep = ".")
#pca.human.genes$pca.Log.Pvalue <- (-log10(pca.human.genes$pca.P.Value) * sign(pca.human.genes$pca.logFC))

#cc.human.genes <- ReadRDSgz("../baseline_lumi/save/toptable.cc.rda")
#colnames(cc.human.genes) <- str_c("cc", colnames(cc.human.genes), sep = ".")
#cc.human.genes$cc.Log.Pvalue <- (-log10(cc.human.genes$cc.P.Value) * sign(cc.human.genes$cc.logFC))

#Read in EBAM
ebam.pca.df <- ReadRDSgz("../baseline_lumi/ebam.pca.df.rda") %>% select(Z.score, Posterior)
ebam.pca.df$Sign.Posterior <- sign(ebam.pca.df$Z.score) * ebam.pca.df$Posterior
colnames(ebam.pca.df) %<>% str_c(".pca")
ebam.pco.df <- ReadRDSgz("../baseline_lumi/ebam.pco.df.rda") %>% select(Z.score, Posterior)
ebam.pco.df$Sign.Posterior <- sign(ebam.pco.df$Z.score) * ebam.pco.df$Posterior
colnames(ebam.pco.df) %<>% str_c(".pco")

all.human.genes <- cbind(ebam.pca.df, ebam.pco.df)
all.human.genes$Symbol <- rownames(all.human.genes)

all.human.reduce <- select(all.human.genes, Symbol, dplyr::contains("Sign.Posterior"))
#human.filtered <- filter(all.human.reduce, Symbol %in% human.only$Symbol)
#human.homology <- left_join(human.filtered, human.only)

#Read in GAA expansion data
gaa.genes <- read.xlsx("../../FRDA project/WGCNA_GAA/gaa.cor.xlsx")
gaa.filtered <- filter(gaa.genes, Symbol %in% human.only$Symbol)
gaa.homology <- left_join(gaa.filtered, human.only)
gaa.homology$gaa.Log.Pvalue <- -log10(gaa.homology$P.value) * sign(gaa.homology$Correlation)

#KIKO
#Read in heart data
kiko.groups <- "KIKO_vs_WT"
kiko.heart <- ReadRDSgz("../../kiko/save/top.object.heart.rda")
colnames(kiko.heart) <- str_c("KIKO_vs_WT", colnames(kiko.heart), sep = ".")
kiko.heart$KIKO_vs_WT.Log.Pvalue <- (-log10(kiko.heart$KIKO_vs_WT.P.Value) * sign(kiko.heart$KIKO_vs_WT.logFC))
kiko.heart$Symbol <- rownames(kiko.heart)

#Only keep genes in HomoloGene
kiko.heart.reduce <- select(kiko.heart, Symbol, KIKO_vs_WT.Log.Pvalue)
kiko.heart.filtered <- filter(kiko.heart.reduce, Symbol %in% mouse.only$Symbol)
kiko.heart.homology <- left_join(kiko.heart.filtered, mouse.only)

#Only keep genes shared between DE and KIKO heart data
kiko.heart.overlap <- intersect(human.homology$HomoloGene.ID, kiko.heart.homology$HomoloGene.ID)
human.kiko.heart <- filter(human.homology, HomoloGene.ID %in% kiko.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.heart.human <- filter(kiko.heart.homology, HomoloGene.ID %in% kiko.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))

#Compute RRHO between DE and KIKO heart data 
kikoh.rrho <- map(human.groups, GetRRHO, "KIKO_vs_WT", human.kiko.heart, kiko.heart.human, "HomoloGene.ID")
kikoh.rrho.logpval <- map(kikoh.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
kikoh.rrho.pval <- exp(-kikoh.rrho.logpval)

#Only keep genes shared between GAA and KIKO heart data
kiko.heart.overlap.gaa <- intersect(gaa.homology$HomoloGene.ID, kiko.heart.homology$HomoloGene.ID)
gaa.kiko.heart <- filter(gaa.homology, HomoloGene.ID %in% kiko.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.heart.gaa <- filter(kiko.heart.homology, HomoloGene.ID %in% kiko.heart.overlap) %>% filter(!duplicated(HomoloGene.ID))

#Compute RRHO between DE and KIKO heart data 
kikohg.rrho <- GetRRHO("KIKO_vs_WT", "gaa", kiko.heart.gaa, gaa.kiko.heart, "HomoloGene.ID")
kikohg.rrho.logpval <- na.omit(kikohg.rrho$hypermat.by) %>% max
kikohg.rrho.pval <- exp(-kikohg.rrho.logpval)

#Read in liver data
kiko.liver <- ReadRDSgz("../../kiko/save/top.object.liver.rda")
colnames(kiko.liver) <- str_c("KIKO_vs_WT", colnames(kiko.liver), sep = ".")
kiko.liver$KIKO_vs_WT.Log.Pvalue <- (-log10(kiko.liver$KIKO_vs_WT.P.Value) * sign(kiko.liver$KIKO_vs_WT.logFC))
kiko.liver$Symbol <- rownames(kiko.liver)

#Only keep genes in HomoloGene
kiko.liver.reduce <- select(kiko.liver, Symbol, KIKO_vs_WT.Log.Pvalue)
kiko.liver.filtered <- filter(kiko.liver.reduce, Symbol %in% mouse.only$Symbol)
kiko.liver.homology <- left_join(kiko.liver.filtered, mouse.only)

#Only keep genes shared between DE and KIKO liver data
kiko.liver.overlap <- intersect(human.homology$HomoloGene.ID, kiko.liver.homology$HomoloGene.ID)
human.kiko.liver <- filter(human.homology, HomoloGene.ID %in% kiko.liver.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.liver.human <- filter(kiko.liver.homology, HomoloGene.ID %in% kiko.liver.overlap) %>% filter(!duplicated(HomoloGene.ID))

#Compute RRHO between DE and KIKO liver data 
kikol.rrho <- map(human.groups, GetRRHO, "KIKO_vs_WT", human.kiko.liver, kiko.liver.human, "HomoloGene.ID")
kikol.rrho.logpval <- map(kikol.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
kikol.rrho.pval <- exp(-kikol.rrho.logpval)

#Only keep genes shared between GAA and KIKO liver data
kiko.liver.overlap.gaa <- intersect(gaa.homology$HomoloGene.ID, kiko.liver.homology$HomoloGene.ID)
gaa.kiko.liver <- filter(gaa.homology, HomoloGene.ID %in% kiko.liver.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.liver.gaa <- filter(kiko.liver.homology, HomoloGene.ID %in% kiko.liver.overlap) %>% filter(!duplicated(HomoloGene.ID))

#Compute RRHO between DE and KIKO liver data 
kikolg.rrho <- GetRRHO("KIKO_vs_WT", "gaa", kiko.liver.gaa, gaa.kiko.liver, "HomoloGene.ID")
kikolg.rrho.logpval <- na.omit(kikolg.rrho$hypermat.by) %>% max
kikolg.rrho.pval <- exp(-kikolg.rrho.logpval)

#Read in muscle data
kiko.muscle <- ReadRDSgz("../../kiko/save/top.object.muscle.rda")
colnames(kiko.muscle) <- str_c("KIKO_vs_WT", colnames(kiko.muscle), sep = ".")
kiko.muscle$KIKO_vs_WT.Log.Pvalue <- (-log10(kiko.muscle$KIKO_vs_WT.P.Value) * sign(kiko.muscle$KIKO_vs_WT.logFC))
kiko.muscle$Symbol <- rownames(kiko.muscle)

#Only keep genes in HomoloGene
kiko.muscle.reduce <- select(kiko.muscle, Symbol, KIKO_vs_WT.Log.Pvalue)
kiko.muscle.filtered <- filter(kiko.muscle.reduce, Symbol %in% mouse.only$Symbol)
kiko.muscle.homology <- left_join(kiko.muscle.filtered, mouse.only)

#Only keep genes shared between DE and KIKO muscle data
kiko.muscle.overlap <- intersect(human.homology$HomoloGene.ID, kiko.muscle.homology$HomoloGene.ID)
human.kiko.muscle <- filter(human.homology, HomoloGene.ID %in% kiko.muscle.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.muscle.human <- filter(kiko.muscle.homology, HomoloGene.ID %in% kiko.muscle.overlap) %>% filter(!duplicated(HomoloGene.ID))

#Compute RRHO between DE and KIKO muscle data 
kikom.rrho <- map(human.groups, GetRRHO, "KIKO_vs_WT", human.kiko.muscle, kiko.muscle.human, "HomoloGene.ID")
kikom.rrho.logpval <- map(kikom.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
kikom.rrho.pval <- exp(-kikom.rrho.logpval)

#Only keep genes shared between GAA and KIKO muscle data
kiko.muscle.overlap.gaa <- intersect(gaa.homology$HomoloGene.ID, kiko.muscle.homology$HomoloGene.ID)
gaa.kiko.muscle <- filter(gaa.homology, HomoloGene.ID %in% kiko.muscle.overlap) %>% filter(!duplicated(HomoloGene.ID))
kiko.muscle.gaa <- filter(kiko.muscle.homology, HomoloGene.ID %in% kiko.muscle.overlap) %>% filter(!duplicated(HomoloGene.ID))

#Compute RRHO between DE and KIKO muscle data 
kikomg.rrho <- GetRRHO("KIKO_vs_WT", "gaa", kiko.muscle.gaa, gaa.kiko.muscle, "HomoloGene.ID")
kikomg.rrho.logpval <- na.omit(kikomg.rrho$hypermat.by) %>% max
kikomg.rrho.pval <- exp(-kikomg.rrho.logpval)

kiko.rrho.pval <- data.frame(KIKO_vs_WT.heart = kikoh.rrho.pval, KIKO_vs_WT.muscle = kikom.rrho.pval, KIKO_vs_WT.liver = kikol.rrho.pval)
kikog.rrho.pval <- c(kikohg.rrho.pval, kikomg.rrho.pval, kikolg.rrho.pval)

#kiki!

#PBMC study
model.pbmc <- ReadRDSgz("../../pbmc/save/model.pbmc.rda")
pbmc.collapse <- ReadRDSgz("../../pbmc/save/pbmc.collapse.rda")

#Only keep genes shared between DE and PBMC data
pbmc.overlap <- intersect(all.human.reduce$Symbol, pbmc.collapse)
human.pbmc <- filter(all.human.reduce, Symbol %in% pbmc.overlap)
pbmc.human <- filter(pbmc.reduce, Symbol %in% pbmc.overlap)

#Read in PBMC data
pbmc.groups <- str_c(c("pca", "pco", "cc"), "pbmc", sep = ".")
pco.pbmc.genes <- ReadRDSgz("../../pbmc/save/top.object.pco.rda")
colnames(pco.pbmc.genes) <- str_c("pco.pbmc", colnames(pco.pbmc.genes), sep = ".")
pco.pbmc.genes$pco.pbmc.Log.Pvalue <- (-log10(pco.pbmc.genes$pco.pbmc.P.Value) * sign(pco.pbmc.genes$pco.pbmc.logFC))

pca.pbmc.genes <- ReadRDSgz("../../pbmc/save/top.object.pca.rda")
colnames(pca.pbmc.genes) <- str_c("pca.pbmc", colnames(pca.pbmc.genes), sep = ".")
pca.pbmc.genes$pca.pbmc.Log.Pvalue <- (-log10(pca.pbmc.genes$pca.pbmc.P.Value) * sign(pca.pbmc.genes$pca.pbmc.logFC))

cc.pbmc.genes <- ReadRDSgz("../../pbmc/save/top.object.cc.rda")
colnames(cc.pbmc.genes) <- str_c("cc.pbmc", colnames(cc.pbmc.genes), sep = ".")
cc.pbmc.genes$cc.pbmc.Log.Pvalue <- (-log10(cc.pbmc.genes$cc.pbmc.P.Value) * sign(cc.pbmc.genes$cc.pbmc.logFC))

all.pbmc.genes <- cbind(pca.pbmc.genes, pco.pbmc.genes, cc.pbmc.genes)
all.pbmc.genes$Symbol <- rownames(all.pbmc.genes)

pbmc.reduce <- select(all.pbmc.genes, Symbol, dplyr::contains("Log.Pvalue"))

#Compute RRHO between DE and PBMC data
pbmch.rrho <- map(pbmc.groups, mkchain( map(human.groups, GetRRHO, .,  human.pbmc, pbmc.human, "Symbol"))) %>% flatten
pbmch.rrho.logpval <- map(pbmch.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
pbmch.rrho.pval <- exp(-pbmch.rrho.logpval)
dim(pbmch.rrho.pval) <- c(3,3)
rownames(pbmch.rrho.pval) <- human.groups
colnames(pbmch.rrho.pval) <- pbmc.groups

#Only keep genes shared between GAA and PBMC data
pbmcg.overlap <- intersect(gaa.genes$Symbol, pbmc.reduce$Symbol)
gaa.pbmc <- filter(gaa.genes, Symbol %in% pbmcg.overlap)
gaa.pbmc$gaa.Log.Pvalue <- -(log10(gaa.pbmc$P.value)) * sign(gaa.pbmc$Correlation)
pbmc.gaa <- filter(pbmc.reduce, Symbol %in% pbmcg.overlap)

#Compute RRHO between GAA and PBMC data
pbmcg.rrho <- map(pbmc.groups, GetRRHO, "gaa", pbmc.gaa, gaa.pbmc, "Symbol")
pbmcg.rrho.logpval <- map(pbmcg.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
pbmcg.rrho.pval <- exp(-pbmcg.rrho.logpval)

#IPSC
#Read in IPSC data
ipsc.groups <- "Healthy.vs.FRDA"
ipsc <- ReadRDSgz("../../ipsc/save/top.object.frdah")
ipsc.reduce <- select(ipsc, logFC, P.Value)
ipsc.reduce$Log.Pvalue <- -(log10(ipsc.reduce$P.Value)) * sign(ipsc.reduce$logFC)
colnames(ipsc.reduce) <- str_c("ipsc", colnames(ipsc.reduce), sep = ".")
ipsc.reduce$Symbol <- rownames(ipsc.reduce)

#Only keep genes shared between DE and IPSC data
ipsc.overlap <- intersect(all.human.reduce$Symbol, ipsc.reduce$Symbol)
human.ipsc <- filter(all.human.reduce, Symbol %in% ipsc.overlap)
ipsc.human <- filter(ipsc.reduce, Symbol %in% ipsc.overlap)

#Compute RRHO between DE and IPSC data
ipsch.rrho <- map(human.groups, GetRRHO, "ipsc", human.ipsc, ipsc.human, "Symbol")
ipsch.rrho.logpval <- map(ipsch.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
ipsch.rrho.pval <- exp(-ipsch.rrho.logpval)

#Only keep genes shared between GAA and IPSC data
ipsc.overlap.gaa <- intersect(gaa.genes$Symbol, ipsc.reduce$Symbol)
gaa.ipsc <- filter(gaa.genes, Symbol %in% ipsc.overlap)
gaa.ipsc$gaa.Log.Pvalue <- -(log10(gaa.ipsc$P.value)) * sign(gaa.ipsc$Correlation)
ipsc.gaa <- filter(ipsc.reduce, Symbol %in% ipsc.overlap)

#Compute RRHO between GAA and IPSC data
ipscg.rrho <- GetRRHO("gaa", "ipsc", gaa.ipsc, ipsc.gaa, "Symbol")
subset1 <- select(gaa.ipsc, Symbol, gaa.Log.Pvalue)
subset2 <- select(ipsc.gaa, Symbol, ipsc.Log.Pvalue)
#Must have outputdir, labels and plots = TRUE
rrho.out <- RRHO(subset1, subset2, alternative = "enrichment", stepsize = 50, plots = TRUE, labels = c("gaa", "ipsc"), outputdir = "./gaa_ipsc", BY = TRUE)
ipscg.rrho.logpval <- na.omit(ipscg.rrho$hypermat.by) %>% max
ipscg.rrho.pval <- exp(-ipscg.rrho.logpval)

gaa.ipsc.sort <- arrange(gaa.ipsc, desc(gaa.Log.Pvalue))
gaa.ipsc.sort$Index <- 1:nrow(gaa.ipsc.sort)
write.xlsx(gaa.ipsc.sort, "./gaa.ipsc.xlsx")
ipsc.gaa.sort <- arrange(ipsc.gaa, desc(ipsc.Log.Pvalue))
ipsc.gaa.sort$Index <- 1:nrow(ipsc.gaa.sort)
write.xlsx(ipsc.gaa.sort, "./ipsc.gaa.xlsx")

#Van Houten
#Read in Van Houten data
vanhouten.groups <- "pco"
vanhouten <- ReadRDSgz("../../VanHauten/baseline/save/top.object.pco.rda")
vanhouten.reduce <- select(vanhouten, logFC, P.Value)
vanhouten.reduce$Log.Pvalue <- -(log10(vanhouten.reduce$P.Value)) * sign(vanhouten.reduce$logFC)
colnames(vanhouten.reduce) <- str_c("vanhouten", colnames(vanhouten.reduce), sep = ".")
vanhouten.reduce$Symbol <- vanhouten$Symbol
vanhouten.reduce %<>% filter(!duplicated(Symbol))

#Only keep genes shared between DE and Van Houten data
vanhouten.overlap <- intersect(all.human.reduce$Symbol, vanhouten.reduce$Symbol)
human.vanhouten <- filter(all.human.reduce, Symbol %in% vanhouten.overlap)
vanhouten.human <- filter(vanhouten.reduce, Symbol %in% vanhouten.overlap)

#Compute RRHO between DE and Van Houten data
vanhoutenh.rrho <- map(human.groups, GetRRHO, "vanhouten", human.vanhouten, vanhouten.human, "Symbol")
vanhoutenh.rrho.logpval <- map(vanhoutenh.rrho, getElement, "hypermat.by") %>% map_dbl(max) 
vanhoutenh.rrho.pval <- exp(-vanhoutenh.rrho.logpval)

#Only keep genes shared between GAA and Van Houten data
vanhouten.overlap.gaa <- intersect(gaa.genes$Symbol, vanhouten.reduce$Symbol)
gaa.vanhouten <- filter(gaa.genes, Symbol %in% vanhouten.overlap)
colnames(gaa.vanhouten)[3:5] <- str_c("gaa", colnames(gaa.vanhouten)[3:5], sep = ".")
gaa.vanhouten$gaa.Log.Pvalue <- -(log10(gaa.vanhouten$gaa.P.value)) * sign(gaa.vanhouten$gaa.Correlation)
vanhouten.gaa <- filter(vanhouten.reduce, Symbol %in% vanhouten.overlap)

#Compute RRHO between DE and Van Houten data
vanhouteng.rrho <- GetRRHO("gaa", "vanhouten", gaa.vanhouten, vanhouten.gaa, "Symbol")
vanhouteng.rrho.logpval <- na.omit(vanhouteng.rrho$hypermat.by) %>% max
vanhouteng.rrho.pval <- exp(-vanhouteng.rrho.logpval)

#Combine all results
pval.rrho <- data.frame(kiko.rrho.pval, pbmch.rrho.pval, FRDA.vs.Healthy.ipsc = ipsch.rrho.pval, FRDA.vs.Healthy.vanhouten = vanhoutenh.rrho.pval)

gaa.pval.rrho <- c(kikog.rrho.pval, pbmch.rrho.pval, ipscg.rrho.pval, vanhouteng.rrho.pval)
pval.rrho.all <- rbind(pval.rrho, gaa = gaa.pval.rrho)

pval.adjust.rrho <- map(pval.rrho.all, p.adjust, method = "fdr", n = (nrow(pval.rrho.all) * ncol(pval.rrho.all))) %>% reduce(cbind) %>% data.frame
colnames(pval.adjust.rrho) <- colnames(pval.rrho.all)
rownames(pval.adjust.rrho) <- rownames(pval.rrho.all)
pval.adjust.rrho$Human.Comparison <- c("Carrier vs. Control", "Patient vs. Carrier", "Patient vs. Control", "GAA1 Correlation")
pval.adjust.rrho$Human.Comparison %<>% factor(levels = pval.adjust.rrho$Human.Comparison)
pval.adjust.plot <- gather(pval.adjust.rrho, Data.Comparison, P.value, -Human.Comparison)
pval.adjust.plot$P.value %<>% signif(3)
pval.adjust.plot$Data.Comparison %<>% factor %>% revalue(c(KIKO_vs_WT.heart = "GSE15848 (Heart)", KIKO_vs_WT.liver = "GSE15848 (Liver)", KIKO_vs_WT.muscle = "GSE15848 (Muscle)", pco.pbmc = "GSE30933 (Patient vs. Control)", pca.pbmc = "GSE30933 (Patient vs. Carrier)", cc.pbmc = "GSE30933 (Carrier vs. Control)", FRDA.vs.Healthy.ipsc = "GSE65399", FRDA.vs.Healthy.vanhouten = "GSE11204"))
pval.adjust.plot$Data.Comparison %<>% factor(levels = unique(pval.adjust.plot$Data.Comparison))

quantile.range <- quantile(pval.adjust.plot$P.value, probs = seq(0, 1, 0.05))
color.palette <- colorRampPalette(c("#FFFFFF", "#df8640"))(length(quantile.range) - 1)
pval.adjust.plot$P.value.quantile <- findInterval(pval.adjust.plot$P.value, quantile.range, all.inside = TRUE)

CairoPDF("overlap.heatmap", width = 8, height = 4, bg = "transparent")
p <- ggplot(pval.adjust.plot, aes(Data.Comparison, Human.Comparison, label = P.value, fill = factor(P.value.quantile))) + geom_raster() + geom_text() + theme_bw()
p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(panel.background = element_blank(), panel.border = element_blank())
p <- p + theme(plot.background = element_blank())
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
p <- p + scale_fill_manual(values = color.palette)
print(p)
dev.off()

 
