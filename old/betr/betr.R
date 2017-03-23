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
library(limma)

#Longitudinal analysis
library(betr)

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

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(vadr)

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

lumi.exprs.collapse <- readRDS.gz("../dtw/save/lumi.collapse.rda")# %>% t
lumi.patient <- readRDS.gz("../ica/save/lumi.collapse.rda")

patient.repgrp <- rep(1:(ncol(lumi.patient)/4), each = 4)
patient.timegrp <- rep(1:4, (ncol(lumi.patient)/4))
patient.out <- betr(lumi.patient, timepoint = pData(lumi.patient)$Sample.Num, replicate = as.integer(factor(pData(lumi.patient)$PIDN)), twoCondition = FALSE)
saveRDS.gz(patient.out, "./save/patient.out.rda")

patient.genes <- data.frame(Symbol = names(patient.out), Probability = patient.out) %>% arrange(desc(Probability))
write.xlsx(patient.genes, "./patient.out.xlsx")

lumi.carrier <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Car")]
carrier.repgrp <- rep(1:(ncol(lumi.carrier)/4), each = 4)
carrier.timegrp <- rep(1:4, (ncol(lumi.carrier)/4))
carrier.out <- betr(lumi.carrier, timepoint = carrier.timegrp, replicate = carrier.repgrp, twoCondition = FALSE)
saveRDS.gz(carrier.out, "./save/carrier.out.rda")

carrier.genes <- data.frame(Symbol = rownames(lumi.exprs.collapse), Probability = carrier.out) %>% arrange(desc(Probability)) 
write.xlsx(carrier.genes[1:500,], "./carrier.out.xlsx")

lumi.control <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Con")]
control.repgrp <- rep(1:(ncol(lumi.control)/4), each = 4)
control.timegrp <- rep(1:4, (ncol(lumi.control)/4))
control.out <- betr(lumi.control, timepoint = control.timegrp, replicate = control.repgrp, twoCondition = FALSE)
saveRDS.gz(control.out, "./save/control.out.rda")

control.genes <- data.frame(Symbol = names(control.out), Probability = control.out) %>% arrange(desc(Probability))
write.xlsx(control.genes[1:500,], "./control.out.xlsx")

lumi.pca <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Pat|Car")]
pca.repgrp <- rep(1:(ncol(lumi.pca)/4), each = 4)
pca.timegrp <- rep(1:4, (ncol(lumi.pca)/4))
pca.condgrp <- str_detect(colnames(lumi.pca), "Pat") %>% as.numeric
pca.out <- betr(exprs(lumi.pca), cond = factor(lumi.pca$Status), timepoint = factor(lumi.pca$Sample.Num), replicate = factor(lumi.pca$PIDN), twoCondition = TRUE, verbose = TRUE)
saveRDS.gz(pca.out, "./save/pca.out.rda")

pca.genes <- data.frame(Symbol = names(pca.out), Probability = pca.out) %>% arrange(desc(Probability))
write.xlsx(pca.genes, "./pca.out.xlsx")

lumi.pco <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Pat|Con")]
pco.repgrp <- rep(1:(ncol(lumi.pco)/4), each = 4)
pco.timegrp <- rep(1:4, (ncol(lumi.pco)/4))
pco.condgrp <- str_detect(colnames(lumi.pco), "Pat") %>% as.numeric
pco.out <- betr(lumi.pco, cond = pco.condgrp, timepoint = pco.timegrp, replicate = pco.repgrp, twoCondition = TRUE)
saveRDS.gz(pco.out, "./save/pco.out.rda")

pco.genes <- data.frame(Symbol = names(pco.out), Probability = pco.out) %>% arrange(desc(Probability))
write.xlsx(pco.genes, "./pco.out.xlsx")

lumi.cc <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Car|Con")]
cc.repgrp <- rep(1:(ncol(lumi.cc)/4), each = 4)
cc.timegrp <- rep(1:4, (ncol(lumi.cc)/4))
cc.condgrp <- str_detect(colnames(lumi.cc), "Car") %>% as.numeric
cc.out <- betr(lumi.cc, cond = cc.condgrp, timepoint = cc.timegrp, replicate = cc.repgrp, twoCondition = TRUE)
saveRDS.gz(cc.out, "./save/cc.out.rda")

cc.genes <- data.frame(Symbol = names(cc.out), Probability = cc.out) %>% arrange(desc(Probability))
write.xlsx(cc.genes, "./cc.out.xlsx")

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort
