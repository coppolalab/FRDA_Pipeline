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

lumi.exprs.collapse <- readRDS.gz("../dtw/save/lumi.exprs.collapse.rda") %>% t
lumi.patient <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Pat")]

patient.reps <- rep(ncol(lumi.patient)/4, nrow(lumi.patient))
patient.repgrp <- rep(1:(ncol(lumi.patient)/4), each = 4)
patient.timegrp <- rep(1:4, (ncol(lumi.patient)/4))
patient.out <- mb.long(lumi.patient, method = "1D", type = "robust", times = 4, reps = patient.reps, rep.grp = patient.repgrp, time.grp = patient.timegrp)
saveRDS.gz(patient.out, "./save/patient.out.rda")

patient.genes <- data.frame(Symbol = rownames(lumi.exprs.collapse), Hotelling = patient.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% filter(Hotelling > 10)
write.xlsx(patient.genes, "./patient.out.xlsx")

lumi.carrier <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Car")]
carrier.reps <- rep(ncol(lumi.carrier)/4, nrow(lumi.carrier))
carrier.repgrp <- rep(1:(ncol(lumi.carrier)/4), each = 4)
carrier.timegrp <- rep(1:4, (ncol(lumi.carrier)/4))
carrier.out <- mb.long(lumi.carrier, method = "1D", type = "robust", times = 4, reps = carrier.reps, rep.grp = carrier.repgrp, time.grp = carrier.timegrp)
saveRDS.gz(carrier.out, "./save/carrier.out.rda")

carrier.genes <- data.frame(Symbol = rownames(lumi.exprs.collapse), Hotelling = carrier.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% filter(Hotelling > 13)
write.xlsx(carrier.genes, "./carrier.out.xlsx")

lumi.control <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Con")]
control.reps <- rep(ncol(lumi.control)/4, nrow(lumi.control))
control.repgrp <- rep(1:(ncol(lumi.control)/4), each = 4)
control.timegrp <- rep(1:4, (ncol(lumi.control)/4))
control.out <- mb.long(lumi.control, method = "1D", type = "robust", times = 4, reps = control.reps, rep.grp = control.repgrp, time.grp = control.timegrp)
saveRDS.gz(control.out, "./save/control.out.rda")

control.genes <- data.frame(Symbol = rownames(lumi.exprs.collapse), Hotelling = control.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% filter(Hotelling > 15)
write.xlsx(control.genes, "./control.out.xlsx")

lumi.pca <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Pat|Car")]
pca.repgrp <- rep(1:(ncol(lumi.pca)/4), each = 4)
pca.timegrp <- rep(1:4, (ncol(lumi.pca)/4))
pca.condgrp <- str_detect(colnames(lumi.pca), "Pat") %>% as.numeric
pca.reps <- cbind(rep(length(which(pca.condgrp == 0))/4, nrow(lumi.pca)), rep(length(which(pca.condgrp == 1))/4, nrow(lumi.pca)))
pca.out <- mb.long(lumi.pca, method = "2D", type = "robust", times = 4, reps = pca.reps, rep.grp = pca.repgrp, time.grp = pca.timegrp, condition.grp = pca.condgrp)
saveRDS.gz(pca.out, "./save/pca.out.rda")

pca.genes <- data.frame(Symbol = rownames(lumi.exprs.collapse), Hotelling = pca.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% filter(Hotelling > 10)
write.xlsx(pca.genes, "./pca.out.xlsx")

lumi.pco <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Pat|Con")]
pco.repgrp <- rep(1:(ncol(lumi.pco)/4), each = 4)
pco.timegrp <- rep(1:4, (ncol(lumi.pco)/4))
pco.condgrp <- str_detect(colnames(lumi.pco), "Pat") %>% as.numeric
pco.reps <- cbind(rep(length(which(pco.condgrp == 0))/4, nrow(lumi.pco)), rep(length(which(pco.condgrp == 1))/4, nrow(lumi.pco)))
pco.out <- mb.long(lumi.pco, method = "2D", type = "robust", times = 4, reps = pco.reps, rep.grp = pco.repgrp, time.grp = pco.timegrp, condition.grp = pco.condgrp)
saveRDS.gz(pco.out, "./save/pco.out.rda")

pco.genes <- data.frame(Symbol = rownames(lumi.exprs.collapse), Hotelling = pco.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% filter(Hotelling > 10)
write.xlsx(pco.genes, "./pco.out.xlsx")

lumi.cc <- lumi.exprs.collapse[,str_detect(colnames(lumi.exprs.collapse), "Car|Con")]
cc.repgrp <- rep(1:(ncol(lumi.cc)/4), each = 4)
cc.timegrp <- rep(1:4, (ncol(lumi.cc)/4))
cc.condgrp <- str_detect(colnames(lumi.cc), "Car") %>% as.numeric
cc.reps <- cbind(rep(length(which(cc.condgrp == 0))/4, nrow(lumi.cc)), rep(length(which(cc.condgrp == 1))/4, nrow(lumi.cc)))
cc.out <- mb.long(lumi.cc, method = "2D", type = "robust", times = 4, reps = cc.reps, rep.grp = cc.repgrp, time.grp = cc.timegrp, condition.grp = cc.condgrp)
saveRDS.gz(cc.out, "./save/cc.out.rda")

cc.genes <- data.frame(Symbol = rownames(lumi.exprs.collapse), Hotelling = cc.out$HotellingT2) %>% arrange(desc(Hotelling)) %>% filter(Hotelling > 14)
write.xlsx(cc.genes, "./cc.out.xlsx")

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort
