#Functional programming
library(magrittr)
library(purrr)
library(functional)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)

#Reading and writing tables
library(readr)
library(openxlsx)

#Classification
library(MASS)
library(randomForest)
library(rpart)
library(supclust)
library(class)
library(e1071)
library(pamr)
library(glmnet)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(parallel)

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

expr.PEER <- readRDS.gz(file = "../baseline/save/intensities1.peer.rda")
targets <- readRDS.gz(file = "../baseline/save/targets1.rmreps.rda") 
load(file = "../WGCNA/save/annot.reduce.rda")
status <- model.matrix( ~ 0 + targets$Status)
colnames(status) <- c("Carrier", "Control", "Patient")

up.patient <- read.xlsx("../baseline/enrichr/fdr/patient_vs_carrier_up.xlsx")
down.patient <- read.xlsx("../baseline/enrichr/fdr/patient_vs_carrier_down.xlsx")
expr.pat.up <- select(data.frame(t(expr.PEER)), one_of(up.patient$Probe_Id))
expr.pat.down <- select(data.frame(t(expr.PEER)), one_of(down.patient$Probe_Id))

rf.output.up <- randomForest(expr.pat.up, targets$Status)
rf.up.genes <- data.frame(Probe_Id = rownames(rf.output$importance), Importance = rf.output$importance) %>% join(annot.reduce) %>% arrange(MeanDecreaseGini)
rf.output.down <- randomForest(expr.pat.down, targets$Status)
rf.down.genes <- data.frame(Probe_Id = rownames(rf.output.down$importance), Importance = rf.output.down$importance) %>% join(annot.reduce) %>% arrange(MeanDecreaseGini)

expr.df.up <- mutate(expr.pat.up, Status = targets$Status)
lda.up <- sda.ranking(t(expr.PEER), targets$Status)
lda.df <- data.frame(Probe_Id = rownames(lda.up), Ranking = lda.up) %>% join(annot.reduce) %>% sort(desc(Ranking))
lda.1000 <- slice(lda.up, 1:1000)
