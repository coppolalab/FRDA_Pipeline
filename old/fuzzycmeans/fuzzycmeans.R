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
library(Mfuzz)

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

csel.repeat <- function(eset, m, crange, max.runs, csel.mat = matrix(), count.runs = 0)
{
    csel.out <- cselection(eset, m = m, crange = crange, repeats = 5, visu = FALSE)
    if (count.runs == 0)
    {
        csel.return <- csel.out
    }
    else
    {
        csel.return <- rbind(csel.mat, csel.out)
    }
    count.runs <- count.runs + 1
    if (count.runs < max.runs)
    {
        csel.repeat(eset, m, crange, max.runs, csel.return, count.runs)
    }
    else
    {
        return(csel.return)
    }
}

dmin.repeat <- function(eset, m, crange, max.runs, dmin.mat = vector(), count.runs = 0)
{
    dmin.out <- Dmin(eset, m = m, crange = crange, repeats = 5, visu = FALSE)
    if (count.runs == 0)
    {
        dmin.return <- dmin.out
    }
    else
    {
        dmin.return <- rbind(dmin.mat, dmin.out)
    }
    count.runs <- count.runs + 1
    if (count.runs < max.runs)
    {
        dmin.repeat(eset, m, crange, max.runs, dmin.return, count.runs)
    }
    else
    {
        return(dmin.return)
    }
}

match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|"))

intensities.mean <- readRDS.gz("../dtw/save/intensities.means.rda") 
intensities.patient <- intensities.mean$Patient
rownames(intensities.patient) <- intensities.patient$Symbol
intensities.patient %<>% select(-Symbol)

patient.betr.genes <- read.xlsx("../betr/patient.out.xlsx") %>% select(Symbol)
patient.timecourse.genes <- read.xlsx("../timecourse/patient.out.xlsx") %>% select(Symbol)
longnet.patient <- readRDS.gz("../longnet/save/all.patient.rda")
ica.patient <- readRDS.gz("../longnet/save/ica.patient.rda")
cluster.genes <- c(unlist(patient.betr.genes), unlist(patient.timecourse.genes)) %>% unique %>% match.exact

lumi.expr.collapse <- readRDS.gz("../dtw/save/lumi.exprs.collapse.rda") %>% t
expr.patients <- lumi.expr.collapse[,str_detect(colnames(lumi.expr.collapse), "Pat")]

split.factor <- rep(1:(ncol(expr.patients)/4), each = 4)
expr.cluster <- expr.patients[grepl(cluster.genes, rownames(expr.patients)),]
cluster.gene.expr.list <- split(data.frame(t(expr.cluster)), split.factor) %>% map(t)

patient.genes <- match.exact(ica.patient)
patient.cluster <- intensities.patient[grepl(patient.genes, rownames(intensities.patient)),] 
#patient.cluster <- intensities.patient[match(longnet.patient, rownames(intensities.patient)),] 
patient.eset <- ExpressionSet(assayData = as.matrix(patient.cluster)) %>% standardise

m.estimate <- mestimate(patient.eset) 

csel.runs <- csel.repeat(patient.eset, m = m.estimate, crange = 4:20, max.runs = 5)
csel.ratio <- (4:20) - colMeans(csel.runs) 

dmin.runs <- dmin.repeat(patient.eset, m = m.estimate, crange = seq(4, 40, 4), max.runs = 5)

seed.mfuzz <- function(eset, c.num, m, mfuzz.list = list(), iter.count = 0)
{
   if (iter.count < 250) 
   {
        mfuzz.new <- mfuzz(eset, c = c.num, m = m)
        iter.count <- iter.count + 1
        print(iter.count)
        mfuzz.add <- mfuzz.new$membership
        colnames(mfuzz.add) <- paste("X", 1:c.num, sep = "")
        mfuzz.list[[iter.count]] <- mfuzz.add
        seed.mfuzz(eset = eset, c.num = c.num, m = m, mfuzz.list = mfuzz.list, iter.count = iter.count)
   }
   else
   {
        return(mfuzz.list)
   }
}

rename.columns <- function(mfuzz.object)
{
    center.colnames <- colnames(mfuzz.object$centers)[1]
    patient.name <- paste(str_split(center.colnames[1], "_")[[1]][1:2], collapse = "_")
    mfuzz.membership <- mfuzz.object$membership
    colnames(mfuzz.membership) <- paste(patient.name, colnames(mfuzz.membership), sep = "_")
    return(mfuzz.membership)
}

map.cluster <- map(cluster.gene.expr.list, Compose(ExpressionSet, standardise)) 
map.mestimate <- map_dbl(map.cluster, mestimate)
nclust <- 7
map.mfuzz <- map(map.cluster, mfuzz, c = nclust, m = map.mestimate[1])
map.renamed <- map(map.mfuzz, rename.columns) %>% reduce(cbind) %>% t
collapse.membership <- collapseRows(map.renamed, rep(1:nclust, 38), rownames(map.renamed), methodFunction = colMedians)$datETcollapsed %>% t %>% data.frame
module.membership <- apply(collapse.membership, 1, which.max)
collapse.membership$cluster <- module.membership
cluster.members <- select(collapse.membership, cluster)

summary(factor(module.membership))

cluster.patient7 <- seed.mfuzz(eset = patient.eset, c.num = 7, m = m.estimate)
median.patient7 <- melt(cluster.patient7) %>% dcast(Var1 ~ Var2, median)

cluster.patient13 <- seed.mfuzz(eset = patient.eset, c.num = 13, m = m.estimate)
median.patient13 <- melt(cluster.patient13) %>% dcast(Var1 ~ Var2, median)

cluster.patient14 <- seed.mfuzz(eset = patient.eset, c.num = 14, m = m.estimate)
median.patient14 <- melt(cluster.patient14) %>% dcast(Var1 ~ Var2, median)

cluster.patient <- mfuzz(patient.eset, c = 7, m = m.estimate)
mfuzz.plot(patient.eset, cl = cluster.patient, mfrow = c(4,4), time.labels = 1:4)
#cluster.filter <- (cluster.patient$membership > 0.4) %>% apply(1, any)
#cluster.filtered <- cluster.patient$membership[cluster.filter,]
#patient.members <- apply(cluster.patient$membership, 1, which.max)
patient.df <- data.frame(Symbol = rownames(cluster.members), Cluster = cluster.members)
write.xlsx(patient.df, "./patient.cmeans.xlsx")

partcoef.patient <- partcoef(patient.eset)

