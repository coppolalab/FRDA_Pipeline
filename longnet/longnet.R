#String operations
library(stringr)

#Reading and writing tables
library(readr)
library(openxlsx)

#For plotting
library(ggplot2)

#For microarray stuff
library(lumi)
library(lumiHumanAll.db)

#Network Analysis
library(minerva)
library(minet)
library(moments)
library(Rgraphviz)
library(flashClust)
library(WGCNA)
enableWGCNAThreads()

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(Cairo)
library(reshape2)

#Data arrangement
library(dplyr)
library(tidyr)

source("../common_functions.R")

enrichr.submit <- function(index, full.df, enrichr.terms, use.weights)
{
    dataset <- full.df[[index]]
    dir.create(file.path("./enrichr", index), showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), enrichr.wkbk, enrichr.data, index)
}

enrichr.wkbk <- function(subindex, full.df, index)
{
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    dir.create(file.path("./enrichr", index))
    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

stringdb.submit <- function(module.color, module.df)
{
    module.genes <- module.df[[module.color]]
    get.stringdb(module.genes, module.color, "./stringdb")
}
intensities.means <- readRDS.gz("../dtw/save/intensities.means.rda")
intensities.ica.genes <- readRDS.gz("../ica/save/intensities.ica.genes.rda")
dtw.cc <- readRDS.gz("../dtw/save/dtw.cc.rda")
dtw.pca <- readRDS.gz("../dtw/save/dtw.pca.rda")
dtw.pco <- readRDS.gz("../dtw/save/dtw.pco.rda")
fdata <- readRDS.gz("../dtw/save/fdata.rda")

genes <- list()
#intensities.patient <- intensities.means$Patient
intensities.patient <- readRDS.gz("../ica/save/expr.means.rda")
ica.patient <- intensities.ica.genes$Patient %>% map(select, Symbol) %>% map(Compose(unlist, as.character)) %>% reduce(c) #%>% unique
dtw.patient <- unique(c(dtw.pca[1:300,]$Symbol, dtw.pco[1:300,]$Symbol)) #%>% str_replace("\\-", "\\.")
dtw.patient <- dtw.patient[dtw.patient %in% rownames(intensities.patient)]
genes$Patient <- unique(c(ica.patient, dtw.patient))
saveRDS.gz(ica.patient, "./save/ica.patient.rda")
saveRDS.gz(genes$Patient, "./save/all.patient.rda")

intensities.carrier <- intensities.means$Carrier
ica.carrier <- intensities.ica.genes$Carrier %>% map(select, Symbol) %>% map(Compose(unlist, as.character)) %>% reduce(c) %>% unique
dtw.carrier <- unique(c(dtw.pca[1:300,]$Symbol, dtw.cc[1:300,]$Symbol)) %>% str_replace("\\-", "\\.")
genes$Carrier <- unique(c(ica.carrier, dtw.carrier))
saveRDS.gz(all.carrier, "./save/all.carrier.rda")

intensities.control <- intensities.means$Control
ica.control <- intensities.ica.genes$Control %>% map(select, Symbol) %>% map(Compose(unlist, as.character)) %>% reduce(c) %>% unique
dtw.control <- unique(c(dtw.cc[1:2000,]$Symbol, dtw.pco[1:2000,]$Symbol)) %>% str_replace("\\-", "\\.")
genes$Control <- unique(c(ica.control, dtw.control))
saveRDS.gz(all.control, "./save/all.control.rda")

test.name <- rownames(intensities.patient) %>% toupper %>% str_replace_all("\\-", "\\.")
test.name2 <- toupper(genes$Patient)

mi.expr <- intensities.patient[match(genes$Patient, rownames(intensities.patient)),]
mi.expr.random <- apply(mi.expr, 1, sample, 4) 

#mi.empirical best
#mi.empirical R^2 tie, steeper slope
mi.mrnet <- minet(t(mi.expr), method = "mrnet", estimator = "mi.shrink", disc = "equalwidth")
TOM.PEER <- TOMsimilarity(mi.mrnet, verbose = 5)
dissimilarity.TOM <- 1 - TOM.PEER

CairoPDF("scalefree", height = 6, width = 6)
scaleFreePlot(TOM.PEER) #Meh
dev.off()

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")

CairoPDF("./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 20
#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, cutHeight = 0.995, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF("./gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

symbol.vector <- rownames(mi.expr) %>% str_replace("\\.", "\\-")
module.expr.out <- data.frame(Symbol = symbol.vector, module.color = dynamic.colors)
write.xlsx(module.expr.out, "./module_status.xlsx", sep = ".")
return(module.expr.out)

mi.networks <- lapply(names(genes), gen.mi, intensities.means, genes)
mi.network <- gen.mi("Patient", intensities.means, genes)

module.patients <- read.xlsx("./module_status.xlsx")
module.symbols <- split(module.patients, factor(module.patients$module.color))

source("../../code/GO/enrichr.R")
map(names(module.symbols), stringdb.submit, module.symbols)
brown.stringdb <- get.stringdb(module.symbols$brown, "brown")
blue.stringdb <- get.stringdb(module.symbols$blue, "blue", edge.threshold = 700)
red.stringdb <- get.stringdb(module.symbols$red, "red", edge.threshold = 700)

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "Humancyc_2016", "NCI-Nature_2016", "Panther_2016") 

trap1 <- map(names(module.symbols)[-c(1:3)], enrichr.submit, module.symbols, enrichr.terms, FALSE)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#black GO forget it
black.gobiol.file <- "./enrichr/black/black_GO_Biological_Process_2015.xlsx"
black.gobiol <- read.xlsx(black.gobiol.file)
black.gobiol$Database <- "GO Biological Process"
black.gobiol$Num.Genes <- map(black.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.gobiol.filter <- filter(black.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(black.gobiol.filter, module.symbols$black$Symbol, file_path_sans_ext(black.gobiol.file))

black.reactome.file <- "./enrichr/black/black_Reactome_2016.xlsx"
black.reactome <- read.xlsx(black.reactome.file)
black.reactome$Database <- "Reactome"
black.reactome$Num.Genes <- map(black.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
black.reactome.filter <- filter(black.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(black.reactome.filter, module.symbols$black$Symbol, file_path_sans_ext(black.reactome.file))

#blue GO
blue.gobiol.file <- "./enrichr/blue/blue_GO_Biological_Process_2015.xlsx"
blue.gobiol <- read.xlsx(blue.gobiol.file)
blue.gobiol$Database <- "GO Biological Process"
blue.gobiol$Num.Genes <- map(blue.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
blue.gobiol.filter <- filter(blue.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(blue.gobiol.filter, module.symbols$blue$Symbol, file_path_sans_ext(blue.gobiol.file))
blue.gobiol.final <- slice(blue.gobiol, c(4, 57, 66, 77))

blue.reactome.file <- "./enrichr/blue/blue_Reactome_2016.xlsx"
blue.reactome <- read.xlsx(blue.reactome.file)
blue.reactome$Database <- "Reactome"
blue.reactome$Num.Genes <- map(blue.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
blue.reactome.filter <- filter(blue.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(blue.reactome.filter, module.symbols$blue$Symbol, file_path_sans_ext(blue.reactome.file))
blue.reactome.final <- slice(blue.reactome, 33)

blue.enrichr <- rbind(blue.gobiol.final, blue.reactome.final)
gen.enrichrplot(blue.enrichr, "blue.enrichr")

#brown GO ??? maybe
brown.gobiol.file <- "./enrichr/brown/brown_GO_Biological_Process_2015.xlsx"
brown.gobiol <- read.xlsx(brown.gobiol.file)
brown.gobiol$Database <- "GO Biological Process"
brown.gobiol$Num.Genes <- map(brown.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
brown.gobiol.filter <- filter(brown.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(brown.gobiol.filter, module.symbols$brown$Symbol, file_path_sans_ext(brown.gobiol.file))

#green GO
green.gobiol.file <- "./enrichr/green/green_GO_Biological_Process_2015.xlsx"
green.gobiol <- read.xlsx(green.gobiol.file)
green.gobiol$Database <- "GO Biological Process"
green.gobiol$Num.Genes <- map(green.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
green.gobiol.filter <- filter(green.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(green.gobiol.filter, module.symbols$green$Symbol, file_path_sans_ext(green.gobiol.file))

#red GO
red.gobiol.file <- "./enrichr/red/red_GO_Biological_Process_2015.xlsx"
red.gobiol <- read.xlsx(red.gobiol.file)
red.gobiol$Database <- "GO Biological Process"
red.gobiol$Num.Genes <- map(red.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
red.gobiol.filter <- filter(red.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(red.gobiol.filter, module.symbols$red$Symbol, file_path_sans_ext(red.gobiol.file))

red.reactome.file <- "./enrichr/red/red_Reactome_2016.xlsx"
red.reactome <- read.xlsx(red.reactome.file)
red.reactome$Database <- "Reactome"
red.reactome$Num.Genes <- map(red.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
red.reactome.filter <- filter(red.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(red.reactome.filter, module.symbols$red$Symbol, file_path_sans_ext(red.reactome.file))

#turquoise GO
turquoise.gobiol.file <- "./enrichr/turquoise/turquoise_GO_Biological_Process_2015.xlsx"
turquoise.gobiol <- read.xlsx(turquoise.gobiol.file)
turquoise.gobiol$Database <- "GO Biological Process"
turquoise.gobiol$Num.Genes <- map(turquoise.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
turquoise.gobiol.filter <- filter(turquoise.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(turquoise.gobiol.filter, module.symbols$turquoise$Symbol, file_path_sans_ext(turquoise.gobiol.file))
turquoise.gobiol.final <- slice(turquoise.gobiol, c(1, 32, 37, 94, 103))

turquoise.reactome.file <- "./enrichr/turquoise/turquoise_Reactome_2016.xlsx"
turquoise.reactome <- read.xlsx(turquoise.reactome.file)
turquoise.reactome$Database <- "Reactome"
turquoise.reactome$Num.Genes <- map(turquoise.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
turquoise.reactome.filter <- filter(turquoise.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(turquoise.reactome.filter, module.symbols$turquoise$Symbol, file_path_sans_ext(turquoise.reactome.file))
turquoise.reactome.final <- slice(turquoise.reactome, 28)

turquoise.enrichr <- rbind(turquoise.gobiol.final, turquoise.reactome.final)
gen.enrichrplot(turquoise.enrichr, "turquoise.enrichr")

#yellow
yellow.gobiol.file <- "./enrichr/yellow/yellow_GO_Biological_Process_2015.xlsx"
yellow.gobiol <- read.xlsx(yellow.gobiol.file)
yellow.gobiol$Database <- "GO Biological Process"
yellow.gobiol$Num.Genes <- map(yellow.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
yellow.gobiol.filter <- filter(yellow.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(yellow.gobiol.filter, module.symbols$yellow$Symbol, file_path_sans_ext(yellow.gobiol.file))

