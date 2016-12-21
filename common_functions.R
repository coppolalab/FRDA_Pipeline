library(irr)
library(tools)
library(heatmap.plus)
library(WGCNA)
library(igraph)
library(TeachingDemos)
library(readr)
enableWGCNAThreads()

biogrid.ppi <- read_tsv("~/Documents/_Research/code/BIOGRID-ORGANISM-3.4.136.tab2/BIOGRID-ORGANISM-Homo_sapiens-3.4.136.tab2.txt") %>% data.frame
inweb.ppi <- read_tsv("~/Documents/_Research/code/InWeb3_HC_NonRed.txt", col_names = FALSE)
iw.hugo <- read_tsv("~/Documents/_Research/code/IWtoHugo.txt", col_names = FALSE)

SaveRDSgz <- function(object,file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

ReadRDSgz <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

GetKappa <- function(term.current, all.terms) {
    map(all.terms, cbind, term.current) %>% map(kappa2) %>% map_dbl(getElement, "value")
}

GetKappaCluster <- function(filename, enrichr.output, gene.names, sig.measure = "P.value") {
    num.genes <- length(gene.names)
    enrichr.list <- map(enrichr.output$Genes, str_split, ",") %>% map(getElement, 1) 
    enrichr.match <- map(enrichr.list, is.element, el = toupper(gene.names)) %>% reduce(rbind) %>% t
    rownames(enrichr.match) <- toupper(gene.names)
    colnames(enrichr.match) <- enrichr.output$Term
    enrichr.match.df <- data.frame(enrichr.match)

    enrichr.kappa <- map(enrichr.match.df, GetKappa, enrichr.match.df) %>% reduce(rbind)
    enrichr.kappa[enrichr.kappa < 0] <- 0

    rownames(enrichr.kappa) <- colnames(enrichr.kappa) <- enrichr.output$Term

    CairoPDF(str_c(filename, "heatmap", sep = "."), width = 30, height = 30)
    heatmap.plus(enrichr.kappa, col = heat.colors(40), symm = TRUE, margins = c(20,20))
    dev.off()

    kappa.dist <- dist(enrichr.kappa, method = "manhattan")
    kappa.clust <- hclust(kappa.dist, method = "average")

    CairoPDF(str_c(filename, "clust", sep = "."), height = 30, width = 30)
    plot(kappa.clust)
    dev.off()

    kappa.modules <- cutreeDynamic(kappa.clust, minClusterSize = 2, method = "tree")
    kappa.modules.df <- data.frame(Term = rownames(enrichr.kappa), Module = kappa.modules)

    enrichr.output$Module <- kappa.modules
    if (sig.measure == "Bayes.Factor"){
        enrichr.output %<>% arrange_(str_c("desc(", sig.measure, ")")) 
        enrichr.output$Index <- 1:nrow(enrichr.output)
    }
    else{
        enrichr.output %<>% arrange_(sig.measure) 
    }
    enrichr.output %<>% select_("Index", "Term", sig.measure, "Num.Genes", "Module", "Genes") %>% arrange(Module)
    
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = enrichr.output, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 1, cols = c(1, 3:ncol(enrichr.output)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    
    saveWorkbook(wb, str_c(filename, "table.xlsx", sep = "."), overwrite = TRUE) 
}

GetPPI <- function(gene.list) {
    biogrid.ppi.reduce <- select(biogrid.ppi, contains("Official"))

    colnames(inweb.ppi) <- c("Interactor.A", "Interactor.B")
    colnames(iw.hugo) <- c("IW.ID", "HUGO")

    ppi.id <- filter(iw.hugo, is.element(HUGO, gene.list))$IW.ID
    ppi.inweb.interactors <- filter(inweb.ppi, is.element(Interactor.A, ppi.id)) %>% filter(is.element(Interactor.B, ppi.id))
    ppi.inweb.symbols <- merge(ppi.inweb.interactors, iw.hugo, by.x = "Interactor.A", by.y = "IW.ID") %>% merge(iw.hugo, by.x = "Interactor.B", by.y = "IW.ID")
    ppi.inweb.final <- select(ppi.inweb.symbols, contains("HUGO"))
    colnames(ppi.inweb.final) <- c("Symbol.A", "Symbol.B")

    ppi.biogrid <- filter(biogrid.ppi.reduce, Official.Symbol.Interactor.A %in% gene.list) %>% filter(Official.Symbol.Interactor.B %in% gene.list)
    colnames(ppi.biogrid) <- c("Symbol.A", "Symbol.B")

    ppi.combined <- rbind(ppi.inweb.final, ppi.biogrid)
    ppi.combined$unique <- str_c(ppi.combined$Symbol.A, ppi.combined$Symbol.B, sep = ".")
    ppi.unique <- filter(ppi.combined, !duplicated(unique)) %>% select(-unique)
    ppi.self <- apply(ppi.unique, 1, reduce, identical)
    ppi.unique.final <- filter(ppi.unique, !ppi.self)
    return(ppi.unique.final)
}

#PPI Plot
PlotPPI <- function(adjacency.expr, gene.list, ppi.edges, filename, plot = TRUE, prune = FALSE, expr = TRUE, clust.keep = 1, plot.width = 7, plot.height = 7, vertex.size = 35) {
    gene.list %<>% str_split(",") %>% extract2(1)
    expr.adjacency <- adjacency.expr[gene.list,gene.list]
    ppi.adjacency <- matrix(0, ncol = ncol(expr.adjacency), nrow = nrow(expr.adjacency), dimnames = list(rownames(expr.adjacency), colnames(expr.adjacency)))
    ppi.adjacency[as.matrix(ppi.edges)] <- 1
    final.adjacency <- expr.adjacency * ppi.adjacency
    final.mins <- apply(final.adjacency, 2, min)
    final.maxs <- apply(final.adjacency, 2, max)
    final.scaled <- scale(final.adjacency, center = final.mins, scale = final.maxs - final.mins)
    final.scaled[is.nan(final.scaled)] <- 0

    final.igraph <- graph_from_adjacency_matrix(expr.adjacency, mode = "undirected", weighted = TRUE, diag = FALSE)
    ppi.igraph <- graph_from_adjacency_matrix(final.scaled, mode = "undirected", weighted = TRUE, diag = FALSE)

    if (prune == TRUE) {
        #num.edges <- map(1:vcount(ppi.igraph), incident, graph = ppi.igraph) %>% map_dbl(length)

        final.igraph <- delete.vertices(final.igraph, which(clusters(ppi.igraph)$membership != clust.keep))
        ppi.igraph <- delete.vertices(ppi.igraph, which(clusters(ppi.igraph)$membership != clust.keep))
    }

    ppi.colors <- rainbow(length(unique(clusters(ppi.igraph)$membership)))
    final.colors <- ppi.colors[clusters(ppi.igraph)$membership]

    final.edges <- attr(E(final.igraph), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_")
    final.df <- data.frame(Edges = final.edges, Weight = edge_attr(final.igraph, "weight"))
    if (sum(final.scaled) > 0) {
        ppi.edges <- attr(E(ppi.igraph), "vnames") %>% str_split("\\|") %>% map(sort) %>% reduce(rbind) %>% apply(1, paste, collapse = "_")
        ppi.df <- data.frame(Edges = ppi.edges, Weight = edge_attr(ppi.igraph, "weight"))

        final.filter <- is.element(final.df$Edges, ppi.df$Edges)
        final.df[final.filter,]$Weight <- ppi.df$Weight 
        final.df[!final.filter,]$Weight <- 0
        final.df$Color <- "#dddddd99"
        final.df[final.filter,]$Color <- "#0000FF99"
    } else {
        final.df <- list(Color = "#0000FF99", Weight = 0)
    }

    if (plot == TRUE) {
        CairoPDF(str_c(filename, "_ppi"), width = plot.width, height = plot.height)
        par(mar=c(0,0,0,0) + 0.5)
        plot.igraph(final.igraph, layout = layout_nicely(final.igraph), vertex.size = vertex.size, vertex.label.degree = pi/2, vertex.label.font = 2, vertex.label.color = "black", edge.color = final.df$Color, edge.width = 5*final.df$Weight)
        dev.off()

        if (expr == TRUE){
            CairoPDF(str_c(filename, "_expr"), width = plot.width, height = plot.height)
            par(mar=c(0,0,0,0) + 0.5)
            plot.igraph(final.igraph, layout = layout_nicely(final.igraph), vertex.size = vertex.size, vertex.label.degree = pi/2, vertex.label.font = 2, vertex.label.color = "black", edge.width = 20*edge_attr(final.igraph, "weight"))
            dev.off()
        }
    }

    return(ppi.igraph)
}
