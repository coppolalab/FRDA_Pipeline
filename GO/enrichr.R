library(httr)
library(readr)
library(openxlsx)
library(plyr)
library(dplyr)
library(stringr)
library(jsonlite)
library(purrr)
library(reshape2)
library(magrittr)

get.enrichrdata <- function(database, gene.df, use.weights)
{
    mainurl <- "http://amp.pharm.mssm.edu/Enrichr"
    if (use.weights == TRUE)
    {
        gene.list <- select(gene.df, Symbol, Weight)
        gene.list.combined <- paste(gene.list$Symbol, gene.list$Weight, sep = ",") 
        gene.list.format <- paste(gene.list.combined, collapse = "\n")
    }
    else
    {
        gene.list <- select(gene.df, Symbol) %>% as.matrix %>% as.vector
        gene.list.format <- paste(gene.list, collapse = "\n")
    }

    post.request <- POST(url = paste(mainurl, "addList", sep = "/"), body = list(list = gene.list.format, description = ""))
    userlist <- content(post.request, "text") %>% str_extract("[0-9]+")
    get.request <- GET(url = paste(mainurl, "enrich", sep = "/"), query = list(backgroundType = database, userListId = userlist))
    get.content <- content(get.request)[[1]]

    content.df <- lapply(get.content, reshapedata) %>% reduce(rbind) %>% data.frame
    if (ncol(content.df) == 7)
    {
        colnames(content.df) <- c("Index", "GO.Term", "P.value", "Z.score", "Combined.Score", "Genes", "Adj.P.value")
        content.df$P.value %<>% as.numeric
        content.df$Z.score %<>% as.numeric
        content.df$Combined.Score %<>% as.numeric
        content.df$Adj.P.value %<>% as.numeric
        content.df$Index %<>% as.numeric

        content.df$GO.Term %<>% as.character
        content.df$Genes %<>% as.character

        content.df %<>% select(Index:P.value, Adj.P.value, Z.score:Genes)
        return(content.df)
    }    
    else
    {
        print(paste(database, "returned no results"))
        return(NA)
    }
}

reshapedata <- function(orig.list)
{
    orig.list[[6]] %<>% paste(collapse = ",")
    return(orig.list)
}
