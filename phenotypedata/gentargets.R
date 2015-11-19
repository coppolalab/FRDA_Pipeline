#Git test
library(magrittr)
library(plyr)
library(dplyr)
library(openxlsx)
library(readr)
library(stringr)
library(purrr)
library(functional)
library(lubridate)
library(R.utils)
library(reshape2)

#Find sample ranges which skip
find.skips <- function(data.vector)
{
    if (sum(data.vector) != length(data.vector))
    {
        return(TRUE)
    }
    else
    {
        return(FALSE)
    }
}

#Find samples ranges which do not start at 1
find.baseline <- function(data.vector)
{
    if (min(data.vector) > 1)
    {
        return(TRUE)
    }
    else
    {
        return(FALSE)
    }
}

#Phenotype data
my_select <- function(trained)
{
    n_fmts <- nchar(gsub("[^%]", "", names(trained))) + grepl("%y", names(trained))*1.5
    names(trained[ which.max(n_fmts) ])
}

get.families <- function(row.vector)
{
    raw.string <- c(row.vector["PIDN"], row.vector["Children"], row.vector["Parents"] , row.vector["Sibling"], row.vector["Nephew.Niece"], row.vector["Grandchild"])
    raw.string <- raw.string[!is.na(raw.string)]
    raw.string %<>% str_replace(" ", "")
    final.string <- lapply(raw.string, str_split, ",") %>% unlist
}

add.regex <- function(string.vector, one.PIDN)
{
    regex.vector <- paste("^", string.vector, sep = "") %>% paste("$", sep = "") 
    regex.key <- paste(regex.vector, collapse = "|")
    PIDN.index <- str_detect(one.PIDN, regex.key) %>% which
    return(PIDN.index)
}

get.blocks <- function(PIDN, multiple.melt)
{
    PIDN.regex <- paste("^", PIDN, "$", sep = "")
    PIDN.df <- filter(multiple.melt, grepl(PIDN.regex, value))
    return(PIDN.df$L1)
}

join.dups <- function(dup.vector, multiple.melt)
{
    regex.vector <- paste("^", dup.vector, sep = "") %>% paste("$", sep = "") %>% paste(collapse = "|")
    dup.df <- filter(multiple.melt, grepl(regex.vector, L1))
    final.PIDNs <- unique(dup.df$value) %>% as.character
    return(final.PIDNs)
}

#Read in phenotype data
chop.1 <- read.xlsx("./CURRENT_chop_rna.xlsx")
chop.1$dob %<>% str_replace("age ", "") 
fixdob <- filter(chop.1, grepl("@", dob)) %>% select(dob) %>% as.matrix %>% apply(1, str_replace_all, ' ', '') %>% str_split("@") %>% reduce(rbind) %>% as.matrix 
colnames(fixdob) <- c("Age", "Date")
adjusteddob <- parse_date_time(fixdob[,2], "mdy", select_formats = my_select) 
year(adjusteddob) <- year(adjusteddob) - as.numeric(fixdob[,1])
chop.1$dob[str_detect(chop.1$dob, '@') & !is.na(chop.1$dob)] <- as.character(as.Date(adjusteddob))

chop.2 <- read.xlsx("./Current_Request 2_Copy of chop_missing.xlsx")
chop.2$dob %<>% str_replace("age ", "") 
fixdob.2 <- filter(chop.2, grepl("@", dob)) %>% select(dob) %>% as.matrix %>% apply(1, str_replace_all, ' ', '') %>% str_split("@") %>% reduce(rbind) %>% as.matrix 
adjusteddob.2 <- parse_date_time(fixdob.2[,2], "mdy", select_formats = my_select) 
year(adjusteddob.2) <- year(adjusteddob.2) - as.numeric(fixdob.2[,1])
chop.2$dob[str_detect(chop.2$dob, '@') & !is.na(chop.2$dob)] <- as.character(as.Date(adjusteddob.2))
chop.2$dob %<>% parse_date_time(c("mdy", "ymd")) %>% as.Date %>% as.character
chop.2$Status <- str_replace(chop.2$categorized_diagnosis, "FRDA_", "")

chop.3 <- read.xlsx("./newpatients_chop_3.xlsx")
chop.3$Status %<>% capitalize
chop.3$Date.Drawn %<>% parse_date_time("mdy") %>% as.Date #%>% as.character
chop.3$DOB <- chop.3$Date.Drawn
year(chop.3$DOB) <- year(chop.3$DOB) - as.numeric(chop.3$Age)
chop.3$DOB %<>% as.character
chop.3.reduce <- select(chop.3, PIDN, Sex, Status, DOB, Age.of.onset, Repeat.1, Repeat.2)
colnames(chop.3.reduce)[5:7] <- c("Onset", "GAA1", "GAA2")
write.xlsx("chop.3.reduce.xlsx")

ucla.1 <- read.xlsx("./ucla_rna updated- w UCLA 083115 update.xlsx")
ucla.1$dob[ucla.1$dob == "MISSING"] <- NA
ucla.1$dob <- parse_date_time(ucla.1$dob, c("mdy", "ymd")) %>% as.Date %>% as.character

ucla.2 <- read.xlsx("./Copy of UCLA DATA FOR JEN 9.25.15.xlsx")
ucla.2$dob[ucla.2$dob == "MISSING"] <- NA
ucla.2$dob <- parse_date_time(ucla.2$dob, c("mdy", "ymd")) %>% as.Date %>% as.character

ucla.3 <- read.xlsx("./newpatients_ucla oct 2015.xlsx")
ucla.3$DOB <- parse_date_time(ucla.3$DOB, c("mdy", "ymd")) %>% as.Date %>% as.character 
ucla.3.reduce <- select(ucla.3, PIDN, Sex, Status, DOB, Onset, GAA1, GAA2)

#chop.all <- rbind(select(chop.1, -Order), select(chop.2, -Notes))
chop.all <- read.xlsx("./Phenotype Data 2015-11-13.xlsx", detectDates = FALSE)
chop.all$DOB %<>% parse_date_time("mdy") %>% as.Date %>% as.character

full.list <- apply(chop.all, 1, get.families)
full.lengths <- sapply(full.list, length)
multiple.PIDN <- full.list[full.lengths > 1]
one.PIDN <- full.list[full.lengths < 2] %>% reduce(c)
indices <- lapply(multiple.PIDN, add.regex, one.PIDN) %>% reduce(c) %>% unique
one.PIDN.filter <- one.PIDN[-indices]

names(multiple.PIDN) <- 1:length(multiple.PIDN)
multiple.melt <- melt(multiple.PIDN)
dup.melt <- duplicated(multiple.melt$value) | duplicated(multiple.melt$value, fromLast = TRUE)
dup.PIDN <- multiple.melt[dup.melt,]$value %>% unique
dup.blocks <- lapply(dup.PIDN, get.blocks, multiple.melt)
joined.blocks <- lapply(dup.blocks, join.dups, multiple.melt)
joined.blocks <- joined.blocks[!duplicated(joined.blocks)]
joined.blocks[[16]] <- c(joined.blocks[[16]], joined.blocks[[17]]) %>% unique
joined.blocks <- joined.blocks[-17]
dup.df <- melt(joined.blocks)
final.block <- max(dup.df$L1) + 1

clean.key <- unlist(dup.blocks) %>% unique %>% as.numeric
multiple.clean <- multiple.PIDN[-clean.key]
names(multiple.clean) <- final.block:(final.block + length(multiple.clean) - 1)
multiple.clean.df <- melt(multiple.clean)

multiple.all <- rbind(dup.df, multiple.clean.df)
colnames(multiple.all) <- c("PIDN", "Family")
multiple.all$Family %<>% as.numeric
final.block2 <- max(multiple.all$Family) + 1
one.df <- data.frame(PIDN = one.PIDN.filter, Family = final.block2:(final.block2 + length(one.PIDN.filter) - 1))
families.df <- rbind(multiple.all, one.df)

chop.families <- join(chop.all, families.df) %>% select(PIDN:GAA2, Family)
final.block3 <- max(chop.families$Family) + 1

chop.3.reduce$Family <- final.block3:(final.block3 + nrow(chop.3.reduce) - 1)
chop.families.all <- rbind(chop.families, chop.3.reduce)
chop.families.all$Site <- "CHOP"
final.block4 <- max(chop.families.all$Family) + 1

ucla.rbind <- rbind(ucla.1, ucla.2) %>% select(PIDN, sex, Status, dob, onset, gaa1, gaa2) 
colnames(ucla.rbind) <- c("PIDN", "Sex", "Status", "DOB", "Onset", "GAA1", "GAA2")
ucla.all <- rbind(ucla.rbind, ucla.3.reduce)
ucla.all$PIDN %<>% str_replace("FA_", "")
ucla.all$Family <- final.block4:(final.block4 + nrow(ucla.all) - 1)
ucla.all$Site <- "UCLA"

subjects.reduce <- rbind(chop.families.all, ucla.all)
save(subjects.reduce, file = "subjects.reduce.rda")

#Read in microarray targets
targets <- read.xlsx("./targets_chop_upd.xlsx") 
targets.key <- select(targets, SampleName, FA_PatientID, Label, Slide, Batch, SCGC.Code)
colnames(targets.key) <- c("Sample.Name", "PIDN", "Sample.Num", "Slide.ID", "Batch", "SCGC.Code")
targets.key$Sample.Num[targets.key$Sample.Num == "ng"] <- 1

#Read in new targets sheet (2014-278).  Select and rename columns as before.  
targets.new <- read_csv("./newtargets.csv")
targets.new.key <- select(targets.new, External.ID, PIDN, Sample.Num, SCGC.Slide.ID2)
colnames(targets.new.key) <- c("Sample.Name", "PIDN", "Sample.Num", "Slide.ID")
targets.new.key$PIDN %<>% str_replace_all("FA_", "")
targets.new.key$Batch <- 17
targets.new.key$SCGC.Code <- "2014-278"

#Read in newer targets sheet (2015-9185)
targets.newer <- read.xlsx("./2015-9185A_Sample Key-FA.xlsx")
targets.newer$Slide.ID <- paste(targets.newer$general.array, targets.newer$genexstripe.controling.stripe, sep = "_")
targets.newer$FA_PatientID %<>% str_replace("FA_", "")
targets.newer$Sample.Num <- str_split(targets.newer$External.ID, "_") %>% laply(`[`, 3)
targets.newer.key <- select(targets.newer, External.ID, FA_PatientID, Sample.Num, Slide.ID)
targets.newer.key$Batch <- 18
targets.newer.key$SCGC.Code <- "2015-9185"
colnames(targets.newer.key) <- c("Sample.Name", "PIDN", "Sample.Num", "Slide.ID", "Batch", "SCGC.Code")
targets.key.all <- rbind(targets.key, targets.new.key, targets.newer.key)

#Add date drawn so that age at draw can be calculated
#Had to manually add 1 missing RIN
targets.dates.full <- read_csv("../phenotypedata/dan_allrna_20151104T200216.csv")
targets.dates.full$PIDN %<>% str_replace_all("FA_", "")
targets.dates <- select(targets.dates.full, PIDN, RIN, one_of(c("RNA ID", "Date Drawn","Date Received", "Sample Num")))# %>% arrange(PIDN)
colnames(targets.dates) %<>% str_replace(" ", ".")

targets.final <- join(targets.key.all, targets.dates)
targets.final$Array.Type <- "Human HT-12 v 4.0"
reps.conditions <- targets.final$Sample.Num == "1r"
reps.names <- targets.final[reps.conditions,]$Sample.Name %>% str_replace("r", "") %>% paste(collapse = "|")
targets.final[reps.conditions,]$Date.Received <- filter(targets.final, grepl(reps.names, Sample.Name) & Sample.Num != "1r")$Date.Received
targets.final[reps.conditions,]$Date.Drawn <- filter(targets.final, grepl(reps.names, Sample.Name) & Sample.Num != "1r")$Date.Drawn
targets.final[reps.conditions,]$RIN <- filter(targets.final, grepl(reps.names, Sample.Name) & Sample.Num != "1r")$RIN
targets.final[reps.conditions,]$RNA.ID <- filter(targets.final, grepl(reps.names, Sample.Name) & Sample.Num != "1r")$RNA.ID

#Fix some PIDNs which changed
replace.PIDNs <- data.frame(Orig = c("6237", "6247", "6236", "6122", "6415", "6418"), New = c("4756", "4757", "4763", "4257", "4759", "206"))
replace.PIDNs.key <- paste(replace.PIDNs$Orig, collapse = "|")
old.PIDNs <- targets.final[grepl(replace.PIDNs.key, targets.final$PIDN),]$PIDN
targets.final[grepl(replace.PIDNs.key, targets.final$PIDN),]$PIDN <- as.character(replace.PIDNs[match(old.PIDNs, replace.PIDNs$Orig),]$New)
targets.final %<>% arrange(PIDN, Date.Drawn)

#Fix samples ranges which skip or do not start at 1.  Temporarily remove replicates because they break detection functions!
targets.noreps <- filter(targets.final, Sample.Num != "1r")
skips <- by(targets.noreps, targets.noreps$PIDN, select, Sample.Num) %>% laply(Compose(as.numeric, sort, diff, find.skips))
baselines <- by(targets.noreps, targets.noreps$PIDN, select, Sample.Num) %>% laply(Compose(as.numeric, find.baseline))
skip.ids <- sort(unique(targets.noreps$PIDN))[skips]
baseline.ids <- sort(unique(targets.noreps$PIDN))[baselines]
all.ids <- unique(c(skip.ids, baseline.ids))

targets.fix <- filter(targets.final, PIDN %in% all.ids & Sample.Num != "1r")
fixed.samplenums <- by(targets.fix, targets.fix$PIDN, select, Sample.Num) %>% llply(function(x){ return(1:length(x)) }) %>% reduce(c) %>% as.character #Use a proper lambda!
targets.final %<>% arrange(PIDN)
targets.final[targets.final$PIDN %in% all.ids & targets.final$Sample.Num != "1r",]$Sample.Num <- fixed.samplenums

subjects.gaa <- select(subjects.reduce, PIDN, GAA1, GAA2) %>% filter(!is.na(GAA1)) %>% filter(!is.na(GAA2))
repeat1 <- apply(select(subjects.gaa, GAA1:GAA2), 1, Compose(as.numeric, min))
repeat2 <- apply(select(subjects.gaa, GAA1:GAA2), 1, Compose(as.numeric, max))
subjects.reduce[!is.na(subjects.reduce$GAA1) & !is.na(subjects.reduce$GAA2), ]$GAA1 <- as.character(repeat1)
subjects.reduce[!is.na(subjects.reduce$GAA1) & !is.na(subjects.reduce$GAA2), ]$GAA2 <- as.character(repeat2)
subjects.reduce$GAA1 %<>% as.numeric
subjects.reduce$GAA2 %<>% as.numeric

hell.conditions <- subjects.reduce$GAA1 < 40 & !is.na(subjects.reduce$GAA1)
subjects.reduce[hell.conditions,]$GAA1 <- NA
subjects.reduce[is.na(subjects.reduce$GAA1),]$GAA1 <- subjects.reduce[is.na(subjects.reduce$GAA1),]$GAA2

#subjects.out <- filter(subjects.all, Site == "CHOP") %>% select(FA_PatientID:dob, onset, gaa1, gaa2)
#colnames(subjects.out) <- c("PIDN", "Sex", "Status", "DOB", "Onset", "GAA1", "GAA2")
#write.xlsx(subjects.out, "subjects.out.xlsx")

targets.final <- join(targets.final, subjects.reduce)
targets.final$Sample.Num %<>% factor
targets.final$Batch %<>% factor
#targets.final$Sex %<>% factor
targets.final$Status %<>% factor
targets.final$Sample.Name %<>% str_replace_all(" ", "")

needs.samplenum <- str_split(targets.final$Sample.Name, "_") %>% llply(length) %>% reduce(c) %>% laply(`==`, 2) %>% which
targets.final[needs.samplenum,]$Sample.Name %<>% paste(targets.final[needs.samplenum,]$Sample.Num, sep = "_")

targets.final$Sex %<>% toupper %>% factor
targets.final$Sample.Name %<>% str_split("_") %>% llply(function(x) { return(x[-length(x)]) }) %>% llply(paste, collapse = "_") %>% reduce(c)
targets.final$Sample.Name %<>% paste(targets.final$Sample.Num, str_sub(targets.final$Status, 1, 3), sep = "_") 
targets.final$Draw.Age <- as.Date(targets.final$Date.Drawn) - as.Date(targets.final$DOB)

save(targets.final, file = "./targets.final.rda")
