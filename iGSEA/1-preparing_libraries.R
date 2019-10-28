##############################################################
## This code prepares some library files (Chaperome list, UniProt, UniProt2EntrezGeneID conversion table) for iGSEA
##############################################################

setwd ("~/Google Drive/Project/AD/code submission/7-iGSEA/") ## Put the current working directory in which the codes were placed

## The chaperome members were prepared based on the work from Breheme M. et al, Cell Rep, 2014 (PUBMED: 25437566), Table S2 (1-s2.0-S2211124714008250-mmc3.xls)
list.chaperome.brehme <- read.csv("./library/chaperome_list.csv", stringsAsFactors = F)

## The UniProt library were download on 12/20/2017 form UniProt.org
uniprot.human <- read.csv ("./library/uniprot-reviewed_human.1220.2017_extended.tab", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)

## Prepare a mapping table between UniProt ID to EntrezGene ID
library(UniProt.ws)
up <- UniProt.ws(taxId=9606)
uniprot2entrezID <- select(up, uniprot$Entry, "ENTREZ_GENE")
write.csv (uniprot2entrezID, "./library/uniprot2entrezID.csv")
    
    ## Add some manual mappings
    uniprot2entrezID <- read.csv("./library/uniprot2entrezID.csv")
    uniprot2entrezID.complementary <- read.csv("./library/uniprot2entrezID_complementary.csv")
    uniprot2entrezID <- rbind (uniprot2entrezID, uniprot2entrezID.complementary)

list.chaperome.brehme.uniprot <- merge (list.chaperome.brehme, uniprot2entrezID, by.x="EntrezID",by.y="ENTREZ_GENE",all.x=TRUE,sort=F)
list.chaperome.brehme.uniprot <- list.chaperome.brehme.uniprot[!is.na(list.chaperome.brehme.uniprot$UNIPROTKB),]
list.chaperome.brehme.uniprot <- merge (list.chaperome.brehme.uniprot, uniprot.human[,c("Entry","Entry.name")], by.x="UNIPROTKB",by.y="Entry",all.x=TRUE,sort=F)
write.csv (list.chaperome.brehme.uniprot, "./library/list.chaperome.brehme.uniprot.csv")
save (list.chaperome.brehme.uniprot, file="./library/list.chaperome.brehme.uniprot.Rdata")
