################################################
#### Tested with R 3.6.1 under OS X 10.13.6 on an iMac (Late 2015)
#### R packages:
#### 'clusterProfiler','org.Hs.eg.db'
################################################

################################################
#### Loading PPI databases (BioGrid+Intact) ####
################################################
setwd ("...") ## Specify the working directory
load (file="./library/ppi_2018_0508.Rdata")

################################################
#### Loading Chaperome ID                   ####
################################################
load (file="./library/list.chaperome.brehme.uniprot.Rdata")

################################################
#### Loading Chaperomics Data (with DE analyses)
################################################

load ("./input/ms1.exp1-3.uniprot.results.Rdata") ## Human AD
data.human.ad <- ms1.uniprot.results.all
remove (ms1.uniprot.results.all)

## Human AD up/down (AD>ND / AD<ND)
input <- data.human.ad
data.human.ad.up   <- input [input$Raw.p <= 0.25 & input$Raw.FC > 1,]
data.human.ad.down <- input [input$Raw.p <= 0.25 & input$Raw.FC < 1,]
remove (input)

uniprot <- read.csv ("./library/uniprot-reviewed_human.1220.2017_extended.tab", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)

################################################
#### Converting UniProtAC to EntrezID
################################################
#library(UniProt.ws)
#up <- UniProt.ws(taxId=9606)
#uniprot2entrezID <- select(up, uniprot$Entry, c("ENTRY-NAME","ENTREZ_GENE"),"UNIPROTKB")
#write.csv (uniprot2entrezID, "~/Google Drive/R/UniProt/uniprot2entrezID_2.csv")
#uniprot2entrezID <- read.csv("~/Google Drive/R/UniProt/uniprot2entrezID.csv")

#uniprot2entrezID <- read.csv ("~/Google Drive/R/UniProt/2018_05_08/uniprot-reviewed_human.0508.2018.slim.tab", quote = "", header=TRUE, sep="\t", stringsAsFactors = F)
#uniprot2entrezID <- subset (uniprot2entrezID, select=c(Entry, Cross.reference..GeneID.))
#names (uniprot2entrezID) <- c("UNIPROTKB","ENTREZ_GENE")
#library (splitstackshape)
#uniprot2entrezID <- concat.split.multiple(uniprot2entrezID, split.cols = "ENTREZ_GENE",seps=";",direction="long")
#write.csv (uniprot2entrezID, "~/Google Drive/R/UniProt/uniprot2entrezID.csv",row.names = F,quote=F)

uniprot2entrezID <- read.csv("./library/uniprot2entrezID.csv")

attach.entregene <- function (x,uniprotkb_col) {
  output <- merge (x, uniprot2entrezID, by.x=uniprotkb_col,by.y="UNIPROTKB",all.x=TRUE,sort=F)
  output <- output [!is.na (output$ENTREZ_GENE),]
  return (output)
}

uniprot.entry.name2entreID <- function (x) {
  input <- data.frame (ID=x,stringsAsFactors = F)
  output <- merge (input, uniprot[,c("Entry.name","Entry")],by.x="ID",by.y="Entry.name",all.x=TRUE,sort=F)
  output <- merge (output, uniprot2entrezID, by.x="Entry",by.y="UNIPROTKB",all.x=TRUE,sort=F)
  return (output$ENTREZ_GENE)
}

################################################
#### Preparing input file
################################################

input <- data.human.ad
list.chaperome.input <- input$Entry.name [input$Entry %in% list.chaperome.brehme.uniprot$UNIPROTKB]

################################################
#### Preapring iGSEA (GO Enrichment module)  
################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")
#biocLite("org.Hs.eg.db")

library ("clusterProfiler")
library ("org.Hs.eg.db")

iGSEA <- function (input, list.chaperome.input, ppi, pAdjustMethod, ont, pvalueCutoff, outputfile) {
             
               for (i in 1:length (list.chaperome.input)) {
                if (i==1) {df.total<-NULL}
                inp <- list.chaperome.input[i]
                
                ppi.selected.1 <- ppi[ppi$PARTICIPANT_A_Entry.name %in% inp,]
                ppi.selected.2 <- ppi[ppi$PARTICIPANT_B_Entry.name %in% inp,]
                
                interactors <- c(ppi.selected.1$PARTICIPANT_B_Entry.name, ppi.selected.2$PARTICIPANT_A_Entry.name)
                interactors <- interactors [interactors %in% input$Entry.name]
                
                interactors.entrezID <- uniprot.entry.name2entreID(interactors)
                interactors.entrezID <- interactors.entrezID[!is.na (interactors.entrezID)]
                
                print (paste0("Chaperome_ID:",inp," ",i,"/",length (list.chaperome.input)))
                if (length(interactors.entrezID)>0) {
                  go.enrichment <- enrichGO (interactors.entrezID, 'org.Hs.eg.db',keyType = "ENTREZID", pAdjustMethod = "BH", ont="BP",readable = T,pvalueCutoff = 1)
                  df <- as.data.frame(go.enrichment)
                  if (nrow (df)>0) {df$Chaperome_ID <- inp}
                  if (is.null(df.total)) {df.total <- df}
                  else {df.total <- rbind (df.total, df)}
                  print (df)
                }
                
               }
              return (df.total)
              write.table (df.total, outputfile, row.names = F, quote=F,sep="\t")
}

input <- data.human.ad.up
input <- attach.entregene(input, "Entry")
df.total.up <- iGSEA (input, list.chaperome.input, ppi, "BH", "BP", 1, "df.total.up.txt") ## Run iGSEA with the DE proteins (AD>ND)

input <- data.human.ad.down
input <- attach.entregene(input, "Entry")
df.total.down <- iGSEA (input, list.chaperome.input, ppi, "BH", "BP", 1, "df.total.down.txt")  ## Run iGSEA with the DE proteins (AD<ND)

save (df.total.up, df.total.down,file="./output/iGSEA.results.Rdata")

## Combine data together
p.cutoff <- 0.001
df.total.up.full <- df.total.up [(df.total.up$p.adjust < p.cutoff),]
df.total.down.full <- df.total.down [(df.total.down$p.adjust < p.cutoff),]
log10(df.total.up.full$p.adjust) -> df.total.up.full$log10.p.adjust
log10(df.total.down.full$p.adjust) -> df.total.down.full$log10.p.adjust
df.total.up.full$Status <- "Up"
df.total.down.full$Status <- "Down"
df.total.up.down.full <- rbind (df.total.up.full, df.total.down.full)

## Write output
write.table (df.total.up.down.full, "./output/df.total.up.down.full.txt",sep="\t",row.names = F, quote=F) ## Full raw output of iGSEA

################################################
#### SynGO
#### Synpase related gene set was prepared based on the SynGO database (updated till Jan 2018)
#### (The Synapse Gene Ontology and Annotation Initiative, http://www.geneontology.org/page/syngo-synapse-biology) 
#### with exclusion of three overly vague GO terms: “protein binding”, “plasma membrane” and “integral component of plasma membrane”. 
#### Related SynGO terms were included in Table S4
################################################

go.SynaptomeDB <- read.csv("./library/synaptomeDB.txt",sep="\t",header=FALSE)
load (file="./output/iGSEA.results.Rdata")
p.cutoff <- 0.1

df.total.up.synGO <- df.total.up [(df.total.up$ID %in% go.SynaptomeDB$V5) & (df.total.up$p.adjust < p.cutoff),]
df.total.down.synGO <- df.total.down [(df.total.down$ID %in% go.SynaptomeDB$V5) & (df.total.down$p.adjust < p.cutoff),]

log10(df.total.up.synGO$p.adjust) -> df.total.up.synGO$log10.p.adjust
log10(df.total.down.synGO$p.adjust) -> df.total.down.synGO$log10.p.adjust

df.total.up.synGO$Status <- "Up"
df.total.down.synGO$Status <- "Down"
df.total.up.down.synGO <- rbind (df.total.up.synGO, df.total.down.synGO)

write.table (df.total.up.synGO, "./output/df.total.up.syngo.txt",sep="\t",row.names = F, quote=F)
write.table (df.total.down.synGO, "./output/df.total.down.syngo.txt",sep="\t",row.names = F, quote=F)
write.table (df.total.up.down.synGO, "./output/df.total.up.down.synGO.txt",sep="\t",row.names = F, quote=F) ## iGSEA results with only SynGO related GO termes, this table was further imported into Cytoscape to generate the iGSEA figure (See Method, section "Interactomic gene-set enrichment analysis (iGESA)")
