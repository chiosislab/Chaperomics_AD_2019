######################################################
## Loading library
######################################################
library(plyr) 
library(splitstackshape)

######################################################
## Loading data
######################################################
setwd ("~/...") ## Specify the working directory

## Loading BioGrid (H.Sapiens)
## The raw data of BioGrid of can be downloaded from https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.4.160/
BioGrid.raw <- read.csv ("./library/BioGrid/3.4.160/BIOGRID-ORGANISM-Homo_sapiens-3.4.160.tab2.txt", colClasses=as.vector (sapply("character", function (x) rep(x,24))), header=TRUE, sep="\t", stringsAsFactors = F)

# Loading Entrez_Gene to UniProtAC
uniprot <-  read.csv ("./library/uniprot-reviewed_human.0508.2018.slim.tab", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)
mapping_uniprot2GeneID <- concat.split.multiple(data = uniprot, split.cols = "Cross.reference..GeneID.", seps=";", direction = "long") ## Ignore the warning
mapping_uniprot2GeneID <- subset (mapping_uniprot2GeneID, select=c("Entry","Cross.reference..GeneID."))

uniprot.current <- read.csv ("~/Google Drive/Project/AD/libraries/uniprot/uniprot-reviewed_human.1220.2017_extended.tab", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)
mapping_uniprot2GeneID <- merge (mapping_uniprot2GeneID, uniprot.current[,c("Entry","Entry.name")],all.x=TRUE,sort=F)
######################################################
## Mapping BioGrid (adding UniProt ID)
######################################################

BioGrid.UniProtID <- merge (BioGrid.raw, mapping_uniprot2GeneID, by.x="Entrez.Gene.Interactor.A", by.y="Cross.reference..GeneID.", all.x=TRUE,sort=F)
names(BioGrid.UniProtID)[names(BioGrid.UniProtID)=="Entry"] <- "PARTICIPANT_A_Entry"
names(BioGrid.UniProtID)[names(BioGrid.UniProtID)=="Entry.name"] <- "PARTICIPANT_A_Entry.name"

BioGrid.UniProtID <- merge (BioGrid.UniProtID, mapping_uniprot2GeneID, by.x="Entrez.Gene.Interactor.B", by.y="Cross.reference..GeneID.", all.x=TRUE,sort=F)
names(BioGrid.UniProtID)[names(BioGrid.UniProtID)=="Entry"] <- "PARTICIPANT_B_Entry"
names(BioGrid.UniProtID)[names(BioGrid.UniProtID)=="Entry.name"] <- "PARTICIPANT_B_Entry.name"

BioGrid.UniProtID.unmapped <- BioGrid.UniProtID[(is.na(BioGrid.UniProtID$PARTICIPANT_A_Entry) | is.na(BioGrid.UniProtID$PARTICIPANT_B_Entry)),] # Record un-UniProt Entries
BioGrid.UniProtID <- BioGrid.UniProtID[!((is.na(BioGrid.UniProtID$PARTICIPANT_A_Entry) | is.na(BioGrid.UniProtID$PARTICIPANT_B_Entry))),] # Remove un-UniProt Entries

write.csv (BioGrid.UniProtID.unmapped,"BioGrid.UniProtID.unmapped.csv")
BioGrid.slim <- data.frame (PARTICIPANT_A_Entry=BioGrid.UniProtID$PARTICIPANT_A_Entry, PARTICIPANT_B_Entry=BioGrid.UniProtID$PARTICIPANT_B_Entry, PARTICIPANT_A_Entry.name=BioGrid.UniProtID$PARTICIPANT_A_Entry.name, PARTICIPANT_B_Entry.name=BioGrid.UniProtID$PARTICIPANT_B_Entry.name, Source="BioGrid", Experiment=BioGrid.UniProtID$Experimental.System, Literature=BioGrid.UniProtID$Pubmed.ID, Interaction=BioGrid.UniProtID$Experimental.System.Type, stringsAsFactors = F)
#write.csv (BioGrid.slim,"BioGrid.csv")


######################################################
# Loading IntAct
######################################################

## The raw data of IntAct of can be downloaded from ftp://ftp.ebi.ac.uk/pub/databases/intact/2018-05-11
## The original file raw file (all.zip) contains PPIs from all species, we extracted only human PPIs (in Linux command line)
IntAct.raw <- read.csv (unz("./IntAct/05.08.2018/intact_human.zip", "intact_human.txt"), header=F, sep="\t",stringsAsFactors = F)

IntAct.raw2 <- IntAct.raw
IntAct.raw2 <- IntAct.raw2[!(grepl ("intact|ensembl",IntAct.raw2[,1]) | grepl ("intact|ensembl",IntAct.raw2[,2])),] # Remove un-UniProt Entries

ac1 <- gsub ("uniprotkb:","",IntAct.raw2[,1])
ac1 <- gsub ("-[A-Za-z0-9_]*","",ac1)
ac2 <- gsub ("uniprotkb:","",IntAct.raw2[,2])
ac2 <- gsub ("-[A-Za-z0-9_]*","",ac2)

IntAct.raw3 <- data.frame (PARTICIPANT_A_Entry=ac1, PARTICIPANT_B_Entry=ac2, Source="IntAct",Experiment=IntAct.raw2[,7], Literature=IntAct.raw2[,9], Interaction=IntAct.raw2[,12], stringsAsFactors = F)
IntAct.raw3 <- IntAct.raw3 [!IntAct.raw3$PARTICIPANT_A_Entry=="",]
IntAct.raw3 <- IntAct.raw3 [!IntAct.raw3$PARTICIPANT_B_Entry=="",]

IntAct.raw3 <- merge (IntAct.raw3, mapping_uniprot2GeneID, by.x="PARTICIPANT_A_Entry", by.y="Entry", all.x=TRUE)
names(IntAct.raw3)[names(IntAct.raw3)=="Entry.name"] <- "PARTICIPANT_A_Entry.name"
IntAct.raw3 <- merge (IntAct.raw3, mapping_uniprot2GeneID, by.x="PARTICIPANT_B_Entry", by.y="Entry", all.x=TRUE)
names(IntAct.raw3)[names(IntAct.raw3)=="Entry.name"] <- "PARTICIPANT_B_Entry.name"

IntAct.slim <- data.frame (PARTICIPANT_A_Entry=IntAct.raw3$PARTICIPANT_A_Entry, PARTICIPANT_B_Entry=IntAct.raw3$PARTICIPANT_B_Entry, PARTICIPANT_A_Entry.name=IntAct.raw3$PARTICIPANT_A_Entry.name, PARTICIPANT_B_Entry.name=IntAct.raw3$PARTICIPANT_B_Entry.name, Source=IntAct.raw3$Source,Experiment=IntAct.raw3$Experiment, Literature=IntAct.raw3$Literature, Interaction=IntAct.raw3$Interaction, stringsAsFactors = F)
IntAct.slim <- IntAct.slim[!(is.na(IntAct.slim$PARTICIPANT_A_Entry.name) | is.na(IntAct.slim$PARTICIPANT_B_Entry.name)),] # Remove un-UniProt Entries


######################################################
# Merging BioGrid and IntAct
######################################################
interactions <- rbind (BioGrid.slim, IntAct.slim)
unique (interactions$Interaction)

###### Filter #######
interactions <- interactions[!interactions$Interaction == "genetic",] ## Discard "genetic" interactions
interactions <- interactions[!interactions$Interaction == "psi-mi:MI:0208(genetic interaction)",] ## Discard "genetic" interactions
sel.delete.experiment <- c("Co-localization","genetic interference","Synthetic Rescue","Synthetic Growth Defect","Synthetic Lethality")
sel.delete.interaction <- c("colocalization","genetic interaction")
interactions <- interactions [! grepl (paste0(sel.delete.experiment,collapse = "|"),interactions$Experiment), ]
interactions <- interactions [! grepl (paste0 (sel.delete.interaction, collapse="|"), interactions$Interaction),]

ppi <- interactions
str (ppi)
save (ppi,file="./library/ppi_2018_0508.Rdata")
