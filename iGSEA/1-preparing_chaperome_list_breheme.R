setwd ("./library/")
list.chaperome.brehme <- read.csv("1-chaperome_list.csv", stringsAsFactors = F)
uniprot.human <- read.csv ("uniprot-reviewed_human.1220.2017_extended.tab", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)
#uniprot.mouse <- read.csv ("~/Google Drive/Project/AD/libraries/uniprot/uniprot-reviewed_mouse.1220.2017_extended.tab", quote="", header=TRUE, sep="\t", stringsAsFactors = FALSE)

## PREPARE A UniProt2EntrezID table

#library(UniProt.ws)
#up <- UniProt.ws(taxId=9606)
#uniprot2entrezID <- select(up, uniprot$Entry, "ENTREZ_GENE")
#write.csv (uniprot2entrezID, "~/Google Drive/R/UniProt/uniprot2entrezID.csv")

uniprot2entrezID <- read.csv("uniprot2entrezID.csv",row.names = "X")
uniprot2entrezID.complementary <- read.csv("uniprot2entrezID_complementary.csv")
uniprot2entrezID <- rbind (uniprot2entrezID, uniprot2entrezID.complementary)
list.chaperome.brehme.uniprot <- merge (list.chaperome.brehme, uniprot2entrezID, by.x="EntrezID",by.y="ENTREZ_GENE",all.x=TRUE,sort=F)
list.chaperome.brehme.uniprot <- list.chaperome.brehme.uniprot[!is.na(list.chaperome.brehme.uniprot$UNIPROTKB),]
list.chaperome.brehme.uniprot <- merge (list.chaperome.brehme.uniprot, uniprot.human[,c("Entry","Entry.name")], by.x="UNIPROTKB",by.y="Entry",all.x=TRUE,sort=F)
write.csv (list.chaperome.brehme.uniprot, "./library/list.chaperome.brehme.uniprot.csv")
save (list.chaperome.brehme.uniprot, file="./library/list.chaperome.brehme.uniprot.Rdata")
