################################################
#### Loading PPI databases (BioGrid+Intact) ####
################################################
load (file="./library/ppi_2018_0508.Rdata")
#write.csv (ppi, "ppi.csv", row.names = F, quote=F)
################################################
#### Loading Chaperome ID                   ####
################################################
load (file="./library/list.chaperome.brehme.uniprot.Rdata")

################################################
#### Loading Data                           ####
################################################
load ("./proteomics_input/ms1.uniprot.results.Rdata")
input <- ms1.uniprot.results.all
remove (ms1.uniprot.results.all)
## Human AD up/down
data.up   <- input [input$LFQ.B.p <= 0.1 & input$LFQ.B.FC > 1,]
data.down <- input [input$LFQ.B.p <= 0.1 & input$LFQ.B.FC < 1,]

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

list.chaperome.input <- input$Entry.name [input$Entry %in% list.chaperome.brehme.uniprot$UNIPROTKB]

################################################
#### Preapring iGSEA (GO Enrichment module)  
################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")
#biocLite("org.Hs.eg.db")

library ("clusterProfiler")
library ("org.Hs.eg.db")


iGSEA.v2 <- function (input, list.chaperome.input, ppi, pAdjustMethod, ont, pvalueCutoff, outputfile) {
             
               for (i in 1:length (list.chaperome.input)) {
                if (i==1) {df.total<-NULL
                           interactors <- NULL
                           df <- NULL}
                inp <- list.chaperome.input[i]
                
                ppi.selected.1 <- ppi[ppi$PARTICIPANT_A_Entry.name %in% inp,]
                ppi.selected.2 <- ppi[ppi$PARTICIPANT_B_Entry.name %in% inp,]
                
                interactors <- c(ppi.selected.1$PARTICIPANT_B_Entry, ppi.selected.2$PARTICIPANT_A_Entry) # Switch from UniProt Entry name to Entry
                interactors <- interactors [interactors %in% input]
                
                #interactors.entrezID <- uniprot.entry.name2entreID(interactors)
                #interactors.entrezID <- interactors.entrezID[!is.na (interactors.entrezID)]
                
                print (paste0("Chaperome_ID:",inp," ",i,"/",length (list.chaperome.input)))
                print (paste0("Interactors_ID:",interactors))
                if (length(interactors)>0) {
                  go.enrichment <- enrichGO (interactors, 'org.Hs.eg.db',keyType = "UNIPROT", pAdjustMethod = "BH", ont="BP",readable = T,pvalueCutoff = 0.1)
                  ## Note: Depending on the version of ClusterProfiler, the argment keytype could be either "keyType" or "keytype". "KeyType" was used on iMac-29''
                  df <- as.data.frame(go.enrichment)
                  if (nrow (df)>0) {df$Chaperome_ID <- inp}
                  if (is.null(df.total)) {df.total <- df}
                  else {df.total <- rbind (df.total, df)}
                }
                print (df)
               }
              return (df.total)
              write.table (df.total, outputfile, row.names = F, quote=F,sep="\t")
}

setwd("./LFQ.B")
input.up <- as.character(data.up$Entry)
iGSEA.up <- iGSEA.v2 (input.up, list.chaperome.input, ppi, "BH", "BP", 1, "df.total.up.txt")
save (iGSEA.up,file="iGSEA.up.results.Rdata")

input.down <- as.character(data.down$Entry)
iGSEA.down <- iGSEA.v2 (input.down, list.chaperome.input, ppi, "BH", "BP", 1, "df.total.down.txt")
save (iGSEA.down,file="iGSEA.down.results.Rdata")

load (file="iGSEA.down.results.Rdata")

save (iGSEA.up, iGSEA.down,file="iGSEA.results.Rdata")


iGSEA.up$Status <- "Up"
iGSEA.down$Status <- "Down"

iGSEA.total <- rbind (iGSEA.up, iGSEA.down)

p.cutoff <- 0.001

log10 (iGSEA.total$p.adjust) -> iGSEA.total$log10.p.adjust

iGSEA.selected <- iGSEA.total [iGSEA.total$p.adjust <= p.cutoff, ]
df.sidetable <- data.frame(name=unique(iGSEA.selected$Chaperome_ID), Description=unique (iGSEA.selected$Chaperome_ID))

write.table (iGSEA.selected, "iGSEA.selected.txt",sep="\t",row.names = F, quote=F)
write.table (df.sidetable, "sidetable.txt",sep="\t",row.names = F, quote=F)
