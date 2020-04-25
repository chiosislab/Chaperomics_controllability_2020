load (file="~/Google Drive/Project/Joshi/MiaPaca_vs_468_2019/3-iGSEA/LFQ.A/iGSEA.up.results.Rdata")
load (file="~/Google Drive/Project/Joshi/MiaPaca_vs_468_2019/3-iGSEA/LFQ.A/iGSEA.down.results.Rdata")
iGSEA.up$Status <- "Up.LFQ.A"
iGSEA.down$Status <- "Down.LFQ.A"
iGSEA.LFQ.A.total <- rbind (iGSEA.up, iGSEA.down)

load (file="~/Google Drive/Project/Joshi/MiaPaca_vs_468_2019/3-iGSEA/LFQ.B/iGSEA.up.results.Rdata")
load (file="~/Google Drive/Project/Joshi/MiaPaca_vs_468_2019/3-iGSEA/LFQ.B/iGSEA.down.results.Rdata")
iGSEA.up$Status <- "Up.LFQ.B"
iGSEA.down$Status <- "Down.LFQ.B"
iGSEA.LFQ.B.total <- rbind (iGSEA.up, iGSEA.down)

iGSEA.LFQ.A.B.total <- rbind (iGSEA.LFQ.A.total, iGSEA.LFQ.B.total)

p.cutoff <- 0.001

log10 (iGSEA.LFQ.A.B.total$p.adjust) -> iGSEA.LFQ.A.B.total$log10.p.adjust

iGSEA.selected <- iGSEA.LFQ.A.B.total [iGSEA.LFQ.A.B.total$p.adjust <= p.cutoff, ]
df.sidetable <- data.frame(name=unique(iGSEA.selected$Chaperome_ID), Description=unique (iGSEA.selected$Chaperome_ID))

setwd(".")
write.table (iGSEA.selected, "iGSEA.selected.combined.txt",sep="\t",row.names = F, quote=F)
write.table (df.sidetable, "sidetable.combined.txt",sep="\t",row.names = F, quote=F)

