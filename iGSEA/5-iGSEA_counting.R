iGSEA <- read.csv("iGSEA.selected.combined.txt",sep = "\t", stringsAsFactors = F)

chaperome <- unique (iGSEA$Chaperome_ID)

f <- NULL
for (i in 1:length (chaperome)) {
  df <- iGSEA [iGSEA$Chaperome_ID==chaperome[i],]
  t <- data.frame (table (df$Status))
  a <- data.frame (Chaperome_ID=chaperome[i], t)
  f <- rbind (f, a)
}
library(dplyr)

mat <- reshape(f,timevar = c("Var1"),idvar = c("Chaperome_ID"), direction="wide")
mat [is.na(mat)] <- 0

row.names (mat) <- mat$Chaperome_ID
mat <- subset (mat,select=-Chaperome_ID)  
mat.LFQ.B <- mat [,c("Freq.Up.LFQ.B", "Freq.Down.LFQ.B")]
mat.LFQ.A <- mat [,c("Freq.Up.LFQ.A", "Freq.Down.LFQ.A")]
mat.LFQ.combined <- mat [,c("Freq.Up.LFQ.B", "Freq.Down.LFQ.B", "Freq.Up.LFQ.A", "Freq.Down.LFQ.A")]
write.csv (mat.LFQ.combined, "iGSEA-Counting.csv")

library (ComplexHeatmap)
library (circlize)
inp <- as.matrix (mat.LFQ.combined)

pdf ("iGSEA-Counting.pdf", width = 3.5, height=0.15*nrow (inp), useDingbats = FALSE)
Heatmap(inp, col = colorRamp2(c(min(inp),5,50, max (inp)), c("white", "lightblue","yellow", "red")), cluster_columns = FALSE)
dev.off()

colSums(mat.LFQ.combined)