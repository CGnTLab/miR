########################
### Batch Correction ###
########################
## Batch correction in RNA-Seq data involves mitigating technical variations or batch effects
## introduced during sample processing. Methods like ComBat or surrogate variable analysis are 
## employed to adjust gene expression values, ensuring more accurate comparisons and reliable 
## downstream analyses by minimizing batch-related variability.
library(sva)

df<-read.csv("input.csv")
batch<-df$Batch
batchCorrected<-ComBat(df,batch = batch, mod = NULL)
mat<-as.data.frame(batchCorrected)
row.names(mat1)<-row
write.csv(mat1,"batchCorrected_file.csv")
