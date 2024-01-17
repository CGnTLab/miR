########################
### Batch Correction ###
########################

library(sva)

df<-read.csv("input.csv")
batch<-df$Batch
batchCorrected<-ComBat(df,batch = batch, mod = NULL)
mat<-as.data.frame(batchCorrected)
row.names(mat1)<-row
write.csv(mat1,"batchCorrected_file.csv")
