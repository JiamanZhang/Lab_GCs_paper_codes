library(stats)
library(base)
args<-commandArgs(TRUE)
inp<-args[1]
gcount<-args[2]
chr<-args[3]
oup<-args[4]
#setwd("/Lustre01/tangqianzi/data/HiCmerge/results/GOM/")
#data<-read.table(gzfile("oe_matrix.GOM.100000.15.txt.gz"),sep="\t",header=F)
data<-read.table(inp,sep="\t",header=F)
mymatrix<-data

result<-prcomp(t(mymatrix),center=TRUE)
x1<-as.vector(result$x[,1])

data1<-read.table(file=gcount,header=F)
newdata<-data1[data1$V1==chr,4]

if (cor(newdata,x1)>0){
myPC1<-x1
print (cor(newdata,x1))
} else {
myPC1<-(-x1)
print (cor(newdata,x1))
}

result<-cbind(data1[data1$V1==chr,c(1:3)],myPC1)
write.table(result,file=oup,sep="\t",quote=F,row.names=F,col.names=F)
