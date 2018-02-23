#!/usr/bin/Rscript
args<-commandArgs(TRUE)
x = read.csv(args, header = TRUE, sep = " ", quote = "\"", dec = ".", fill = TRUE)
pdf(paste(args, ".pdf"))
y <- as.numeric(unlist(x[1]))      # turn into numbers
hist(y,breaks=100)
dev.off()



# x = read.csv("dist.D.Sorted.MovingAvU.bed.CTCFDS11167ex7500", header = TRUE, sep = " ", quote = "\"", dec = ".", fill = TRUE)
# yy <- as.numeric(unlist(y[1]))
#hist(xx,breaks=100,col=rgb(1,0,0,0.5))
#> hist(zz,breaks=100,col=rgb(0,0,1,1),add=T)
#> hist(yy,breaks=100,col=rgb(0,1,0,1),add=T)

