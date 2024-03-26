args<-commandArgs(T)
dat = read.table(args[1],sep="\t",header=F)

pois_lambda = mean(dat$V4)

dat$V5 <- sapply(dat$V4,function(x)ppois(x,pois_lambda) )
dat$V6 = 1-dat$V5
write.table(dat,paste0(args[1],"poisson_test.txt"),row.names=F,col.names=F,quote=F)
