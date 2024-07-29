args<-commandArgs(T)

dat = read.table(args[1],sep="\t",header=F)
dat$pVal <- dbinom(dat$V4,(dat$V4+dat$V5), as.numeric(args[2]))
write.table(dat,paste0(args[1],"binomial_test.cov"),row.names=F,sep="\t",quote=F)

system(paste0("gzip ",args[1],"binomial_test.cov"))

