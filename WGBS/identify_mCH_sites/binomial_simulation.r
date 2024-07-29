args<-commandArgs(T)

dat = read.table(args[1],sep="\t")
dat$cov=dat$V4+dat$V5
dat$simulation = rbinom(size=dat$cov, n=nrow(dat), prob=as.numeric(args[2]))
dat$un = dat$cov - dat$simulation
dat = dat[, c("V1","V2","V3","simulation","un","V6","V7")]
write.table(dat,paste0(args[1],".total.simu.cov"),quote=F,col.names=F,sep="\t",row.names=F)
