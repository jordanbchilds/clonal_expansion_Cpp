wt = read.delim("ipu_wldCount.txt",sep="\t",stringsAsFactors=FALSE,header=FALSE)
mut = read.delim("ipu_mntCount.txt",sep="\t",stringsAsFactors=FALSE,header=FALSE)

# Seems to be a column of NA at the right of these files
imax = 121
Tmax = 120.0

wt = wt[,1:imax]
mut = mut[,1:imax]
t = Tmax*(0:(imax-1))/(imax-1)

cn = wt+mut
mutload = mut/cn

Nmax = 141132
png("MutationLoad.png",type="cairo-png",width=480*5,height=480*5,pointsize=12*5)
plot(NULL,xlim=range(t),ylim=c(0.0,1.0),xlab="Time (y)",ylab="Mutation load")
for(i in 1:Nmax) points(t,mutload[i,],type="l",col=rgb(0,0,0,0.01))
dev.off()

png("CopyNumber.png",type="cairo-png",width=480*5,height=480*5,pointsize=12*5)
plot(NULL,xlim=range(t),ylim=c(0.0,max(cn)),xlab="Time (y)",ylab="Total Copy Number")
for(i in 1:Nmax) points(t,cn[i,],type="l",col=rgb(0,0,0,0.01))
dev.off()

Nbr = 750
br = seq(0,Nbr)/Nbr
mids = (br[2:length(br)]+br[1:(length(br)-1)])/2
res = matrix(nrow=length(t),ncol=length(mids))
for(i in 1:length(t)){
 res[i,] = hist(mutload[,i], br, plot = FALSE)$density
}
png("Rainbow.png",width=2000,height=2000)
 image(res[5:121,],col = rainbow(10000))
dev.off()
pdf("Densities.pdf")
for(i in 1:length(t)){
plot(mids,res[i,],type="l",main=paste("t =",t[i]),xlab="Mutation load",ylab="Density")
}
dev.off()