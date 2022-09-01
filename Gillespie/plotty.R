wld_raw = read.table("./OUTPUT/rand_ICNe3/ipu_wldCount.txt")
mnt_raw = read.table("./OUTPUT/rand_ICNe3/ipu_mntCount.txt")

imax = 121
Tmax = 120.0

wld_ts = wld_raw[,1:imax]
mnt_ts = mnt_raw[,1:imax]

t = Tmax*(0:(imax-1))/(imax-1)

copy_num = wld_ts+mnt_ts
mut_load = mnt_ts/copy_num

Nmax = 141132

png("./FIGURE/rand_ICNe3.png",width=2*480*5,height=480*5,pointsize=12*5)
{
  par(mfrow=c(1,2), mar=c(6,6,3,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  
  plot(NULL,xlim=range(t),ylim=c(0.0,1.0),
       xlab="Time (years)",ylab="Mutation load")
  for(i in 1:Nmax){
    points(t,mut_load[i,],type="l",col=rgb(0,0,0,0.01))
  }
  
  plot(NULL,xlim=range(t),ylim=c(0.0,2000.0),
       xlab="Time (years)",ylab="Copy Number")
  for(i in 1:Nmax) points(t,copy_num[i,],type="l",col=rgb(0,0,0,0.01))

}
dev.off()







