
myBlack = function(alpha) { rgb(0,0,0,alpha)}

Nsim = 1e3
ml = read.table("./syn_ml.txt", sep="\n")[[1]]
cn = read.table("./syn_cn.txt", sep="\n")[[1]]

ts = 0:120

pdf("./synData.pdf", width=12, height=9)
{
  par(mfrow=c(1,2), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  
  plot(NA, xlim=c(0,120), ylim=c(0,1), 
       main="Mutation Load", xlab="Time (years)", ylab="")
  for(i in 0:120){
    points(i+rnorm(Nsim), ml[i*1000+(1:1000)], pch=20, col=myBlack(0.01))
  }

  plot(NA, xlim=c(0,120), ylim=c(0,2000), 
       main="Copy Number", xlab="Time (years)", ylab="")
  for(i in 0:120){
    points(i+rnorm(Nsim), cn[i*1000+(1:1000)], pch=20, col=myBlack(0.01))
  }
}
dev.off()
