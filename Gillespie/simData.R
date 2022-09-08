# gen some fake data
myBlack = function(alpha){ rgb(0,0,0, alpha) }
Ndata = 1e4
times = c(25,55,65)

mut_load = matrix(NA, nrow=Ndata, ncol=length(times))
mut_load[,1] = rnorm(Ndata, 0.4, 0.03)
mut_load[,2] = rnorm(Ndata, 0.8, 0.03)
mut_load[,3] = rnorm(Ndata, 0.9, 0.02)

copy_num = matrix(NA, nrow=Ndata, ncol=length(times))
copy_num[,1] = rnorm(Ndata, 1000, 30)
copy_num[,2] = rnorm(Ndata, 950, 50)
copy_num[,3] = rnorm(Ndata, 1050, 60)

pdf("./FIGURE/simulated_data.pdf", width=12, height=7)
{
  par(mfrow=c(1,2), mar=c(6,6,6,3), cex.main=2, cex.lab=2, cex.axis=1.5)
  plot(NA, xlim=c(0,120), ylim=c(0,1), 
       main="", xlab="Time (years)", ylab="Mutation Load")
  for(i in 1:length(times)){
    points(times[i]+rnorm(Ndata,0,0.5), mut_load[,i], pch=20, cex=0.7, col=myBlack(0.02))
  }
  
  plot(NA, xlim=c(0,120), ylim=c(0,2000), 
       main="", xlab="Time (years)", ylab="Copy Number")
  for(i in 1:length(times)){
    points(times[i]+rnorm(Ndata,0,0.5), copy_num[,i], pch=20, cex=0.7, col=myBlack(0.02))
  }
title(main="Simulated Data", line=-2, outer=TRUE)
}
dev.off()

dir.create("../simulated_data", showWarnings = FALSE)

write.csv(mut_load, file="../simulated_data/mut_load.csv", row.names=FALSE)
write.csv(copy_num, file="../simulated_data/sopy_num.csv", row.names=FALSE)


# priors
curve(dnorm(x,2.64e-3,5e-4), from=1e-3, to=9e-3, 
      ylab="Density", xlab="Reaction Rate", main="My Prior")
curve(dnorm(x,2e-3,5e-4), from=1e-3, to=9e-3, 
      ylab="Density", xlab="Control Parameter", main="My Prior")
hist(rnorm(1e5, 1000,100), breaks=500:1500, col=myBlack(0.05),border=myBlack(0.05),
     freq=FALSE, xlab="Initial Copy Number", main="My Prior")
curve(dunif(x,0.2,0.6), from=0, to=1, 
      ylab="Density", xlab="Initial Mutation Load", main="My Prior")












