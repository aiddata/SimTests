#Visualizations for SciClone outputs
#Number of files in directory
f <- list.files("/mnt/sc/geoMatch_testing/Results",pattern="*.csv")
iterations <- length(f)

for(i in 1:iterations)
{
  if(i == 1)
  {
    record.dataframe <- read.csv(paste("/mnt/sc/geoMatch_testing/Results/",f[i],sep=""))
  }
  else
  {
    record.dataframe[i,] <- read.csv(paste("/mnt/sc/geoMatch_testing/Results/",f[i],sep=""))
  }
}





plot(xlim=c(0,3),ylim=c(0,3),record.dataframe$avg.spill, 
     record.dataframe$abs.error.geoMatch, 
     col=rgb(1,0,0,alpha=0.5), pch=3, cex=0.5,
     main="Absolute error in Treatment Estimates",
     ylab="Error",
     xlab="Spillover")

lines(lowess(record.dataframe$avg.spill, 
             record.dataframe$abs.error.geoMatch), 
      col=rgb(1,0,0), pch=3)

lines(lowess(record.dataframe$avg.spill, 
             record.dataframe$abs.error.matchit), 
      col=rgb(0,0,1), pch=4)
points(record.dataframe$avg.spill, 
       record.dataframe$abs.error.matchit, 
       col=rgb(0,0,1,alpha=0.5), pch=4, cex=0.5)

lines(lowess(record.dataframe$avg.spill, 
             record.dataframe$abs.error.lm), 
      col=rgb(0,1,0), pch=4)
points(record.dataframe$avg.spill, 
       record.dataframe$abs.error.lm, 
       col=rgb(0,1,0,alpha=0.5), pch=4, cex=0.5)

legend("topleft",
       cex = 0.65,
       legend=c("geoMatch","Baseline LM", "Baseline MatchIt"), 
       pch=c(pch = 3, pch=4, pch=4),
       col=c(col="red", col="green", col="blue"), title = "Legend")
