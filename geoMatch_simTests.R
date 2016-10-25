library(devtools)
install_local('/home/aiddata/Desktop/Github/geoMatch')
library(geoMatch)

source("spatial.simulation.R")

iterations <- 4

record.dataframe <- data.frame(n.points=round(runif(iterations,min=500,max=500)),
                               cov.Decay = runif(iterations, min=25, max=25),
                               se.Decay = runif(iterations, min=1, max=1),
                               t.Decay = runif(iterations, min=0, max=50),
                               sim.T.e = runif(iterations, min=0.0, max=0.0),
                               T.percent = runif(iterations, min=0.10, max=0.10),
                               sim.Y.scale = runif(iterations, min=0.0, max=0.0),
                               Theta = runif(iterations, min=1.0, max=1.0),
                               sim.Y.cov.effect = runif(iterations, min=0.0, max=0.0),
                               sim.Y.het.effect = runif(iterations, min=0.0, max=0.0),
                               sim.Y.e = runif(iterations, min=0.0, max=0.0),
                               spill.mag = runif(iterations, min=0.0, max=1.0),
                               avg.spill = NA,
                               est.T.lm = NA,
                               est.T.matchit = NA,
                               est.T.geoMatch = NA,
                               abs.error.lm = NA,
                               abs.error.matchit = NA,
                               abs.error.geoMatch = NA,
                               mean.error.lm = NA,
                               mean.error.matchit = NA,
                               mean.error.geoMatch = NA,
                               true.T.ie = NA)
for(i in 1:iterations)
  {
    print(i)
    it.spdf <- spatial.simulation(n.points = record.dataframe["n.points"][i,], 
                        cov.Decay = record.dataframe["cov.Decay"][i,], 
                        se.Decay = record.dataframe["se.Decay"][i,],
                        t.Decay = record.dataframe["t.Decay"][i,],
                        sim.T.e= record.dataframe["sim.T.e"][i,], 
                        T.percent = record.dataframe["T.percent"][i,],
                        sim.Y.scale = record.dataframe["sim.Y.scale"][i,],
                        Theta = record.dataframe["Theta"][i,],
                        sim.Y.cov.effect = record.dataframe["sim.Y.cov.effect"][i,],
                        sim.Y.het.effect = record.dataframe["sim.Y.het.effect"][i,],
                        sim.Y.e = record.dataframe["sim.Y.e"][i,],
                        spill.mag = record.dataframe["spill.mag"][i,])
    it.geoMatch <- geoMatch(T ~ X, 
             method = "nearest", 
             caliper=0.25, 
             data = it.spdf, 
             outcome.variable="Y", 
             outcome.suffix="_adjusted")
    
    it.Match <- matchit(T ~ X, 
                        method = "nearest", 
                        caliper=0.25, 
                        data = it.spdf@data)
    
    lm.model <- lm(Y ~ T + X, data=it.spdf@data)
    matchit.model <- lm(Y ~ T + X, data=match.data(it.Match))
    geo.model <- lm(Y_adjusted ~ T + X, data=match.data(it.geoMatch))

    record.dataframe["mean.error.lm"][i,] <- mean((lm.model$coefficients["T"] - it.spdf@data[it.spdf@data$T == 1,]$ie.spill))
    record.dataframe["mean.error.matchit"][i,] <- mean((matchit.model$coefficients["T"] - it.spdf@data[it.spdf@data$T == 1,]$ie.spill))
    record.dataframe["mean.error.geoMatch"][i,] <- mean((geo.model$coefficients["T"] - it.spdf@data[it.spdf@data$T == 1,]$ie.spill))
    
    record.dataframe["abs.error.lm"][i,] <- mean(abs((lm.model$coefficients["T"] - it.spdf@data[it.spdf@data$T == 1,]$ie.spill)))
    record.dataframe["abs.error.matchit"][i,] <- mean(abs((matchit.model$coefficients["T"] - it.spdf@data[it.spdf@data$T == 1,]$ie.spill)))
    record.dataframe["abs.error.geoMatch"][i,] <- mean(abs((geo.model$coefficients["T"] - it.spdf@data[it.spdf@data$T == 1,]$ie.spill)))
    
    
    record.dataframe["est.T.lm"][i,] <- lm.model$coefficients["T"]
    record.dataframe["est.T.matchit"][i,] <- matchit.model$coefficients["T"]
    record.dataframe["est.T.geoMatch"][i,] <- geo.model$coefficients["T"]
    record.dataframe["true.T.ie"][i,] <- mean(it.spdf@data[it.spdf@data$T == 1,]["ie.spill"][[1]])
    record.dataframe["avg.spill"][i,] <- mean(it.spdf@data[it.spdf@data$T == 0,]["t.spill"][[1]])
}





plot(ylim=c(-10,10),record.dataframe$avg.spill, 
     record.dataframe$abs.error.geoMatch, 
     col=rgb(1,0,0,alpha=0.5), pch=3, cex=0.5,
     main="Relative absolute error in Treatment Estimates",
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

legend("bottomright",
       cex = 0.65,
       legend=c("geoMatch","Baseline LM", "Baseline MatchIt"), 
       pch=c(pch = 3, pch=4, pch=4),
       col=c(col="red", col="green", col="blue"), title = "Legend")
