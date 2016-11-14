#SciClone Edition
#library(devtools)
#install_github('itpir/geoMatch')
library(geoMatch)

source("/sciclone/home00/geogdan/geoMatch_testing/SimTests/spatial.simulation.R")

#For SciClone - LEAVE THIS AT 1!  Iterations are handled in python.
iterations <- 1

record.dataframe <- data.frame(n.points=round(runif(iterations,min=200,max=200)),
                               cov.Decay = runif(iterations, min=0.01, max=0.01),
                               se.Decay = runif(iterations, min=0.01, max=0.01),
                               t.Decay = runif(iterations, min=50, max=50),
                               sim.T.e = runif(iterations, min=0.0, max=0.0),
                               T.percent = runif(iterations, min=0.5, max=0.5),
                               sim.Y.scale = runif(iterations, min=0.0, max=0.0),
                               Theta = runif(iterations, min=1.0, max=1.0),
                               sim.Y.cov.effect = runif(iterations, min=0.0, max=0.0),
                               sim.Y.het.effect = runif(iterations, min=0.0, max=0.0),
                               sim.Y.e = runif(iterations, min=0.0, max=0.0),
                               spill.mag = runif(iterations, min=0.0, max=2.0),
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
                               true.T.ie = NA,
                               converge = NA)
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
             outcome.suffix="_adjusted",
             max.it=10,
             report = TRUE)
    
    record.dataframe$converge <- it.geoMatch$optim

    
    it.Match <- matchit(T ~ X, 
                        method = "nearest", 
                        caliper=0.25, 
                        data = it.spdf@data)
    
    lm.model <- lm(Y ~ T + X, data=it.spdf@data)
    matchit.model <- lm(Y ~ T + X, data=match.data(it.Match))
    geo.model <- lm(Y_adjusted ~ T + X, data=match.data(it.geoMatch))

    record.dataframe["mean.error.lm"][i,] <- mean((lm.model$coefficients["T"] - mean(it.spdf@data[it.spdf@data$T == 1,]["ie.spill"][[1]])))
    record.dataframe["mean.error.matchit"][i,] <- mean((matchit.model$coefficients["T"] - mean(it.spdf@data[it.spdf@data$T == 1,]["ie.spill"][[1]])))
    record.dataframe["mean.error.geoMatch"][i,] <- mean((geo.model$coefficients["T"] - mean(it.spdf@data[it.spdf@data$T == 1,]["ie.spill"][[1]])))
    
    record.dataframe["abs.error.lm"][i,] <- mean(abs((lm.model$coefficients["T"] - it.spdf@data[it.spdf@data$T == 1,]["ie.spill"][[1]])))
    record.dataframe["abs.error.matchit"][i,] <- mean(abs((matchit.model$coefficients["T"] - it.spdf@data[it.spdf@data$T == 1,]["ie.spill"][[1]])))
    record.dataframe["abs.error.geoMatch"][i,] <- mean(abs((geo.model$coefficients["T"] - it.spdf@data[it.spdf@data$T == 1,]["ie.spill"][[1]])))
    
    
    record.dataframe["est.T.lm"][i,] <- lm.model$coefficients["T"]
    record.dataframe["est.T.matchit"][i,] <- matchit.model$coefficients["T"]
    record.dataframe["est.T.geoMatch"][i,] <- geo.model$coefficients["T"]
    record.dataframe["true.T.ie"][i,] <- mean(it.spdf@data[it.spdf@data$T == 1,]["ie.spill"][[1]])
    record.dataframe["avg.spill"][i,] <- mean(it.spdf@data[it.spdf@data$T == 0,]["t.spill"][[1]])
}

write.csv(record.dataframe, paste('/sciclone/home00/geogdan/geoMatch_testing/Results/out_',date(),runif(1),'.csv',sep=""))
