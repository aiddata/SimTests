library(sp)
library(gstat)

#Script to simulate a Spatial Data Frame given a set of inputs.
spatial.simulation <- function(n.points,
                               cov.Decay, se.Decay, t.Decay,
                               sim.T.e, 
                               T.percent, 
                               sim.Y.scale,
                               Theta,
                               sim.Y.cov.effect,
                               sim.Y.het.effect,
                               sim.Y.e,
                               t.Thresh,
                               spill.mag)
{
 
  #generate random coordinates
  random.longitude <- runif(n.points, -5.0, 5.0)
  random.latitude <- runif(n.points, -5.0, 5.0)
  
  # create dataframe
  spdf <- data.frame(id = 1:n.points,
                     longitude = random.longitude,
                     latitude = random.latitude)
  
  # convert to spatial points dataframe
  coordinates(spdf) <- c("longitude", "latitude")
  proj4string(spdf) <- CRS("+proj=longlat +datum=WGS84")

  #Iteratively simulates values into this spatial field following spatial
  #distance-decay models (spherical assumed).
  #(A) A covariate (X)
  #(B) A spatial error field (e)
  
  var.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                       model=vgm(psill=1, model="Sph", range=cov.Decay),
                       nmax=10)
  var.sim <- predict(var.g.dummy, newdata=spdf, nsim=1)
  var.sim$sim1 <- var.sim$sim1 / max(var.sim$sim1)
  spdf@data$X <- var.sim$sim1

  
  var.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                       model=vgm(psill=1, model="Sph", range=se.Decay),
                       nmax=10)
  var.sim <- predict(var.g.dummy, newdata=spdf, nsim=1)
  var.sim$sim1 <- var.sim$sim1 / max(var.sim$sim1)
  spdf@data$e <- var.sim$sim1

  
  #Each unit of observation (i) is assigned a sim probability
  #to receive treatment, calculated by:
  #sim.T.prob = Xi + sim.T.e
  #where sim.T.e is a parameterized error function representing random error
  #introduced into the Xi / T relationship.
  e.vec <- runif(n.points, 
                 min(spdf@data$X),
                 max(spdf@data$X)) * sim.T.e
  
  spdf@data$sim.T.prob <- spdf@data$X + e.vec 
  min.vec <- (spdf@data$sim.T.prob + abs(min(spdf@data$sim.T.prob)))
  max.vec <- (max(spdf@data$sim.T.prob) + abs(min(spdf@data$sim.T.prob)))
  spdf@data$sim.T.prob <- min.vec / max.vec 
    
  #Assign T based on the probabilities generated
  spdf@data$T <- 0
  spdf@data$T[sample.int(n=n.points, size=(T.percent * n.points), prob=spdf@data$sim.T.prob)] <- 1
    
  
  #Once treatment is defined, an outcome measure is generated
  #for estimation according to a linear model:
  #Yi = sim.Y.scale + (Theta * Ti) + (sim.Y.cov.effect * Xi) + (sim.Y.het.effect * (Ti * Xi)) + sim.Y.e
  #where sim.Y.scale is a universal scalar, theta the main treatment effect,
  #sim.Y.cov.effect the effect of covariate X, and sim.Y.het.effect the effect of heterogeneity
  #in impacts across covariate X.  sim.Y.e represents simulation noise.
  
  spdf@data$Y <- sim.Y.scale + 
    (Theta * spdf@data$T) + 
    (sim.Y.cov.effect * spdf@data$X) +
    (sim.Y.het.effect * (spdf@data$T * spdf@data$X))
  

  spdf@data$Y <- spdf@data$Y + runif(n.points, 
                                min(spdf@data$Y),
                                max(spdf@data$Y)) * sim.Y.e
  
  #For each unit of observation Ti = 1, the total effect (Theta + sim.Y.het.effect)
  #is calculated.  
  
  spdf@data$ie.nospill <- (Theta * spdf@data$T) + 
    (sim.Y.het.effect * (spdf@data$T * spdf@data$X))
  
  
  #For each C, spillover is calculated from all Ti values based on spatial distances.
  #The average of Ti values is used to estimate spillover, based on the closest t.Thresh percent
  #of treated cases (according to a spherical distance decay function).
  spdf@data$t.spill <- 0
  t.vals <- spdf@data[spdf@data$T == 1,]["ie.nospill"]
  ct.dists <- spDists(x=spdf[spdf@data$T == 0,], y=spdf[spdf@data$T == 1,])
  for (i in 1:length(spdf[spdf@data$T == 0,]))
  {

     t.weights <- 1-(((3/2) * (ct.dists[i,] / t.Decay) - (1/2) * (ct.dists[i,] / t.Decay)^3))
     
     #select t.Thresh closest observations; ties are included (so total can exceed this value).
     max.t.dist <- max(head(sort(ct.dists[i,]), t.Thresh))
     closest.vector <- (ct.dists[i,] <= max.t.dist)
     
     spdf@data[spdf@data$T == 0,][i,]["t.spill"] <- mean(t.weights[closest.vector] * t.vals[closest.vector,]) * spill.mag
    
  }
  

  #Final impact estimate is generated, including spillover effects:
  spdf@data$ie.spill <- spdf@data$t.spill + spdf@data$ie.nospill
  
  
  
  return(spdf)
}





