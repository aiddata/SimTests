source("spatial.simulation.R")
#Example function
#n.points - number of points
#*.Decay - range at which spatial autocorrelation is gone (kilomers, relative to approx 550x550km area)
#sim.T.e - Error multiplier.  Random values between min(x) and max(X) multiplied by this are added to X.
#T.percent - percent of observations which are assigned as treated.
#sim.Y.scale - intercept for outcome calculation
#Theta - treatment effect (main)
#sim.Y.cov.effect - ancillary variable effect on outcome
#sim.Y.het.effect - ancillary variable effect on outcome in treated areas (interaction)
#sim.Y.e - Error multiplier. Random values between min(y) and max(Y) multiplied by this are added to Y in a second stage.
#t.Thresh - number of nearby cases that contribute to spillover averaging.
#spill.mag - Multiplier on spillover effect
test <- spatial.simulation(n.points = 500, 
                           cov.Decay = 750, se.Decay = 500, t.Decay = 500,
                           sim.T.e=0.0, 
                           T.percent = 0.15,
                           sim.Y.scale = 0.0,
                           Theta = 1.0,
                           sim.Y.cov.effect = 0.1,
                           sim.Y.het.effect = 0.0,
                           sim.Y.e = 0.0,
                           t.Thresh = 10,
                           spill.mag = 0.5)



