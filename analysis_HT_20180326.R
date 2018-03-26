#################
#               #
#  HT INSIGHT   #
#  26 MAR 2018  # 
#               #
#################



library(MASS)
library(lme4)

rm(list = ls())

setwd("/Users/nsircar/Dropbox/Papers/HTAnalysis")

dat <- read.csv("/Users/nsircar/Dropbox/Papers/HTAnalysis/Renomination20092014.csv", stringsAsFactors = F)

logmoveable <- log10(dat$Moveable)
elected <- as.numeric(dat$Position2009 == 1)
crim <- as.numeric(dat$CriminalCase)
pg <- as.numeric(dat$PostGraduate)
party <- as.numeric(as.factor(paste(dat$Party2009, dat$State)))
elected <- as.numeric(dat$Position2009 == 1)

keepval <- which(logmoveable > 0 & !is.na(logmoveable)) ## Subset of data to analyze - 1020 candidates


## Renomination Rates

table(dat$ReNominated[keepval]) ## Table of Renomination
mean(dat$ReNominated[keepval]) ## Rate of Renomination

mean(dat$ReNominated[keepval][elected[keepval] == 1]) ## Elected Renomination Rate
mean(dat$ReNominated[keepval][elected[keepval] == 0]) ## Unelected Renomination Rate

## Measuring Renomination in the Same Party

sameparty <- as.numeric(dat$Party2009[keepval] == dat$Party2014[keepval])[dat$ReNominated[keepval] == 1]
sameparty[dat$Party2009[keepval][dat$ReNominated[keepval] == 1] == "AUDF" & dat$Party2014[keepval][dat$ReNominated[keepval] == 1] == "AIUDF"] <- 1

1 - mean(sameparty, na.rm = T) ## Rate of Party Change

1 - mean(sameparty[elected[keepval][dat$ReNominated[keepval] == 1] == 1], na.rm = T) ## For Elected
1 - mean(sameparty[elected[keepval][dat$ReNominated[keepval] == 1] == 0], na.rm = T) ## For Unelected


## Run model and simulate coefficients

fit <- glmer(dat$ReNominated[keepval]~logmoveable[keepval]+crim[keepval]+pg[keepval] + (1|party[keepval]), family = binomial(link = "logit")) ## Core Model for Piece

summary(fit)


fit2 <- glmer(dat$ReNominated[keepval]~logmoveable[keepval]+crim[keepval]+pg[keepval] + elected[keepval] + (1|party[keepval]), family = binomial(link = "logit")) ## Control for whether elected in 2009
summary(fit2) ## Wealth and Criminality Effects seem to be driven by electability

coef.sims <- mvrnorm(1000, summary(fit)$coef[,1], summary(fit)$vcov)  ## Simulation from Original Model

## Holding at the Mean for other predictors for predicted values

crimval <- cbind(rep(1,2), rep(mean(logmoveable[keepval]), 2), c(0,1), rep(mean(pg[keepval]), 2)) ## X values for Crim Case Simulation -- Mean moveable assets, mean PG completion
pgval <- cbind(rep(1,2), rep(mean(logmoveable[keepval]), 2), rep(mean(crim[keepval]), 2),  c(0,1)) ## X values for PG degree Simulation -- Mean moveable assets, mean criminality



invlogit <- function(x) 1/(1 + exp(-x))  ## Inverse Logit function


## Making Figure 1

vals <- c(apply(invlogit((pgval %*% t(coef.sims))), 1, mean), apply(invlogit((crimval %*% t(coef.sims))), 1, mean)) ## Predicted Values from Simulations
vals.ci <- cbind(apply(invlogit((pgval %*% t(coef.sims))), 1, quantile, c(.05, .95)), apply(invlogit((crimval %*% t(coef.sims))), 1, quantile, c(.05, .95))) ## 90% intervals from simulations

pdf("/Users/nsircar/Dropbox/Papers/HTanalysis/renom_fig1.pdf", width=8, height=6)

par(mar = rep(0,4), mai = rep(0,4))
plot(0,0, xlim = c(-.85,7), ylim = c(-.08, .6), ann = F, axes = F, type = "n")
x.coords <- barplot(vals, col = c(rgb(1,125/255,125/255), rgb(1,125/255,125/255), rgb(102/255,178/255,1), rgb(102/255,178/255,1)), space = c(.5,.5,1,.5), axes = F, add = T )
axis(2, at=seq(0,.6,.1), labels = seq(0,.6,.1), pos = 0)

cats <- c("Graduate\nor Less", "Post\nGraduate", "Not Major\nCriminal", "Major\nCriminal")

cols <- c(rgb(204/255,0,0), rgb(204/255,0,0), rgb(0,0,204/255), rgb(0,0,204/255))
for (i in 1:4){
  
  text(x.coords[i],-.05, cats[i], cex = 1.2)
  arrows(x.coords[i], vals.ci[1,i], x.coords[i], vals.ci[2,i], col = cols[i], code = 3, angle = 90, length = .07, lwd = 2)
  text(x.coords[i]+.3, vals[i], round(vals[i], 2), pos = 3, cex = 1.5)
  
  
}

text(-.83, .3, "Estimated Probability of Renomination", srt=90, cex = 1.2)

dev.off()




xval <- logmoveable[keepval][log10(1000000) <= logmoveable[keepval] & log10(5000000) >= logmoveable[keepval]]
xpred <- cbind(rep(1, length(xval)), xval, rep(mean(crim[keepval]), length(xval)), rep(mean(pg[keepval]), length(xval)))

mean(colMeans(invlogit((xpred %*% t(coef.sims)))))  ## Mean Renomination btw 10 and 50 Lakh

xval <- logmoveable[keepval][log10(10000000) <= logmoveable[keepval] ]
xpred <- cbind(rep(1, length(xval)), xval, rep(mean(crim[keepval]), length(xval)), rep(mean(pg[keepval]), length(xval)))

mean(colMeans(invlogit((xpred %*% t(coef.sims))))) ## Mean Renomination Crorepati


## Making Figure 2

cutpoints <- quantile(logmoveable[keepval], seq(0,1,.05)) ## Generate Cutpoints for "Bins" every 5 Percentile Points


vec <- seq(5,8, .01) ## Values over which Simulated (Log) Moveable Assets Displayed

assetval <- cbind(rep(1, length(vec)), vec, rep(mean(crim[keepval]), length(vec)), rep(mean(pg[keepval]), length(vec))) ## Mean levels of PG degree completion and criminality

X <- cbind(rep(1, length(dat$ReNominated)), logmoveable, crim, pg)[keepval, ] ## Predictor Matrix

renom <- dat$ReNominated

xbin5 <- ybin5 <- error <- rep(NA, length(cutpoints) - 1)  ## xbin5 = Binned X value; ybin5 = Binned Y value
for (i in 2:length(cutpoints)){
  xbin5[i-1] <- mean(cutpoints[(i-1):i])  
  
  
  xval <- logmoveable[keepval][(logmoveable[keepval] > cutpoints[i-1]) & (logmoveable[keepval] <= cutpoints[i]) ]
  yval <- renom[keepval][(logmoveable[keepval] > cutpoints[i-1]) & (logmoveable[keepval] <= cutpoints[i]) ]
  estval <- invlogit(X[(logmoveable[keepval] > cutpoints[i-1]) & (logmoveable[keepval] <= cutpoints[i]), ] %*% fixef(fit))
  
  error[i-1] <- mean(yval - estval)
  
  fitval <- mean(invlogit(assetval[(assetval[,2] > cutpoints[i-1]) & (assetval[,2] <= cutpoints[i]), ] %*% fixef(fit)))
  
  ybin5[i-1] <- fitval + error[i-1]  
}  


vals <- apply(invlogit((assetval %*% t(coef.sims))), 1, mean) ## Predicted Value from Simulations
vals.ci <- apply(invlogit((assetval %*% t(coef.sims))), 1, quantile, c(.05, .95)) ## 90% Intervals from Simulations

colsline <- rgb(153/255,255/255,153/255)
colsfill <- rgb(0,102/255,0)

pdf("/Users/nsircar/Dropbox/Papers/HTanalysis/renom_fig2.pdf", width=8, height=6)

par(mar = rep(0,4), mai = rep(0,4))
plot(0,0, xlim = c(min(vec) - .4, max(vec)+.1), ylim = c(-.1, 0.5), ann =F, axes = F, type = "n")


polygon(c(vec, rev(vec)), c(vals.ci[1,], rev(vals.ci[2,])), col = colsline,border = NA   )  ## Colored Interval Band
points(vec, vals, lwd = 2.5, type = "l", col = colsfill)  ## Predicted Value Curve
points(xbin5[2:19], ybin5[2:19], pch=19, col = colsfill, cex = 1.5) ## Plotting Points for Binned X and Y values

axis(1, at = seq(5, 8, 1), labels = c("1 Lakh", "10 Lakh", "1 Crore", "10 Crore"),  pos=0)  ## Axis on a Log (base 10) scale
axis(2, at = seq(0,.5,.1), pos = 5)

text(6.5, -.09, "Reported Candidate Asset Wealth (Moveable)", cex = 1.2)
text(4.65, .25, "Estimated Probability of Renomination", cex = 1.2, srt = 90)

dev.off()
