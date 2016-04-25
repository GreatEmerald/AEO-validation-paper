# Main script for AOE: Multivariate Analysis paper
install.packages("inspectr", repos="http://R-Forge.R-project.org")
library(pls)
library(inspectr)
source("optisim.r")

setwd("/home/dainius/Documents/Advanced Earth Observation/AEO-validation-paper/script/")

# Load data (assemble from 2 CSVs)
DATA1 = read.csv("soilNL_test.csv", header=TRUE);
DATA2 = read.csv("soilNL_training.csv", header=TRUE);
names(DATA2)[4] = "sqrt_OM"
DATA = rbind(DATA1, DATA2)
rm(DATA1)
rm(DATA2)

# Reorder by sample ID. This is because the samples were ordered by SOM and
# files were split using interleaved sampling, whereas a usual database would
# have them ordered by whichever sample was taken first.
DATA = DATA[order(DATA$Sample, DATA$RowID),]

# Make into an inspectr SpectraDataFrame
sDATA = DATA
spectra(sDATA) = id ~ pH + Moisture + OM ~ 420:2500

# Plot validation plots
# Wrapper function for simpler plsr
fplsr = function(yvar, xstart=7, data, ...)
{
    xmat = as.matrix(data[xstart:ncol(data)])
    ymat = as.numeric(unlist(data[yvar]))
    return(plsr(ymat~xmat, ...))
}
yvar = "pH"

# Leave-one-out (equivalent to validation="CV", segments=nrow(DATA))
loo = fplsr(yvar, data=DATA, ncomp=20, validation = "LOO")
validationplot(loo, legendpos="topright") 

# Random: 2 groups, one compared with the other
cv.rand = fplsr(yvar, data=DATA, ncomp=20, validation = "CV", segments=2, segment.type="random")
validationplot(cv.rand, legendpos="topright", main="Random sampling")

# Consecutive: 2 groups, one compared with the other
cv.cons = fplsr(yvar, data=DATA, ncomp=20, validation = "CV", segments=2, segment.type="consecutive")
validationplot(cv.cons, legendpos="topright", main="Consecutive sampling")

# Interleaved: 2 groups, one compared with the other
cv.intl = fplsr(yvar, data=DATA, ncomp=20, validation = "CV", segments=2, segment.type="interleaved")
validationplot(cv.intl, legendpos="topright", main="Interleaved sampling")

# Random: 3 groups, one compared with the other
cv.rand3 = fplsr(yvar, data=DATA, ncomp=20, validation = "CV", segments=3, segment.type="random")
validationplot(cv.rand3, legendpos="topright", main="Random sampling")

# Consecutive: 3 groups, one compared with the other
cv.cons3 = fplsr(yvar, data=DATA, ncomp=20, validation = "CV", segments=3, segment.type="consecutive")
validationplot(cv.cons3, legendpos="topright", main="Consecutive sampling")

# Interleaved: 3 groups, one compared with the other
cv.intl3 = fplsr(yvar, data=DATA, ncomp=20, validation = "CV", segments=3, segment.type="interleaved")
validationplot(cv.intl3, legendpos="topright", main="Interleaved sampling")

# Monte Carlo: 100 random 2-group sampling
mc.rmsep = data.frame()
for(i in 1:100)
{
    cv.mc = fplsr(yvar, data=DATA, ncomp=20, validation = "CV", segments=2, segment.type="random")
    mc.rmsep = rbind(mc.rmsep, data.frame(t(RMSEP(cv.mc)$val[2,1,])))
}
mc.cvadj = colMeans(mc.rmsep)
plot(mc.cvadj ~ as.numeric(0:(length(mc.cvadj)-1)), type="l")

# Monte Carlo: 100 random 3-group sampling
mc.rmsep = data.frame()
for(i in 1:100)
{
    cv.mc = fplsr(yvar, data=DATA, ncomp=20, validation = "CV", segments=3, segment.type="random")
    mc.rmsep = rbind(mc.rmsep, data.frame(t(RMSEP(cv.mc)$val[2,1,])))
}
mc.cvadj3 = colMeans(mc.rmsep)
plot(mc.cvadj3 ~ as.numeric(0:(length(mc.cvadj3)-1)), col=8)

# Showcase of the effect of increasing group number
validationplot(loo, main="")
rgramp = c("black", "magenta", "orange", "cyan", "green", "yellow", "grey", "blue")
for(i in 2^(1:8))
{
    iter.intl = fplsr(yvar, data=DATA, ncomp=20, validation = "CV", segments=i, segment.type="interleaved")
    lines(RMSEP(iter.intl)$val[2,1,] ~ RMSEP(iter.intl)$comp, col=rgramp[log(i, 2)])
}
lines(RMSEP(loo)$val[2,1,] ~ RMSEP(loo)$comp, col="red", lwd=3, lty=3)
legend("topright", c("k=2", "k=4", "k=8", "k=16", "k=32", "k=64", "k=128", "k=256", "LOO"), col=c(rgramp,"red"), lty=c(1,1,1,1,1,1,1,1,3), lwd=c(1,1,1,1,1,1,1,1,3))

# Plot all 2-group adjCV results in one large graph
validationplot(loo, estimate="adjCV", main="") 
# Random
lines(RMSEP(cv.rand)$val[2,1,] ~ RMSEP(cv.rand)$comp, col=2)
# Consecutive
lines(RMSEP(cv.cons)$val[2,1,] ~ RMSEP(cv.cons)$comp, col=3)
# Interleaved
lines(RMSEP(cv.intl)$val[2,1,] ~ RMSEP(cv.intl)$comp, col=4)
# Monte Carlo
lines(mc.cvadj ~ as.numeric(0:(length(mc.cvadj)-1)), col=5)
legend("topright", c("Leave-one-out", "Random", "Consecutive", "Interleaved", "Monte Carlo, 100 runs"), col=1:5, lty=c(1,1,1,1,1))

# For reference:
# CV
#lines(RMSEP(cv.intl3)$val[1,1,] ~ RMSEP(cv.intl3)$comp)
# AdjCV
#lines(RMSEP(cv.intl3)$val[2,1,] ~ RMSEP(cv.intl3)$comp)

## Holdout validation
# Number of components will be 16 (alternatively, use 6 or 12)
compnr = 16
# Graph of training set size vs RMSE 

# Random sampling (without replacement)
sampleresponse = data.frame()
for(i in floor((seq(4.5, 16, 0.5))^2))
{
    samp = sample(nrow(DATA), i)
    tset = DATA[samp,]
    vset = DATA[-samp,]
    vspec = vset[7:ncol(vset)]
    rsm = fplsr(yvar, data=tset, ncomp=compnr, validation = "none")
    pred.tst <- predict(rsm, comps = 1:compnr, newdata = as.matrix(vspec))
    RMSE.xval <- sqrt(mean((vset$pH - pred.tst)^2))
    sampleresponse = rbind(sampleresponse, data.frame(i=i, RMSE=RMSE.xval))
    #plot (pred.tst ~ vset$pH, xlab="Observed (validation)", ylab="Predicted (training)", main="pH observed vs predicted")
    #abline(0, 1)
}
plot(sampleresponse$RMSE~sampleresponse$i, type="l")

# Random sampling (with replacement) AKA bootstrapping:
# For example, if we have 20 samples, we pretend to have more (30) by duplicating some.
# The validation set then uses all of the samples we didn't use for training for best estimates.
# X values are "effective sample size", which is an inflated version of the real one.
# Pretty much requires Monte Carlo treatment!
sampleresponse = data.frame()
for(i in floor((seq(4.5, 16, 0.5))^2))
{
    samp = sample(nrow(DATA), i, replace=TRUE)
    tset = DATA[samp,]
    vset = DATA[-samp,]
    vspec = vset[7:ncol(vset)]
    rsm = fplsr(yvar, data=tset, ncomp=compnr, validation = "none")
    pred.tst <- predict(rsm, comps = 1:compnr, newdata = as.matrix(vspec))
    RMSE.xval <- sqrt(mean((vset$pH - pred.tst)^2))
    sampleresponse = rbind(sampleresponse, data.frame(i=i, RMSE=RMSE.xval))
    #plot (pred.tst ~ vset$pH, xlab="Observed (validation)", ylab="Predicted (training)", main="pH observed vs predicted")
    #abline(0, 1)
}
plot(sampleresponse$RMSE~sampleresponse$i, type="l")

# Stratified sampling: can we tell soil types for non-parametric sampling?

hist(DATA$pH, freq=FALSE, breaks=16)
# Else, parametric optimisation can be applied (like kenstone, but bins instead of samples)
# This can use k-means clustering to determine strata, but how many? 2 (peaks in pH)? 16 (components)? By pH, or other variables?
# Leave for later

# Repeat the above in Monte Carlo

# Kennard-Stone algorithm: select the spectra with the largest differences
# Rather slow! It returns the same things every time, so just running it
# with the highest setting once is enough.
sampleresponse = data.frame()
sampfull = kenstone(sDATA, floor(16^2))
for(i in floor((seq(4.5, 16, 0.5))^2))
{
    samp = sampfull[1:i]
    tset = DATA[samp,]
    vset = DATA[-samp,]
    vspec = vset[7:ncol(vset)]
    rsm = fplsr(yvar, data=tset, ncomp=compnr, validation = "none")
    pred.tst <- predict(rsm, comps = 1:compnr, newdata = as.matrix(vspec))
    RMSE.xval <- sqrt(mean((vset$pH - pred.tst)^2))
    sampleresponse = rbind(sampleresponse, data.frame(i=i, RMSE=RMSE.xval))
    #plot (pred.tst ~ vset$pH, xlab="Observed (validation)", ylab="Predicted (training)", main="pH observed vs predicted")
    #abline(0, 1)
}
plot(sampleresponse$RMSE~sampleresponse$i, type="l", xlab="Samples in training set", ylab="RMSE")

# OptiSim: a quicker version of Kennard-Stone, using k=4 as per paper
# This is random subsets of 4 samples, so technically Monte Carlo would apply,
# but the results shouldn't vary enough for it to be worth the time
# (and definitely not if it exceeds Kennard-Stone time of calculation)
sampleresponse = data.frame()
sampfull = optisim(sDATA, 16^2, 4)
for(i in floor((seq(4.5, 16, 0.5))^2))
{
    # Rather odd way of converting factors to IDs, since optisim() returns selected SpectraDataFrame,
    # whereas we need the non-selected ones too
    samp = as.numeric(levels(sampfull@id[,1]))[sampfull@id[,1]][1:i]
    tset = DATA[samp,]
    vset = DATA[-samp,]
    vspec = vset[7:ncol(vset)]
    rsm = fplsr(yvar, data=tset, ncomp=compnr, validation = "none")
    pred.tst <- predict(rsm, comps = 1:compnr, newdata = as.matrix(vspec))
    RMSE.xval <- sqrt(mean((vset$pH - pred.tst)^2))
    sampleresponse = rbind(sampleresponse, data.frame(i=i, RMSE=RMSE.xval))
    #plot (pred.tst ~ vset$pH, xlab="Observed (validation)", ylab="Predicted (training)", main="pH observed vs predicted")
    #abline(0, 1)
}
plot(sampleresponse$RMSE~sampleresponse$i, type="l", xlab="Samples in training set", ylab="RMSE")
