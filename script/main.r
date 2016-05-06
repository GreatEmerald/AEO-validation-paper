# Main script for AEO: Multivariate Analysis paper
# Copyright (C) 2016 Dainius Masiliunas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

install.packages("inspectr", repos="http://R-Forge.R-project.org")
library(pls)
library(inspectr)

setwd("/home/dainius/Documents/Advanced Earth Observation/AEO-validation-paper/script/")
source("optisim.r")

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
spectra(sDATA) = id ~ pH ~ 420:2500

# Harmonise the names in the australia dataset
AUSTRALIA = australia
names(AUSTRALIA)[names(AUSTRALIA)=="ph"]="pH"

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

# Semi-generic function for running all the different sampling strategies
# Returns a data.frame with relevant statistics (RMSE, R2, etc.) for the amount of samples
# type defines which validation strategy to use: random, bootstrap, stratified, none (if samples are provided)
# strata defines on which values of the response variable to break the samples into strata, only 1 supported right now
# mc.runs defines how many times to run the Monte Carlo analysis (1 means no MC)
# validation.set overrides the validation set used (when unset, non-selected samples are used)
hvalidation = function(type = "none", samples = c(), strata = c(), mc.runs=1, data=DATA, validation.set=NULL, ncomp=compnr, static.seed=TRUE)
{
    if(static.seed)
        set.seed(430001)
    sampleresponse = data.frame()
    for(i in floor((seq(4.5, 16, 0.5))^2))
    {
        rundata = data.frame()
        for(run in 1:mc.runs)
        {
            if(type == "random")
                samp = sample(nrow(data), i)
            if(type == "bootstrap")
                samp = sample(nrow(data), i, replace=TRUE)
            if(type == "stratified" && is.numeric(strata))
            {
                samp1 = sample(as.numeric(rownames(data[data$pH<strata[[1]],])), round(i/2))
                samp2 = sample(as.numeric(rownames(data[data$pH>=strata[[1]],])), round(i/2))
                samp = c(samp1, samp2)
            }
            if(type == "none" && is.numeric(samples))
                samp = samples[1:i]
            tset = data[samp,]
            if(is.data.frame(validation.set))
                vset = validation.set
            else
                vset = data[-samp,]
            vspec = subset(vset, select=c(X420:X1849, X1951:X2500))
            model = fplsr(yvar, data=tset, ncomp=ncomp, validation = "none")
            pred = predict(model, comps = 1:ncomp, newdata = as.matrix(vspec))
            pred.RMSE = sqrt(mean((vset$pH - pred)^2))
            pred.RPD = sd(vset$pH)/pred.RMSE
            rundata = rbind(rundata, data.frame(samples=i, RMSE=pred.RMSE, RPD=pred.RPD))
        }
        sampleresponse = rbind(sampleresponse, colMeans(rundata))
        names(sampleresponse) = names(rundata)
    }
    return(sampleresponse)
}

# Random sampling (without replacement)
hv.rand = hvalidation(type="random")
plot(hv.rand$RMSE~hv.rand$sample, type="l")

# Monte Carlo equivalent
hv.randmc = hvalidation(type="random", mc.runs=100)
plot(hv.randmc$RMSE~hv.randmc$sample, type="l")

# Random sampling (with replacement) AKA bootstrapping:
# For example, if we have 20 samples, we pretend to have more (30) by duplicating some.
# The validation set then uses all of the samples we didn't use for training for best estimates.
# X values are "effective sample size", which is an inflated version of the real one.
# Pretty much requires Monte Carlo treatment!
hv.boot = hvalidation(type="bootstrap")
plot(hv.boot$RMSE~hv.boot$sample, type="l")

# Monte Carlo equivalent
hv.bootmc = hvalidation(type="bootstrap", mc.runs=100)
plot(hv.bootmc$RMSE~hv.bootmc$sample, type="l")

# Stratified sampling: can we tell soil types for non-parametric sampling?
hist(DATA$pH, freq=FALSE, breaks=16)
# Else, parametric optimisation can be applied (like kenstone, but bins instead of samples)
# This can use k-means clustering to determine strata, but how many? 2 (peaks in pH)? 16 (components)? By pH, or other variables?
# Go with an acid/alkaline split for now: median is close to 7 and gives balanced results without bootstrap
hv.strt = hvalidation(type="stratified", strata=c(median(DATA$pH)))

# Repeat the above in Monte Carlo
hv.strtmc = hvalidation(type="stratified", strata=c(median(DATA$pH)), mc.runs=100)

# Kennard-Stone algorithm: select the spectra with the largest differences
# Rather slow! It returns the same things every time, so just running it
# with the highest setting once is enough.
ks.full = kenstone(sDATA, 16^2)
hv.kstn = hvalidation(type="none", samples=ks.full)
plot(hv.kstn$RMSE~hv.kstn$sample, type="l")

# OptiSim: a quicker version of Kennard-Stone, using k=4 as per paper
# This is random subsets of 4 samples, so technically Monte Carlo would apply,
# but the results shouldn't vary enough for it to be worth the time
# (and definitely not if it exceeds Kennard-Stone time of calculation)
os.full = optisim(sDATA, 16^2, 4)
# Factors to numbers; ought to be just as.numeric, but oh well...
os.samples = as.numeric(levels(os.full@id[,1]))[os.full@id[,1]]
hv.osim = hvalidation(type="none", samples=os.samples)
plot(hv.osim$RMSE~hv.osim$sample, type="l")

# One large plot
plot(hv.rand$RMSE~hv.rand$sample, type="l", xlab="Samples in training set", ylab="RMSE", col=1,
    ylim=c(min(min(hv.rand$RMSE), min(hv.boot$RMSE), min(hv.kstn$RMSE), min(hv.osim$RMSE)),
        max(max(hv.rand$RMSE), max(hv.boot$RMSE), max(hv.kstn$RMSE), max(hv.osim$RMSE))))
lines(hv.randmc$RMSE~hv.randmc$sample, col=2)
lines(hv.boot$RMSE~hv.boot$sample, col=7)
lines(hv.bootmc$RMSE~hv.bootmc$sample, col=3)
lines(hv.strt$RMSE~hv.strt$sample, col=8)
lines(hv.strtmc$RMSE~hv.strtmc$sample, col=4)
lines(hv.kstn$RMSE~hv.kstn$sample, col=5)
lines(hv.osim$RMSE~hv.osim$sample, col=6)
legend("topright", c("Random", "Random Monte Carlo", "Bootstrap", "Bootstrap Monte Carlo",
    "Stratified random", "Stratified random Monte Carlo", "Kennard-Stone", "Optisim, k=4"), col=c(1,2,7,3,8,4,5,6), lty=c(1,1,1,1,1,1,1,1))

# RPD
plot(hv.rand$RPD~hv.rand$sample, type="l", xlab="Samples in training set", ylab="RPD", col=1,
    ylim=c(min(min(hv.rand$RPD), min(hv.boot$RPD), min(hv.kstn$RPD), min(hv.osim$RPD)),
        max(max(hv.rand$RPD), max(hv.boot$RPD), max(hv.kstn$RPD), max(hv.osim$RPD))))
lines(hv.randmc$RPD~hv.randmc$sample, col=2)
lines(hv.boot$RPD~hv.boot$sample, col=3)
lines(hv.bootmc$RPD~hv.bootmc$sample, col=4)
lines(hv.kstn$RPD~hv.kstn$sample, col=5)
lines(hv.osim$RPD~hv.osim$sample, col=6)

## Rerun with australia as validation
av.rand = hvalidation(type="random", validation.set=AUSTRALIA)
av.randmc = hvalidation(type="random", mc.runs=100, validation.set=AUSTRALIA)
av.boot = hvalidation(type="bootstrap", validation.set=AUSTRALIA)
av.bootmc = hvalidation(type="bootstrap", mc.runs=100, validation.set=AUSTRALIA)
av.strt = hvalidation(type="stratified", strata=c(median(DATA$pH)), validation.set=AUSTRALIA)
av.strtmc = hvalidation(type="stratified", strata=c(median(DATA$pH)), mc.runs=100, validation.set=AUSTRALIA)
# ks.full has already been calculated, no need to redo that
av.kstn = hvalidation(type="none", samples=ks.full, validation.set=AUSTRALIA)
# No need to recalculate os.samples either
av.osim = hvalidation(type="none", samples=os.samples, validation.set=AUSTRALIA)

plot(av.rand$RMSE~av.rand$sample, type="l", xlab="Samples in training set", ylab="RMSE", col=1,
    ylim=c(min(min(av.rand$RMSE), min(av.boot$RMSE), min(av.kstn$RMSE), min(av.osim$RMSE)),
        max(max(av.rand$RMSE), max(av.boot$RMSE), max(av.kstn$RMSE), max(av.osim$RMSE))))
lines(av.randmc$RMSE~av.randmc$sample, col=2)
lines(av.boot$RMSE~av.boot$sample, col=7)
lines(av.bootmc$RMSE~av.bootmc$sample, col=3)
lines(av.strt$RMSE~av.strt$sample, col=8)
lines(av.strtmc$RMSE~av.strtmc$sample, col=4)
lines(av.kstn$RMSE~av.kstn$sample, col=5)
lines(av.osim$RMSE~av.osim$sample, col=6)
legend("topright", c("Random", "Random Monte Carlo", "Bootstrap", "Bootstrap Monte Carlo",
    "Stratified random", "Stratified random Monte Carlo", "Kennard-Stone", "Optisim, k=4"), col=c(1,2,7,3,8,4,5,6), lty=c(1,1,1,1,1,1,1,1))

plot(av.rand$RPD~hv.rand$sample, type="l", xlab="Samples in training set", ylab="RPD", col=1,
    ylim=c(min(min(av.rand$RPD), min(av.boot$RPD), min(av.kstn$RPD), min(av.osim$RPD)),
        max(max(av.rand$RPD), max(av.boot$RPD), max(av.kstn$RPD), max(av.osim$RPD))))
lines(av.randmc$RPD~av.randmc$sample, col=2)
lines(av.boot$RPD~av.boot$sample, col=7)
lines(av.bootmc$RPD~av.bootmc$sample, col=3)
lines(av.strt$RPD~av.strt$sample, col=8)
lines(av.strtmc$RPD~av.strtmc$sample, col=4)
lines(av.kstn$RPD~av.kstn$sample, col=5)
lines(av.osim$RPD~av.osim$sample, col=6)

## Retry everything with 6 components
hv6.rand = hvalidation(type="random", ncomp=6)
hv6.randmc = hvalidation(type="random", mc.runs=100, ncomp=6)
hv6.boot = hvalidation(type="bootstrap", ncomp=6)
hv6.bootmc = hvalidation(type="bootstrap", mc.runs=100, ncomp=6)
hv6.strt = hvalidation(type="stratified", strata=c(median(DATA$pH)), ncomp=6)
hv6.strtmc = hvalidation(type="stratified", strata=c(median(DATA$pH)), mc.runs=100, ncomp=6)
hv6.kstn = hvalidation(type="none", samples=ks.full, ncomp=6)
hv6.osim = hvalidation(type="none", samples=os.samples, ncomp=6)

plot(hv6.rand$RMSE~hv6.rand$sample, type="l", xlab="Samples in training set", ylab="RMSE", col=1,
    ylim=c(min(min(hv6.rand$RMSE), min(hv6.boot$RMSE), min(hv6.kstn$RMSE), min(hv6.osim$RMSE)),
        max(max(hv6.rand$RMSE), max(hv6.boot$RMSE), max(hv6.kstn$RMSE), max(hv6.osim$RMSE))))
lines(hv6.randmc$RMSE~hv6.randmc$sample, col=2)
lines(hv6.boot$RMSE~hv6.boot$sample, col=7)
lines(hv6.bootmc$RMSE~hv6.bootmc$sample, col=3)
lines(hv6.strt$RMSE~hv6.strt$sample, col=8)
lines(hv6.strtmc$RMSE~hv6.strtmc$sample, col=4)
lines(hv6.kstn$RMSE~hv6.kstn$sample, col=5)
lines(hv6.osim$RMSE~hv6.osim$sample, col=6)
legend("topright", c("Random", "Random Monte Carlo", "Bootstrap", "Bootstrap Monte Carlo",
    "Stratified random", "Stratified random Monte Carlo", "Kennard-Stone", "Optisim, k=4"), col=c(1,2,7,3,8,4,5,6), lty=c(1,1,1,1,1,1,1,1))

# RPD
plot(hv6.rand$RPD~hv6.rand$sample, type="l", xlab="Samples in training set", ylab="RPD", col=1,
    ylim=c(min(min(hv6.rand$RPD), min(hv6.boot$RPD), min(hv6.kstn$RPD), min(hv6.osim$RPD)),
        max(max(hv6.rand$RPD), max(hv6.boot$RPD), max(hv6.kstn$RPD), max(hv6.osim$RPD))))
lines(hv6.randmc$RPD~hv6.randmc$sample, col=2)
lines(hv6.boot$RPD~hv6.boot$sample, col=7)
lines(hv6.bootmc$RPD~hv6.bootmc$sample, col=3)
lines(hv6.strt$RPD~hv6.strt$sample, col=8)
lines(hv6.strtmc$RPD~hv6.strtmc$sample, col=4)
lines(hv6.kstn$RPD~hv6.kstn$sample, col=5)
lines(hv6.osim$RPD~hv6.osim$sample, col=6)

av6.rand = hvalidation(type="random", validation.set=AUSTRALIA, ncomp=6)
av6.randmc = hvalidation(type="random", mc.runs=100, validation.set=AUSTRALIA, ncomp=6)
av6.boot = hvalidation(type="bootstrap", validation.set=AUSTRALIA, ncomp=6)
av6.bootmc = hvalidation(type="bootstrap", mc.runs=100, validation.set=AUSTRALIA, ncomp=6)
av6.strt = hvalidation(type="stratified", strata=c(median(DATA$pH)), validation.set=AUSTRALIA, ncomp=6)
av6.strtmc = hvalidation(type="stratified", strata=c(median(DATA$pH)), mc.runs=100, validation.set=AUSTRALIA, ncomp=6)
# ks.full has already been calculated, no need to redo that
av6.kstn = hvalidation(type="none", samples=ks.full, validation.set=AUSTRALIA, ncomp=6)
# No need to recalculate os.samples either
av6.osim = hvalidation(type="none", samples=os.samples, validation.set=AUSTRALIA, ncomp=6)

plot(av6.rand$RMSE~av6.rand$sample, type="l", xlab="Samples in training set", ylab="RMSE", col=1,
    ylim=c(min(min(av6.rand$RMSE), min(av6.boot$RMSE), min(av6.kstn$RMSE), min(av6.osim$RMSE)),
        max(max(av6.rand$RMSE), max(av6.boot$RMSE), max(av6.kstn$RMSE), max(av6.osim$RMSE))))
lines(av6.randmc$RMSE~av6.randmc$sample, col=2)
lines(av6.boot$RMSE~av6.boot$sample, col=7)
lines(av6.bootmc$RMSE~av6.bootmc$sample, col=3)
lines(av6.strt$RMSE~av6.strt$sample, col=8)
lines(av6.strtmc$RMSE~av6.strtmc$sample, col=4)
lines(av6.kstn$RMSE~av6.kstn$sample, col=5)
lines(av6.osim$RMSE~av6.osim$sample, col=6)
legend("topright", c("Random", "Random Monte Carlo", "Bootstrap", "Bootstrap Monte Carlo",
    "Stratified random", "Stratified random Monte Carlo", "Kennard-Stone", "Optisim, k=4"), col=c(1,2,7,3,8,4,5,6), lty=c(1,1,1,1,1,1,1,1))

plot(av6.rand$RPD~av6.rand$sample, type="l", xlab="Samples in training set", ylab="RPD", col=1,
    ylim=c(min(min(av6.rand$RPD), min(av6.boot$RPD), min(av6.kstn$RPD), min(av.osim$RPD)),
        max(max(av6.rand$RPD), max(av6.boot$RPD), max(av6.kstn$RPD), max(av.osim$RPD))))
lines(av6.randmc$RPD~av6.randmc$sample, col=2)
lines(av6.boot$RPD~av6.boot$sample, col=3)
lines(av6.bootmc$RPD~av6.bootmc$sample, col=4)
lines(av6.strt$RPD~av6.strt$sample, col=5)
lines(av6.strtmc$RPD~av6.strtmc$sample, col=6)
lines(av6.kstn$RPD~av6.kstn$sample, col=7)
lines(av6.osim$RPD~av6.osim$sample, col=8)
