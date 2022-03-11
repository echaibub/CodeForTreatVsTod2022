
require(synapser)
synLogin()

## Read csv file with the subset of the mPower data
## needed to reproduce all the analyses in this paper
## to get access to the data please go over the qualification 
## process described in:
## https://www.synapse.org/#!Synapse:syn4993293/wiki/247860
##
dat <- read.csv(synGet("syn27655357")$path)
dat$createdOn <- as.POSIXct(dat$createdOn)

## Get data from participant id "cdb18d7b-6793-48c0-8085-63c2a71ce5a9"
##
id <- "cdb18d7b-6793-48c0-8085-63c2a71ce5a9"
sdat <- dat[dat$healthCode == id,]

## Order the data according to the dates 
##
sdat <- sdat[order(sdat$createdOn),]

## Get the indexes of activities performed before and after taking medication
##
idx.before <- which(sdat$medTimepoint == "Immediately before Parkinson medication")
idx.after <- which(sdat$medTimepoint == "Just after Parkinson medication (at your best)") 


##############################
## Generate Fig 1
##############################

my.cex <- 1.2
cl <- 1.5
ca <- 1.5
cm <- 1.3
cl2 <- 1.7

mat <- matrix(c(1, 3,
                2, 3), 2, 2, byrow = TRUE)

figpath <- ""

#pdf(paste(figpath, "confounding_by_tod_new.pdf", sep = ""), height = 5, width = 9)
par(mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.75, 0))
layout(mat)
plot(tod ~ createdOn, data = sdat, type = "l", ylim = c(0, 24),
     ylab = "time-of-the-day", xlab = "date", main = id,
     cex.lab = cl, cex.axis = ca, cex.main = cm, col = "lightgray")
points(sdat$createdOn[idx.before], sdat$tod[idx.before], col = "red",
       pch = 20, cex = my.cex)
points(sdat$createdOn[idx.after], sdat$tod[idx.after], col = "blue",
       pch = 20, cex = my.cex)
mtext("(a)", side = 3, at = as.numeric(sdat$createdOn[1])-1500000, line = 0.75, cex = cl)
####
plot(numberTaps ~ createdOn, data = sdat, type = "l", ylim = c(30, 90),
     ylab = "number of taps", xlab = "date", main = id,
     cex.lab = cl, cex.axis = ca, cex.main = cm, col = "lightgray")
points(sdat$createdOn[idx.before], sdat$numberTaps[idx.before], col = "red",
       pch = 20, cex = my.cex)
points(sdat$createdOn[idx.after], sdat$numberTaps[idx.after], col = "blue",
       pch = 20, cex = my.cex)
mtext("(b)", side = 3, at = as.numeric(sdat$createdOn[1])-1500000, line = 0.75, cex = cl)
####
fit <- lm(numberTaps ~ tod, data = sdat)
plot(numberTaps ~ tod, data = sdat, type = "n", ylim = c(30, 90),
     ylab = "number of taps", xlab = "time-of-the-day", main = id,
     cex.lab = cl, cex.axis = ca, cex.main = cm)
abline(a = fit$coefficients[1], b = fit$coefficients[2])
points(sdat$tod[idx.before], sdat$numberTaps[idx.before], col = "red",
       pch = 20, cex = my.cex)
points(sdat$tod[idx.after], sdat$numberTaps[idx.after], col = "blue",
       pch = 20, cex = my.cex)
legend("topright", legend = c("before", "after"), text.col = c("red", "blue"), cex = cl2, bty = "n")
mtext("(c)", side = 3, at = 4.3, line = 0.75, cex = cl)
#dev.off()


