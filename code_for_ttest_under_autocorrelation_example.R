
####################################
## code for generating Fig 5
####################################

RunSimulationsUnderNull <- function(nSim, n, ar) {  
  ttestPvals <- matrix(NA, nSim, 2)
  colnames(ttestPvals) <- c("randomizedDesign", "pairedDesign") 
  acEstimates <- rep(NA, nSim) 
  
  for (i in seq(nSim)) {
    cat(i, "\n")
    
    ## simulate AR(1) time series data under the null
    y <- as.vector(arima.sim(n = n, model = list(ar = ar)))
    
    ## estimate the autocorrelation at lag 1
    acEstimates[i] <- acf(y, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
    
    ## sample binary variables representing the before/after labels
    g1 <- sample(c(0, 1), n, replace = TRUE) ## random design
    g2 <- rep(c(0, 1), n/2) ## paired design
    
    ## convert  the binary labels to character labels
    gg1 <- c("before", "after")[g1 + 1]
    gg2 <- c("before", "after")[g2 + 1]
    
    ## get t-test p-values
    ttestPvals[i, 1] <- t.test(y ~ gg1, var.equal = TRUE)$p.value
    ttestPvals[i, 2] <- t.test(y ~ gg2, var.equal = TRUE)$p.value
  } 
  
  list(ttestPvals = ttestPvals,
       acEstimates = acEstimates)
}


######################################################
######################################################
######################################################


## simulations that generate the outputs in panels j to r
nSim <- 10000
n <- 60
set.seed(123456789)
sim0 <- RunSimulationsUnderNull(nSim, n, ar = 0)
simP <- RunSimulationsUnderNull(nSim, n, ar = 0.95)
simN <- RunSimulationsUnderNull(nSim, n, ar = -0.95)
#save(sim0, simP, simN, file = "randomized_vs_paired_designs.RData", compress = TRUE)
#load("randomized_vs_paired_designs.RData")


## synthetic data examples used for panels a to i

n <- 60

gp <- rep(c("b", "a"), n/2)
set.seed(1234567)
gr <- sample(c("b", "a"), n, replace = TRUE)


set.seed(12345)
y1 <- 10 + as.vector(arima.sim(n = n, model = list(ar = 0)))

set.seed(123456)
y2 <- 10 + as.vector(arima.sim(n = n, model = list(ar = 0.95)))

set.seed(1234567)
y3 <- 10 + as.vector(arima.sim(n = n, model = list(ar = -0.95)))


## negative autocorrelation
t.test(y3 ~ gp, paired = FALSE, var.equal = FALSE) ## paired case
t.test(y3 ~ gr, paired = FALSE, var.equal = FALSE) ## random case


## positive autocorrelation
t.test(y2 ~ gp, paired = FALSE, var.equal = FALSE) ## paired case
t.test(y2 ~ gr, paired = FALSE, var.equal = FALSE) ## random case



mat <- matrix(c(1, 1, 2, 2, 2, 3, 3, 3, 
                4, 4, 5, 5, 5, 6, 6, 6,
                7, 7, 8, 8, 8, 9, 9, 9,
                10, 10, 11, 11, 11, 12, 12, 12,
                13, 13, 14, 14, 14, 15, 15, 15,
                16, 16, 17, 17, 17, 18, 18, 18), 6, 8, byrow = TRUE)


ca <- 1.2
myline <- 0.5
cl <- 1.4
ca2 <- 1.4
cm <- 1.6
mydensi <- 100
cd <- 1.75
nc <- 30


figpath <- ""

## generates Fig 5
#pdf(file = paste0(figpath, "ttest_examples_new.pdf"), width = 10, height = 10)
par(mar = c(3, 3, 2, 0.5), mgp = c(1.5, 0.5, 0))
layout(mat)
#####
acf(y3, main = "")
mtext("(a)", side = 3, at = 0, cex = ca, line = myline)
##
plot(y3, type = "n", xlab = "day", xaxt = "n", main = "negative autocorrelation, paired", ylab = "outcome",
     cex.lab = cl, cex.axis = ca2, cex.main = cm)
axis(side = 1, at = seq(1.5, 60.5, by = 2), labels = seq(30))
rect(xleft = seq(0.5, 59.5, by = 4), xright = seq(2.5, 60.5, by = 4), ybottom = 0, ytop = 20, 
     density = mydensi, col = "lightgrey")
lines(y3, type = "l")
points(which(gp == "b"), y3[which(gp == "b")], col = "red", pch = 20, cex = cd)
points(which(gp == "a"), y3[which(gp == "a")], col = "blue", pch = 20, cex = cd)
mtext("(b)", side = 3, at = 0, cex = ca, line = myline)
##
plot(y3, type = "n", xlab = "day", xaxt = "n", main = "negative autocorrelation, random", ylab = "outcome",
     cex.lab = cl, cex.axis = ca2, cex.main = cm)
axis(side = 1, at = seq(1, 60, by = 1), labels = seq(60))
rect(xleft = seq(0.5, 59.5, by = 2), xright = seq(1.5, 60.5, by = 2), ybottom = 0, ytop = 20, 
     density = mydensi, col = "lightgrey")
lines(y3, type = "l")
points(which(gr == "b"), y3[which(gr == "b")], col = "red", pch = 20, cex = cd)
points(which(gr == "a"), y3[which(gr == "a")], col = "blue", pch = 20, cex = cd)
mtext("(c)", side = 3, at = 0, cex = ca, line = myline)
#####
acf(y2, main = "")
mtext("(d)", side = 3, at = 0, cex = ca, line = myline)
##
plot(y2, type = "n", xlab = "day", xaxt = "n", main = "positive autocorrelation, paired", ylab = "outcome",
     cex.lab = cl, cex.axis = ca2, cex.main = cm)
axis(side = 1, at = seq(1.5, 60.5, by = 2), labels = seq(30))
rect(xleft = seq(0.5, 59.5, by = 4), xright = seq(2.5, 60.5, by = 4), ybottom = 0, ytop = 20, 
     density = mydensi, col = "lightgrey")
lines(y2, type = "l")
points(which(gp == "b"), y2[which(gp == "b")], col = "red", pch = 20, cex = cd)
points(which(gp == "a"), y2[which(gp == "a")], col = "blue", pch = 20, cex = cd)
mtext("(e)", side = 3, at = 0, cex = ca, line = myline)
##
plot(y2, type = "n", xlab = "day", xaxt = "n", main = "positive autocorrelation, random", ylab = "outcome",
     cex.lab = cl, cex.axis = ca2, cex.main = cm)
axis(side = 1, at = seq(1, 60, by = 1), labels = seq(60))
rect(xleft = seq(0.5, 59.5, by = 2), xright = seq(1.5, 60.5, by = 2), ybottom = 0, ytop = 20, 
     density = mydensi, col = "lightgrey")
lines(y2, type = "l")
points(which(gr == "b"), y2[which(gr == "b")], col = "red", pch = 20, cex = cd)
points(which(gr == "a"), y2[which(gr == "a")], col = "blue", pch = 20, cex = cd)
mtext("(f)", side = 3, at = 0, cex = ca, line = myline)
#####
acf(y1, main = "")
mtext("(g)", side = 3, at = 0, cex = ca, line = myline)
##
plot(y1, type = "n", xlab = "day", xaxt = "n", main = "no autocorrelation, paired", ylab = "outcome",
     cex.lab = cl, cex.axis = ca2, cex.main = cm)
axis(side = 1, at = seq(1.5, 60.5, by = 2), labels = seq(30))
rect(xleft = seq(0.5, 59.5, by = 4), xright = seq(2.5, 60.5, by = 4), ybottom = 0, ytop = 20, 
     density = mydensi, col = "lightgrey")
lines(y1, type = "l")
points(which(gp == "b"), y1[which(gp == "b")], col = "red", pch = 20, cex = cd)
points(which(gp == "a"), y1[which(gp == "a")], col = "blue", pch = 20, cex = cd)
mtext("(h)", side = 3, at = 0, cex = ca, line = myline)
##
plot(y1, type = "n", xlab = "day", xaxt = "n", main = "no autocorrelation, random", ylab = "outcome",
     cex.lab = cl, cex.axis = ca2, cex.main = cm)
axis(side = 1, at = seq(1, 60, by = 1), labels = seq(60))
rect(xleft = seq(0.5, 59.5, by = 2), xright = seq(1.5, 60.5, by = 2), ybottom = 0, ytop = 20, 
     density = mydensi, col = "lightgrey")
lines(y1, type = "l")
points(which(gr == "b"), y1[which(gr == "b")], col = "red", pch = 20, cex = cd)
points(which(gr == "a"), y1[which(gr == "a")], col = "blue", pch = 20, cex = cd)
mtext("(i)", side = 3, at = 0, cex = ca, line = myline)
####
####
####
par(mar = c(3, 3, 3, 0.5), mgp = c(1.5, 0.5, 0))
hist(simN$acEstimates, probability = TRUE, main = "negative autocorr.", 
     xlab = "autocorrelation (lag = 1)", cex.axis = ca2, cex.lab = cl, cex.main = cm)
mtext("(j)", side = 3, at = -1, cex = ca, line = myline)
##
hist(simN$ttestPvals[, "pairedDesign"], probability = TRUE, xlim = c(0, 1),
     main = "negative autocorr., paired", xlab = "p-values", nclass = nc, 
     cex.axis = ca2, cex.lab = cl, cex.main = cm)
mtext("(k)", side = 3, at = 0, cex = ca, line = myline)
##
hist(simN$ttestPvals[, "randomizedDesign"], probability = TRUE, xlim = c(0, 1), 
     main = "negative autocorr., random", xlab = "p-values", nclass = nc, 
     cex.axis = ca2, cex.lab = cl, cex.main = cm)
mtext("(l)", side = 3, at = 0, cex = ca, line = myline)
####
hist(simP$acEstimates, probability = TRUE, main = "positive autocorr.", 
     xlab = "autocorrelation (lag = 1)", cex.axis = ca2, cex.lab = cl, cex.main = cm)
mtext("(m)", side = 3, at = 0.2, cex = ca, line = myline)
##
hist(simP$ttestPvals[, "pairedDesign"], probability = TRUE, xlim = c(0, 1),
     main = "positive autocorr., paired", xlab = "p-values", nclass = 15, 
     cex.axis = ca2, cex.lab = cl, cex.main = cm)
mtext("(n)", side = 3, at = 0, cex = ca, line = myline)
##
hist(simP$ttestPvals[, "randomizedDesign"], probability = TRUE, xlim = c(0, 1), 
     main = "positive autocorr., random", xlab = "p-values", nclass = nc, 
     cex.axis = ca2, cex.lab = cl, cex.main = cm)
mtext("(o)", side = 3, at = 0, cex = ca, line = myline)
####
hist(sim0$acEstimates, probability = TRUE, main = "no autocorrelation", 
     xlab = "autocorrelation (lag = 1)", cex.axis = ca2, cex.lab = cl, cex.main = cm)
mtext("(p)", side = 3, at = -0.5, cex = ca, line = myline)
##
hist(sim0$ttestPvals[, "pairedDesign"], probability = TRUE, xlim = c(0, 1),
     main = "no autocorrelation, paired", xlab = "p-values", nclass = nc, 
     cex.axis = ca2, cex.lab = cl, cex.main = cm)
mtext("(q)", side = 3, at = 0, cex = ca, line = myline)
##
hist(sim0$ttestPvals[, "randomizedDesign"], probability = TRUE, xlim = c(0, 1), 
     main = "no autocorrelation, random", xlab = "p-values", nclass = nc, 
     cex.axis = ca2, cex.lab = cl, cex.main = cm)
mtext("(r)", side = 3, at = 0, cex = ca, line = myline)
#dev.off()

