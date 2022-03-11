
##########################################
## code for generating Fig 6
##########################################

library(sandwich)
library(lmtest)
library(forecast)


SimNullDataAR1 <- function(n, 
                           beta1,
                           design = c("randomized", "paired")) {
  y <- rep(NA, n)
  e <- rnorm(n)
  x <- switch(design,
              randomized = sample(c("b", "a"), n, replace = TRUE),
              paired = rep(c("b", "a"), n/2))
  y[1] <- e[1]
  for (i in 2:n) {
    y[i] <- beta1 * y[i-1] + e[i]
  }
  
  data.frame(y, x)
}


RunSimulationsUnderNull <- function(nSim, D, design, myseeds) {  
  Pvals <- matrix(NA, nSim, 3)
  colnames(Pvals) <- c("lm", "Newey-West", "arima") 
  Betas <- matrix(NA, nSim, 3)
  colnames(Betas) <- c("true", "lm", "arima") 
  acEstimates <- rep(NA, nSim)
  arimaModel <- matrix(NA, nSim, 3)
  colnames(arimaModel) <- c("p", "d", "q")
  
  for (i in seq(nSim)) {
    cat(i, "\n")
    n <- D[i, "n"]
    beta1 <- D[i, "beta1"]
    set.seed(myseeds[i])
    dat <- SimNullDataAR1(n, 
                          beta1,
                          design)
    
    ## estimate the autocorrelation at lag 1
    acEstimates[i] <- acf(dat$y, lag.max = 1, plot = FALSE)$acf[2, 1, 1]
    
    ## get t-test p-values
    fit1 <- lm(y ~ x, data = dat)
    su <- summary(fit1)
    Pvals[i, 1] <- su$coef[2, 4]
    Pvals[i, 2] <- coeftest(fit1, vcov = NeweyWest)[2, 4]
    
    xreg <- rep(NA, n)
    xreg[dat$x == "b"] <- 0
    xreg[dat$x == "a"] <- 1
    fit2 <- AutoArimaTest(x = dat$y, xreg = xreg)
    Pvals[i, 3] <- fit2$pval  
    
    arimaModel[i, "p"] <- fit2$p
    arimaModel[i, "d"] <- fit2$d
    arimaModel[i, "q"] <- fit2$q   
    
    Betas[i, 1] <- 0
    Betas[i, 2] <- su$coef[2, 1]
    Betas[i, 3] <- fit2$est
  } 
  
  list(Pvals = Pvals, acEstimates = acEstimates, arimaModel = arimaModel)
}


AutoArimaTest <- function(x, xreg, ic = "bic") {
  aux <- auto.arima(x = x, xreg = xreg, ic = ic, approximation = TRUE)
  auxOrder <- arimaorder(aux)
  p <- auxOrder[1]
  d <- auxOrder[2]
  q <- auxOrder[3]
  est <- aux$coef["xreg"]
  tstat <- abs(est)/sqrt(aux$var.coef["xreg", "xreg"])
  pval <- 2 * pt(tstat, df = length(x) - 2, lower.tail = FALSE)
  
  list(pval = pval, p = p, d = d, q = q, est = est)
}


TransformDesign <- function(x, par.ranges, is.discrete) {
  for (i in 1:ncol(x)) {
    par.range <- par.ranges[[i]]
    par.min <- min(par.range)
    par.max <- max(par.range)
    x[, i] <- x[, i] * (par.max - par.min) + par.min
    if (is.discrete[[i]]) {
      x[, i] <- round(x[, i])
    }
  }
  x
}


GetEmpiricalPower <- function(x, alphaGrid = seq(0.001, 1, by = 0.001)) {
  EmpiricalPower <- function(x, alphaGrid) {
    x <- x[!is.na(x)]
    ntests <- length(x)
    nalpha <- length(alphaGrid)
    out <- rep(NA, nalpha)
    for (i in seq(nalpha)) {
      out[i] <- sum(x <= alphaGrid[i])/ntests
    }
    out
  }
  lapply(x, EmpiricalPower, alphaGrid)
}


######################################
######################################
######################################

nSim <- 10000

set.seed(12345)
UD <- matrix(runif(nSim * 2), nSim, 2)


######################
## positive only
######################

par.ranges <- list(n = c(15, 400),  
                   beta1 = c(0, 0.9))
is.discrete <- c(TRUE, FALSE)
D <- TransformDesign(UD, par.ranges, is.discrete)
colnames(D) <- names(par.ranges)
D[, "n"] <- 2*D[, "n"] ## make sure n is even for paired design

set.seed(12345)
myseeds <- sample(10000:100000, nSim, replace = FALSE)

sim.N.R.arma10.P <- RunSimulationsUnderNull(nSim, D, design = "randomized", myseeds = myseeds)
sim.N.P.arma10.P <- RunSimulationsUnderNull(nSim, D, design = "paired", myseeds = myseeds)


######################
## negative only
######################

par.ranges <- list(n = c(15, 400),  
                   beta1 = c(-0.9, 0))
is.discrete <- c(TRUE, FALSE)
D <- TransformDesign(UD, par.ranges, is.discrete)
colnames(D) <- names(par.ranges)
D[, "n"] <- 2*D[, "n"] ## make sure n is even for paired design

set.seed(12345)
myseeds <- sample(10000:100000, nSim, replace = FALSE)

sim.N.R.arma10.N <- RunSimulationsUnderNull(nSim, D, design = "randomized", myseeds = myseeds)
sim.N.P.arma10.N <- RunSimulationsUnderNull(nSim, D, design = "paired", myseeds = myseeds)



######################
## no autocorrelation
######################

par.ranges <- list(n = c(15, 400),  
                   beta1 = c(0, 0))
is.discrete <- c(TRUE, FALSE)
D <- TransformDesign(UD, par.ranges, is.discrete)
colnames(D) <- names(par.ranges)
D[, "n"] <- 2*D[, "n"] ## make sure n is even for paired design

set.seed(12345)
myseeds <- sample(10000:100000, nSim, replace = FALSE)

sim.N.R.arma10.0 <- RunSimulationsUnderNull(nSim, D, design = "randomized", myseeds = myseeds)
sim.N.P.arma10.0 <- RunSimulationsUnderNull(nSim, D, design = "paired", myseeds = myseeds)




#####################################################
#####################################################
#####################################################


GetEmpiricalPower <- function(x, alphaGrid = seq(0.001, 1, by = 0.001)) {
  EmpiricalPower <- function(x, alphaGrid) {
    x <- x[!is.na(x)]
    ntests <- length(x)
    nalpha <- length(alphaGrid)
    out <- rep(NA, nalpha)
    for (i in seq(nalpha)) {
      out[i] <- sum(x <= alphaGrid[i])/ntests
    }
    out
  }
  lapply(x, EmpiricalPower, alphaGrid)
}




alphaGrid <- seq(0.001, 1, by = 0.001)
mylwd <- 1
cexleg <- 1.5
cexl <- 1.5
cl <- 1.4
ca <- 1.4
cm <- 1.5

figpath <- ""

## Fig 6
#pdf(paste(figpath, "null_sim_study_lm_nw_arima_plos_new.pdf", sep = ""), width = 8, height = 6.5)
par(mfrow = c(2, 3), mar = c(4, 3, 3, 1), mgp = c(2, 0.5, 0))
mymain <- "negative autocor., paired"
sim <- sim.N.P.arma10.N
ep <- GetEmpiricalPower(x = list(lm = sim$Pvals[, 1], 
                                 nw = sim$Pvals[, 2], 
                                 ar = sim$Pvals[, 3]), alphaGrid)
plot(alphaGrid, ep[[1]], type = "n", ylim = c(0, 1),
     ylab = "empirical type I error", xlab = "nominal significance level", cex.lab = cl,
     cex.main = cm, cex.axis = ca, main = mymain)
abline(a = 0, b = 1, col = "grey", lwd = 2)
lines(alphaGrid, ep[[1]], col = "darkgreen", lwd = mylwd)
lines(alphaGrid, ep[[2]], col = "blue", lwd = mylwd)
lines(alphaGrid, ep[[3]], col = "red", lwd = mylwd)
legend("bottomright", legend = c("linear regression", "Newey-West HAC", "ARIMA"), 
       text.col = c("darkgreen", "blue", "red"), bty = "n", cex = cexleg)
mtext("(a)", side = 3, at = -0.1, line = 1, cex = cl)
#############
mymain <- "positive autocor., paired"
sim <- sim.N.P.arma10.P
ep <- GetEmpiricalPower(x = list(lm = sim$Pvals[, 1], 
                                 nw = sim$Pvals[, 2], 
                                 ar = sim$Pvals[, 3]), alphaGrid)
plot(alphaGrid, ep[[1]], type = "n", ylim = c(0, 1),
     ylab = "empirical type I error", xlab = "nominal significance level", cex.lab = cl,
     cex.main = cm, cex.axis = ca, main = mymain)
abline(a = 0, b = 1, col = "grey", lwd = 2)
lines(alphaGrid, ep[[1]], col = "darkgreen", lwd = mylwd)
lines(alphaGrid, ep[[2]], col = "blue", lwd = mylwd)
lines(alphaGrid, ep[[3]], col = "red", lwd = mylwd)
legend("topleft", legend = c("linear regression", "Newey-West HAC", "ARIMA"), 
       text.col = c("darkgreen", "blue", "red"), bty = "n", cex = cexleg)
mtext("(b)", side = 3, at = -0.1, line = 1, cex = cl)
#############
mymain <- "no autocor., paired"
sim <- sim.N.P.arma10.0
ep <- GetEmpiricalPower(x = list(lm = sim$Pvals[, 1], 
                                 nw = sim$Pvals[, 2], 
                                 ar = sim$Pvals[, 3]), alphaGrid)
plot(alphaGrid, ep[[1]], type = "n", ylim = c(0, 1),
     ylab = "empirical type I error", xlab = "nominal significance level", cex.lab = cl,
     cex.main = cm, cex.axis = ca, main = mymain)
abline(a = 0, b = 1, col = "grey", lwd = 2)
lines(alphaGrid, ep[[1]], col = "darkgreen", lwd = mylwd)
lines(alphaGrid, ep[[2]], col = "blue", lwd = mylwd)
lines(alphaGrid, ep[[3]], col = "red", lwd = mylwd)
legend("bottomright", legend = c("linear regression", "Newey-West HAC", "ARIMA"), 
       text.col = c("darkgreen", "blue", "red"), bty = "n", cex = cexleg)
mtext("(c)", side = 3, at = -0.1, line = 1, cex = cl)
#############
mymain <- "negative autocor., random"
sim <- sim.N.R.arma10.N
ep <- GetEmpiricalPower(x = list(lm = sim$Pvals[, 1], 
                                 nw = sim$Pvals[, 2], 
                                 ar = sim$Pvals[, 3]), alphaGrid)
plot(alphaGrid, ep[[1]], type = "n", ylim = c(0, 1),
     ylab = "empirical type I error", xlab = "nominal significance level", cex.lab = cl,
     cex.main = cm, cex.axis = ca, main = mymain)
abline(a = 0, b = 1, col = "grey", lwd = 2)
lines(alphaGrid, ep[[1]], col = "darkgreen", lwd = mylwd)
lines(alphaGrid, ep[[2]], col = "blue", lwd = mylwd)
lines(alphaGrid, ep[[3]], col = "red", lwd = mylwd)
legend("bottomright", legend = c("linear regression", "Newey-West HAC", "ARIMA"), 
       text.col = c("darkgreen", "blue", "red"), bty = "n", cex = cexleg)
mtext("(d)", side = 3, at = -0.1, line = 1, cex = cl)
#############
mymain <- "positive autocor., random"
sim <- sim.N.R.arma10.P
ep <- GetEmpiricalPower(x = list(lm = sim$Pvals[, 1], 
                                 nw = sim$Pvals[, 2], 
                                 ar = sim$Pvals[, 3]), alphaGrid)
plot(alphaGrid, ep[[1]], type = "n", ylim = c(0, 1),
     ylab = "empirical type I error", xlab = "nominal significance level", cex.lab = cl,
     cex.main = cm, cex.axis = ca, main = mymain)
abline(a = 0, b = 1, col = "grey", lwd = 2)
lines(alphaGrid, ep[[1]], col = "darkgreen", lwd = mylwd)
lines(alphaGrid, ep[[2]], col = "blue", lwd = mylwd)
lines(alphaGrid, ep[[3]], col = "red", lwd = mylwd)
legend("bottomright", legend = c("linear regression", "Newey-West HAC", "ARIMA"), 
       text.col = c("darkgreen", "blue", "red"), bty = "n", cex = cexleg)
mtext("(e)", side = 3, at = -0.1, line = 1, cex = cl)
#############
mymain <- "no autocor., random"
sim <- sim.N.R.arma10.0
ep <- GetEmpiricalPower(x = list(lm = sim$Pvals[, 1], 
                                 nw = sim$Pvals[, 2], 
                                 ar = sim$Pvals[, 3]), alphaGrid)
plot(alphaGrid, ep[[1]], type = "n", ylim = c(0, 1),
     ylab = "empirical type I error", xlab = "nominal significance level", cex.lab = cl,
     cex.main = cm, cex.axis = ca, main = mymain)
abline(a = 0, b = 1, col = "grey", lwd = 2)
lines(alphaGrid, ep[[1]], col = "darkgreen", lwd = mylwd)
lines(alphaGrid, ep[[2]], col = "blue", lwd = mylwd)
lines(alphaGrid, ep[[3]], col = "red", lwd = mylwd)
legend("bottomright", legend = c("linear regression", "Newey-West HAC", "ARIMA"), 
       text.col = c("darkgreen", "blue", "red"), bty = "n", cex = cexleg)
mtext("(f)", side = 3, at = -0.1, line = 1, cex = cl)
#############
#dev.off()

