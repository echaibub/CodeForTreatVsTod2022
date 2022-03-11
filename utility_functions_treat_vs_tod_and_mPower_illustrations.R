
library(pROC)
library(forecast)
library(gplots)
library(sandwich)
library(lmtest)


################################################
## treat vs tod functions
################################################

SummarizeEvidence <- function(x, thr = 0.05) {
  evidence <- apply(x, 1, EvidenceRule, thr)
  
  data.frame(x, evidence, check.names = FALSE)
}


EvidenceRule <- function(pvals, thr) {
  ## 1: a(T,X)=0, 2: a(Y,X)=0, 3: a(Y,T)=0, 4: a(Y,X|T)=0, 5: a(Y,T|X)=0
  evidence <- NA
  if (sum(is.na(pvals)) == 0) {
    if(pvals[2] <= thr & pvals[3] <= thr & pvals[4] <= thr & pvals[5] <= thr) {
      evidence <- "treat and tod effects"
    }
    if(pvals[2] <= thr & pvals[3] <= thr & pvals[4] <= thr & pvals[5] > thr) {
      evidence <- "treat effect"
    }
    if(pvals[2] <= thr & pvals[3] > thr & pvals[4] <= thr & pvals[5] > thr) {
      evidence <- "treat effect"
    }
    if(pvals[2] <= thr & pvals[3] <= thr & pvals[4] > thr & pvals[5] <= thr) {
      evidence <- "tod effect"
    }
    if(pvals[2] > thr & pvals[3] <= thr & pvals[4] > thr & pvals[5] <= thr) {
      evidence <- "tod effect"
    }
    if(pvals[2] <= thr & pvals[3] > thr & pvals[4] > thr & pvals[5] > thr) {
      evidence <- "inconsistent treat effect"
    }
    if(pvals[2] <= thr & pvals[3] <= thr & pvals[4] > thr & pvals[5] > thr) {
      evidence <- "inconsistent treat and tod effect"
    }
    if(pvals[2] > thr & pvals[3] <= thr & pvals[4] > thr & pvals[5] > thr) {
      evidence <- "inconsistent tod effect"
    }
  }
  
  evidence
}


#####################################
## UI-tests functions
#####################################

ParticipantTreatVsTodEffectsLmNw <- function(dat, id, featNames) {
  dat <- dat[dat$healthCode == id,]
  nfeat <- length(featNames)
  
  lmBetas <- data.frame(matrix(NA, nfeat, 5))
  rownames(lmBetas) <- featNames
  colnames(lmBetas) <- c("b.TX", "b.YX", "b.YT", "b.YX|T", "b.YT|X")  
  lmPvals <- data.frame(matrix(NA, nfeat, 5))
  rownames(lmPvals) <- featNames
  colnames(lmPvals) <- c("a(T,X)=0", "a(Y,X)=0", "a(Y,T)=0", "a(Y,X|T)=0", "a(Y,T|X)=0") 
  nwPvals <- data.frame(matrix(NA, nfeat, 5))
  rownames(nwPvals) <- featNames
  colnames(nwPvals) <- c("a(T,X)=0", "a(Y,X)=0", "a(Y,T)=0", "a(Y,X|T)=0", "a(Y,T|X)=0") 
  
  lmResAcf <- data.frame(matrix(NA, nfeat, 4))
  rownames(lmResAcf) <- featNames
  colnames(lmResAcf) <- c("lag1_Y~X", "lag1_Y~X+T",  "pval_Y~X", "pval_Y~X+T") 
  
  x <- dat$medTimepoint
  tod <- dat$tod 
  fit1 <- lm(tod ~ x)
  su1 <- summary(fit1)$coefficients 
  lmBetas[, "b.TX"] <- su1[2, 1]
  lmPvals[, "a(T,X)=0"] <- su1[2, 4]
  nwPvals[, "a(T,X)=0"] <- coeftest(fit1, vcov = NeweyWest)[2, 4]
  for (i in seq(nfeat)) {
    #cat(i, "\\n")
    y <- dat[, featNames[i]]
    ## lm fits
    fit2 <- lm(y ~ x)
    fit3 <- lm(y ~ tod)
    fit4 <- lm(y ~ x + tod)    
    su2 <- summary(fit2)$coefficients
    su3 <- summary(fit3)$coefficients   
    su4 <- summary(fit4)$coefficients    
    lmBetas[i, "b.YX"] <- su2[2, 1]
    lmBetas[i, "b.YT"] <- su3[2, 1]   
    lmBetas[i, "b.YX|T"] <- su4[2, 1]
    lmBetas[i, "b.YT|X"] <- su4[3, 1]
    lmPvals[i, "a(Y,X)=0"] <- su2[2, 4]
    lmPvals[i, "a(Y,T)=0"] <- su3[2, 4]    
    lmPvals[i, "a(Y,X|T)=0"] <- su4[2, 4]
    lmPvals[i, "a(Y,T|X)=0"] <- su4[3, 4]
    ## lm with Newey-West covariance estimation
    nwPvals[i, "a(Y,X)=0"] <- coeftest(fit2, vcov = NeweyWest)[2, 4]
    nwPvals[i, "a(Y,T)=0"] <- coeftest(fit3, vcov = NeweyWest)[2, 4]
    nwPvals[i, "a(Y,X|T)=0"] <- coeftest(fit4, vcov = NeweyWest)[2, 4]
    nwPvals[i, "a(Y,T|X)=0"] <- coeftest(fit4, vcov = NeweyWest)[3, 4]
    
    ## get residuals from lm fits
    res2 <- fit2$resid
    lmResAcf[i, "lag1_Y~X"] <- acf(res2, plot = FALSE)$acf[2]
    lmResAcf[i, "pval_Y~X"] <- Box.test(res2, lag = 1, type = "Ljung-Box", fitdf = 0)$p.value
    res4 <- fit4$resid
    lmResAcf[i, "lag1_Y~X+T"] <- acf(res4, plot = FALSE)$acf[2]
    lmResAcf[i, "pval_Y~X+T"] <- Box.test(res4, lag = 1, type = "Ljung-Box", fitdf = 0)$p.value
  }
  
  list(lmBetas = lmBetas, 
       lmPvals = lmPvals, 
       nwPvals = nwPvals,
       lmResAcf = lmResAcf)
}



ParticipantTreatVsTodEffectsArima <- function(dat, id, featNames) {
  arimaAnalysisTest1 <- function(y, xreg, lag = 5, ic = "bic", max.p = 10, max.q = 10) {
    aux <- auto.arima(x = y, xreg = xreg, ic = ic, approximation = TRUE, max.p = max.p, max.q = max.q)
    auxOrder <- arimaorder(aux)
    p <- auxOrder[1]
    d <- auxOrder[2]
    q <- auxOrder[3]
    eff <- aux$coef["xreg"]
    zstat <- abs(eff)/sqrt(aux$var.coef["xreg", "xreg"])
    pval <- 2 * pnorm(zstat, lower.tail = FALSE)
    list(eff = eff, pval = pval, p = p, d = d, q = q, autoArimaFit = aux)
  }
  arimaAnalysisTest2 <- function(y, xreg, lag = 5, ic = "bic", max.p = 10, max.q = 10) {
    aux <- auto.arima(x = y, xreg = xreg, ic = ic, approximation = TRUE, max.p = max.p, max.q = max.q)
    auxOrder <- arimaorder(aux)
    p <- auxOrder[1]
    d <- auxOrder[2]
    q <- auxOrder[3]
    ####
    eff1 <- aux$coef["xreg1"]
    zstat1 <- abs(eff1)/sqrt(aux$var.coef["xreg1", "xreg1"])
    pval1 <- 2 * pnorm(zstat1, lower.tail = FALSE)
    ####
    eff2 <- aux$coef["xreg2"]
    zstat2 <- abs(eff2)/sqrt(aux$var.coef["xreg2", "xreg2"])
    pval2 <- 2 * pnorm(zstat2, lower.tail = FALSE)  
    
    list(eff1 = eff1, pval1 = pval1, eff2 = eff2, pval2 = pval2, 
         p = p, d = d, q = q, autoArimaFit = aux)
  }
  dat <- dat[dat$healthCode == id,]
  nfeat <- length(featNames)
  
  arimaBetas <- data.frame(matrix(NA, nfeat, 5))
  rownames(arimaBetas) <- featNames
  colnames(arimaBetas) <- c("b.TX", "b.YX", "b.YT", "b.YX|T", "b.YT|X")  
  arimaPvals <- data.frame(matrix(NA, nfeat, 5))
  rownames(arimaPvals) <- featNames
  colnames(arimaPvals) <- c("a(T,X)=0", "a(Y,X)=0", "a(Y,T)=0", "a(Y,X|T)=0", "a(Y,T|X)=0") 
  arimaResAcf <- data.frame(matrix(NA, nfeat, 4))
  rownames(arimaResAcf) <- featNames
  colnames(arimaResAcf) <- c("lag1_Y~X", "lag1_Y~X+T",  "pval_Y~X", "pval_Y~X+T")
  
  arimaModel <- data.frame(matrix(NA, nfeat, 6))
  rownames(arimaModel) <- featNames
  colnames(arimaModel) <- c("p_Y~X", "d_Y~X", "q_Y~X", "p_Y~X+T", "d_Y~X+T", "q_Y~X+T") 
  
  x <- dat$medTimepoint
  tod <- dat$tod 
  
  afit1 <- try(arimaAnalysisTest1(y = tod, xreg = as.numeric(x) - 1, 
                                  lag = lag, ic = "bic", max.p = 10, max.q = 10), silent = TRUE)
  if (!inherits(afit1, "try-error")) {
    arimaBetas[, "b.TX"] <- afit1$eff
    arimaPvals[, "a(T,X)=0"] <- afit1$pval 
  }
  
  for (i in seq(nfeat)) {
    #cat(i, "\n")
    y <- dat[, featNames[i]]
    
    ## arima fits
    afit2 <- try(arimaAnalysisTest1(y = y, xreg = as.numeric(x) - 1), silent = TRUE)
    if (!inherits(afit2, "try-error")) {
      arimaBetas[i, "b.YX"] <- afit2$eff
      arimaPvals[i, "a(Y,X)=0"] <- afit2$pval
      arimaModel[i, "p_Y~X"] <- afit2$p
      arimaModel[i, "d_Y~X"] <- afit2$d
      arimaModel[i, "q_Y~X"] <- afit2$q 
      res2 <- afit2$autoArimaFit$resid
      res2 <- res2[!is.na(res2)]
      arimaResAcf[i, "lag1_Y~X"] <- acf(res2, plot = FALSE)$acf[2]
      arimaResAcf[i, "pval_Y~X"] <- Box.test(res2, lag = 1, type = "Ljung-Box", fitdf = 0)$p.value    
    }    
    afit3 <- try(arimaAnalysisTest1(y = y, xreg = tod), silent = TRUE)  
    if (!inherits(afit3, "try-error")) {
      arimaBetas[i, "b.YT"] <- afit3$eff 
      arimaPvals[i, "a(Y,T)=0"] <- afit3$pval 
    }
    xreg <- cbind(as.numeric(x) - 1, tod)
    colnames(xreg) <- NULL
    afit4 <- try(arimaAnalysisTest2(y = y, xreg = xreg), silent = TRUE)
    if (!inherits(afit4, "try-error")) {
      arimaBetas[i, "b.YX|T"] <- afit4$eff1
      arimaBetas[i, "b.YT|X"] <- afit4$eff2 
      arimaPvals[i, "a(Y,X|T)=0"] <- afit4$pval1
      arimaPvals[i, "a(Y,T|X)=0"] <- afit4$pval2
      arimaModel[i, "p_Y~X+T"] <- afit4$p
      arimaModel[i, "d_Y~X+T"] <- afit4$d
      arimaModel[i, "q_Y~X+T"] <- afit4$q  
      res4 <- afit4$autoArimaFit$resid
      res4 <- res4[!is.na(res4)]
      arimaResAcf[i, "lag1_Y~X+T"] <- acf(res4, plot = FALSE)$acf[2]
      arimaResAcf[i, "pval_Y~X+T"] <- Box.test(res4, lag = 1, type = "Ljung-Box", fitdf = 0)$p.value
    }
  }
  
  list(arimaBetas = arimaBetas,
       arimaPvals = arimaPvals, 
       arimaResAcf = arimaResAcf,
       arimaModel = arimaModel)
}


RunTreatmentVsTodEffectsLmNw <- function(dat, featNames) {
  ids <- as.character(unique(dat$healthCode))
  nids <- length(ids) 
  lmPvals <- vector(mode = "list", length = nids)
  names(lmPvals) <- ids
  nwPvals <- lmPvals
  lmBetas <- lmPvals
  lmResAcf <- lmPvals
  for (i in seq(nids)) {
    cat("participant: ", i, "\n")
    aux <- ParticipantTreatVsTodEffectsLmNw(dat, ids[i], featNames)
    lmBetas[[i]] <- aux$lmBetas
    lmPvals[[i]] <- aux$lmPvals
    nwPvals[[i]] <- aux$nwPvals
    lmResAcf[[i]] <- aux$lmResAcf
  }
  
  list(lmBetas = lmBetas, lmPvals = lmPvals, nwPvals = nwPvals, lmResAcf = lmResAcf)
}


RunTreatmentVsTodEffectsArima <- function(dat, featNames) {
  ids <- as.character(unique(dat$healthCode))
  nids <- length(ids) 
  arimaPvals <- vector(mode = "list", length = nids)
  names(arimaPvals) <- ids
  arimaBetas <- arimaPvals
  arimaResAcf <- arimaPvals
  arimaModel <- arimaPvals
  for (i in seq(nids)) {
    cat("participant: ", i, "\n")
    aux <- ParticipantTreatVsTodEffectsArima(dat, ids[i], featNames)
    arimaBetas[[i]] <- aux$arimaBetas
    arimaPvals[[i]] <- aux$arimaPvals
    arimaResAcf[[i]] <- aux$arimaResAcf
    arimaModel[[i]] <- aux$arimaModel
  }
  
  list(arimaBetas = arimaBetas, 
       arimaPvals = arimaPvals,        
       arimaResAcf = arimaResAcf,
       arimaModel = arimaModel)
}


ParticipantEvidence <- function(x, mtMethod = "BH", thr = 0.05, numericToo = TRUE) {
  nfeat <- nrow(x)
  ax <- PvalAdjustment(x, mtMethod, byrow = TRUE)
  ev <- as.character(SummarizeEvidence(ax, thr)[, "evidence"])
  evN <- NULL
  if (numericToo) {
    evN <- matrix(NA, nfeat, 1)
    rownames(evN) <- rownames(x)
    evN[ev == "tod effect"] <- -10
    evN[ev == "treat effect"] <- 10
    evN[ev == "treat and tod effects"] <- 0  
    evN[ev == "inconsistent treat effect"] <- 5 
    evN[ev == "inconsistent tod effect"] <- -5
    evN[ev == "inconsistent treat and tod effect"] <- 1 
  }
  
  list(ev = ev, evN = evN)
}


ParticipantPutativeTreatPval <- function(x, mtMethod1 = "none", mtMethod2 = "BH", thr = 0.05) {
  ev <- ParticipantEvidence(x, mtMethod1, thr, FALSE)$ev
  ## only adjust for tod if we detect a consistent 
  ## tod (or tod and treat) effect
  ev[ev == "treat effect"] <- NA 
  ev[ev == "inconsistent treat effect"] <- NA 
  ev[ev == "inconsistent treat and tod effect"] <- NA 
  ev[ev == "inconsistent tod effect"] <- NA 
  nfeat <- length(ev)
  pvals <- rep(NA, nfeat)
  names(pvals) <- names(ev)
  idxMarg <- which(is.na(ev))
  idxCond <- which(!is.na(ev))
  pvals[idxMarg] <- x[idxMarg, "a(Y,X)=0"]
  pvals[idxCond] <- x[idxCond, "a(Y,X|T)=0"]
  pvals <- p.adjust(pvals, method = mtMethod2)
  
  min(pvals)
}


ParticipantPutativeTodPval <- function(x, mtMethod1 = "none", mtMethod2 = "BH", thr = 0.05) {
  ev <- ParticipantEvidence(x, mtMethod1, thr, FALSE)$ev
  ## only adjust for treat if we detect a consistent 
  ## treat (or tod and treat) effect
  ev[ev == "tod effect"] <- NA 
  ev[ev == "inconsistent tod effect"] <- NA 
  ev[ev == "inconsistent treat and tod effect"] <- NA 
  ev[ev == "inconsistent treat effect"] <- NA 
  nfeat <- length(ev)
  pvals <- rep(NA, nfeat)
  names(pvals) <- names(ev)
  idxMarg <- which(is.na(ev))
  idxCond <- which(!is.na(ev))
  pvals[idxMarg] <- x[idxMarg, "a(Y,T)=0"]
  pvals[idxCond] <- x[idxCond, "a(Y,T|X)=0"]
  pvals <- p.adjust(pvals, method = mtMethod2)
  
  min(pvals)
}


PvalAdjustment <- function(x, mtMethod, byrow = TRUE) {
  if (byrow) {
    ax <- t(apply(x, 1, p.adjust, mtMethod))
  }
  else {
    ax <- apply(x, 2, p.adjust, mtMethod)
  }
  
  ax
}


ParticipantAcf <- function(pvals, acfs, mtMethod1, thr, 
                           margName = "lag1_Y~X", 
                           condName = "lag1_Y~X+T") {
  ev <- ParticipantEvidence(pvals, mtMethod1, thr, FALSE)$ev
  ## only adjust for tod if we detect a consistent 
  ## tod (or tod and treat) effect
  ev[ev == "treat effect"] <- NA 
  ev[ev == "inconsistent treat effect"] <- NA 
  ev[ev == "inconsistent treat and tod effect"] <- NA 
  ev[ev == "inconsistent tod effect"] <- NA 
  nfeat <- length(ev)
  selacfs <- rep(NA, nfeat)
  sellbpvals <- rep(NA, nfeat)
  names(selacfs) <- names(ev)
  names(sellbpvals) <- names(ev)
  idxMarg <- which(is.na(ev))
  idxCond <- which(!is.na(ev))
  selacfs[idxMarg] <- acfs[idxMarg, "lag1_Y~X"]
  selacfs[idxCond] <- acfs[idxCond, "lag1_Y~X+T"]
  sellbpvals[idxMarg] <- acfs[idxMarg, "pval_Y~X"]
  sellbpvals[idxCond] <- acfs[idxCond, "pval_Y~X+T"]
  
  list(selacfs = selacfs, sellbpvals = sellbpvals)
}


OrganizeAcfs <- function(pvalsL, acfsL, mtMethod1 = "BH", thr = 0.05) {
  ids <- names(pvalsL)
  nids <- length(ids)
  nfeat <- nrow(pvalsL[[1]])
  acfM <- matrix(NA, nfeat, nids)
  colnames(acfM) <- ids
  rownames(acfM) <- rownames(pvalsL[[1]])
  lboxM <- acfM
  for (i in seq(nids)) {
    aux <- ParticipantAcf(pvalsL[[i]], acfsL[[i]], mtMethod1, thr) 
    acfM[,i] <- aux$selacfs
    lboxM[,i] <- aux$sellbpvals
  }
  
  list(acfM = acfM, lboxM = lboxM)
}


GenerateEvidenceMatrix <- function(xL, thr = 0.05, numericToo = TRUE) {
  AdjustedParticipantEvidence <- function(x, nids, nfeat, thr = 0.05, numericToo = TRUE) {
    ax <- PvalAdjustmentFeaturesAndHealthCode(x, nids, nfeat)
    ev <- as.character(SummarizeEvidence(ax, thr)[, "evidence"])
    evN <- NULL
    if (numericToo) {
      evN <- matrix(NA, nfeat, 1)
      rownames(evN) <- rownames(x)
      evN[ev == "tod effect"] <- -10
      evN[ev == "treat effect"] <- 10
      evN[ev == "treat and tod effects"] <- 0  
      evN[ev == "inconsistent treat effect"] <- 5 
      evN[ev == "inconsistent tod effect"] <- -5
      evN[ev == "inconsistent treat and tod effect"] <- 1 
    }
    list(ev = ev, evN = evN)
  }
  PvalAdjustmentFeaturesAndHealthCode <- function(x, nids, nfeat) {
    ntests <- nids * nfeat
    x <- x * ntests
    x[x > 1] <- 1
    x
  }
  ids <- names(xL)
  nids <- length(ids)
  nfeat <- nrow(xL[[1]])
  ev <- matrix(NA, nfeat, nids)
  rownames(ev) <- rownames(xL[[1]])
  colnames(ev) <- ids
  if (numericToo) {
    evN <- ev
    for (i in seq(nids)) {
      aux <- AdjustedParticipantEvidence(xL[[i]], nids, nfeat, thr, numericToo)
      ev[, i] <- aux$ev
      evN[, i] <- aux$evN
    }
  }
  else {
    for (i in seq(nids)) {
      ev[, i] <- AdjustedParticipantEvidence(xL[[i]], nids, nfeat, thr, numericToo)$ev
    }    
  }
  
  list(ev = ev, evN = evN)
}


UItestsPutativeTreat <- function(xL, mtMethod1 = "none", mtMethod2 = "BH", thr = 0.05) {
  ids <- names(xL)
  nids <- length(ids)
  nfeat <- nrow(xL[[1]])
  pvals <- matrix(NA, nids, 1)
  rownames(pvals) <- ids
  colnames(pvals) <- "pvalUItest"
  for (i in seq(nids)) {
    pvals[i, 1] <- ParticipantPutativeTreatPval(xL[[i]], mtMethod1, mtMethod2, thr)
  }
  
  pvals
}


UItestsPutativeTod <- function(xL, mtMethod1 = "none", mtMethod2 = "BH", thr = 0.05) {
  ids <- names(xL)
  nids <- length(ids)
  nfeat <- nrow(xL[[1]])
  pvals <- matrix(NA, nids, 1)
  rownames(pvals) <- ids
  colnames(pvals) <- "pvalUItest"
  for (i in seq(nids)) {
    pvals[i, 1] <- ParticipantPutativeTodPval(xL[[i]], mtMethod1, mtMethod2, thr)
  }
  
  pvals
}


GetTreatPvalMatrix <- function(pvalList, mtMethod1 = "BH", thr = 0.05) {
  ids <- names(pvalList)
  nids <- length(ids)
  nfeat <- nrow(pvalList[[1]])
  pvalM <- matrix(NA, nids, nfeat)
  rownames(pvalM) <- ids
  colnames(pvalM) <- rownames(pvalList[[1]])
  for (i in seq(nids)) {
    ax <- PvalAdjustment(pvalList[[i]], mtMethod1, byrow = TRUE)
    ev <- as.character(SummarizeEvidence(ax, thr)[, "evidence"])
    ev[ev == "treat effect"] <- NA 
    ev[ev == "inconsistent treat effect"] <- NA 
    ev[ev == "inconsistent treat and tod effect"] <- NA 
    ev[ev == "inconsistent tod effect"] <- NA 
    idxMarg <- which(is.na(ev))
    idxCond <- which(!is.na(ev))
    pvalM[i, idxMarg] <- pvalList[[i]][idxMarg, "a(Y,X)=0"]
    pvalM[i, idxCond] <- pvalList[[i]][idxCond, "a(Y,X|T)=0"]
  }
  
  pvalM
}


GetTodPvalMatrix <- function(pvalList, mtMethod1 = "BH", thr = 0.05) {
  ids <- names(pvalList)
  nids <- length(ids)
  nfeat <- nrow(pvalList[[1]])
  pvalM <- matrix(NA, nids, nfeat)
  rownames(pvalM) <- ids
  colnames(pvalM) <- rownames(pvalList[[1]])
  for (i in seq(nids)) {
    ax <- PvalAdjustment(pvalList[[i]], mtMethod1, byrow = TRUE)
    ev <- as.character(SummarizeEvidence(ax, thr)[, "evidence"])
    ev[ev == "tod effect"] <- NA 
    ev[ev == "inconsistent tod effect"] <- NA 
    ev[ev == "inconsistent treat and tod effect"] <- NA 
    ev[ev == "inconsistent treat effect"] <- NA 
    idxMarg <- which(is.na(ev))
    idxCond <- which(!is.na(ev))
    pvalM[i, idxMarg] <- pvalList[[i]][idxMarg, "a(Y,T)=0"]
    pvalM[i, idxCond] <- pvalList[[i]][idxCond, "a(Y,T|X)=0"]
  }
  
  pvalM
}


############################################
## additional functions
############################################

GetDataForNof1 <- function(dat, beforeThr, afterThr) {
  CountBeforeAfterMedicationPerParticipant <- function(dat) {
    participantIds <- unique(dat$healthCode)
    nParticipants <- length(participantIds)
    out <- data.frame(matrix(NA, nParticipants, 4))
    colnames(out) <- c("healthCode", "n", "nBefore", "nAfter")
    out[, 1] <- participantIds
    for (i in seq(nParticipants)) {
      #cat(i, "\n")
      pdat <- dat[which(dat$healthCode == participantIds[i]),]
      aux <- as.character(pdat$medTimepoint)
      out[i, "n"] <- nrow(pdat)
      out[i, "nBefore"] <- sum(aux == "Immediately before Parkinson medication", na.rm = TRUE)
      out[i, "nAfter"] <- sum(aux == "Just after Parkinson medication (at your best)", na.rm = TRUE)
    }
    out[order(out[, "nBefore"] + out[, "nAfter"], decreasing = TRUE),]
  }
  GetParticipants <- function(x, beforeThr, afterThr) {
    idx <- which(x$nBefore >= beforeThr & x$nAfter >= afterThr)
    x$healthCode[idx]
  }
  ## Get data from (self-reported) Parkinson's patients
  dat <- dat[dat$PD == TRUE,]
  ## Keep only the "before medication" and "after medication" data
  ## (drops the "at another time" data, as well as, the data missing 
  ## the medTimepoint response)
  idxBefore <- which(dat$medTimepoint == "Immediately before Parkinson medication")
  idxAfter <- which(dat$medTimepoint == "Just after Parkinson medication (at your best)")
  dat <- dat[sort(c(idxBefore, idxAfter)),]
  ## Make sure that there are only two levels for 
  ## the medTimepoint factor
  dat$medTimepoint <- factor(dat$medTimepoint)
  ## Count the number of activity tasks performed before and
  ## after medication
  countBA <- CountBeforeAfterMedicationPerParticipant(dat)
  ## Select participants that performed at least "beforeThr" 
  ## activity tasks before medication, and at least "afterThr" 
  ## tasks after medication
  tokeep <- GetParticipants(x = countBA, beforeThr = beforeThr, afterThr = afterThr)
  dat <- dat[dat$healthCode %in% tokeep,]  
  ## Replace infinite values by NA
  dat[sapply(dat, is.infinite)] <- NA
  
  dat
}


IncludeUTCandLocalTimeVariables <- function(dat, tzDat) {
  GetDayAndTimeOfDay <- function(x) {
    n <- length(x)
    tod <- rep(NA, n)
    day <- rep(NA, n)
    idx <- which(!is.na(x))
    x <- x[!is.na(x)]
    x <- as.character(x)
    n <- length(x)
    aux <- unlist(strsplit(x, " "))
    if (length(aux) == (2*n)) {
      auxDay <- aux[seq(1, 2*n-1, by = 2)]
      auxTime <- aux[seq(2, 2*n, by = 2)]
      auxTime <- as.numeric(unlist(strsplit(auxTime, ":")))
      ih <- seq(1, 3*n, by = 3)
      im <- seq(2, 3*n, by = 3)
      is <- seq(3, 3*n, by = 3)
      secs <- 3600 * auxTime[ih] + 60 * auxTime[im] + auxTime[is]
      tod[idx] <- secs/3600
      day[idx] <- unlist(auxDay) 
    }
    if (length(aux) == (3*n)) {
      auxDay <- aux[seq(1, 3*n-2, by = 3)]
      auxTime <- aux[seq(2, 3*n-1, by = 3)]
      auxTime <- as.numeric(unlist(strsplit(auxTime, ":")))
      ih <- seq(1, 3*n, by = 3)
      im <- seq(2, 3*n, by = 3)
      is <- seq(3, 3*n, by = 3)
      secs <- 3600 * auxTime[ih] + 60 * auxTime[im] + auxTime[is]
      tod[idx] <- secs/3600
      day[idx] <- unlist(auxDay) 
    }
    list(tod = tod, day = day)
  }
  ## rename the original createdOn time stamp
  ## as the createdOnUTC
  dat$createdOnUTC <- as.POSIXct(dat$createdOn, tz = "UTC")
  ## create day and time of the day variables in UTC
  auxTime <- GetDayAndTimeOfDay(dat$createdOnUTC)
  dat$todUTC <- auxTime$tod
  dat$dayUTC <- as.POSIXct(auxTime$day, tz = "UTC")
  ## get the day and tod values in local time
  ids <- as.character(unique(dat$healthCode))
  aux <- tzDat[match(ids, tzDat$healthCode), c("healthCode", "UTC_offset")]
  utcOffset <- rep(NA, nrow(dat))
  for (i in seq(length(ids))) {
    utcOffset[which(dat$healthCode %in% aux[i, "healthCode"])] <- aux[i, "UTC_offset"]
  }
  tod <- dat$todUTC + utcOffset
  ## negative values correspond to tasks done later in the
  ## day, that show up as done earlier in the next day in
  ## UTC time
  nextday <- which(tod < 0)  
  ## fix the local tod
  tod[nextday] <- 24 + tod[nextday]
  ## fix the local day (go one day back)
  day <- dat$dayUTC
  day[nextday] <- day[nextday] - (24*60*60)
  ## include the local time variables 
  ## (it still prints as UTC, though)
  dat$tod <- tod
  dat$day <- day
  dat$createdOn <- as.POSIXct(day + tod * 3600)
  dat 
}


TransformFeatures <- function(dat, featNames) {
  NormalTrans <- function(x) {
    n <- sum(!is.na(x))
    r <- rank(x, na.last = "keep")
    qnorm((r - 0.5)/n)
  }
  ids <- as.character(unique(dat$healthCode))
  nids <- length(ids)
  nfeat <- length(featNames)
  for (i in seq(nids)) {
    idx <- which(dat$healthCode == ids[i])
    for (j in seq(nfeat)) {
      dat[idx, featNames[j]] <- NormalTrans(dat[idx, featNames[j]])
    }   
  }
  
  dat
}


LoessDetrendedFeatures <- function(dat, featNames) {
  participantIds <- unique(dat$healthCode)
  nParticipants <- length(participantIds)
  featColumns <- match(featNames, colnames(dat))
  for (i in seq(nParticipants)) {
    #cat(i, "\n")
    idx1 <- which(dat$healthCode == participantIds[i])
    pdat <- dat[idx1,]
    o <- order(pdat$createdOn)
    pdat <- pdat[o,]
    xaxis <- seq(nrow(pdat))
    for (j in featColumns) {
      aux <- loess(pdat[, j] ~ xaxis)
      pdat[as.numeric(names(aux$residuals)), j] <- aux$residuals + median(pdat[, j], na.rm = TRUE)
    }    
    dat[idx1[o],] <- pdat
  }
  dat
}


GetPairnessScore <- function(dat, participants, sortResults = TRUE) {
  PairnessScore <- function(dat, participant) {
    sdat <- dat[dat$healthCode == participant,]
    xx <- as.data.frame.matrix(table(sdat$day, sdat$medTimepoint))
    xx[xx != 0] <- 1
    aux <- apply(xx, 1, sum)
    sum(aux == 2)/length(aux)  
  }
  n <- length(participants)
  pairedScore <- matrix(NA, n, 1, dimnames = list(participants, "score"))
  for (i in seq(n)) {
    pairedScore[i, 1] <- PairnessScore(dat, participant = participants[i])
  }
  if (sortResults) {
    pairedScore <- pairedScore[order(pairedScore[, 1], decreasing = TRUE),, drop = FALSE]
  }
  
  pairedScore
}


GetProportions <- function(uiTre, uiTod, thr = 0.05, r = 2) {
  out <- matrix(NA, 5, 3)
  rownames(out) <- c("treat", "tod", "treat_alone", "tod_alone", "treat_and_tod")
  colnames(out) <- c("lm", "nw", "arima")
  n <- nrow(uiTre)
  out[1, 1] <- round(sum(uiTre[, 1] <= thr)/n, r)
  out[2, 1] <- round(sum(uiTod[, 1] <= thr)/n, r)
  out[3, 1] <- round(sum(uiTre[, 1] <= thr & uiTod[, 1] > thr)/n, r)
  out[4, 1] <- round(sum(uiTod[, 1] <= thr & uiTre[, 1] > thr)/n, r)
  out[5, 1] <- round(sum(uiTre[, 1] <= thr & uiTod[, 1] <= thr)/n, r)
  ####
  out[1, 2] <- round(sum(uiTre[, 2] <= thr)/n, r)
  out[2, 2] <- round(sum(uiTod[, 2] <= thr)/n, r)
  out[3, 2] <- round(sum(uiTre[, 2] <= thr & uiTod[, 2] > thr)/n, r)
  out[4, 2] <- round(sum(uiTod[, 2] <= thr & uiTre[, 2] > thr)/n, r)
  out[5, 2] <- round(sum(uiTre[, 2] <= thr & uiTod[, 2] <= thr)/n, r)
  ####
  out[1, 3] <- round(sum(uiTre[, 3] <= thr)/n, r)
  out[2, 3] <- round(sum(uiTod[, 3] <= thr)/n, r)
  out[3, 3] <- round(sum(uiTre[, 3] <= thr & uiTod[, 3] > thr)/n, r)
  out[4, 3] <- round(sum(uiTod[, 3] <= thr & uiTre[, 3] > thr)/n, r)
  out[5, 3] <- round(sum(uiTre[, 3] <= thr & uiTod[, 3] <= thr)/n, r)
  
  out
}


PadjustMatrix <- function(x, mtMethod = "BH") {
  nr <- nrow(x)
  nc <- ncol(x)
  apvals <- matrix(p.adjust(as.vector(x), method = mtMethod), nr, nc)
  dimnames(apvals) <- dimnames(x)
  
  apvals
}


MyImage <- function(x, mar, featCex = 0.5, partCex = 0.5, line = 0.3, main = "") {
  par(mar = mar)
  nfeat <- nrow(x)
  nids <- ncol(x)
  image(t(x[rev(seq(nfeat)),]), col = redgreen(12), xaxt = "n", yaxt = "n")
  mtext("features", side = 2, at = 0.5, las = 0, cex = featCex, line = line)
  mtext("participants", side = 1, at = 0.5, las = 0, cex = partCex, line = line)
  mtext(main, side = 3, at = 0.5, las = 1, cex = 1, line = 0.5)
}

