#' ---
#' title: "The prevalence of trauma, PTSD, and other mental health issues: evidence from a representative Slovak sample"
#' author: "Ivan Ropovik"
#' date: "`r Sys.Date()`"
#' output:
#'    html_document:
#'       toc: true
#'       toc_float: true
#'       code_folding: show
#'       fig_retina: 2
#' always_allow_html: yes
#' ---
#+ setup, include=FALSE, message=FALSE
# libraries and settings --------------------------------------------------
#knitr::opts_chunk$set(echo=FALSE, warning = FALSE)

rm(list = ls())

# install required R libraries if not installed already
list.of.packages <- c("tidyverse", "magrittr", "devtools", "overlapping", "entropy", "metap", "MBESS", "CarletonStats", "goftest", "qqtest", "beepr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load required libraries
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)

# source
source("functions.R")

# read in the data
dat <- read_delim("data/data.csv", ";", escape_double = FALSE, col_types = cols(paper = col_character(), sample = col_integer(), variable = col_integer(), group = col_integer(), nClusters = col_integer(), n = col_integer(), mean = col_double(), sd = col_double(), SE = col_double(), decimalM = col_integer(), decimalSD = col_integer(), discreteItemsCount = col_integer(), nPerCellAssumed = col_integer(), notes = col_character()), na = "NA", trim_ws = TRUE)
# compute SD from SE
dat <- dat %>%
  mutate(sd = ifelse(is.na(sd) & !is.na(SE), SE*sqrt(n), sd),
         varType = 1,
         studyID = paste(paper, "/", sample, sep = ""),
         trialID = as.numeric(paper) + as.numeric(sample)/100) %>%
  group_by(trialID) %>% mutate(id = cur_group_id()) %>% ungroup()

data <- dat %>% select(id, variable, group, n, mean, sd, decimalM, decimalSD, varType, studyID) %>% as.data.frame()
data <- data[1:100,]

# calculate p-values using Monte Carlo method
set.seed(1)
pMC <- sim_distr(m = 10000, data = data, plot_flag = F)
# nizke p <- priliz rozdielne, vysoke p, priliz rovnake

# calculate p-values using regular ANOVA method
pAnova <- get_anova_p_vals(data)

# find the matrix value closer to .5
ps <- matrix(nrow = dim(pMC)[1], ncol = dim(pMC)[2])
for(i in 1:dim(pMC)[1]){
  for(j in 1:dim(pMC)[2]){
ps[i,j] <- ifelse(is.na(pMC[i,j]), NA, c(pMC[i,j], pAnova[i,j])[which.min(abs(c(pMC[i,j], pAnova[i,j]) - .5))])
  }
}
rownames(ps) <- rownames(pMC)
ps


# exclude significant p-values and and rescale the non-significant to (0, 1)
ps <- ifelse(ps < .05, NA, (ps-0.05)/0.95)

# Analysis based on p-values for individual within-study variables --------

# choose a single p-value for every independent study and permute
nIterations <- 1000
pSets <- matrix(nrow = nrow(ps), ncol = nIterations)
for(i in 1:nIterations){
  pSets[,i] <- apply(ps, 1, function(x){sample(x[!is.na(x)], size = 1)})
}
rownames(pSets) <- rownames(pMC)

# combine the permuted sets of p-values using Stouffer's method
pSetsStouffer <- apply(pSets, 2, function(x){
  1 - pnorm(sum(sapply(x, qnorm))/sqrt(length(x)))
})

# average over the iterations and compute the sd
(meanStouffer <- c("Mean Stouffer's p over iterations" = mean(pSetsStouffer),
                  "SD of Stouffer's p over iterations" = sd(pSetsStouffer)))

# Wald's z-test for Stouffer's method
(waldsStouffer <- c("zWaldStouffer" = zStouffer <- (meanStouffer[1] - .5)/(meanStouffer[2]/sqrt(nIterations)),
                    "pWaldStouffer" = 2*pnorm(abs(as.numeric(zStouffer)), lower.tail = F)))

# combine the permuted sets of p-values using Fisher's method
pSetsFisher <- apply(pSets, 2, function(x){
  pchisq(-2 * sum(log(x)), 2*length(x))
})

# average over the iterations and compute the sd
(meanFisher <- c("Mean Stouffer's p over iterations" = mean(pSetsFisher),
                  "SD of Stouffer's p over iterations" = sd(pSetsFisher)))

# Wald's z-test for Fisher's method
(waldsFisher <- c("zFisher" = zFisher <- (meanFisher[1] - .5)/(meanFisher[2]/sqrt(nIterations)),
                  "pFisher" = 2*pnorm(abs(as.numeric(zFisher)), lower.tail = F)))

# Combining the p-values into a study-level p-value -----------------------
assumedIntercorr <- .2
pStudyStouffer <- NA
pStudyFisher <- NA
for(i in 1:nrow(ps)){# urobit aj korekciu na zaklade effective number of tests
  corMatrix <- matrix(assumedIntercorr, nrow = sum(!is.na(ps[i,])), ncol = sum(!is.na(ps[i,])))
  diag(corMatrix) <- 1
  pStudyStouffer[i] <- stouffer(1 - as.numeric(na.omit(ps[i,])),
                                R = mvnconv(corMatrix, target = "p", cov2cor = T, side = 2),
                                adjust = "empirical", size = 100000)$p
  pStudyFisher[i] <- fisher(1 - as.numeric(na.omit(ps[i,])),
                                R = mvnconv(corMatrix, target = "p", cov2cor = T, side = 2),
                                adjust = "empirical", size = 100000)$p
}
pStudyStouffer
pStudyFisher

distrTests <- list(
  "K-S for Stouffer's p" = ks.test(pStudyStouffer,"punif", 0, 1),
  "K-S for Fisher's p" = ks.test(pStudyFisher,"punif", 0, 1),
  "A-D for Stouffer's p" = ad.test(pStudyStouffer, "punif", 0, 1),
  "A-D for Fisher's p" = ad.test(pStudyFisher, "punif", 0, 1)
)

qqtest::qqtest(pStudyStouffer, dist = "uniform", legend = T, xlim = c(0, 1), ylim = c(0, 1),
               xAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               yAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               bty = "n", drawPercentiles = T, nreps = 10000)


# Distribution overlap ----------------------------------------------------

overlapOut <- NA
for(i in 1:100){
  y <- runif(10000, 0, 1)
  overlapOut[i] <- overlap(list(pStudyStouffer,y), boundaries = list(from = 0, to = 1), plot = F)$OV
}
(overlapEst <- c("Overlap estimate mean" = mean(overlapOut),
                 "Overlap estimate SD" = sd(overlapOut)))

# Overlap density plot
plotData <- as.data.frame(cbind("pvalues" = c(pStudyStouffer, y), "Distribution" = c(rep(1, length(pStudyStouffer)), rep(2, length(y)))))
plotData$Distribution <- as.factor(plotData$Distribution)
ggplot(data = plotData, aes(x = pvalues, group = Distribution, fill = Distribution)) +
  geom_density(adjust = 0.5, alpha = .5) + theme_classic()

# lack of variability in SDs
nSim <- 10 # define the number of simulations
theoreticalSDSD <- theoreticalRangeSD <- empiricalSDSD <- empiricalRangeSD <- data.frame(NA)
listTheoreticalSDSD <- listTheoreticalRangeSD <- list(NA)
set.seed(1)
for(n in 1:nSim){
  for(i in 1:max(data$id)){
    for(j in 1:max(data$variable)){
      d <- data %>% filter(id == i, variable == j) # subset rows for each variable within each independent sample
      pooledSD <- mean(d$sd) # compute pooled SD (assuming the null)
      empiricalSDSD[i,j] <- sd(d$sd) # compute sd of sd
      empiricalRangeSD[i,j] <- diff(range(d$sd)) # compute range of sd
      theoreticalSDSD[i,j] <- d %>% rowwise() %>% transmute(sd(rnorm(n, mean, pooledSD))) %>% unlist() %>% sd()
      theoreticalRangeSD[i,j] <- d %>% rowwise() %>% transmute(sd(rnorm(n, mean, pooledSD))) %>% unlist() %>% range() %>% diff()
    }
  }
  listTheoreticalSDSD[[n]] <- theoreticalSDSD
  listTheoreticalRangeSD[[n]] <- theoreticalRangeSD %>% mutate_if(is.numeric, list(~na_if(., -Inf)))
}

# average over lists of sdsds and compute the probability of observing equal or smaller degree of variability in SDs
# higher probability = smaller than expected variability
excessSmallSDs <-  Reduce("+", lapply(listTheoreticalSDSD, function(x){x >= empiricalSDSD}))/nSim
excessSmallRanges <-  Reduce("+", lapply(listTheoreticalRangeSD, function(x){x >= empiricalRangeSD}))/nSim
rownames(excessSmallSDs) <- rownames(excessSmallRanges) <- rownames(pMC)
excessSmallSDs
excessSmallRanges

# find the smaller p-value for the SD and Ranges matrices (picking the less extreme p-value)
selectRangeSD <- matrix(nrow = dim(excessSmallSDs)[1], ncol = dim(excessSmallSDs)[2])
for(i in 1:dim(excessSmallSDs)[1]){
  for(j in 1:dim(excessSmallSDs)[2]){
    selectRangeSD[i,j] <- ifelse(is.na(excessSmallSDs[i,j]), NA, c(excessSmallSDs[i,j], excessSmallRanges[i,j])[which.max(abs(c(excessSmallSDs[i,j], excessSmallRanges[i,j]) - 1))])
  }
}
rownames(selectRangeSD) <- rownames(pMC)
selectRangeSD


######################################
# Graveyard
######################################


hist(pMC[pMC>.05], breaks = 100)
hist(pC[pC>.05], breaks = 100)

view(round(as.matrix(pMC, dimnames = NULL), 3))
view(round(as.matrix(pAnova, dimnames = NULL), 3))

allmetap(beckerp, method = "all")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EmpiricalBrownsMethod")
library(EmpiricalBrownsMethod)

# the call to combine p-values
empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=TRUE)
