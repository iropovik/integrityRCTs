#' ---
#' title: "The integrity of data underlying randomized controlled studies of psychological and educational interventions"
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

readRDS("workspace.RDS")
rm(list = ls())

# install required R libraries if not installed already
if (!requireNamespace("drat", quietly = TRUE)) install.packages("drat")
drat::addRepo("daqana")

list.of.packages <- c("tidyverse", "magrittr", "devtools", "overlapping", "entropy", "metap", "poolr", "MBESS", "CarletonStats", "goftest", "qqtest", "dqsample", "beepr")
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
  mutate(sd = ifelse(is.na(sd) & !is.na(SE), round(SE*sqrt(n), decimalSD), sd),
         varType = 1,
         studyID = paste(paper, "/", sample, sep = ""),
         trialID = as.numeric(paper) + as.numeric(sample)/100,
         rowNo = 1:nrow(dat)) %>%
  group_by(trialID) %>% mutate(id = cur_group_id()) %>% ungroup()

data <- dat %>% select(id, variable, group, n, mean, sd, decimalM, decimalSD, varType, studyID) %>%  as.data.frame()

# calculate p-values using Monte Carlo method
set.seed(1)
# ps <- sim_distr(m = 10000, data = data, plot_flag = F) #6.65115min
ps <- readRDS("ps.RDS") # nizke p <- priliz rozdielne, vysoke p, priliz rovnake

# Compensate for the lack of precision in MC calculation in p-values
ps <- ifelse(as.matrix(ps) == 0, .000005, ifelse(as.matrix(ps) == 1, .999995, as.matrix(ps)))

# how many significant
paste(round(table(as.matrix(ps) < .05)[2]/sum(table(as.matrix(ps) < .05))*100, 3), "%")

# Analysis based on p-values for individual within-sample variables --------

# choose a single p-value for every independent sample and permute
set.seed(1)
nIterations <- 100000

pSets <- matrix(nrow = nrow(ps), ncol = nIterations)
for(i in 1:nIterations){  # excluding significant p-values and and rescaling the non-significant p-values to (0, 1)
  pSets[,i] <- apply(ifelse(as.matrix(ps) < .05, NA, (as.matrix(ps)-0.05)/0.95), 1,
                     function(x){sample(x[!is.na(x)], size = 1)})
}
rownames(pSets) <- rownames(ps)

# Compensate for the lack of precision in MC calculation in p-values
pSets <- ifelse(as.matrix(pSets) == 0, .00001, ifelse(as.matrix(pSets) == 1, .99999, as.matrix(pSets)))

# saveRDS(pSets, "pSets.RDS")
# pSets <- readRDS("pSets.RDS")

# combine the permuted sets of p-values using Stouffer's method
# low combined p-value means too small baseline difference, inconsistent with random sampling
pSetsStouffer <- apply(pSets, 2, function(x){
  1 - pnorm(sum(sapply(x, qnorm), na.rm = T)/sqrt(length(x[!is.na(x)])))
})

# average over the iterations and compute the sd
(meanStouffer <- c("Mean Stouffer's p over iterations" = mean(pSetsStouffer),
                  "SD of Stouffer's p over iterations" = sd(pSetsStouffer)))

# Wald's z-test for Stouffer's method
(waldsStouffer <- c("zWaldStouffer" = zStouffer <- (meanStouffer[1] - .5)/(meanStouffer[2]/sqrt(nIterations)),
                    "pWaldStouffer" = 2*pnorm(abs(as.numeric(zStouffer)), lower.tail = F)))

ggplot(mapping = aes(pSetsStouffer, fill = pSetsStouffer)) + geom_density(adjust = 0.7, fill = "#4b88ad", alpha = .5) + theme_classic()

# combine the permuted sets of p-values using Fisher's method
pSetsFisher <- apply(pSets, 2, function(x){
  pchisq(-2 * sum(log(x), na.rm = T), 2*length(x[!is.na(x)]))
})

# average over the iterations and compute the sd
(meanFisher <- c("Mean Fisher's p over iterations" = mean(pSetsFisher),
                 "SD of Fisher's p over iterations" = sd(pSetsFisher)))

# Wald's z-test for Fisher's method
(waldsFisher <- c("zFisher" = zFisher <- (meanFisher[1] - .5)/(meanFisher[2]/sqrt(nIterations)),
                  "pFisher" = 2*pnorm(abs(as.numeric(zFisher)), lower.tail = F)))

ggplot(mapping = aes(pSetsFisher, fill = pSetsFisher)) + geom_density(adjust = 0.7, fill = "#4b88ad", alpha = .5) + theme_classic()

# show that the distribution of stouffer's p is uniform under uniform distribution of p-values
# yy <- NA
# for(i in 1:100){
#   xx <- runif(10000, 0, 1)
#   yy[i] <- 1 - pnorm(sum(sapply(xx, qnorm), na.rm = T)/sqrt(length(xx[!is.na(xx)])))
# }
# hist(yy)

# Combining the p-values into a study-level p-value -----------------------
psCombined <- list()
for(i in gsub("/.*", "", rownames(ps))){
  psCombined[[i]] <- as.data.frame(ps) %>% rownames_to_column() %>% filter(grepl(i, rowname, fixed = TRUE)) %>% select(-rowname) %>% unlist()
}
psCombined <- lapply(psCombined, function(x){as.numeric(x[!is.na(x)])})

assumedIntercorr <- .3 # urobit aj korekciu na zaklade effective number of tests
pStudyStouffer <- list(NA)
pStudyFisher <- list(NA)
for(i in 1:length(psCombined)){
  corMatrix <- matrix(assumedIntercorr, nrow = length(psCombined[[i]]), ncol = length(psCombined[[i]]))
  diag(corMatrix) <- 1
  pStudyStouffer[[i]] <- tryCatch(stouffer(1 - psCombined[[i]],
                                           R = mvnconv(corMatrix, target = "p", cov2cor = T, side = 2),
                                           adjust = "empirical", size = 10000)$p, error = function(e) NULL)
  pStudyFisher[[i]] <- tryCatch(fisher(1 - psCombined[[i]],
                                       R = mvnconv(corMatrix, target = "p", cov2cor = T, side = 2),
                                       adjust = "empirical", size = 10000)$p, error = function(e) NULL)
}
names(pStudyStouffer) <- names(pStudyFisher) <- rownames(psCombined)
pStudyStouffer <- unlist(pStudyStouffer)
pStudyFisher <- unlist(pStudyFisher)

list("Probabilities by Stouffer's method" =
       list("<.001" = prop.table(table(pStudyStouffer<.001))[2]*100,
          "<.01" = prop.table(table(pStudyStouffer<.01))[2]*100,
          "<.05" = prop.table(table(pStudyStouffer<.05))[2]*100),
     "Probabilities by Fisher's method" =
       list("<.001" = prop.table(table(pStudyFisher<.001))[2]*100,
          "<.01" = prop.table(table(pStudyFisher<.01))[2]*100,
          "<.05" = prop.table(table(pStudyFisher<.05))[2]*100))

(distrTests <- list(
  "K-S for Stouffer's p" = ks.test(pStudyStouffer,"punif", 0, 1),
  "K-S for Fisher's p" = ks.test(pStudyFisher,"punif", 0, 1),
  "A-D for Stouffer's p" = ad.test(pStudyStouffer, "punif", 0, 1),
  "A-D for Fisher's p" = ad.test(pStudyFisher, "punif", 0, 1)
))

qqtest(pStudyStouffer, dist = "uniform", legend = T, xlim = c(0, 1), ylim = c(0, 1),
               xAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               yAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               bty = "n", nreps = 100000, cex = .5)

qqtest(pStudyFisher, dist = "uniform", legend = T, xlim = c(0, 1), ylim = c(0, 1),
               xAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               yAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               bty = "n", nreps = 100000, cex = .5)

# Distribution overlap ----------------------------------------------------
set.seed(1)
y <- runif(1000000, 0, 1)
overlap(list(pStudyStouffer,y), boundaries = list(from = 0, to = 1), plot = F)$OV

# Overlap density plot
plotData <- as.data.frame(cbind("pvalues" = c(pStudyStouffer, y), "Distribution" = c(rep(1, length(pStudyStouffer)), rep(2, length(y)))))
plotData$Distribution <- as.factor(plotData$Distribution)
ggplot(data = plotData, aes(x = pvalues, group = Distribution, fill = Distribution)) +
  geom_density(alpha = .5) + theme_classic()

# Lack of variability in SDs ----------------------------------------------
# To do: speed up
# rnormSdFun <- function(df){
#   sd(rnorm(n = df[["n"]],mean =  df[["mean"]], sd = pooledSD))
# }

# set.seed(1)
# d %>% rowwise() %>% transmute(sd(rnorm(n, mean, pooledSD))) %>% unlist() %>% sd()
# set.seed(1)
# d %>% rowwise() %>% rnormSdFun()

nSim <- 1000 # define the number of simulations
theoreticalSDSD <- theoreticalRangeSD <- empiricalSDSD <- empiricalRangeSD <- data.frame(NA)
listTheoreticalSDSD <- listTheoreticalRangeSD <- list(NA)
set.seed(1)
system.time(for(n in 1:nSim){
              for(i in 1:max(data$id)){
                for(j in 1:max(data$variable)){
                  d <- data %>% filter(id == i, variable == j) # subset rows for each variable within each independent sample
                  pooledSD <- mean(d$sd) # compute pooled SD (assuming the null)
                  empiricalSDSD[i,j] <- sd(d$sd) # compute sd of sd
                  theoreticalSDSD[i,j] <- d %>% rowwise() %>% transmute(sd(rnorm(n, mean, pooledSD))) %>% unlist() %>% sd()
                }
              }
              listTheoreticalSDSD[[n]] <- theoreticalSDSD
            }
)

set.seed(1)
for(n in 1:nSim){
  for(i in 1:max(data$id)){
    for(j in 1:max(data$variable)){
      d <- data %>% filter(id == i, variable == j) # subset rows for each variable within each independent sample
      pooledSD <- mean(d$sd) # compute pooled SD (assuming the null)
      empiricalRangeSD[i,j] <- diff(range(d$sd)) # compute range of sd
      theoreticalRangeSD[i,j] <- d %>% rowwise() %>% transmute(sd(rnorm(n, mean, pooledSD))) %>% unlist() %>% range() %>% diff()
    }
  }
  listTheoreticalRangeSD[[n]] <- theoreticalRangeSD %>% mutate_if(is.numeric, list(~na_if(., -Inf))) # remove -Inf values
}

# average over lists of sdsds and compute the probability of observing equal or smaller degree of variability in SDs
# higher probability = smaller than expected variability
excessSmallSDs <-  Reduce("+", lapply(listTheoreticalSDSD, function(x){x >= empiricalSDSD}))/nSim
excessSmallRanges <-  Reduce("+", lapply(listTheoreticalRangeSD, function(x){x >= empiricalRangeSD}))/nSim
rownames(excessSmallSDs) <- rownames(excessSmallRanges) <- rownames(ps)
excessSmallSDs
excessSmallRanges

# find the smaller p-value for the SD and Ranges matrices (picking the less extreme p-value)
# insuffVariability <- matrix(nrow = dim(excessSmallSDs)[1], ncol = dim(excessSmallSDs)[2])
# for(i in 1:dim(excessSmallSDs)[1]){
#   for(j in 1:dim(excessSmallSDs)[2]){
#     insuffVariability[i,j] <- ifelse(is.na(excessSmallSDs[i,j]), NA, c(excessSmallSDs[i,j], excessSmallRanges[i,j])[which.max(abs(c(excessSmallSDs[i,j], excessSmallRanges[i,j]) - 1))])
#   }
# }
# rownames(insuffVariability) <- rownames(ps)
# insuffVariability

# Average into trial-level probabilities
# studyInsuffVariability <- apply(insuffVariability, 1, mean, na.rm = T)
# names(studyInsuffVariability) <- rownames(ps)
#
# studyInsuffVariability
# pStudyStouffer

# Compensate for the lack of precision in MC calculation in p-values

insuffVariability <- excessSmallSDs #delete if analysis based on ranges is carried out
insuffVariability <- as.data.frame(ifelse(as.matrix(insuffVariability) == 0, .0005, ifelse(as.matrix(insuffVariability) == 1, .9995, as.matrix(insuffVariability))))

insuffVariabilityTrial <- list()
for(i in gsub("/.*", "", rownames(ps))){
  insuffVariabilityTrial[[i]] <- insuffVariability %>% rownames_to_column() %>% filter(grepl(i, rowname, fixed = TRUE)) %>% select(-rowname) %>% unlist()
}
insuffVariabilityTrial <- lapply(insuffVariabilityTrial, function(x){as.numeric(x[!is.na(x)])})

sdStudyStouffer <- list(NA)
sdStudyFisher <- list(NA)
for(i in 1:length(insuffVariabilityTrial)){
  corMatrix <- matrix(assumedIntercorr, nrow = length(insuffVariabilityTrial[[i]]), ncol = length(insuffVariabilityTrial[[i]]))
  diag(corMatrix) <- 1
  sdStudyStouffer[[i]] <- tryCatch(stouffer(1 - insuffVariabilityTrial[[i]],
                                           R = mvnconv(corMatrix, target = "p", cov2cor = T, side = 2),
                                           adjust = "empirical", size = 100000)$p, error = function(e) NULL)
  sdStudyFisher[[i]] <- tryCatch(fisher(1 - insuffVariabilityTrial[[i]],
                                       R = mvnconv(corMatrix, target = "p", cov2cor = T, side = 2),
                                       adjust = "empirical", size = 100000)$p, error = function(e) NULL)
}
names(sdStudyStouffer) <- names(sdStudyFisher) <- rownames(insuffVariabilityTrial)
sdStudyStouffer <- unlist(sdStudyStouffer)
sdStudyFisher <- unlist(sdStudyFisher)

# read RDS objects
sdStudyStouffer <- readRDS("sdStudyStouffer.RDS")
sdStudyFisher <- readRDS("sdStudyFisher.RDS")

mSdCor <- .3
pSdStouffer <- apply(cbind(pStudyFisher,sdStudyFisher), 1, function(x){
  stouffer(x, R = mvnconv(matrix(c(1, mSdCor, mSdCor, 1), nrow = 2, ncol = 2),
                          target = "p", cov2cor = T, side = 2), adjust = "empirical", size = 10000)$p})
pSdFisher <- apply(cbind(pStudyFisher,sdStudyFisher), 1, function(x){
  fisher(x, R = mvnconv(matrix(c(1, mSdCor, mSdCor, 1), nrow = 2, ncol = 2),
                        target = "p", cov2cor = T, side = 2), adjust = "empirical", size = 10000)$p})

list("Probabilities by Stouffer's method" =
       list("<.001" = prop.table(table(pSdStouffer<.001))[2]*100,
            "<.01" = prop.table(table(pSdStouffer<.01))[2]*100,
            "<.05" = prop.table(table(pSdStouffer<.05))[2]*100),
     "Probabilities by Fisher's method" =
       list("<.001" = prop.table(table(pSdFisher<.001))[2]*100,
            "<.01" = prop.table(table(pSdFisher<.01))[2]*100,
            "<.05" = prop.table(table(pSdFisher<.05))[2]*100))

qqtest(pSdStouffer, dist = "uniform", legend = T, xlim = c(0, 1), ylim = c(0, 1),
               xAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               yAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               bty = "n", nreps = 100000, cex = .5)

qqtest(pSdFisher, dist = "uniform", legend = T, xlim = c(0, 1), ylim = c(0, 1),
               xAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               yAxisProbs = c(0.2, 0.4, 0.6, 0.8, 1),
               bty = "n", nreps = 100000, cex = .5)

# GRIM & GRIMMER inconsistencies ------------------------------------------

grimAndGrimmer(dat)
# proportion of GRIM (M) inconsistencies
paste(round((table(dat$grim)[1]/table(!is.na(dat$discreteItemsCount))[2])*100, 2), "%")
# proportion of GRIMMER (SD) inconsistencies
paste(round(((table(dat$grimmer)[1] - table(dat$grim)[1])/table(!is.na(dat$discreteItemsCount))[2])*100, 2), "%")
# proportion of GRIM or GRIMMER (M | SD) inconsistencies
paste(round((table(dat$grimmer)[1]/table(!is.na(dat$discreteItemsCount))[2])*100, 2), "%")

######################################
# Graveyard
######################################
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

allmetap(beckerp, method = "all")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EmpiricalBrownsMethod")
library(EmpiricalBrownsMethod)

# the call to combine p-values
empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=TRUE)
