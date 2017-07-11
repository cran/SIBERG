### R code from vignette source 'SIBER.Rnw'

###################################################
### code chunk number 1: SIBER.Rnw:42-44
###################################################
options(width=80)
options(continue=' ')


###################################################
### code chunk number 2: loadLibrary
###################################################
library(SIBERG)


###################################################
### code chunk number 3: SIBER.Rnw:85-90
###################################################
set.seed(1000)
N <- 100 # sample size
G <- 200 # number of simulated genes
# RNAseq count data simulated from NB model with mean 1000, dispersion=0.2
Dat <-  matrix(rnbinom(G*N, mu=1000, size=1/0.2), nrow=G) 


###################################################
### code chunk number 4: SIBER.Rnw:96-97
###################################################
SIBER(y=Dat[1, ], model='LN')


###################################################
### code chunk number 5: SIBER.Rnw:102-103
###################################################
SIBER(y=Dat[1, ], model='NB')


###################################################
### code chunk number 6: SIBER.Rnw:108-109
###################################################
SIBER(y=Dat[1, ], model='GP')


###################################################
### code chunk number 7: SIBER.Rnw:115-116
###################################################
SIBER(y=log(Dat[1, ]+1), model='NL')


###################################################
### code chunk number 8: SIBER.Rnw:135-137
###################################################
library(edgeR)
TMM <- calcNormFactors(Dat, method='TMM')


###################################################
### code chunk number 9: SIBER.Rnw:154-155
###################################################
SIBER(y=Dat[1, ], d=1/TMM, model='LN')


###################################################
### code chunk number 10: SIBER.Rnw:227-236
###################################################
data(simDat)
ind <- 1
# true parameter generating the simulated data
parList$LN[ind, ]
# fit by E model
fitLN(y=dataList$LN[ind, ], base=exp(1), eps=1, model='E')
# fit by V model. 
fitLN(y=dataList$LN[ind, ], base=exp(1), eps=1, model='V')



###################################################
### code chunk number 11: SIBER.Rnw:242-252
###################################################
ind <- 5 # 0-inflated gene
# true parameter generating the simulated data
parList$LN[ind, ]
# fit by E model. 0-inflated model is disabled by setting zeroPercentThr=1.
# the result is biased. 
fitLN(y=dataList$LN[ind, ], base=exp(1), eps=1, model='E', zeroPercentThr=1)
# fit by 0-inflated model. 0-inflated model overrides the E model since percentage
# of observed zero counts exceeds the threshold.
fitLN(y=dataList$LN[ind, ], base=exp(1), eps=1, model='E', zeroPercentThr=0.2)



###################################################
### code chunk number 12: sessionInfo
###################################################
getwd()
sessionInfo()


