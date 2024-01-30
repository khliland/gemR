# MS.GEM.02 Kristians  script design are group ms sex ageLP Jkn Perm.r

#setwd("H:/A_R/R/R data/R GEM Method publ/MS")
#setwd("D:/R GEM Method publ/MS")

setwd("C:/Users/ellen.fargestad/Documents/R copy.C/R MS rep")
getwd()

library(gemR)
library(permute)
library(gprofiler2)
library(corrplot)

rm(list = ls())

#####################################################
# Load data 
#####################################################

#load("./R copy.C/R MS rep/WS/GEM loaded and saved after PLS-GEM.RData")
# save(My.Array.1,file="./Data/My.Array.1.RData")

load(file="./Data/My.Array.1.RData")

##############################################
# Function
##############################################

run.ci <- function(CI){
  CI = as.data.frame(CI)
  id1 <- which(CI$Left > 0);n1 <-rownames(CI)[id1]
  id2 <- which(CI$Right< 0);n2 <- rownames(CI)[id2]
  n.select <- c(n1,n2)
  out <- list(n.select, n1, n2)
  return(out)
}

check.error <- function(pls.mod,ncomp,y){
  # Argument: pls model
  # output:   gives the proportion of errornes classification
  error <- vector()
  
  n <- length(y)
  for(i in 1:ncomp){
    y.predicted <- pls.mod$validation$pred[,1,i]
    y.pred.sign <- sign(y.predicted)
    error[i]    <- sum(sign(y)!=sign(y.predicted))/n
  }
  error <- c(0.5,error)
  return(error) 
}

propVar <- function(object){
  vars <- object$sdev^2
  vars/sum(vars)
}

##############################################
# Load
##############################################

#load(file='H:/A_R/R/R data/R MS publ.SR/Data/My.Array.1.Rdata')
#names(My.Array.1)
Design    <- My.Array.1$Design;head(Design)
clinics   <- My.Array.1$clinics;clinics[1:2,]
D         <- Design
M         <- as.matrix(as.data.frame(unlist(My.Array.1$M.1)))
colnames(M)<- substr(colnames(M),3,20) # omitt x that is added the protein names;colnames(M)[1:10]
proteins  <- M;rownames(proteins)<- rownames(Design)
colnames(proteins)[1:10]

group     <- as.numeric(My.Array.1$factorA)
ms        <- as.numeric(My.Array.1$factorB)
sex       <- rep(NA,each=dim(D)[1]);names(sex)<- rownames(D)
sex       <- as.numeric(sex)
id        <- which(clinics$Sex=='M');sex[id]<- -1
id        <- which(clinics$Sex=='F');sex[id]<-  1
Sex.FM    <- clinics$Sex
ageLP     <- clinics$Age_at_LP
group.12  <- (group+3)/2
ms.12     <- (ms+3)/2
sex.12    <- (sex+3)/2

###################
## colors and pch  
###################

# in plots of factor groups and ms:  
my.cls.s.group  <- c('blue','red')[group.12];# ctrl- MS
my.pch.s.ms     <- c(2,15)[ms.12]

# in plots of factor sex:  
my.cls.s.sex    <- c('skyblue','deeppink')[sex.12]

# in plots of ageLP: Make groups of age for the colorcoding in the plots
# may be simplified....
ageLP.gr <- quantile(ageLP, probs = seq(0, 1, 0.25))
id1        <- which(ageLP<ageLP.gr[2])
id2        <- which(ageLP<ageLP.gr[3])
id3        <- which(ageLP<ageLP.gr[4])
id.q1 <- id1;                         sort(ageLP[id.q1])
id.q2 <- setdiff(id2,id1);            sort(ageLP[id.q2])
id.q3 <- setdiff(id3,c(id1,id2));     sort(ageLP[id.q3])
id.q4 <- setdiff(c(1:dim(M)[1]),c(id1,id2,id3)); sort(ageLP[id.q4])
length(c(id.q1,id.q2,id.q3,id.q4))

qu.gr <- rep(NA,each=dim(M)[1])
qu.gr[id.q1]<- 1
qu.gr[id.q2]<- 2
qu.gr[id.q3]<- 3
qu.gr[id.q4]<- 4
cls.4 <- c('brown','green4','skyblue','black')
my.cls.s.ageLP <- cls.4 [qu.gr]
my.pch.s.ageLP <- c()


# Organise data as data.frame.
###################################

MS.data <- data.frame(proteins = I(proteins),
                      group = group,
                      ms = ms,
                      sex = sex,
                      ageLP = ageLP)

# Permuted data 
#################
id        <- shuffle(group);  group.p <- group[id]
id        <- shuffle(ms);     ms.p    <- ms[id]
id        <- shuffle(sex);    sex.p   <- sex[id]
id        <- shuffle(ageLP);  ageLP.p <- ageLP[id]

MS.data.p    <- data.frame(
  proteins =I(proteins),
  group.p  = group.p, 
  ms.p     = ms.p,
  sex.p    = sex.p,
  ageLP.p  = ageLP.p)

# MS.data <- data.frame(proteins = I(proteins),
#                       group = factor(group),
#                       ms = factor(ms),
#                       sex = factor(sex),
#                       ageLP = ageLP)


# MS.data.p    <- data.frame(
#     proteins =I(proteins),
#     group.p  = factor(group.p), 
#     ms.p     = factor(ms.p),
#     sex.p    = factor(sex.p),
#     ageLP.p  = factor(ageLP.p))

# save to reproduce the results 
save(MS.data.p,file='./Data/MS.data.p.RData')
load(file='./Data/MS.data.p.RData')

##########################
## GEM
##########################

# Step 1 - GLM:
MS.gem <- GEM(proteins ~ ms + group + sex + ageLP, data = MS.data)

# Matrices of effects and ER values can be extracted for custom analyses
E.group       <- MS.gem$effects$group
E.ms          <- MS.gem$effects$ms
E.sex         <- MS.gem$effects$sex
E.ageLP       <- MS.gem$effects$ageLP
ER.group      <- MS.gem$ER.values$group
ER.ms         <- MS.gem$ER.values$ms
ER.sex        <- MS.gem$ER.values$sex
ER.ageLP      <- MS.gem$ER.values$ageLP

# Step 2 - PLS
ncomp <- 10
# PLS analysis 

# 'group' effect
################
y                 <- group
MS.pls.group      <- pls(MS.gem, 'group', ncomp, validation = "LOO",
                 jackknife = TRUE, df.used = 3)
#colSums(MS.pls.group$classes == as.numeric(MS.data$group))
error.group     <- check.error(MS.pls.group,ncomp=ncomp,y);error.group
ncomp.opt         <- 2
MS.coef.group     <- coef(MS.pls.group)
MS.scores.group   <- scores(MS.pls.group)
MS.loadings.group <- loadings(MS.pls.group)
MS.jack.group     <- MS.pls.group$jack[,1,ncomp.opt]
MS.signif.group   <- names(which(MS.jack.group < 0.05))
cls.v             <- rep(1,each=dim(proteins)[2]);names(cls.v)<- colnames(proteins)
cls.v[MS.signif.group]<- 2
cls.v             <- as.numeric(cls.v)
my.cls.v.group    <- c('gray70','black')[cls.v]
my.pch.v.group    <- c(1,19)[cls.v]

# 'ms' effect
###############
y         <- ms
MS.pls.ms <- pls(MS.gem, 'ms', ncomp, validation = "LOO",
                 jackknife = TRUE, df.used = 3)
#colSums(MS.pls.ms$classes == as.numeric(MS.data$ms))
error.ms     <- check.error(pls.mod=MS.pls.ms,ncomp=ncomp,y=y);error.ms
ncomp.opt      <- 2
MS.coef.ms     <- coef(MS.pls.ms)
MS.scores.ms   <- scores(MS.pls.ms)
MS.loadings.ms <- loadings(MS.pls.ms)
MS.jack.ms     <- MS.pls.ms$jack[,1,ncomp.opt]
MS.signif.ms   <- which(MS.jack.ms < 0.05)
cls.v[MS.signif.ms]<- 2
cls.v             <- as.numeric(cls.v)
my.cls.v.ms    <- c('gray70','black')[cls.v]
my.pch.v.ms    <- c(1,19)[cls.v]

# 'sex' effect
################
y          <- sex
MS.pls.sex <- pls(MS.gem, 'sex', ncomp, validation =  "LOO",
                 jackknife = TRUE, df.used = 3)
#colSums(MS.pls.sex$classes == as.numeric(MS.data$sex))
error.sex     <- check.error(pls.mod=MS.pls.sex,ncomp=ncomp,y=y);error.sex
ncomp.opt      <- 2
MS.coef.sex     <- coef(MS.pls.sex)
MS.scores.sex   <- scores(MS.pls.sex)
MS.loadings.sex <- loadings(MS.pls.sex)
MS.jack.sex     <- MS.pls.sex$jack[,1,ncomp.opt]
MS.signif.sex   <- which(MS.jack.sex < 0.05)
cls.v             <- rep(1,each=dim(proteins)[2]);names(cls.v)<- colnames(proteins)
cls.v[MS.signif.sex]<- 2
cls.v             <- as.numeric(cls.v)
my.cls.v.sex    <- c('gray70','black')[cls.v]
my.pch.v.sex    <- c(1,19)[cls.v]

# 'ageLP' effect
################
MS.pls.ageLP <- pls(MS.gem, 'ageLP', ncomp, validation = "LOO",
                  jackknife = TRUE, df.used = 3)
#colSums(MS.pls.sex$classes == as.numeric(MS.data$ageLP))
ncomp.opt         <- 2
MS.coef.ageLP     <- coef(MS.pls.ageLP)
MS.scores.ageLP   <- scores(MS.pls.ageLP)
MS.loadings.ageLP <- loadings(MS.pls.ageLP)
MS.jack.ageLP     <- MS.pls.ageLP$jack[,1,ncomp.opt]
MS.signif.ageLP   <- which(MS.jack.ageLP < 0.05)
cls.v             <- rep(1,each=dim(proteins)[2]);names(cls.v)<- colnames(proteins)
cls.v[MS.signif.ageLP]<- 2
cls.v             <- as.numeric(cls.v)
my.cls.v.ageLP    <- c('gray70','black')[cls.v]
my.pch.v.ageLP    <- c(1,19)[cls.v]
ageLP.predicted   <- MS.pls.ageLP$fitted.values[,,2]

#####################################################################
# PLS permuted data
#####################################################################

# Step 1 - GLM:
MS.gem.p <- GEM(proteins ~ group.p + ms.p + sex.p + ageLP.p, data = MS.data.p)

# goup
################
my.ncomp    <- 10
y           <- group.p
pls.mod  <- pls(MS.gem.p, 'group.p', ncomp= my.ncomp, validation='LOO', df.used = 3)
error.group.p   <- check.error(pls.mod,ncomp=my.ncomp,y);error.group.p
#colSums(MS.pls.group.p$classes == as.numeric(MS.data.p$group.p))

# ms
#############
my.ncomp    <- 10
y           <- ms
pls.mod     <- pls(MS.gem.p, 'ms.p', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
error.ms.p     <- check.error(pls.mod,ncomp=my.ncomp,y);error.ms.p

# sex
#############
my.ncomp    <- 10
y           <- sex
pls.mod     <- pls(MS.gem.p, 'sex.p', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
error.sex.p     <- check.error(pls.mod,ncomp=my.ncomp,y);error.sex.p

# ageLP
############
my.ncomp    <- 10
y           <- ageLP.p
pls.mod     <- pls(MS.gem.p, 'ageLP.p', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
ageLP.p.predicted <- pls.mod$fitted.values[,,2]


##################
# PLS plots
##################

my.ylim.error <- c(0,0.6)

oldpar <- par(mfrow = c(4,3), mar = c(4,4,2,1))


# group
dataset <- 'PLS-GEM group'
scoreplot(MS.pls.group, main ="", panel.first = abline(h=0, v=0, col="gray"),
          col = my.cls.s.group, pch = my.pch.s.ms)
title(paste0(letters[1],") ",dataset, " - Scores" ), adj = 0,line = 1)
loadingplot(MS.pls.group, scatter = TRUE,col=my.cls.v.group, pch=my.pch.v.group,
            main = "",panel.first = abline(h=0, v=0, col="gray"))
title(paste0(letters[2],") ",dataset, " - Kristians Loadings" ), adj = 0,line = 1)

# det beste er om "Classification accuracy" kan starte p?? 0. Kristian: kan du endre til det? 
# Men dette funker ikke n??r jeg ikke har definert factorene som factorer
# corrplot(MS.pls.ms, main = "Correlation loadings", label = "names")# f??r feilmelding
# plot(colSums(MS.pls.group$classes == as.numeric(MS.data$group)),
#      ylim = c(0,100),
#      ylab = "# correct", xlab = "# components",panel.first = grid(),
#      main = "Classification accuracy")
plot(error.group,ylim=my.ylim.error)
lines(error.group,ylim=my.ylim.error)
lines(error.group.p,ylim=my.ylim.error,col='blue')
abline(h=0.5,lty=2,col='gray')
title(paste0(letters[3],") ",dataset, " - Error" ), adj = 0,line = 1)

# ms
dataset <- 'PLS-GEM ms'
scoreplot(MS.pls.ms, main = "", panel.first = abline(h=0, v=0, col="gray"),
          col = my.cls.s.group, pch = my.pch.s.ms)
title(paste0(letters[4],") ",dataset, " - Scores" ), adj = 0,line = 1)
loadingplot(MS.pls.ms, scatter = TRUE,col=my.cls.v.ms, pch=my.pch.v.ms,
            main = "",panel.first = abline(h=0, v=0, col="gray"))
title(paste0(letters[5],") ",dataset, " - Loadings" ), adj = 0,line = 1)
plot(error.ms,ylim=my.ylim.error)
lines(error.ms,ylim=my.ylim.error)
lines(error.ms.p,ylim=my.ylim.error,col='blue')
abline(h=0.5,lty=2,col='gray')
title(paste0(letters[6],") ",dataset, " - Error" ), adj = 0,line = 1)

# sex
dataset <- 'PLS-GEM sex'
scoreplot(MS.pls.sex, main = "", panel.first = abline(h=0, v=0, col="gray"),
          col = my.cls.s.sex, pch = my.pch.s.ms)
title(paste0(letters[4],") ",dataset, " - Scores" ), adj = 0,line = 1)
loadingplot(MS.pls.sex, scatter = TRUE,col=my.cls.v.sex, pch=my.pch.v.sex,
            main = "",panel.first = abline(h=0, v=0, col="gray"))
title(paste0(letters[5],") ",dataset, " - Loadings" ), adj = 0,line = 1)

plot(error.sex,ylim=my.ylim.error)
lines(error.sex,ylim=my.ylim.error)
lines(error.sex.p,ylim=my.ylim.error,col='blue')
abline(h=0.5,lty=2,col='gray')
title(paste0(letters[6],") ",dataset, " - Error" ), adj = 0,line = 1)

# ageLP
dataset <- 'PLS-GEM ageLP'
ncomp.opt.ageLP <- 3
scoreplot(MS.pls.ageLP, main = "", panel.first = abline(h=0, v=0, col="gray"),
          col = my.cls.s.ageLP, pch = my.pch.s.ms)
title(paste0(letters[7],") ",dataset, " - Scores" ), adj = 0,line = 1)
loadingplot(MS.pls.ageLP, scatter = TRUE,col=my.cls.v.ageLP, pch=my.pch.v.ageLP,
            main = "",panel.first = abline(h=0, v=0, col="gray"))
title(paste0(letters[8],") ",dataset, " - Loadings" ), adj = 0,line = 1)

ageLP.predicted <- MS.pls.ageLP$fitted.values[,,ncomp.opt.ageLP]
plot(ageLP,ageLP.predicted,col=my.cls.s.ageLP,pch=my.pch.s.ms)
title(paste0(letters[9],") ",dataset, " - Pred vs observed" ), adj = 0,line = 1)

