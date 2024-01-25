# MS.GEM.01 Run A B C quant1 Jkn Perm.r

#setwd("H:/A_R/R/R data/R GEM Method publ/MS")
#setwd("D:/R GEM Method publ/MS")

setwd("C:/Users/ellen.fargestad/Documents/R copy.C/R MS rep")
getwd()

rm(list = ls())
library(gemR)
library(permute)
library(gprofiler2)

l           <- 0.05

#####################################################
# Load data
#####################################################

#load("./R copy.C/R MS rep/WS/GEM loaded and saved after PLS-GEM.RData")
# save(My.Array.1,file="./Data/My.Array.1.RData")

load(file="My.Array.1.RData")

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

##################################################################
# cols and points
#################################################################

cls.2 <- c('blue','red')
points.2 <- c(2,15)

##############################################
# Load
##############################################

#load(file='H:/A_R/R/R data/R MS publ.SR/Data/My.Array.1.Rdata')
#names(My.Array.1)
Design.1  <- My.Array.1$Design;head(Design.1)
clinics.1 <- My.Array.1$clinics;clinics.1[1:2,]
D.1       <- Design.1
M.1       <- My.Array.1$M.1    # 779 all proteins
M.1.c     <- My.Array.1$M.1.c  # 357 common proteins of the two data sets
M         <- M.1;dim(M)
M         <- as.matrix(as.data.frame(unlist(M)))
colnames(M)<- substr(colnames(M),3,20) # omitt x that is added the protein names
colnames(M)[1:10]

factorA   <- as.numeric(My.Array.1$factorA)
factorB   <- as.numeric(My.Array.1$factorB)
factorC   <- rep(NA,each=dim(D.1)[1]);names(factorC)<- rownames(D.1)
factorC   <- as.numeric(factorC)
id        <- which(clinics.1$Sex=='M');factorC[id]<- -1
id        <- which(clinics.1$Sex=='F');factorC[id]<-  1
Sex       <- clinics.1$Sex
quant1    <- clinics.1$Age_at_LP;length(quant1)
A.12      <- (factorA+3)/2
B.12      <- (factorB+3)/2
C.12      <- (factorC+3)/2
D         <- data.frame(factorA,factorB,factorC,A.12,B.12,C.12,quant1)
rownames(D)<- rownames(M)
my.labels.sex <- C.12
my.cls.sex    <- c('skyblue','deeppink')[C.12]

quant.1.gr <- quantile(quant1, probs = seq(0, 1, 0.25))
id1        <- which(quant1<quant.1.gr[2])
id2        <- which(quant1<quant.1.gr[3])
id3        <- which(quant1<quant.1.gr[4])
id.q1 <- id1;                         sort(quant1[id.q1])
id.q2 <- setdiff(id2,id1);            sort(quant1[id.q2])
id.q3 <- setdiff(id3,c(id1,id2));     sort(quant1[id.q3])
id.q4 <- setdiff(c(1:dim(M)[1]),c(id1,id2,id3)); sort(quant1[id.q4])
length(c(id.q1,id.q2,id.q3,id.q4))

qu.gr <- rep(NA,each=dim(M)[1])
qu.gr[id.q1]<- 1
qu.gr[id.q2]<- 2
qu.gr[id.q3]<- 3
qu.gr[id.q4]<- 4

cls.4 <- c('brown','green4','skyblue','black')

#cls.4 <- palette("Classic Tableau")[1:4]
my.cls.qu.gr <- cls.4 [qu.gr]

D.all <- data.frame(as.vector(factorA), as.vector(factorB),as.vector(factorC),as.vector(quant1),as.vector(qu.gr),as.vector(qu.gr),as.matrix(D.1),Sex)
colnames(D.all)[1:5] <- c('factorA','factorB','factorC','quant1','qu.gr')
head(D.all)
dim(D.all)

#write.table(D.all,file='./Data/D.all.txt',sep='\t')


# Original data  Add sex and quant
####################################
My.Array    <- data.frame(
  M       =I(as.matrix(M)),
  D       =I(as.matrix(D)),
  factorA = factorA, # cluster ID
  factorB = factorB, # MS
  factorC = factorC, # Sex
  quant1 = as.numeric(quant1), # Age_at_LP
  qu.gr= factor(qu.gr))

# Permuted data
#################
id        <- shuffle(factorA);factorA.p <- factorA[id]
id        <- shuffle(factorB);factorB.p <- factorB[id]
id        <- shuffle(factorC);factorC.p <- factorC[id]
id        <- shuffle(quant1);quant1.p   <- quant1[id]
A.12.p    <- (factorA.p+3)/2
B.12.p    <- (factorB.p+3)/2
C.12.p    <- (factorC.p+3)/2
# #

My.Array.p    <- data.frame(
    M       =I(as.matrix(M)),
    factorA.p = factorA.p,
    factorB.p = factorB.p,
    factorC.p = factorC.p,
    quant1.p = quant1.p)

Permuted.data.1 <- list(My.Array.p,factorA.p,factorB.p,factorC.p,quant1.p)
#save(Permuted.data.1,file='./Data/Permuted.data.1.RData')
load(file='./Data/Permuted.data.1.RData')

My.Array.p <- Permuted.data.1[[1]]
factorA.p  <- Permuted.data.1[2][[1]]
factorB.p  <- Permuted.data.1[3][[1]]
factorC.p  <- Permuted.data.1[4][[1]]
quant1.p   <- Permuted.data.1[5][[1]]

##################################################################
# GEM data set 1,all proteins, to make eED and ED
##################################################################

# original data - original model
# gem.o       <- GEM(M ~ factorA*factorB, data = My.Array)

# original data - new model
gem        <- GEM(M ~ factorA+factorB+factorC+quant1, data = My.Array)
er.factorA <- gem$ER.values$factorA
er.factorB <- gem$ER.values$factorB
er.factorC <- gem$ER.values$factorC
er.quant1  <- gem$ER.values$quant1
eff.factorA <- gem$effects$factorA
eff.factorB <- gem$effects$factorB
eff.factorC <- gem$effects$factorC
eff.quant1  <- gem$effects$quant1
gem.res     <- gem$residuals


# original data, knockout factorA
er.knockout.A <- knock.out(gem, 'factorA')

# permuted data
gem.p     <- GEM(M ~ factorA.p+factorB.p+factorC.p+quant1.p, data = My.Array.p)

#####################################################################
# PLS original data
#####################################################################

# Model.1 (Mosleth et al. 2021): A*B and Model.2: A B C quant

## Model.1
####################

# factor A (Cluster ID)
###########################
my.ncomp    <- 8
y           <- factorA
pls.mod     <- pls(gem, 'factorA', ncomp= my.ncomp, validation = "LOO", df.used = 2)
error.A     <- check.error(pls.mod,ncomp=my.ncomp,y);error.A
PLS.scores.A    <- pls.mod$scores
PLS.loadings.A  <- pls.mod$loadings

# Jackknifed comp 5 (lowest error)
my.ncomp    <- 2
pls.mod     <- pls(gem, 'factorA', ncomp= my.ncomp, validation = "LOO", jackknife = TRUE, df.used = 3)
pls.mod.A   <- pls.mod
p.jkn       <- pls.mod$jack[,1,my.ncomp];
id          <- which(p.jkn<l)
n.s.jkn.A.comp5   <- colnames(M)[id];length(n.s.jkn.A.comp5)

# Jackknifed comp 2
my.ncomp    <- 2
pls.mod     <- pls(gem, 'factorA', ncomp= my.ncomp, validation = "LOO", jackknife = TRUE, df.used = 3)
pls.mod.A   <- pls.mod
p.jkn       <- pls.mod$jack[,1,my.ncomp];
id          <- which(p.jkn<l)
n.s.jkn.A.comp2   <- colnames(M)[id];length(n.s.jkn.A.comp2)
test              <- intersect(n.s.jkn.A.comp2,n.s.jkn.A.comp5);length(test)

n.s.jkn.A   <- n.s.jkn.A.comp2
coeff.A     <- pls.mod$coefficients[,1,my.ncomp]
id          <- which(coeff.A<0);neg.A <- names(coeff.A)[id]
id          <- which(coeff.A>0);pos.A <- names(coeff.A)[id]
n.s.jkn.A.neg <- intersect(n.s.jkn.A,neg.A)
n.s.jkn.A.pos <- intersect(n.s.jkn.A,pos.A)


# factor B (MS)
##################
my.ncomp    <- 8
y           <- factorB
pls.mod     <- pls(gem, 'factorB', ncomp= my.ncomp, validation = "LOO", df.used = 3)
error.B     <- check.error(pls.mod,ncomp=my.ncomp,y);error.B
PLS.scores.B    <- pls.mod$scores
PLS.loadings.B  <- pls.mod$loadings

# Jackknife
my.ncomp    <- 6
pls.mod     <- pls(gem, 'factorB', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
pls.mod.B   <- pls.mod
p.jkn       <- pls.mod$jack[,1,my.ncomp]
id          <- which(p.jkn<l)
n.s.jkn.B.comp6   <- colnames(M)[id];length(n.s.jkn.B.comp6)

my.ncomp    <- 2
pls.mod     <- pls(gem, 'factorB', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
pls.mod.B   <- pls.mod
p.jkn       <- pls.mod$jack[,1,my.ncomp]
id          <- which(p.jkn<l)
n.s.jkn.B.comp2   <- colnames(M)[id];length(n.s.jkn.B.comp2)
test              <- setdiff(n.s.jkn.B.comp6,n.s.jkn.B.comp2);length(test)

n.s.jkn.B   <- n.s.jkn.B.comp2
coeff.B     <- pls.mod$coefficients[,1,my.ncomp]
id          <- which(coeff.B<0);neg.B <- names(coeff.B)[id]
id          <- which(coeff.B>0);pos.B <- names(coeff.B)[id]
n.s.jkn.B.neg <- intersect(n.s.jkn.B,neg.B)
n.s.jkn.B.pos <- intersect(n.s.jkn.B,pos.B)

#test <- gconvert(n.s.jkn.B.comp2)[,1:6];id <- which(substr(test$description,1,5)=="compl");test[id,5]
#test <- gconvert(n.s.jkn.B.comp6)[,1:6];id <- which(substr(test$description,1,5)=="compl");test[id,5]

# factor C (Sex)
################
my.ncomp    <- 8
y           <- factorC
pls.mod     <- pls(gem, 'factorC', ncomp= my.ncomp, validation = "LOO", df.used = 3)
error.C     <- check.error(pls.mod,ncomp=my.ncomp,y);error.C
PLS.scores.C    <- pls.mod$scores
PLS.loadings.C  <- pls.mod$loadings

# Jackknife
my.ncomp    <- 3
y           <- factorC
pls.mod     <- pls(gem, 'factorC', ncomp= my.ncomp, validation = "LOO",jackknife = TRUE, df.used = 3)
pls.mod.C   <- pls.mod
p.jkn       <- pls.mod$jack[,1,my.ncomp]
id          <- which(p.jkn<l)
n.s.jkn.C.comp3   <- colnames(M)[id];length(n.s.jkn.C.comp3)

my.ncomp    <- 2
y           <- factorC
pls.mod     <- pls(gem, 'factorC', ncomp= my.ncomp, validation = "LOO",jackknife = TRUE, df.used = 3)
pls.mod.C   <- pls.mod
p.jkn       <- pls.mod$jack[,1,my.ncomp]
id          <- which(p.jkn<l)
n.s.jkn.C.comp2   <- colnames(M)[id];length(n.s.jkn.C.comp2)
test              <- intersect(n.s.jkn.C.comp2,n.s.jkn.C.comp3);length(test)

n.s.jkn.C   <- n.s.jkn.C.comp2
coeff.C     <- pls.mod$coefficients[,1,my.ncomp]
id          <- which(coeff.C<0);neg.C <- names(coeff.C)[id]
id          <- which(coeff.C>0);pos.C <- names(coeff.C)[id]
n.s.jkn.C.neg <- intersect(n.s.jkn.C,neg.C)
n.s.jkn.C.pos <- intersect(n.s.jkn.C,pos.C)

# quant1 (Age at LP)
##########################
my.ncomp    <- 8
y           <- quant1
pls.mod     <- pls(gem, 'quant1', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
pls.mod.quant1 <- pls.mod
# summary(pls.mod.quant1)
# Kristian jeg skrev dette inn selv, hvor finner jeg disse tallene i pls.mod
quant.1.RMSEP.CV    <- c(8.807,8.486,7.394,6.695,6.928,6.970,7.298,7.283,7.252)
quant.1.RMSEP.adj.CV<- c(8.807,8.486,7.372,6.681,6.901,6.944,7.266,7.252,7.221)

PLS.scores.quant1    <- pls.mod$scores
PLS.loadings.quant1  <- pls.mod$loadings

# Jackknifed coefficient P-values
my.ncomp      <- 3
pls.mod       <- pls(gem, 'quant1', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
p.jkn         <- pls.mod$jack[,1,my.ncomp]
id            <- which(p.jkn<l); length(id)
n.s.jkn.quant1.comp3   <- colnames(M)[id]

my.ncomp      <- 2
pls.mod       <- pls(gem, 'quant1', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
p.jkn         <- pls.mod$jack[,1,my.ncomp]
id            <- which(p.jkn<l); length(id)
n.s.jkn.quant1.comp2   <- colnames(M)[id]
test <- intersect(n.s.jkn.quant1.comp2,n.s.jkn.quant1.comp3);length(test)

n.s.jkn.quant1   <- n.s.jkn.quant1.comp2
coeff.quant1     <- pls.mod$coefficients[,1,my.ncomp]
residuals.quant1 <- pls.mod$residuals[,1,my.ncomp]
id <- which(coeff.quant1<0);neg.q <- names(coeff.quant1)[id]
id <- which(coeff.quant1>0);pos.q <- names(coeff.quant1)[id]
n.s.jkn.quant1.neg <- intersect(n.s.jkn.quant1,neg.q)
n.s.jkn.quant1.pos <- intersect(n.s.jkn.quant1,pos.q)

#####################################################################
# PLS permuted data
#####################################################################

# factor A.p
################
my.ncomp    <- 8
y           <- factorA.p
pls.mod     <- pls(gem.p, 'factorA.p', ncomp= my.ncomp, validation='LOO', df.used = 3)
error.A.p     <- check.error(pls.mod,ncomp=my.ncomp,y);error.A.p

# factor B.p
#############
my.ncomp    <- 8
y           <- factorB.p
pls.mod     <- pls(gem, 'factorB', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
error.B.p     <- check.error(pls.mod,ncomp=my.ncomp,y);error.B.p

# factor C.p
#############
my.ncomp    <- 8
y           <- factorC.p
pls.mod     <- pls(gem, 'factorC', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
error.C.p     <- check.error(pls.mod,ncomp=my.ncomp,y);error.C.p

# quant1
############
my.ncomp    <- 8
y           <- quant1.p
pls.mod     <- pls(gem.p, 'quant1.p', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 3)
pls.mod.quant1.p <- pls.mod
# summary(pls.mod.quant1.p)
# Kristian jeg skrev dette inn selv, hvor finner jeg disse tallene i pls.mod
quant.1.p.RMSEP.CV    <- as.numeric(c(8.807,8.932,9.696,9.769,9.338,9.666,9.146,9.346,9.199))
quant.1.p.RMSEP.adj.CV <- c(8.807,8.933,9.679,9.741,9.301,9.623,9.099,9.297,9.154)


################################################
## Plot results of original vs permuted model
################################################

pdf('./Figures/GEM.PLS Prediction results of data and permuted.data.pdf')

my.cls.model <- 'black'
my.cls.model.p <- 'blue'

par(mfrow=c(3,2))

dataset <- 'Group ID'
my.xlim <- c(1,12);my.ylim <- c(0,0.8)
plot(error.A,col='white',xlab='no of PLS factors',ylab='error rate',xlim=my.xlim,ylim=my.ylim)
lines(error.A,col=my.cls.model,xlim=my.xlim,ylim=my.ylim)
lines(error.A.p,col=my.cls.model.p,xlim=my.xlim,ylim=my.ylim,lty=2)
abline(h=0.5,lty=3)
text(10.2,0.085,label='Data',col=my.cls.model)
text(10.6,0.64,label='Permuted',col=my.cls.model.p)
title(paste0(letters[1],") Error rate ",dataset), adj = 0,line = 1)

dataset <- 'MS'
my.xlim <- c(1,12);my.ylim <- c(0,0.8)
plot(error.B,col='white',xlab='no of PLS factors',ylab='error rate',xlim=my.xlim,ylim=my.ylim)
lines(error.B,col=my.cls.model,xlim=my.xlim,ylim=my.ylim)
lines(error.B.p,col=my.cls.model.p,xlim=my.xlim,ylim=my.ylim,lty=2)
abline(h=0.5,lty=3)
text(10.2,0.085,label='Data',col=my.cls.model)
text(10.7,0.45,label='Permuted',col=my.cls.model.p)
title(paste0(letters[2],") Error rate ",dataset), adj = 0,line = 1)

dataset <- 'Sex'
my.xlim <- c(1,12);my.ylim <- c(0,0.8)
plot(error.C,col='white',xlab='no of PLS factors',ylab='error rate',xlim=my.xlim,ylim=my.ylim)
lines(error.C,col=my.cls.model,xlim=my.xlim,ylim=my.ylim)
lines(error.B.p,col=my.cls.model.p,xlim=my.xlim,ylim=my.ylim,lty=2)
abline(h=0.5,lty=3)
text(10.2,0.2,label='Data',col=my.cls.model)
text(10.7,0.45,label='Permuted',col=my.cls.model.p)
title(paste0(letters[3],") Error rate ",dataset), adj = 0,line = 1)

dataset <- 'Age at LP'
my.xlim <- c(1,12);my.ylim <- c(6,10)
plot(quant.1.p.RMSEP.adj.CV,col='white',xlim=my.xlim,ylim=my.ylim,xlab='no of PLS factors',ylab='RMSEP')
lines(quant.1.RMSEP.CV,xlim=my.xlim,ylim=my.ylim,col=my.cls.model)
#lines(quant.1.RMSEP.adj.CV,xlim=my.xlim,ylim=my.ylim,col=my.cls.model,lty=3)
lines(quant.1.p.RMSEP.CV,xlim=my.xlim,ylim=my.ylim,col=my.cls.model.p,lty=2)
#lines(quant.1.p.RMSEP.adj.CV,xlim=my.xlim,ylim=my.ylim,col=my.cls.model.p,lty=3)
text(10.5,7.2,label='Data')
text(10.7,9.2,label='Permuted',col='blue')
title(paste0(letters[4],") RMSEP ",dataset), adj = 0,line = 1)

dev.off()


################
# PLS Plots
################

pdf('./Figures/GEM.PLS scores and loadings A (Cluster), B(MS), C (Sex) and quant1(Age) cohort 1.pdf')

comps <- c(1,2)

par(mfrow=c(4,2))
my.cls.s    <- cls.2[D$A.12];# ctrl- MS
my.pch.s    <- points.2[D$B.12]

# factor A
dataset <- '\nGroup ID'
plot(PLS.scores.A[,comps],col=my.cls.s,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[1],") Scores ",dataset), adj = 0,line = 1)

my.xlim <- c(1.2*min(PLS.loadings.A[,1]),1.2*max(PLS.loadings.A[,1]))
my.ylim <- c(1.2*min(PLS.loadings.A[,2]),1.2*max(PLS.loadings.A[,2]))
v           <- rep(1,each=dim(PLS.loadings.A)[1]);length(v)
names(v)    <- rownames(PLS.loadings.A)
v[n.s.jkn.A] <- 2
my.cls.v    <- c('gray65','black')[v]
my.pch.v    <- c(1,15,17,19)[v]
plot(PLS.loadings.A[,comps],col=my.cls.v,pch=my.pch.v,
     xlim=my.xlim,ylim=my.ylim,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[2],") Loadings ",dataset), adj = 0,line = 1)

# factor B
dataset <- '\nMS'
plot(PLS.scores.B[,comps],col=my.cls.s,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[3],") Scores ",dataset), adj = 0,line = 1)

my.xlim <- c(1.2*min(PLS.loadings.B[,1]),1.2*max(PLS.loadings.B[,1]))
my.ylim <- c(1.2*min(PLS.loadings.B[,2]),1.2*max(PLS.loadings.B[,2]))
v           <- rep(1,each=dim(PLS.loadings.B)[1]);length(v)
names(v)    <- rownames(PLS.loadings.B)
v[n.s.jkn.B] <- 2
my.cls.v    <- c('gray65','black')[v]
my.pch.v    <- c(1,15,17,19)[v]
plot(PLS.loadings.B[,comps],col=my.cls.v,pch=my.pch.v,
     xlim=my.xlim,ylim=my.ylim,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[4],") Loadings ",dataset), adj = 0,line = 1)

dataset <- '\nSex'
plot(PLS.scores.C[,comps],col='white',
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
text(PLS.scores.C[,comps],col=my.cls.sex,labels=my.labels.sex)
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[5],") Scores ",dataset), adj = 0,line = 1)

my.xlim <- c(1.2*min(PLS.loadings.C[,1]),1.2*max(PLS.loadings.C[,1]))
my.ylim <- c(1.2*min(PLS.loadings.C[,2]),1.2*max(PLS.loadings.C[,2]))
v           <- rep(1,each=dim(PLS.loadings.C)[1]);length(v)
names(v)    <- rownames(PLS.loadings.C)
v[n.s.jkn.C] <- 2
my.cls.v    <- c('gray65','black')[v]
my.pch.v    <- c(1,15,17,19)[v]
plot(PLS.loadings.C[,comps],col=my.cls.v,pch=my.pch.v,
     xlim=my.xlim,ylim=my.ylim,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[6],") Loadings ",dataset), adj = 0,line = 1)


# quant1
dataset <- '\nAge at LP'
plot(PLS.scores.quant1[,comps],col='white',
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
text(PLS.scores.quant1[,comps],col=my.cls.qu.gr,labels=quant1)
title(paste0(letters[7],") Scores ",dataset), adj = 0,line = 1)

my.xlim <- c(1.2*min(PLS.loadings.quant1[,1]),1.2*max(PLS.loadings.quant1[,1]))
my.ylim <- c(1.2*min(PLS.loadings.quant1[,2]),1.2*max(PLS.loadings.quant1[,2]))
v           <- rep(1,each=dim(PLS.loadings.quant1)[1]);length(v)
names(v)    <- rownames(PLS.loadings.quant1)
v[n.s.jkn.quant1] <- 2
my.cls.v    <- c('gray65','black')[v]
my.pch.v    <- c(1,15,17,19)[v]
plot(PLS.loadings.quant1[,comps],col=my.cls.v,pch=my.pch.v,
     xlim=my.xlim,ylim=my.ylim,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[8],") Loadings ",dataset), adj = 0,line = 1)

dev.off()

#############################################################
## Find name and desc of significant proteins
#############################################################

# consistent for factor A withint both levels of factor B and vice versa

# factor A
##############################

# A pos
N.sel       <- n.s.jkn.A.pos;length(N.sel)
GN          <- gconvert(N.sel)
my.gr.names <- c('Gr1.nonMS','Gr1.MS','Gr2.nonMS','Gr2.MS')
Data.sel    <- M[,N.sel]
gr4         <- D.1$gr4
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
cons.31             <- which(Means.gr4[3,]>Means.gr4[1,]);
cons.42             <- which(Means.gr4[4,]>Means.gr4[2,]);
names.consistent    <- intersect(names(cons.31),names(cons.42));length(names.consistent)
y                   <- cbind(names.consistent,names.consistent);colnames(y)<- c('input','input.copy')
x                   <- GN
xy                  <- merge(x,y,by.x='input',by.y='input');dim(x);dim(xy)
GN.consistent       <- xy
description         <- sub('\\[.*', '', GN.consistent$description)

GN.s.jkn.consistent.A.pos   <- data.frame(GN.consistent[,c('input','target','name')],description)
unique.N.consistent         <- unique(xy$input)
Means.gr4.consistent.A.pos  <- Means.gr4[,unique.N.consistent];
M.consistent.A.pos          <- M[,unique.N.consistent];
dim(GN.s.jkn.consistent.A.pos)
dim(Means.gr4.consistent.A.pos)


# A neg
N.sel       <- n.s.jkn.A.neg;length(N.sel)
GN          <- gconvert(N.sel)
my.gr.names <- c('Gr1.nonMS','Gr1.MS','Gr2.nonMS','Gr2.MS')
Data.sel    <- M[,N.sel]
gr4         <- D.1$gr4
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
cons.31             <- which(Means.gr4[3,]<Means.gr4[1,]);
cons.42             <- which(Means.gr4[4,]<Means.gr4[2,]);
names.consistent    <- intersect(names(cons.31),names(cons.42));length(names.consistent)
y                   <- cbind(names.consistent,names.consistent);colnames(y)<- c('input','input.copy')
x                   <- GN
xy                  <- merge(x,y,by.x='input',by.y='input');dim(x);dim(xy)
GN.consistent       <- xy
description         <- sub('\\[.*', '', GN.consistent$description)

GN.s.jkn.consistent.A.neg   <- data.frame(GN.consistent[,c('input','target','name')],description)
unique.N.consistent         <- unique(xy$input)
Means.gr4.consistent.A.neg  <- Means.gr4[,unique.N.consistent];
M.consistent.A.neg          <- M[,unique.N.consistent];
dim(GN.s.jkn.consistent.A.neg)
dim(Means.gr4.consistent.A.neg)


# factor B
##############################

# B pos
N.sel       <- n.s.jkn.B.pos;length(N.sel)
GN          <- gconvert(N.sel)
my.gr.names <- c('Gr1.nonMS','Gr1.MS','Gr2.nonMS','Gr2.MS')
Data.sel    <- M[,N.sel]
gr4         <- D.1$gr4
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
cons.21             <- which(Means.gr4[2,]>Means.gr4[1,]);#cons.12 # MS - nonD within Cluster 1
cons.43             <- which(Means.gr4[4,]>Means.gr4[3,]);#cons.43 # MS - nonD within Cluster 2
names.consistent    <- intersect(names(cons.21),names(cons.43));length(names.consistent)
y                   <- cbind(names.consistent,names.consistent);colnames(y)<- c('input','input.copy')
x                   <- GN
xy                  <- merge(x,y,by.x='input',by.y='input');dim(x);dim(xy)
GN.consistent       <- xy
description         <- sub('\\[.*', '', GN.consistent$description)

GN.s.jkn.consistent.B.pos   <- data.frame(GN.consistent[,c('input','target','name')],description)
unique.N.consistent         <- unique(xy$input)
Means.gr4.consistent.B.pos  <- Means.gr4[,unique.N.consistent];
M.consistent.B.pos          <- M[,unique.N.consistent];
dim(GN.s.jkn.consistent.B.pos)
dim(Means.gr4.consistent.B.pos)


# B neg
N.sel       <- n.s.jkn.B.neg;length(N.sel)
GN          <- gconvert(N.sel)
my.gr.names <- c('Cl1.nonMS','Cl1.MS','Cl2.nonMS','Cl2.MS')
Data.sel    <- M[,N.sel]
gr4         <- D.1$gr4
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
cons.21             <- which(Means.gr4[2,]<Means.gr4[1,]);#cons.12 # MS - nonD within Cluster 1
cons.43             <- which(Means.gr4[4,]<Means.gr4[3,]);#cons.43 # MS - nonD within Cluster 2
names.consistent    <- intersect(names(cons.21),names(cons.43));length(names.consistent)
y                   <- cbind(names.consistent,names.consistent);colnames(y)<- c('input','input.copy')
x                   <- GN
xy                  <- merge(x,y,by.x='input',by.y='input');dim(x);dim(xy)
GN.consistent       <- xy
description         <- sub('\\[.*', '', GN.consistent$description)
GN.s.jkn.consistent.B.neg   <- data.frame(GN.consistent[,c('input','target','name')],description)
unique.N.consistent         <- unique(xy$input)
Means.gr4.consistent.B.neg  <- Means.gr4[,unique.N.consistent];
M.consistent.B.neg          <- M[,unique.N.consistent];
dim(GN.s.jkn.consistent.B.neg)
dim(Means.gr4.consistent.B.neg)


# sign for B and for C
########################################

# B.consistent.pos C.pos
N.sel           <- intersect(GN.s.jkn.consistent.B.pos$input,n.s.jkn.C.pos);length(N.sel)
GN              <- gconvert(N.sel)
description     <- sub('\\[.*', '', GN$description)
GN.s.jkn.B.pos.and.C.pos  <- data.frame(GN[,c('input','target','name')],description)
Means.gr4.consistent.B.pos.C.pos <- Means.gr4.consistent.B.pos[,GN.s.jkn.B.pos.and.C.pos$input]
colnames(Means.gr4.consistent.B.pos.C.pos) <- GN.s.jkn.B.pos.and.C.pos$name

# B.consistent.pos C.neg  - o
# N.sel           <- intersect(GN.s.jkn.consistent.B.pos$input,n.s.jkn.C.neg);length(N.sel)
# GN              <- gconvert(N.sel)
# description     <- sub('\\[.*', '', GN$description)
# GN.s.jkn.B.pos.and.C.neg  <- data.frame(GN[,c('input','target','name')],description)

# B.consistent.neg C.neg
N.sel           <- intersect(GN.s.jkn.consistent.B.neg$input,n.s.jkn.C.neg);length(N.sel)
GN              <- gconvert(N.sel)
description     <- sub('\\[.*', '', GN$description)
GN.s.jkn.B.neg.and.C.neg  <- data.frame(GN[,c('input','target','name')],description)
Means.gr4.consistent.B.neg.C.neg <- Means.gr4.consistent.B.neg[,GN.s.jkn.B.neg.and.C.neg$input]
colnames(Means.gr4.consistent.B.neg.C.neg) <- GN.s.jkn.B.neg.and.C.neg$name

# B.consistent.neg C.pos
N.sel           <- intersect(GN.s.jkn.consistent.B.neg$input,n.s.jkn.C.pos);length(N.sel)
GN              <- gconvert(N.sel)
description     <- sub('\\[.*', '', GN$description)
GN.s.jkn.B.neg.and.C.pos  <- data.frame(GN[,c('input','target','name')],description)
Means.gr4.consistent.B.neg.C.pos <- Means.gr4.consistent.B.neg[,GN.s.jkn.B.neg.and.C.pos$input]
colnames(Means.gr4.consistent.B.neg.C.pos) <- GN.s.jkn.B.neg.and.C.pos$name

# B.consistent.neg C.neg
N.sel           <- intersect(GN.s.jkn.consistent.B.neg$input,n.s.jkn.C.neg);length(N.sel)
GN              <- gconvert(N.sel)
description     <- sub('\\[.*', '', GN$description)
GN.s.jkn.B.neg.and.C.neg  <- data.frame(GN[,c('input','target','name')],description)
Means.gr4.consistent.B.neg.C.neg <- Means.gr4.consistent.B.neg[,GN.s.jkn.B.neg.and.C.neg$input]
colnames(Means.gr4.consistent.B.neg.C.neg) <- GN.s.jkn.B.neg.and.C.neg$name

# D quant
N.sel           <- n.s.jkn.quant1.neg
GN              <- gconvert(N.sel)
description     <- sub('\\[.*', '', GN$description)
GN.s.jkn.D.neg  <- data.frame(GN[,c('input','target','name')],description)

N.sel           <- n.s.jkn.quant1.pos
GN              <- gconvert(N.sel)
description     <- sub('\\[.*', '', GN$description)
GN.s.jkn.D.pos  <- data.frame(GN[,c('input','target','name')],description)


# factor A C and B
##############################

# B pos.C.pos
#################
N.sel       <- GN.s.jkn.B.pos.and.C.pos$input;length(N.sel)
P.sel       <- GN.s.jkn.B.pos.and.C.pos$name;length(P.sel)
GN          <- gconvert(N.sel)
my.gr.names <- c('Gr1.nonMS','Gr1.MS','Gr2.nonMS','Gr2.MS')
id.F        <- which(clinics.1$Sex=='F');id.F
id.M        <- which(clinics.1$Sex=='M');id.M

# Females
Data.sel    <- M[id.F,N.sel]
gr4         <- D.1[id.F,'gr4']
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
colnames(Means.gr4) <- P.sel
F.Means.gr4.B.pos.and.C.pos <-Means.gr4

# Males
Data.sel    <- M[id.M,N.sel]
gr4         <- D.1[id.M,'gr4']
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- Data.sel[id,] # only one
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
colnames(Means.gr4) <- P.sel
M.Means.gr4.B.pos.and.C.pos <-Means.gr4

# B.neg.C.pos
##################
N.sel       <- GN.s.jkn.B.neg.and.C.pos$input;length(N.sel)
P.sel       <- GN.s.jkn.B.neg.and.C.pos$name;length(P.sel)
GN          <- gconvert(N.sel)
my.gr.names <- c('Gr1.nonMS','Gr1.MS','Gr2.nonMS','Gr2.MS')
id.F        <- which(clinics.1$Sex=='F');id.F
id.M        <- which(clinics.1$Sex=='M');id.M

# Females
Data.sel    <- M[id.F,N.sel]
gr4         <- D.1[id.F,'gr4']
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
colnames(Means.gr4) <- P.sel
F.Means.gr4.B.neg.and.C.pos <-Means.gr4

# Males
Data.sel    <- M[id.M,N.sel]
gr4         <- D.1[id.M,'gr4']
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- Data.sel[id,] # only one
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
colnames(Means.gr4) <- P.sel
M.Means.gr4.B.neg.and.C.pos <-Means.gr4

# B.neg.C.neg
##################
N.sel       <- GN.s.jkn.B.neg.and.C.neg$input;length(N.sel)
P.sel       <- GN.s.jkn.B.neg.and.C.neg$name;length(P.sel)
GN          <- gconvert(N.sel)
my.gr.names <- c('Gr1.nonMS','Gr1.MS','Gr2.nonMS','Gr2.MS')

#write.table(GN.s.jkn.B.neg.and.C.neg,file='./Results/GEM/GN.s.jkn.B.neg.and.C.neg.txt',sep='\t')

# Females
Data.sel    <- M[id.F,N.sel]
gr4         <- D.1[id.F,'gr4']
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
colnames(Means.gr4) <- P.sel
F.Means.gr4.B.neg.and.C.neg <-Means.gr4

# Males
Data.sel    <- M[id.M,N.sel]
gr4         <- D.1[id.M,'gr4']
id          <- which(gr4==1);mean.gr1 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==2);mean.gr2 <- Data.sel[id,] # only one
id          <- which(gr4==3);mean.gr3 <- apply(Data.sel[id,],2,mean)
id          <- which(gr4==4);mean.gr4 <- apply(Data.sel[id,],2,mean)
Means.gr4   <- rbind(mean.gr1,mean.gr2,mean.gr3,mean.gr4)
rownames(Means.gr4) <- my.gr.names
colnames(Means.gr4) <- P.sel
M.Means.gr4.B.neg.and.C.neg <-Means.gr4


#############################
# save to table
#############################

dim(GN.s.jkn.B.neg.and.C.neg)
dim(GN.s.jkn.B.neg.and.C.pos)
dim(GN.s.jkn.B.pos.and.C.pos)
#dim(GN.s.jkn.B.pos.and.C.neg)

write.table(GN.s.jkn.B.neg.and.C.neg,file='./Results/GN.s.jkn.B.neg.and.C.neg.txt',sep='\t')
write.table(GN.s.jkn.B.neg.and.C.pos,file='./Results/GN.s.jkn.B.neg.and.C.pos.txt',sep='\t')
write.table(GN.s.jkn.B.pos.and.C.pos,file='./Results/GN.s.jkn.B.pos.and.C.pos.txt',sep='\t')

write.table(GN.s.jkn.D.pos,file='./Results/GN.s.jkn.D.pos.txt',sep='\t')
write.table(GN.s.jkn.D.neg,file='./Results/GN.s.jkn.D.neg.txt',sep='\t')



#############################
# PCA of selected proteins
#############################

pdf('./Figures/PCA.original data after GEM.PLS.jkn. A (Cluster), B (MS), C (Sex) and quant1(Age) cohort 1.pdf')

par(mfrow=c(4,2))

# factor A
Inn             <- M[,n.s.jkn.A]
pca.mod         <- prcomp(Inn,scale=FALSE) # all data are already scales
PCA.scores.A    <- pca.mod$x
PCA.loadings.A  <- pca.mod$rotation
expl.var        <- round(propVar(pca.mod)*100,digits = 0);expl.var
dataset         <- 'factor A (Cluster ID)'

comps <- c(1,2)
plot(PCA.scores.A[,comps],col=my.cls.s,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[1],") Scores ",dataset), adj = 0,line = 1)

plot(PCA.loadings.A[,comps],
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
#text(PCA.loadings.A[,comps],labels=colnames(Inn),cex=0.8)
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[2],") Loadings ",dataset), adj = 0,line = 1)


# factor B
Inn             <- M[,n.s.jkn.B]
pca.mod         <- prcomp(Inn,scale=FALSE) # all data are already scales
PCA.scores.B    <- pca.mod$x
PCA.loadings.B  <- pca.mod$rotation
expl.var        <- round(propVar(pca.mod)*100,digits = 0);expl.var
dataset         <- 'factor B (MS)'

comps <- c(1,2)
plot(PCA.scores.B[,comps],col=my.cls.s,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[1],") Scores ",dataset), adj = 0,line = 1)

plot(PCA.loadings.B[,comps],
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[2],") Loadings ",dataset), adj = 0,line = 1)

# factor C
Inn             <- M[,n.s.jkn.C]
pca.mod         <- prcomp(Inn,scale=FALSE) # all data are already scales
PCA.scores.C    <- pca.mod$x
PCA.loadings.C  <- pca.mod$rotation
expl.var        <- round(propVar(pca.mod)*100,digits = 0);expl.var
dataset         <- 'factor C (Sex)'

comps <- c(1,2)
plot(PCA.scores.C[,comps],col=my.cls.sex,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
text(PCA.scores.C[,comps],col=my.cls.sex,labels=my.labels.sex)
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[1],") Scores ",dataset), adj = 0,line = 1)

plot(PCA.loadings.C[,comps],
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[2],") Loadings ",dataset), adj = 0,line = 1)

# quant1
# factor B
Inn             <- M[,n.s.jkn.quant1]
pca.mod         <- prcomp(Inn,scale=FALSE) # all data are already scales
PCA.scores.quant1    <- pca.mod$x
PCA.loadings.quant1  <- pca.mod$rotation
expl.var        <- round(propVar(pca.mod)*100,digits = 0);expl.var
dataset         <- 'quant1 (Age)'

comps <- c(1,2)
plot(PCA.scores.quant1[,comps],col=my.cls.qu.gr,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[3],") Scores ",dataset), adj = 0,line = 1)

plot(PCA.loadings.quant1[,comps],
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title(paste0(letters[4],") Loadings ",dataset), adj = 0,line = 1)

dev.off()



######################################
## Plot selected proteins
######################################

#GN.sel.B <- rbind(GN.s.jkn.B.neg,GN.s.jkn.B.pos)

pdf('./Figures/Line plots of Means - Complement proteins, and Im kappa and heavy.pdf')

par(mfcol=c(2,2))
my.xlim <- c(0,6);my.ylim <- c(-1,2)
my.cex.main <- 1.1
gr4 <- D.1$gr4;length(gr4)

# Complements
GN        <- GN.s.jkn.consistent.B.neg
id        <- which(substr(GN$description,1,5)=='compl');id
GN.compl  <- GN[id,]
Inn       <- Means.gr4.consistent.B.neg[,GN.compl$input]
ER.sel      <- er.factorB[,colnames(Inn)]
ER.sel.mean <- aggregate(ER.sel,list(gr4),mean)[,-1];
rownames(ER.sel.mean)<- rownames(Inn)
k         <- dim(Inn)[2];k

plot(Inn[1:2,1],xlab='',xaxt='n',ylab='exprs.level',
     xlim=my.xlim,ylim=my.ylim,col='white')
axis(1,at=c(1,2,4,5),my.gr.names,las=3)
for (i in 1:k){lines(Inn[1:2,i],col='blue')}
for (i in 1:k){segments(x0=4, y0=Inn[3,i], x1 = 5, y1 = Inn[4,i],col='red')}
title(paste0(letters[1],") Complement proteins \nOriginal centred data\nsignifikantly neg for factor B (MS) "), adj = 0,line = 1,cex.main=my.cex.main)

Inn <- ER.sel.mean
plot(Inn[1:2,1],xlab='',xaxt='n',ylab='exprs.level',
     xlim=my.xlim,ylim=my.ylim,col='white')
axis(1,at=c(1,2,4,5),my.gr.names,las=3)
for (i in 1:k){lines(Inn[1:2,i],col='blue')}
for (i in 1:k){segments(x0=4, y0=Inn[3,i], x1 = 5, y1 = Inn[4,i],col='red')}
title(paste0(letters[2],") Complement proteins \nER values of MS\nsignifikantly neg for factor B (MS) "), adj = 0,line = 1,cex.main=my.cex.main)

dev.off()

