# MS.CIS 01 running script for SR-Publication.r

#setwd("H:/A_R/R/R data/R a run ELMO/E.MOD/MS/Externe data/Stoop 2017")

setwd("H:/A_R/R/R data/R MS publ.SR")
getwd()

#setwd("C:/Users/ellen.fargestad/Documents/R copy/MS/R MS publ.SR")
#getwd()

library(corrplot)
library(ER)
library(pls)
library(plsVarSel)
library(readxl)
library(corrplot)
library(ffmanova)
library(gprofiler2)
library(corrgram)

rm(list = ls())

####################
## Func
####################
propVar <- function(object){
  vars <- object$sdev^2
  vars/sum(vars)
}

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
  return(error) 
}

##############################################
# Set color and points
##############################################

cls.2     <- c('blue','red')
cls.3     <- c('blue','blue','red')
points.2  <- c(2,15)
points.3  <- c(1,2,15)

##############################################
# Get Data
##############################################

#load("H:/A_R/R/R data/R a run ELMO/E.MOD/MS/Results/G.all.RData")
#G.all[1:4,] . Here is not all, 
#dim(G.all)

Imp <- as.data.frame(read_excel("Data/Info2R. from Table.1.ALL proteins with.gene.names.xlsx", sheet = "PG.2.R"))
rownames(Imp)<- Imp[,1]
PG.id <- Imp
rm(Imp)

PG.id[1:3,]

out <- sort(substr(PG.id$names.short,1,6))
out

## Cohort 1
##################

# Patient data incl cons IgG
Imp <- as.data.frame(read_excel("Rawdata/Table S1. Cohort 1. Patient information.xlsx", sheet = "PatientData2R"))
rownames(Imp) <- Imp[,1]
DF.patients <- Imp

# Karim's data
Imp.K  <- as.data.frame(read_excel("Rawdata/Cohort 1. Table S2. PCA loadings (scores displayed as figs).xlsx", sheet = "(c) MS.loadings"))
n.K.ms <- Imp.K[,1]
n.K.sub4 <- c('P01834', 'P01625', 'P01871', 'P51693')
n.K.ms9  <-c('P06331', 'P01591', 'Q15782', 'P01857', 'Q99983', 'P02763', 'P00739', 'Q96ID5', 'Q9Y5Y7') 

# Elisa

# Proteome
Imp <- as.data.frame(read_excel("RawData/Cohort 1. Design.xlsx"))
rownames(Imp)<- Imp[,1]
D.1 <- Imp[,-1]
rm(Imp)

# CNS normalised data from Opsahl et al 2016
# CNS median normaliserte, log2 transformed. imputed 
# Statistisk organisering av data fxr protein profilering:, from Jill 2. februar 2016 kl. 15:12
# 1) log2
# 2) imputere missing values (setter inn protein intensitet basert pe normalfordeling). 
# 3) z-scorer proteiner: trekker protein intensiteten fra average intensitet.

Imp <- as.data.frame(
  read_excel("Rawdata/Cohort 1. Til Ellen_Eksoprt fra perseus_log2_imput_zscore.xlsx", sheet = "log2_impu"))
idc1  <- which(colnames(Imp)=='nonMS1')
idc2  <- which(colnames(Imp)=='MS37')
Temp  <- Imp[,idc1:idc2]
rownames(Temp)<-Imp$Accession
M.1.temp <- t(Temp);

# Protein info from the same file as the proteins
M.1.info      <- Imp[,(idc2+1):dim(Imp)[2]]
rownames(M.1.info)<- M.1.info$Accession
rm(Imp)

# Patient info
Imp <- as.data.frame(read_excel("Rawdata/Cohort 1. pmic12254-sup-0002-table1 Patient info.xlsx",sheet = "individual patient info"))
dim(Imp)
rownames(Imp)<- Imp[,1]
C.1.temp <- Imp[,-1]
rm(Imp)

# Reorder according to the order of the Design matrix
M.1 <- scale(M.1.temp[rownames(D.1),])
M.1.notscaled <- M.1.temp[rownames(D.1),]
clinics.1 <- C.1.temp[rownames(D.1),]
rm(C.1.temp);rm(M.1.temp);rm(Temp);rm(idc1);rm(idc2)

## Cohort 2
##################

# Design and Clinical data
Imp <- as.data.frame(read_excel("Rawdata/Cohort 2. Design, Clinics and Peptides From prca1891-sup-0002-suppmat2.xlsx", sheet = "D.Clinics"))
rownames(Imp)<- Imp[,1]
# reorder
# temp1 <- substr(rownames(Imp),3,3);temp1
# temp2 <- Imp$Conversion_to_CDMS;temp2
# temp3 <- Imp$Gender;temp3
# temp123 <- paste0(temp1,temp2,temp3)
# o <- order(temp123);o
# temp123[o]

C.2.temp <- Imp[,-1]

# Peptid-protein data from Stoop 2017 - Here: take means over several reads of the same protiens 
Imp <- as.data.frame(read_excel("Rawdata/Cohort 2. Design, Clinics and Peptides From prca1891-sup-0002-suppmat2.xlsx", sheet = "Peptides"))

#Imp <- as.data.frame(read_excel("H:/A_R/R/R data/R a run ELMO/E.MOD/MS/Externe data/Stoop 2017/Data/prca1891-sup-0002-suppmat2 . Stoop 2017 CIS vs CTRL.Ed.Ellen.xlsx", sheet = "ed peptides"))
#Imp[1:4,1:5]
idc1        <- which(colnames(Imp)=="CIS.10195");idc1
idc2        <- which(colnames(Imp)=="CTRL.A136");idc2
Ass.name    <- Imp[,1]
Data        <- Imp[,idc1:dim(Imp)[2]]
Info.2      <- Imp[,1:(idc1-1)];Info.2[1:3,]
Ass.Data    <- data.frame(Ass.name,Data)
Agg.Data    <- aggregate(Ass.Data[,-1], by = list(Ass.Data$Ass.name),mean, na.rm=TRUE)
rownames(Agg.Data)<- Agg.Data[,1]
DF2.with.NA <- Agg.Data[,-1];
DF2         <- DF2.with.NA[complete.cases(DF2.with.NA),]
dim(DF2);dim(DF2.with.NA)   # 51 out of 469 proteins have NA and are omitted in the file DF2
DF2.log     <- t(log2(DF2)) # log 2 transformation. there is no small numbers #min(Proteins.set2)
M.2.temp    <- as.matrix(DF2.log[rownames(C.2.temp),]) # order as Desing and clinical data

# Make design of new data with factors A and B
v         <- rep(NA, each=dim(C.2.temp)[1]);names(v) <- rownames(C.2.temp);length(v)
A         <- v; B <- v;names(A) 
A         <- as.numeric(A);B <- as.numeric(B)
A.12      <- A;B.12<- B
idr.ctrl  <- which(substr(rownames(C.2.temp),1,3)=='CTR');idr.ctrl
idr.cis   <- which(substr(rownames(C.2.temp),1,3)=='CIS');idr.cis
B[idr.ctrl] <- -1;B.12[idr.ctrl] <- 1;
B[idr.cis]  <-  1;B.12[idr.cis]  <- 2;

idr.neg   <- which(C.2.temp$Mutliple_oligoclonal_bands=='Negative');length(idr.neg)
idr.pos   <- which(C.2.temp$Mutliple_oligoclonal_bands=='Positive');length(idr.pos)
idr.unk   <- which(C.2.temp$Mutliple_oligoclonal_bands=='Unknown'); length(idr.unk)
A[1:length(A)]  <- -1;A.12[1:length(A)]   <- 1
A[idr.neg]      <- -1;A.12[idr.neg] <- 1;
A[idr.pos]      <-  1;A.12[idr.pos] <- 2;
A[idr.unk]      <-  NA;A.12[idr.unk] <- NA;

D.2.temp        <- data.frame(A,B,A.12,B.12)
rownames(D.2.temp)<- rownames(C.2.temp)
id <- order(paste0(D.2.temp$A,D.2.temp$B));names.order <- rownames(D.2.temp)[id]

clinics.2     <- C.2.temp[names.order,]
D.2           <- D.2.temp[names.order,]
M.2           <- scale(M.2.temp[names.order,])
M.2.notscaled <- M.2.temp[names.order,]

dim(D.2)
id      <- which(D.1$A.12==1);n1.a1 <- rownames(D.1)[id];length(id)
id      <- which(D.1$A.12==2);n1.a2 <- rownames(D.1)[id];length(id)
M.1.a1  <- M.1[n1.a1,]
M.1.a2  <- M.1[n1.a2,]
dim(M.1.a1)

# in Dataset 2. Stoop, use only a1 (subgr 1), CIS without Ig bands
id2.a1  <- which(D.2$A.12==1);length(id2.a1)
id2.a2  <- which(D.2$A.12==2)
M.2.a1  <- scale(M.2[id2.a1,]) # scale within a1
M.2.a1.notscales <- M.2.notscaled[id2.a1,]
M.2.a2  <- scale(M.2[id2.a2,]) # scale within a1
M.2.a2.notscales <- M.2.notscaled[id2.a2,]

C.2.a1 <- clinics.2[id2.a1,]
D.2.a1 <- D.2[id2.a1,]
D.2.a2 <- D.2[id2.a2,]
dim(D.2.a1)

rm(Imp);rm(DF2.with.NA);rm(Agg.Data);rm(DF2);rm(Data);rm(Ass.Data);rm(DF2.log)
rm(D.2.temp);rm(C.2.temp);rm(M.2.temp)

# common proteins of M.1 and M.2
common.proteins <- intersect(colnames(M.1),colnames(M.2))

##############################################
# Cohort 1, box plot of consIgG
#########################################################
# we do not have IgG numbers of non-neurological patients

DF.noNA <- na.omit(DF.patients)
dim(DF.noNA)
id <- which(DF.noNA$`IgG concentration_g.L`>0.5)
DF.noNA.2 <- DF.noNA[-id,]

cons.IgG  <- DF.noNA.2$`IgG concentration_g.L`
gr        <- DF.noNA.2$gr.4
my.labels <- c('1.C','1.MS','2.C','2.MS')
my.gr     <- my.labels[gr]
my.ylim   <- c(-0.05,0.22)

boxplot(cons.IgG~my.gr,col=c('blue','blue','red','red'),
        xlab='',xaxt='n',density=10,
        ylim=my.ylim)
#axis(1,at=c(1,2,3,4),labels =my.labels )
abline(v=2.5,lty=2,col='gray')
text(1.5,-0.03,labels='group 1',col='blue')
text(3.5,-0.03,labels='group 2',col='red')
text(1,-0.015,labels='Controls',col='blue')
text(2,-0.015,labels='MS',   col='blue')
text(3,-0.015,labels='Controls',col='red')
text(4,-0.015,labels='MS',   col='red')

# We can not calculate:
# IgG index was calculated as the CSF-plasma concentration quotient for IgG divided by the CSF-plasma concentration quotient for albumin (QIgG/ Qalb) 

##############################################
# Proteins selected by Opsahl
##############################################

Imp   <- as.data.frame(read_excel("Rawdata/Opsahl.pmic12254-sup-0003-table2.xlsx",sheet = "Table S2.p0.05"))
rownames(Imp)   <- Imp[,1]
Opsahl.Tab.S2   <- Imp
n.Opsahl.Tab.S2 <- rownames(Opsahl.Tab.S2)
Temp            <- PG.id[rownames(Opsahl.Tab.S2),];dim(Temp)

# Opsahl Fig.S1 proteins involved in neurogenesis and axogenesis
n.Opsahl.Fig.S1<- c('Q9Y6N8','P55290','P43146','P54764','Q15375','Q92876',
                    'P32004','Q92823','P04156','O75093','Q8IW52','O94991',
                    'P04216','Q9P2S2')
Temp  <- PG.id[n.Opsahl.Fig.S1,]


####################################################################
## ffmanova ignoring groups across the data- similar to Opsahl
####################################################################

# l       <- 0.05
# R       <- M.1 ;dim(R)
# B       <- D.1$B
# my.nSim <- 9999
# my.data     <- data.frame(B,R);dim(my.data)
# FF.res.1    <- ffmanova(R ~ B, data = my.data,stand=TRUE, nSim = my.nSim, verbose=TRUE)
# FF.res      <- FF.res.1
# FF.res
# p.Raw       <- FF.res$pRaw
# p.FDR       <- FF.res$pAdjFDR
# p.values.1  <- rbind(p.Raw,p.FDR);dim(p.values.1)
# rownames(p.values.1)<- c('int','p.Raw.B','int','p.FDR.B')
# id.raw      <- which(p.values.1[2,]<l);length(id.raw)
# id.fdr      <- which(p.values.1[4,]<l);length(id.fdr)
# n.fdr       <- colnames(p.values.1)[id.fdr]
# n.raw       <- colnames(p.values.1)[id.raw]

# of the 14 proteins in Opsahl fig S2, 11 is sign by RAW p-values, 3 by FDR adj pvalues


##############################################
# Make arrays
##############################################

M.1.c <- M.1[,common.proteins]
M.2.c <- M.2[,common.proteins]
M.2.c.a1 <- (M.2.c[id2.a1,])
M.2.c.a2 <- (M.2.c[id2.a2,])

# Make arrays
My.Array.1 <- data.frame(
  Design = I(D.1),
  factorA = D.1$A,
  factorB = D.1$B,
  clinics = I(clinics.1),
  M.1     = I(M.1),# M.log2.imputed, subtracted protein intensity from average intensity, scales
  M.1.c   = I(M.1.c))
names(My.Array.1)

My.Array.2 <- data.frame(
  Design = I(D.2),
  factorA = D.2$A,
  factorB = D.2$B,
  clinics = I(clinics.2),
  M.2     = I(M.2),# M.log2.scaled
  M.2.c   = I(M.2.c))
names(My.Array.2)

# Cohort 1 Subset. compare MS with healthy within gr 1
idr10       <- which(D.1$gr==10);idr10 # within gr 1, without neuroligical disorders
idr12       <- which(D.1$gr==12);idr12 # within gr 1, MS
M.1.subset.1  <- M.1[c(idr10,idr12),]
D.1.subset.1  <- D.1[c(idr10,idr12),]
factorB.subset.1 <- D.1.subset.1$B

My.Array.1.subset.1  <- data.frame(
  factorB.subset.1 = factorB.subset.1, 
  M.1.subset.1 = M.1.subset.1) # is scales

# Cohort 1 Subset. compare MS with patients with outher neurological disorders within gr 1
idr11       <- which(D.1$gr==11);idr11 # within gr 1, without neuroligical disorders
idr12       <- which(D.1$gr==12);idr12 # within gr 1, MS
M.1.subset.2  <- M.1[c(idr11,idr12),]
D.1.subset.2  <- D.1[c(idr11,idr12),]
factorB.subset.2 <- D.1.subset.2$B

My.Array.1.subset.2  <- data.frame(
  factorB.subset = factorB.subset.2, 
  M.1.subset.2 = M.1.subset.2) # is scales

# Cohort 2 Subset. compare only converters of CIS vs ctrl
id        <- which(clinics.2$Conversion_to_CDMS=='Yes');n.conv <- rownames(clinics.2)[id]
id        <- which(D.2$B==-1);n.ctrl <- rownames(D.2)[id]
n.ctrl.conv <- c(n.ctrl,n.conv)
id        <- which(D.2$A==-1);n.A1   <- rownames(D.2)[id]
n.subset  <- intersect(n.A1,n.ctrl.conv)
D.2.subset   <- D.2[n.subset,];dim(D.2.subset)
M.2.c.subset <- M.2.c[n.subset,]

My.Array.2.c.subset  <- data.frame(
  factorA.subset=D.2.subset$A,  
  factorB.subset=D.2.subset$B, 
  M.2.subset=I(M.2.c.subset))

# only females
####################
id <- which(clinics.1$Sex=='F');    n.1.females <- rownames(clinics.1)[id]
id <- which(clinics.2$Gender=='F'); n.2.females <- rownames(clinics.2)[id]
n.2.a1.females <- intersect(n.2.females,rownames(M.2.a1))
length(n.2.females)
length(n.2.a1.females)

M.1.f       <- M.1[n.1.females,]
M.1.c.f     <- M.1.c[n.1.females,]
D.1.f       <- D.1[n.1.females,]
M.2.c.a1.f  <- M.2.c.a1[n.2.a1.females,]
D.2.a1.f    <- D.2[n.2.a1.females,]

My.Array.1.f  <- data.frame(
  factorA=D.1.f$A,  
  factorB=D.1.f$B, 
  M.1.f=I(M.1.f),
  M.1.c.f=I(M.1.c.f))

# only males
##################
id <- which(clinics.1$Sex=='M');    n.1.males <- rownames(clinics.1)[id]
id <- which(clinics.2$Gender=='M'); n.2.males <- rownames(clinics.2)[id]
n.2.a1.males <- intersect(n.2.males,rownames(M.2.a1))
length(n.2.a1.males)

M.1.m       <- M.1[n.1.males,]
M.1.c.m     <- M.1.c[n.1.males,]
D.1.m       <- D.1[n.1.males,]
M.2.c.a1.m  <- M.2.c.a1[n.2.a1.males,]
D.2.a1.m    <- D.2.a1[n.2.a1.males,]

dim(M.1.c.m);dim(D.1.m)
dim(M.2.c.a1.m);dim(D.2.a1.m)

My.Array.1.m  <- data.frame(
  factorA=D.1.m$A,  
  factorB=D.1.m$B, 
  M.1.f=I(M.1.m),
  M.1.c.f=I(M.1.c.m))

##############################################
# Save and start here
##############################################

# save.image("C:/Users/ellen.fargestad/Documents/R copy/MS/R MS publ.SR/Data/MS data loaded.RData")
# setwd("C:/Users/ellen.fargestad/Documents/R copy/MS/R MS publ.SR")

save.image("H:/A_R/R/R data/R MS publ.SR/Worksheets/MS data loaded.RData")

setwd("H:/A_R/R/R data/R MS publ.SR")
getwd()

library(corrplot)
library(ER)
library(pls)
library(plsVarSel)
library(readxl)
library(corrplot)
library(ffmanova)
library(gprofiler2)
library(corrgram)

load("H:/A_R/R/R data/R MS publ.SR/Worksheets/MS data loaded.RData")


##############################################
# Data Analysis
##############################################

#######################
# Conf Intervall
#######################

# Dataset 1
my.array <- My.Array.1

# across groups
# Compare MS and non-MS patients within cluster 1 and 2
CI.1.B    <- with(my.array, confints(M.1[factorB == -1,], M.1[factorB ==  1,]))


# within each factor
# Compare MS and non-MS patients within cluster 1 and 2
CI.1.B.wicla.A1    <- with(my.array, confints(M.1[factorB == -1 & factorA == -1,],
                                              M.1[factorB ==  1 & factorA == -1,]))
CI.1.B.wicla.A2    <- with(my.array, confints(M.1[factorB == -1 & factorA ==  1,],
                                              M.1[factorB ==  1 & factorA ==  1,]))

# Compare Cluster within MS and within nonMS
CI.1.A.wicla.B1    <- with(my.array, confints(M.1[factorB == -1 & factorA == -1,],
                                              M.1[factorB == -1 & factorA ==  1,]))
CI.1.A.wicla.B2    <- with(my.array, confints(M.1[factorB ==  1 & factorA == -1,],
                                              M.1[factorB ==  1 & factorA ==  1,]))
# Find sign
CI.1.s.A.wicla.B1   <- run.ci(CI=CI.1.A.wicla.B1)
CI.1.s.A.wicla.B2   <- run.ci(CI=CI.1.A.wicla.B2)
CI.1.s.B.wicla.A1   <- run.ci(CI=CI.1.B.wicla.A1)
CI.1.s.B.wicla.A2   <- run.ci(CI=CI.1.B.wicla.A2)

n.1.s.ci.AinB1.pos   <- CI.1.s.A.wicla.B1[[2]];n.1.s.ci.AinB1.neg   <- CI.1.s.A.wicla.B1[[3]]
n.1.s.ci.AinB2.pos   <- CI.1.s.A.wicla.B2[[2]];n.1.s.ci.AinB2.neg   <- CI.1.s.A.wicla.B2[[3]]
n.1.s.ci.BinA1.pos   <- CI.1.s.B.wicla.A1[[2]];n.1.s.ci.BinA1.neg   <- CI.1.s.B.wicla.A1[[3]]
n.1.s.ci.BinA2.pos   <- CI.1.s.B.wicla.A2[[2]];n.1.s.ci.BinA2.neg   <- CI.1.s.B.wicla.A2[[3]]

# sign in both class levels
n.1.s.ci.A.wicla.pos <- intersect(n.1.s.ci.AinB1.pos,n.1.s.ci.AinB2.pos)
n.1.s.ci.A.wicla.neg <- intersect(n.1.s.ci.AinB1.neg,n.1.s.ci.AinB2.neg)
n.1.s.ci.B.wicla.pos <- intersect(n.1.s.ci.BinA1.pos,n.1.s.ci.BinA2.pos)
n.1.s.ci.B.wicla.neg <- intersect(n.1.s.ci.BinA1.neg,n.1.s.ci.BinA2.neg)
n.1.s.ci.A.wicla     <- c(n.1.s.ci.A.wicla.pos,n.1.s.ci.A.wicla.neg)
n.1.s.ci.B.wicla     <- c(n.1.s.ci.B.wicla.pos,n.1.s.ci.B.wicla.neg)

length(n.1.s.ci.A.wicla) # sign for A within both B levels, with consistent pattern
length(n.1.s.ci.B.wicla) # sign for B within both A levels, with consistent pattern
gconvert(n.1.s.ci.B.wicla)

# within each subgr 1 subset: compare MS without non-neurol
# In subgr 1 is Non-neuro vs ctrl: gr=10, MS=gr12
my.array              <- My.Array.1.subset.1
CI.1.B.subset.1         <- with(my.array, confints(M.1.subset.1[factorB.subset.1 == -1,],M.1.subset.1[factorB.subset.1 ==  1,]))
CI.1.s.conf.B.subset.1  <- run.ci(CI=CI.1.B.subset.1)
n.1.s.ci.B.subset.pos.1 <- CI.1.s.conf.B.subset.1[[2]]
n.1.s.ci.B.subset.neg.1 <- CI.1.s.conf.B.subset.1[[3]]
n.1.s.ci.B.subset.neg.and.subset.1 <- intersect(n.1.s.ci.B.wicla.neg,n.1.s.ci.B.subset.neg.1);
length(n.1.s.ci.B.subset.neg.and.subset.1)

# within each subgr 1 subset: compare MS with non-neurol
# In subgr 1 is Non-neuro vs ctrl: gr=10, MS=gr12
my.array              <- My.Array.1.subset.2;str(my.array)
CI.1.B.subset.2         <- with(my.array, confints(M.1.subset.2[factorB.subset.2 == -1,],M.1.subset.2[factorB.subset.2 ==  1,]))
head(CI.1.B.subset.2)
CI.1.s.conf.B.subset.2  <- run.ci(CI=CI.1.B.subset.2)
n.1.s.ci.B.subset.pos.2 <- CI.1.s.conf.B.subset.2[[2]]
n.1.s.ci.B.subset.neg.2 <- CI.1.s.conf.B.subset.2[[3]]
n.1.s.ci.B.subset.neg.and.subset.2 <- intersect(n.1.s.ci.B.wicla.neg,n.1.s.ci.B.subset.neg.2);
length(n.1.s.ci.B.subset.neg.and.subset.2)

# Cohort 2
##############
# among the significant, which are present also in cohort 2
n.1.s.ci.B.wicla.c <- intersect(n.1.s.ci.B.wicla.neg,colnames(M.2.c));
length(n.1.s.ci.B.wicla.c)

my.array <- My.Array.2 ;names(my.array)
CI.2.B.wicla.A1  <- with(my.array,confints(M.2.c[factorB == -1 & factorA == -1,],
                                           M.2.c[factorB ==  1 & factorA == -1,]))
CI.2.s.B.wicla.A1 <- run.ci(CI=CI.2.B.wicla.A1)
n.2.s.ci.B.pos    <- CI.2.s.B.wicla.A1[[2]]
n.2.s.ci.B.neg    <- CI.2.s.B.wicla.A1[[3]]
n.s.ci.B.neg.12   <- intersect(n.1.s.ci.B.wicla.neg,n.2.s.ci.B.neg)# sign among common 
length(n.s.ci.B.neg.12)

my.array                 <- My.Array.2.c.subset ;names(my.array)
CI.2.B.wicla.A1.subset   <- with(my.array,confints(M.2.subset[factorB.subset == -1, ],M.2.subset[factorB.subset ==  1, ]))
CI.2.s.B.wicla.A1.subset <- run.ci(CI=CI.2.B.wicla.A1.subset)
n.2.s.ci.B.neg.subset    <- CI.2.s.B.wicla.A1.subset[[3]]
n.12.s.ci.B.neg.subset   <- intersect(n.1.s.ci.B.wicla.neg,n.2.s.ci.B.neg.subset)# sign among common 
length(n.12.s.ci.B.neg.subset)
gconvert(n.12.s.ci.B.neg.subset)

##################################################################
# ER data set 1,all proteins, to make eED and ED
##################################################################

# ER
my.array <- My.Array.1
y         <- D.1$B;length(y)
er.1      <- ER(M.1 ~ factorA*factorB, data = my.array)
er        <- er.1
ER.set1.A <- as.data.frame(unclass(er$ER.values$factorA))
ER.set1.B <- as.data.frame(unclass(er$ER.values$factorB))
ER.set1.AB<- as.data.frame(unclass(er$ER.values$`factorA:factorB`))
Eff.set1.A    <- as.data.frame(unclass(er$effects$factorA))
Eff.set1.B    <- as.data.frame(unclass(er$effects$factorB))
Eff.set1.AB   <- as.data.frame(unclass(er$effects$`factorA:factorB`))
ER.res        <- as.data.frame(unclass(er$residuals))

# eED
temp.A  <- aggregate(ER.set1.B,by=list(D.1$A.12),mean)[,-1]
temp.B  <- aggregate(ER.set1.B,by=list(D.1$B.12),mean)[,-1]
eED.A   <- temp.A[2,]-temp.A[1,]
eED.B   <- temp.B[2,]-temp.B[1,]
eED.At  <-  t(eED.A)
eED.Bt  <-  t(eED.B)
DF1.ed  <- cbind(eED.At,eED.Bt);colnames(DF1.ed)<- c('eED.A','eED.B')

##################################################################
# ER data set 1c, common
##################################################################

getwd()
save(My.Array.1,file='H:/A_R/R/R data/R MS publ.SR/Data/My.Array.1.Rdata')

# Set1 common 
my.array    <- My.Array.1
y           <- D.1$B;length(y)
er.1c       <- ER(M.1.c ~ factorA*factorB, data = my.array)
er          <- er.1c
ER.set1c.A    <- as.data.frame(unclass(er$ER.values$factorA))
ER.set1c.B    <- as.data.frame(unclass(er$ER.values$factorB))
ER.set1c.AB   <- as.data.frame(unclass(er$ER.values$`factorA:factorB`))
Eff.set1c.A   <- as.data.frame(unclass(er$effects$factorA))
Eff.set1c.B   <- as.data.frame(unclass(er$effects$factorB))
Eff.set1c.AB  <- as.data.frame(unclass(er$effects$`factorA:factorB`))
temp.A      <- aggregate(ER.set1c.A,by=list(D.1$A.12),mean)[,-1]
temp.B      <- aggregate(ER.set1c.B,by=list(D.1$B.12),mean)[,-1]
set1c.eED.A <- temp.A[2,]-temp.A[1,]
set1c.eED.B <- temp.B[2,]-temp.B[1,]
rm(temp.A);rm(temp.B)

# Set2 common, we analyse only patient gr 1 as no CIS are in patient gr 2 
# set 2
temp.B      <- aggregate(M.2.c.a1,by=list(D.2.a1$B.12),mean)[,-1]
set2c.ED.B  <- temp.B[2,]-temp.B[1,]
dim(set1c.eED.B);
dim(set2c.ED.B)

# set 1 and 2
set12c.ED   <- rbind(set1c.eED.B,set2c.ED.B)
set12c.ED.t <- t(set12c.ED);
colnames(set12c.ED.t)<- c('eED.set1c','ED.set2c')
DF1and2.ed <- set12c.ED.t
dim(DF1and2.ed)
DF1and2.ed[1:4,]

#my.ylim <- c(-1.5,1.5)

#####################################################################
# PLSDA within Data set - common proteins
#####################################################################

l     <- 0.05
er    <- er.1c
protein.ass <- colnames(M.1.c)

# set 1 factor A
my.ncomp    <- 6
my.cls.s    <- cls.3[D.1$A.3]
my.pch.s    <- points.3[D.1$B.3]
y           <- D.1$A
pls.mod     <- pls(er, 'factorA', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 1)
scores.1.A    <- pls.mod$scores
loadings.1.A  <- pls.mod$loadings
coeff.1.A.pc1 <- pls.mod$coefficients[,1,1]
residuals.1.A <- pls.mod$residuals[,1,1]
predicted.1.A <- pls.mod$fitted.values[,1,1]

# Jackknifed coefficient P-values
my.ncomp    <- 1
pls.mod     <- pls(er, 'factorA', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 1)
p.jkn.1.A   <- pls.mod$jack[,1,my.ncomp];length(p.jkn.1.A) # select sign only by first PC
id          <- which(p.jkn.1.A<l); length(id)
n.s.1pc    <- protein.ass[id];length(n.s.1pc)
n.s.jkn.1.A.1pc  <- n.s.1pc

my.ncomp    <- 2
pls.mod     <- pls(er, 'factorA', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 1)
p.jkn.1.A   <- pls.mod$jack[,1,my.ncomp];length(p.jkn.1.A) # select sign only by first PC
id          <- which(p.jkn.1.A<l); length(id)
n.s.2pc    <- protein.ass[id];length(n.s.2pc)
n.s.jkn.1.A.2pc  <- n.s.2pc
check.error(pls.mod,ncomp=my.ncomp,y)

# n.s.jkn.1.A.1pc.2pc <- unique(n.s.jkn.1.A.1pc,n.s.jkn.1.A.2pc)
n.s.jkn.1.A.1pc.2pc <- unique(c(n.s.jkn.1.A.1pc,n.s.jkn.1.A.2pc))
length(n.s.jkn.1.A.1pc.2pc)

# set 1 factor B
my.ncomp    <- 5
y           <- D.1$B
pls.mod     <- pls(er, 'factorB', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 1)
scores.1.B    <- pls.mod$scores
loadings.1.B  <- pls.mod$loadings
coeff.1.B.pc1 <- pls.mod$coefficients[,1,1]
coeff.1.B.pc2 <- pls.mod$coefficients[,1,2]
residuals.1.B.1pc <- pls.mod$residuals[,1,1]
predicted.1.B.1pc <- pls.mod$fitted.values[,1,1]
residuals.1.B.2pc <- pls.mod$residuals[,1,2]
predicted.1.B.2pc <- pls.mod$fitted.values[,1,2]
check.error(pls.mod,ncomp=my.ncomp,y)

# Jackknifed coefficient P-values (sorted)
my.ncomp    <- 4
pls.mod     <- pls(er, 'factorB', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 1)
p.jkn.1.B   <- pls.mod$jack[,1,my.ncomp];length(p.jkn.1.B) # select sign only by first PC
id          <- which(p.jkn.1.B<l); length(id)
n.s.jkn.1.B.4pc  <- protein.ass[id]

my.ncomp    <- 3
pls.mod     <- pls(er, 'factorB', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 1)
p.jkn.1.B   <- pls.mod$jack[,1,my.ncomp];length(p.jkn.1.B) # select sign only by first PC
id          <- which(p.jkn.1.B<l); length(id)
n.s.jkn.1.B.3pc  <- protein.ass[id]

my.ncomp    <- 2
pls.mod     <- pls(er, 'factorB', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 1)
p.jkn.1.B   <- pls.mod$jack[,1,my.ncomp];length(p.jkn.1.B) # select sign only by first PC
id          <- which(p.jkn.1.B<l); length(id)
n.s.jkn.1.B.2pc  <- protein.ass[id]

my.ncomp    <- 1
pls.mod     <- pls(er, 'factorB', ncomp= my.ncomp, validation='LOO',jackknife = TRUE, df.used = 1)
p.jkn.1.B   <- pls.mod$jack[,1,my.ncomp];length(p.jkn.1.B) # select sign only by first PC
id          <- which(p.jkn.1.B<l); length(id)
n.s.jkn.1.B.1pc  <- protein.ass[id]

n.s.jkn.1.B.1234pc <- unique(c(n.s.jkn.1.B.1pc,n.s.jkn.1.B.2pc,n.s.jkn.1.B.3pc,n.s.jkn.1.B.4pc))

# set2
###########
my.ncomp  <- 2
my.cls.s  <- cls.2[D.2.a1$A.12]
my.pch.s  <- points.2[D.2.a1$B.12]
y         <- as.numeric(D.2.a1$B)
x         <- as.matrix(M.2.a1);dim(x)
my.data   <- data.frame(y = x, x = x)
pls.mod   <- plsr(y~x, data = my.data, jackknife=TRUE, ncomp = my.ncomp, validation = "LOO")
check.error(pls.mod,ncomp=my.ncomp,y)
scores.2.B    <- pls.mod$scores
loadings.2.B  <- pls.mod$loadings
coeff.2.B.pc1 <- pls.mod$coefficients[,1,1]
coeff.2.B.pc2 <- pls.mod$coefficients[,1,2]
residuals.2.B.1pc <- pls.mod$residuals[,1,1]
predicted.2.B.1pc <- pls.mod$fitted.values[,1,1]
residuals.2.B.2pc <- pls.mod$residuals[,1,2]
predicted.2.B.2pc <- pls.mod$fitted.values[,1,2]

my.ncomp  <- 2
pls.mod   <- plsr(y~x, data = my.data, jackknife=TRUE, ncomp = my.ncomp, validation = "LOO")
jkn       <- jack.test(pls.mod, ncomp = my.ncomp, use.mean = TRUE)
p.jkn     <- jkn$pvalues
id        <- which(p.jkn<l);length(id)
n.s.jkn.2.B.2pc <- colnames(x)[id]

my.ncomp  <- 1
pls.mod   <- plsr(y~x, data = my.data, jackknife=TRUE, ncomp = my.ncomp, validation = "LOO")
jkn       <- jack.test(pls.mod, ncomp = my.ncomp, use.mean = TRUE)
p.jkn     <- jkn$pvalues
id        <- which(p.jkn<l);length(id)
n.s.jkn.2.B.1pc <- colnames(x)[id]

n.s.jkn.2.B.1pc.2pc <- unique(c(n.s.jkn.2.B.1pc,n.s.jkn.2.B.2pc))

# select proteins for factor A
###############################
N.s.A.jkn.set1    <- n.s.jkn.1.A.1pc;length(N.s.A.jkn.set1)

# pos/neg for coeff 
id <- which(coeff.1.A.pc1>0); n.coeff.1.pos <- names(coeff.1.A.pc1)[id] 
id <- which(coeff.1.A.pc1<0); n.coeff.1.neg <- names(coeff.1.A.pc1)[id] 

# pos for coeff at 2pc in both sets
N.s.A.jkn.set1.pos <- intersect(N.s.A.jkn.set1,n.coeff.1.pos)
N.s.A.jkn.set1.neg <- intersect(N.s.A.jkn.set1,n.coeff.1.neg)

length(N.s.A.jkn.set1.pos)
length(N.s.A.jkn.set1.neg)



# select proteins for factor B - use up to optimal for each cohort
####################################################################

# include pc up to optimal in both studies, use this
N.s.set1    <- unique(c(n.s.jkn.1.B.1pc,n.s.jkn.1.B.2pc,n.s.jkn.1.B.3pc,n.s.jkn.1.B.4pc));length(N.s.set1)
N.s.set2    <- unique(c(n.s.jkn.2.B.1pc,n.s.jkn.2.B.2pc));length(N.s.set2)
N.s.B.jkn.set1and2  <- intersect(N.s.set1,N.s.set2);length(N.s.B.jkn.set1and2)
#N.sel2     <- N.s.B.jkn.set1and2

# neg for coeff at 2pc in both sets
id <- which(coeff.1.B.pc2<0); n.coeff.1.neg <- names(coeff.1.B.pc2)[id] 
id <- which(coeff.2.B.pc2<0); n.coeff.2.neg <- names(coeff.2.B.pc2)[id] 
N.coeff.1and2.neg <- intersect(n.coeff.1.neg,n.coeff.2.neg)

# pos for coeff at 2pc in both sets
id <- which(coeff.1.B.pc2>0); n.coeff.1.pos <- names(coeff.1.B.pc2)[id] 
id <- which(coeff.2.B.pc2>0); n.coeff.2.pos <- names(coeff.2.B.pc2)[id] 
N.coeff.1and2.pos <- intersect(n.coeff.1.pos,n.coeff.2.pos)

N.s.B.jkn.set1and2.pos <- intersect(N.s.B.jkn.set1and2,N.coeff.1and2.pos)
N.s.B.jkn.set1and2.neg <- intersect(N.s.B.jkn.set1and2,N.coeff.1and2.neg)
N.s.B.jkn.set1and2.consistent <- unique(c(N.s.B.jkn.set1and2.pos,N.s.B.jkn.set1and2.neg))
length(N.s.B.jkn.set1and2.consistent) # selected in both sets, consistent pattern

sort(PG.id[N.s.B.jkn.set1and2.consistent,3])

##################################################################
# Combine the data, include all proteins 
##################################################################

# add gr4 to D.2.a1
gr4       <- D.2.a1$B.12 # all are a1 so gr4 and B.12 is the same for D.2
D.2.a1    <- cbind(D.2.a1,gr4);dim(D.2.a1)
common.D  <- intersect(colnames(D.1),colnames(D.2.a1))

ER.set1c.B  <- as.matrix(ER.set1c.B)
N.sel       <- intersect(colnames(ER.set1c.B),colnames(M.2.a1));length(N.sel)
M.12 <- as.matrix(rbind(ER.set1c.B[,N.sel],M.2.a1[,N.sel]));dim(M.12)
D.12 <- rbind(D.1[,common.D],D.2.a1[,common.D])   ;dim(D.12)
id1  <- which(substr(rownames(D.12),1,1)=='n')
id2  <- which(substr(rownames(D.12),1,1)=='M')
id3  <- which(substr(rownames(D.12),3,3)=='S')
id4  <- which(substr(rownames(D.12),3,3)=='R')
Id.set1 <- sort(c(id1,id2))
Id.set2 <- sort(c(id3,id4))
C             <-  rep(NA,each=dim(D.12)[1]);names(C)<- rownames(D.12) 
C.12          <-  C
C[Id.set1]    <- -1
C[Id.set2]    <-  1
C.12[Id.set1] <-  1
C.12[Id.set2] <-  2
Gender.1      <- matrix(NA,nrow=dim(D.1)[1],ncol=2);dim(D.1);dim(Gender.1)
Gender.2      <- matrix(NA,nrow=dim(D.2.a1)[1],ncol=2);dim(D.2);dim(Gender.2)
rownames(Gender.1)            <- rownames(D.1);   rownames(Gender.2)<- rownames(D.2.a1)
colnames(Gender.1)            <- c('G','G.extra');colnames(Gender.2)<- c('G','G.extra')
Gender.1[n.1.females,1:2]     <- 'f';Gender.1[n.1.males,1:2]<- 'm'
Gender.2[n.2.a1.females,1:2]  <- 'f';Gender.2[n.2.a1.males,1:2]<- 'm'
Gender  <- rbind(Gender.1,Gender.2);dim(Gender)
G       <- Gender[,1]
DCG.12          <- data.frame(D.12,C,C.12,G);D.12[1:3,]
factorB       <- as.factor(D.12$B)
factorC       <- as.factor(unname(C))
dim(D.12)

My.Array.12 <- data.frame(D.12 = I(DCG.12), 
                          factorB = factorB, 
                          factorC = factorC,
                          G = G,
                          M.12 = I(M.12))

N.sel <- c(N.s.B.jkn.set1and2.consistent,n.s.ci.B.neg.12);length(N.sel)
Inn          <- M.12[,N.sel]

N.inn1     <- n.s.ci.B.neg.12;length(N.inn1)
N.inn2     <- setdiff(colnames(Inn),n.s.ci.B.neg.12)

Inn1  <- Inn[,N.inn1] # sign by CI
Inn2  <- Inn[,N.inn2] # additional proteins sign by Jkn

my.array <- My.Array.12
er.12      <- ER(M.12 ~ factorC*factorB, data = my.array)
er        <- er.12
ER.set12.C <- as.data.frame(unclass(er$ER.values$factorC))
ER.set12.B <- as.data.frame(unclass(er$ER.values$factorB))
ER.set12.CB<- as.data.frame(unclass(er$ER.values$`factorC:factorB`))

rm(Gender.1);rm(Gender.2);rm(factorB);rm(factorC)

#######################################
## combine all data a1 and a2
#######################################

common.D.all <- intersect(colnames(D.1),colnames(D.2));common.D.all
D.12.all <- rbind(D.1[,common.D.all],D.2[,common.D.all])   ;dim(D.12.all)
id1  <- which(substr(rownames(D.12.all),1,1)=='n')
id2  <- which(substr(rownames(D.12.all),1,1)=='M')
id3  <- which(substr(rownames(D.12.all),3,3)=='S')
id4  <- which(substr(rownames(D.12.all),3,3)=='R')
Id.set1 <- sort(c(id1,id2))
Id.set2 <- sort(c(id3,id4))
C             <-  rep(NA,each=dim(D.12.all)[1]);names(C)<- rownames(D.12.all) 
C.12          <-  C
C[Id.set1]    <- -1
C[Id.set2]    <-  1
C.12[Id.set1] <-  1
C.12[Id.set2] <-  2
Gender.1      <- matrix(NA,nrow=dim(D.1)[1],ncol=2);dim(D.1);dim(Gender.1)
Gender.2      <- matrix(NA,nrow=dim(D.2)[1],ncol=2);dim(D.2);dim(Gender.2)
rownames(Gender.1)            <- rownames(D.1);   rownames(Gender.2)<- rownames(D.2)
colnames(Gender.1)            <- c('G','G.extra');colnames(Gender.2)<- c('G','G.extra')
Gender.1[n.1.females,1:2]     <- 'f';Gender.1[n.1.males,1:2]<- 'm'
Gender.2[n.2.females,1:2]     <- 'f';Gender.2[n.2.males,1:2]<- 'm'
Gender  <- rbind(Gender.1,Gender.2);dim(Gender)
G       <- Gender[,1]
DCG.12.all          <- data.frame(D.12.all,C,C.12,G);DCG.12.all[1:3,]
factorB       <- as.factor(D.12$B)
factorC       <- as.factor(unname(C))

# A1
D             <- DCG.12.all
id            <- which(D$A.12==1);length(id)
D.12.A1       <- D[id,]
Temp          <- rbind(M.1.c,M.2);dim(Temp)
M.12.A1       <- Temp[id,]
factorB       <- D.12.A1$B
factorC       <- D.12.A1$C

My.Array.12.A1 <- data.frame(D.12 = I(D.12.A1), 
                             factorB = factorB, 
                             factorC = factorC,
                             M.12.A1 = I(M.12.A1))

# A2
id <- which(D$A.12==2);length(id)
D.12.A2       <- D[id,]
Temp          <- rbind(M.1.c,M.2);dim(Temp)
M.12.A2       <- Temp[id,]
factorB       <- D.12.A2$B
factorC       <- D.12.A2$C


My.Array.12.A2 <- data.frame(D.12 = I(D.12.A2), 
                          factorB = factorB, 
                          factorC = factorC,
                          M.12.A2 = I(M.12.A2))

#######################################
## CI across the data
#######################################

# across ER.set 1 combined with set 2.a1
my.array      <- My.Array.12;names(my.array)
CI.12.B       <- with(my.array,confints(M.12[factorB == -1, ],M.12[factorB ==  1, ]))
CI.12.B.s     <- run.ci(CI=CI.12.B)
n.s.CI.12.B   <- CI.12.B.s[[3]]
length(n.s.CI.12.B)


# across both A1
my.array      <- My.Array.12.A1;names(my.array)
CI.12.B.A1    <- with(my.array,confints(M.12.A1[factorB == -1, ],M.12.A1[factorB ==  1, ]))
CI.12.B.s.A1  <- run.ci(CI=CI.12.B.A1)
n.s.CI.12.B.A1   <- CI.12.B.s.A1[[3]]
length(n.s.CI.12.B.A1)


# across both A2
my.array      <- My.Array.12.A2;names(my.array)
CI.12.B.A2    <- with(my.array,confints(M.12.A2[factorB == -1, ],M.12.A2[factorB ==  1, ]))
CI.12.B.s.A2  <- run.ci(CI=CI.12.B.A2)
n.s.CI.12.B.A2   <- CI.12.B.s.A2[[3]]
length(n.s.CI.12.B.A2)


n.s.CI.12.within.A1.A2 <- intersect(n.s.CI.12.B.A1,n.s.CI.12.B.A2);length(n.s.CI.12.within.A1.A2)
sort(PG.id[n.s.CI.12.within.A1.A2,3])

test <- intersect(n.s.CI.12.within.A1.A2,n.1.s.ci.B.wicla);length(test)
sort(PG.id[test,3])

test <- intersect(n.s.CI.12.within.A1.A2,N.s.B.jkn.set1and2.consistent);length(test)
sort(PG.id[test,3])



#######################################
## ffmanova across the data
#######################################

id <- which(D.12$C.12==1);n.set1 <- rownames(D.12)[id]
id <- which(D.12$C.12==2);n.set2 <- rownames(D.12)[id]

# data set 1 and 2 combined
###############################

l       <- 0.05
R       <- M.12 ;dim(R)
B       <- DCG.12$B
C       <- DCG.12$C
my.nSim <- 9999
my.data     <- data.frame(B,C,R);dim(my.data)
#FF.res.12   <- ffmanova(R ~ B+C, data = my.data,stand=TRUE, nSim = my.nSim, verbose=TRUE)
FF.res.12   <- ffmanova(R ~ B*C, data = my.data,stand=TRUE, nSim = my.nSim, verbose=TRUE)
FF.res.12
p.Raw       <- FF.res$pRaw[-1,]
p.FDR       <- FF.res$pAdjFDR[-1,]
p.values.12 <- rbind(p.Raw,p.FDR);dim(p.values.12)
FF.res      <- FF.res.12
p.FDR       <- FF.res$pAdjFDR[-1,]
rownames(p.values.12)<- c('p.Raw.B','p.Raw.C','p.Raw.BC','p.FDR.B','p.FDR.C','p.FDR.BC')

n.B.sign    <- names(which(p.FDR[1,]<l));length(n.B.sign);#sort(PG.id[n.B.sign,3])
n.BC.sign   <- names(which(p.FDR[3,]<l));length(n.BC.sign);#sort(PG.id[n.BC.sign,3])

# sign by different levels
nB.0.05     <- names(which(p.FDR[1,n.B.sign]<0.05))
nB.0.01     <- names(which(p.FDR[1,n.B.sign]<0.01))
nB.0.001    <- names(which(p.FDR[1,n.B.sign]<0.001))
n.FDR.12.B.sign <- n.B.sign

nB.only.0.001 <- nB.0.001
nB.only.0.01 <- setdiff(nB.0.01,nB.0.001)
nB.only.0.05 <- setdiff(nB.0.05,nB.0.01)
length(nB.only.0.001)+length(nB.only.0.01)+length(nB.only.0.05)

# combined with Jkn
nB.0.05.jkn <- intersect(nB.only.0.05,N.s.B.jkn.set1and2.consistent);length(nB.0.05.jkn)
nB.0.01.jkn <- intersect(nB.only.0.01,N.s.B.jkn.set1and2.consistent);length(nB.0.01.jkn)
nB.0.001.jkn <- intersect(nB.only.0.001,N.s.B.jkn.set1and2.consistent);length(nB.0.001.jkn)
length(nB.0.05.jkn)+length(nB.0.01.jkn)+length(nB.0.001.jkn)

# sort(PG.id[nB.0.05.jkn,3])
# sort(PG.id[nB.0.01.jkn,3])
# sort(PG.id[nB.0.001.jkn,3])
# PG.id[nB.0.05,]# Q92876, KLK6
# PG.id['Q92876',]
# PG.id['Q8N3J6',]
# PG.id['Q15582',]

length(n.FDR.12.B.sign)
length(n.s.CI.12.B)
n.CIFDR.B.sign <- intersect(n.FDR.12.B.sign,n.s.CI.12.B)
length(n.CIFDR.B.sign)

# sign by pls consistent
length(N.s.B.jkn.set1and2.consistent)

test <- intersect(N.s.B.jkn.set1and2.consistent,n.s.CI.12.B);length(test)
test <- intersect(N.s.B.jkn.set1and2.consistent,n.FDR.B.sign);length(test)
test <- intersect(N.s.B.jkn.set1and2.consistent,n.CIFDR.B.sign);length(test)


########################
## Compare with Opsahl
########################

# Sign by Opsahl.FigS1.and by ANOVA.FDR of M.12
#n.test <- intersect(nB.0.05.jkn,n.Opsahl.Fig.S1);length(n.test)

# Sign by Opsahl.TablS2.and by ANOVA.FDR of M.12
#n.test <- intersect(nB.0.05.jkn,n.Opsahl.Tab.S2);length(test)

# Sign by Opsahl.TablS2.and by PLS common
# n.test <- intersect(N.s.B.jkn.set1and2.consistent,n.Opsahl.Tab.S2);length(n.test)
# PG.id[n.test,]

# Sign by Opsahl.TablS2.and by PLS common
# n.test <- intersect(n.s.jkn.1.B.1234pc,n.Opsahl.Tab.S2);length(n.test)
# PG.id[n.test,]


####################################
# PCA Karim - Data set 1
####################################

my.cls.s    <- cls.3[D.1$A.3];# ctrl.non.neur - ctrl.other.neuro - MS
my.pch.s    <- points.3[D.1$B.3]

# Karims list of 9 proteins Inn MS
length(colnames(M.1))
Inn       <- M.1[,n.K.ms9];dim(Inn);Inn[1:3,]

# males
TEMP1     <- t(Inn)
TEMP2     <- data.frame(PG.id[rownames(TEMP1),],TEMP1);dim(TEMP2)
o         <- order(TEMP2$name)
TEMP3     <- TEMP2[o,];
rownames(TEMP3)<- TEMP3$name
idc1      <- which(colnames(TEMP3)=='nonMS1')
TEMP4     <- TEMP3[,idc1:dim(TEMP3)[2]]
Inn       <- t(TEMP4)

pca.mod   <- prcomp(Inn,scale=FALSE) # all data are already scales
scores    <- pca.mod$x
loadings  <- pca.mod$rotation
expl.var  <- round(propVar(pca.mod)*100,digits = 0);expl.var

##################################################
# Fig 1. PCA 9 proteins
##################################################


pdf('./Figures/Fig 1 PCA of 9 proteins selected for MS versus controls.pdf')

par(mfrow=c(2,2))
plot(scores[,1:2],col=my.cls.s,pch=my.pch.s,
     xlab=paste0('PC1 ',expl.var[1],' %'),
     ylab=paste0('PC2 ',expl.var[2],' %'))
abline(h=0,v=0,lty=1,col='gray50')
title("(a)", adj = 0,line = 1)

my.xlim <- c(-0.25,0.5)
plot(loadings[,1:2],col='white',cex=0.00001,
     xlim=my.xlim,
     xlab=paste0('PC1 ',expl.var[1],' %'),
     ylab=paste0('PC2 ',expl.var[2],' %'))
text(loadings[,1:2],labels=rownames(loadings),cex=0.8,
     xlim=my.xlim)
abline(h=0,v=0,lty=1,col='gray50')
title("(b)", adj = 0,line = 1)

dev.off()

###############################################
# Fig S1 1D plot of proteins selected by Karim
###############################################

pdf('./Figures/Fig S1 1D plot of proteins. MS vs controls.pdf')

reorder.col <- c("IGHG1","IGHV4-34","OMD","JCHAIN","CHI3L2","ORM1","LYVE1","HPR","IGSF21")   
Inn <- Inn[,reorder.col]
Means.inn <- aggregate(Inn,by=list(D.1$gr4),mean)[D.1$gr4,-1]

par(mfrow=c(3,3))
for (i in 1:dim(Inn)[2]){
  plot(Inn[,i],col=my.cls.s,pch=my.pch.s,
       cex.axis =1.2,cex.lab=1.4,
       xlab='participants',ylab='z-scores',main=colnames(Inn)[i])
  abline(v=63,lty=2,col='gray50')
  abline(h=0,lty=1,col='gray50')
  lines(Means.inn[,i])
}

dev.off()

#########################
# Fig 2 IgGs Barplots 
#########################

#N.sel.A <- n.1.s.ci.A.wicla;length(N.sel.A)
#N.sel.B <- n.1.s.ci.B.wicla;length(N.sel.B)

pdf('./figures/Fig 2 Cohort 1. Bar plot Ig.sign P.gr.pdf')

my.cex.text <- 1.2
my.ylim <-  c(-1.6,1.6)

# A
######
# IgG
id         <- which(substr(PG.id$names.short,1,3)=='Ig ')
n.Ig       <- rownames(PG.id)[id]
N.sel      <- intersect(n.Ig,n.1.s.ci.A.wicla);length(N.sel)
gr.names.4 <- vector()
Inn.M      <-  M.1[,N.sel];dim(Inn.M)
Inn.D     <- D.1
my.cls.ci <- c('blue','red')[c(rep(1,each=110),rep(2,each=120))]
id <- which(D.1$gr4==1);gr.names.4[id] <- '1.A1.B1'
id <- which(D.1$gr4==2);gr.names.4[id] <- '2.A1.B2'
id <- which(D.1$gr4==3);gr.names.4[id] <- '3.A2.B1'
id <- which(D.1$gr4==4);gr.names.4[id] <- '4.A2.B2'
Means.sel   <- aggregate(Inn.M,by=list(gr.names.4),mean)[,-1] 
D.gr        <- aggregate(Inn.D,by=list(gr.names.4),mean)
Means.sel.t <- t(as.matrix(Means.sel));dim(Means.sel.t)
colnames(Means.sel.t)<- c('A1.B1.N','A1.B2','A2.B1.N','A2.B2')
Inn  <- Means.sel.t;dim(Inn)
Means.data.Ig  <- Inn
barplot(c(Inn[,1],c(0,0,0,0,0,0),Inn[,2],c(0,0,0,0,0,0,0,0),Inn[,3],
          c(0,0,0,0,0,0),Inn[,4]),
        col=my.cls.ci,xaxt='n',border=NA, 
        ylim = my.ylim,ylab='z-scores')
abline(v=135,lty=2,col='gray50')
text(30,-1.25,labels='1 C',   col='blue',cex=my.cex.text)
text(90,-1.25,labels='1 MS',    col='blue',cex=my.cex.text)
text(180,-1.25,labels='2 C',  col='red',cex=my.cex.text)
text(240,-1.25,labels='2 MS',   col='red',cex=my.cex.text)

dev.off()



##################################################
# Fig 3 Bar plot of selected for MS by Jkn 5 groups
##################################################

pdf('./figures/Fig 3 Cohort 1. Bar plot. Protein sign for MS vs ctrl.pdf')
par(mfrow=c(2,2))
N.sel       <- n.1.s.ci.B.wicla.neg
my.cls.ci   <- c('blue','red')[c(rep(1,each=60),rep(2,each=60))]
my.ylim <-  c(-1.6,1.6)

# Data and ER
N.sel       <- n.1.s.ci.B.wicla.neg
my.cls.ci   <- c('blue','red')[c(rep(1,each=60),rep(2,each=60))]
my.ylim <-  c(-1.6,1.6)

Inn.M       <-  M.1[,N.sel];dim(Inn.M)
Means.sel   <- aggregate(Inn.M,by=list(gr.names.4),mean)[,-1] 
D.gr        <- aggregate(Inn.D,by=list(gr.names.4),mean)
Means.sel.t <- t(as.matrix(Means.sel));dim(Means.sel.t)
colnames(Means.sel.t)<- c('A1.B1.N','A1.B2','A2.B1.N','A2.B2')
Inn         <- Means.sel.t
Means.data.sign.MS  <- Inn
barplot(c(Inn[,1],c(0,0,0,0,0,0),Inn[,2],c(0,0,0,0,0,0,0,0),Inn[,3],
          c(0,0,0,0,0,0),Inn[,4]),
        col=my.cls.ci,xaxt='n',border=NA, 
        ylim = my.ylim,ylab='z-scores')
abline(v=70,lty=2,col='gray50')
text(15,-1.25,labels='1 C',   col='blue',cex=my.cex.text)
text(50,-1.25,labels='1 MS',    col='blue',cex=my.cex.text)
text(90,-1.25,labels='2 C',  col='red',cex=my.cex.text)
text(125,-1.25,labels='2 MS',   col='red',cex=my.cex.text)
title("(a)", adj = 0,line = 1)

# Er.A
Inn.ER.A    <-  ER.set1.A[,N.sel];dim(Inn)
Inn.ER.A.sel<-  Inn.ER.A[,N.sel];dim(Inn.M)
Means.sel   <- aggregate(Inn.ER.A.sel,by=list(gr.names.4),mean)[,-1] 
D.gr        <- aggregate(Inn.D,by=list(gr.names.4),mean)
Means.sel.t <- t(as.matrix(Means.sel));dim(Means.sel.t)
colnames(Means.sel.t)<- c('A1.B1.N','A1.B2','A2.B1.N','A2.B2')
Inn         <- Means.sel.t
Means.ER.A.sign.MS  <- Inn
barplot(c(Inn[,1],c(0,0,0,0,0,0),Inn[,2],c(0,0,0,0,0,0,0,0),Inn[,3],
          c(0,0,0,0,0,0),Inn[,4]),
        col=my.cls.ci,xaxt='n',border=NA, 
        ylim = my.ylim,ylab='z-scores')
abline(v=70,lty=2,col='gray50')
text(15,-1.25,labels='1 C',   col='blue',cex=my.cex.text)
text(50,-1.25,labels='1 MS',    col='blue',cex=my.cex.text)
text(90,-1.25,labels='2 C',  col='red',cex=my.cex.text)
text(125,-1.25,labels='2 MS',   col='red',cex=my.cex.text)
title("(b)", adj = 0,line = 1)

# ER.B
Inn.ER.B    <-  ER.set1.B[,N.sel];dim(Inn)
Inn.ER.B.sel<-  Inn.ER.B[,N.sel];dim(Inn.M)
Means.sel   <- aggregate(Inn.ER.B.sel,by=list(gr.names.4),mean)[,-1] 
D.gr        <- aggregate(Inn.D,by=list(gr.names.4),mean)
Means.sel.t <- t(as.matrix(Means.sel));dim(Means.sel.t)
colnames(Means.sel.t)<- c('A1.B1.N','A1.B2','A2.B1.N','A2.B2')
Inn         <- Means.sel.t
Means.ER.B.sign.MS  <- Inn
barplot(c(Inn[,1],c(0,0,0,0,0,0),Inn[,2],c(0,0,0,0,0,0,0,0),Inn[,3],
          c(0,0,0,0,0,0),Inn[,4]),
        col=my.cls.ci,xaxt='n',border=NA, 
        ylim = my.ylim,ylab='z-scores')
abline(v=70,lty=2,col='gray50')
text(15,-1.25,labels='1 C',   col='blue',cex=my.cex.text)
text(50,-1.25,labels='1 MS',    col='blue',cex=my.cex.text)
text(90,-1.25,labels='2 C',  col='red',cex=my.cex.text)
text(125,-1.25,labels='2 MS',   col='red',cex=my.cex.text)
title("(c)", adj = 0,line = 1)

dev.off()

#######################################
## Fig 4 PLS plots, common proteins
#######################################

pdf('./Figures/Fig 4 PLS-DA  A and B cohort 1 and B cohort 2.pdf')

par(mfrow=c(3,2))

# scores 

# set 1 factor A
head(D.1)
D.1$A.3
my.cls.s    <- cls.3[D.1$A.3];# ctrl.non.neur - ctrl.other.neuro - MS
my.pch.s    <- points.3[D.1$B.3]
#my.cls.s    <- cls[D.1$A.12];#my.cls.s
#my.pch.s    <- points[D.1$B.12]
plot(scores.1.A[,1:2],col=my.cls.s,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title("(a)", adj = 0,line = 1)

v           <- rep(1,each=dim(loadings.1.A)[1]);length(v)
names(v)    <- rownames(loadings.1.A)
v[n.s.jkn.1.A.1pc.2pc] <- 2
my.cls.v    <- c('gray65','black')[v]
my.pch.v    <- c(1,15,17,19)[v]

plot(loadings.1.A[,1:2],col=my.cls.v,pch=my.pch.v,
     xlim=c(-0.15,0.15),ylim=c(-0.15,0.15),
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title("(b)", adj = 0,line = 1)

#par(mfrow=c(2,2))
# set 1 factor B
my.cls.s    <- cls.3[D.1$A.3];# ctrl.non.neur - ctrl.other.neuro - MS
my.pch.s    <- points.3[D.1$B.3]
#my.cls.s    <- cls[D.1$A.12]
#my.pch.s    <- points[D.1$B.12]
plot(scores.1.B[,1:2],col=my.cls.s,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title("(c)", adj = 0,line = 1)

v           <- rep(1,each=dim(loadings.1.B)[1]);length(v)
names(v)    <- rownames(loadings.1.A)
#v[N.s.set1] <- 2
v[N.s.B.jkn.set1and2.consistent]<- 3
my.cls.v    <- c('gray65','gray50','black')[v]
my.pch.v    <- c(1,19,15)[v]

plot(loadings.1.B[,1:2],col=my.cls.v,pch=my.pch.v,
     xlim=c(-0.15,0.15),ylim=c(-0.15,0.15),
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title("(d)", adj = 0,line = 1)

my.cls.s    <- cls.2[D.2$A.12]
my.pch.s    <- points.2[D.2$B.12]
plot(scores.2.B[,1:2],col=my.cls.s,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
#     labels=Names.set2)
abline(v=0,h=0,lty=2,col='gray50')
title("(e)", adj = 0,line = 1)

v           <- rep(1,each=dim(loadings.1.B)[1]);length(v)
names(v)    <- rownames(loadings.1.A)
#v[N.s.set2] <- 2
v[N.s.B.jkn.set1and2.consistent]<- 3
my.cls.v    <- c('gray65','gray50','black')[v]
my.pch.v    <- c(1,19,15)[v]

plot(loadings.2.B[,1:2],col=my.cls.v,pch=my.pch.v,
     xlim=c(-0.11,0.11),ylim=c(-0.23,0.23),
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title("(f)", adj = 0,line = 1)

dev.off()


###########################################
# Fig S2 Cohort 1. CI with protein ASS.id
###########################################

pdf('./Figures/Fig S2 Cohort 1 Plot of CI. sign B.pdf')

my.ylim <- c(-2,2)
my.cex <- 0.6
# Cohort 1. CI for B within A1
par(mfrow=c(2,2))
N.sel     <- n.1.s.ci.B.wicla.neg;length(N.sel)
CI        <- CI.1.B.wicla.A1[N.sel,];
TEMP      <- data.frame(CI,PG.id[N.sel,])
o         <- order(TEMP$name)
CI.names  <- TEMP[o,];CI.names[1:3,]
rownames(CI.names)<- CI.names$Gene.names
Inn <- CI.names;colnames(Inn)
plot(Inn$A,xaxt='n',
     xlab='',ylab='z-scores',
     ylim=my.ylim)
abline(h=0,lty=2,col='gray20')
lines(Inn$A,ylim=my.ylim)
lines(Inn$Left,ylim=my.ylim,col='gray50')
lines(Inn$Right,ylim=my.ylim,col='gray50')
axis(1,at=1:length(CI$A),labels = Inn$name,las=2,cex.axis=my.cex)
title("(a)", adj = 0,line = 1)

# Cohort 1. CI for B within A2
N.sel     <- n.1.s.ci.B.wicla.neg;length(N.sel)
CI        <- CI.1.B.wicla.A2[N.sel,];CI[1:3,]
TEMP      <- data.frame(CI,PG.id[N.sel,]);TEMP[1:3,]
o         <- order(TEMP$name)
CI.names  <- TEMP[o,];CI.names[1:3,]
rownames(CI.names)<- CI.names$Gene.names
Inn <- CI.names;colnames(Inn)
plot(Inn$A,xaxt='n',
     xlab='',ylab='z-scores',
     ylim=my.ylim)
abline(h=0,lty=2,col='gray20')
lines(Inn$A,ylim=my.ylim)
lines(Inn$Left,ylim=my.ylim,col='gray50')
lines(Inn$Right,ylim=my.ylim,col='gray50')
axis(1,at=1:length(CI$A),labels = Inn$name,las=2,cex.axis=my.cex)
title("(b)", adj = 0,line = 1)

N.sel     <- n.1.s.ci.B.wicla.neg;length(N.sel)
CI        <- CI.1.B.subset[N.sel,];CI[1:3,]
TEMP      <- data.frame(CI,PG.id[N.sel,]);TEMP[1:3,]
o         <- order(TEMP$name)
CI.names  <- TEMP[o,];CI.names[1:3,]
rownames(CI.names)<- CI.names$Gene.names
Inn <- CI.names;colnames(Inn)
plot(Inn$A,xaxt='n',
     xlab='',ylab='z-scores',
     ylim=my.ylim)
abline(h=0,lty=2,col='gray20')
lines(Inn$A,ylim=my.ylim)
lines(Inn$Left,ylim=my.ylim,col='gray50')
lines(Inn$Right,ylim=my.ylim,col='gray50')
axis(1,at=1:length(CI$A),labels = Inn$name,las=2,cex.axis=my.cex)
title("(c)", adj = 0,line = 1)

dev.off()


##################################################################
# Fig S4. Dataset 2 P01834 
##################################################################
M.2 <- M.2.c

pdf('./Figures/Fig S4 Cohort 2 P01834 IGKC.pdf')
par(mfrow=c(2,2))
my.ylim <- c(-3,4)
my.cls.s  <- cls.2[D.2.a1$A.12]
my.pch.s  <- points.2[D.2.a1$B.12]
plot(M.2[id2.a1,'P01834'],col=my.cls.s,
     xlab='Participants',ylab='z-scores',
     pch=my.pch.s,ylim=my.ylim)
abline(h=0,lty=1,col='gray50')
title("(a)", adj = 0,line = 1)

my.cls.s  <- cls.2[D.2.a2$A.12]
my.pch.s  <- points.2[D.2.a2$B.12]
plot(M.2[id2.a2,'P01834'],col=my.cls.s,
     xlab='Participants',ylab='z-scores',
     pch=my.pch.s,ylim=my.ylim)
abline(h=0,lty=1,col='gray50')
title("(b)", adj = 0,line = 1)
dev.off()

###########################################
# Fig S5 Cohort 2. CI with protein ASS.id
###########################################

# Cohort 2. CI for B within A1
pdf('./Figures/Fig S5 Cohort 2 Plot of CI. sign B.pdf')
par(mfrow=c(2,2))
N.sel     <- n.s.ci.B.neg.12;length(N.sel)
CI        <- CI.2.B.wicla.A1[N.sel,];
TEMP      <- data.frame(CI,PG.id[N.sel,])
o         <- order(TEMP$name)
CI.names  <- TEMP[o,];CI.names[1:3,]
rownames(CI.names)<- CI.names$Gene.names
Inn       <- CI.names;colnames(Inn)
plot(Inn$A,xaxt='n',
     xlab='',ylab='z-scores',
     ylim=my.ylim)
abline(h=0,lty=2,col='gray20')
lines(Inn$A,ylim=my.ylim)
lines(Inn$Left,ylim=my.ylim,col='gray50')
lines(Inn$Right,ylim=my.ylim,col='gray50')
axis(1,at=1:length(CI$A),labels = Inn$name,las=2,cex.axis=my.cex)
title("(a)", adj = 0,line = 1)

# Cohort 2. CI for B for subset, contr vs CIs.converters
N.sel     <- n.s.ci.B.neg.12;length(N.sel)
CI        <- CI.2.B.wicla.A1.subset[N.sel,];
TEMP      <- data.frame(CI,PG.id[N.sel,])
o         <- order(TEMP$name)
CI.names  <- TEMP[o,];CI.names[1:3,]
rownames(CI.names)<- CI.names$Gene.names
Inn <- CI.names;
plot(Inn$A,xaxt='n',
     xlab='',ylab='z-scores',
     ylim=my.ylim)
abline(h=0,lty=2,col='gray20')
lines(Inn$A,ylim=my.ylim)
lines(Inn$Left,ylim=my.ylim,col='gray50')
lines(Inn$Right,ylim=my.ylim,col='gray50')
axis(1,at=1:length(CI$A),labels = Inn$name,las=2,cex.axis=my.cex)
title("(b)", adj = 0,line = 1)
dev.off()

#############################
## Fig S7 E+R demonstration
#############################


pdf('./figures/Fig S7 Ill ER Compl B.pdf')

er        <- er.1
my.cex.main <- 0.0000001
my.cls.s  <- cls.3[D.1$A.3];# ctrl.non.neur - ctrl.other.neuro - MS
my.pch.s  <- points.3[D.1$B.3]
Inn       <- M.1
D         <- data.frame(D.1,my.cls.s,my.pch.s)

er.A  <- ER.set1.A;dim(er.A)
er.B  <- ER.set1.B
er.AB <- ER.set1.AB

eff.A  <- Eff.set1.A
eff.B  <- Eff.set1.B
eff.AB <- Eff.set1.AB
res    <- ER.res

n.sel <- "P00751";gconvert(n.sel) # CFB
#n.sel <- "P01834" # Ig kappa C chain

par(mfrow=c(2,2))
my.ylim     <-  c(-3,3)
my.cex.axis <- 1.2
my.cex.lab  <- 1.2
dim(Inn)
y.sel     <- which(colnames(Inn)==n.sel);y.sel

Inn.means <- aggregate(Inn,by=list(D$gr4),mean)[D$gr4,n.sel]
plot(Inn[,n.sel],col=my.cls.s,pch=my.pch.s,
     xlab='Participants',ylab='z-scores',
     ylim=my.ylim,cex.axis =my.cex.axis,
     cex.lab=my.cex.lab)
lines(Inn.means,lwd = 2)
abline(h=0,lty=1,col='gray50')
abline(v=62.5,lty=2,col='gray50')
title("(a)", adj = 0,line = 1)

plot(res[,n.sel],col=my.cls.s,pch=my.pch.s,
     xlab='Participants',ylab='Residuals',
     ylim=my.ylim,cex.axis =my.cex.axis,
     cex.lab=my.cex.lab)
abline(h=0,lty=1,col='gray50')
abline(v=64.5,lty=2,col='gray50')
title("(b)", adj = 0,line = 1)

plot(ER.set1.A[,n.sel],col=my.cls.s,pch=my.pch.s,
     xlab='Participants',ylab='ER values of patient gr.',
     ylim=my.ylim,cex.axis =my.cex.axis,
     cex.lab=my.cex.lab)
abline(h=0,lty=1,col='gray50')
abline(v=64.5,lty=2,col='gray50')
lines(eff.A[,n.sel],lty=1,lwd = 2)
title("(c)", adj = 0,line = 1)

plot(ER.set1.B[,n.sel],col=my.cls.s,pch=my.pch.s,
     xlab='Participants',ylab='ER values of disease category',
     ylim=my.ylim,cex.axis =my.cex.axis,
     cex.lab=my.cex.lab)
abline(h=0,lty=1,col='gray50')
abline(v=64.5,lty=2,col='gray50')
lines(eff.B[,n.sel],lty=1,lwd = 2)
title("(d)", adj = 0,line = 1)

dev.off()


###############################################
# Fig S9 Normal probability plots of PLS models
###############################################

pdf('./Figures/Fig S9 Normal probability plots of PLS models a b c.pdf')

par(mfrow=c(2,3))
qqnorm(residuals.1.A,
       ylab="Standardized Residuals", xlab="Normal Scores",main='')
qqline(residuals.1.A) 
title("(a)", adj = 0,line = 1)

qqnorm(residuals.1.B.2pc,
       ylab="Standardized Residuals", xlab="Normal Scores",main='')
qqline(residuals.1.B.2pc) 
title("(b)", adj = 0,line = 1)

qqnorm(residuals.2.B.2pc,
       ylab="Standardized Residuals", xlab="Normal Scores",main='')
qqline(residuals.2.B.2pc) 
title("(c)", adj = 0,line = 1)

dev.off()



##########################################
## Fig S12 Corr across both data sets
##########################################
N.sel <- N.s.B.jkn.set1and2.consistent;length(N.sel)

# Original data
pdf('./Figures/Fig S12 Cohort 12. Original data. Corr plot sign Jkn .pdf')
Inn           <- M.12[,N.sel]
GN.sel        <- PG.id[colnames(Inn),]
G.sel         <- GN.sel$name
colnames(Inn) <- G.sel
o             <- order(G.sel)
Inn.o         <- Inn[,o];
corr.all      <- cor(Inn.o)
Inn.2corr     <- Inn
sort(colnames(Inn.2corr))

corrgram(Inn.2corr, order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="")
dev.off()


Corr.all <- cor(Inn.2corr);Corr.all[1:5,1:5]
idc      <- which(colnames(Corr.all)=='CFB')
idr      <- which(Corr.all[,idc]>0.6);g.high.corr.CFB <- rownames(Corr.all)[idr]
n.high.corr.CFB <- GN.sel[1:3,]
n.sel    <- which(substr(colnames(Corr.all),1,1)=='C')
Temp     <- Corr.all[n.sel,n.sel]
Compl    <- Temp[1:7,1:7]

# After ER modelling of the combined data
# pdf('./Figures/Fig S14 Cohort 12.ER.values. Corr plot sign Jkn .pdf')
# Inn           <- ER.set12.B[,N.sel]
# GN.sel        <- PG.id[colnames(Inn),]
# G.sel         <- GN.sel$name
# colnames(Inn) <- G.sel
# o             <- order(G.sel)
# Inn.o         <- Inn[,o];
# corr.all      <- cor(Inn.o)
# Inn.2corr     <- Inn
# 
# corrgram(Inn, order=TRUE, lower.panel=panel.shade,
#          cex.labels = 0.4,
#          upper.panel=panel.pie, text.panel=panel.txt,
#          main="")
# dev.off()


###########################################################
## Fig S13 Corr across both data sets selected proteins
###########################################################

C3    <- Inn.2corr[,'C3'] 
TF    <- Inn.2corr[,'TF']
BTD   <- Inn.2corr[,'BTD']
FSTL1 <- Inn.2corr[,'FSTL1']
CFB   <- Inn.2corr[,'CFB']
CFH   <- Inn.2corr[,'CFH']
RBP4  <- Inn.2corr[,'RBP4']
IGFBP4<- Inn.2corr[,'IGFBP4']
MDH1  <- Inn.2corr[,'MDH1']
APOH  <- Inn.2corr[,'APOH']
SAA4  <- Inn.2corr[,'SAA4']
RELN  <- Inn.2corr[,'RELN']

pdf('./Figures/Fig S13. Cohort 1 and 2 Corr Selected proteins.pdf')
my.cls.s.12 <- c('gray','black')[D.12$B.12]
my.pch.s.12 <- c(15,17)[DCG.12$C.12]

par(mfrow=c(2,2))
r  <- round(cor(C3,CFB),digits = 2);r
plot(CFB,C3,col=my.cls.s.12,pch=my.pch.s.12)
abline(lm(CFB~C3),lty=2)
text(2.4,2.77,labels=paste('r = ',r))
title("(a)", adj = 0,line = 1)

r  <- round(cor(C3,CFH),digits = 2);r
plot(CFH,C3,col=my.cls.s.12,pch=my.pch.s.12)
abline(lm(CFH~C3),lty=2)
text(2.27,2.65,labels=paste('r = ',r))
title("(b)", adj = 0,line = 1)

r  <- round(cor(C3,TF),digits = 2);r
plot(TF,C3,col=my.cls.s.12,pch=my.pch.s.12)
abline(lm(C3~TF),lty=2)
text(2.3,2.65,labels=paste('r = ',r))
title("(c)", adj = 0,line = 1)

r  <- round(cor(TF,RBP4),digits = 2);r
plot(TF,RBP4,col=my.cls.s.12,pch=my.pch.s.12)
abline(lm(RBP4~TF),lty=2)
text(2.35,2.65,labels=paste('r = ',r))
title("(d)", adj = 0,line = 1)

dev.off()

r  <- round(cor(C3,RBP4),digits = 2);r
plot(RBP4,C3,col=my.cls.s.12,pch=my.pch.s.12)
abline(lm(C3~RBP4),lty=2)
text(2.55,2.6,labels=paste('r = ',r))
title("(d)", adj = 0,line = 1)



###########################
# Means of means
###########################

# Means.Means.data.Ig  <- apply(Means.data.Ig,2,mean);Means.Means.data.Ig
# Means.Means.ER.A.Ig  <- apply(Means.ER.A.Ig,2,mean);Means.Means.ER.A.Ig
# Means.Means.ER.B.Ig  <- apply(Means.ER.B.Ig,2,mean);Means.Means.ER.B.Ig
# Means.Means.data.sign.MS  <- apply(Means.data.sign.MS,2,mean);Means.Means.data.sign.MS
# Means.Means.ER.A.sign.MS  <- apply(Means.ER.A.sign.MS,2,mean);Means.Means.ER.A.sign.MS
# Means.Means.ER.B.sign.MS  <- apply(Means.ER.B.sign.MS,2,mean);Means.Means.ER.B.sign.MS
# par(mfrow=c(2,3))
# my.ylim <- c(-1.5,1.5)
# barplot(Means.Means.data.Ig,
#         ylim = my.ylim,
#         col=c('blue','blue','red','red'))
# barplot(Means.Means.ER.A.Ig,
#         ylim = my.ylim,
#         col=c('blue','blue','red','red'))
# barplot(Means.Means.ER.B.Ig,
#         ylim = my.ylim,
#         col=c('blue','blue','red','red'))
# barplot(Means.Means.data.sign.MS,
#         ylim = my.ylim,
#         col=c('blue','blue','red','red'))
# barplot(Means.Means.ER.A.sign.MS,
#         ylim = my.ylim,
#         col=c('blue','blue','red','red'))
# barplot(Means.Means.ER.B.sign.MS,
#         ylim = my.ylim,
#         col=c('blue','blue','red','red'))


#######################################
## Fig 4 PLS plots, common proteins
#######################################

pdf('./Figures/Fig 4 PLS-DA  A and B cohort 1 and B cohort 2.pdf')

par(mfrow=c(3,2))


par(mfrow=c(2,2))
# set 1 factor B
my.cls.s    <- cls.3[D.1$A.3];# ctrl.non.neur - ctrl.other.neuro - MS
my.pch.s    <- points.3[D.1$B.3]
#my.cls.s    <- cls[D.1$A.12]
#my.pch.s    <- points[D.1$B.12]
plot(scores.1.B.sel[,1:2],col=my.cls.s,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title("(a)", adj = 0,line = 1)

plot(loadings.1.B.sel[,1:2],
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title("(b)", adj = 0,line = 1)

my.cls.s    <- cls.2[D.2$A.12]
my.pch.s    <- points.2[D.2$B.12]
plot(scores.2.B.sel[,1:2],col=my.cls.s,pch=my.pch.s,
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
#     labels=Names.set2)
abline(v=0,h=0,lty=2,col='gray50')
title("(c)", adj = 0,line = 1)

plot(loadings.2.B.sel[,1:2],
     xlab= 'PLS factor 1',
     ylab= 'PLS factor 2')
abline(v=0,h=0,lty=2,col='gray50')
title("(d)", adj = 0,line = 1)

dev.off()



###########################################
## CI.12 
###########################################

pdf('./Figures/Fig set 12 set1.ER.set2 a.CI b.gender.pdf')

my.ylim <- c(-2,2)
my.cex <- 0.6
# Cohort 1. CI for B within A1
par(mfrow=c(2,1))

N.sel     <- N.s.B.jkn.set1and2.consistent;length(N.sel)
CI        <- CI.12.B[N.sel,];
TEMP      <- data.frame(CI,PG.id[N.sel,])
o         <- order(TEMP$name)
CI.names  <- TEMP[o,];rownames(CI.names)<- CI.names$Gene.names
Inn       <- CI.names;colnames(Inn)
plot(Inn$A,xaxt='n',
     xlab='',ylab='z-scores',
     ylim=my.ylim)
abline(h=0,lty=2,col='gray20')
lines(Inn$A,ylim=my.ylim,lwd=2)
lines(Inn$Left,ylim=my.ylim,col='gray50')
lines(Inn$Right,ylim=my.ylim,col='gray50')
axis(1,at=1:length(CI$A),labels = Inn$name,las=2,cex.axis=my.cex)
title("(a)", adj = 0,line = 1)

# Plot means of females, and men.r

# set 12, females
N.sel       <- N.s.B.jkn.set1and2.consistent;length(N.sel)
my.array    <- My.Array.12;names(my.array)
y           <- DCG.12$B;length(y)
er.12       <- ER(M.12 ~ factorC*factorB, data = my.array)
er          <- er.12
ER.12.C     <- as.data.frame(unclass(er$ER.values$factorC))
ER.12.B     <- as.data.frame(unclass(er$ER.values$factorB))
ER.12.CB    <- as.data.frame(unclass(er$ER.values$`factorC:factorB`))
id          <- which(DCG.12$G=='f');n.g.f <- rownames(DCG.12)[id]
id          <- which(DCG.12$G=='m');n.g.m <- rownames(DCG.12)[id]
Inn.f       <- ER.12.B[n.g.f,N.sel] # selected in both sets, consistent pattern
Inn.m       <- ER.12.B[n.g.m,N.sel] # selected in both sets, consistent pattern
DCG.12.f    <- DCG.12[n.g.f,]
DCG.12.m    <- DCG.12[n.g.m,]
dim(Inn.f);dim(Inn.m)

Means.ER.12.B <- t(aggregate(ER.12.B[,N.sel],by=list(DCG.12$B.12),mean)[,-1]);dim(Means.ER.12.B)
Means.f       <- t(aggregate(Inn.f,by=list(DCG.12.f$B.12),mean)[,-1])
Means.m       <- t(aggregate(Inn.m,by=list(DCG.12.m$B.12),mean)[,-1])
Diff.ER.12.B  <- Means.ER.12.B[,2]- Means.ER.12.B[,1]
Diff.f        <- Means.f[,2]-Means.f[,1]
Diff.m        <- Means.m[,2]-Means.m[,1]
Diff      <- cbind(Diff.ER.12.B,Diff.f,Diff.m);dim(Diff)
TEMP2     <- data.frame(PG.id[rownames(Diff),],Diff)
o         <- order(TEMP2$name)
TEMP3     <- TEMP2[o,];
rownames(TEMP3) <- TEMP3$name
Inn             <- TEMP3;

plot(Inn$Diff.ER.12.B,ylim=my.ylim,xaxt='n',
     xlab='',ylab='z-scores')
lines(Inn$Diff.ER.12.B,col='black',ylim=my.ylim,lwd=2)
lines(Inn$Diff.f,col='blue',ylim=my.ylim)
lines(Inn$Diff.m,col='orange',ylim=my.ylim)
abline(h=0,lty=2,col='gray30')
axis(1,at=1:length(Diff.f),labels = rownames(Inn),las=2,cex.axis=my.cex)
title("(b)", adj = 0,line = 1)

dev.off()



##################################
# CI of combined datasets for A2
##################################

# sign for PLS for M.12 and by CI within A1 A2, by directly combinding the cohort
N.sel     <- intersect(n.s.CI.12.B.A2,N.s.B.jkn.set1and2.consistent);length(N.sel)

N.sel     <- n.s.CI.12.within.A1.A2;length(N.sel)


pdf('./Figures/Fig set 12 A1 and A2. directly combined for both sets common for A1 A2.pdf')

my.ylim <- c(-2,2)
my.cex <- 0.6
par(mfrow=c(2,1))

# A1
my.array      <- My.Array.12.A1;names(my.array)
CI.12.B.A1    <- with(my.array,confints(M.12.A1[factorB == -1, ],M.12.A1[factorB ==  1, ]))
CI.12.B.s.A1  <- run.ci(CI=CI.12.B.A1)
n.s.CI.12.B.A1   <- CI.12.B.s.A1[[3]]
CI        <- CI.12.B.A1[N.sel,];
TEMP      <- data.frame(CI,PG.id[N.sel,])
o         <- order(TEMP$name)
CI.names  <- TEMP[o,];rownames(CI.names)<- CI.names$Gene.names
Inn       <- CI.names;colnames(Inn)
plot(Inn$A,xaxt='n',
     xlab='',ylab='z-scores',
     ylim=my.ylim)
abline(h=0,lty=2,col='gray20')
lines(Inn$A,ylim=my.ylim,lwd=2)
lines(Inn$Left,ylim=my.ylim,col='gray50')
lines(Inn$Right,ylim=my.ylim,col='gray50')
axis(1,at=1:length(CI$A),labels = Inn$name,las=2,cex.axis=my.cex)
title("(a)", adj = 0,line = 1)

# A2
my.array      <- My.Array.12.A2;names(my.array)
CI.12.B.A2    <- with(my.array,confints(M.12.A2[factorB == -1, ],M.12.A2[factorB ==  1, ]))
CI.12.B.s.A2  <- run.ci(CI=CI.12.B.A2)
n.s.CI.12.B.A2   <- CI.12.B.s.A2[[3]]
CI        <- CI.12.B.A2[N.sel,];
TEMP      <- data.frame(CI,PG.id[N.sel,])
o         <- order(TEMP$name)
CI.names  <- TEMP[o,];rownames(CI.names)<- CI.names$Gene.names
Inn       <- CI.names;colnames(Inn)
plot(Inn$A,xaxt='n',
     xlab='',ylab='z-scores',
     ylim=my.ylim)
abline(h=0,lty=2,col='gray20')
lines(Inn$A,ylim=my.ylim,lwd=2)
lines(Inn$Left,ylim=my.ylim,col='gray50')
lines(Inn$Right,ylim=my.ylim,col='gray50')
axis(1,at=1:length(CI$A),labels = Inn$name,las=2,cex.axis=my.cex)
title("(b)", adj = 0,line = 1)
dev.off()

##################################
# Save worksheet
##################################
save.image("H:/A_R/R/R data/R MS publ.SR/Worksheets/MS data loaded and analysed.RData")


##################################
# Corr-pattern selected proteins
##################################


# Original data
pdf('./Figures/Fig New. Cohort 1. corrgram of 24 proteins.pdf')
GN.sel <- gconvert(n.1.s.ci.B.wicla)[,2:6] # 24 proteins neg for MS in cohort 1
N.sel <- GN.sel$input
G.sel <- GN.sel$name
x     <- GN.sel$description
P.sel <- sub('\\[.*', '', x)


Inn           <- ER.set1.B[,N.sel];dim(Inn)
colnames(Inn) <- G.sel
o             <- order(G.sel)
Inn.o         <- Inn[,o];colnames(Inn.o)
corr.all      <- cor(Inn.o)
Inn.2corr     <- Inn.o
colnames(Inn.2corr)

corrgram(Inn.2corr, order=FALSE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="")
dev.off()

Inn.o.means.2.gr <- t(aggregate(Inn.o,by=list(D.1$B.12),mean))  
Inn.o.means.2.gr

Table.24.proteins.set1 <- data.frame(N.sel[o],G.sel[o],P.sel[o],Inn.o.means.2.gr[-1,])
write.table(Table.24.proteins.set1,file='./Results/Table.24.proteins.set1.txt',sep='\t')


Inn             <- M.1[,N.sel]
Inn.means.4.gr  <- aggregate(Inn,by=list(D.1$gr4),mean)

par(mfrow=c(3,4))
n.sel <- 'P06681'
g.sel           <- gconvert(n.sel)[,5]
inn             <- Inn.means.4.gr[,n.sel]
barplot(inn,col=c('blue','red','blue','red'))
abline(v=2.5,lty=2,col='gray')
title(paste0('(a) ', g.sel), adj = 0,line = 1)

n.sel           <- 'P00751'
g.sel           <- gconvert(n.sel)[,5]
inn             <- Inn.means.4.gr[,n.sel]
barplot(inn,col=c('blue','red','blue','red'))
abline(v=2.5,lty=2,col='gray')
title(paste0('(b) ', g.sel), adj = 0,line = 1)

n.sel           <- 'Q9NZP8'
g.sel           <- gconvert(n.sel)[,5]
inn             <- Inn.means.4.gr[,n.sel]
barplot(inn,col=c('blue','red','blue','red'))
abline(v=2.5,lty=2,col='gray')
title(paste0('(c) ', g.sel), adj = 0,line = 1)

n.sel           <- 'P05156'
g.sel           <- gconvert(n.sel)[,5]
inn             <- Inn.means.4.gr[,n.sel]
barplot(inn,col=c('blue','red','blue','red'))
abline(v=2.5,lty=2,col='gray')
title(paste0('(d) ', g.sel), adj = 0,line = 1)

n.sel           <- 'P02787'
g.sel           <- gconvert(n.sel)[,5]
inn             <- Inn.means.4.gr[,n.sel]
barplot(inn,col=c('blue','red','blue','red'))
abline(v=2.5,lty=2,col='gray')
title(paste0('(e) ', g.sel), adj = 0,line = 1)

n.sel           <- 'P02774'
g.sel           <- gconvert(n.sel)[,5]
inn             <- Inn.means.4.gr[,n.sel]
barplot(inn,col=c('blue','red','blue','red'))
abline(v=2.5,lty=2,col='gray')
title(paste0('(f) ', g.sel), adj = 0,line = 1)

n.sel           <- 'P02647'
g.sel           <- gconvert(n.sel)[,5]
inn             <- Inn.means.4.gr[,n.sel]
barplot(inn,col=c('blue','red','blue','red'))
abline(v=2.5,lty=2,col='gray')
title(paste0('(g) ', g.sel), adj = 0,line = 1)

n.sel           <- 'O75882'
g.sel           <- gconvert(n.sel)[,5]
inn             <- Inn.means.4.gr[,n.sel]
barplot(inn,col=c('blue','red','blue','red'))
abline(v=2.5,lty=2,col='gray')
title(paste0('(h) ', g.sel), adj = 0,line = 1)


#################
# PCA
#################

# N.sel       <- N.s.B.jkn.set1and2.consistent;length(N.sel)
# M.1.c.a1    <- M.1.a1[,colnames(M.1.c)]
# M.1.c.a2    <- M.1.a2[,colnames(M.1.c)]
# Inn.all   <- rbind(M.1.c.a1,M.2.c.a1);dim(Inn.all)
# Inn       <- Inn.all[,N.sel];dim(Inn)
# pca.mod   <- prcomp(Inn,scale=TRUE)
# scores    <- pca.mod$x
# loadings  <- pca.mod$rotation
# expl.var  <- round(propVar(pca.mod)*100,digits = 0);expl.var[1:8]
# D         <- DCG.12[rownames(Inn),];dim(D)
# my.cls.s  <- cls.2[D$B.12];my.cls.s
# my.pch.s  <- points.2[D$C.12];my.pch.s
# plot(scores[,1],
#      main=paste0('PCA scores'),
#      xlab=paste0('Participants'),
#      ylab=paste0('PC1:',expl.var[1],'%'),
#      col=my.cls.s,pch=my.pch.s)
# abline(h=0,v=0,lty=1,col='gray50')
# title("(a)", adj = 0,line = 1)
# plot(loadings[,1],
#      main=paste0('PCA loadings'),
#      xlab=paste0('PC1:',expl.var[1],'%'),
#      ylab=paste0('PC2:',expl.var[2],'%'))
# abline(h=0,v=0,lty=1,col='gray50')
# title("(b)", adj = 0,line = 1)

# Inn.all   <- rbind(M.1.c.a2,M.2.c.a2);dim(Inn.all)
# Inn       <- Inn.all[,N.sel];rownames(Inn)
# pca.mod   <- prcomp(Inn,scale=TRUE)
# scores    <- pca.mod$x
# loadings  <- pca.mod$rotation
# expl.var  <- round(propVar(pca.mod)*100,digits = 0);expl.var[1:8]
# D         <- DCG.12[rownames(Inn),];rownames(DCG.12)
# my.cls.s  <- cls.2[D$B.12];my.cls.s
# my.pch.s  <- points.2[D$C.12];my.pch.s
# plot(scores[,1],
#      main=paste0('PCA scores'),
#      xlab=paste0('PC1:',expl.var[1],'%'),
#      ylab=paste0('PC2:',expl.var[2],'%'),
#      col=my.cls.s,pch=my.pch.s)
# abline(h=0,v=0,lty=1,col='gray50')
# title("(c)", adj = 0,line = 1)
# plot(loadings[,1],
#      main=paste0('PCA loadings'),
#      xlab=paste0('PC1:',expl.var[1],'%'),
#      ylab=paste0('PC2:',expl.var[2],'%'))
# abline(h=0,v=0,lty=1,col='gray50')
# title("(d)", adj = 0,line = 1)

