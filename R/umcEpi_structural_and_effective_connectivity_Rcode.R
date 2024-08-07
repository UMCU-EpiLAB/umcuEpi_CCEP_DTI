# STReEF05_compare_networks_R
# inter-modal similarity and network topography comparison between structural (derived from DWI) and effective (derived from SPES) networks

# author: Susanne Jelsma & Dorien van Blooijs
# date: May 2022

# Be aware! This file has a twin written in matlab code and to execute the code correctly, run first the file: 'umcuEpi_structural_and_effective_connectivity' ,and run the sections of this file when indicated.

# calculate the jaccard index, construct a linear multilevel model

-----------------------------------------------------------------------
SECTION 1 : calculate the jaccard index
-----------------------------------------------------------------------
# install packages (only needed if you run this file for the first time)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("qvalue"))

install.packages("jaccard",dependencies = TRUE)
install.packages("Rcpp",dependencies = TRUE)
install.packages("dplyr",dependencies = TRUE)
install.packages("magrittr",dependencies = TRUE)
install.packages("R.matlab",dependencies = TRUE)
install.packages("xlsx")

# check if the package is present in the list of available packages (only needed if you run this file for the first time)
av <- available.packages(filters=list())
av[av[, "Package"] == jaccard]

# compute expected jaccard similarity coefficient under independence
library(jaccard)
library(xlsx)

i <- 0
expected <- integer(1) # the expected Jaccard Index
score <- integer(1) # the Jaccard Index
infos <- data.frame(matrix(ncol=5,nrow=1)) # present the statisctics of the Jaccard Index in a table
colnames(infos) <- c('statistics','pvalue', 'expectation','accuracy','error.type')
patients <- c('RESP0703','RESP0749', 'RESP0753','RESP0778','RESP0856','RESP0882','RESP0905','RESP0928','RESP0968','RESP0978','RESP0984','RESP0992','RESP0993')
for (pat in patients)
{ 
i <- i+1
datapath_dwi <- sprintf("/Fridge/users/susanne/derivatives/STReEF/pipeline_all_patients/dwi_matlab/sub-%s/SC_matrix_dwi",pat) # see personalDataPath_mrtrix file for the exact dir (**)
datapath_ccep <- sprintf("/Fridge/users/susanne/derivatives/STReEF/pipeline_all_patients/dwi_matlab/sub-%s/SC_matrix_ccep",pat)

#datapath_dwi <- sprintf("/Fridge/users/susanne/derivatives/dwi_matlab2/sub-%s/SC_matrix_dwi",pat) # see personalDataPath_mrtrix file for the exact dir (**)
#datapath_ccep <- sprintf("/Fridge/users/susanne/derivatives/dwi_matlab2/sub-%s/SC_matrix_ccep",pat)
dwi <- read.table(datapath_dwi) # structural network
spes <- read.table(datapath_ccep) # effective network
dwi2 <- as.matrix(dwi)
spes2 <- as.matrix(spes)
score[i] <- jaccard(dwi2,spes2,center =FALSE,px=NULL,py=NULL) # the Jaccard Index
expected[i] <- jaccard.ev(dwi2,spes2,px=NULL,py=NULL) # the expected Jaccard Index
k = c(1,6,8,9,10,2,3,4,5,7,11,12,13);
score[k]
expected[k]
test <- jaccard.test.mca(dwi2,spes2,accuracy=1e-05,verbose=TRUE) 
infos[i,] <- test} # statistical information about the Jaccard Index

tc = infos[,1]
mean(tc)
jaccard_statistiek = t.test(tc,alternative= 'g')

# present the Jaccard Index per electrode type (only if you load in all 13 patients)
grid <- c(1,6,8,9,10);
seeg <- c(2,3,4,5,7,11,12,13);
grid_info <- infos[grid,]
seeg_info <- infos[seeg,]

# save the results in one xlsx file with different sheets
write.xlsx(infos,file ="*dwi_matlab*jaccard_index.xlsx")
write.xlsx(grid_info,file ="/*dwi_matlab*/jaccard_index.xlsx",sheetName='grid',append=TRUE)
write.xlsx(seeg_info,file ="*dwi_matlab*/jaccard_index.xlsx",sheetName='seeg',append=TRUE)
write.xlsx(score,file ="*dwi_matlab*/jaccard_index.xlsx",sheetName='score',append=TRUE)
write.xlsx(expected,file ="*dwi_matlab*/jaccard_index.xlsx",sheetName='expected',append=TRUE)

-----------------------------------------------------------------------
SECTION 2: construct a linear multilevel model
-----------------------------------------------------------------------
# install packages (only needed if you run this file for the first time)
install.packages("swirl")
install.packages("lrstat")
install.packages("parameters")
install.packages("lme4")

install.packages("sjstats",dependencies = TRUE)
install.packages("car",dependencies = TRUE)
install.packages("carData",dependencies = TRUE)
install.packages(c("backports","effects","ggplot2","interactions","lme4","lmerTest","psych","plyr"))

library(backports)     # to revive the isFALSE() function for sim_slopes()
library(effects)       # for probing interactions
library(ggplot2)       # for data visualization
library(interactions)  # for probing/plotting interactions
library(lme4)          # for multilevel models
library(lmerTest)      # for p-values
library(psych)         # for describing the data
library(plyr)          # for data manipulation
library("car")
library("swirl")
library("carData")
library(lrstat)

library(parameters)
library(lme4)

# read in the matrix file 
dir = '/Fridge/users/susanne/derivatives/STReEF/pipeline_all_patients/dwi_matlab/'
path = paste(dir,'data_long_LMM.csv',sep = '')
data_long = read.csv(file=path,header=TRUE)

# construct a linear multilevel model with backward elimination of possible predictors. The possible predictors at patient level were the effective network characteristic, the node proximity, the volume of the structural electrode contact areas, and the nodes in the seizure onset zone (SOZ). We added each variable as fixed effect to the model. The variable with the highest p-value was removed at each step until all variables had a p-value<0.05. 

# intercept only model
interceptonlymodel = lmer(formula = SCD ~ 1 + (1|subj) , data = data_long) 
summary(interceptonlymodel) #to get parameter estimates.

# ICC
RandomEffects = as.data.frame(VarCorr(interceptonlymodel))
ICC_between = RandomEffects[1,4]/(RandomEffects[1,4]+RandomEffects[2,4]) 
ICC_between

# model with all predictors and 11 subjects degree (Degree effective networks) + Node proximity + volume(Volume electrode areas) + epi(SOZ nodes) 
data_long_11 = subset(data_long,!(subj == 9|subj==13))

model1 = lmer(formula = SCD ~ 1 + ECD + NP + VEA + SOZ +(1|subj), 
               data    = data_long_11) 
summary(model1)
model1_par = parameters::parameters(model1)
t1 = model1_par$t
p1 = model1_par$Parameter

model2 <- lmer(formula = SCD ~ 1 + ECD + NP + SOZ + (1|subj), 
               data    = data_long_11) # without VEA
summary(model2)
model2_par = parameters::parameters(model2)


model3 <- lmer(formula = SCD ~ 1 + ECD + NP + (1|subj), #final model 11 pt
               data    = data_long_11)
summary(model3)
model3_par = parameters::parameters(model3)

model4 <- lmer(formula = SCD ~ 1 + ECD + NP + (1|subj), #final model all pt
               data    = data_long)
summary(model4)

model4_par = parameters::parameters(model4)

# statistica; treshold t-value with p-value 0.01 and df around 670
tresh = qt(0.05,670)

data = data.frame('Statistical Treshold'= tresh, 'IC' = ICC_between, 'model1 par' = p1, 'model1 t' = t1, 'final_model_par'=model4_par)

path = paste(dir,'LLMdata_R.csv',sep = '')
write.csv(data, path)







