#STReEF05_Rmodel
# This R code is developed for the manuscript 'Structural and
# Effective brain connectivity in focal epilepsy' by Jelsma et al. 

# author: Susanne Jelsma
# date: October 2021

# calculate the jaccard index, construct a linear multilevel model
-----------------------------------------------------------------------
 # SECTION 1 : calculate the jaccard index
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


# compute expected jaccard similarity coefficient under independence
library(jaccard)
library(xlsx)
library(R.matlab)

dir = '/blabla/shareData_STReEF/derivatives/'

i = 0
expected = integer(1) # the expected Jaccard Index
score = integer(1) # the Jaccard Index
infos = data.frame(matrix(ncol=5,nrow=1)) # present the statisctics of the Jaccard Index in a table
colnames(infos) = c('statistics','pvalue', 'expectation','accuracy','error.type')
patients = c('STREEF01','STREEF02','STREEF03','STREEF04','STREEF05','STREEF06','STREEF07','STREEF08','STREEF09','STREEF10','STREEF11','STREEF12','STREEF13')

for (pat in patients)
{ 
  i = i+1
  datapath_SC = paste(dir, sprintf("sub-%s/sub-%s_ses-1_Structural_Connectivity.mat",pat,pat),sep = '')
  datapath_EC = paste(dir, sprintf("sub-%s/sub-%s_ses-1_Effective_Connectivity.mat",pat,pat),sep = '')
  
  SC = readMat(datapath_SC) # structural network
  EC = readMat(datapath_EC) # effective network
  SC_matrix = SC$SC.matrix
  EC_matrix = EC$EC.matrix
  
  score[i] <- jaccard(SC_matrix,EC_matrix,center =FALSE,px=NULL,py=NULL) # the Jaccard Index
  expected[i] <- jaccard.ev(SC_matrix,EC_matrix,px=NULL,py=NULL) # the expected Jaccard Index
  test = jaccard.test.mca(SC_matrix,EC_matrix,accuracy=1e-05,verbose=TRUE) 
  infos[i,] = test} # statistical information about the Jaccard Index

jaccard_statistiek = t.test(infos[,1],alternative= 'g')

p_adjust = p.adjust(infos[,2], method = 'fdr', n = length(score))
reject = p_adjust<0.05

-----------------------------------------------------------------------
  #SECTION 2: construct a linear multilevel model
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
dir = '/blabla/shareData_STReEF/derivatives/'
path = paste(dir,'input_LMM_model.csv',sep = '')
data_long = read.csv(file=path,header=TRUE)

# construct a linear multilevel model with backward elimination of possible predictors. The possible predictors at patient level were the effective network characteristic, the node proximity, the volume of the structural electrode contact areas, and the nodes in the seizure onset zone (SOZ). We added each variable as fixed effect to the model. The variable with the highest p-value was removed at each step until all variables had a p-value<0.05. 
df = data.frame(matrix(nrow = 5, ncol = 14))
names(df) = c("lmm1_predictor", "lmm1_t_value", "lmm2_predictor", "lmm2_t_value", "lmm3_predictor", "lmm3_t_value", "lmm_final_predictor", "lmm_final_t_value", "lmm_final_coef", "lmm_final_CI_low", "lmm_final_CI_high", "lmm_final_p" , 'statistical_tresh', 'ICC')

# intercept only model
interceptonlymodel = lmer(formula = SCD ~ 1 + (1|subj) , data = data_long) 
summary(interceptonlymodel) #to get parameter estimates.

# ICC
RandomEffects = as.data.frame(VarCorr(interceptonlymodel))
df$ICC[1] = RandomEffects[1,4]/(RandomEffects[1,4]+RandomEffects[2,4]) 

# model with all predictors and 11 subjects degree (Degree effective networks) + Node proximity + volume(Volume electrode areas) + epi(SOZ nodes) 
data_long_11 = subset(data_long,!(subj == 9|subj==13))

model1 = lmer(formula = SCD ~ 1 + ECD + NP + SOZ + VEA + (1|subj), 
              data    = data_long_11) 
summary(model1)
model1_par = parameters::parameters(model1)
df$lmm1_predictor = model1_par$Parameter[1:5]
df$lmm1_t_value = model1_par$t[1:5]


model2 = lmer(formula = SCD ~ 1 + ECD + NP + SOZ + (1|subj), 
               data    = data_long_11) # without VEA
summary(model2)
model2_par = parameters::parameters(model2)
df$lmm2_predictor[1:4] = model2_par$Parameter[1:4]
df$lmm2_t_value[1:4] = model2_par$t[1:4]

model3 = lmer(formula = SCD ~ 1 + ECD + NP + (1|subj), #final model 11 pt
               data    = data_long_11)
summary(model3)
model3_par = parameters::parameters(model3)
df$lmm3_predictor[1:3] = model3_par$Parameter[1:3]
df$lmm3_t_value[1:3] = model3_par$t[1:3]

model4 = lmer(formula = SCD ~ 1 + ECD + NP + (1|subj), #final model all pt
               data    = data_long)
summary(model4)

model4_par = parameters::parameters(model4)
df$lmm_final_predictor[1:3] = model4_par$Parameter[1:3]
df$lmm_final_t_value[1:3] = model4_par$t[1:3]
df$lmm_final_coef[1:3] = model4_par$Coefficient[1:3]
df$lmm_final_CI_low[1:3] = model4_par$CI_low[1:3]
df$lmm_final_CI_high[1:3] = model4_par$CI_high[1:3]
df$lmm_final_p[1:3] = model4_par$p[1:3]


# statistica; treshold t-value with p-value 0.01 and df around 670
df$statistical_tresh[1] = qt(0.05,670)

path = paste(dir,'output_LMM_model.csv',sep = '')
write.csv(df, path)
