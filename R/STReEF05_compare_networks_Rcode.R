# STReEF05_compare_networks_R
# inter-modal similarity and network topography comparison between structural (derived from DWI) and effective (derived from SPES) networks

# author: Susanne Jelsma & Dorien van Blooijs
# date: May 2022

# Be aware! This file has a twin written in matlab code and to execute the code correctly, run first the file: ' STReEF05_compare_networks' ,and run the sections of this file when indicated.

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
# example jaccard
set.seed(1234)
x = rbinom(100,1,.5)
y = rbinom(100,1,.5)
jaccard.test.mca(x,y,accuracy = 1e-05)


# compute expected jaccard similarity coefficient under independence
library(jaccard)
library(xlsx)

i <- 0
expected <- integer(13) # the expected Jaccard Index
score <- integer(13) # the Jaccard Index
infos <- data.frame(matrix(ncol=5,nrow=13)) # present the statisctics of the Jaccard Index in a table
colnames(infos) <- c('statistics','pvalue', 'expectation','accuracy','error.type')
patients <- c('RESP0703','RESP0749','RESP0753','RESP0778','RESP0856','RESP0882','RESP0905','RESP0928','RESP0968','RESP0978','RESP0984','RESP0992','RESP0993')
for (pat in patients)
{ 
i <- i+1
datapath_dwi <- sprintf("*dwi_matlab*/sub-%s/SC_matrix_dwi",pat) # see personalDataPath_mrtrix file for the exact dir (**)
datapath_ccep <- sprintf("*dwi_matlab*/sub-%s/SC_matrix_ccep",pat)
dwi <- read.table(datapath_dwi) # structural network
spes <- read.table(datapath_ccep) # effective network
dwi2 <- as.matrix(dwi)
spes2 <- as.matrix(spes)
score[i] <- jaccard(dwi2,spes2,center =FALSE,px=NULL,py=NULL) # the Jaccard Index
expected[i] <- jaccard.ev(dwi2,spes2,px=NULL,py=NULL) # the expected Jaccard Index
test <- jaccard.test.mca(dwi2,spes2,accuracy=1e-05,verbose=TRUE) 
infos[i,] <- test} # statistical information about the Jaccard Index

# present the Jaccard Index per electrode type
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
install.packages("sjstats",dependencies = TRUE)
install.packages("car",dependencies = TRUE)
library("car")

# read in the matrix file 
data_long <- read.csv(file="*R*/data_long.csv",header=TRUE)
describe(data_long$degreeS)

# construct a linear multilevel model with backward elimination of possible predictors. The possible predictors at patient level were the effective network characteristic, the node proximity, the volume of the structural electrode contact areas, and the nodes in the seizure onset zone (SOZ). We added each variable as fixed effect to the model. The variable with the highest p-value was removed at each step until all variables had a p-value<0.05. 

# intercept only model
interceptonlymodel <- lmer(formula = degreeSC ~ 1 + (1|subj) , data = data_long) 
summary(interceptonlymodel) #to get paramater estimates.

# ICC
RandomEffects <- as.data.frame(VarCorr(interceptonlymodel))
ICC_between <- RandomEffects[1,4]/(RandomEffects[1,4]+RandomEffects[2,4]) 
ICC_between

# model with all predictors degree (Degree effective networks) + distance (Node proximity) + volume(Volume electrode areas) + epi(SOZ nodes) 

model1 <- lmer(formula = degreeSC ~ 1 + degreeEC + distance + volume + epi +(1|subj), 
               data    = data_long)
summary(model1)
RandomEffects <- as.data.frame(VarCorr(model1))
ICC_unexplained <- RandomEffects[1,4]/(RandomEffects[1,4]+RandomEffects[2,4]) 
ICC_unexplained

#model without volume predictor  degree (Degree effective networks) + distance (Node proximity) + epi(SOZ nodes) 
model2 <- lmer(formula = degreeSC ~ 1 + degreeEC + distance*epi + epi + (1|subj), 
               data    = data_long)
summary(model2)
RandomEffects <- as.data.frame(VarCorr(model2))
ICC_unexplained <- RandomEffects[1,4]/(RandomEffects[1,4]+RandomEffects[2,4]) 
ICC_unexplained

vif(model2)


#model without volume and SOZ nodes predictor  degree (Degree effective networks) + distance (Node proximity)  
model3 <- lmer(formula = degreeSC ~ 1 + degreeEC*epi + distance + (1|subj), 
               data    = data_long)
model3 <- lmer(formula = degreeSC ~ 1 + degreeEC + distance + (1|subj), 
               data    = data_long)
summary(model3)
RandomEffects <- as.data.frame(VarCorr(model3))
ICC_unexplained <- RandomEffects[1,4]/(RandomEffects[1,4]+RandomEffects[2,4]) 
ICC_unexplained

confint(model3)
vif(model3)

library(ggplot2)
p <- ggplot(data_long, aes(x = degreeEC, y = degreeSC, colour = subj)) +
  geom_point(size=3) +
  geom_line(aes(y = predict(model3)),size=1) 
print(p)


# check if the data was suitable for the statistics performed
## pearson
model6 = cor(data_long$degreeSC,data_long$degreeEC,method="pearson")
model6 = cor(data_long$degreeSC[data_long$subj == 703],data_long$degreeEC[data_long$subj == 703],method="spearman")

regrline <- lm(data_long$degreeSC[data_long$subj == 856] ~ data_long$degreeEC[data_long$subj == 856], data_long)
regrline <- lm(data_long$degreeSC[data_long$subj == 703] ~ data_long$distance[data_long$subj == 703], data_long)

qqnorm(resid(regrline)) 
qqline(resid(regrline), col = "red") # add a perfect fit line

## normality
model6 = regrline
plot(fitted(model6), resid(model6, type = "pearson"))# this will create the plot
abline(0,0, col="red")

qqnorm(resid(model6)) 
qqline(resid(model6), col = "red") # add a perfect fit line

qqnorm(ranef(model6)$subj[,1] )
qqline(ranef(model6)$subj[,1], col = "red")

## plot
plot(degreeSC~degreeEC,data_long)
plot(degreeSC[data_long$subj == 978]~degreeEC[data_long$subj == 978],data_long) # 1 patient
plot(jitter(degreeSC,4)~degreeEC,data_long)
regrline <- lm(degreeSC ~ degreeEC, data_long)
abline(regrline, lwd=3, col='red')
summary(regrline)

cor(data_long$degreeSC,data_long$degreeEC,method="spearman")
cor(data_long$degreeSC[data_long$subj == 978],data_long$degreeEC[data_long$subj == 978],method="spearman")


