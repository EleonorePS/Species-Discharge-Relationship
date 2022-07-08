#code info-------------------------------------------------

# author: Eleonore Pierrat - July 2022
# citation : Global modelling of water consumption impacts on freshwater fish biodiversity in LCA, submitted
# this code performed the K10 cross validation to select the SDR model (data preprocessing done elsewhere) 

#import libraries-------------------------

#for SDR fit
library(ggplot2)
library(devtools)
library(psych)
library(DescTools)
library(dplyr)
library(car)
library(stats)
library(rpart)
library(openxlsx)
library(purrr)
library(caret)
library(MuMIn)
library(jtools)
#install_github("vqv/ggbiplot")
#install_github("vbarbarossa/valerioUtils")
library(valerioUtils)
library(raster)
library(lme4)
library(openxlsx)
#for sensitivity
library(splitTools)
library(ape)
library(spdep)

search()
rm(list=ls())


windowsFonts(Font=windowsFont("Arial"))

#set working directory------------------

getwd()
setwd("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/")
#results will be saved there

#Import files -------------------------

#preprocessed dataset used for fitting the SDR
d_pcrglob<- read.csv("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/output/dataset_SDR_no_filter_20220316.csv")

#prerpocessed dataset for sensitivity analysis
d_sensitivity<- read.csv("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/output/dataset_SDR_GSIM_20211006.csv")#recalc

#observed discharge and species richness for sensitivity
d<- read.csv("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/schipper_barbarossa_2021/input_tab_divAREA_ep.csv")

#equivalent keys between pcrglob and GSIM
joint<- read.xlsx("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/output/joint.xlsx")

d_out = read.xlsx("PhD Water Footprint/GLAM/EQ - new indicator/code/idlist_sdr.xlsx")#ok




#functions--------------------------------
results_model <- function(m, check){
  #this function calculates the predicted values of the model using a chack dataset and exports coefficients comparing real values and simulated values
  check_vals <- predict(m, check)
  
  results<- data.frame(
    BIC = BIC(m),#minimize
    kge = KGE(log(check$SR_tot),check_vals)[4],#the closest to 1 the better
    CorrCoef= KGE(log(check$SR_tot),check_vals)[1],
    BiasRatio= KGE(log(check$SR_tot),check_vals)[2],
    RelativeVariability= KGE(log(check$SR_tot),check_vals)[3],
    RootMeanSquareError = RMSE(log(check$SR_tot),check_vals),#minimize
    R2= R2(log(check$SR_tot),check_vals),#maximize
    MeanStandardError = MSE(log(check$SR_tot),check_vals)#,#
  )
  #extract model parameters
  coeffs <-  data.frame(t(m$coefficients))
  
  return(list(results= results,coeffs= coeffs))
}



results_sim_obs <- function(sim,obs){
  
  results_crosscheck<- data.frame(
  kge = KGE(obs,sim)[4],#the closest to 1 the better
  CorrCoef= KGE(obs,sim)[1],
  BiasRatio= KGE(obs,sim)[2],
  RelativeVariability= KGE(obs,sim)[3],
  RootMeanSquareError = RMSE(obs,sim),#minimize
  R2= R2(obs,sim),#maximize
  MeanStandardError = MSE(obs,sim))
  #print(results_crosscheck, digits =3)
  return(results_crosscheck)
}

best_model<-function(m){
  #dredging 
  options(na.action = 'na.fail')
  
  dd= dredge(m, beta="none", evaluate = 1, rank = 'BIC', 
             subset =! ("log(q)"&& "q"), fixed = "log(q)")#remove q keep only log(q)
  
  return(list(models = dd, best= get.models(dd, 1)[[1]]))
  
}

compare <- function(a,b){
  dif<-a-b
  c<- r.squared(a,b)
  print(plot(a,b)+  abline(a=0, b=1))
  paste0("r2 = ", c)
}



#data preparation---------------------------------------------------------------------

d_pcrglob <- na.omit(d_pcrglob)

#convert to factors - habitat realm ecoregion climate 5 clamate 30 idbasin pcrglob
d_pcrglob$ecoregion.f<-factor(d_pcrglob$ecoregion)
d_pcrglob$realm.f<-factor(d_pcrglob$realm)
d_pcrglob$habitat.f<-factor(d_pcrglob$habitat)
d_pcrglob$climate5.f<-factor(d_pcrglob$climate5)
d_pcrglob$id_basin_pcrglob.f<-factor(d_pcrglob$id_basin_pcrglob)#need to add that"""


#exclude num variables corresponding to factors
d_pcrglob.factors <- subset(d_pcrglob, select= c( 'realm.f', 'habitat.f', 'climate5.f', 'id_basin_pcrglob.f'))
d_pcrglob.vals <-subset(d_pcrglob, select= c('SR_tot', 'q', 'prec', 'temp', 'area', 'prec_delta','temp_delta', 'ti','slope','elevation'))
d_pcrglob.vals$logA <- log(d_pcrglob$area)
d_pcrglob.vals$logq <- log(d_pcrglob$q)
d_pcrglob.vals$logprec <- log(d_pcrglob$prec)
d_pcrglob.vals$native <- (d_pcrglob$SR_tot - d_pcrglob$SR_exo)/d_pcrglob$SR_tot*100
d_pcrglob_subset<-data.frame(d_pcrglob.vals,d_pcrglob.factors)




#sdr curve
#ggplot(subset(d_pcrglob_subset,(q>=0)&(area>2550)), aes(log(q),log(SR_tot)))+#best model
#  geom_point()+
#  theme_bw()+
#  theme(
#    panel.grid.major = element_blank(), 
#    panel.grid.minor = element_blank())+
  #geom_abline(slope= 1, intercept = 0)
  #ylim(0,10)+
  #xlim(0,10)+
  #facet_wrap(~ climate5.f)

#q distribution study
hist(log(d_pcrglob.vals$q,10))
plot(seq(0, 1, 0.1), quantile(log(d_pcrglob.vals$q,10),probs = seq(0, 1, 0.1)),  type = "b", xlab = c("percentile"))
quantile(d_pcrglob.vals$q,probs = seq(0, 1, 0.1))

#sr distribution study
hist(d_pcrglob.vals$SR_tot)
plot(seq(0, 1, 0.1), quantile(log(d_pcrglob.vals$SR_tot,10),probs = seq(0, 1, 0.1)),  type = "b", xlab = c("percentile"))
quantile(d_pcrglob.vals$SR_tot,probs = seq(0, 1, 0.1))

#native species
hist(d_pcrglob$SR_exo/d_pcrglob$SR_tot*100)
quantile(d_pcrglob.vals$native,probs = seq(0, 1, 0.1))

#climate and habitat distribution
hist(d_pcrglob$climate5)
hist(d_pcrglob$habitat, label=TRUE)

##test correlation between continuous vars
subset(d_pcrglob.vals, select= c('SR_tot', 'q', 'prec', 'temp', 'area', 'prec_delta','temp_delta', 'ti','slope','elevation'))%>%
         as.numeric(.)%>%
         cor(.,y=., method=c("spearman"), use="everything")


#correlation between variables (unprocessed)
cor(d_pcrglob.vals[-c(1)],d_pcrglob.vals[-c(1)], use = "everything", method = "spearman")

d_pcrglob.vals[-c(1)]%>%
  subset(.,(q>0)&(area>2550))%>%
  scale(., center=TRUE, scale=TRUE)%>%
  cor(.,., use = "complete.obs", method = "spearman")
        

table <- table(as.numeric(d_pcrglob.factors$realm.f), as.numeric(d_pcrglob.factors$habitat.f))
table
chisq.test(table)

table <- table(as.numeric(d_pcrglob.factors$climate5.f), as.numeric(d_pcrglob.factors$habitat.f))
table
chisq.test(table)


#remove correlated variables (with logQ)
d_pcrglob_subset<-subset(d_pcrglob_subset, select = -c(prec,logprec))


#preprocess dataset----------------------------------------------------------------

#basin filters 
d_pcrglob_subset<-subset(d_pcrglob_subset,(q>0)&(area>2550))

#figures of climate and habitat distributions
table(as.numeric(d_pcrglob_subset$habitat))
ggplot(d_pcrglob_subset, aes(d_pcrglob_subset$habitat))+
  geom_histogram(stat="count")+
  theme_bw()

table(as.numeric(d_pcrglob_subset$climate))
ggplot(d_pcrglob_subset, aes(d_pcrglob_subset$climate))+
  geom_histogram(stat="count")+
  theme_bw()

#center and scale if needed
d_pcrglob.vals.scale<-scale(subset(d_pcrglob_subset, select =-c(SR_tot,climate5.f, realm.f, habitat.f,id_basin_pcrglob.f)), center= TRUE, scale =TRUE)#not centered to avoid negative values
d_pcrglob.factors.scale <-subset(d_pcrglob_subset, select = c(climate5.f, realm.f, habitat.f, id_basin_pcrglob.f))
d_pcrglob_preproc<-data.frame(d_pcrglob.vals.scale,d_pcrglob.factors.scale)
d_pcrglob_preproc$SR_tot<-d_pcrglob_subset$SR_tot#we want to predice SR not centered and scaled SR

#run this line for non centered and non scale values!
d_pcrglob_preproc<-d_pcrglob_subset


##test correlation between continuous and categorical vars
kruskal.test(logq~climate5.f, data = d_pcrglob_preproc)#climate  explains significantly the variance in logq
kruskal.test(d_pcrglob_subset$logq~ d_pcrglob_subset$habitat.f)#habitat explains significantly the variance in loq
kruskal.test(temp~climate5.f, data = d_pcrglob_preproc)#climate explains significantly the variance in temperature
kruskal.test(elevation~climate5.f, data = d_pcrglob_preproc)
kruskal.test(temp~habitat.f, data = d_pcrglob_preproc)
kruskal.test(elevation~habitat.f, data = d_pcrglob_preproc)
kruskal.test(area~habitat.f, data =d_pcrglob_preproc)
kruskal.test(logq~habitat.f, data = d_pcrglob_preproc)
kruskal.test(prec_delta~habitat.f, data = d_pcrglob_preproc)
kruskal.test(temp~habitat.f, data =d_pcrglob_preproc)









#fit regression SDR with cross validation------------------------------------------------------------


#create 10 partitions for which the largest set is used for training i-e- 75% of the dataset
set.seed(5000)#first 2011, 2103, 5000


#K10 folds definitions---------------------------
N=10#number of splits

#stratified folds that maintain the major habitat distribution
folds<-createFolds(d_pcrglob_preproc$habitat.f, k=10, list = FALSE)

d<- d_pcrglob_preproc
d$folds<-folds

#M1 global SDR model-------------------------------------------------------------------
results <- vector("list",length = N)
coeffs <- vector("list",length = N)
model <-vector("list",length = N)


for (i in 1:N){
  #check<-d[-folds[,i],]
  #train<-d[folds[,i],]
  check <-subset(d, folds==i)
  train<-subset(d,folds=!i)
  
  #M1  = glm(SR_tot ~ logq,family=gaussian(link='log'), data=train)
  M1  = glm(log(SR_tot) ~ logq, family=gaussian(link='identity'), data=train)
  #M1  = glm(SR_tot ~ logq,family=poisson(link='log'), data=train)
  
  
  
  results <- rbind(results, results_model(M1,check)$results)
  coeffs<- rbind(coeffs,results_model(M1,check)$coeffs)
  model<-rbind(model, M1)#it does not work here
  
}


resultsM1<-results
modelM1<-model[-1]
bestM1<-modelM1[which.min(resultsM1$BIC)]

#find min BIC to select best model



#M2 global SDR+ ---------------------------------------------------------------


results <- vector("list",length = N)
coeffs <- vector("list",length = N)
model<-vector("list",length = N)

for (i in 1:N){
  
  #check<-d[-folds[,i],]
  #train<-d[folds[,i],]
  check <-subset(d, folds==i)
  train<-subset(d,folds=!i)
  check<-subset(check,select = c('SR_tot','logq','elevation', 'prec_delta', 'temp','temp_delta','slope','ti',"logA"))#exclude A and logA + categories
  train<-subset(train,select = c('SR_tot','logq','elevation', 'prec_delta', 'temp','temp_delta','slope','ti',"logA"))
  
  
  #M2 =glm(SR_tot ~ logq + ., family= gaussian(link = 'log'), data=train)
  M2 =glm(log(SR_tot) ~ logq + ., family= gaussian(link = 'identity'), data=train)
  #M2 =glm(SR_tot ~ logq + . ,family= poisson(link='log'), data=train)
  
  m<-M2
  
  options(na.action = 'na.fail')
  dd= dredge(m, beta="none", evaluate = 1, rank = 'BIC', 
             subset = with(+logq) && (!(('slope'&&'ti')|
                                        #('elevation'&&'slope')|
                                        #('elevation'&&'ti')|
                                          #('area'&&'logA')|
                                          ('temp'&&'temp_delta')
                                          #('temp'&&'prec_delta')
                                          #('temp'&&'ti'))
                                           ))
             )#remove q keep only log(q)

  bestM<- get.models(dd, 1)[[1]]

  results <- rbind(results,results_model(bestM,check)$results)
  model<-rbind(model,bestM)

}

resultsM2<-results
modelM2<-model[-1,]
bestM2<-modelM2[which.min(resultsM2$BIC),]

#M3 major habitat type model-----------------------------------------------------


results <- vector("list",length = N)
coeffs <- vector("list",length = N)
model <- vector("list",length = N)

for (i in 1:N){

  check <-subset(d, folds==i)
  train<-subset(d,folds=!i)
  
  #check<-subset(check,select = c('SR_tot','q','climate30.f','climate5.f','habitat.f','realm.f'))
  #train<-subset(train,select = c('SR_tot','q','climate30.f','climate5.f','habitat.f','realm.f'))
  check<-subset(check,select = c('SR_tot','q','logq','habitat.f','realm.f', "logA"))
  train<-subset(train,select = c('SR_tot','q','logq','habitat.f','realm.f',"logA"))

  
  #M3 <-glm(SR_tot ~ logq +., family = gaussian(link = 'log'), data=train)
  M3 <-glm(log(SR_tot) ~ logq +., family = gaussian(link = 'identity'), data=train)
  #M3<-glm(SR_tot ~ logq + . ,family= poisson(link='log'), data=train)
  m<- M3
  
  
  #when testing GLMs
  options(na.action = 'na.fail')
  dd= dredge(m, beta="none", evaluate = 1, rank = 'BIC', 
             subset =with(+logq) && (!(("logq"&& "q")|("habitat.f"&& "realm.f"))))#remove q keep only log(q)
  bestM<- get.models(dd, 1)[[1]]
  results <- rbind(results,results_model(bestM,check)$results)
  model<-rbind(model,bestM)

 
   #when testing lmer
  
  #M3 <-lmer(log(SR_tot) ~ logq + (1|habitat.f), REML=FALSE, data=train)
  #m<- M3
  #options(na.action = 'na.fail')
  #dd= dredge(m, beta="none", evaluate = 1, rank = 'BIC')
  #bestM<- get.models(dd, 1)[[1]]
  #results <- rbind(results,results_model(bestM, check))
  
  }

resultsM3<-results
modelM3<-model[-1,]
bestM3<-modelM3[which.min(resultsM3$BIC),]
  
#M4 climate and elevation model------------------------------------------------------------------------


results <- vector("list",length = N)
coeffs <- vector("list",length = N)
model <- vector("list",length = N)

for (i in 1:N){

  check <-subset(d, folds==i)
  train<-subset(d,folds=!i)
  
  #check<-subset(check,select = c('SR_tot','q','logq','elevation','climate5.f'))
  #train<-subset(train,select = c('SR_tot','q','logq','elevation','climate5.f'))
  check<-subset(check,select = c('SR_tot','q','logq','elevation','climate5.f','habitat.f', "logA"))
  train<-subset(train,select = c('SR_tot','q','logq','elevation','climate5.f','habitat.f', "logA"))

    #using GLM
 M4 =glm(log(SR_tot) ~ logq + .,  family = gaussian(link='identity'),data=train)
 #M4 =glm(SR_tot ~ logq + .,  family = gaussian(link='log'),data=train)
 #M4<-glm(SR_tot ~ logq + . ,family= poisson(link='log'), data=train)
 m<- M4
  
  options(na.action = 'na.fail')
  #when testing GLM use lines below
  dd= dredge(m, beta="none", evaluate = 1, rank = 'BIC', 
             subset =with(+logq) && (! (("logq"&& "q")|('climate5.f'&&'habitat.f'))))#remove q keep only log(q)
  bestM<- get.models(dd, 1)[[1]]
  results <- rbind(results,results_model(bestM,check)$results)
  model<-rbind(model,bestM)
  
  
  #using lmer
  #M4 <-lmer(log(SR_tot) ~ logq + elevation+ (1|climate5.f), REML=FALSE, data=train)
  #m<- M4
  #options(na.action = 'na.fail')
  #dd= dredge(m, beta="none", evaluate = 1, rank = 'BIC') 
  #bestM<- get.models(dd, 1)[[1]]
  #results <- rbind(results,results_model(bestM,check))
  
}

resultsM4<-results
modelM4<-model[-1,]
bestM4<-modelM4[which.min(resultsM4$BIC),]


#M5-------------------------------------------------------------------


#only ecoregions were removed
results <- vector("list",length = N)
coeffs <- vector("list",length = N)
model<- vector("list",length = N)

for (i in 1:N){

  check <-subset(d, folds==i)
  train<-subset(d,folds=!i)
  check<-subset(check,select = -c(id_basin_pcrglob.f, area, folds, native))
  train<-subset(train,select = -c(id_basin_pcrglob.f, area, folds, native))
  
  #remoidbasin
  
  #M5 =glm(SR_tot ~ logq + .,family = gaussian(link= 'log'), data=train)
  M5 =glm(log(SR_tot) ~ logq + .,family = gaussian(link= 'identity'), data=train)
  #M5<-glm(SR_tot ~ logq + .,family = poisson(link= 'log'), data=train)
  
  m<- M5
  
  options(na.action = 'na.fail')
  dd= dredge(m, beta="none", evaluate = 1, rank = 'BIC', #remove q keep only log(q)
             subset =with(+logq) && (! (("logq"&& "q")|
                                          #('elevation'&&'slope')|
                                          #('elevation'&&'ti')|
                                          ('slope'&&'ti')|
                                          #('area'&&'logA')|
                                          ('temp'&&'temp_delta')|
                                          #('prec_delta'&&'temp_delta')|
                                          #('prec_delta'&&'temp')|
                                          #('temp'&&'ti')|
                                          ('temp'&&'climate5.f')|
                                          ('elevation'&&'habitat.f')|
                                          ('climate5.f'&&'habitat.f')|
                                          ('habitat.f'&&'realm.f')|
                                          ('climate5.f'&&'realm.f'))
                                     )
             )
  
  
  bestM<- get.models(dd, 1)[[1]]
  results <- rbind(results,results_model(bestM,check)$results)
  model<-rbind(model,bestM)
}

resultsM5<-results
modelM5<-model[-1,]
bestM5<-modelM5[which.min(resultsM5$BIC),]




#models comparisons--------------------------------------------------------------



model_comparison_mean <- as.data.frame(
                  rbind(map_dbl(resultsM1,mean),
                          map_dbl(resultsM2,mean),
                          map_dbl(resultsM3,mean),
                          map_dbl(resultsM4,mean),
                          map_dbl(resultsM5,mean)))
row.names(model_comparison_mean)= c('M1','M2','M3','M4','M5')
print(model_comparison_mean, digits =3)
#the problem in MSE comes from the basins with SR=1



model_comparison_sd <- as.data.frame(
  rbind(map_dbl(resultsM1,sd),
        map_dbl(resultsM2,sd),
        map_dbl(resultsM3,sd),
        map_dbl(resultsM4,sd),
        map_dbl(resultsM5,sd)))

row.names(model_comparison_sd)= c('M1','M2','M3','M4','M5')
model_comparison_sd



#table of the equations resulting from the CV
model_equations<-as.data.frame(
  rbind(M1$formula,
        bestM2$formula,
        bestM3$formula,
        bestM4$formula,
        bestM5$formula))
row.names(model_equations)= c('M1','M2','M3','M4','M5')
model_equations



#best model selection---------------------------------------------------------------------

i<-which.min(resultsM5$BIC)
check<-subset(d, folds ==i)
train<-subset(d, folds !=i)
M5 =glm(bestM5$formula, family = gaussian(link = identity),data=train)
summary(M5)
coeffs5<-data.frame(CI=confint(M5, level = 0.95, digits = 3), estimate = data.frame(M5$coefficients))
print(coeffs5, digits = 2)

i<-which.min(resultsM2$BIC)
check<-subset(d, folds ==i)
train<-subset(d, folds !=i)
M2 =glm(bestM2$formula, family = gaussian(link = identity),data=train)
summary(M2)
coeffs2<-data.frame(CI=confint(M2, level = 0.95, digits = 3), estimate = data.frame(M2$coefficients))
print(coeffs2, digits = 2)


i<-which.min(resultsM4$BIC)
check<-subset(d, folds ==i)
train<-subset(d, folds !=i)
M4 =glm(bestM4$formula, family = gaussian(link = identity),data=train)
summary(M4)
coeffs4<-data.frame(CI=confint(M4, level = 0.95, digits = 2), estimate = data.frame(M4$coefficients))
print(coeffs4, digits = 2)


i<-which.min(resultsM1$BIC)
check<-subset(d, folds ==i)
train<-subset(d, folds !=i)
M1 =glm(bestM1$call, data=train)
summary(M1)
coeffs1<-data.frame(CI=confint(M1, level = 0.95, digits = 2), estimate = data.frame(M1$coefficients))
print(coeffs1, digits = 2)







#dominance analysis ----------------------------------------------------------
library(dominanceanalysis)

dominance <-dominanceAnalysis(M4)
print(dominance)
summary(dominance)


contributionByLevel(dominance)
averageContribution(dominance)






#verify GLM hypothesis for M5 ----------------------------------------------------------

best_M <- M5
summary(best_M)

#verify normal distribution of residuals, independence, linear regression, homogeneity
par(mfrow = c(2, 2))
plot(best_M, which = 1:4)
par(mfrow = c(1,1))

plot(M5$fitted.values, best_M$residuals)


sim50<-predict(best_M, d_pcrglob_preproc)#complete dataset
test50<-cbind(d_pcrglob_subset, sim50)






#verify GLM hypothesis for M4 ----------------------------------------------------------

best_M <- M4
summary(best_M)#not the right model


#verify normal distribution of residuals and homogeneity of variance
par(mfrow = c(2, 2))
plot(best_M, which = 1:4)
par(mfrow = c(1,1))


sim40<-predict(best_M, d_pcrglob_preproc)
obs40<-log(d_pcrglob_preproc$SR_tot)

test40<-cbind(d_pcrglob_preproc, sim40)

ggplot(test40, aes(log(SR_tot),sim40, colour = climate5.f))+#best model
  geom_point()+
  geom_abline(slope= 1, intercept = 0)
  





#sensitivity with Observed data-----------------------------

#this is the preprocessed dataset with python
d_sensitivity<- read.csv("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/output/dataset_SDR_GSIM_20211006.csv")#recalc
d_sensitivity <- d_sensitivity[,-1]
d_sensitivity$BAS <- d_sensitivity$id_basin_gsim

#join preprocessed data 
joint<- read.xlsx("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/output/joint.xlsx")
d_out <- merge(joint, d_sensitivity, by = 'BAS')#basin that has both pcrglob and gsim basin


#select latitude
gsim_meta<- read.csv("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/GSIM_catchment_characteristics.csv")
gsim_meta<-subset(gsim_meta, select=c('gsim.no', 'lat.new'))
gsim_meta$ID <- gsim_meta$gsim.no

#select basins with latitude +/- 42 (validity of xenopoulos SDR)
d_out <- merge(d_out, gsim_meta, by = 'ID')
d_out<- subset(d_out, (lat.new > -42.0000)&(lat.new < 42.0000))#n=362

#add species richness and observed discharge to preprocessed dataset
d<- read.xlsx("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/schipper_barbarossa_2021/input_tab_divAREA_ep.xlsx")
d<-subset(d,select= c('SR_tot', 'Q_MEAN', 'PREC_PRES', 'TEMP_PRES', 'AREA', 'PREC_DELTA','TEMP_DELTA', 'TI','SLOPE','ELEVATION','ID'))
d_out1 <- merge(d_out, d, by = 'ID')#basin that has both pcrglob and gsim basin


#prepare data
d <- d_out1#select valerios values 
#convert to factors - habitat realm ecoregion climate 5 clamate 30 idbasin pcrglob
d$ecoregion.f<-factor(d$ecoregion)
d$realm.f<-factor(d$realm)
d$habitat.f<-factor(d$habitat)
d$climate5.f<-factor(d$climate5)
d$climate30.f<-factor(d$climate30)
d$BAS.f<-factor(d$BAS)#need to add that"""

#exclude basins that are not within the validity of the fitted SDR from CV
d<-subset(d, (AREA > 2550)&(Q_MEAN>0))

#exclude num variables corresponding to factors
d.factors <- subset(d, select= c( 'realm.f', 'habitat.f', 'climate5.f', 'BAS.f'))
d.vals <-subset(d, select= c('SR_tot', 'Q_MEAN', 'prec', 'temp', 'area', 'prec_delta','temp_delta', 'ti','slope','elevation'))#my values
d.vals$logA <- log(d$area)
d.vals$logq <- log(d$Q_MEAN)

#data from Schipper and Barbarossa 2022
#d.vals <-subset(d, select= c('SR_tot', 'Q_MEAN', 'PREC_PRES', 'TEMP_PRES', 'AREA', 'PREC_DELTA','TEMP_DELTA', 'TI','SLOPE','ELEVATION'))#my values
#colnames(d.vals) <- c('SR_tot','q','prec', 'temp','area','prec_delta','temp_delta', 'ti', 'slope', 'elevation' )
#d.vals$logA <- log(d$AREA)
#d.vals$logq <- log(d$Q_MEAN)

d.vals.scale<-scale(subset(d.vals, select =-c(SR_tot)), center= FALSE, scale =FALSE)#select center and scale option like in for fitting the sdr above!
d.factors.scale <-subset(d.factors, select = -c(BAS.f))
d_preproc<-data.frame(d.vals.scale,d.factors.scale, d$SR_tot)


#sensitivity test

#M4
sim4<-predict(M4, d_preproc)
test4 = na.omit(cbind(d_preproc, sim4))

ggplot(test4, aes(log(d.SR_tot),sim4, colour = climate5.f))+#best model
  geom_point()+
  ylim(0,8)+
  xlim(0,8)+
  geom_abline(slope= 1, intercept = 0)
a<-results_sim_obs(test4$sim4, log(test4$d.SR_tot))

#M1
sim1<-predict(M1, d_preproc)
test1 = na.omit(cbind(d_preproc, sim1))
ggplot(test1, aes(log(d.SR_tot),sim1))+#best model
  geom_abline(slope= 1, intercept = 0)+
  geom_point()


e<-results_sim_obs(test1$sim1, log(test1$d.SR_tot))


#xenopoulos
sim0<-0.4*log(d$Q_MEAN)+4.2
test0 = na.omit(cbind(d_preproc, sim0))

ggplot(test0, aes(log(d.SR_tot),sim0))+#best model
  geom_point()+
  ylim(0,8)+
  xlim(0,8)+
  geom_abline(slope= 1, intercept = 0)

c<-results_sim_obs(test0$sim0, log(test0$d.SR_tot))

#table with results of sensitivity
results_cross_check <- rbind(a,e,c)
row.names(results_cross_check)= c('M4','M1','Xenopoulos')
print(results_cross_check, digits = 2)





#moran I test --------------------------------------------

#library(raster)

# define file path
shpFile <- "C:/Users/easpi/Documents/PhD Water Footprint/GLAM/EQ - new indicator/qgis/basin_5min_pcrglob_pseudo-mercator.shp"

# load shapefile
shp <- shapefile(shpFile)

# convert SpatialPointsDataFrame to data.frame
shpDF <- as.data.frame(shp)
shp_subset <- shpDF[shp$ID %in% d_pcrglob_preproc$id_basin_pcrglob,]

centroids <- as.data.frame(coordinates(shp))
centroids$ID <- shpDF$ID

#select rows where we have data
shp1<-shp[shp$ID %in% d_pcrglob_preproc$id_basin_pcrglob,]

centroids_subset <- centroids[centroids$ID %in% d_pcrglob_preproc$id_basin_pcrglob,]
str(centroids_subset)

plot(shp1)
points(centroids_subset[,c(1,2)], col='blue')

#calculate inverse of distance between catchments
dist_mat <- dist(centroids_subset[,1:2], method = 'euclidian', upper=1)
dist_mat <- dist_mat / 1000 # conversion in km
inv_dist_mat <- dist_mat^(-1)#element wise power
str(inv_dist_mat)
inv_dist_mat<-as.matrix(inv_dist_mat)


#perform test Moran I
library(ape)
#https://mgimond.github.io/simple_moransI_example/
best_M<-M4
sim<-predict(best_M, d_pcrglob_preproc)
obs<-log(d_pcrglob_preproc$SR_tot)
residuals <- obs-sim

MI <- Moran.I(residuals, inv_dist_mat, scaled= TRUE)#WEIGHTED BY INVERSE OF DISTANCE MATRIX
as.data.frame(MI)

library(spdep)
#other test with spdep
dist_mat1 <- dnearneigh(coordinates(shp1), 0, 1000000)#1000 km limit
lw<- nb2listw(dist_mat1, style = 'W', zero.policy = TRUE)#the list is complete, select only included polygons
#style W means that the links are summed per polygon

moran.test(residuals, lw, zero.policy=TRUE)#H0= random distribution of the residuals
moran.mc(residuals, lw, zero.policy=TRUE, nsim=599)
moran.plot(residuals, lw, zero.policy=TRUE)



#block validation for area---------------

library(splitTools)

N=10#number of splits

#select fold type
#folds<-create_folds(d_pcrglob_preproc$climate5.f,k=10, type=c("block"))
#folds<-create_folds(d_pcrglob_preproc$climate5.f,k=10, type=c("stratified"))
folds<-create_folds(d_pcrglob_preproc$area,k=10, type=c("stratified"))

d<- d_pcrglob_preproc


#model 4 training on stratified set
results <- vector("list",length = N)
coeffs <- vector("list",length = N)
model <- vector("list",length = N)

for (i in 1:N){
  check<-d[-folds[[i]],]
  train<-d[folds[[i]],]
  #check <-subset(d, folds==i)
  #train<-subset(d,folds=!i)
  
  check<-subset(check,select = c('SR_tot','q','logq','elevation','climate5.f','habitat.f',"logA"))
  train<-subset(train,select = c('SR_tot','q','logq','elevation','climate5.f','habitat.f',"logA"))
  
  #using GLM
  M4 =glm(log(SR_tot) ~ logq + .,  family = gaussian(link='identity'),data=train)
  #M4 =glm(SR_tot ~ logq + .,  family = gaussian(link='log'),data=train)
  #M4<-glm(SR_tot ~ logq + . ,family= poisson(link='log'), data=train)
  m<- M4

  options(na.action = 'na.fail')
  #when testing GLM use lines below
  dd= dredge(m, beta="none", evaluate = 1, rank = 'BIC', 
             subset =with(+logq) && (! (("logq"&& "q")|('climate5.f'&&'habitat.f'))))#remove q keep only log(q)
  bestM<- get.models(dd, 1)[[1]]
  results <- rbind(results,results_model(bestM,check)$results)
  model<-rbind(model,bestM)
  
}

resultsM4<-results
modelM4<-model[-1,]
bestM4<-modelM4[which.min(resultsM4$BIC),]

#show CV results
resultsM4%>%
  map_dbl(.,mean)%>%
  rbind(.)%>%
  as.data.frame(.)%>%
  arrange(.,BIC)%>%
  print(.)
#select best model
i<-which.min(resultsM4$BIC)
check<-d[-folds[[i]],]
train<-d[folds[[i]],]
M4bis =glm(bestM4$formula, family = gaussian(link = identity),data=train)
summary(M4bis)
coeffs4bis<-data.frame(CI=confint(M4bis, level = 0.95, digits = 2), estimate = data.frame(M4bis$coefficients))
#coefficients sensitivity
print(coeffs4bis, digits = 2)
#coefficients CV above
coeffs4

#export selected model------------------
#save results for not centered and not scaled results
write.xlsx(as.data.frame(coeffs4), "output/SDR_regression.xlsx", overwrite = FALSE, rowNames=TRUE, colNames=TRUE)



