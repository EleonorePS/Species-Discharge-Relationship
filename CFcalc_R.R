
#import libraries-------------------------


library(ggplot2)
library(stats)
library(openxlsx)
library(purrr)
rm(list=ls())



windowsFonts(Font=windowsFont("Arial"))

#set working directory------------------

getwd()
setwd("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/")
#PUT INFUT FILES HERE


#import files --------------------------------------

#this is the preprocessed dataset with python
d_pcrglob<- read.xlsx("output/dataset_SDR_no_filter_20220316_climatecorrected.xlsx")

#this is the selected modelfrom the cross validation
#sdr<- read.xlsx("output/SDR_regression_2022_07_06.xlsx")#ok
sdr<- read.xlsx("output/SDR_regression_2022_07_06_newm4.xlsx")#ok

#global extinction probabilities
gep<-read.xlsx("GEP_fish_all.xlsx")#does not change

#depletion factors
dfdt<- read.xlsx("depletion_factors_smooth_all_mix-1960-2000_20220316.xlsx")




#-PATCH  FUNCTION------------------
patch<-function(d,vals,col){
  d$isna <- is.na(d[col])
  
  s<-as.data.frame(vals[[col]])
  
  f<-which((d$isna==TRUE)&(is.na(d$climate5)))
  d[f,col]<- median(na.omit(d[[col]]))#median where there are values
  
  f<-which((d$isna==TRUE)&(d$climate5 == 0))
  d[f,col]<- median(na.omit(d[[col]]))
  
  for (i in 1:5){
    f<-which((d$isna==TRUE)&(d$climate5 == i)) 
    if (length(f)==0) next
    d[f,col]<- s[i,2]
    next
  }
  colnames(d)[length(colnames(d))]<-paste0('isna',col,collapse="")
  return(d)
}
#clean data-------------------------------
d_pcrglob_subset<-subset(d_pcrglob, area>2550 & q>0)#-20318 - 13618 basins 
d_pcrglob_subset<-na.omit(d_pcrglob_subset)
str(d_pcrglob_subset)



#MARGINAL EFFECT FACTOR --------------

a<-sdr[8,]#logq  coeffs not scaled
a

#marginal EF
#EF_marginal<-subset(d_pcrglob_subset,select = c('q','id_basin_pcrglob'))
EF_marginal_m<-a[['M4.coefficients']]/d_pcrglob_subset$q/(3600*24*365.25)# in PDF.m^-3.yr, warning! the model must fit non centered non scaled data 
EF_marginal2.5<-a[['CI.2.5..']]/d_pcrglob_subset$q/(3600*24*365.25)# in PDF.m^-3.yr, warning! the model must be fit n the not cetered and ot scaled vars
EF_marginal97.5<-a[['CI.97.5..']]/d_pcrglob_subset$q/(3600*24*365.25)# in PDF.m^-3.yr, warning! the model must be fit n the not cetered and ot scaled vars

EF_marginal <-data.frame(EF_marginal_m,EF_marginal2.5,EF_marginal97.5)
EF_marginal$id_basin_pcrglob <- d_pcrglob_subset$id_basin_pcrglob
EF_marginal$climate5 <- d_pcrglob_subset$climate5
EF_marginal<- na.omit(EF_marginal)

#distribution table
table<-map(EF_marginal[,1:3],aggregate,by = EF_marginal['climate5'],FUN=median)
table<-data.frame(table$EF_marginal_m$x, table$EF_marginal2.5$x, table$EF_marginal97.5$x)
colnames(table)<- c('EF marginal estimate', 'CI2.5','CI97.5')
row.names(table)<- c('Tropical','Dry','Temperate','Continental', 'Polar')
print(table, digits = 3)

table1<-map(EF_marginal[,1:3],quantile,c(0,0.25,0.5,0.75,1))
print(data.frame(table1),digits = 3)
data.frame(table1[1])


#violin plot---------
lab<- c("Tropical","Arid","Temperate","Cold","Polar")
ggplot(EF_marginal, aes(climate5, log(EF_marginal_m,10), fill = as.factor(climate5)))+#best model
  geom_violin()+
  #geom_boxplot(outlier.size = 0.1, varwidth = 1)+
  labs(fill = "Climate")+
  xlab("")+
  ylab('log10 EFm')+
  scale_fill_manual(#values = c("Tropical" = "dark blue", "Arid" = "orange", "Temperate" = "darkgreen", "Cold" = "lightblue", "Polar"="deeppink"),
    values = c("1" = "dark blue", "2" = "orange", "3" = "forestgreen", "4" = "darkslategray1", "5"="deeppink"),
    labels = c("Tropical","Arid","Temperate","Cold","Polar"))+
  ylim(c(-15,0))+
  theme_bw()+
  theme(
    text = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
  )


ggsave("output/violinplot_marginal.png", width = 4, height = 4)



#AVERAGE EFFECT FACTORS-----------------------------------------------

a<-sdr[8,]#logq  coeffs not scaled

EF_average<-subset(d_pcrglob_subset,select = c('qnat','q','climate5','id_basin_pcrglob'))
colnames(EF_average)<- c('qnatural','qhuman','climate5','id_basin_pcrglob')
EF_average<-subset(EF_average, qnatural>qhuman)#loose 3592-2785
EF_average$dif<-EF_average$qnatural-EF_average$qhuman


EF_average$EF_average_m<-(1-(EF_average$qhuman/EF_average$qnatural)^a[['M4.coefficients']])/((EF_average$qnatural-EF_average$qhuman)*365.25*3600*24)
EF_average$EF_average2.5<-(1-(EF_average$qhuman/EF_average$qnatural)^a[['CI.2.5..']])/((EF_average$qnatural - EF_average$qhuman)*365.25*3600*24)
EF_average$EF_average97.5<-(1-(EF_average$qhuman/EF_average$qnatural)^a[['CI.97.5..']])/((EF_average$qnatural - EF_average$qhuman)*365.25*3600*24)
EF_average <- na.omit(EF_average)#remove no data. rests 2785 
str(EF_average)


#table
table<-map(EF_average[,6:8],quantile,c(0,0.25,0.5,0.75,1))
print(data.frame(table),digits = 3)
data.frame(table[1])

#violin plot
ggplot(EF_average, aes(climate5, log(EF_average_m,10), fill = as.factor(climate5)))+#best model
  geom_violin()+
  #geom_boxplot(outlier.size = 0.1, varwidth = 1)+
  labs(fill = "Climate")+
  xlab("")+
  ylab('log10 EFa')+
  
  scale_fill_manual(#values = c("Tropical" = "dark blue", "Arid" = "orange", "Temperate" = "darkgreen", "Cold" = "lightblue", "Polar"="deeppink"),
    values = c("1" = "dark blue", "2" = "orange", "3" = "forestgreen", "4" = "darkslategray1", "5"="deeppink"),
    labels = c("Tropical","Arid","Temperate","Cold","Polar"))+
  ylim(c(-15,-0))+
  theme_bw()+
  theme(
    text = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    
  )

ggsave("output/EQ - new indicator/results/violinplot_avg.png", width = 4, height =4)







#compare marginal and average EF-------------
EF<-merge(EF_marginal, EF_average,by='id_basin_pcrglob')

test=subset(EF,EF$EF_marginal_m<=EF$EF_average_m)
test$ratio<-test$qhuman/test$qnatural
str(test)#2744lines for marginal < average and 34 toherwise

EF$ratio <-EF$EF_marginal_m/EF$EF_average_m

ggplot(EF, aes(qhuman,EF_marginal_m))+
  geom_point()+
  #geom_point(EF_average_m)+
  scale_y_log10()

ggplot(EF, aes(qhuman,EF_average_m))+
  geom_point()+
  scale_y_log10()


print(quantile(EF$ratio, seq(0,1,0.1),na.rm = TRUE), digits = 3)
mean(EF$ratio)




#PATCHING EFFECT FACTORS ----------------------------------------------------------------



df_merge_mgl<-merge.data.frame(subset(d_pcrglob, select =c('id_basin_pcrglob','climate5')),EF_marginal[,1:4],by='id_basin_pcrglob', all = TRUE)


#calc average values per climate zones
s<-map(subset(EF_marginal, select = c('EF_marginal_m', 'EF_marginal2.5', 'EF_marginal97.5')), aggregate, EF_marginal['climate5'], FUN=max)
str(s)#insertvalues where there are NA

df_m_patch <- patch(df_merge_mgl,s, 'EF_marginal_m')
df_m_patch <- patch(df_m_patch,s, 'EF_marginal2.5')
df_m_patch <- patch(df_m_patch,s, 'EF_marginal97.5')


s1<-map(subset(EF_average, select = c('EF_average_m', 'EF_average2.5', 'EF_average97.5')), aggregate, EF_average['climate5'], FUN=max)
str(s1)#insertvalues where there are NA

df_merge_avg<-merge.data.frame(subset(d_pcrglob, select =c('id_basin_pcrglob','climate5')),EF_average[,c(4,6,7,8)],by = 'id_basin_pcrglob', all.x = TRUE)#only basins with data
df_a_patch <- patch(df_merge_avg,s1, 'EF_average_m')
df_a_patch <- patch(df_a_patch,s1,'EF_average2.5')
df_a_patch <- patch(df_a_patch,s1,'EF_average97.5')


#EXPORT EFFECT FACTORS-------------------------------------------------
ex1<-subset(df_m_patch, select = c('id_basin_pcrglob','EF_marginal_m', 'EF_marginal2.5','EF_marginal97.5','isnaEF_marginal_m'))
colnames(ex1)[5]<-c('isna_marginal')

ex2<-subset(df_a_patch, select = c('id_basin_pcrglob','EF_average_m', 'EF_average2.5','EF_average97.5','isnaEF_average_m','climate5'))
colnames(ex2)[5]<-c('isna_average')
export = merge(ex1,ex2,by=c('id_basin_pcrglob'))

write.xlsx(export, "output/EF_patch_07082022_newm4.xlsx", overwrite = TRUE)#ok

#EF[EF$id_basin_pcrglob==17313,"EF_marginal_m"]



#PATCH FATE FACTORS-----------------------------------------------


dfdt<- read.csv("C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code/output/depletion_factors_smooth_all_mix-1960-2000_20220708.csv")


df_merge<-merge.data.frame(subset(d_pcrglob, select =c('id_basin_pcrglob','climate5')),dfdt,by = 'id_basin_pcrglob', all.x = TRUE)#only basins with data
dfdt<-merge.data.frame(subset(d_pcrglob, select =c('id_basin_pcrglob','climate5')),dfdt,by = 'id_basin_pcrglob', all.y = TRUE)


s0<-map(na.omit(subset(dfdt, exclude = c('id_basin_pcrglob','climate5'))), 
       aggregate, 
       dfdt['climate5'],
       FUN=median)#median value


dfdt_patch<-patch(df_merge,s0,'DF_q')
dfdt_patch<-patch(dfdt_patch,s0,'DF_gws')
#dfdt_patch<-patch(dfdt_patch,s0,'DT_gws')
dfdt_patch<-patch(dfdt_patch,s0,'DF_et')
#dfdt_patch<-patch(dfdt_patch,s0,'DT_et')
dfdt_patch<-patch(dfdt_patch,s0,'DF_sm')
#dfdt_patch<-patch(dfdt_patch,s0,'DT_sm')

write.xlsx(dfdt_patch, "output/DFDT_patch_20220708.xlsx",overwrite = TRUE)#ok'

dfdt_patch<- read.xlsx("output/DFDT_patch_20220708.xlsx")






#CALCULATE CHARACTERIZATION FACTORS------------------


FF<-subset(dfdt_patch,select=c('id_basin_pcrglob', 'DF_q','isnaDF_q'))
colnames(FF)[3]<-c('isna_FF')#patched value

FF$FF<--FF$DF_q
for (i in seq(nrow(FF))) if (FF[i,'DF_q']<0) FF[i,'FF']=-FF[i,'DF_q'] else FF[i,'FF']="No discharge depletion"
#No discharge depletion = (FF<0)
nrow(subset(FF, (FF=="No discharge depletion")&(isna_FF==FALSE)))/nrow(subset(FF, isna_FF==FALSE))
#8% of river basins

FF<-subset(FF, select = c('id_basin_pcrglob', 'FF',"isna_FF"))
FF$FF<-as.numeric(FF$FF)

EF<-export#patched effect factors file

CF<-merge(FF, EF, by='id_basin_pcrglob')
CF<-merge(CF, gep, by.x='id_basin_pcrglob',by.y="basins_5min_pcrglobwb_adjusted",all.x = TRUE)

CF$CF_marginal_m<-CF$FF*CF$EF_marginal_m#NA means that the FF is <0
CF$CF_marginal2.5<-CF$FF*CF$EF_marginal2.5
CF$CF_marginal97.5<-CF$FF*CF$EF_marginal97.5

CF$CF_average_m<-CF$FF*CF$EF_average_m#NA means that the FF is <0
CF$CF_average2.5<-CF$FF*CF$EF_average2.5
CF$CF_average97.5<-CF$FF*CF$EF_average97.5

CF$CF_gep_marginal_m<-CF$CF_marginal_m*CF$GEP_freshwater_fish#NA means that the FF is <0
CF$CF_gep_marginal2.5<-CF$FF*CF$EF_marginal2.5*CF$GEP_freshwater_fish
CF$CF_gep_marginal97.5<-CF$FF*CF$EF_marginal97.5*CF$GEP_freshwater_fish

CF$CF_gep_average_m<-CF$FF*CF$EF_average_m*CF$GEP_freshwater_fish#NA means that the FF is <0
CF$CF_gep_average2.5<-CF$FF*CF$EF_average2.5*CF$GEP_freshwater_fish
CF$CF_gep_average97.5<-CF$FF*CF$EF_average97.5*CF$GEP_freshwater_fish

CF[CF$id_basin_pcrglob==17313,"EF_marginal_m"]

#DISTRIBUTION
#marginal
CF1=na.omit(CF)
CF1=CF1[(CF1$isna_FF==FALSE)&(CF1$isna_marginal==FALSE),]
table<-quantile(CF1$CF_marginal_m,c(0.0,0.25,0.5,0.75,1), na.rm = TRUE)
print(data.frame(table),digits = 3)

CF1=CF1[CF1$CF_gep_marginal_m!=0,]
table<-quantile(CF1$CF_gep_marginal_m,c(0.0,0.25,0.5,0.75,1), na.rm = TRUE)
print(data.frame(table),digits = 3)





#average
CF1=na.omit(CF)
CF1=CF1[(CF1$isna_FF==FALSE)&(CF1$isna_average==FALSE),]
table<-quantile(CF1$CF_average_m,c(0.0,0.25,0.5,0.75,1), na.rm = TRUE)
print(data.frame(table),digits = 3)

CF1=CF1[(CF1$isna_FF==FALSE)&(CF1$isna_average==FALSE)&(CF1$CF_gep_average_m!=0),]
table<-quantile(CF1$CF_gep_average_m,c(0.0,0.25,0.5,0.75,1), na.rm = TRUE)
print(data.frame(table),digits = 3)

CF2=CF[(CF$isna_FF==FALSE)&(CF$isna_average==FALSE),]
table<-quantile(CF2$CF_average_m,c(0.0,0.25,0.5,0.75,1), na.rm = TRUE)
print(data.frame(table),digits = 3)



#all values

table<-quantile(CF$CF_marginal_m,c(0.0,0.25,0.5,0.75,1), na.rm = TRUE)
print(data.frame(table),digits = 3)




#EXPORT CHARACERIZATION FACTORS -----------------------------


write.xlsx(CF, "output/CF_patch_newm4_20220708_forstats.xlsx", keepNA=TRUE, overwrite = TRUE)#ok'

#tables SI



tab=rbind(CF[CF$id_basin_pcrglob==(17313),"EF_marginal_m"],
          CF[CF$id_basin_pcrglob==(22098),"EF_marginal_m"],
          CF[CF$id_basin_pcrglob==(16558),"EF_marginal_m"],
          CF[CF$id_basin_pcrglob==(14990),"EF_marginal_m"],
          CF[CF$id_basin_pcrglob==(24212),"EF_marginal_m"],
          CF[CF$id_basin_pcrglob==(20745),"EF_marginal_m"])
row.names(tab)=c("Ganges river",
                 "Nile river",
                 "Yangtze river",
                 "Euphrates river",
                 "Amazon river",
                 "Niger river")
colnames(tab)=c('EF_marginal_m')
print(tab)

tab1=rbind(CF[CF$id_basin_pcrglob==(17313),c("CF_marginal_m","CF_gep_marginal_m",'FF')],
          CF[CF$id_basin_pcrglob==(19039),c("CF_marginal_m","CF_gep_marginal_m",'FF')],
          CF[CF$id_basin_pcrglob==(13884),c("CF_marginal_m","CF_gep_marginal_m",'FF')],
          CF[CF$id_basin_pcrglob==(17644),c("CF_marginal_m","CF_gep_marginal_m",'FF')],
          CF[CF$id_basin_pcrglob==(14842),c("CF_marginal_m","CF_gep_marginal_m",'FF')],
          CF[CF$id_basin_pcrglob==(14842),c("CF_marginal_m","CF_gep_marginal_m",'FF')])
tab1
row.names(tab1)=c("Ganges river",
                 "Godavari river",
                 "Yellow river",
                 "Pearl river",
                 "Red river",
                 "Arkansas river")
colnames(tab1)=c('CF_marginal_regional','CF_marginal_global','FF')
print(tab1)


