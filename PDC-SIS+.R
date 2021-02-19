gt <- read.csv("E:/Deep learning/Log GTCA/SIS-PDC+/NewGTlog10.csv")

gtm<- read.csv("E:/Deep learning/Log GTCA/gt/gt matachecd/GTstationarylog10Matched.csv")
gdp <- df[,"GDP"]

data <- df[,3:83]
library(dplyr)
## apply Lag1 to Lag3
PredictLag3 <- as.data.frame(data) %>%
  mutate_all(list(lag1 = ~lag(.), lag2 = ~lag(., 2))) %>%
  select(gtools:: mixedorder(names(.)))
Xwl <- PredictLag3[4:201,]

## alternative way to lag3
#library(data.table)
#PredictLag3 <- as.data.frame(data) %>%
  #mutate_all(list(lag1 = ~lag(.), lag2 = ~lag(., 2), lag3 = ~lag(., 3))) %>%
 #select(gtools:: mixedorder(names(.)))
#PredictLag3 <- setDT(data)[, paste0(rep(names(data), each = 2), "Lag", 1:3) := shift(.SD, n = 1:3)]
#setcolorder(data, order(sub("Lag\\d+", "", names(data))))

#source('funs.r')
library(energy)
library(tsDyn)
library(mgcv)
library(mnormt)
library(tseries)

## sample size
n = 201

## total number of series
d=  81
## lags to consider as variables for each series
r = 3

## used for burn in
n2 = n
library(Hmisc)

## create response and its 3 lags
Y0 <- as.data.frame(gdp)
Lag1 <- as.data.frame(Lag(gdp,1)) 
Lag2 <- as.data.frame(Lag(gdp,2))
Lag3 <- as.data.frame(Lag(gdp,3))

# Remove the missing values
Y0 = Y0[4:201,]
Y1 = Lag1[4:201,]
Y2 = Lag2[4:201,]
Y3 = Lag3[4:201,]


###############################

###distance correlation screening: DC-SIS
dc =rep(0,d*r)
for(j in 1:(d*r)){
  dc[j] = dcor(Y0,Xwl[,j])
}
dc.sort <- sort(abs(dc), method= "sh", index=TRUE, decreasing=TRUE)

dc.ind <- dc.sort$ix


## PDC-SIS and threshold
thresh1= ceil(sqrt()) 
thresh = 24


##method 1: for each lag, condition on the previous lag for the variable in consideration
## three lags of response are always part of the conidtioning set in

pdc = rep(0,d*r)
for(j in 1:d){
  pdc[3*j-2] =  pdcor(Y0, Xwl[,3*j-2], cbind(Y1,Y2,Y3))
  pdc[3*j-1] =  pdcor(Y0, Xwl[,3*j-1], cbind(Xwl[,3*j-2], Y1,Y2,Y3))
  pdc[3*j] = pdcor(Y0, Xwl[,3*j], cbind(Xwl[,3*j-2], Xwl[,3*j-1], Y1,Y2,Y3))
}

pdc.sort <- sort(abs(pdc), method= "sh", index=TRUE, decreasing=TRUE)
pdc.ind <- pdc.sort$ix
indices = pdc.ind[1:thresh]
screened.set.sis = Xwl[,indices] 

write.csv(screened.set.sis,file = "gt24.csv")

### PDC-SIS+
## Use block bootstrap to get threshold for PDC-SIS
k=5
null=c()

##BB <- as.data.frame(null)

for(i in 1:k){
  Yfull=tsbootstrap(Y0,nb=30)
  Ynew=Yfull[4:(n+3)]
  Ynew1=Yfull[3:(n+2)]
  Ynew2=Yfull[2:(n+1)]
  Ynew3=Yfull[1:n]
  for(j in 1:d){
    null[j+(i-1)*d]=pdcor(Ynew, Xwl[,3*j-2], cbind(Ynew1,Ynew2,Ynew3))
  }
}

Condset=cbind(Y1,Y2,Y3)
allind=NULL

## choose threshold as 99th percentile 
pdc.rind <- quantile(null,.99)

## choose threshold  
thresh = 24


i=1
pdci=c()
while(i <= r) {##loop through each lag
  pdc2 = rep(0,d)
  if(i==1)
  {
    for(j in 1:d){
      pdci[3*j-3+i] =  pdcor(Y0, Xwl[,3*j-3+i], Condset)
      pdc2[j]=pdci[3*j-3+i]
      
    }
  }
  if(i==2) {
    for(j in 1:d){
      pdci[3*j-3+i] =  pdcor(Y0, Xwl[,3*j-3+i], cbind(Condset,Xwl[,3*j-3+i+1]))
      pdc2[j]=pdci[3*j-3+i]
      
    }
  }
  if(i==3) {
    for(j in 1:d){
      pdci[3*j-3+i] =  pdcor(Y0, Xwl[,3*j-3+i], cbind(Condset,Xwl[,3*j-3+i-1],Xwl[,3*j-3+i-2]))
      pdc2[j]=pdci[3*j-3+i]
      
    }
  }
  indi = which(abs(pdc2)>max(abs(pdc.rind)))
  pdc.sort2 <- sort(abs(pdc2), method= "sh", index=TRUE, decreasing=TRUE)
  ## choose a cap of sqrt(n) for how many variables we can add to our conditioning set
  indi1=pdc.sort2$ix[1:16]
  indi2=intersect(indi,indi1)
  allind = c(allind, 3*(indi2)-3+i)
  Condset = cbind(Condset, Xwl[,3*(indi2)-3+i])
  i=i+1
}

pdc.sort <- sort(abs(pdci), method= "sh", index=TRUE, decreasing=TRUE)
pdc.ind <- pdc.sort$ix
indices = pdc.ind[1:thresh]
screened.set.sis1 = Xwl[,indices] 
data1 = as.data.frame(screened.set.sis1)

write.csv(screened.set.sis1,file = "gt_pdc_24.csv")
###########################################

