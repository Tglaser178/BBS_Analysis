##############################################################################
# Analysis of  bbs data to get trends by BCR
# Created by T. Glaser April 09 2020
# Based on code sent by Jim Sarroco for WIWA BBS analysis (wiwa_network_bbs_model1.R)
# which  was in turn based on code from J. Sauer
##############################################################################
rm(list=ls())

library("ff") # to read in large BBS data set
library("rgeos")
library("rgdal")
library("maptools")
library("sp")
library("jagsUI")
library('tidyverse')

SPECIES="Wilson's Warbler"
YEAR=c(1967:2018)
OUTFILE="data/outBBS/BBS_WIWA_BCR.RData"



out<-list()
#-----------------------------------------------------------------------------
# 1 - Read in and process BBS data. 
#-----------------------------------------------------------------------------

### read in filtered data for species
bbs.dat <- read.csv("./data/BBS_WIWA_BCR.filtered.csv")
#Create unique values for filtered BCRs
BCR=unique(bbs.dat$BCR)

### Create unique route field 
bbs.dat$Route <- factor(do.call(paste, list(bbs.dat$BCR, bbs.dat$Route, sep=".")))  
# select years of interest
#bbs.dat <- bbs.dat[bbs.dat$Year>(YEAR[1]-1) & bbs.dat$Year<(YEAR[2]+1),]

# Save number of routes to "out" list
out$nRoutes<-length(unique(bbs.dat$Route))


### Observer data
bbs.dat$Rte <- factor(do.call(paste, list(bbs.dat$State, bbs.dat$Route, sep=".")))
# Create Observer ID based on route x observer
bbs.dat$obs <- factor(do.call(paste, list(bbs.dat$Obs, bbs.dat$Rte, sep="."))) 
# create first year indicator variable
bbs.dat$fy <- ifelse(duplicated(bbs.dat$ObsN), 0, 1) 


# fill in zeros for years where route counted with no observations of species
bbs.dat$SpeciesTotal[is.na(bbs.dat$SpeciesTotal)] <- 0 

# Extract Coordinates
bbs.locs <- bbs.dat[,c("Longitude","Latitude")] 

#-----------------------------------------------------------------------------
# 2 - Read in and process spatial data, match to BBS data
#-----------------------------------------------------------------------------

### Read in regions = BCR shapefile
REGIONS <- readOGR(dsn="./data/bcr_terrestrial_shape",layer="BCR_Terrestrial_master_International")
#regions.proj=proj4string(REGIONS)
#project to Albers equal area so areas calculated will be in sq meters
REGIONS<-spTransform(REGIONS,CRS("+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5"))

# Apply area to each filtered BCR
bbs.BCR=BCR
BCR=data.frame(BCR=bbs.BCR, AreaBCR=sapply(1:length(BCR),function(i){gArea(subset(REGIONS,BCR==BCR[i]),byid=T)}))

# Merge BCR areas with bbs dataframe
BCR$BCR=as.factor(bbs.BCR)
AreaBCR.df=as.data.frame(tapply(BCR$AreaBCR,BCR$BCR,sum,simplify=T))
names(AreaBCR.df)<-c("AreaBCR")
AreaBCR.df$AreaBCR<-row.names(AreaBCR.df)
BCR=merge(BCR,AreaBCR.df,by.x="BCR")
bbs.dat=merge(bbs.dat,BCR,by.x="BCR")
bbs.dat$BCR<-factor(bbs.dat$BCR,levels=as.character(unique(bbs.dat$BCR)))  # remove strata from which ther are no observations


#-----------------------------------------------------------------------------
#  Write BBS model to file. Code provided by J. Sauer. (via Jim Sarocco)
#-----------------------------------------------------------------------------

source("LinkSauerModel.R")

#-----------------------------------------------------------------------------
# 4 - prep/bundle data for bbs trend model, set inits, run model
#-----------------------------------------------------------------------------

### Set Parameters for BBS Model

count <- bbs.dat$SpeciesTotal
ncounts <- length(bbs.dat$SpeciesTotal)
obser <- bbs.dat$ObsN
obser <- as.numeric(factor(bbs.dat$ObsN))
nobservers <- length(unique(obser))   
firstyr <- bbs.dat$fy
year <- as.numeric(factor(bbs.dat$Year))
nyears <- length(unique(year))
strat <- as.numeric(factor(bbs.dat$BCR))    
nstrata <- length(unique(strat))

# Set area weight for BCRs
aw <- unique(subset(bbs.dat, select = c(BCR, AreaBCR)))
aw <- aw[order(aw$BCR),]
areaweight <- aw$AreaBCR

# calculate z weights 
rte.all <- read.csv("./data/routes.csv") 
rte.all$Rte <- factor(do.call(paste, list(rte.all$BCR, rte.all$Route, sep=".")))
#rte.all <- rte.all[rte.all$Year>(YEAR[1]-1) & rte.all$Year<(YEAR[2]+1),]

rte.all=rte.all[sapply(rte.all$BCR,function(x){any(x==BCR$BCR)}),]
BCR.df=data.frame(BCR=as.factor(bbs.BCR))
rte.all=merge(rte.all,BCR.df,by="BCR")   

Rte <- read.csv("./data/routes.csv") 
rte.sum <- aggregate(Rte~BCR, rte.all, length) # use BCR scale here
names(rte.sum)[2] <- "tot.rtes"
spec.rte.sum <- aggregate(Rte~BCR, unique(subset(bbs.dat, select=c(Rte, BCR))), length)
names(spec.rte.sum)[2] <- "detec.rtes"
wts <- merge(spec.rte.sum, rte.sum)
wts$nonzeroweight <- wts$detec.rtes/wts$tot.rtes
wts<-wts[order(wts$BCR),]
nonzeroweight <- wts$nonzeroweight


### bundle data:
nYears<-length(YEAR)-YEAR[1]
jags.data <- list(count=count, year = year, obser=obser,  nyears=nyears, firstyr=firstyr, ncounts=ncounts, strat=strat,  
                  nobservers=nobservers, nstrata=nstrata, areaweight=areaweight, nonzeroweight=nonzeroweight, fixedyear=round(nYears/2))

### bundle data:
nYears<-YEAR[2]-YEAR[1]
jags.data <- list(count=count, year = year, obser=obser,  nyears=nyears, firstyr=firstyr, ncounts=ncounts, strat=strat,  
                  nobservers=nobservers, nstrata=nstrata, areaweight=areaweight, nonzeroweight=nonzeroweight, fixedyear=round(nYears/2))

# Initial values
inits <- function(){
  list(tauyear=rep(1,nstrata),taunoise=1,tauobs=1,beta=rep(0,nstrata),strata=rep(0,nstrata),eta=0)
}  

# Parameters monitored
parameters <- c("eta", "N", "sdnoise", "sdobs", "CompIndex", "Bbar")

# Set MCMC to desired settings
ni <- 50
nt <- 3
nb <- 10
nc <- 1

print("Calling JAGS")
bbs.out <- jags(jags.data, inits, parameters, "bbs_model_13.txt",n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
bbs.BCR=levels(bbs.dat$BCR)
bbs.years=YEAR
bbs.areas<-as.data.frame(aw,row.names=1:nrow(aw))
save(bbs.out,bbs.BCR,bbs.areas,bbs.years,file=OUTFILE)


# years<-seq(1967,2018)
# reg.counts<-matrix(NA,length(years),11)
# reg.ci.025<-matrix(NA,length(years),11)
# reg.ci.975<-matrix(NA,length(years),11)
# for(i in 1:length(years)) reg.counts[i,]<-apply(bbs.out$sims.list$n[,,i],2,mean)
# for(i in 1:length(years)) reg.ci.025[i,]<-apply(bbs.out$sims.list$n[,,i],2,quantile, probs=c(0.025))
# for(i in 1:length(years)) reg.ci.975[i,]<-apply(bbs.out$sims.list$n[,,i],2,quantile, probs=c(0.975))



# reg.trends <- 100*((bbs.out$sims.list$n[,,17]/bbs.out$sims.list$n[,,1])^(1/16)-1)
# reg.tr.mn <- apply(reg.trends, 2, mean)
# reg.tr.ci <- apply(reg.trends, 2, quantile, probs=c(0.025, 0.975))
# 
# 
# reg.tr.df <- data.frame(Reg = levels(bbs.dat$Reg), mn = reg.tr.mn, t(reg.tr.ci))
