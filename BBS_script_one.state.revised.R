rm(list=ls())
library('tidyverse')


read.bbs.dat.state<- function(fname,AOU){
  
  
  ## Read in full dataset with each state included ##
  
  bbs.dat<-read.csv(fname)  
  myStateNum = bbs.dat$StateNum[1]
  
  ## select only regular bbs data; i.e., no double observers, replication ##
  bbs.dat <- bbs.dat[bbs.dat$RPID%in%"101",]
  
  #filter out everything not speces of interest
  bbs.dat <- bbs.dat %>% filter(AOU==AOU.interest)
  
  ## Read in weather/observer data ##
  Weather<-read.csv("Weather/weather.csv")
  Weather <- Weather[Weather$RPID%in%"101",]
  
  ## Join both the counts and weather/observer data ##
  df <- left_join(Weather,bbs.dat,by=c("CountryNum"="CountryNum","StateNum"="StateNum","Route"="Route","Year"="Year"))
  
  ## Read in route data and join with the dataframe ##
  routes <- read.csv("routes.csv")
  df_2 <- left_join(df,routes,by=c("CountryNum"="CountryNum","StateNum"="StateNum","Route"="Route"))
  
  return(df_2)
  
}


## Specify Species of Interest ##
AOU.interest=6850

## Run each state/province through the function created ##
Alabama<- read.bbs.dat.state(fname="Alabama/Alabama.csv",AOU=AOU.interest)
Alaska<- read.bbs.dat.state(fname="Alaska/Alaska.csv",AOU=AOU.interest)
Alberta<- read.bbs.dat.state(fname="Alberta/Alberta.csv",AOU=AOU.interest)
Arizona<- read.bbs.dat.state(fname="Arizona/Arizona.csv",AOU=AOU.interest)
Arkansas<- read.bbs.dat.state(fname="Arkansa/Arkansa.csv",AOU=AOU.interest)
British_Colombia<- read.bbs.dat.state(fname="BritCol/BritCol.csv",AOU=AOU.interest)
California<- read.bbs.dat.state(fname="Califor/Califor.csv",AOU=AOU.interest)
Colorado<- read.bbs.dat.state(fname="Colorad/Colorad.csv",AOU=AOU.interest)
Connecticut<- read.bbs.dat.state(fname="Connect/Connect.csv",AOU=AOU.interest)
Delaware<- read.bbs.dat.state(fname="Delawar/Delawar.csv",AOU=AOU.interest)
Florida<- read.bbs.dat.state(fname="Florida/Florida.csv",AOU=AOU.interest)
Georgia<- read.bbs.dat.state(fname="Georgia/Georgia.csv",AOU=AOU.interest)
Idaho<- read.bbs.dat.state(fname="Idaho/Idaho.csv",AOU=AOU.interest)
Illinois<- read.bbs.dat.state(fname="Illinoi/Illinoi.csv",AOU=AOU.interest)
Indiana<- read.bbs.dat.state(fname="Indiana/Indiana.csv",AOU=AOU.interest)
Iowa<- read.bbs.dat.state(fname="Iowa/Iowa.csv",AOU=AOU.interest)
Kansas<- read.bbs.dat.state(fname="Kansas/Kansas.csv",AOU=AOU.interest)
Kentucky<- read.bbs.dat.state(fname="Kentuck/Kentuck.csv",AOU=AOU.interest)
Louisiana<- read.bbs.dat.state(fname="Louisia/Louisia.csv",AOU=AOU.interest)
Maine<- read.bbs.dat.state(fname="Maine/Maine.csv",AOU=AOU.interest)
Manitoba<- read.bbs.dat.state(fname="Manitob/Manitob.csv",AOU=AOU.interest)
Maryland<- read.bbs.dat.state(fname="Marylan/Marylan.csv",AOU=AOU.interest)
Massachusetts<- read.bbs.dat.state(fname="Massach/Massach.csv",AOU=AOU.interest)
Michigan<- read.bbs.dat.state(fname="Michiga/Michiga.csv",AOU=AOU.interest)
Minnesota<- read.bbs.dat.state(fname="Minneso/Minneso.csv",AOU=AOU.interest)
Mississippi<- read.bbs.dat.state(fname="Mississ/Mississ.csv",AOU=AOU.interest)
Missouri<- read.bbs.dat.state(fname="Missour/Missour.csv",AOU=AOU.interest)
Montana<- read.bbs.dat.state(fname="Montana/Montana.csv",AOU=AOU.interest)
New_Brunswick<- read.bbs.dat.state(fname="NBrunsw/NBrunsw.csv",AOU=AOU.interest)
North_Carolina<- read.bbs.dat.state(fname="NCaroli/NCaroli.csv",AOU=AOU.interest)
North_Dakota<- read.bbs.dat.state(fname="NDakota/NDakota.csv",AOU=AOU.interest)
Nebraska<- read.bbs.dat.state(fname="Nebrask/Nebrask.csv",AOU=AOU.interest)
Nevada<- read.bbs.dat.state(fname="Nevada/Nevada.csv",AOU=AOU.interest)
Newfoundland<- read.bbs.dat.state(fname="Newfoun/Newfoun.csv",AOU=AOU.interest)
New_Hampshire<- read.bbs.dat.state(fname="NHampsh/NHampsh.csv",AOU=AOU.interest)
New_Jersey<- read.bbs.dat.state(fname="NJersey/NJersey.csv",AOU=AOU.interest)
New_Mexico<- read.bbs.dat.state(fname="NMexico/NMexico.csv",AOU=AOU.interest)
Nova_Scotia<- read.bbs.dat.state(fname="NovaSco/NovaSco.csv",AOU=AOU.interest)
Nunavut<- read.bbs.dat.state(fname="Nunavut/Nunavut.csv",AOU=AOU.interest)
NWTerri<- read.bbs.dat.state(fname="NWTerri/NWTerri.csv",AOU=AOU.interest)
New_York<- read.bbs.dat.state(fname="NYork/NYork.csv",AOU=AOU.interest)
Ohio<- read.bbs.dat.state(fname="Ohio/Ohio.csv",AOU=AOU.interest)
Oklahoma<- read.bbs.dat.state(fname="Oklahom/Oklahom.csv",AOU=AOU.interest)
Ontario<- read.bbs.dat.state(fname="Ontario/Ontario.csv",AOU=AOU.interest)
Oregon<- read.bbs.dat.state(fname="Oregon/Oregon.csv",AOU=AOU.interest)
PEI<- read.bbs.dat.state(fname="PEI/PEI.csv",AOU=AOU.interest)
Pennsylvania<- read.bbs.dat.state(fname="Pennsyl/Pennsyl.csv",AOU=AOU.interest)
Quebec<- read.bbs.dat.state(fname="Quebec/Quebec.csv",AOU=AOU.interest)
Rhode_Island<- read.bbs.dat.state(fname="RhodeIs/RhodeIs.csv",AOU=AOU.interest)
Saskatchewan<- read.bbs.dat.state(fname="Saskatc/Saskatc.csv",AOU=AOU.interest)
South_Carolina<- read.bbs.dat.state(fname="SCaroli/SCaroli.csv",AOU=AOU.interest)
South_Dakota<- read.bbs.dat.state(fname="SDakota/SDakota.csv",AOU=AOU.interest)
Tennessee<- read.bbs.dat.state(fname="Tenness/Tenness.csv",AOU=AOU.interest)
Texas<- read.bbs.dat.state(fname="Texas/Texas.csv",AOU=AOU.interest)
Utah<- read.bbs.dat.state(fname="Utah/Utah.csv",AOU=AOU.interest)
Vermont<- read.bbs.dat.state(fname="Vermont/Vermont.csv",AOU=AOU.interest)
Virginia<- read.bbs.dat.state(fname="Virgini/Virgini.csv",AOU=AOU.interest)
West_Virginia<- read.bbs.dat.state(fname="W_Virgi/W_Virgi.csv",AOU=AOU.interest)
Washington<- read.bbs.dat.state(fname="Washing/Washing.csv",AOU=AOU.interest)
Wisconsin<- read.bbs.dat.state(fname="Wiscons/Wiscons.csv",AOU=AOU.interest)
Wyoming<- read.bbs.dat.state(fname="Wyoming/Wyoming.csv",AOU=AOU.interest)
Yukon<- read.bbs.dat.state(fname="Yukon/Yukon.csv",AOU=AOU.interest)


Alabama$AOU<-AOU.interest
Alabama$SpeciesTotal[is.na(Alabama$SpeciesTotal)]<-0
Alaska$AOU<-AOU.interest
Alaska$SpeciesTotal[is.na(Alaska$SpeciesTotal)]<-0
Alberta$AOU<-AOU.interest
Alberta$SpeciesTotal[is.na(Alberta$SpeciesTotal)]<-0
Arizona$AOU<-AOU.interest
Arizona$SpeciesTotal[is.na(Arizona$SpeciesTotal)]<-0
Arkansas$AOU<-AOU.interest
Arkansas$SpeciesTotal[is.na(Arkansas$SpeciesTotal)]<-0
British_Colombia$AOU<-AOU.interest
British_Colombia$SpeciesTotal[is.na(British_Colombia$SpeciesTotal)]<-0
California$AOU<-AOU.interest
California$SpeciesTotal[is.na(California$SpeciesTotal)]<-0
Colorado$AOU<-AOU.interest
Colorado$SpeciesTotal[is.na(Colorado$SpeciesTotal)]<-0
Connecticut$AOU<-AOU.interest
Connecticut$SpeciesTotal[is.na(Connecticut$SpeciesTotal)]<-0
Delaware$AOU<-AOU.interest
Delaware$SpeciesTotal[is.na(Delaware$SpeciesTotal)]<-0
Florida$AOU<-AOU.interest
Florida$SpeciesTotal[is.na(Florida$SpeciesTotal)]<-0
Georgia$AOU<-AOU.interest
Georgia$SpeciesTotal[is.na(Georgia$SpeciesTotal)]<-0
Idaho$AOU<-AOU.interest
Idaho$SpeciesTotal[is.na(Idaho$SpeciesTotal)]<-0
Illinois$AOU<-AOU.interest
Illinois$SpeciesTotal[is.na(Illinois$SpeciesTotal)]<-0
Indiana$AOU<-AOU.interest
Indiana$SpeciesTotal[is.na(Indiana$SpeciesTotal)]<-0
Iowa$AOU<-AOU.interest
Iowa$SpeciesTotal[is.na(Iowa$SpeciesTotal)]<-0
Kansas$AOU<-AOU.interest
Kansas$SpeciesTotal[is.na(Kansas$SpeciesTotal)]<-0
Kentucky$AOU<-AOU.interest
Kentucky$SpeciesTotal[is.na(Kentucky$SpeciesTotal)]<-0
Louisiana$AOU<-AOU.interest
Louisiana$SpeciesTotal[is.na(Louisiana$SpeciesTotal)]<-0
Maine$AOU<-AOU.interest
Maine$SpeciesTotal[is.na(Maine$SpeciesTotal)]<-0
Manitoba$AOU<-AOU.interest
Manitoba$SpeciesTotal[is.na(Manitoba$SpeciesTotal)]<-0
Maryland$AOU<-AOU.interest
Maryland$SpeciesTotal[is.na(Maryland$SpeciesTotal)]<-0
Massachusetts$AOU<-AOU.interest
Massachusetts$SpeciesTotal[is.na(Massachusetts$SpeciesTotal)]<-0
Michigan$AOU<-AOU.interest
Michigan$SpeciesTotal[is.na(Michigan$SpeciesTotal)]<-0
Minnesota$AOU<-AOU.interest
Minnesota$SpeciesTotal[is.na(Minnesota$SpeciesTotal)]<-0
Mississippi$AOU<-AOU.interest
Mississippi$SpeciesTotal[is.na(Mississippi$SpeciesTotal)]<-0
Missouri$AOU<-AOU.interest
Missouri$SpeciesTotal[is.na(Missouri$SpeciesTotal)]<-0
Montana$AOU<-AOU.interest
Montana$SpeciesTotal[is.na(Montana$SpeciesTotal)]<-0
Nebraska$AOU<-AOU.interest
Nebraska$SpeciesTotal[is.na(Nebraska$SpeciesTotal)]<-0
Nevada$AOU<-AOU.interest
Nevada$SpeciesTotal[is.na(Nevada$SpeciesTotal)]<-0
New_Brunswick$AOU<-AOU.interest
New_Brunswick$SpeciesTotal[is.na(New_Brunswick$SpeciesTotal)]<-0
New_Hampshire$AOU<-AOU.interest
New_Hampshire$SpeciesTotal[is.na(New_Hampshire$SpeciesTotal)]<-0
New_Jersey$AOU<-AOU.interest
New_Jersey$SpeciesTotal[is.na(New_Jersey$SpeciesTotal)]<-0
New_Mexico$AOU<-AOU.interest
New_Mexico$SpeciesTotal[is.na(New_Mexico$SpeciesTotal)]<-0
New_York$AOU<-AOU.interest
New_York$SpeciesTotal[is.na(New_York$SpeciesTotal)]<-0
Newfoundland$AOU<-AOU.interest
Newfoundland$SpeciesTotal[is.na(Newfoundland$SpeciesTotal)]<-0
North_Carolina$AOU<-AOU.interest
North_Carolina$SpeciesTotal[is.na(North_Carolina$SpeciesTotal)]<-0
North_Dakota$AOU<-AOU.interest
North_Dakota$SpeciesTotal[is.na(North_Dakota$SpeciesTotal)]<-0
Nova_Scotia$AOU<-AOU.interest
Nova_Scotia$SpeciesTotal[is.na(Nova_Scotia$SpeciesTotal)]<-0
Nunavut$AOU<-AOU.interest
Nunavut$SpeciesTotal[is.na(Nunavut$SpeciesTotal)]<-0
NWTerri$AOU<-AOU.interest
NWTerri$SpeciesTotal[is.na(NWTerri$SpeciesTotal)]<-0
Ohio$AOU<-AOU.interest
Ohio$SpeciesTotal[is.na(Ohio$SpeciesTotal)]<-0
Oklahoma$AOU<-AOU.interest
Oklahoma$SpeciesTotal[is.na(Oklahoma$SpeciesTotal)]<-0
Ontario$AOU<-AOU.interest
Ontario$SpeciesTotal[is.na(Ontario$SpeciesTotal)]<-0
Oregon$AOU<-AOU.interest
Oregon$SpeciesTotal[is.na(Oregon$SpeciesTotal)]<-0
PEI$AOU<-AOU.interest
PEI$SpeciesTotal[is.na(PEI$SpeciesTotal)]<-0
Pennsylvania$AOU<-AOU.interest
Pennsylvania$SpeciesTotal[is.na(Pennsylvania$SpeciesTotal)]<-0
Quebec$AOU<-AOU.interest
Quebec$SpeciesTotal[is.na(Quebec$SpeciesTotal)]<-0
Rhode_Island$AOU<-AOU.interest
Rhode_Island$SpeciesTotal[is.na(Rhode_Island$SpeciesTotal)]<-0
Saskatchewan$AOU<-AOU.interest
Saskatchewan$SpeciesTotal[is.na(Saskatchewan$SpeciesTotal)]<-0
South_Carolina$AOU<-AOU.interest
South_Carolina$SpeciesTotal[is.na(South_Carolina$SpeciesTotal)]<-0
South_Dakota$AOU<-AOU.interest
South_Dakota$SpeciesTotal[is.na(South_Dakota$SpeciesTotal)]<-0
Tennessee$AOU<-AOU.interest
Tennessee$SpeciesTotal[is.na(Tennessee$SpeciesTotal)]<-0
Texas$AOU<-AOU.interest
Texas$SpeciesTotal[is.na(Texas$SpeciesTotal)]<-0
Utah$AOU<-AOU.interest
Utah$SpeciesTotal[is.na(Utah$SpeciesTotal)]<-0
Vermont$AOU<-AOU.interest
Vermont$SpeciesTotal[is.na(Vermont$SpeciesTotal)]<-0
Virginia$AOU<-AOU.interest
Virginia$SpeciesTotal[is.na(Virginia$SpeciesTotal)]<-0
Washington$AOU<-AOU.interest
Washington$SpeciesTotal[is.na(Washington$SpeciesTotal)]<-0
West_Virginia$AOU<-AOU.interest
West_Virginia$SpeciesTotal[is.na(West_Virginia$SpeciesTotal)]<-0
Wisconsin$AOU<-AOU.interest
Wisconsin$SpeciesTotal[is.na(Wisconsin$SpeciesTotal)]<-0
Wyoming$AOU<-AOU.interest
Wyoming$SpeciesTotal[is.na(Wyoming$SpeciesTotal)]<-0
Yukon$AOU<-AOU.interest
Yukon$SpeciesTotal[is.na(Yukon$SpeciesTotal)]<-0

## bind all states together ##
bbs.dat=rbind(Alabama,Alaska,Alberta,Arizona,Arkansas,British_Colombia,California,Colorado,Connecticut,Delaware,Florida,Georgia,Idaho,Illinois,Indiana,Iowa,Kansas,Kentucky,Louisiana,Maine,Manitoba,Maryland,Massachusetts,Michigan,Minnesota,Mississippi,Missouri,Montana,New_Brunswick,North_Carolina,North_Dakota,Nebraska,Nevada,Newfoundland,New_Hampshire,New_Jersey,New_Mexico,Nova_Scotia,Nunavut,NWTerri,New_York,Ohio,Oklahoma,Ontario,Oregon,PEI,Pennsylvania,Quebec,Rhode_Island,Saskatchewan,South_Carolina,South_Dakota,Tennessee,Texas,Utah,Vermont,Virginia,West_Virginia,Washington,Wisconsin,Wyoming,Yukon) 


## Clear environment of individual states ##
rm(Alabama,Alaska,Alberta,Arizona,Arkansas,British_Colombia,California,Colorado,Connecticut,Delaware,Florida,Georgia,Idaho,Illinois,Indiana,Iowa,Kansas,Kentucky,Louisiana,Maine,Manitoba,Maryland,Massachusetts,Michigan,Minnesota,Mississippi,Missouri,Montana,New_Brunswick,North_Carolina,North_Dakota,Nebraska,Nevada,Newfoundland,New_Hampshire,New_Jersey,New_Mexico,Nova_Scotia,Nunavut,NWTerri,New_York,Ohio,Oklahoma,Ontario,Oregon,PEI,Pennsylvania,Quebec,Rhode_Island,Saskatchewan,South_Carolina,South_Dakota,Tennessee,Texas,Utah,Vermont,Virginia,West_Virginia,Washington,Wisconsin,Wyoming,Yukon) 



## Group by BCR (Only includes BCRs with total counts > 100 for the given time range) ##

bbs.dat.summarized <- bbs.dat %>%  
  group_by(BCR) %>% 
  summarize(n_routes=n(),SpeciesTotal=sum(SpeciesTotal,na.rm=TRUE))%>%
  filter(SpeciesTotal > 100)

bcr.with.species<-bbs.dat.summarized$BCR
bcr.with.species

## Once relevant BCRs identified, filter data to include only those BCRs of interest ##
bbs.dat <- bbs.dat %>% filter(BCR %in% bcr.with.species)


## Save Data for Species to New CSV File ##
#write.csv(bbs.dat,file="BBS_(Species)_state.bind.csv")
