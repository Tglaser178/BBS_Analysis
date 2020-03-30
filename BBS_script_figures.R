library(tidyverse)

# Read in Filtered Data for Species of Interest
#bbs.dat<- read.csv("BBS_(Species)_BCR.filtered.csv")
bbs.dat<- read.csv("BBS_WIWA_BCR.filtered.csv")


# Create new dataframe that summarizes the data by BCR and Year, and produces a new column that averages the number of birds for each route

df_Species <- bbs.dat %>% 
  group_by(BCR,Year) %>% 
  summarize(n_routes=n(),SpeciesTotal=sum(SpeciesTotal,na.rm=TRUE))%>%
mutate(bpr=SpeciesTotal/n_routes)

# Convert BCR variable to a factor for figures
df_Species$BCR<-as_factor(df_Species$BCR)


#Figures
library(ggplot2)

# Average Abundance of Counts (for species of interest) per each BCR
ggplot(df_Species, aes(x=BCR,y=bpr)) +
  geom_bar(stat="identity") +
  labs(x="BCR",y="Average Abundance per Route (1967-2018)")+
  ggtitle("Average Abundance by Bird Conservation Region")

# Average Abundance of Counts (for species of interest) per year for each BCR
ggplot(df_Species, aes(x=Year,y=bpr,group=BCR,color=BCR))+
  geom_line()+
  scale_x_discrete(limits=c(1970,1980,1990,2000,2010,2018),labels=c("1970","1980","1990","2000","2010","2018"))+
  labs(x="Year",y="Average Abundance per Route (1967-2018)")+
  ggtitle("Average Abundance per BCR Route vs. Year")




