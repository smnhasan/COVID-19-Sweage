library(tidyverse)
library(rstatix)
library(ggpubr)
library(extrafont)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(lme4)

setwd("E:\\ResearchProject\\Aminul\\Sweage")


nayeem.data <- read.csv("Boxplot2.csv", 
                        header = TRUE)
nayeem.data$CNLog100 <- nayeem.data$CNLog*100
nayeem.data$Copy.Number
nayeem.data$nTemp <- as.numeric(nayeem.data$Temp)
nayeem.data$nPH <- as.numeric(nayeem.data$PH)
data <- nayeem.data[(nayeem.data$Gene=="ORF1ab"),]

options(scipen = 999)
nayeem.data$VLLog
library(glmmTMB)
fit <- glmmTMB(VLLog ~ Site + Region  + nTemp + ( 1| Month),  na.action=na.omit,  data = data)

summary(fit)

data <- nayeem.data[(nayeem.data$Gene=="N"),]

library(glmmTMB)
fit <- glmmTMB(VLLog ~ Site + Region  + nTemp + ( 1| Month),  na.action=na.omit,  data = data)

summary(fit)

library(RColorBrewer)

################Gene*Month=Site#############

p <- nayeem.data %>% 
  filter(Site %in% c("Drainage","Without drainage")) %>%
  ggplot(aes(x=Gene,y=Ct.Value, fill=Site)) + facet_wrap(~Month) + 
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Ct (number of cycles)", x = "Gene") +  scale_fill_brewer(palette="Accent")+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p



p <- nayeem.data %>% 
  filter(Site %in% c("Drainage","Without drainage")) %>%
  ggplot(aes(x=Region,y=Ct.Value, fill=Site)) + facet_wrap(~Gene) + 
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Ct (number of cycles)", x = "Sectors") +  scale_fill_brewer(palette="Dark2")+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p

nayeem.data$VLLog
p <- nayeem.data %>% 
  filter(Site %in% c("Drainage","Without drainage")) %>%
  ggplot(aes(x=Region,y=VLLog, fill=Site)) + facet_wrap(~Gene) + 
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Viral load (Log10 genome copies/100 mL)", x = "Sectors") +  scale_fill_brewer(palette="Set1")+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p


nayeem.data$CNLog
p <- nayeem.data %>% 
  filter(Site %in% c("Drainage","Without drainage")) %>%
  ggplot(aes(x=Region,y=CNLog, fill=Site)) +facet_wrap(~Gene) + 
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Copy Number (Log10 genome copies/100 mL)", x = "Sectors") +  scale_fill_brewer(palette="Set2")+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p

########1a.	Disparity urban (DRAINAGE) and rural (WITHOUT DRAINAGE)###########

nayeem.data.market <- nayeem.data[nayeem.data$Station=="Market Place" | nayeem.data$Station=="Bus Stand/ Rail station " | nayeem.data$Station=="Communal Pond/River",]

p <- nayeem.data.market %>% 
  filter(Site %in% c("Drainage","Without drainage")) %>%
  ggplot(aes(x=Region,y=Ct.Value, fill=Gene)) + facet_wrap(~Station) + 
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Ct (number of cycles)", x = "Sectors")  +  scale_fill_brewer(palette="Paired")+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p

library(wesanderson)
nayeem.data.bus_stand <- nayeem.data[nayeem.data$Station=="Bus Stand/ Rail station ",]

p <- nayeem.data %>% 
  filter(Station %in% c("Market Place","Isolation Center Drain Effluent")) %>%
  ggplot(aes(x=Region,y=Ct.Value, fill=Station)) + facet_wrap(~Gene) + 
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Ct (number of cycles)", x = "Sectors") +  scale_fill_manual(values = wes_palette("BottleRocket1", n = 3))+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p

nayeem.data.river <- nayeem.data[nayeem.data$Station=="Communal Pond/River",]

p <- nayeem.data.river %>% 
  filter(Site %in% c("Drainage","Without drainage")) %>%
  ggplot(aes(x=Region,y=Ct.Value, fill=Gene)) +
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Ct (number of cycles)", x = "Sectors") +  scale_fill_manual(values = wes_palette("BottleRocket2", n = 3))+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p




nayeem.data.iso <- nayeem.data[nayeem.data$Station=="Isolation Center Drain Effluent" | nayeem.data$Station=="Medical College Drain Wastage" | nayeem.data$Station=="City Drain Effluent",]

p <- nayeem.data.iso %>% 
  filter(Site %in% c("Drainage","Without drainage")) %>%
  ggplot(aes(x=Region,y=Ct.Value, fill=Gene)) + facet_wrap(~Station) + 
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Ct (number of cycles)", x = "Sectors") +  scale_fill_manual(values = wes_palette("Royal1", n = 3))+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p


nayeem.data.hos <- nayeem.data[nayeem.data$Station=="Sadar Hospital Drain Effluent" ,]

p <- nayeem.data.hos %>% 
  filter(Site %in% c("Drainage","Without drainage")) %>%
  ggplot(aes(x=Region,y=Ct.Value, fill=Gene)) +
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Ct (number of cycles)", x = "Sectors") +  scale_fill_manual(values = wes_palette("Darjeeling2", n = 3))+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p

nayeem.data.city <- nayeem.data[nayeem.data$Station=="City Drain Effluent",]

p <- nayeem.data.city %>% 
  filter(Site %in% c("Drainage","Without drainage")) %>%
  ggplot(aes(x=Region,y=Ct.Value, fill=Gene)) +
  geom_boxplot() + geom_jitter(alpha=0.1) +labs(y="Ct (number of cycles)", x = "Sectors") +  scale_fill_manual(values = wes_palette("Royal2", n = 3))+
  theme(axis.title=element_text(size=15), legend.position="bottom"
        ,axis.text.y = element_text(colour="grey20",size=12)
        ,axis.text.x = element_text(colour="grey20",size=12)
        ,legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))
p <- p 
p

library(psych)
###Mean Difference###
#Mymensing Sylhet Rajshahi Rangpur Gaibandha Dhaka Barisal Chittagong  Khulna    
#Cox's Bazar Comilla Kishorgonj Brahmobaria Habigonj
#ORF1ab N
#December January
nayeem.data$Region
#data <- nayeem.data[(nayeem.data$Region=="Padma"),]
data <- nayeem.data[(nayeem.data$Gene=="N"),]
data_i <- data[(data$Station=="Isolation Center Drain Effluent"),]
data_m <- data[(data$Station=="Market Place"),]
describe(data_i$Ct.Value)
describe(data_m$Ct.Value)
t.test(data_i$Ct.Value[1:28], data_m$Ct.Value, paired = T, var.eq = F)


#Mymensing Sylhet Rajshahi Rangpur Gaibandha Dhaka Barisal Chittagong  Khulna    
#Cox's Bazar Comilla Kishorgonj Brahmobaria Habigonj

data <- nayeem.data[(nayeem.data$=="Brahmobaria"),]
data <- data[(data$Gene=="N"),]
data <- data[(data$Month=="January"),]
data <- data[(data$Station=="Isolation Center Drain Effluent" | data$Station=="Sadar Hospital Drain Effluent" | data$Station=="City Drain Effluent"),] 
data
describe(data$Ct.Value)


nayeem.data$