################Covid-19 gene prevalance in Sweage###############
#                     Mohammad Nayeem Hasan                     #
#################################################################

library(MASS)
require(foreign)
require(ggplot2)
require(maptools)
library(tidyverse)
library(betareg)
library(car)
library(gapminder)
library(dplyr)

setwd("E:\\ResearchProject\\Aminul\\Sweage")

################Descriptive#############
library(ggpubr)
library(pastecs)
options(scipen = 999)
COVID <- read.csv("Data_corr.csv")

print(head(COVID, n=3))

stat.desc(COVID)

################Correlation#############
cor.test(COVID$Ph,COVID$Temperature)

ggplot(COVID, aes(x=Temperature, y=Ph, add = "reg.line")) + 
  geom_point()+ geom_smooth(method=lm) +  
  stat_cor(label.x = 24, label.y = 11.4) +
  stat_regline_equation(label.x = 24, label.y = 11.7) +
  labs(y="pH", x = "Temperature")


cor.test(COVID$Ph,COVID$Cy5.1.)
cor.test(COVID$Ph,COVID$FAM.1.)
cor.test(COVID$Ph,COVID$Rox.1.)

cor.test(COVID$Temperature,COVID$Cy5.1.)

ggplot(COVID, aes(x=Temperature, y=Cy5.1., add = "reg.line")) + 
  geom_point()+ geom_smooth(method=lm) +  
  stat_cor(label.x = 30, label.y = 26.2) +
  stat_regline_equation(label.x = 30, label.y = 26.7) +
  labs(y="Cy5(1)", x = "Temperature")

cor.test(COVID$Temperature,COVID$FAM.1.)

ggplot(COVID, aes(x=Temperature, y=FAM.1., add = "reg.line")) + 
  geom_point()+ geom_smooth(method=lm) +  
  stat_cor(label.x = 28, label.y = 26.2) +
  stat_regline_equation(label.x = 28, label.y = 26.7) +
  labs(y="FAM(1)", x = "Temperature")

cor.test(COVID$Temperature,COVID$Rox.1.)

ggplot(COVID, aes(x=Temperature, y=Rox.1., add = "reg.line")) + 
  geom_point()+ geom_smooth(method=lm) +  
  stat_cor(label.x = 30, label.y = 26.2) +
  stat_regline_equation(label.x = 30, label.y = 26.7) +
  labs(y="Rox(1)", x = "Temperature")


cor.test(COVID$Cy5.1.,COVID$FAM.1.)
cor.test(COVID$Cy5.1.,COVID$Rox.1.)

cor.test(COVID$FAM.1.,COVID$Rox.1.)



cor.test(COVID$Ph,COVID$Temperature)
cor.test(COVID$Ph,COVID$Cy5.2.)
cor.test(COVID$Ph,COVID$FAM.2.)
cor.test(COVID$Ph,COVID$Rox.2.)

cor.test(COVID$Temperature,COVID$Cy5.2.)
cor.test(COVID$Temperature,COVID$FAM.2.)
cor.test(COVID$Temperature,COVID$Rox.2.)

cor.test(COVID$Cy5.2.,COVID$FAM.2.)
cor.test(COVID$Cy5.2.,COVID$Rox.2.)

cor.test(COVID$FAM.2.,COVID$Rox.2.)

################Paired t-test#############
print(head(COVID, n=3))
t.test(COVID$Cy5.1., COVID$Cy5.2., paired = TRUE, alternative = "two.sided")
t.test(COVID$FAM.1., COVID$FAM.2., paired = TRUE, alternative = "two.sided")
t.test(COVID$Rox.1., COVID$Rox.2., paired = TRUE, alternative = "two.sided")

################Division wise gene#############

tiff("Box.tiff", units="in", width=6, height=4, res=300)

COVID <- read.csv("data_new_1.csv")

print(head(COVID, n=3))


p <- COVID %>% 
  filter(Gene %in% c("Cy5(1)","FAM(1)","Rox(1)")) %>%
  ggplot(aes(x=Division,y=ï..Ct, fill=Gene)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) +labs(y="Ct (number of cycles)", x = "Division") + theme(legend.position = "none",axis.text.y = element_text(colour="grey20",size=8,face="bold"))

p
COVID <- read.csv("data_new_2.csv")

print(head(COVID, n=3))


q <- COVID %>% 
  filter(Gene %in% c("IC","ORF1ab","N")) %>%
  ggplot(aes(x=Division,y=ï..Ct, fill=Gene)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2)+coord_flip()+labs(y="Ct (number of cycles)", x = "") + theme(axis.text.y = element_text(colour="grey20",size=8,face="bold"))
q
require(gridExtra)
gridExtra::grid.arrange(p,q, ncol=2, widths = 2:3)
dev.off()

COVID <- read.csv("data_new_pH.csv")

print(head(COVID, n=3))


COVID %>% 
  ggplot(aes(x=Division,y=ï..pH)) +
  geom_boxplot(fill = 'Gray', color = 'Black') + geom_jitter(width=0.1,alpha=0.2)+coord_flip()+labs(y="pH", x = "Division")

###########Station wise############
tiff("BoxD.tiff", units="in", width=6, height=4, res=300)
COVID <- read.csv("Dhaka_1.csv")

p <- ggplot(COVID, aes(x=ï..Station, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "Station") + theme(axis.title=element_text(size=15)
                                                                       ,axis.text.y = element_text(colour="grey20",size=12,face="bold")
                                                                       ,legend.title = element_text(color = "black", size = 14),
                                                                       legend.text = element_text(color = "black", size = 14))
p
COVID <- read.csv("Dhaka_2.csv")

q <- ggplot(COVID, aes(x=ï..Station, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "Station") + theme(axis.title=element_text(size=15)
                                                                       ,axis.text.y = element_text(colour="grey20",size=12,face="bold")
                                                                       ,legend.title = element_text(color = "black", size = 14),
                                                                       legend.text = element_text(color = "black", size = 14))
q
require(gridExtra)
gridExtra::grid.arrange(p,q)
dev.off()

#########Station wise Rohingya#######
tiff("BoxR.tiff", units="in", width=6, height=4, res=300)
COVID <- read.csv("Rohingya_1.csv")

p <- ggplot(COVID, aes(x=ï..Station, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "Station") + theme(axis.title=element_text(size=15)
                                                                       ,axis.text.y = element_text(colour="grey20",size=12,face="bold")
                                                                       ,legend.title = element_text(color = "black", size = 14),
                                                                       legend.text = element_text(color = "black", size = 14))
p
COVID <- read.csv("Rohingya_2.csv")

q <- ggplot(COVID, aes(x=ï..Station, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "Station") + theme(axis.title=element_text(size=15)
                                                                       ,axis.text.y = element_text(colour="grey20",size=12,face="bold")
                                                                       ,legend.title = element_text(color = "black", size = 14),
                                                                        legend.text = element_text(color = "black", size = 14))
q
require(gridExtra)
gridExtra::grid.arrange(p,q)
dev.off()
#########district wise Isolation#######
tiff("BoxI.tiff", units="in", width=6, height=4, res=300)
COVID <- read.csv("Isolation_1.csv")

p <- ggplot(COVID, aes(x=ï..Division, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "District") + theme(axis.text.y = element_text(colour="grey20",size=5,face="bold"))

p
COVID <- read.csv("Isolation_2.csv")

q <- ggplot(COVID, aes(x=ï..Division, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "District") + theme(axis.text.y = element_text(colour="grey20",size=5,face="bold"))
q
require(gridExtra)
gridExtra::grid.arrange(p,q)
dev.off()
#########district wise Sadar Hospital #######
tiff("BoxS.tiff", units="in", width=6, height=4, res=300)
COVID <- read.csv("Sadar_Hospital_1.csv")

p <- ggplot(COVID, aes(x=ï..Division, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "District") + theme(axis.text.y = element_text(colour="grey20",size=5,face="bold"))


COVID <- read.csv("Sadar_Hospital_2.csv")

q <- ggplot(COVID, aes(x=ï..Division, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "District") + theme(axis.text.y = element_text(colour="grey20",size=5,face="bold"))

require(gridExtra)
gridExtra::grid.arrange(p,q)
dev.off()
#########district wise Medical College Drain Wastage  #######
tiff("BoxC.tiff", units="in", width=6, height=4, res=300)
COVID <- read.csv("Community_Drain_1.csv")

p <- ggplot(COVID, aes(x=ï..Division, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "District") + theme(axis.text.y = element_text(colour="grey20",size=5,face="bold"))


COVID <- read.csv("Community_Drain_2.csv")

q <- ggplot(COVID, aes(x=ï..Division, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "District") + theme(axis.text.y = element_text(colour="grey20",size=5,face="bold"))

require(gridExtra)
gridExtra::grid.arrange(p,q)
dev.off()

#########district wise Medical College Drain Wastage  #######

COVID <- read.csv("City_Drain_1.csv")

ggplot(COVID, aes(x=ï..Division, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "District")


COVID <- read.csv("City_Drain_2.csv")

ggplot(COVID, aes(x=ï..Division, y=Ct, fill=Gene)) + 
  geom_bar(stat="identity", width=0.4, position=position_dodge(width=0.8)) +
  coord_flip() +labs(y="Ct (number of cycles)", x = "District")


##Heatmap
data <- read.csv("Heatmap.csv")
################heatmap that is working##################33
require(graphics); 
require(grDevices)

names(data)[names(data) == "ï..Confirmed_cases"] <- "Confirmed cases"


number <- cbind(data[,1], data[,2], data[,3], data[,4], data[,5], data[,6])
#biocLite("gplots")
x <- scale(number)
library("gplots")
hmcols<-colorRampPalette(c("Red","yellow"))(75)
heatmap.2(x,
          col=hmcols,
          scale="none",
          key=TRUE,
          symkey=FALSE,
          dendrogram="both",
          density.info="none",
          trace="none",
          cexRow= 0.6,
          cexCol = 0.8,
          labRow = rownames(data),
          labCol = colnames(data),
          margins=c(7, 5)
)

