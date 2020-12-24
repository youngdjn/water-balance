library(ggplot2)
library(viridis)
library(gridExtra)
library(cowplot)
library(Cairo)
library(dplyr)

# Scenario: aspect

d <- read.csv("data/aspect.csv",header=TRUE,stringsAsFactors=FALSE)
d$aspect <- reorder(d$aspect,c(1,2,3,1,2,3))
d$Crop.coefficient <- reorder(d$Crop.coefficient,c(2,2,2,1,1,1))

# Plot PET
pa <- ggplot(d,aes(Crop.coefficient,PET)) +
  geom_bar(aes(fill=aspect),position="dodge",stat="identity") +
  scale_y_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  geom_hline(yintercept=50,colour="blue",linetype=2,size=0.5) +
  theme(text = element_text(size=20)) +
  scale_x_discrete(drop=FALSE) +
  xlab("") +
  ylab("Water quantity (e.g., mm)") +
  ggtitle("PET") +
  scale_fill_viridis(name="Aspect",discrete=TRUE, begin=0.05, end=0.95) +
  theme_bw(8) +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = 0.5))

# Plot AET
pb <- ggplot(d,aes(Crop.coefficient,AET)) +
  geom_bar(aes(fill=aspect),position="dodge",stat="identity") +
  scale_y_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  geom_hline(yintercept=50,colour="blue",linetype=2,size=0.5) +
  theme(text = element_text(size=20)) +
  scale_x_discrete(drop=FALSE) +
  xlab("PET coefficient") +
  ylab("") +
  ggtitle("AET") +
  scale_fill_viridis(name="Aspect",discrete=TRUE, begin=0.05, end=0.95) +
  theme_bw(8) +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = 0.5))

# Plot Deficit
pc <- ggplot(d,aes(Crop.coefficient,Deficit)) +
  geom_bar(aes(fill=aspect),position="dodge",stat="identity") +
  scale_y_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100)) +
  theme(text = element_text(size=20)) +
  scale_x_discrete(drop=FALSE) +
  xlab("") +
  ylab("") +
  ggtitle("CWD") +
  scale_fill_viridis(name="Solar exposure",discrete=TRUE, begin=0.05, end=0.95) +
  geom_hline(aes(color="Water input",yintercept = water),linetype=2,size=0.5) +
  scale_color_manual(name="",values = c("blue"),) +
  theme_bw(8) +
  theme(plot.title = element_text(hjust = 0.5))

Cairo(file="figures/Fig1.png",width=1200,height=500,res=200)
plot_grid(pa,pb,pc,nrow=1,rel_widths=c(1,1,1.55))
dev.off()

# 
# # scenario: water input
# 
# 
# d <- read.csv("water.csv",header=TRUE,stringsAsFactors=FALSE)
# #d$water <- reorder(d$water,c(1,2,3,1,2,3))
# d$water.in <- as.numeric(d$water)
# d$water <- as.factor(d$water.input)
# d$water <- reorder(d$water,c(1,2,3,1,2,3))
# 
# d$scenario <- reorder(d$scenario,c(1,1,1,2,2,2))
# 
# # Plot PET
# ggplot(d,aes(scenario,PET,fill=water)) +
#   geom_bar(position=position_dodge(),stat="identity") +
#   geom_errorbar(position=position_dodge(),aes(y=water.in,ymax=water.in,ymin=water.in),width=0.9,colour="blue",linetype=2,size=1) +
#   scale_y_continuous(limits=c(0,8),breaks=c(0,2,4,6,8)) +
#   scale_fill_brewer(palette="Set2") + theme(text = element_text(size=20)) +
#   scale_x_discrete(drop=FALSE) +
#   xlab("Scenario") +
#   ggtitle("PET")
# 
# # Plot AET
# ggplot(d,aes(scenario,AET,fill=water)) + geom_bar(position=position_dodge(),stat="identity") + geom_errorbar(position=position_dodge(),aes(y=water.in,ymax=water.in,ymin=water.in),width=0.9,colour="blue",linetype=2,size=1) + scale_y_continuous(limits=c(0,8),breaks=c(0,2,4,6,8)) + scale_fill_brewer(palette="Set2") + theme(text = element_text(size=20)) + scale_x_discrete(drop=FALSE) + xlab("Scenario")
# 
# # Plot Deficit
# ggplot(d,aes(scenario,Deficit,fill=water)) + geom_bar(position=position_dodge(),stat="identity") + scale_y_continuous(limits=c(0,8),breaks=c(0,2,4,6,8)) + scale_fill_brewer(palette="Set2") + theme(text = element_text(size=20)) + scale_x_discrete(drop=FALSE) + xlab("Scenario")
# 
# 
# library(gridExtra)
# 
# grid.arrange