library(ggplot2)
library(gridExtra)

######################
# DISCOVERY COHORT
######################
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

proteom   <- read.csv(file="proteom.csv")
gadelivery <- read.csv(file="gadelivery.csv")
T          <- merge(proteom, gadelivery, by="individual")
colnames(T)[which(names(T)=="Gestational.age.at.delivery")] <-"gadelivery" 
T          <- T[order(T$PE, T$gadelivery, decreasing=TRUE),]

proteom1    <- T[!duplicated(T$individual),]
proteom1$ID <- c(1:nrow(proteom1))
map         <- proteom1[,c("individual", "ID")]
T           <- merge (T, map, by="individual")
T           <- T[order(T$PE, T$gadelivery, decreasing=TRUE),]

uniqueID = data.frame(pID = unique(T$ID))

  

T$colorcode <- c("#56B4E9")
T$colorcode[T$PE==0] <- c("#999999")
  
    ggplot(T)+
    geom_rect(aes(xmin = -1, xmax = 43, ymin = -1, ymax = 35), fill = 'gray91') +
    geom_segment(data = uniqueID, aes(x = 0, xend = 42, y = pID, yend = pID), color = 'white', size=0.8) +
geom_segment(aes(x = 12, y = 1, xend = 12, yend = 34), linetype = 'dashed', color = 'white', size = 1) +
geom_segment(aes(x = 26, y = 1, xend = 26, yend = 34), linetype = 'dashed', color = 'white', size = 1) +
geom_segment(aes(x = 37, y = 1, xend = 37, yend = 34), linetype = 'dashed', color = 'white', size = 1) +
     geom_point(aes(x=ga, y=ID, color = colorcode), size = 3) +
    scale_color_identity(guide = "legend", labels = c("PE", "Control"))+
      geom_point(aes(x=gadelivery, y=ID,   color = colorcode,), size = 4, shape = 25)  +
    labs(x="Gestational age (weeks)", y="Patients") +
    ggtitle("Sample Collection Times: Cohort 1") +
     scale_x_continuous(expand = c(0, 0))+
 scale_y_continuous(expand = c(0, 0))+
    geom_text(aes(label = '0', x=0, y=-0.1), size = 4) +
  geom_text(aes(label = '12', x=12, y=-0.1), size = 4) +
  geom_text(aes(label = '26', x=26, y=-0.1), size = 4) +
  geom_text(aes(label = '37', x=37, y=-0.1), size = 4) +
  theme(axis.title.x = element_text(size=18)) +
  theme(plot.title = element_text(size = 18))+
   theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) + theme(legend.position = "none") 
  
  

######################
# VALIDATION COHORT
######################
urine      <- read.csv(file="/Users/Ivana/WORKSPACE/PAPEROMICS/SUBMISSIONNATURECOMMUNICATIONS/CODE/urine.csv")
patients16 <- read.csv(file="/Users/Ivana/WORKSPACE/PAPEROMICS/SUBMISSIONNATURECOMMUNICATIONS/CODE/patients16.csv")
Q <- merge(urine, patients16, by="individual")
colnames(Q)[which(names(Q)=="Gestational.age.at.delivery")] <-"gadelivery" 
colnames(Q)[which(names(Q)=="gest_age")] <-"ga" 
Q <- Q[order(Q$PE, Q$gadelivery, decreasing=TRUE),]

urine1 <- Q[!duplicated(Q$individual),]
urine1$ID <- c(1:nrow(urine1))
map <- urine1[,c("individual", "ID")]
Q <- merge (Q, map, by="individual")
Q <- Q[order(Q$PE, Q$gadelivery, decreasing=TRUE),]

uniqueID = data.frame(pID = unique(Q$ID))

Q$Patient <- c("#56B4E9")
Q$Patient[Q$PE==0] <- c("#999999")
  
 
ggplot(Q)+
geom_rect(aes(xmin = -1, xmax = 43, ymin = -1, ymax = 17), fill = 'gray91') +
geom_segment(data = uniqueID, aes(x = 0, xend = 42, y = pID, yend = pID), color = 'white', size=0.8) +
geom_segment(aes(x = 12, y = 1, xend = 12, yend = 16), linetype = 'dashed', color = 'white', size = 1) +
geom_segment(aes(x = 26, y = 1, xend = 26, yend = 16), linetype = 'dashed', color = 'white', size = 1) +
geom_segment(aes(x = 37, y = 1, xend = 37, yend = 16), linetype = 'dashed', color = 'white', size = 1)  +
geom_point(aes(x=ga, y=ID, color = Patient), size = 3) +
scale_color_identity(guide = "legend", labels = c("PE", "Control"))+
geom_point(aes(x=gadelivery, y=ID,   color = Patient), size = 4, shape = 25)  +
labs(x="Gestational age (weeks)", y="Patients") +
ggtitle("Sample Collection Times: Cohort 2") +
scale_x_continuous(expand = c(0, 0))+
scale_y_continuous(expand = c(0, 0))+
  geom_text(aes(label = '0', x=0, y=-0.1), size = 4) +
  geom_text(aes(label = '12', x=12, y=-0.1), size = 4) +
  geom_text(aes(label = '26', x=26, y=-0.1), size = 4) +
  geom_text(aes(label = '37', x=37, y=-0.1), size = 4) +
  theme(axis.title.x = element_text(size=18)) +
  theme(plot.title = element_text(size = 18))+
   theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +  theme(legend.text=element_text(size=12)) 
  
  




  