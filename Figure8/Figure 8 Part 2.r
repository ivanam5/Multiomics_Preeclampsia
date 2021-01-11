install.packages('Rtsne')
library('Rtsne')
library('ggplot2')

mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

cytof <- read.csv('cytof_edited.csv')
cytof <- cytof[820:ncol(cytof)]
cytof

cytof <- cytof[,c(55,152,214,36,225,163,176,147)]
cytof

for(i in 1:ncol(cytof)){
    sum <- 0
    num <- 0
    for(j in 1:nrow(cytof)){
        if(!is.na(cytof[j,i])){
            sum = sum + cytof[j,i]
            num = num + 1
        }
    }
    avg <- sum/num
    for(j in 1:nrow(cytof)){
        if(is.na(cytof[j,i])){
           cytof[j,i] <- avg
        }
    }
}

proteome <- read.csv("C:/Users/andyt/Downloads/Datasets/proteome_edited.csv")
proteome <- proteome[7:ncol(proteome)]
proteome

proteins <- c('APCS','APOB','CCL23.1','HIPK3','SELE','SPARCL1','SELL','VEGFA.1','PRSS2','LEP','ROR1','IL24')
proteome <- proteome[,proteins]
proteome

ctrl_rows <- c(1,2,3,4,5,6,7,8,9,15,16,17,18,19,25,26,27,28,29,36,37)
pree_rows <- c(10,11,12,13,14,20,21,22,23,24,30,31,32,33,34,35,38,39,40,41,42,43,44,45,46,47)

proteome_ctrl <- proteome[1,]
for(i in 2:length(ctrl_rows)){
    proteome_ctrl <- rbind(proteome_ctrl,proteome[ctrl_rows[i],])
}
proteome_ctrl <- model.matrix(~., data=proteome_ctrl)
proteome_ctrl <- proteome_ctrl[,-1]

proteome_pree <- proteome[10,]
for(i in 2:length(pree_rows)){
    proteome_pree <- rbind(proteome_pree,proteome[pree_rows[i],])
}

proteome_pree <- model.matrix(~., data=proteome_pree)
proteome_pree <- proteome_pree[,-1]

cytof_ctrl <- cytof[1,]
for(i in 2:length(ctrl_rows)){
    cytof_ctrl <- rbind(cytof_ctrl,cytof[ctrl_rows[i],])
}

cytof_ctrl <- model.matrix(~., data=cytof_ctrl)
cytof_ctrl <- cytof_ctrl[,-1]

cytof_pree <- cytof[10,]
for(i in 2:length(pree_rows)){
    cytof_pree <- rbind(cytof_pree,cytof[pree_rows[i],])
}

cytof_pree <- model.matrix(~., data=cytof_pree)
cytof_pree <- cytof_pree[,-1]

proteome_name <- c()
cytof_name <- c()
pvals_pree <- c()
for(i in 1:ncol(proteome_pree)){
    for(j in 1:ncol(cytof_pree)){
        res <- cor.test(proteome_pree[,i],cytof_pree[,j], alternative = 'two.sided')
        pvals_pree <- c(pvals_pree,res$p.value)
        proteome_name <- c(proteome_name,colnames(proteome_pree)[i])
        cytof_name <- c(cytof_name,colnames(cytof_pree)[j])
    }
}

proteome_name2 <- c()
cytof_name2 <- c()
pvals_ctrl <- c()
for(i in 1:ncol(proteome_ctrl)){
    for(j in 1:ncol(cytof_ctrl)){
        res <- cor.test(proteome_ctrl[,i],cytof_ctrl[,j], alternative = 'two.sided')
        pvals_ctrl <- c(pvals_ctrl,res$p.value)
        proteome_name2 <- c(proteome_name2,colnames(proteome_ctrl)[i])
        cytof_name2 <- c(cytof_name2,colnames(cytof_ctrl)[j])
    }
}

pvaldf <- data.frame(pree=pvals_pree,ctrl=pvals_ctrl)
pvaldf$index <- as.numeric(1:nrow(pvaldf))
pvaldf

labels <- c()
for(i in 1:nrow(pvaldf)){
    labels <- c(labels,0)
}
# LEP
labels[77] <- 1
labels[75] <- 2
labels[74] <- 3
labels[73] <- 4
labels[78] <- 5
labels[79] <- 6
labels[80] <- 7

#SELL
labels[52] <- 8
labels[55] <- 9
labels[49] <- 10
labels[50] <- 11
labels[56] <- 12

#VEGFA
labels[57] <- 13
labels[58] <- 14
labels[59] <- 15
labels[60] <- 16
labels[61] <- 17
labels[64] <- 18

#APCS
labels[7] <- 19
labels[5] <- 20

#HIPK3
labels[26] <- 21
labels[31] <- 22
labels[29] <- 23

#PRSS2
labels[68] <- 24
labels[67] <- 25
labels[70] <- 26

#CCL23
labels[22] <- 27

ggplot(pvaldf, aes(x=-log10(pree), y=-log10(ctrl))) +
geom_point(size = 2.5, color='#23ab24') + 
geom_vline(xintercept = 1.3,color="grey", linetype="dashed") + geom_hline(yintercept = 1.3, color="grey", linetype="dashed") +
geom_text(label=ifelse(labels == 10,labels,""), size = 8, nudge_x = -0.1, nudge_y = 0.075) +
geom_text(label=ifelse(labels == 7,labels,""), size = 8, nudge_x = 0.075, nudge_y = -0.075) +
geom_text(label=ifelse(labels == 13,labels,""), size = 8, nudge_x = 0.075, nudge_y = -0.075) +
geom_text(label=ifelse(labels == 17,labels,""), size = 8, nudge_x = 0.1, nudge_y = -0.05) +
geom_text(label=ifelse(labels == 16,labels,""), size = 8, nudge_x = -0.1, nudge_y = 0.075) +
geom_text(label=ifelse(labels == 15,labels,""), size = 8, nudge_x = -0.1, nudge_y = 0.075) +
geom_text(label=ifelse(labels == 23,labels,""), size = 8, nudge_x = -0.0875, nudge_y = -0.0875) +
geom_text(label=ifelse(labels == 27,labels,""), size = 8, nudge_x = -0.0875, nudge_y = 0.075) +
geom_text(label=ifelse(labels == 26,labels,""), size = 8, nudge_x = -0.1, nudge_y = -0.075) +
geom_text(label=ifelse(labels == 25,labels,""), size = 8, nudge_x = 0.125, nudge_y = -0.025) +
geom_text(label=ifelse(labels == 21,labels,""), size = 8, nudge_x = 0.0825, nudge_y = -0.0825) +
geom_text(label=ifelse(labels == 14,labels,""), size = 8, nudge_x = 0.0625, nudge_y = 0.0825) +
geom_text(label=ifelse(labels == 3,labels,""), size = 8, nudge_x = 0.075, nudge_y = 0.05) +
geom_text(label=ifelse(labels > 0 & labels < 20 & labels != 3 & labels != 7  & labels != 14 & labels != 13 & labels != 10 & labels != 17 & labels != 16 & labels != 15,labels,""), size = 8, nudge_x = 0.075, nudge_y = 0.075) +
geom_text(label=ifelse(labels > 0 & labels > 19 & labels != 25 & labels != 23 & labels != 27 & labels != 26 & labels != 21,labels,""), size = 8, nudge_x = 0.125, nudge_y = 0.1) +
theme_bw() +  theme(axis.line = element_line(colour = "black")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.title = element_blank()) +theme(panel.border = element_blank())+ theme(panel.background = element_rect(fill = "transparent")) +theme(plot.background = element_rect(fill = "transparent", color = NA)) + theme(legend.background = element_rect(fill="transparent")) +
theme(axis.title=element_text(size=25)) + theme(axis.text=element_text(size=25))

ggsave("C:/Users/andyt/Downloads/cytofproteomeplotplain_gest.png", device='png')


