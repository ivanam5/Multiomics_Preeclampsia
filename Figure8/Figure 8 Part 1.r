install.packages('Rtsne')
library('Rtsne')
library('ggplot2')

mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

cytof_original <- read.csv('immune_dataFINAL.csv')
T1 <- c(1,4,6,8,10,14,16,20,24,30,34,38,42,46,49,51,55,64,66)
T2 <- c(3,5,7,9,12,15,18,22,25,32,36,40,44,48,50,53,56,65,68)
cytof_ge <- t(cytof_original[,2])
cytof_original <- t(cytof_original[819:ncol(cytof_original)])
f <- (cytof_original[,T2[1]] - cytof_original[,T1[1]])/(cytof_ge[,T2[1]] - cytof_ge[,T1[1]])
cytof_data <- f
for(i in 2:length(T1)){
    f <- (cytof_original[,T2[i]] - cytof_original[,T1[i]])/(cytof_ge[,T2[i]] - cytof_ge[,T1[i]])
    cytof_data <- cbind(cytof_data,f)
}
cytof_data <- t(as.data.frame(cytof_data))
cytof_data <- as.data.frame(cytof_data)
cytof_data['Individual'] <- as.data.frame(c(10523,10524,10541,10547,10634,10639,10712,10721,10722,10727,10732,10733,10736,10803,10804,10809,10816,19017,19600))
cytof_data

for(i in 1:ncol(cytof_data)){
    sum <- 0
    num <- 0
    for(j in 1:nrow(cytof_data)){
        if(!is.na(cytof_data[j,i])){
            sum = sum + cytof_data[j,i]
            num = num + 1
        }
    }
    avg <- sum/num
    print(avg)
    for(j in 1:nrow(cytof_data)){
        if(is.na(cytof_data[j,i])){
           cytof_data[j,i] <- avg
        }
    }
}

score <- read.csv('C:/Users/andyt/Downloads/SCORE2.csv')
scores <- score$s_improvednnls
index <- c(scores[14],scores[17],scores[25],scores[31],scores[42],scores[45],scores[48],scores[53],scores[56],scores[59],scores[62],scores[65],scores[67],scores[70],scores[73],scores[75],scores[78],scores[81],scores[86])
cytof_data$score <- index
cytof_data

pvals <- c()
for(i in 1:371){
    res <- cor.test(cytof_data[,i],cytof_data$score, alternative = 'two.sided')
    pvals <- c(pvals,res$p.value)
}
pvals

df <- read.csv('C:/Users/andyt/Downloads/cytofclusters.csv')

mask <- c()
for(i in 1:nrow(df)){
    mask <- c(mask,0)
}
mask[55] <- 2
mask[152] <- 2
mask[214] <- 2
mask[36] <- 2
mask[225] <- 2
mask[163] <- 2
mask[176] <- 2
mask[147] <- 1

mask[221] <- 1
mask[47] <- 1
mask[176] <- 2
mask[45] <- 1
mask[143] <- 1
df$mask <- mask
df$X <- NULL
df

ggplot(as.data.frame(df),aes(x = V1, y = V2)) + 
geom_point(colour = "#23ab24", size = ifelse(mask == 0 & -log10(pvals) > 1,-log10(pvals),0)) + 
geom_point(colour = "#23ab24", size = ifelse(mask == 0 & -log10(pvals) < 1,1,0)) + 
geom_point(fill = "#23ab24", colour = 'black', size = ifelse(mask == 1 & -log10(pvals) > 1,-log10(pvals),-1), pch = 21) + 
geom_point(fill = "#23ab24", colour = 'black', size = ifelse(mask == 1 & -log10(pvals) < 1,1,-1), pch = 21) + 
geom_point(fill = "#F57C00", colour = 'black', size = ifelse(mask == 2 & -log10(pvals) > 1,-log10(pvals),-1), pch = 21) + 
geom_point(fill = "#F57C00", colour = 'black', size = ifelse(mask == 2 & -log10(pvals) < 1,1,-1), pch = 21) + 
theme_bw() + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.title = element_blank()) + 
theme(panel.background = element_rect(fill = "transparent")) +
theme(plot.background = element_rect(fill = "transparent", color = NA)) +
theme(legend.background = element_rect(fill="transparent")) + theme(legend.text=element_text(size=11)) + theme_void() + 
theme(legend.position = "none")

ggsave("C:/Users/andyt/Downloads/figure13.png", device = 'png', dpi = 300)


