
library(gridExtra)
library(glmnetUtils)
library(glmnet)
library(ROCR)
library(pROC)
library(broom)


mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)
patientsalldata6       <- read.csv( file="patientsalldata6.csv")
patientsalldata6       <- patientsalldata6$x
metabolom_urine_data1  <- read.csv(file="metabolom_urine_data1.csv")
metabolom_urine_data2  <- read.csv( file="metabolom_urine_data2.csv")
patients16             <- read.csv(file="patients16.csv")


##################################
# REDUCED MODEL - 10 COVARIATES
################################
#############
# INITIALIZE 
#############
set.seed(300)
`%notin%` <- Negate(`%in%`)
# Take 33 patients 
datatrain <- metabolom_urine_data2[metabolom_urine_data2$individual %in% patientsalldata6,]
nrow(datatrain[!duplicated(datatrain$individual),]) #33

##################################
# FIT EN MODEL ON 33 PATIENTS - 10 COVARIATES
################################
pe_train             <- datatrain$PE
datatrain$PE         <- NULL
datatrain$individual <- NULL
datatrain$gest_age   <- NULL
datatrain$Gestational.age.at.delivery <- NULL

X <- model.matrix(~., data=datatrain) 
X <- X[,-1]  
model  <- cv.glmnet(X, pe_train, family="binomial", type.measure = "auc",  alpha=.9, nfolds=4)
tmodel <- tidy(coef(model, s = "lambda.min"))
tmodel <- tmodel[-1,] #remove intercept

 # Choose 10 covariates
 tmodel$absvalue <- abs(tmodel$value)
 tmodel <- tmodel[order(tmodel$absvalue, decreasing=TRUE),]
 tmodel10 <- tmodel[c(1:10),]
 tmodel10 <-tmodel10[!is.na(tmodel10$column),]

# Fit the model with selected covariates only
seldatatrain <- subset(datatrain, select = c(tmodel10$row))
X <- model.matrix(~., data=seldatatrain) 
X <- X[,-1]  
model10 <-  cv.glmnet(X, pe_train, family="binomial", type.measure = "auc",   alpha=.9, nfolds=4)

################################
# TEST MODEL ON 16 PATIENTS
################################
datatest <- metabolom_urine_data1[metabolom_urine_data1$individual %notin% patientsalldata6,]
S <- datatest[,c("individual", "PE")]
S$improved <- NA
datatest$PE          <- NULL
datatest$individual  <- NULL
datatest$gest_age    <- NULL
datatest$Gestational.age.at.delivery<- NULL

# Calculate scores 
seldatatest <- subset(datatest, select = c(tmodel10$row))
Xtest <- model.matrix(~., data=seldatatest) 
Xtest <- Xtest[,-1] 
predict16 <- predict(model10, Xtest, s = "lambda.min",  type='response', alpha=0.9)
S$score     <- as.numeric(predict16)  


p16 <- patients16$individual
for (l in 1:length(p16)){	
	patient <- p16[c(l)]
	cc<-S$score[which(S$individual==patient)]
S$improved[which(S$individual==patient)]  <- sum(cc)/length(cc)
}

roc(S$PE,S$improved, direction="<", ci=TRUE) #0.8735
validatedroc <- roc(S$PE,S$improved, direction="<", ci=TRUE) 
plot(validatedroc, print.auc=TRUE, col="darkblue", lwd=3, cex.lab=1.5, main="AUC of the prediction model on the validation cohort")
validatedroc 


