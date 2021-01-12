library(glmnetUtils)
library(nnls)
library(glmnet)
library(pROC)
library(ROCR)
library(broom)
library(MASS)
#####################################
# Data
proteom_data           <- read.csv(file="proteom_data.csv")
metabolom_urine_data2  <- read.csv( file="metabolom_urine_data2.csv")
patientsalldata6       <- read.csv( file="patientsalldata6.csv")
foldidAll              <- read.csv(file="foldidAll.csv")
proteom_data$foldid           <-foldidAll 
metabolom_urine_data2$foldid  <-foldidAll 
patientsalldata6              <- patientsalldata6$x

#############################################
 # URINE METABOLOME MODEL WITH 10 COVARIATES
##############################################
set.seed(300)
# Initialize 
S               <- proteom_data[,c("SampleId", "individual", "PE")]
S$s_metaburine  <- 0


for (i in 1:length(patientsalldata6)){	
   leaveoutpatient <- patientsalldata6[c(i)]
   print(i)   
######################################
# 3. TRAIN SEPARATE MODELS for each j
######################################
# Metabolome urine data   
######################################
    datatrain       <-  metabolom_urine_data2[!(metabolom_urine_data2$individual == leaveoutpatient),]
    pe_train        <-  metabolom_urine_data2$PE[!(metabolom_urine_data2$individual== leaveoutpatient)]
    datatest        <-  metabolom_urine_data2[metabolom_urine_data2$individual == leaveoutpatient,]
    datatrain$PE    <- NULL
    datatest$PE          <- NULL
    datatrain$individual <- NULL
    datatest$individual  <- NULL
    fold=datatrain$foldid
    datatrain$foldid<-NULL
    datatest$foldid<-NULL

##################
# FIT WITH EN ONLY
##################
X <- model.matrix(~., data=datatrain) 
X <- X[,-1]  


fit_metaburine <- cv.glmnet(X, pe_train, family="binomial", type.measure = "auc", foldid=fold,  alpha=.9, nfolds=4)
t_modelpe <- tidy(coef(fit_metaburine, s = "lambda.min"))
t_modelpe <- t_modelpe[-1,] #remove intercept
 
 # Choose 10 covariates
 t_modelpe$absvalue <- abs(t_modelpe$value)
 t_modelpe <- t_modelpe[order(t_modelpe$absvalue, decreasing=TRUE),]
 t_modelpe10 <- t_modelpe[c(1:10),]
 t_modelpe10 <-t_modelpe10[!is.na(t_modelpe10$column),]


#########################################################################################
# Run EN only with selected covariates to obtain final prediction model
#########################################################################################
# Put in right format
######################
# Fit the model with selected covariates only
seldatatrain <- subset(datatrain, select = c(t_modelpe10$row))
X <- model.matrix(~., data=seldatatrain) 
X <- X[,-1]  
finalmodel_metaburine <-  cv.glmnet(X, pe_train, family="binomial", type.measure = "auc", foldid=fold,  alpha=.9, nfolds=4)

# Calculate scores for patient j
seldatatest <- subset(datatest, select = c(t_modelpe10$row))
Xtest <- model.matrix(~., data=seldatatest) 
Xtest <- Xtest[,-1]
s_metabolome_urine <- predict(finalmodel_metaburine, Xtest, s = "lambda.min",  type='response', alpha=0.9)

 S$s_metaburine[S$individual==leaveoutpatient] <-  s_metabolome_urine
 
 } # for j
 
 S$s_metaburineimproved<-0
for (l in 1:length(patientsalldata6)){	
	patient <- patientsalldata6[c(l)]
	cc<-S$s_metaburine[which(S$individual==patient)]
	#print(patient)
	#print(cc)
S$s_metaburineimproved[which(S$individual==patient)]  <- sum(cc)/length(cc)
}
 

roc10metaburine <- roc(S$PE,S$s_metaburineimproved, ci=TRUE)
roc10metaburine

    
###########################################
# SAME FOR PROTEOME: MODEL WITH 10 COVARIATES
############################################
set.seed(300)
# Initialize 
S               <- proteom_data[,c("SampleId", "individual", "PE")]
S$s_proteome    <- 0


for (i in 1:length(patientsalldata6)){	
   leaveoutpatient <- patientsalldata6[c(i)]
   print(i)   
######################################
# 3. TRAIN SEPARATE MODELS for each j
######################################
# Proteome data   
######################################
    datatrain       <-  proteom_data[!(proteom_data$individual == leaveoutpatient),]
    pe_train        <-  proteom_data$PE[!(proteom_data$individual== leaveoutpatient)]
    datatest        <-  proteom_data[proteom_data$individual == leaveoutpatient,]
    datatrain$PE         <- NULL
    datatest$PE          <- NULL
    datatrain$individual <- NULL
    datatest$individual  <- NULL
    fold=datatrain$foldid
    datatrain$foldid     <-NULL
    datatest$foldid      <-NULL

##################
# FIT WITH EN ONLY
##################
X <- model.matrix(~., data=datatrain) 
X <- X[,-1]  

fit_proteome <- cv.glmnet(X, pe_train, family="binomial", type.measure = "auc", foldid=fold,  alpha=.9, nfolds=4)
t_modelpe <- tidy(coef(fit_proteome, s = "lambda.min"))
t_modelpe <- t_modelpe[-1,] #remove intercept
 
 
 # Choose 10 covariates
 t_modelpe$absvalue <- abs(t_modelpe$value)
 t_modelpe <- t_modelpe[order(t_modelpe$absvalue, decreasing=TRUE),]
 t_modelpe10 <- t_modelpe[c(1:10),]
 t_modelpe10 <-t_modelpe10[!is.na(t_modelpe10$column),]


#########################################################################################
# Run logistic regression only with selected covariates to obtain final prediction model
#########################################################################################
# Put in right format
######################
# Fit the model with selected covariates only
seldatatrain <- subset(datatrain, select = c(t_modelpe10$row))
X <- model.matrix(~., data=seldatatrain) 
X <- X[,-1]  
finalmodel_proteome <-  cv.glmnet(X, pe_train, family="binomial", type.measure = "auc", foldid=fold,  alpha=.9, nfolds=4)

# Calculate scores for patient j
seldatatest <- subset(datatest, select = c(t_modelpe10$row))
Xtest <- model.matrix(~., data=seldatatest) 
Xtest <- Xtest[,-1]
s_proteome <- predict(finalmodel_proteome, Xtest, s = "lambda.min",  type='response', alpha=0.9)

 S$s_proteome[S$individual==leaveoutpatient] <-  s_proteome
 
 } # for j
 
 S$s_proteomeimproved<-0
for (l in 1:length(patientsalldata6)){	
	patient <- patientsalldata6[c(l)]
	cc<-S$s_proteome[which(S$individual==patient)]
	#print(patient)
	#print(cc)
S$s_proteomeimproved[which(S$individual==patient)]  <- sum(cc)/length(cc)
}
 

roc10proteome <- roc(S$PE,S$s_proteomeimproved, ci=TRUE, quiet=TRUE, direction="<")
roc10proteome


##############
#FIG 3B
###############
plot(roc10metaburine, print.auc=FALSE, col="#3F51B5",  lwd=3, cex.lab=1.8,  xlab="False Positive Rate", ylab="True Positive Rate")
lines(roc10proteome, col="#F57C00")
legend(0.85,0.15, legend=c("Urine metabolome (AUC=0.88)", "Proteome (AUC=0.83)"), col = c("#3F51B5", "#F57C00"), lty = 1.5, cex = 1.6)
   
################# END ##########################

