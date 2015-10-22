library(caret)

setwd("~/Documents/My Projects/Predict HIV Progression")


hiv_train <- read.csv("training_data.csv", header=TRUE)

hiv_train$Resp <- factor(hiv_train$Resp)

alaA <- c(c("GCT","GCC","GCA"), c("GCG","GCN"), c("GCU")) # Alanine
argR <- c(c("CGT","CGC","CGA","CGG","AGA","AGG"), 
          c("CGN","MGR"), c("CGU")) # Arginine
asnN <- c(c("AAT","AAC"), c("AAY")) # Asparagine
aspD <- c(c("GAT","GAC"), c("GAY"), c("GAU")) # Aspartic acid
cysC <- c(c("TGT","TGC"), c("TGY"), c("UGU", "UGC")) # Cysteine
ginQ <- c(c("CAA","CAG"), c("CAR")) # Glutamine
gluE <- c(c("GAA","GAG"), c("GAR")) # Glutamic acid
glyG <- c(c("GGT","GGC","GGA","GGG"), c("GGN"), c("GGU")) # Glycine
hisH <- c(c("CAT","CAC"), c("CAY"), c("CAU")) # Histidine
iieI <- c(c("ATT","ATC","ATA"), c("ATH"), c("AUU", "AUC", "AUA")) # Isoleucine
startC <- c(c("ATG")) # start1
leuL <- c(c("TTA","TTG","CTT","CTC","CTA","CTG"), c("YTR","CTN"), 
          c("UUA","UUG","CUU","CUC","CUA","CUG")) # Leucine
lysK <- c(c("AAA","AAG"), c("AAR")) # Lysine
metM <- c(c("ATG"), c("AUG")) # Methionine
pheF <- c(c("TTT","TTC"), c("TTY"), c("UUU","UUC")) # Phenylalanine
proP <- c(c("CCT","CCC","CCA","CCG"), c("CCN"), c("CCU")) # Proline
serS <- c(c("TCT","TCC","TCA","TCG","AGT","AGC"), c("TCN","AGY"), 
          c("UCU","UCC","UCA","UCG")) # Serine
thrT <- c(c("ACT","ACC","ACA","ACG"), c("ACN"), c("ACU")) # Threonine     
trpW <- c(c("TGG"), c("UGG")) # Tryptophan   
tyrY <- c(c("TAT","TAC"), c("TAY"), c("UAU","UAC")) # Tyrosine
valV <- c(c("GTT","GTC","GTA","GTG"), c("GTN"), c("GUU","GUC","GUA","GUG")) # Valine
stopC <- c(c("TAA","TGA","TAG"), c("TAR","TRA"), c("UAA","UAG","UGA")) # Stopp1


alanineRT <- 0; arginineRT <- 0; asparagineRT <- 0;  
asparticAcidRT <- 0; cysteineRT <- 0; glutamineRT <- 0;
glutaminAcidRT <- 0; glycineRT <- 0; histidineRT <- 0; 
isoleucineRT <- 0; leucineRT <- 0; lysineRT <- 0; 
methionineRT <- 0; phenylalanineRT <- 0; prolineRT <- 0; 
serineRT <- 0; threonineRT <- 0; tryptophanRT <- 0;
tyrosineRT <- 0; valineRT <- 0; startCodonRT <- 0; 
stopCodonRT <- 0;
alaninePR <- 0; argininePR <- 0; asparaginePR <- 0; 
asparticAcidPR <- 0; cysteinePR <- 0; glutaminePR <- 0;
glutaminAcidPR <- 0; glycinePR <- 0; histidinePR <- 0; 
isoleucinePR <- 0; leucinePR <- 0; lysinePR <- 0; 
methioninePR <- 0; phenylalaninePR <- 0; prolinePR <- 0; 
serinePR <- 0; threoninePR <- 0; tryptophanPR <- 0;
tyrosinePR <- 0; valinePR <- 0; startCodonPR <- 0; 
stopCodonPR <- 0;

hiv_train$alanineRT <- 0; hiv_train$arginineRT <- 0;  hiv_train$asparagineRT <- 0; 
hiv_train$asparticAcidRT <- 0; hiv_train$cysteineRT <- 0; hiv_train$glutamineRT <- 0;
hiv_train$glutaminAcidRT <- 0; hiv_train$glycineRT <- 0;hiv_train$histidineRT <- 0; 
hiv_train$isoleucineRT <- 0; hiv_train$leucineRT <- 0; hiv_train$lysineRT <- 0; 
hiv_train$methionineRT <- 0; hiv_train$phenylalanineRT <- 0; hiv_train$prolineRT <- 0; 
hiv_train$serineRT <- 0; hiv_train$threonineRT <- 0; hiv_train$tryptophanRT <- 0; 
hiv_train$tyrosineRT <- 0; hiv_train$valineRT <- 0; hiv_train$startCodonRT <- 0; 
hiv_train$stopCodonRT <- 0;
hiv_train$alaninePR <- 0; hiv_train$argininePR <- 0; hiv_train$asparaginePR <- 0; 
hiv_train$asparticAcidPR <- 0; hiv_train$cysteinePR <- 0; hiv_train$glutaminePR <- 0;
hiv_train$glutaminAcidPR <- 0; hiv_train$glycinePR <- 0;hiv_train$histidinePR <- 0; 
hiv_train$isoleucinePR <- 0; hiv_train$leucinePR <- 0; hiv_train$lysinePR <- 0; 
hiv_train$methioninePR <- 0;  hiv_train$phenylalaninePR <- 0; hiv_train$prolinePR <- 0; 
hiv_train$serinePR <- 0;  hiv_train$threoninePR <- 0; hiv_train$tryptophanPR <- 0; 
hiv_train$tyrosinePR <- 0;  hiv_train$valinePR <- 0; hiv_train$startCodonPR <- 0; 
hiv_train$stopCodonPR <- 0;

hiv_train$totalSeqRT <- 0; hiv_train$totalCodonRT <- 0;
hiv_train$totalSeqPR <- 0; hiv_train$totalCodonPR <- 0;

# source("dnaRT.R")
# debug(dnaRT)
# hiv_train1 <- dnaRT(hiv_train)
# undebug(dnaRT)
# 
# source("dnaPR.R")
# debug(dnaPR)
# hiv_train2 <- dnaPR(hiv_train1)
# undebug(dnaPR)


source("dnaRT.R")
hiv_train1 <- dnaRT(hiv_train)
write.table(hiv_train1, "hiv_train1.csv", sep=",")

source("dnaPR.R")
hiv_train2 <- dnaPR(hiv_train1)
write.table(hiv_train2, "hiv_train2.csv", sep=",")


hiv_train2$alanine <- hiv_train2$alanineRT + hiv_train2$alaninePR
hiv_train2$arginine <- hiv_train2$arginineRT + hiv_train2$argininePR
hiv_train2$asparagine  <- hiv_train2$asparagineRT + hiv_train2$asparaginePR
hiv_train2$asparticAcid <- hiv_train2$asparticAcidRT + hiv_train2$asparticAcidPR
hiv_train2$cysteine <- hiv_train2$cysteineRT + hiv_train2$cysteinePR
hiv_train2$glutamine <- hiv_train2$glutamineRT + hiv_train2$glutaminePR
hiv_train2$glutaminAcid <- hiv_train2$glutaminAcidRT + hiv_train2$glutaminAcidPR
hiv_train2$glycine <- hiv_train2$glycineRT + hiv_train2$glycinePR 
hiv_train2$histidine <- hiv_train2$histidineRT + hiv_train2$histidinePR 
hiv_train2$isoleucine <- hiv_train2$isoleucineRT + hiv_train2$isoleucinePR 
hiv_train2$leucine <- hiv_train2$leucineRT + hiv_train2$leucinePR 
hiv_train2$lysine <- hiv_train2$lysineRT + hiv_train2$lysinePR 
hiv_train2$methionine <- hiv_train2$methionineRT + hiv_train2$methioninePR 
hiv_train2$phenylalanine <- hiv_train2$phenylalanineRT + hiv_train2$phenylalaninePR 
hiv_train2$proline <- hiv_train2$prolineRT + hiv_train2$prolinePR 
hiv_train2$serine <- hiv_train2$serineRT + hiv_train2$serinePR 
hiv_train2$threonine <- hiv_train2$threonineRT + hiv_train2$threoninePR 
hiv_train2$tryptophan <- hiv_train2$tryptophanRT + hiv_train2$tryptophanPR 
hiv_train2$tyrosine <- hiv_train2$tyrosineRT + hiv_train2$tyrosinePR 
hiv_train2$valine <- hiv_train2$valineRT + hiv_train2$valinePR 
hiv_train2$startCodon <- hiv_train2$startCodonRT + hiv_train2$startCodonPR 
hiv_train2$stopCodon <- hiv_train2$stopCodonRT + hiv_train2$stopCodonPR
hiv_train2$totalSeq <- hiv_train2$totalSeqRT + hiv_train2$totalSeqPR
hiv_train2$totalCondon <- hiv_train2$totalCodonRT + hiv_train2$totalCodonPR


hiv_train3 <- hiv_train2

hiv_train3$PatientID <- NULL
hiv_train3$PR.Seq <- NULL
hiv_train3$RT.Seq <- NULL


table(hiv_train3$Resp)

# qplot() based on Resp vs 3 variables

# -------------------------------------------------------
# partition the data set into training and testing
# -------------------------------------------------------

set.seed(1234) # you will need to shuffle it
inTrain <- createDataPartition(y=hiv_train3$Resp, p=0.7, list=FALSE)
training <- hiv_train3[inTrain,] # Training data - The first 70%?
testing <- hiv_train3[-inTrain,] # Testing data - The rest of 30%?

dim(training)
dim(testing)

(dim(testing)[1] / dim(training)[1]) * 100 # This is not always 100 percent



# -------------------------------------------------------
# Full initial model using GLM 
# -------------------------------------------------------


# I think maybe the caret package has train function directly
fit_full_glm_1 <- train(Resp ~ ., method="glm", data=training)
summary(fit_full_glm_1)

# Alternatively the full model using glm
fit_full_glm_2 <- glm(Resp ~ . , data=training, family="binomial")
summary(fit_full_glm_2)

# Both approaches generate exact same model
AIC(fit_full_glm_1)
AIC(fit_full_glm_2)


# -------------------------------------------------------
# Use step to come up with the best model from GLM
# -------------------------------------------------------


library(MASS)

# Same best model
fit_full_glm_1_best <- step(fit_full_glm_1, direction = "both")
summary(fit_full_glm_1_best)

fit_full_glm_2_best <- step(fit_full_glm_2, direction = "both")
summary(fit_full_glm_2_best)




# -------------------------------------------------------
# Use step to come up with the best model from tree
# Tree doesn't seem to work at all in this program
# -------------------------------------------------------

fit_full_rf_best <- step(fit_full_rf, direction = "both")
summary(fit_full_rf_best)

fit_full_rpart_best <- step(fit_full_rpart, direction = "both")
summary(fit_full_rpart_best)




      
#----------------------------------------------------------
# Sum of condon from best GLM
#----------------------------------------------------------

fit_glm_1_sum <- glm(Resp ~ alanine + arginine + asparagine + asparticAcid + 
                   cysteine + glutamine + glutaminAcid + glycine + 
                   histidine + isoleucine + leucine + lysine + 
                   methionine + phenylalanine + proline + serine + 
                   threonine + tryptophan + tyrosine + valine + 
                   stopCodon + startCodon + totalSeq + totalCondon, 
               data=training, family="binomial")
summary(fit_glm_1_sum)


fit_glm_1_sum_best <- step(fit_glm_1_sum, direction = "both")
summary(fit_glm_1_sum_best)


#----------------------------------------------------------
# PR sequence from best GLM
#----------------------------------------------------------

fit_glm_1_PR <- glm(Resp ~ alaninePR + argininePR + asparaginePR + 
                  asparticAcidPR + cysteinePR + glutaminePR + 
                  glutaminAcidPR + glycinePR + histidinePR + 
                  isoleucinePR + leucinePR + lysinePR + 
                  methioninePR + phenylalaninePR + prolinePR + 
                  serinePR + threoninePR + tryptophanPR + 
                  tyrosinePR + valinePR + stopCodonPR + 
                  startCodonPR + totalSeqPR + totalCodonPR, 
              data=training, family="binomial")
summary(fit_glm_1_PR)

fit_glm_1_PR_best <- step(fit_glm_1_PR, direction="both")
summary(fit_glm_1_PR_best)



#----------------------------------------------------------
# RT sequence from best GLM
#----------------------------------------------------------

fit_glm_1_RT <- glm(Resp ~ alanineRT + arginineRT + asparagineRT + 
                  asparticAcidRT + cysteineRT + glutamineRT + 
                  glutaminAcidRT + glycineRT + histidineRT + 
                  isoleucineRT + leucineRT + lysineRT + 
                  methionineRT + phenylalanineRT + prolineRT + 
                  serineRT + threonineRT + tryptophanRT + 
                  tyrosineRT + valineRT + stopCodonRT + 
                  startCodonRT + totalSeqRT + totalCodonRT, 
              data=training, family="binomial")
summary(fit_glm_1_RT)

fit_glm_1_RT_best <- step(fit_glm_1_RT, direction="both")
summary(fit_glm_1_RT_best)




AIC(fit_full_glm_1)
AIC(fit_full_glm_1_best) # The smallest AIC, we should pick this one

AIC(fit_full_glm_2)
AIC(fit_full_glm_2_best) # The smallest AIC, we should pick this one


AIC(fit_glm_1_sum)
AIC(fit_glm_1_sum_best)
AIC(fit_glm_1_PR)
AIC(fit_glm_1_PR_best)
AIC(fit_glm_1_RT)
AIC(fit_glm_1_RT_best)


#---------------------------------------------#
#       Generalized Regression Model
#---------------------------------------------#
#       Predicting on the Testing set
#---------------------------------------------#

fit_full_glm_2_best_test <- predict(fit_full_glm_2_best, newdata=testing, type="response")
fit_full_glm_2_best_test


table(testing$Resp, fit_full_glm_2_best_test > 0.5)

testing$Resp2 <- 0
for (i in 1:nrow(testing)){
    if (fit_full_glm_2_best_test[i] > 0.5) {
        testing$Resp2[i] <- 1
    }
}


testing$Resp2 <- factor(testing$Resp2)


testing$Res <- FALSE
for (i in 1:nrow(testing)){
    if (testing$Resp[i] == testing$Resp2[i]) {
        testing$Res[i] <- TRUE
    }
}

sum(testing$Res)
testing$Resp
testing$Resp2


233/299

299-233

# anova(fit_full, fit_full_best, 
#       fit_PR, fit_PR_best, 
#       fit_RT, fit_RT_best, 
#       fit_sum, fit_sum_best)
# 
# anova(fit_full_best, fit_PR_best, fit_RT_best, fit_sum_best)

# anova(fit_full_best, test="Chisq")

confint(fit_full_best)
confint.default(fit_full_best)



# plot(predict(fit_full_best), resid(fit_full_best), pch = ".")
# 
# par(mfrow = c(2,2))
# plot(fit_full_best)
# dev.off()





# Find the random forest to use here?
# we will compare random forest to glm



# -------------------------------------------------------
# Regression and Classification Tree from a full initial model (failure)
# -------------------------------------------------------

# fit_full_rf <-  train(Resp ~ ., method="rf", data=training, 
#                       trControl=trainControl(method="cv"), number=3)
# print(fit_full_rf$finalModel)
# plot(fit_full_rf$finalModel, uniform = TRUE, main="Classification Tree")      


# fit_full_rpart <- train(Resp ~ ., method="rpart", data=training)
# print(fit_full_rpart$finalModel)
# plot(fit_full_rpart$finalModel, uniform = TRUE, main="Classification Tree")      


#---------------------------------------------#
#               Decision Tree
#---------------------------------------------#
#   Maybe you can stick with the full model
#---------------------------------------------#


# It doesn't produce a tree, just a root
fit_glm_1_sum <- train(Resp ~ alanine + arginine + asparagine + asparticAcid + 
                         cysteine + glutamine + glutaminAcid + glycine + 
                         histidine + isoleucine + leucine + lysine + 
                         methionine + phenylalanine + proline + serine + 
                         threonine + tryptophan + tyrosine + valine + 
                         stopCodon + startCodon + totalSeq + totalCondon, 
                     method = "rpart", data=training)
summary(fit_glm_1_sum)
print(fit_glm_1_sum$finalModel)
plot(fit_glm_1_sum$finalModel, uniform = TRUE, main="Classification Tree")      


# You should only use this, as it produces less errors I think

library(rattle)
library(rpart.plot)


set.seed(123)
fit_full_rpart <- rpart(Resp ~ ., data=training, method="class")
print(fit_full_rpart)
rpart.plot(fit_full_rpart)
fancyRpartPlot(fit_full_rpart)

rsq.rpart(fit_full_rpart)

tmp <- printcp(fit_full_rpart) 

rsq.val <- 1-tmp[,c(3,4)] ; rsq.val
rsq.val.1 <- 1-tmp[,c(3,4,5)] ; rsq.val.1
rsq.val.1[,3] + rsq.val.1[,1] < rsq.val.1[,2]

plotcp(fit_full_rpart)

##########################################
# You rely on this for error or accuracy
##########################################

summary(tree(fit_full_rpart))

# You need to specify type=class (as a factor I think)
fit_full_rpart_test <- predict(fit_full_rpart, newdata=testing, type="class")
print(fit_full_rpart_test)
summary(fit_full_rpart_test)

cm_full_rpart <- confusionMatrix(fit_full_rpart_test, testing$Resp)
cm_full_rpart

(216 + 17) / 299 = 0.7792642


# I don't think we will need this, it is the same as the table

plot(cm_full_rpart$table, col = cm_full_rpart$byClass, 
     main = paste("Decision Tree Conusion Matrix: Accuray =", 
                  round(cm_full_rpart$overall['Accuray'], 4)))


testing$Resp_fit_full_rpart_test <- 0;


fit_full_rpart_test[,2]


for (i in 1:nrow(testing)) {
    if (testing[i]$row.names == fit_full_rpart_test[i]$row.names) {
        
    }    
}







table(testing$Resp, fit_full_glm_2_best_test > 0.5)



#########################################
# It may not produce better results
#########################################


fit_sum_rpart <- rpart(Resp ~ alanine + arginine + asparagine + asparticAcid + 
                           cysteine + glutamine + glutaminAcid + glycine + 
                           histidine + isoleucine + leucine + lysine + 
                           methionine + phenylalanine + proline + serine + 
                           threonine + tryptophan + tyrosine + valine + 
                           stopCodon + startCodon + totalSeq + totalCondon, 
                       data=training, method="class")
print(fit_sum_rpart)


rsq.rpart(fit_sum_rpart)


tmp <- printcp(fit_sum_rpart) 

rsq.val <- 1-tmp[,c(3,4)] ; rsq.val
rsq.val.1 <- 1-tmp[,c(3,4,5)] ; rsq.val.1
rsq.val.1[,3] + rsq.val.1[,1] < rsq.val.1[,2]

plotcp(fit_sum_rpart)
install.packages("tree")
library(tree)
summary(tree(fit_sum_rpart))


rpart.plot(fit_sum_rpart)
fancyRpartPlot(fit_sum_rpart)



## -----------------------------
#        Random Forest
## -----------------------------

set.seed(11234)

# This seems to have taken a long time to execute this.

fit_full_rm_0 <- randomForest(Resp ~ . , data=training, family="binomial")
fit_ful_rm_0
summary(fit_full_rm_0)

# Alternatively, use caret package for this
library(caret)
fit_full_rm <- train(Resp ~. , data=training, method="rf", prox=TRUE)
fit_full_rm
summary(fit_full_rm)

getTree(fit_full_rm$finalModel, k=2)


# ---- Predict Random Forest ------
fit_full_rm_0_test <- predict(fit_full_rm_0, testing, type="class")
cm_rm_0 <- confusionMatrix(fit_full_rm_0_test, testing$Resp)
cm_rm_0

# Accuracy = 0.7826

fit_full_rm_test <- predict(fit_full_rm, testing)
cm_rm <- confusionMatrix(fit_full_rm_test, testing$Resp)
cm_rm

# Accurary = 0.7826

plot(fit_full_rm_0)
getTree(fit_full_rm_0, 1, labelVar=TRUE)


# -------------------------
# Boosted Regression
# -------------------------
fit_full_gbm <- train(Resp ~. , data=training, method="gbm", verbose=FALSE)
print(fit_full_gbm)

fit_full_gbm_test <- predict(fit_full_gbm, newdata=testing)

qplot(fit_full_gbm_test, Resp, data=testing)

cm_gbm <- confusionMatrix(fit_full_gbm_test, testing$Resp)
cm_gbm

# Accuracy = 0.786

plot(cm_gbm)
