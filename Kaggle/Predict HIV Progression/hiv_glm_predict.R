setwd("~/Documents/My Projects/Predict HIV Progression")

hiv_test <- read.csv("test_data.csv", header=TRUE)

hiv_test$Resp <- factor(hiv_test$Resp)

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


hiv_test$alanineRT <- 0; hiv_test$arginineRT <- 0;  hiv_test$asparagineRT <- 0; 
hiv_test$asparticAcidRT <- 0; hiv_test$cysteineRT <- 0; hiv_test$glutamineRT <- 0;
hiv_test$glutaminAcidRT <- 0; hiv_test$glycineRT <- 0;hiv_test$histidineRT <- 0; 
hiv_test$isoleucineRT <- 0; hiv_test$leucineRT <- 0; hiv_test$lysineRT <- 0; 
hiv_test$methionineRT <- 0; hiv_test$phenylalanineRT <- 0; hiv_test$prolineRT <- 0; 
hiv_test$serineRT <- 0; hiv_test$threonineRT <- 0; hiv_test$tryptophanRT <- 0; 
hiv_test$tyrosineRT <- 0; hiv_test$valineRT <- 0; hiv_test$startCodonRT <- 0; 
hiv_test$stopCodonRT <- 0;
hiv_test$alaninePR <- 0; hiv_test$argininePR <- 0; hiv_test$asparaginePR <- 0; 
hiv_test$asparticAcidPR <- 0; hiv_test$cysteinePR <- 0; hiv_test$glutaminePR <- 0;
hiv_test$glutaminAcidPR <- 0; hiv_test$glycinePR <- 0;hiv_test$histidinePR <- 0; 
hiv_test$isoleucinePR <- 0; hiv_test$leucinePR <- 0; hiv_test$lysinePR <- 0; 
hiv_test$methioninePR <- 0;  hiv_test$phenylalaninePR <- 0; hiv_test$prolinePR <- 0; 
hiv_test$serinePR <- 0;  hiv_test$threoninePR <- 0; hiv_test$tryptophanPR <- 0; 
hiv_test$tyrosinePR <- 0;  hiv_test$valinePR <- 0; hiv_test$startCodonPR <- 0; 
hiv_test$stopCodonPR <- 0;
hiv_test$totalSeqRT <- 0; hiv_test$totalCodonRT <- 0;
hiv_test$totalSeqPR <- 0; hiv_test$totalCodonPR <- 0;

source("dnaRT.R")
hiv_test1 <- dnaRT(hiv_test)
write.table(hiv_test1, "hiv_test1.csv", sep=",")

source("dnaPR.R")
hiv_test2 <- dnaPR(hiv_test1)
write.table(hiv_test2, "hiv_test2.csv", sep=",")


hiv_test2$alanine <- hiv_test2$alanineRT + hiv_test2$alaninePR
hiv_test2$arginine <- hiv_test2$arginineRT + hiv_test2$argininePR
hiv_test2$asparagine  <- hiv_test2$asparagineRT + hiv_test2$asparaginePR
hiv_test2$asparticAcid <- hiv_test2$asparticAcidRT + hiv_test2$asparticAcidPR
hiv_test2$cysteine <- hiv_test2$cysteineRT + hiv_test2$cysteinePR
hiv_test2$glutamine <- hiv_test2$glutamineRT + hiv_test2$glutaminePR
hiv_test2$glutaminAcid <- hiv_test2$glutaminAcidRT + hiv_test2$glutaminAcidPR
hiv_test2$glycine <- hiv_test2$glycineRT + hiv_test2$glycinePR 
hiv_test2$histidine <- hiv_test2$histidineRT + hiv_test2$histidinePR 
hiv_test2$isoleucine <- hiv_test2$isoleucineRT + hiv_test2$isoleucinePR 
hiv_test2$leucine <- hiv_test2$leucineRT + hiv_test2$leucinePR 
hiv_test2$lysine <- hiv_test2$lysineRT + hiv_test2$lysinePR 
hiv_test2$methionine <- hiv_test2$methionineRT + hiv_test2$methioninePR 
hiv_test2$phenylalanine <- hiv_test2$phenylalanineRT + hiv_test2$phenylalaninePR 
hiv_test2$proline <- hiv_test2$prolineRT + hiv_test2$prolinePR 
hiv_test2$serine <- hiv_test2$serineRT + hiv_test2$serinePR 
hiv_test2$threonine <- hiv_test2$threonineRT + hiv_test2$threoninePR 
hiv_test2$tryptophan <- hiv_test2$tryptophanRT + hiv_test2$tryptophanPR 
hiv_test2$tyrosine <- hiv_test2$tyrosineRT + hiv_test2$tyrosinePR 
hiv_test2$valine <- hiv_test2$valineRT + hiv_test2$valinePR 
hiv_test2$startCodon <- hiv_test2$startCodonRT + hiv_test2$startCodonPR 
hiv_test2$stopCodon <- hiv_test2$stopCodonRT + hiv_test2$stopCodonPR
hiv_test2$totalSeq <- hiv_test2$totalSeqRT + hiv_test2$totalSeqPR
hiv_test2$totalCondon <- hiv_test2$totalCodonRT + hiv_test2$totalCodonPR


hiv_test3 <- hiv_test2

hiv_test3$PatientID <- NULL
hiv_test3$PR.Seq <- NULL
hiv_test3$RT.Seq <- NULL
# You remove response outcome
hiv_test3$Resp <- NULL


pred <- predict(fit_full_best, newdata=hiv_test3, type='response')
sum(as.integer(pred * 100) > 50)
sum(as.numeric(pred * 100) > 50)
sum((pred* 100) > 50)
sum((pred* 100) > 75)


hiv_test3$Resp <- 0 # Initial Resp
pr <- 85
for (i in 1:nrow(hiv_test3)) {
    if ((pred[i]*100) > pr) {
        hiv_test3$Resp[i] <- 1
    }
}

sum(hiv_test3$Resp)

# summary(pred)

