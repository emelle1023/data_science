dnaPR <- function(hiv) {
    
    for (i in 1:nrow(hiv)) {        
        c <- as.character(hiv$PR.Seq[i]); 
        print(nchar(c))

        hiv$totalSeqPR[i] <- as.integer(nchar(c))
        hiv$totalCodonPR[i] <- as.integer(nchar(c)) /3

        for (j in 1:nchar(c)) {
            if (nchar(c) >= 3) {
                m <- substring(c, 1, 3)
                
                if (m %in% alaA) {
                    alaninePR <- alaninePR + 1; 
                    hiv$alaninePR[i] <- alaninePR
                } else if (m %in% argR) {
                    argininePR <- argininePR + 1; 
                    hiv$argininePR[i] <- argininePR + 1
                } else if (m %in% asnN) {                  
                    asparaginePR <- asparaginePR + 1; 
                    hiv$asparaginePR[i] <- asparaginePR
                } else if (m %in% aspD) {
                    asparticAcidPR <- asparticAcidPR + 1;
                    hiv$asparticAcidPR[i] <- asparticAcidPR
                } else if (m %in% cysC) {
                    cysteinePR <- cysteinePR + 1;
                    hiv$cysteinePR[i] <- cysteinePR                    
                } else if (m %in% ginQ) {
                    glutaminePR <- glutaminePR + 1;
                    hiv$glutaminePR[i] <- glutaminePR
                } else if (m %in% gluE) {
                    glutaminAcidPR <- glutaminAcidPR + 1;
                    hiv$glutaminAcidPR[i] <- glutaminAcidPR
                } else if (m %in% glyG) {
                    glycinePR <- glycinePR + 1;
                    hiv$glycinePR[i] <- glycinePR
                } else if (m %in% hisH) {
                    histidinePR <- histidinePR + 1;
                    hiv$histidinePR[i] <- histidinePR
                } else if (m %in% iieI) {
                    isoleucinePR <- isoleucinePR + 1;
                    hiv$isoleucinePR[i] <- isoleucinePR
                } else if (m %in% leuL) {
                    leucinePR <- leucinePR + 1;
                    hiv$leucinePR[i] <- leucinePR
                } else if (m %in% lysK) {
                    lysinePR <- lysinePR + 1;
                    hiv$lysinePR[i] <- lysinePR
                } else if (m %in% metM) {
                    methioninePR <- methioninePR + 1;
                    hiv$methioninePR[i] <- methioninePR
                } else if (m %in% pheF) {
                    phenylalaninePR <- phenylalaninePR + 1;
                    hiv$phenylalaninePR[i] <- phenylalaninePR
                } else if (m %in% proP) {
                    prolinePR <- prolinePR + 1;
                    hiv$prolinePR[i] <- prolinePR
                } else if (m %in% serS) {
                    serinePR <- serinePR + 1;
                    hiv$serinePR[i] <- serinePR
                } else if (m %in% thrT) {
                    threoninePR <- threoninePR + 1;
                    hiv$threoninePR[i] <- threoninePR
                } else if (m %in% trpW) {
                    tryptophanPR <- tryptophanPR + 1;
                    hiv$tryptophanPR[i] <- tryptophanPR
                } else if (m %in% tyrY) {
                    tyrosinePR <- tyrosinePR + 1;
                    hiv$tyrosinePR[i] <- tyrosinePR
                } else if (m %in% valV) {
                    valinePR <- valinePR + 1;
                    hiv$valinePR[i] <- valinePR
                } else if (m %in% startC) {
                    startCodonPR <- startCodonPR + 1;
                    hiv$startCodonPR[i] <- startCodonPR
                } else if (m %in% stopC) {
                    stopCodonPR <- stopCodonPR + 1;
                    hiv$stopCodonPR[i] <- stopCodonPR
                } else if (m == "") {
                    print("wrong")
                } else {
                    print(paste0("(patent=", i, "), (seq=", j, "), (codon=",m, ")"))
                }
            } else {
                break
            }
            c <- substring(c, 4)
        }    
        alaninePR <- 0; argininePR <- 0; asparaginePR <- 0; 
        asparticAcidPR <- 0; cysteinePR <- 0; glutaminePR <- 0;
        glutaminAcidPR <- 0; glycinePR <- 0; histidinePR <- 0; 
        isoleucinePR <- 0; leucinePR <- 0; lysinePR <- 0; 
        methioninePR <- 0; phenylalaninePR <- 0; prolinePR <- 0; 
        serinePR <- 0; threoninePR <- 0; tryptophanPR <- 0;
        tyrosinePR <- 0; valinePR <- 0; startCodonPR <- 0; 
        stopCodonPR <- 0;
                
    }
    return (hiv)
}