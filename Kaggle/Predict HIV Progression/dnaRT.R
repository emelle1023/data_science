dnaRT <- function(hiv) {
    
    for (i in 1:nrow(hiv)) {    
        c <- as.character(hiv$RT.Seq[i]); 
        print(nchar(c))

        hiv$totalSeqRT[i] <- as.integer(nchar(c))
        hiv$totalCodonRT[i] <- as.integer(nchar(c)) / 3

        for (j in 1:nchar(c)) {
            if (nchar(c) >= 3) {
                m <- substring(c, 1, 3)

                if (m %in% alaA) {
                    alanineRT <- alanineRT + 1; 
                    hiv$alanineRT[i] <- alanineRT
                } else if (m %in% argR) {
                    arginineRT <- arginineRT + 1; 
                    hiv$arginineRT[i] <- arginineRT + 1
                } else if (m %in% asnN) {                  
                    asparagineRT <- asparagineRT + 1; 
                    hiv$asparagineRT[i] <- asparagineRT
                } else if (m %in% aspD) {
                    asparticAcidRT <- asparticAcidRT + 1;
                    hiv$asparticAcidRT[i] <- asparticAcidRT
                } else if (m %in% cysC) {
                    cysteineRT <- cysteineRT + 1;
                    hiv$cysteineRT[i] <- cysteineRT                    
                } else if (m %in% ginQ) {
                    glutamineRT <- glutamineRT + 1;
                    hiv$glutamineRT[i] <- glutamineRT
                } else if (m %in% gluE) {
                    glutaminAcidRT <- glutaminAcidRT + 1;
                    hiv$glutaminAcidRT[i] <- glutaminAcidRT
                } else if (m %in% glyG) {
                    glycineRT <- glycineRT + 1;
                    hiv$glycineRT[i] <- glycineRT
                } else if (m %in% hisH) {
                    histidineRT <- histidineRT + 1;
                    hiv$histidineRT[i] <- histidineRT
                } else if (m %in% iieI) {
                    isoleucineRT <- isoleucineRT + 1;
                    hiv$isoleucineRT[i] <- isoleucineRT
                } else if (m %in% leuL) {
                    leucineRT <- leucineRT + 1;
                    hiv$leucineRT[i] <- leucineRT
                } else if (m %in% lysK) {
                    lysineRT <- lysineRT + 1;
                    hiv$lysineRT[i] <- lysineRT
                } else if (m %in% metM) {
                    methionineRT <- methionineRT + 1;
                    hiv$methionineRT[i] <- methionineRT
                } else if (m %in% pheF) {
                    phenylalanineRT <- phenylalanineRT + 1;
                    hiv$phenylalanineRT[i] <- phenylalanineRT
                } else if (m %in% proP) {
                    prolineRT <- prolineRT + 1;
                    hiv$prolineRT[i] <- prolineRT
                } else if (m %in% serS) {
                    serineRT <- serineRT + 1;
                    hiv$serineRT[i] <- serineRT
                } else if (m %in% thrT) {
                    threonineRT <- threonineRT + 1;
                    hiv$threonineRT[i] <- threonineRT
                } else if (m %in% trpW) {
                    tryptophanRT <- tryptophanRT + 1;
                    hiv$tryptophanRT[i] <- tryptophanRT
                } else if (m %in% tyrY) {
                    tyrosineRT <- tyrosineRT + 1;
                    hiv$tyrosineRT[i] <- tyrosineRT
                } else if (m %in% valV) {
                    valineRT <- valineRT + 1;
                    hiv$valineRT[i] <- valineRT
                } else if (m %in% startC) {
                    startCodonRT <- startCodonRT + 1;
                    hiv$startCodonRT[i] <- startCodonRT
                } else if (m %in% stopC) {
                    stopCodonRT <- stopCodonRT + 1;
                    hiv$stopCodonRT[i] <- stopCodonRT
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
        
        alanineRT <- 0; arginineRT <- 0; asparagineRT <- 0;  
        asparticAcidRT <- 0; cysteineRT <- 0; glutamineRT <- 0;
        glutaminAcidRT <- 0; glycineRT <- 0; histidineRT <- 0; 
        isoleucineRT <- 0; leucineRT <- 0; lysineRT <- 0; 
        methionineRT <- 0; phenylalanineRT <- 0; prolineRT <- 0; 
        serineRT <- 0; threonineRT <- 0; tryptophanRT <- 0;
        tyrosineRT <- 0; valineRT <- 0; startCodonRT <- 0; 
        stopCodonRT <- 0;
    }
    
    return (hiv)
}