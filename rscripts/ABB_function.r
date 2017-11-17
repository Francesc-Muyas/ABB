suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(gamlss))

####### FUNCTIONS
roc.curve <- function(s,TABLE=MERGE_test){
  #S <- predict(ALL_SC,type=MERGE_test,newdata=TABLE)
  Ps <- (S >= s)*1
  
  FP <- as.numeric(sum((Ps==1)*(TABLE$Y==0)))
  TP <- as.numeric(sum((Ps==1)*(TABLE$Y==1)))
  TN <- as.numeric(sum((Ps==0)*(TABLE$Y==0)))
  FN <- as.numeric(sum((Ps==0)*(TABLE$Y==1)))
  
  PPV <- TP/(TP+FP)
  
  return(PPV)
}

BINOM_TEST <- function(X,Y,P){
  #print (c(X,Y,P))
  if (X/Y > P){
    P_VAL <- 2*(1-pbinom(X-1,Y,P,lower.tail = T))
  }
  if (X/Y < P){
    P_VAL <- 2*(1-pbinom(X,Y,P,lower.tail = F))
  }
  if (X/Y == P){
    #P_VAL <- 2*(1-pbinom(X,Y,P,lower.tail = F))-dbinom(X,Y,P)
    P_VAL <- 1
  }
return(P_VAL)
}


ABB_func <- function(SCORE_1,SCORE_2,SCORE_3){
  # Getting glm value
  POS <- data.frame(SCORE_1,SCORE_2,SCORE_3)
  colnames(POS) <- c("SCORE_1", "SCORE_2", "SCORE_3")
  logABB <- predict(ALL_SC,newdata=POS,type="response")
  
  # Getting ABB value
  ABB <- round(min(roc.curve(logABB,MERGE_test)),digits=4)
  return(as.numeric(ABB))
}

DIPLOID_func <- function(ALT,DP,alpha){
  
  ####### PARAMETERS 0/0 GT
  MU_0 <- 0.03270537
  SIGMA_0 <- 0.1452151
  NU_0 <- 1.689303e+01
  TAU_0 <- 1.011730e-14
  
  ####### PARAMETERS 1/1 GT  
  MU_1 <- 0.97263928
  SIGMA_1 <-  0.1364325
  NU_1 <- 2.719154e-15
  TAU_1 <- 4.626247e+00
  
  ####### PARAMETERS 0/1 GT
  P <- 0.5
  
  
  ####### INFO
  #ALT <- c(0,1,10,9,19,NA)
  ALT2 <- ALT[!(is.na(ALT))]
  
  #DP <- c(10,9,20,17,20,NA)
  DP2 <- DP[!(is.na(DP))]
  
  AB <- ALT2/DP2
  AB2 <- (ALT2-1)/DP2
  #print (sum(is.na(AB)*1))
  
  ## P-values for GT 0/0
  P_HOM0 <- as.numeric(lapply(AB2,function(x) ifelse(x < 0, 1, pBEINF(x, MU_0, SIGMA_0, NU_0, TAU_0, lower.tail=FALSE))))

  ## P-values for GT 0/1
  P_heter <- mapply(function(x, y) BINOM_TEST(x,y,P), (ALT2), (DP2))

  ## P-values for GT 1/1
  P_HOM1 <- pBEINF(AB, MU_1, SIGMA_1, NU_1, TAU_1, lower.tail=T)
  
  ## Getting the greatest P value
  HIGH_P <- mapply(function(x, y, z) max(x,y,z), (P_heter), (P_HOM0),(P_HOM1))
  
  ## Calculate the absolute distance taking into account the distribution which each position-sample fits better
  DISTANCE <- mapply(function(x, y, z, w) ifelse(is.na(x) | is.na(y)| is.na(z), NA, if (max(x,y,z) == x){abs(w-P)}else if(max(x,y,z) == y){w}else if(max(x,y,z) == z){abs(w-1)}), (P_heter), (P_HOM0),t(P_HOM1),AB)
  
  ## Matrix with binari results -> 1: significant; 0:non-significant. Taking alpha as threshold
  SIG <- mapply(function(x, y, z) ifelse(is.na(x) | is.na(y)| is.na(z), NA, if(x < alpha/2 & y < alpha & z < alpha){1}else{0}),t(P_heter), t(P_HOM0),t(P_HOM1))
  
  ## SCORE 1 -> mean of the absolute distance to the ideal distribution per position
  SCORE_1 <- mean(DISTANCE)
  
  ## SCORE 2 -> proportion of significant variants in each position in the analyzed samples
  SCORE_2 <- sum(SIG)/length(SIG)
  
  ## SCORE 3 -> mean of highest -log10(probabilities) per position 
  SCORE_3 <- mean(-log10(HIGH_P))
  
  ## Getting ABB score value
  ABB_SCORE <- ABB_func(SCORE_1,SCORE_2,SCORE_3)
  #ABB_SCORE <- 1
  return(c(SCORE_1, SCORE_2, SCORE_3,ABB_SCORE))
}

HAPLOID_func <- function(ALT,DP, alpha){
  
  ####### PARAMETERS 0/0 GT
  MU_0 <- 0.03270537
  SIGMA_0 <- 0.1452151
  NU_0 <- 1.689303e+01
  TAU_0 <- 1.011730e-14
  
  ####### PARAMETERS 1/1 GT  
  MU_1 <- 0.97263928
  SIGMA_1 <-  0.1364325
  NU_1 <- 2.719154e-15
  TAU_1 <- 4.626247e+00
  
  ####### PARAMETERS 0/1 GT
  P <- 0
  
  
  ####### INFO
  #ALT <- c(0,1,10,9,19,NA)
  ALT2 <- ALT[!(is.na(ALT))]
  
  #DP <- c(10,9,20,17,20,NA)
  DP2 <- DP[!(is.na(DP))]
  
  AB <- ALT2/DP2
  AB2 <- (ALT2-1)/DP2
  
  ## P-values for GT 0/0
  P_HOM0 <- as.numeric(lapply(AB2,function(x) ifelse(x < 0, 1, pBEINF(x, MU_0, SIGMA_0, NU_0, TAU_0, lower.tail=FALSE))))
  
  ## P-values for GT 0/1
  P_heter <- 0
    
  ## P-values for GT 1/1
  P_HOM1 <- pBEINF(AB, MU_1, SIGMA_1, NU_1, TAU_1, lower.tail=T)
  
  ## Getting the greatest P value
  HIGH_P <- mapply(function(x, y, z) max(x,y,z), (P_heter), (P_HOM0),(P_HOM1))
  
  ## Calculate the absolute distance taking into account the distribution which each position-sample fits better
  DISTANCE <- mapply(function(x, y, z, w) ifelse(is.na(x) | is.na(y)| is.na(z), NA, if (max(x,y,z) == x){abs(w-P)}else if(max(x,y,z) == y){w}else if(max(x,y,z) == z){abs(w-1)}), (P_heter), (P_HOM0),t(P_HOM1),AB)
  
  ## Matrix with binari results -> 1: significant; 0:non-significant. Taking alpha as threshold
  SIG <- mapply(function(x, y, z) ifelse(is.na(x) | is.na(y)| is.na(z), NA, if(x < alpha/2 & y < alpha & z < alpha){1}else{0}),t(P_heter), t(P_HOM0),t(P_HOM1))
  
  ## SCORE 1 -> mean of the absolute distance to the ideal distribution per position
  SCORE_1 <- mean(DISTANCE)
  
  ## SCORE 2 -> proportion of significant variants in each position in the analyzed samples
  SCORE_2 <- sum(SIG)/length(SIG)
  
  ## SCORE 3 -> mean of highest -log10(probabilities) per position 
  SCORE_3 <- mean(-log10(HIGH_P))
  
  ## Getting ABB score value
  ABB_SCORE <- ABB_func(SCORE_1,SCORE_2,SCORE_3)
  return(c(SCORE_1, SCORE_2, SCORE_3,ABB_SCORE))
}



