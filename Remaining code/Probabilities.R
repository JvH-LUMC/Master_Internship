library(foreign, quietly = TRUE)
library(rms, verbose = FALSE, quietly = TRUE)
library(mctest, verbose = FALSE, quietly = TRUE)
library(mice, quietly = TRUE)
library(glmnet,quietly = TRUE)
library(glmpath,quietly = TRUE)

part1 <- read.spss("~/Documents/OneDrive/Voorbereiding R coding/NewCode/Part1.sav", to.data.frame = TRUE, use.value.labels = TRUE) 
part2 <- read.spss("~/Documents/OneDrive/Voorbereiding R coding/NewCode/Part2.sav", to.data.frame = TRUE, use.value.labels = TRUE)

part1_imp <- mice(part1, m=10, seed=1, print=FALSE) 
part2_imp <- mice(part2, m=10, seed=1, print=FALSE) 

part1_imp1 <- complete(part1_imp, 1)
part1_imp10 <- complete(part1_imp, 1)
for(i in 2:10) {
  part1_imp10 <- rbind(part1_imp10, complete(part1_imp, i))
}
part1_imp1$w <- 1/10
part1_imp10$w <- 1/10

part2_imp1 <- complete(part2_imp, 1)
part2_imp10 <- complete(part2_imp, 1)
for(i in 2:10) {
  part2_imp10 <- rbind(part2_imp10, complete(part2_imp, i))
}
part2_imp1$w <- 1/10
part2_imp10$w <- 1/10

# Loop for predicted probability per row
Prob<- matrix(nrow=887, ncol=8) 
colnames(Prob) <- c("PredProb", "Mc_FailedFemoralApproach", "ArchElongation", "AorticVariant", "InnCca_90", "AngleBif", "NASCET", "Total")  
for (i in 1:dim(part1_imp1)[1]) {
  int <- -3.3512221
  ct1 <- 1.2490924
  ct2 <- 0.1560974
  ct3 <- 0.4053100
  ct4 <- 0.3570322
  ct5 <- 0.3961826
  p.ct1 <- as.numeric(part1_imp1$ArchElongation_dichIII[i])
  p.ct2 <- as.numeric(part1_imp1$AorticVariant_dichC[i])
  p.ct3 <- as.numeric(part1_imp1$InnCca_90orLarger[i])
  p.ct4 <- as.numeric(part1_imp1$AngleFollowingBifurcation[i])
  p.ct5 <- as.numeric(part1_imp1$ICAE_NASCET_70[i])
  
  formula <- int + (p.ct1*ct1) + (p.ct2*ct2) + (p.ct3*ct3) + (p.ct4*ct4) + (p.ct5*ct5)
  
  probability <- exp(formula)/(1+(exp(formula)))
  round(probability, 3)
  Prob[i,1] <-  round(probability*100, 3)
  Prob[i,2] <- part1_imp1$Mc_FailedFemoralApproach[i]
  Prob[i,3] <- p.ct1
  Prob[i,4] <- p.ct2
  Prob[i,5] <- p.ct3
  Prob[i,6] <- p.ct4
  Prob[i,7] <- p.ct5
  Prob[i,8] <- sum(p.ct1+p.ct2+p.ct3+p.ct4+p.ct5)
}
Prob <- as.data.frame(Prob)

# Calculate frequency per predicted probability interval 
sum(Prob[,1]<10)
sum(Prob[,1]>=10 & Prob[,1]<20)
sum(Prob[,1]>=20 & Prob[,1]<30)
sum(Prob[,1]>=30 & Prob[,1]<40)
sum(Prob[,1]>=40 & Prob[,1]<50)
sum(Prob[,1]>=50 & Prob[,1]<60)
sum(Prob[,1]>=60 & Prob[,1]<70)
sum(Prob[,1]>=70 & Prob[,1]<80)
sum(Prob[,1]>=80 & Prob[,1]<90)
sum(Prob[,1]>=90 & Prob[,1]<101)

# Calculate frequency per how many charachteristics in total 
sum(Prob[,8]<=0)
sum(Prob[,8]>=1 & Prob[,8]<2)
sum(Prob[,8]>=2 & Prob[,8]<3)
sum(Prob[,8]>=3 & Prob[,8]<4)
sum(Prob[,8]>=4 & Prob[,8]<5)
sum(Prob[,8]>=5)


install.packages("writexl")
library(writexl)
write_xlsx(Prob, "~/Desktop/prob.xlsx")


# Now for failed femoral approach patients
sum(Prob[,8]<=0 & Prob[,2]>0)
sum(Prob[,8]>=1 & Prob[,8]<2 & Prob[,2]>0)
sum(Prob[,8]>=2 & Prob[,8]<3 & Prob[,2]>0)
sum(Prob[,8]>=3 & Prob[,8]<4 & Prob[,2]>0)
sum(Prob[,8]>=4 & Prob[,8]<5 & Prob[,2]>0)
sum(Prob[,8]>=5 & Prob[,2]>0)