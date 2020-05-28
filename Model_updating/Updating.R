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

# external validation 
mod02 <- lrm(Mc_FailedFemoralApproach ~ ArchElongation_dichIII + AorticVariant_dichC + InnCca_90orLarger + AngleFollowingBifurcation +ICAE_NASCET_70, data = part2_imp10, weight=w, x=T, y=T)
coefDEF <- c(-4.3125508, 1.6299069, 0.9145084, 1.0513471, 0.7295294, 1.2001544)  
lp_10 <- mod02$x %*% coefDEF[-1]
ext <- lrm(Mc_FailedFemoralApproach ~ lp_10, data = part2_imp10, weights = w, x=T, y=T)
validate(ext, B=200, bw=F, type="individual", dxy=T)
rcorr.cens(ext$linear.predictors, ext$y)
boot <- rcorr.cens(ext$linear.predictors, ext$y)
boot[1] -1.96*boot[3]/2
boot[1] +1.96*boot[3]/2

mod02$coefficients
ext$coefficients
# Search for mean lambda over each set
m <- 10 
Results.m <- matrix(nrow=10, ncol=5) 
colnames(Results.m) <- c("i", "lambda","shrinkage","df", "Apparent c" )  
for(i3 in 1:m) {
  data1 <- as.data.frame(complete(part2_imp, i3)) 
  data1$ArchElongation_dichIII <- ifelse(data1$ArchElongation_dichIII == "I or II", 0,1)
  data1$AorticVariant_dichC <- ifelse(data1$AorticVariant_dichC == "A or B", 0,1)
  data1$InnCca_90orLarger <- ifelse(data1$InnCca_90orLarger == "0", 0,1)
  data1$AngleFollowingBifurcation <- ifelse(data1$AngleFollowingBifurcation == "no", 0,1)
  data1$ICAE_NASCET_70 <- ifelse(data1$ICAE_NASCET_70 == "<70", 0,1)
  
  mod02 <- lrm(Mc_FailedFemoralApproach ~ ArchElongation_dichIII + AorticVariant_dichC + InnCca_90orLarger + AngleFollowingBifurcation +ICAE_NASCET_70, data = data1, x=T, y=T)
  coefDEF <- c(-4.3125508, 1.6299069, 0.9145084, 1.0513471, 0.7295294, 1.2001544)  
  lp <- mod02$x %*% coefDEF[-1] 
  lp <- lp*-1
  y2 <- mod02$y
  x2 <- model.matrix(Mc_FailedFemoralApproach ~ ArchElongation_dichIII +  AorticVariant_dichC + InnCca_90orLarger 
                     + AngleFollowingBifurcation + ICAE_NASCET_70, data = data1)[,-1]
  mod02.2 <- lrm(Mc_FailedFemoralApproach ~ ArchElongation_dichIII + AorticVariant_dichC + InnCca_90orLarger + AngleFollowingBifurcation +ICAE_NASCET_70, offset(lp), data = data1, x=T, y=T)
  
  mod01.1 <- lrm(formula = Mc_FailedFemoralApproach ~ ArchElongation_dichIII + AorticVariant_dichC + InnCca_90orLarger + AngleFollowingBifurcation +ICAE_NASCET_70, data = data1, x=T, y=T) 
  x <- model.matrix(Mc_FailedFemoralApproach ~ ArchElongation_dichIII + AorticVariant_dichC + InnCca_90orLarger + AngleFollowingBifurcation +ICAE_NASCET_70, data = data1)[,-1]
  
  set.seed(123)
  cv.glmmod <- cv.glmnet(mod02.2$x, mod02.2$y, alpha=1, family = "binomial") 
  model.L1 <- glmnet(mod02.2$x, mod02.2$y, alpha=1, family = "binomial", lambda = cv.glmmod$lambda.min) 
  lp1 <- x2 %*% coef(model.L1)[-1]
  Results.m[i3,c(1,2)] <- c(i3, cv.glmmod$lambda.min) 
  shrinkage <- 1/coef(lrm(y2~ lp1))
  Results.m[i3,3] <-  shrinkage[2]
  Results.m[i3,4] <- sum(abs(coef(model.L1)) > 0) 
  Results.m[i3,5] <- rcorr.cens(as.numeric(lp1),y2)[2] / 2 + .5 
}
Results.m
LASSO.10imp <- apply(Results.m, 2, mean) 
LASSO.10imp

part2_imp10$ArchElongation_dichIII <- ifelse(part2_imp10$ArchElongation_dichIII == "I or II", 0,1)
part2_imp10$AorticVariant_dichC <- ifelse(part2_imp10$AorticVariant_dichC == "A or B", 0,1)
part2_imp10$InnCca_90orLarger <- ifelse(part2_imp10$InnCca_90orLarger == "0", 0,1)
part2_imp10$AngleFollowingBifurcation <- ifelse(part2_imp10$AngleFollowingBifurcation == "no", 0,1)
part2_imp10$ICAE_NASCET_70 <- ifelse(part2_imp10$ICAE_NASCET_70 == "<70", 0,1)

# Use mean lambda for the lasso model 
mod02DEF <- lrm(Mc_FailedFemoralApproach ~ ArchElongation_dichIII + AorticVariant_dichC + InnCca_90orLarger + AngleFollowingBifurcation +ICAE_NASCET_70, offset(lp_10.2), data = part2_imp10, weight=w, x=T, y=T)
lp_10.2 <- lp_10*-1
mod02DEF$coefficients
mod02$coefficients

# Validate the model before the lasso 
validate(mod02DEF, B=200, bw=F, type="individual", dxy=T) 
rcorr.cens(mod02DEF$linear.predictors, mod02DEF$y)
boot <- rcorr.cens(mod02DEF$linear.predictors, mod02DEF$y)
boot[1] -1.96*boot[3]/2
boot[1] +1.96*boot[3]/2

set.seed(123)
lassomodel <- glmnet(x=mod02DEF$x, y=mod02DEF$y, weights =part2_imp10$w, alpha=1, lambda = LASSO.10imp[2], family="binomial") 
coef(lassomodel)
round(coef(lassomodel)/mod02DEF$coef,3) 
mean(coef(lassomodel)/mod02DEF$coef)
lassocoef <- c(Intercept = -3.3512221, ArchElongation_dichIII =  1.2490924, AorticVariant_dichC = 0.1560974, InnCca_90orLarger = 0.4053100, 
               AngleFollowingBifurcation = 0.3570322, ICAE_NASCET_70 = 0.3961826)
mod02DEF$coefficients <- lassocoef

# Validate the model after the lasso 
validate(mod02DEF, B=200, bw=F, type="individual", dxy=T) 
rcorr.cens(mod02DEF$linear.predictors, mod02DEF$y)
boot <- rcorr.cens(mod02DEF$linear.predictors, mod02DEF$y)
boot[1] -1.96*boot[3]/2
boot[1] +1.96*boot[3]/2

# Calculate odds ratio's
#dd <- datadist(part1_imp10)
#options(datadist="dd")
p1 <- glm(Mc_FailedFemoralApproach ~ ICAE_NASCET_70, data = part1_imp10, weights = w, x=T, y=T, family = binomial(link = "logit"))
exp(cbind(coef(p1), confint(p1)))
