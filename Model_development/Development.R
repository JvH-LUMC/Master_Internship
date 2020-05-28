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

# Search for mean lambda over each imputed set
m <- 10 
Results.m <- matrix(nrow=10, ncol=5) 
colnames(Results.m) <- c("i", "lambda","shrinkage","df", "Apparent c" )  
for(i in 1:m) {
  data1 <- as.data.frame(complete(part1_imp, i))
  data1$ArchElongation_dichIII <- ifelse(data1$ArchElongation_dichIII == "Type I or II", 0,1)
  data1$AorticVariant_dichC <- ifelse(data1$AorticVariant_dichC == "Type A, B, D, other", 0,1)
  data1$InnCca_90orLarger <- ifelse(data1$InnCca_90orLarger == "1", 0,1)
  data1$AngleFollowingBifurcation <- ifelse(data1$AngleFollowingBifurcation == "No angle >= 90 degrees", 0,1)
  data1$ICAE_NASCET_70 <- ifelse(data1$ICAE_NASCET_70 == "1", 0,1)
  data1$M_prev_ht <- ifelse(data1$M_prev_ht == "No", 0,1)
  mod01.1 <- lrm(formula = Mc_FailedFemoralApproach ~ M_age + ArchElongation_dichIII + AngleAaInn_OR_AaCca_dich45 + AorticVariant_dichB + AorticVariant_dichC + InnCca_90orLarger 
                 + ICA_90orLarger + AngleFollowingBifurcation + OrigoStenosis50percent + ICAE_NASCET_70 + ICAI_stenosis50 + M_prev_ht, data = data1, x=T, y=T) 
  x <- model.matrix(Mc_FailedFemoralApproach ~ M_age + ArchElongation_dichIII + AngleAaInn_OR_AaCca_dich45 + AorticVariant_dichB + AorticVariant_dichC + InnCca_90orLarger 
                    + ICA_90orLarger + AngleFollowingBifurcation + OrigoStenosis50percent + ICAE_NASCET_70 + ICAI_stenosis50 + M_prev_ht, data = data1)[,-1]
  set.seed(123)
  cv.glmmod <- cv.glmnet(mod01.1$x, mod01.1$y, alpha=1, family = "binomial") 
  par(mfrow=c(1,1))
  plot(cv.glmmod)
  model.L1 <- glmnet(mod01.1$x, mod01.1$y, alpha=1, family = "binomial", lambda = cv.glmmod$lambda.min) 
  lp1 <- x %*% coef(model.L1)[-1]
  Results.m[i,c(1,2)] <- c(i, cv.glmmod$lambda.min) 
  shrinkage <- 1/coef(lrm(mod01.1$y~ as.numeric(lp1)))
  Results.m[i,3] <-  shrinkage[2] 
  Results.m[i,4] <- sum(abs(coef(model.L1)) > 0) 
  Results.m[i,5] <- rcorr.cens(as.numeric(lp1),mod01.1$y)[2] / 2 + .5
}
Results.m
LASSO.10imp <- apply(Results.m, 2, mean) 
LASSO.10imp

# For supplementary figures lasso
mod01.1 <- lrm(formula = Mc_FailedFemoralApproach ~ M_age + ArchElongation_dichIII + AngleAaInn_OR_AaCca_dich45 + AorticVariant_dichB + AorticVariant_dichC + InnCca_90orLarger 
               + ICA_90orLarger + AngleFollowingBifurcation + OrigoStenosis50percent + ICAE_NASCET_70 + ICAI_stenosis50 + M_prev_ht, data = part1_imp1, x=T, y=T) 
cv.glmmod <- cv.glmnet(mod01.1$x, mod01.1$y, alpha=1, family = "binomial") 
model.L1 <- glmnet(mod01.1$x, mod01.1$y, alpha=1, family = "binomial", lambda = cv.glmmod$lambda.min) 
par(mfrow=c(1,1))
plot(cv.glmmod)

# Boostrap function for lasso predictor selection
boot.LASSO <- function(data=part1_imp10, x=mod01.1$x, y=mod01.1$y, B=10){
  j <- sample(nrow(part1_imp10), replace=T)
  Rmat <- matrix(nrow = B, ncol=5) 
  colnames(Rmat) <- c("lambda", "slope.orig", "D.orig", "slope.test", "D.test")
  yj <- part1_imp10[j,1] 
  for (i in 1:B){
    cv.glmmod <- cv.glmnet(x=mod01.1$x[j,], y=yj, alpha=1, family="binomial") 
    
    Rmat[1] <- c(LASSO.10imp[2]) 
    model.L1j <- glmnet(x=mod01.1$x[j,], y=yj, alpha=1, lambda = LASSO.10imp[2], family = "binomial")
    xj <- model.matrix(Mc_FailedFemoralApproach ~ M_age + ArchElongation_dichIII + AngleAaInn_OR_AaCca_dich45 + AorticVariant_dichB + AorticVariant_dichC + InnCca_90orLarger 
                       + ICA_90orLarger + AngleFollowingBifurcation + OrigoStenosis50percent + ICAE_NASCET_70 + ICAI_stenosis50 + M_prev_ht, data = part1_imp10[j, ])
    x <- model.matrix(Mc_FailedFemoralApproach ~ M_age + ArchElongation_dichIII + AngleAaInn_OR_AaCca_dich45 + AorticVariant_dichB + AorticVariant_dichC + InnCca_90orLarger 
                      + ICA_90orLarger + AngleFollowingBifurcation + OrigoStenosis50percent + ICAE_NASCET_70 + ICAI_stenosis50 + M_prev_ht, data = part1_imp10)
    lp1j <- xj %*% coef(model.L1j) 
    Rmat[2] <- coef(lrm(yj~as.numeric(lp1j)))[-1] 
    Rmat[3] <- rcorr.cens(as.numeric(lp1j), yj)[2]
    lp1 <-  x %*% coef(model.L1j) 
    shrinkage <- coef(lrm(mod01.1$y~as.numeric(lp1)))
    Rmat[4] <- shrinkage[2]
    Rmat[5] <- rcorr.cens(as.numeric(lp1), mod01.1$y)[2] 
  }
  return(Rmat)
}

# Predictor selection with previously defined bootstrap function 
B <- 200 
Results.mat <- matrix(nrow=1, ncol=(4+5)) 
colnames(Results.mat) <- c("D.full", "slope.full", "D.step","slope.step",
                           "lambda","slope.opt", "D.opt", "slope.corr", "D.corr"  )

mod01.1 <- lrm(formula = Mc_FailedFemoralApproach ~ M_age + ArchElongation_dichIII + AngleAaInn_OR_AaCca_dich45 + AorticVariant_dichB + AorticVariant_dichC + InnCca_90orLarger 
               + ICA_90orLarger + AngleFollowingBifurcation + OrigoStenosis50percent + ICAE_NASCET_70 + ICAI_stenosis50 + M_prev_ht, weights = w, data = part1_imp10, x=T, y=T) 
x <- model.matrix(Mc_FailedFemoralApproach ~ M_age + ArchElongation_dichIII + AngleAaInn_OR_AaCca_dich45 + AorticVariant_dichB + AorticVariant_dichC + InnCca_90orLarger 
                  + ICA_90orLarger + AngleFollowingBifurcation + OrigoStenosis50percent + ICAE_NASCET_70 + ICAI_stenosis50 + M_prev_ht, data = part1_imp10)

Results.mat[1:2] <- validate(mod01.1, B=B, maxit=99)[c(1,3), 5] 
Results.mat[3:4] <- validate(mod01.1, B=B, bw=T, rule = "aic", maxit=99)[c(1,3), 5] 

set.seed(123)
cv.glmmod <- cv.glmnet(mod01.1$x, mod01.1$y, weights = part1_imp10$w, alpha=1, family = "binomial") 
model.L1 <- glmnet(x=mod01.1$x, y=mod01.1$y, alpha=1, lambda = LASSO.10imp[2], family = "binomial") 
lp1 <- x %*% coef(model.L1) 
LASSO.means <- apply(boot.LASSO(data=part1_imp10, x=mod01.1$x, y=mod01.1$y, B=B), 2, mean)
Results.mat[5] <- LASSO.means[2]
Results.mat[6] <- LASSO.means[3] - LASSO.means[5] 
Results.mat[7] <- LASSO.means[4] - LASSO.means[6] 
slope <- 1/coef(lrm(mod01.1$y~ as.numeric(lp1)))
Results.mat[8] <- slope[2] - Results.mat[7] 
Results.mat[9] <- rcorr.cens(as.numeric(lp1), mod01.1$y)[2]- Results.mat[8] 
apply(Results.mat, 2, mean) 
apply(Results.mat, 2, mean)[c(2,4,10)]/2 +.5
LASSO.10imp <- apply(Results.mat, 2, mean)  
LASSO.10imp

# Validation of the model before the lasso
mod01DEF <- lrm(Mc_FailedFemoralApproach ~ ArchElongation_dichIII + AorticVariant_dichC + InnCca_90orLarger + AngleFollowingBifurcation +ICAE_NASCET_70, data = part1_imp10, weights = w, x=T, y=T)
set.seed(123)
validate(mod01DEF, B=200, bw=F, type="individual", dxy=T)
rcorr.cens(mod01DEF$linear.predictors, mod01DEF$y)
boot <- rcorr.cens(mod01DEF$linear.predictors, mod01DEF$y)
boot[1] -1.96*boot[3]/2
boot[1] +1.96*boot[3]/2

# Fit lasso model and insert coefficients in lrm
lassomodel <- glmnet(x=mod01DEF$x, y=mod01DEF$y, weights =part1_imp10$w, alpha=1, lambda = LASSO.10imp[2], family="binomial")
coef(lassomodel)
round(coef(lassomodel)/mod01DEF$coef,3) 
mean(coef(lassomodel)/mod01DEF$coef) 
lassocoef <- c(Intercept = -4.3125508, ArchElongation_dichIII = 1.6299069, AorticVariant_dichC = 0.9145084, InnCca_90orLarger = 1.0513471, 
               AngleFollowingBifurcation = 0.7295294, ICAE_NASCET_70 = 1.2001544)
mod01DEF$coefficients <- lassocoef


# Validation of this model after the lasso
validate(mod01DEF, B=200, bw=F, type="individual", dxy=T)
rcorr.cens(mod01DEF$linear.predictors, mod01DEF$y)
boot <- rcorr.cens(mod01DEF$linear.predictors, mod01DEF$y)
boot[1] -1.96*boot[3]/2
boot[1] +1.96*boot[3]/2

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