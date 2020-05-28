#' ---
#'title: "Data manipulations"
#'author: Julia van Hees
#' date: April 5th, 2020
#' output: 
#'      html_document:
#'             toc: true
#'             highlight: default
#'             theme: readable
#' ---
#------------- Introduction -------------------------
#' 
#' **Included in this file:**
#' 
#' * Load part 1 (derivation cohort) data into R and check data quality
#' * Load part 2 (validation cohort) data into R and check data quality
#' * Clinical variable distributions   
#' * Treatment variable distributions
#' * Imaging variable distributions
#' * Vascular Anatomical variable distributions
#' <br>
#' <br>
#------------------------ Load part 1 -----------------------------
#'
#' ## Load part 1 (derivation cohort) data into R and check data quality
#' <br>
library(foreign, quietly = TRUE)
library(rms, verbose = FALSE, quietly = TRUE)
library(mctest, verbose = FALSE, quietly = TRUE)
library(mice, quietly = TRUE)
library(glmnet,quietly = TRUE)
library(glmpath,quietly = TRUE)
library(Hmisc, quietly = TRUE)
#' <br>
part1 <- read.spss("~/Documents/OneDrive/Voorbereiding R coding/NewCode/part1_imp1.sav", to.data.frame = TRUE, use.value.labels = TRUE) 
part2 <- read.spss("~/Documents/OneDrive/Voorbereiding R coding/NewCode/part2_imp1.sav", to.data.frame = TRUE, use.value.labels = TRUE)
#' <br>
part1 <- mice(part1_imp1, m=10, seed=1, print=FALSE) 
part2 <- mice(part2_imp1, m=10, seed=1, print=FALSE) 
#' <br>
part1_imp1 <- complete(part1_imp1, 1)
part2_imp1 <- complete(part2_imp1, 1)
#' <br>
#' <br>
#' 
#------------------------ Clinical -----------------------------
#'
#' ## Clinical variable distributions 
#' <br>
#' The distribution of numeric (continuous) variables was determined first.
#' Start with M_age. Check if normally distributed. If the Shapiro-Wilk test p>0.05, then normality can be assumed.
shapiro.test(part1_imp1$M_age) # w=0.96, p=1.1e-14
#' Also check with a histogram:  
par(mfrow=c(1,1))
hist(part1_imp1$M_age,probability = T, main = "Distribution of M_age" )
#' And a QQplot: 
qqnorm(part1_imp1$M_age, main = "QQ plot for M_age in part 1", pch=19)
qqline(part1_imp1$M_age)
#' Shapiro-Wilk indicates no normality in the data. The data in the histogram seems to be a little skewed to the right. Additionally, the two tails of the QQ plot deviate from the middle line to the same side, this also indicates non-normality in the data. Due to this, I calculated a median with interquartile range (IQR) for distribution of this variable.
library(plyr)
median <- median(part1_imp1$M_age, na.rm = TRUE)
iqr <- IQR(part1_imp1$M_age, na.rm = TRUE)
iqr_low <- quantile(part1_imp1$M_age, na.rm = TRUE)[2]
iqr_high <- quantile(part1_imp1$M_age, na.rm = TRUE)[4]
print(median)
print(iqr_low)
print(iqr_high)
#' <br>
#' Now for M_NIHSS_BL (severity of stroke; value 0 to 4):
#' Shapiro-Wilk:
shapiro.test(part1_imp1$M_NIHSS_BL) # w=0.98, p=5.3e-07
#' histogram:
library(MASS)
par(mfrow=c(1,1))
truehist(part1_imp1$M_NIHSS_BL, main = "Distribution of M_NIHSS_BL", col = "white")
#' And a QQplot: 
qqnorm(part1_imp1$M_NIHSS_BL, main = "QQ plot for M_NIHSS_BL in part 1", pch=19)
qqline(part1_imp1$M_NIHSS_BL)
#' Calculate median and IQR:
median <- median(part1_imp1$M_NIHSS_BL, na.rm = TRUE)
iqr <- IQR(part1_imp1$M_NIHSS_BL, na.rm = TRUE)
iqr_low <- quantile(part1_imp1$M_NIHSS_BL, na.rm = TRUE)[2]
iqr_high <- quantile(part1_imp1$M_NIHSS_BL, na.rm = TRUE)[4]
print(median)
print(iqr_low)
print(iqr_high)
#' <br>
#' Continue on with M_togroin (duration from onset of the stroke to puncture of the groin artery):
shapiro.test(part1_imp1$M_togroin) # w=0.98, p=3.9e-07
#' Histogram:
par(mfrow=c(1,1))
truehist(part1_imp1$M_togroin, main = "Distribution of M_togroin", col = "white")
#' And a QQplot: 
qqnorm(part1_imp1$M_togroin, main = "QQ plot for M_togroin in part 1", pch=19)
qqline(part1_imp1$M_togroin)
#' Calculate the median and IQR values:
median <- median(part1_imp1$M_togroin, na.rm = TRUE)
iqr <- IQR(part1_imp1$M_togroin, na.rm = TRUE)
iqr_low <- quantile(part1_imp1$M_togroin, na.rm = TRUE)[2]
iqr_high <- quantile(part1_imp1$M_togroin, na.rm = TRUE)[4]
print(median)
print(iqr_low)
print(iqr_high)
#' <br>
#' Now M_rr_syst (systolic blood pressure):
shapiro.test(part1_imp1$M_rr_syst) # w=0.99, p=1.4e-05
#' Histogram:
par(mfrow=c(1,1))
truehist(part1_imp1$M_rr_syst, main = "Distribution of M_rr_syst", col = "white")
#' QQ plot:
qqnorm(part1_imp1$M_rr_syst, main = "QQ plot for M_rr_syst in part 1", pch=19)
qqline(part1_imp1$M_rr_syst)
#' Calculate the median and IQR values:
median <- median(part1_imp1$M_rr_syst, na.rm = TRUE)
iqr <- IQR(part1_imp1$M_rr_syst, na.rm = TRUE)
iqr_low <- quantile(part1_imp1$M_rr_syst, na.rm = TRUE)[2]
iqr_high <- quantile(part1_imp1$M_rr_syst, na.rm = TRUE)[4]
print(median)
print(iqr_low)
print(iqr_high)
#' <br>
#' Other clinical variables were factors and hence the total number with percentage was calculated:
#' Sex, female or male 
total.men<- summary(part1_imp1$M_sex)[2]
percentage <- ((summary(part1_imp1$M_sex)[2])/887)*100
print(total.men,0)
print(percentage,0)
#' Hypertension, Yes or No 
total.ht<- summary(part1_imp1$M_prev_ht)[2]
percentage <- ((summary(part1_imp1$M_prev_ht)[2])/887)*100
print(total.ht,0)
print(percentage,0)
#' Diabetes mellitus, Yes or No  
total.dm<- summary(part1_imp1$M_prev_dm)[2]
percentage <- ((summary(part1_imp1$M_prev_dm)[2])/887)*100
print(total.dm,0)
print(percentage,0)
#' Previous ischemic stroke, Yes or No  
total.str<- summary(part1_imp1$M_prev_str)[2]
percentage <- ((summary(part1_imp1$M_prev_str)[2])/887)*100
print(total.str,0)
print(percentage,0)
#' Atrial fibrillation, Yes or No
total.af<- summary(part1_imp1$M_prev_af)[2]
percentage <- ((summary(part1_imp1$M_prev_af)[2])/887)*100
print(total.af,0)
print(percentage,0)
#' Prestroke modified ranking scale, 0, 1, 2 >3
total.0 <- summary(part1_imp1$M_premrs)[1]
percentage <- ((summary(part1_imp1$M_premrs)[1])/887)*100
print(total.0,0) # mrs 0
print(percentage,0)
total.1 <- summary(part1_imp1$M_premrs)[2]
percentage <- ((summary(part1_imp1$M_premrs)[2])/887)*100
print(total.1,0) # mrs 1
print(percentage,0)
total.2 <- summary(part1_imp1$M_premrs)[3]
percentage <- ((summary(part1_imp1$M_premrs)[3])/887)*100
print(total.2,0) # mrs 2
print(percentage,0)
total.35 <- (summary(part1_imp1$M_premrs)[4])+(summary(part1_imp1$M_premrs)[5])+(summary(part1_imp1$M_premrs)[6])
percentage <- (total.35/887)*100
print(total.35,0) # mrs >3
print(percentage,0)
#' <br>
#------------------------ Derivation cohort: Treatment ----------------------------- 
#'
#' ## Treatment variable distributions
#' <br>
#' Now for M_durproc (duration of the procedure):
shapiro.test(part1_imp1$M_durproc) # w=0.95, p=3.2e-15
#' Histogram:
par(mfrow=c(1,1))
truehist(part1_imp1$M_durproc, main = "Distribution of M_durproc", col = "white")
#' QQ plot:
qqnorm(part1_imp1$M_durproc, main = "QQ plot for M_durproc in part 1", pch=19)
qqline(part1_imp1$M_durproc)
#' Calculate the median and IQR values:
median <- median(part1_imp1$M_durproc, na.rm = TRUE)
iqr <- IQR(part1_imp1$M_durproc, na.rm = TRUE)
iqr_low <- quantile(part1_imp1$M_durproc, na.rm = TRUE)[2]
iqr_high <- quantile(part1_imp1$M_durproc, na.rm = TRUE)[4]
print(median)
print(iqr_low)
print(iqr_high)
#' <br>
#' Other treatment variables were factors and hence the total number with percentage was calculated:
#' Failed femoral approach, Yes or No
part1_imp1$Mc_FailedFemoralApproach <- factor(part1_imp1$Mc_FailedFemoralApproach)
total.failed <- summary(part1_imp1$Mc_FailedFemoralApproach)[2]
percentage <- ((summary(part1_imp1$Mc_FailedFemoralApproach)[2])/887)*100
print(total.failed,0)
print(percentage,0)
#' Treatment with IV alteplase, Yes or No
total.1 <- summary(part1_imp1$M_ivtrom)[1]
percentage <- ((summary(part1_imp1$M_ivtrom)[1])/887)*100
print(total.1,0) # ivtrom No
print(percentage,0)
total.2 <- summary(part1_imp1$M_ivtrom)[2]
percentage <- ((summary(part1_imp1$M_ivtrom)[2])/887)*100
print(total.2,0) # ivtrom Yes
print(percentage,0)
#' Post-eTICI, 0, 1, 2(A-C), 3 
total.ptici0<- summary(part1_imp1$M_posttici_c)[1] # postici 0
percentage <- (total.ptici0/887)*100
print(total.ptici0,0)
print(percentage,0)
total.ptici1<- summary(part1_imp1$M_posttici_c)[2] # postici 1
percentage <- (total.ptici1/887)*100
print(total.ptici1,0)
print(percentage,0)
total.ptici2a<- summary(part1_imp1$M_posttici_c)[3] # postici 2A
percentage <- (total.ptici2a/887)*100
print(total.ptici2a,0)
print(percentage,0)
total.ptici2b<- summary(part1_imp1$M_posttici_c)[4] # postici 2B
percentage <- (total.ptici2b/887)*100
print(total.ptici2b,0)
print(percentage,0)
total.ptici2c<- summary(part1_imp1$M_posttici_c)[5] # postici 2C
percentage <- (total.ptici2c/887)*100
print(total.ptici2c,0)
print(percentage,0)
total.ptici3<- summary(part1_imp1$M_posttici_c)[6] # postici 3
percentage <- (total.ptici3/887)*100
print(total.ptici3,0)
print(percentage,0)
#------------------------ Derivation cohort: Imaging -----------------------------
#'
#' ## Imaging variable distributions
#' <br>
#' Continue on with M_ASPECTS_BL (baseline ASPECTS score; determines severity of the stroke by looking at the penetration in several brain areas):
shapiro.test(part1_imp1$M_ASPECTS_BL) # w=0.83, p=2.2e-16
#' Histogram:
par(mfrow=c(1,1))
truehist(part1_imp1$M_ASPECTS_BL, main = "Distribution of M_ASPECTS_BL", col = "white")
#' QQ plot:
qqnorm(part1_imp1$M_ASPECTS_BL, main = "QQ plot for M_ASPECTS_BL in part 1", pch=19)
qqline(part1_imp1$M_ASPECTS_BL)
#' Calculate the median and IQR values:
median <- median(part1_imp1$M_ASPECTS_BL, na.rm = TRUE)
iqr <- IQR(part1_imp1$M_ASPECTS_BL, na.rm = TRUE)
iqr_low <- quantile(part1_imp1$M_ASPECTS_BL, na.rm = TRUE)[2]
iqr_high <- quantile(part1_imp1$M_ASPECTS_BL, na.rm = TRUE)[4]
print(median)
print(iqr_low)
print(iqr_high)
#' <br>
#' Now the last numeric variable in part1_imp1, M_CBS_BL (clot burden score):
shapiro.test(part1_imp1$M_CBS_BL) # w=0.93, p=2.2e-16
#' Histogram:
par(mfrow=c(1,1))
truehist(part1_imp1$M_CBS_BL, main = "Distribution of M_CBS_BL", col = "white")
#' QQ plot:
qqnorm(part1_imp1$M_CBS_BL, main = "QQ plot for M_CBS_BL in part 1", pch=19)
qqline(part1_imp1$M_CBS_BL)
#' Median and IQR:
median <- median(part1_imp1$M_CBS_BL, na.rm = TRUE)
iqr <- IQR(part1_imp1$M_CBS_BL, na.rm = TRUE)
iqr_low <- quantile(part1_imp1$M_CBS_BL, na.rm = TRUE)[2]
iqr_high <- quantile(part1_imp1$M_CBS_BL, na.rm = TRUE)[4]
print(median)
print(iqr_low)
print(iqr_high)
#' <br>
#' Other imaging variables were factors and hence the total number with percentage was calculated:
#' Pre-eTICI, 0, 1, 2A, 2B-3
total.1 <- summary(part1_imp1$Mc_pretici_short)[1]
percentage <- ((summary(part1_imp1$Mc_pretici_short)[1])/887)*100
print(total.1,0) # pretici 0
print(percentage,0)
total.2 <- summary(part1_imp1$Mc_pretici_short)[2]
percentage <- ((summary(part1_imp1$Mc_pretici_short)[2])/887)*100
print(total.2,0) # pretici 1
print(percentage,0)
total.3 <- summary(part1_imp1$Mc_pretici_short)[3]
percentage <- ((summary(part1_imp1$Mc_pretici_short)[3])/887)*100
print(total.3,0) # pretici 2
print(percentage,0)
total.4 <- (summary(part1_imp1$Mc_pretici_short)[4])
percentage <- (total.4/887)*100
print(total.4,0) # pretici >3
print(percentage,0)
print(percentage,0)
#' Occlusion segment, ICA-I, ICA-T, M1, M2-M3
total.1 <- summary(part1_imp1$M_occlsegment_c_short)[1]
percentage <- ((summary(part1_imp1$M_occlsegment_c_short)[1])/887)*100
print(total.1,0) # ICAI
print(percentage,0)
total.2 <- summary(part1_imp1$M_occlsegment_c_short)[2]
percentage <- ((summary(part1_imp1$M_occlsegment_c_short)[2])/887)*100
print(total.2,0) # ICAT
print(percentage,0)
total.3 <- summary(part1_imp1$M_occlsegment_c_short)[3]
percentage <- ((summary(part1_imp1$M_occlsegment_c_short)[3])/887)*100
print(total.3,0) # M1
print(percentage,0)
total.4 <- summary(part1_imp1$M_occlsegment_c_short)[4]
percentage <- ((summary(part1_imp1$M_occlsegment_c_short)[4])/887)*100
print(total.4,0) # Other
print(percentage,0)
#' Collateral score, absent, <50% filling, >50<100 filling, 100% of occluded area
total.1 <- summary(part1_imp1$M_collaterals)[1]
percentage <- ((summary(part1_imp1$M_collaterals)[1])/887)*100
print(total.1,0) # absent
print(percentage,0)
total.2 <- summary(part1_imp1$M_collaterals)[2]
percentage <- ((summary(part1_imp1$M_collaterals)[2])/887)*100
print(total.2,0) # <50% fillin
print(percentage,0)
total.3 <- summary(part1_imp1$M_collaterals)[3]
percentage <- ((summary(part1_imp1$M_collaterals)[3])/887)*100
print(total.3,0) # >50<100 fillin
print(percentage,0)
total.4 <- (summary(part1_imp1$M_collaterals)[4])
percentage <- (total.4/887)*100
print(total.4,0) # 100% of occluded area
print(percentage,0)
#' Prestroke modified ranking scale reversed, 0, 1, 2, >3 
total.0 <- summary(part1_imp1$M_premrs)[1]
percentage <- ((summary(part1_imp1$M_premrs)[1])/887)*100
print(total.0,0) #  0
print(percentage,0)
total.1 <- summary(part1_imp1$M_premrs)[2]
percentage <- ((summary(part1_imp1$M_premrs)[2])/887)*100
print(total.1,0) #  1
print(percentage,0)
total.2 <- summary(part1_imp1$M_premrs)[3]
percentage <- ((summary(part1_imp1$M_premrs)[3])/887)*100
print(total.2,0) #  2
print(percentage,0)
total.35 <- (summary(part1_imp1$M_premrs)[4])+(summary(part1_imp1$M_premrs)[5])+(summary(part1_imp1$M_premrs)[6])
percentage <- (total.35/887)*100
print(total.35,0) #  >3
print(percentage,0)
#------------------------ Derivation cohort: Vascular Anatomy -----------------------------
#'
#' ## Vascular Anatomical variable distributions
#' <br>
#' Now that I have processed the numeric clinical, treatment and imaging variables in part1_imp1, I still have to do the variables that we have measured ourselves, the vascular anatomy variables.
#' I am going to start with AngleAaInn_Or_AaCca
#' Shapiro-Wilk:
shapiro.test(part1_imp1$AngleAaInn_OR_AaCca) # w=0.93, p=2.2e-16
#' Histogram:
par(mfrow=c(1,1))
truehist(part1_imp1$AngleAaInn_OR_AaCca, main = "Distribution of AngleAaInn_OR_AaCca", col = "white")
#' QQ plot:
qqnorm(part1_imp1$AngleAaInn_OR_AaCca, main = "QQ plot for AngleAaInn_OR_AaCca in part 1", pch=19)
qqline(part1_imp1$AngleAaInn_OR_AaCca)
#' Calculate mean and IQR:
median <- median(part1_imp1$AngleAaInn_OR_AaCca, na.rm = TRUE)
iqr <- IQR(part1_imp1$AngleAaInn_OR_AaCca, na.rm = TRUE)
iqr_low <- quantile(part1_imp1$AngleAaInn_OR_AaCca, na.rm = TRUE)[2]
iqr_high <- quantile(part1_imp1$AngleAaInn_OR_AaCca, na.rm = TRUE)[4]
print(median)
print(iqr_low)
print(iqr_high) 
#' <br> 
#' Now AngleCcaIca
#' Shapiro-Wilk:
shapiro.test(part1_imp1$AngleCcaIca) # w=0.97, p=1.1e-11
#' Histogram:
par(mfrow=c(1,1))
truehist(part1_imp1$AngleCcaIca, main = "Distribution of AngleCcaIca", col = "white")
#' QQ plot:
qqnorm(part1_imp1$AngleCcaIca, main = "QQ plot for AngleCcaIca in part 1", pch=19)
qqline(part1_imp1$AngleCcaIca)
#' Determine the median and IQR:
median <- median(part1_imp1$AngleCcaIca, na.rm = TRUE)
iqr <- IQR(part1_imp1$AngleCcaIca, na.rm = TRUE)
iqr_low <- quantile(part1_imp1$AngleCcaIca, na.rm = TRUE)[2]
iqr_high <- quantile(part1_imp1$AngleCcaIca, na.rm = TRUE)[4]
print(median)
print(iqr_low)
print(iqr_high)
#' <br>
#' Continue on with the NASCET score (how much smaller is the lumen of the ICA?)
#' shapiro-Wilk:
shapiro.test(part1_imp1$ICAE_NASCET_Degree) # w=0.61, p=2.2e-16
#' Histogram:
par(mfrow=c(1,1))
truehist(part1_imp1$ICAE_NASCET_Degree, main = "Distribution of ICAE_NASCET_Degree" )
#' QQ plot:
qqnorm(part1_imp1$ICAE_NASCET_Degree, main = "QQ plot for ICAE_NASCET_Degree in part 1", pch=19)
qqline(part1_imp1$ICAE_NASCET_Degree)
#' <br>
#' Other vascular anatomical variables were factors and hence the total number with percentage was calculated:
#' Aortic variant, A, B, C
total.1 <- summary(part1_imp1$AorticVariant)[1]+summary(part1_imp1$AorticVariant)[4]
percentage <- ((summary(part1_imp1$AorticVariant)[1])/887)*100
print(total.1,0) # A+D
print(percentage,0)
total.2 <- summary(part1_imp1$AorticVariant)[2]
percentage <- ((summary(part1_imp1$AorticVariant)[2])/887)*100
print(total.2,0) # B
print(percentage,0)
total.3 <- summary(part1_imp1$AorticVariant)[3]
percentage <- ((summary(part1_imp1$AorticVariant)[3])/887)*100
print(total.3,0) # C
print(percentage,0)
#' Arch elongation, I, II, III
total.1 <- summary(part1_imp1$ArchElongation)[1]
percentage <- ((summary(part1_imp1$ArchElongation)[1])/887)*100
print(total.1,0) # A+D
print(percentage,0)
total.2 <- summary(part1_imp1$ArchElongation)[2]
percentage <- ((summary(part1_imp1$ArchElongation)[2])/887)*100
print(total.2,0) # B
print(percentage,0)
total.3 <- summary(part1_imp1$ArchElongation)[3]
percentage <- ((summary(part1_imp1$ArchElongation)[3])/887)*100
print(total.3,0) # C
print(percentage,0)
#' Inn origin atherosclerosis, not present, yes, both, not applicable
total.1 <- summary(part1_imp1$InnAtherosclerosis)[1]
percentage <- ((summary(part1_imp1$InnAtherosclerosis)[1])/887)*100
print(total.1,0) # not present
print(percentage,0)
total.2 <- summary(part1_imp1$InnAtherosclerosis)[2]
percentage <- ((summary(part1_imp1$InnAtherosclerosis)[2])/887)*100
print(total.2,0) # yes
print(percentage,0)
total.3 <- summary(part1_imp1$InnAtherosclerosis)[3]
percentage <- ((summary(part1_imp1$InnAtherosclerosis)[3])/887)*100
print(total.3,0) # both
print(percentage,0)
total.4 <- summary(part1_imp1$InnAtherosclerosis)[4]
percentage <- ((summary(part1_imp1$InnAtherosclerosis)[4])/887)*100
print(total.4,0) # not applicable
print(percentage,0)
#' CCA origin atherosclerosis, not present, yes, both, not applicable
total.1 <- summary(part1_imp1$CCAOrigin)[1]
percentage <- ((summary(part1_imp1$CCAOrigin)[1])/887)*100
print(total.1,0) # not present
print(percentage,0)
total.2 <- summary(part1_imp1$CCAOrigin)[2]
percentage <- ((summary(part1_imp1$CCAOrigin)[2])/887)*100
print(total.2,0) # yes
print(percentage,0)
total.3 <- summary(part1_imp1$CCAOrigin)[3]
percentage <- ((summary(part1_imp1$CCAOrigin)[3])/887)*100
print(total.3,0) # both
print(percentage,0)
total.4 <- summary(part1_imp1$CCAOrigin)[4]
percentage <- ((summary(part1_imp1$CCAOrigin)[4])/887)*100
print(total.4,0) # not applicable 
print(percentage,0)
#' Inn origin atherosclerosis, not present, yes, both, not applicable 
total.1 <- summary(part1_imp1$InnOrigin)[1]
percentage <- ((summary(part1_imp1$InnOrigin)[1])/887)*100
print(total.1,0) # not present
print(percentage,0)
total.2 <- summary(part1_imp1$InnOrigin)[2]
percentage <- ((summary(part1_imp1$InnOrigin)[2])/887)*100
print(total.2,0) # yes
print(percentage,0)
total.3 <- summary(part1_imp1$InnOrigin)[3]
percentage <- ((summary(part1_imp1$InnOrigin)[3])/887)*100
print(total.3,0) # both
print(percentage,0)
total.4 <- summary(part1_imp1$InnOrigin)[4]
percentage <- ((summary(part1_imp1$InnOrigin)[4])/887)*100
print(total.4,0) # not applicable
print(percentage,0)
#' CCA atherosclerosis, not present, yes spots, yes<50, yes>50
summary(part1_imp1$CcaAherosclerosis)
total.1 <- summary(part1_imp1$CcaAherosclerosis)[1]
percentage <- ((summary(part1_imp1$CcaAherosclerosis)[1])/887)*100
print(total.1,0) # not present
print(percentage,0)
total.2 <- summary(part1_imp1$CcaAherosclerosis)[2]
percentage <- ((summary(part1_imp1$CcaAherosclerosis)[2])/887)*100
print(total.2,0) # yes spots
print(percentage,0)
total.3 <- summary(part1_imp1$CcaAherosclerosis)[3]
percentage <- ((summary(part1_imp1$CcaAherosclerosis)[3])/887)*100
print(total.3,0) # yes <50
print(percentage,0)
total.4 <- summary(part1_imp1$CcaAherosclerosis)[4]
percentage <- ((summary(part1_imp1$CcaAherosclerosis)[4])/887)*100
print(total.4,0) # yes >50 
print(percentage,0)
#' Number of 90 angles in Inn + CCA, 0, 1, 2, >3
total.1 <- summary(part1_imp1$InnCca_nr_90orLarger)[1]
percentage <- ((summary(part1_imp1$InnCca_nr_90orLarger)[1])/887)*100
print(total.1,0) # 0 angles
print(percentage,0)
total.2 <- summary(part1_imp1$InnCca_nr_90orLarger)[2]
percentage <- ((summary(part1_imp1$InnCca_nr_90orLarger)[2])/887)*100
print(total.2,0) # 1 angle
print(percentage,0)
total.3 <- summary(part1_imp1$InnCca_nr_90orLarger)[3]
percentage <- ((summary(part1_imp1$InnCca_nr_90orLarger)[3])/887)*100
print(total.3,0) # 2 angles
print(percentage,0)
total.4 <- summary(part1_imp1$InnCca_nr_90orLarger)[4]
percentage <- ((summary(part1_imp1$InnCca_nr_90orLarger)[4])/887)*100
print(total.4,0) # 3 angles
print(percentage,0)
#' Number of 90 angles in ICA, 0, 1, 2, >3
total.1 <- summary(part1_imp1$ICA_nr_90orLarger)[1]
percentage <- ((summary(part1_imp1$ICA_nr_90orLarger)[1])/887)*100
print(total.1,0) # 0 angles
print(percentage,0)
total.2 <- summary(part1_imp1$ICA_nr_90orLarger)[2]
percentage <- ((summary(part1_imp1$ICA_nr_90orLarger)[2])/887)*100
print(total.2,0) # 1 angle
print(percentage,0)
total.2 <- summary(part1_imp1$ICA_nr_90orLarger)[3]
percentage <- ((summary(part1_imp1$ICA_nr_90orLarger)[3])/887)*100
print(total.2,0) # 1 angle
print(percentage,0)
total.3 <- summary(part1_imp1$ICA_nr_90orLarger)[4]+summary(part1_imp1$ICA_nr_90orLarger)[5]+summary(part1_imp1$ICA_nr_90orLarger)[6]
percentage <- (total.3/887)*100
print(total.3,0) # 2 angles
print(percentage,0)
#' ICAI atherosclerosis, not present, yes spots, yes<50, yes>50
total.1 <- summary(part1_imp1$ICAIAtherosclerosis)[1]
percentage <- ((summary(part1_imp1$ICAIAtherosclerosis)[1])/887)*100
print(total.1,0) # not present
print(percentage,0)
total.2 <- summary(part1_imp1$ICAIAtherosclerosis)[2]
percentage <- ((summary(part1_imp1$ICAIAtherosclerosis)[2])/887)*100
print(total.2,0) # yes spots
print(percentage,0)
total.3 <- summary(part1_imp1$ICAIAtherosclerosis)[3]
percentage <- ((summary(part1_imp1$ICAIAtherosclerosis)[3])/887)*100
print(total.3,0) # yes<50
print(percentage,0)
total.4 <- summary(part1_imp1$ICAIAtherosclerosis)[4]
percentage <- ((summary(part1_imp1$ICAIAtherosclerosis)[4])/887)*100
print(total.4,0) # yes>50
print(percentage,0)
#' ICAI stenosis, <50, >50
total.1 <- summary(part1_imp1$ICAI_stenosis50)[1]
percentage <- ((summary(part1_imp1$ICAI_stenosis50)[1])/887)*100
print(total.1,0) # <50
print(percentage,0)
total.2 <- summary(part1_imp1$ICAI_stenosis50)[2]
percentage <- ((summary(part1_imp1$ICAI_stenosis50)[2])/887)*100
print(total.2,0) # >50
print(percentage,0)
#' ICAE atherosclerosis, not present, yes calcified spots, yes stenosis
total.1 <- summary(part1_imp1$ICAEAtherosclerosis)[1]
percentage <- ((summary(part1_imp1$ICAEAtherosclerosis)[1])/887)*100
print(total.1,0) # not present
print(percentage,0)
total.2 <- summary(part1_imp1$ICAEAtherosclerosis)[2]
percentage <- ((summary(part1_imp1$ICAEAtherosclerosis)[2])/887)*100
print(total.2,0) # yes calcified spots
print(percentage,0)
total.3 <- summary(part1_imp1$ICAEAtherosclerosis)[3]
percentage <- ((summary(part1_imp1$ICAEAtherosclerosis)[3])/887)*100
print(total.3,0) # yes stenosis
print(percentage,0)
#' Origostenosis, not at origin, at origin
total.1 <- summary(part1_imp1$OrigoStenosis50percent)[1]
percentage <- ((summary(part1_imp1$OrigoStenosis50percent)[1])/887)*100
print(total.1,0) # 0 angles
print(percentage,0)
total.2 <- summary(part1_imp1$OrigoStenosis50percent)[2]
percentage <- ((summary(part1_imp1$OrigoStenosis50percent)[2])/887)*100
print(total.2,0) # 0 angles
print(percentage,0)
#' Aortic variant A+B, C
total.AB<- (summary(part1_imp1$AorticVariant)[1])+(summary(part1_imp1$AorticVariant)[2])+(summary(part1_imp1$AorticVariant)[4]) # variant A, B and D
percentage <- (total.AB/887)*100
print(total.AB,0)
print(percentage,0)
total.C<- summary(part1_imp1$AorticVariant)[3] # variant C
percentage <- (total.C/887)*100
print(total.C,0)
print(percentage,0)
#' Arch elongation type I+II, III
total.I<- (summary(part1_imp1$ArchElongation)[1])+(summary(part1_imp1$ArchElongation)[2]) # Type I and II
percentage <- (total.I/887)*100
print(total.I,0)
print(percentage,0)
total.III<- summary(part1_imp1$ArchElongation)[3] # variant C
percentage <- (total.III/887)*100
print(total.III,0)
print(percentage,0)
#' Angle <45 AA and Inn or AA and CCA, Yes , No
total.1 <- summary(part1_imp1$AngleAaInn_OR_AaCca_dich45)[2]
percentage <- ((summary(part1_imp1$AngleAaInn_OR_AaCca_dich45)[2])/887)*100
print(total.1,0) # 0 angles
print(percentage,0)
#' 90 angle in the Inn or CCA, Yes, No 
total.ang<- (summary(part1_imp1$InnCca_90orLarger)[1])
percentage <- (total.ang/887)*100
print(total.ang,0)
print(percentage,0)
#' 90 angle after carotid bifurcation, Yes, No
total.bif<- (summary(part1_imp1$AngleFollowingBifurcation)[2])
percentage <- (total.bif/887)*100
print(total.bif,0)
print(percentage,0)
#' 90 angle in the ICA, Yes, No 
total.1 <- summary(part1_imp1$ICA_90orLarger)[2]
percentage <- ((summary(part1_imp1$ICA_90orLarger)[2])/887)*100
print(total.1,0) # 0 angles
print(percentage,0)
#' NASCET >70, Yes, No 
total.nascet<- (summary(part1_imp1$ICAE_NASCET_70)[2])
percentage <- (total.nascet/887)*100
print(total.nascet,0)
print(percentage,0)
#' <br>
#' These were all the variables for part 1. Same code was repeated for part 2 (code not included here). 
#' <br>
#' End of variable checks and distribution assessment. The above code is used to produce:
#' 
#' * Table 1
#' * Table 2
#'
#' <br>
#' <br> 
#' 
#' 
#' 
#' 
#' 