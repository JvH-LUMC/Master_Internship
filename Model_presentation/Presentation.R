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

# Nomogram
dd <- datadist(part2_imp10);options(datadist = "dd")
par(mar=c(2,2,1,1))
nom <- nomogram(mod02DEF, fun=function(x)1/(1+exp(-x)),lp.at = seq(-5,3, by=.5), maxscale=10, funlabel = "Failure", lp=F, fun.at = c(.001,.01,.05,seq(.1,.9,by=.1),.95,.99,.999))
plot(nom)
print(nom)


# Code function E. Steyerberg for validation plot
val.prob.ci.2 <- function(p, y, logit, group, weights = rep(1, length(y)), normwt = F, pl = T,
                          smooth = c("loess","rcs",F), CL.smooth="fill",CL.BT=F,lty.smooth=1,col.smooth="black",lwd.smooth=1,
                          nr.knots=5,logistic.cal = F,lty.log=1,col.log="black",lwd.log=1, xlab = "Predicted probability", ylab =
                            "Observed proportion", xlim = c(-0.02, 1),ylim = c(-0.15,1), m, g, cuts, emax.lim = c(0, 1),
                          legendloc =  c(0.50 , 0.27), statloc = c(0,.85),dostats=T,cl.level=0.95,method.ci="pepe",roundstats=2,
                          riskdist = "predicted", cex=0.75,cex.leg = 0.75, connect.group =
                            F, connect.smooth = T, g.group = 4, evaluate = 100, nmin = 0, d0lab="0", d1lab="1", cex.d01=0.7,
                          dist.label=0.04, line.bins=-.05, dist.label2=.03, cutoff, las=1, length.seg=1,
                          y.intersp=1,lty.ideal=1,col.ideal="red",lwd.ideal=1, pch.group=2, col.pch="black", cex.group=1, ...)
{
  if(smooth[1]==F){smooth <- "F"}
  smooth <- match.arg(smooth)
  if(!missing(p))
    if(any(!(p>=0 | p<=1))){stop("Probabilities can not be > 1 or < 0.")}
  if(missing(p))
    p <- 1/(1 + exp( - logit))
  else logit <- log(p/(1 - p))
  if(!all(y%in%0:1)){stop("The vector with the binary outcome can only contain the values 0 and 1.")}
  if(length(p) != length(y))
    stop("lengths of p or logit and y do not agree")
  names(p) <- names(y) <- names(logit) <- NULL
  if(!missing(group)) {
    if(length(group) == 1 && is.logical(group) && group)
      group <- rep("", length(y))
    if(!is.factor(group))
      group <- if(is.logical(group) || is.character(group))
        as.factor(group) else cut2(group, g =
                                     g.group)
    names(group) <- NULL
    nma <- !(is.na(p + y + weights) | is.na(group))
    ng <- length(levels(group))
  }
  else {
    nma <- !is.na(p + y + weights)
    ng <- 0
  }
  logit <- logit[nma]
  y <- y[nma]
  p <- p[nma]
  if(ng > 0) {
    group <- group[nma]
    weights <- weights[nma]
    return(val.probg(p, y, group, evaluate, weights, normwt, nmin)
    )
  }
  
  # Sort vector with probabilities
  y     <- y[order(p)]
  logit <- logit[order(p)]
  p     <- p[order(p)]
  
  
  if(length(p)>5000 & smooth=="loess"){warning("Number of observations > 5000, RCS is recommended.",immediate. = T)}
  if(length(p)>1000 & CL.BT==T){warning("Number of observations is > 1000, this could take a while...",immediate. = T)}
  
  
  if(length(unique(p)) == 1) {
    #22Sep94
    P <- mean(y)
    Intc <- log(P/(1 - P))
    n <- length(y)
    D <- -1/n
    L01 <- -2 * sum(y * logit - log(1 + exp(logit)), na.rm = T)
    L.cal <- -2 * sum(y * Intc - log(1 + exp(Intc)), na.rm = T)
    U.chisq <- L01 - L.cal
    U.p <- 1 - pchisq(U.chisq, 1)
    U <- (U.chisq - 1)/n
    Q <- D - U
    
    stats <- c(0, 0.5, 0, D, 0, 1, U, U.chisq, U.p, Q, mean((y - p[
      1])^2), Intc, 0, rep(abs(p[1] - P), 2))
    names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq",
                      "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier",
                      "Intercept", "Slope", "Emax", "Eavg", "ECI")
    return(stats)
  }
  i <- !is.infinite(logit)
  nm <- sum(!i)
  if(nm > 0)
    warning(paste(nm, "observations deleted from logistic calibration due to probs. of 0 or 1"))
  i.2 <- i
  f.or <- lrm(y[i]~logit[i])
  f <- lrm.fit(logit[i], y[i])
  cl.slope <- confint(f,level=cl.level)[2,]
  f2 <-	lrm.fit(offset=logit[i], y=y[i])
  if(f2$fail){
    warning("The lrm function did not converge when computing the calibration intercept!",immediate.=T)
    f2 <- list()
    f2$coef <- NA
    cl.interc <- rep(NA,2)
  }else{
    cl.interc <- confint(f2,level=cl.level)
  }
  stats <- f$stats
  cl.auc <- ci.auc(y,p,cl.level,method.ci)
  
  n <- stats["Obs"]
  predprob <- seq(emax.lim[1], emax.lim[2], by = 0.0005)
  lt <- f$coef[1] + f$coef[2] * log(predprob/(1 - predprob))
  calp <- 1/(1 + exp( - lt))
  emax <- max(abs(predprob - calp))
  if (pl) {
    plot(0.5, 0.5, xlim = xlim, ylim = ylim, type = "n", xlab = xlab,
         ylab = ylab, las=las,...)
    clip(0,1,0,1)
    abline(0, 1, lty = lty.ideal,col=col.ideal,lwd=lwd.ideal)
    do.call("clip", as.list(par()$usr))
    
    
    lt <- lty.ideal
    lw.d <- lwd.ideal
    all.col <- col.ideal
    leg <- "Ideal"
    marks <- -1
    if (logistic.cal) {
      lt <- c(lt, lty.log)
      lw.d <- c(lw.d,lwd.log)
      all.col <- c(all.col,col.log)
      leg <- c(leg, "Logistic calibration")
      marks <- c(marks, -1)
    }
    if(smooth!="F"){all.col <- c(all.col,col.smooth)}
    if (smooth=="loess") {
      #Sm <- lowess(p,y,iter=0)
      Sm <- loess(y~p,degree=2)
      Sm <- data.frame(Sm$x,Sm$fitted); Sm.01 <- Sm
      
      if (connect.smooth==T & CL.smooth!="fill") {
        clip(0,1,0,1)
        lines(Sm, lty = lty.smooth,lwd=lwd.smooth,col=col.smooth)
        do.call("clip", as.list(par()$usr))
        lt <- c(lt, lty.smooth)
        lw.d <- c(lw.d,lwd.smooth)
        marks <- c(marks, -1)
      }else if(connect.smooth==F & CL.smooth!="fill"){
        clip(0,1,0,1)
        points(Sm,col=col.smooth)
        do.call("clip", as.list(par()$usr))
        lt <- c(lt, 0)
        lw.d <- c(lw.d,1)
        marks <- c(marks, 1)
      }
      if(CL.smooth==T | CL.smooth=="fill"){
        to.pred <- seq(min(p),max(p),length=200)
        if(CL.BT==T){
          cat("Bootstrap samples are being generated.\n\n\n")
          
          replicate(2000,BT.samples(y,p,to.pred)) -> res.BT
          apply(res.BT,1,quantile,c(0.025,0.975)) -> CL.BT
          colnames(CL.BT) <- to.pred
          
          if(CL.smooth=="fill"){
            clip(0,1,0,1)
            polygon(x = c(to.pred, rev(to.pred)), y = c(CL.BT[2,],
                                                        rev(CL.BT[1,])),
                    col = rgb(177, 177, 177, 177, maxColorValue = 255), border = NA)
            if (connect.smooth==T) {
              lines(Sm, lty = lty.smooth,lwd=lwd.smooth,col=col.smooth)
              lt <- c(lt, lty.smooth)
              lw.d <- c(lw.d,lwd.smooth)
              marks <- c(marks, -1)
            }else if(connect.smooth==F){
              points(Sm,col=col.smooth)
              lt <- c(lt, 0)
              lw.d <- c(lw.d,1)
              marks <- c(marks, 1)
            }
            do.call("clip", as.list(par()$usr))
            leg <- c(leg, "Flexible calibration (Loess)")
          }else{
            
            clip(0,1,0,1)
            lines(to.pred,CL.BT[1,],lty=2,lwd=1,col=col.smooth);clip(0,1,0,1);lines(to.pred,CL.BT[2,],lty=2,lwd=1,col=col.smooth)
            do.call("clip", as.list(par()$usr))
            leg <- c(leg,"Flexible calibration (Loess)","CL flexible")
            lt <- c(lt,2)
            lw.d <- c(lw.d,1)
            all.col <- c(all.col,col.smooth)
            marks <- c(marks,-1)
          }
          
        }else{
          Sm.0 <- loess(y~p,degree=2)
          predict(Sm.0,type="fitted",se=T) -> cl.loess
          clip(0,1,0,1)
          if(CL.smooth=="fill"){
            polygon(x = c(Sm.0$x, rev(Sm.0$x)), y = c(cl.loess$fit+cl.loess$se.fit*1.96,
                                                      rev(cl.loess$fit-cl.loess$se.fit*1.96)),
                    col = rgb(177, 177, 177, 177, maxColorValue = 255), border = NA)
            if (connect.smooth==T) {
              lines(Sm, lty = lty.smooth,lwd=lwd.smooth,col=col.smooth)
              lt <- c(lt, lty.smooth)
              lw.d <- c(lw.d,lwd.smooth)
              marks <- c(marks, -1)
            }else if(connect.smooth==F){
              points(Sm,col=col.smooth)
              lt <- c(lt, 0)
              lw.d <- c(lw.d,1)
              marks <- c(marks, 1)
            }
            do.call("clip", as.list(par()$usr))
            leg <- c(leg, "Flexible calibration (Loess)")
          }else{
            lines(Sm.0$x,cl.loess$fit+cl.loess$se.fit*1.96,lty=2,lwd=1,col=col.smooth)
            lines(Sm.0$x,cl.loess$fit-cl.loess$se.fit*1.96,lty=2,lwd=1,col=col.smooth)
            do.call("clip", as.list(par()$usr))
            leg <- c(leg,"Flexible calibration (Loess)","CL flexible")
            lt <- c(lt,2)
            lw.d <- c(lw.d,1)
            all.col <- c(all.col,col.smooth)
            marks <- c(marks,-1)
          }
          
        }
        
      }else{
        leg <- c(leg, "Flexible calibration (Loess)")}
      cal.smooth <- approx(Sm.01, xout = p)$y
      eavg <- mean(abs(p - cal.smooth))
      ECI <- mean((p-cal.smooth)^2)*100
    }
    if(smooth=="rcs"){
      par(lwd=lwd.smooth,bty="n",col=col.smooth)
      if(!is.numeric(nr.knots)){stop("Nr.knots must be numeric.")}
      if(nr.knots==5){
        tryCatch(rcspline.plot(p,y,model="logistic",nk=5,show="prob", statloc = "none", noprint=T,
                               add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth),error=function(e){
                                 warning("The number of knots led to estimation problems, nk will be set to 4.",immediate. = T)
                                 tryCatch(rcspline.plot(p,y,model="logistic",nk=4,show="prob", statloc = "none", noprint=T,
                                                        add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth)
                                          ,error=function(e){
                                            warning("Nk 4 also led to estimation problems, nk will be set to 3.",immediate.=T)
                                            rcspline.plot(p,y,model="logistic",nk=3,show="prob", statloc = "none", noprint=T,
                                                          add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p)))
                                                          ,lty=lty.smooth)
                                          })
                               })
      }else if(nr.knots==4){
        tryCatch(rcspline.plot(p,y,model="logistic",nk=4,show="prob", statloc = "none", noprint=T,
                               add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth),error=function(e){
                                 warning("The number of knots led to estimation problems, nk will be set to 3.",immediate.=T)
                                 rcspline.plot(p,y,model="logistic",nk=3,show="prob", statloc = "none", noprint=T,
                                               add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth)
                               })
      }else if(nr.knots==3){
        tryCatch(rcspline.plot(p,y,model="logistic",nk=3,show="prob", statloc = "none", noprint=T,
                               add=T,showknots=F,xrange=c(min(na.omit(p)),max(na.omit(p))),lty=lty.smooth),
                 error=function(e){
                   stop("Nk = 3 led to estimation problems.")
                 })
      }else{stop(paste("Number of knots = ",nr.knots,sep="", ", only 5 >= nk >=3 is allowed."))}
      
      par(lwd=1,bty="o",col="black")
      leg <- c(leg,"Flexible calibration (RCS)","CL flexible")
      lt <- c(lt,lty.smooth,2)
      lw.d <- c(lw.d,rep(lwd.smooth,2))
      all.col <- c(all.col,col.smooth)
      marks <- c(marks,-1,-1)
    }
    if(!missing(m) | !missing(g) | !missing(cuts)) {
      if(!missing(m))
        q <- cut2(p, m = m, levels.mean = T, digits = 7)
      else if(!missing(g))
        q <- cut2(p, g = g, levels.mean = T, digits = 7)
      else if(!missing(cuts))
        q <- cut2(p, cuts = cuts, levels.mean = T, digits = 7)
      means <- as.single(levels(q))
      prop <- tapply(y, q, function(x)mean(x, na.rm = T))
      points(means, prop, pch = pch.group, cex=cex.group, col=col.smooth)
      #18.11.02: CI triangles
      ng	<-tapply(y, q, length)
      og	<-tapply(y, q, sum)
      ob	<-og/ng
      se.ob	<-sqrt(ob*(1-ob)/ng)
      g		<- length(as.single(levels(q)))
      
      for (i in 1:g) lines(c(means[i], means[i]), c(prop[i],min(1,prop[i]+1.96*se.ob[i])), type="l", col=col.pch)
      for (i in 1:g) lines(c(means[i], means[i]), c(prop[i],max(0,prop[i]-1.96*se.ob[i])), type="l", col=col.pch)
      
      if(connect.group) {
        lines(means, prop)
        lt <- c(lt, 1)
        lw.d <- c(lw.d,1)
      }
      else lt <- c(lt, 0)
      lw.d <- c(lw.d,1)
      leg <- c(leg, "Grouped observations")
      all.col <- c(all.col,col.smooth)
      marks <- c(marks, pch.group)
    }
  }
  lr <- stats["Model L.R."]
  p.lr <- stats["P"]
  D <- (lr - 1)/n
  L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm = TRUE)
  U.chisq <- L01 - f$deviance[2]
  p.U <- 1 - pchisq(U.chisq, 2)
  U <- (U.chisq - 2)/n
  Q <- D - U
  Dxy <- stats["Dxy"]
  C <- stats["C"]
  R2 <- stats["R2"]
  B <- sum((p - y)^2)/n
  # ES 15dec08 add Brier scaled
  Bmax  <- mean(y) * (1-mean(y))^2 + (1-mean(y)) * mean(y)^2
  Bscaled <- 1 - B/Bmax
  stats <- c(Dxy, C, R2, D, lr, p.lr, U, U.chisq, p.U, Q, B,
             f2$coef[1], f$coef[2], emax, Bscaled)
  names(stats) <- c("Dxy", "C (ROC)", "R2", "D", "D:Chi-sq",
                    "D:p", "U", "U:Chi-sq", "U:p", "Q", "Brier", "Intercept",
                    "Slope", "Emax", "Brier scaled")
  if(smooth=="loess")
    stats <- c(stats, c(Eavg = eavg),c(ECI = ECI))
  
  # Cut off definition
  if(!missing(cutoff)) {
    arrows(x0=cutoff,y0=.1,x1=cutoff,y1=-0.025,length=.15)
  }
  if(pl) {
    if(min(p)>plogis(-7) | max(p)<plogis(7)){
      
      lrm(y[i.2]~qlogis(p[i.2]))-> lrm.fit.1
      if(logistic.cal)  lines(p[i.2],plogis(lrm.fit.1$linear.predictors),lwd=lwd.log,lty=lty.log,col=col.log)
      
    }else{logit <- seq(-7, 7, length = 200)
    prob <- 1/(1 + exp( - logit))
    pred.prob <- f$coef[1] + f$coef[2] * logit
    pred.prob <- 1/(1 + exp( - pred.prob))
    if(logistic.cal) lines(prob, pred.prob, lty=lty.log,lwd=lwd.log,col=col.log)
    }
    #	pc <- rep(" ", length(lt))
    #	pc[lt==0] <- "."
    lp <- legendloc
    if (!is.logical(lp)) {
      if (!is.list(lp))
        lp <- list(x = lp[1], y = lp[2])
      legend(lp, leg, lty = lt, pch = marks, cex = cex.leg, bty = "n",lwd=lw.d,
             col=all.col,y.intersp = y.intersp)
    }
    if(!is.logical(statloc)) {
      if(dostats[1]==T){
        stats.2 <- paste('Calibration\n',
                         '...intercept: '
                         , sprintf(paste("%.",roundstats,"f",sep=""), stats["Intercept"]), " (",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.interc[1])," to ",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.interc[2]),")",'\n',
                         '...slope: '
                         , sprintf(paste("%.",roundstats,"f",sep=""), stats["Slope"]), " (",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.slope[1])," to ",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.slope[2]),")",'\n',
                         'Discrimination\n',
                         '...c-statistic: '
                         , sprintf(paste("%.",roundstats,"f",sep=""), stats["C (ROC)"]), " (",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.auc[2])," to ",
                         sprintf(paste("%.",roundstats,"f",sep=""),cl.auc[3]),")"
                         , sep = '')
        text(statloc[1], statloc[2],stats.2,pos=4,cex=cex)
        
      }else{
        dostats <- dostats
        leg <- format(names(stats)[dostats])	#constant length
        leg <- paste(leg, ":", format(stats[dostats], digits=roundstats), sep =
                       "")
        if(!is.list(statloc))
          statloc <- list(x = statloc[1], y = statloc[2])
        text(statloc, paste(format(names(stats[dostats])),
                            collapse = "\n"), adj = 0, cex = cex)
        text(statloc$x + (xlim[2]-xlim[1])/3 , statloc$y, paste(
          format(round(stats[dostats], digits=roundstats)), collapse =
            "\n"), adj = 1, cex = cex)
      }
    }
    if(is.character(riskdist)) {
      if(riskdist == "calibrated") {
        x <- f$coef[1] + f$coef[2] * log(p/(1 - p))
        x <- 1/(1 + exp( - x))
        x[p == 0] <- 0
        x[p == 1] <- 1
      }
      else x <- p
      bins <- seq(0, min(1,max(xlim)), length = 101)
      x <- x[x >= 0 & x <= 1]
      #08.04.01,yvon: distribution of predicted prob according to outcome
      f0	<-table(cut(x[y==0],bins))
      f1	<-table(cut(x[y==1],bins))
      j0	<-f0 > 0
      j1	<-f1 > 0
      bins0 <-(bins[-101])[j0]
      bins1 <-(bins[-101])[j1]
      f0	<-f0[j0]
      f1	<-f1[j1]
      maxf <-max(f0,f1)
      f0	<-(0.1*f0)/maxf
      f1	<-(0.1*f1)/maxf
      
      segments(bins1,line.bins,bins1,length.seg*f1+line.bins)
      segments(bins0,line.bins,bins0,length.seg*-f0+line.bins)
      lines(c(min(bins0,bins1)-0.01,max(bins0,bins1)+0.01),c(line.bins,line.bins))
      text(max(bins0,bins1)+dist.label,line.bins+dist.label2,d1lab,cex=cex.d01)
      text(max(bins0,bins1)+dist.label,line.bins-dist.label2,d0lab,cex=cex.d01)
      
    }
  }
  if(dostats==T){
    cat(paste("\n\n A ",cl.level*100,
              "% confidence interval is given for the calibration intercept, calibration slope and c-statistic. \n\n",
              sep=""))}
  stats
}


getAUCmw <- function(x, y){
  
  xy <- expand.grid(x, y)
  mean(ifelse(xy[,2] < xy[,1], 1, (ifelse(xy[,2] == xy[,1], 1/2, 0))))
  
}

#get the variance V_2(\theta) from method 4 of Newcombe 2006
getVARmw  <- function(theta, m, n){
  nstar <- mstar <- (m+n)/2  - 1
  theta *
    (1-theta) *
    (1 + nstar*(1-theta)/(2-theta) + mstar*theta/(1+theta)) /
    (m*n)
  
}

#note that here x < y
#m is the length of x and n is the length of y
getS2mw <- function(x, y, m, n){
  
  xi <- sort(x)
  yj <- sort(y)
  rxy <- rank(c(xi, yj))
  Ri <- rxy[1:m]
  Sj <- rxy[m+1:n]
  S102 <- 1/((m-1)*n^2) * (sum((Ri-1:m)^2) - m*(mean(Ri)-(m+1)/2)^2)
  S012 <- 1/((n-1)*m^2) * (sum((Sj-1:n)^2) - n*(mean(Sj)-(n+1)/2)^2)
  S2 <- (m*S012 + n*S102) / (m + n)
  S2
  
}


#get the CI by method 5 of Newcombe 2006
getAUCCImwu  <- function(hat.theta, zalpha, m, n){
  
  nstar <- mstar <- (m+n)/2  - 1
  
  a <- -1 - nstar - mstar
  b <- 1 + 2*mstar
  c <- 2 + nstar
  d <- 1 + 2*hat.theta
  ht2 <- hat.theta^2
  e <- 2 - 2*hat.theta - ht2
  f <- ht2 - 4*hat.theta
  g <- 2*ht2
  
  z2 <- zalpha^2
  mn <- m*n
  z5 <- -mn  + z2*a
  z4 <- mn*d - z2*(a-b)
  z3 <- mn*e - z2*(b-c)
  z2 <- mn*f - z2*c
  z1 <- mn*g
  
  roots <- polyroot(c(z1, z2, z3, z4, z5))
  real <- Re(roots[sapply(1:4, function(i) all.equal(Im(roots[i]), 0))])
  real <- real[real > 0 & real < 1]
  ci <- sort(real)
  if(length(ci) > 2)
    warning("There are three roots meet the requirement when computing
            the confidence interval.")
  else{
    if(length(real) == 1){
      if(real <= hat.theta){
        ci <- c(real, 1)
      }
      else{
        ci <- c(0, real)
      }
    }
    if(length(real) == 0){
      ci <- c(0, 1)
    }
  }
  ci
}



#get the CI by method 5 of Newcombe 2006
auc.mw.newcombe <- function(x, y, alpha){
  
  point <- getAUCmw(x, y)
  nx <- length(x)
  ny <- length(y)
  zalpha <- qnorm(1-alpha/2)
  ci <- getAUCCImwu(point, zalpha, ny, nx)
  
  c(point, ci)
}

auc.mw.zhou <- function(x, y, alpha){
  if(max(y) < min(x)){
    c(1, 1, 1)
  }
  else{
    point <- getAUCmw(x, y)
    nx <- length(x)
    ny <- length(y)
    zalpha <- qnorm(1-alpha/2)
    varHatTheta <- getVARmw(point,nx, ny)
    
    Z <- 1/2 * log((1+point)/(1-point))
    varZ <- 4 / (1-point^2)^2 * varHatTheta
    LL <- Z - zalpha*sqrt(varZ)
    UL <- Z + zalpha*sqrt(varZ)
    ci <- c((exp(2*LL) - 1)/ (exp(2*LL) + 1), (exp(2*UL) - 1)/ (exp(2*UL) + 1))
    
    c(point,ci)
  }
}

auc.mw.pepe <- function(x, y, alpha){
  
  if(max(y) < min(x)){
    c(1, 1, 1)
  }
  else{
    point <- getAUCmw(x, y)
    nx <- length(x)
    ny <- length(y)
    zalpha <- qnorm(1-alpha/2)
    
    varHatTheta <- (nx+ny) * getS2mw(y, x, ny, nx) / (nx*ny)
    #varHatTheta <- getVARmw(point,nx, ny)
    LL <- log(point/(1-point)) - zalpha*sqrt(varHatTheta)/(point*(1-point))
    UL <- log(point/(1-point)) + zalpha*sqrt(varHatTheta)/(point*(1-point))
    ci <- c(exp(LL) / (1+exp(LL)), exp(UL) / (1+exp(UL)))
    
    c(point,ci)
  }
}

auc.mw.delong <- function(x, y, alpha){
  
  point <- getAUCmw(x, y)
  nx <- length(x)
  ny <- length(y)
  zalpha <- qnorm(1-alpha/2)
  
  D10 <- sapply(1:ny, function(i)
    mean(ifelse(x > y[i], 1, ifelse(x == y[i], 1/2, 0))))
  D01 <- sapply(1:nx, function(i)
    mean(ifelse(x[i] > y, 1, ifelse(x[i] == y, 1/2, 0))))
  varDhatTheta <- 1/(ny*(ny-1))*sum((D10-point)^2) +
    1/(nx*(nx-1))*sum((D01-point)^2)
  ci <- c(point - zalpha*sqrt(varDhatTheta),
          point + zalpha*sqrt(varDhatTheta))
  
  c(point, ci)
}

mw.jackknife <- function(x, y){
  nx <- length(x)
  ny <- length(y)
  n <- nx + ny
  
  hatThetaPartial <- rep(0, n)
  
  for(i in 1:nx){
    hatThetaPartial[i] <- getAUCmw(x[-i], y)
  }
  for(i in 1:ny){
    hatThetaPartial[i+nx] <- getAUCmw(x, y[-i])
  }
  
  hatThetaPartial
  
}

auc.mw.jackknife <- function(x, y, alpha){
  
  nx <- length(x)
  ny <- length(y)
  n <- nx + ny
  
  hatTheta <- getAUCmw(x, y)
  
  hatThetaPseudo <- rep(0, n)
  
  hatThetaPartial <- mw.jackknife(x, y)
  
  for(i in 1:n){
    hatThetaPseudo[i] <- n*hatTheta - (n-1)*hatThetaPartial[i]
  }
  
  point <- mean(hatThetaPseudo)
  ST2 <- mean((hatThetaPseudo - point)^2) / (n-1)
  ST <- sqrt(ST2)
  z.alpha2 <- qt(1 - alpha/2, df=n-1)
  ci <- c(point - z.alpha2*ST, point + z.alpha2*ST)
  
  c(point, ci)
}

auc.mw.boot <- function(x, y, alpha, nboot=1000, method){
  if(max(y) < min(x)){
    c(1, 1, 1)
  }
  else{
    nx <- length(x)
    ny <- length(y)
    
    point <- getAUCmw(x, y)
    index.x <- matrix(sample.int(nx, size = nx*nboot, replace = TRUE),
                      nboot, nx)
    index.y <- matrix(sample.int(ny, size = ny*nboot, replace = TRUE),
                      nboot, ny)
    mw.boot <- sapply(1:nboot, function(i) getAUCmw(x[index.x[i,]],
                                                    y[index.y[i,]]))
    if(method=="P"){
      ci <- as.vector(quantile(mw.boot, c(alpha/2, 1-alpha/2), type=6))
    }
    else{
      hatZ0 <- qnorm(mean(mw.boot < point))
      
      partial <- mw.jackknife(x, y)
      mpartial <- mean(partial)
      hatA <- sum((mpartial - partial)^3) /
        (6 * (sum((mpartial - partial)^2))^(3/2))
      
      alpha1 <- pnorm(hatZ0 + (hatZ0 + qnorm(alpha/2)) /
                        (1 - hatA*(hatZ0 + qnorm(alpha/2))))
      alpha2 <- pnorm(hatZ0 + (hatZ0 + qnorm(1-alpha/2)) /
                        (1 - hatA*(hatZ0 + qnorm(1-alpha/2))))
      
      ci <- as.vector(quantile(mw.boot, c(alpha1, alpha2), type=6))
    }
    c(point, ci)
  }
}

auc.nonpara.mw <- function(x, y, conf.level=0.95,
                           method=c("newcombe", "pepe", "delong", "jackknife", "bootstrapP", "bootstrapBCa"),
                           nboot){
  
  alpha <- 1 - conf.level
  method <- match.arg(method)
  estimate <- switch(method,
                     newcombe=auc.mw.newcombe(x, y, alpha),
                     pepe=auc.mw.pepe(x, y, alpha),
                     delong=auc.mw.delong(x, y, alpha),
                     jackknife=auc.mw.jackknife(x, y, alpha),
                     bootstrapP=auc.mw.boot(x, y, alpha, nboot, method="P"),
                     bootstrapBCa=auc.mw.boot(x, y, alpha, nboot, method="BCa"))
  estimate
}

cat("AUC(xb.hat,y)", "\n")

AUC <- function(xb.hat,y){
  n<-length(xb.hat)
  n1<-sum(y)
  mean.rank <- mean(rank(xb.hat)[y == 1])
  AUC<-(mean.rank - (n1 + 1)/2)/(n - n1)
  return(AUC)
}



ci.auc <- function(crit,pred,conf.level=0.95,method="pepe"){
  
  tmp <- cbind.data.frame(crit=crit,pred=pred)
  # healthy
  nondis <- tmp[which(tmp$crit==0),]
  
  # disease
  dis <- tmp[which(tmp$crit==1),]
  
  # ci auc
  if(!grepl("bootstrap",method)){
    result <- auc.nonpara.mw(dis$pred,nondis$pred,conf.level,method)
  }else if(grepl("bootstrap",method)){
    cat("Bootstrap-based methods are not supported by this package. Method will be set to 'pepe'. \n\n")
    result <- auc.nonpara.mw(dis$pred,nondis$pred,conf.level,method="pepe")
  }
  
  return(result)
}


# ------------------------------ validation plot ---------------------------------------

# the val.prob.ci.2 function is obtained from E. Steyerberg 

par(mfrow=c(1,1))
val.prob.ci.2(logit=mod02DEF$linear.predictors, y=mod02DEF$y, smooth = "rcs", col.smooth = "cornflowerblue", lwd.smooth=2, CL.smooth = "fill", dostats = TRUE, logistic.cal = TRUE, riskdist = "predicted",
              xlim = c(-0.05,.4), ylim=c(-0.15,.8), connect.smooth = TRUE, statloc = c(-.05,.65), pl=TRUE, xlab="Predicted probability according to updated model"
              , ylab="Observed proportion in the validation cohort", lwd.ideal=1, col.ideal="black", legendloc = c(-.05,.5), lty.smooth = 1, col.log="orange", cex.d01 = 1)
title("Predictive performance after model updating")
