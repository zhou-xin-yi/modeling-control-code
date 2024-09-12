### R code from vignette source 'mle_2012.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: opts
###################################################
x11(width=10,height=6) ## hack for figure margins too small error
options(continue=" ")


###################################################
### code chunk number 2: double (eval = FALSE)
###################################################
## options(SweaveHooks=list( fig=function() par(mfrow=c(1,2), bty="l", las=1) ))


###################################################
### code chunk number 3: quad (eval = FALSE)
###################################################
## options(SweaveHooks=list( fig=function() { par(mfrow=c(2,2), bty="l", las=1)} ))


###################################################
### code chunk number 4: single
###################################################
## DM original: 5,13,4,11
options(SweaveHooks=list( fig=function() { par(mar=c(5,10,4,10)+0.1, bty="l", las=1) }))


###################################################
### code chunk number 5: readgdat
###################################################
dat <- read.table("gophertortoise.txt",header=TRUE)


###################################################
### code chunk number 6: calcs1
###################################################
tt <- with(dat,table(size,status))
tt0 <- with(subset(dat,size>=230 & size<240),table(status))
neg0 <- tt0[1]
pos0 <- tt0[2]
tot0 <- sum(tt0)
with(subset(dat,size>=230 & size<240),table(Sex,status))
with(subset(dat,size>=230 & size<240),table(Sex,status))


###################################################
### code chunk number 7: binomfig1
###################################################
xvec = 0:tot0
plot(xvec,dbinom(xvec,size=tot0,prob=0.84),
     xlab="Number seropositive (x)",ylab="Probability")
abline(v=pos0)


###################################################
### code chunk number 8: mle_2012.Rnw:178-179
###################################################
extremes <- qbinom(c(0.01,0.99),prob=0.84,size=tot0)


###################################################
### code chunk number 9: binomfig2
###################################################
plot(xvec,dbinom(xvec,size=tot0,prob=0.84,log=TRUE),
     xlab="Number seropositive (x)",ylab="Log probability")
abline(v=pos0)


###################################################
### code chunk number 10: mle_2012.Rnw:198-199 (eval = FALSE)
###################################################
## plot(xvec,dbinom(xvec,size=tot0,prob=0.84),log="y")


###################################################
### code chunk number 11: likcurve
###################################################
options(SweaveHooks=list( fig=function() par(mfrow=c(1,2), bty="l", las=1) ))
xlabstr <- "Per capita probability seropositive (p)" 
curve(dbinom(pos0,prob=x,size=tot0),
      from=0,to=1,
      xlab=xlabstr,ylab="Likelihood")
abline(v=pos0/tot0,lty=2)
mtext(side=3,at=pos0/tot0,expression(hat(p)))
abline(h=dbinom(pos0,size=tot0,prob=pos0/tot0),lty=2)
curve(-dbinom(pos0,prob=x,size=tot0,log=TRUE),
      from=0,to=1,ylim=c(0,30),
      xlab=xlabstr,ylab="Negative log-likelihood")
abline(v=pos0/tot0,lty=2)
mtext(side=3,at=pos0/tot0,expression(hat(p)))
abline(h=-dbinom(pos0,size=tot0,prob=pos0/tot0,log=TRUE),lty=2)


###################################################
### code chunk number 12: aggreg
###################################################
library(plyr)
sizetab <- ddply(dat,"size",
                function(x) c(tot=nrow(x),pos=sum(x$status)))


###################################################
### code chunk number 13: mle_2012.Rnw:319-320
###################################################
(prob0 <- with(sizetab,sum(pos)/sum(tot)))


###################################################
### code chunk number 14: comblik
###################################################
with(sizetab,-sum(dbinom(pos,size=tot,prob=prob0,log=TRUE)))


###################################################
### code chunk number 15: likfun
###################################################
NLLfun_binom0 <- function(prob,dat=sizetab) {
  with(dat,-sum(dbinom(pos,size=tot,prob=prob,log=TRUE)))
}


###################################################
### code chunk number 16: mle_2012.Rnw:339-341
###################################################
pvec <- seq(0,1,length=201)
NLLvec <- sapply(pvec,NLLfun_binom0)


###################################################
### code chunk number 17: vectorfun (eval = FALSE)
###################################################
## ## slick vectorization stuff to hurt your brain if you look in here ...
## vNLLfun <- Vectorize(NLLfun_binom0,vectorize.args="prob")
## curve(vNLLfun_binom0(x),from=0,to=1,ylim=c(10,400))
## tdat <- data.frame(tot=31,pos=26)
## abline(v=prob0,lty=2)
## curve(vNLLfun_binom0(x,dat=tdat),add=TRUE,col=2)
## abline(v=pos0/tot0,lty=2,col=2)
## minNLL0 <- min(NLLvec)
## minNLL1 <- NLLfun_binom0(pos0/tot0,tdat)


###################################################
### code chunk number 18: mle_2012.Rnw:357-358
###################################################
## DM original: 5,13,4,11
options(SweaveHooks=list( fig=function() { par(mar=c(5,10,4,10)+0.1, bty="l", las=1) }))


###################################################
### code chunk number 19: likcurve2
###################################################
plot(pvec,NLLvec,type="l",
     xlab=xlabstr,ylab="Negative log likelihood")


###################################################
### code chunk number 20: exercise1 (eval = FALSE)
###################################################
## NLLfun_binom0 <- function(prob,dat=sizetab) {
##    with(dat,-sum(dbinom(pos,size=tot,prob=prob,log=TRUE)))
##  }
## pvec <- seq(0,1,length=201)
## NLLvec <- sapply(pvec,NLLfun_binom0)
## plot(pvec,NLLvec,type="l",
##       xlab=xlabstr,ylab="Negative log likelihood",ylim=c(0,800))
## NLLfun_binom0 <- function(prob,dat=sizetab[50:59,]) {
##    with(dat,-sum(dbinom(pos,size=tot,prob=prob,log=TRUE)))
##  }
## pvec <- seq(0,1,length=201)
## NLLvec2 <- sapply(pvec,NLLfun_binom0)
## lines(pvec,NLLvec2,col=2)
## 
## plot(pvec,NLLvec-min(NLLvec),type="l",
##       xlab=xlabstr,ylab="Adjusted Negative log likelihood",ylim=c(0,200))
## lines(pvec,NLLvec2-min(NLLvec2),col=2)


###################################################
### code chunk number 21: optim-call
###################################################
m0<-optim(par=05,fn=NLLfun_binom0,method="Brent",lower=0,upper=1)
m0


###################################################
### code chunk number 22: mle_2012.Rnw:457-459
###################################################
library(bbmle)
(m1 <- mle2(NLLfun_binom0,start=list(prob=0.5)))


###################################################
### code chunk number 23: mlefit2
###################################################
m1F <- mle2(pos~dbinom(prob=prob,size=tot),
            start=list(prob=0.5),data=sizetab)


###################################################
### code chunk number 24: constplot
###################################################
with(sizetab,plot(size,pos/tot,cex=tot,
                 xlab="size",ylab="status/fraction seropositive"))
abline(h=prob0,col=2)


###################################################
### code chunk number 25: mle_2012.Rnw:547-552
###################################################
## does bobyqa work?
#library(optimx)
#mle2(pos~dbinom(prob=prob,size=tot),start=list(prob=0.5),data=sizetab,
#     optimizer="optimx",method="bobyqa",lower=0,upper=1)
#mle2(pos~dbinom(prob=plogis(logitprob),size=tot),start=list(logitprob=qlogis(0.5)),data=sizetab)


###################################################
### code chunk number 26: likcurve2B
###################################################
plot(pvec,NLLvec,type="l",
     xlab=xlabstr,ylab="Negative log likelihood",ylim=c(114,120),
     xlim=c(0.6,0.9))
abline(h=c(min(NLLvec),min(NLLvec)+1.92),col=2,lty=2)
cc <- confint(m1,quietly=TRUE)
## ccq <- confint(m1,method="quad")
## need quadratic approx
## curve(min(NLLvec)+((x-coef(m1))/(stdEr(m1)))^2/2,add=TRUE,col=4)
abline(v=cc)
## abline(v=ccq,col=4,lty=2)
## legend("topright",c("profile","quad. approx"),
##       lty=1,col=c(1,4),xpd=NA)


###################################################
### code chunk number 27: uppertail
###################################################
(s <- 2*(c(logLik(m1))-(-NLLfun_binom0(0.5))))  ## test statistic


###################################################
### code chunk number 28: pval
###################################################
pchisq(s,df=1,lower.tail=FALSE)


###################################################
### code chunk number 29: cfun
###################################################
cfun <- function(x,min=0.001,max=0.999) { 
  ifelse(x<min,min,ifelse(x>max,max,x))
}


###################################################
### code chunk number 30: m1L
###################################################
(m1L <- mle2(pos~dbinom(prob=cfun(a+b*(size-200)),
                       size=tot),start=list(a=0.5,b=0.001),data=sizetab))


###################################################
### code chunk number 31: pred1
###################################################
pp1 <- data.frame(size=70:310,tot=1)
pp1$pred <- predict(m1L,pp1)


###################################################
### code chunk number 32: linpredplot
###################################################
with(sizetab,plot(size,pos/tot,cex=tot,
                 xlab="size",ylab="status/fraction seropositive"))
with(pp1,lines(size,pred,col=2))


###################################################
### code chunk number 33: mle_2012.Rnw:694-702
###################################################
## we could define the negative log-likelihood function:
NLLfun_binom_lin <- function(a,b,dat=sizetab) {
  with(dat,
       {
         prob <-cfun(a+b*(size-200))
         -sum(dbinom(pos,size=tot,prob=prob,log=TRUE))
     })
}


###################################################
### code chunk number 34: mle_2012.Rnw:717-725
###################################################
## hack to get useRaster=TRUE without having to explain it
library(emdbook)
cc <- curve3d(m1L@minuslogl(x,y),
              xlim=c(0.45,0.65),ylim=c(0.002,0.006),
              sys3d="image", xlab="intercept",ylab="slope",
              useRaster=TRUE)
contour(cc$x,cc$y,cc$z,add=TRUE)
points(coef(m1L)[1],coef(m1L)[2],pch=16)


###################################################
### code chunk number 35: mle_2012.Rnw:728-734 (eval = FALSE)
###################################################
## library(emdbook)
## cc <- curve3d(m1L@minuslogl(x,y),
##               xlim=c(0.45,0.65),ylim=c(0.002,0.006),
##               sys3d="image", xlab="intercept",ylab="slope")
## contour(cc$x,cc$y,cc$z,add=TRUE)
## points(coef(m1L)[1],coef(m1L)[2],pch=16)


###################################################
### code chunk number 36: mle_2012.Rnw:748-749
###################################################
options(SweaveHooks=list( fig=function() par(mfrow=c(1,2), bty="l", las=1) ))


###################################################
### code chunk number 37: prof2
###################################################
plot(profile(m1L))


###################################################
### code chunk number 38: mle_2012.Rnw:758-759
###################################################
confint(m1L)


###################################################
### code chunk number 39: mle_2012.Rnw:771-773
###################################################
m1L2 <- mle2(pos~dbinom(prob=plogis(a+b*(size-200)),
                       size=tot),start=list(a=0.5,b=0.001),data=sizetab)


###################################################
### code chunk number 40: mle_2012.Rnw:779-793
###################################################
## DM original: 5,13,4,11
options(SweaveHooks=list( fig=function() { par(mar=c(5,10,4,10)+0.1, bty="l", las=1) }))
NLLfun_binom_logist <- function(a,b,dat=sizetab) {
  with(dat,
       {
         prob <-plogis(a+b*(size-200))
         -sum(dbinom(pos,size=tot,prob=prob,log=TRUE))
     })
}
cc2 <- curve3d(m1L2@minuslogl(x,y),xlim=c(-1,1),ylim=c(0.01,0.05),
              sys3d="image",useRaster=TRUE,ann=FALSE)
mtext(side=1,line=2.5,"intercept",cex=1.5)
mtext(side=2,line=3.5,"slope",las=0,cex=1.5)
contour(cc2$x,cc2$y,cc2$z,add=TRUE)
points(coef(m1L2)[1],coef(m1L2)[2],pch=16)


###################################################
### code chunk number 41: mle_2012.Rnw:798-799
###################################################
options(SweaveHooks=list( fig=function() par(mfrow=c(1,2), bty="l", las=1) ))


###################################################
### code chunk number 42: prof2L
###################################################
plot(profile(m1L2))


###################################################
### code chunk number 43: mle_2012.Rnw:806-808
###################################################
pp2 <- pp1
pp2$pred <- predict(m1L2,pp2)


###################################################
### code chunk number 44: logistpredplot
###################################################
## DM original: 5,13,4,11
options(SweaveHooks=list( fig=function() { par(mar=c(5,10,4,10)+0.1, bty="l", las=1) }))
with(sizetab,plot(size,pos/tot,cex=tot,
                 xlab="size",ylab="status/fraction seropositive"))
with(pp1,lines(size,pred,col=2))
with(pp2,lines(size,pred,col=4))
legend(200,0.3,c("linear","logistic"),col=c(2,4),lty=1)


###################################################
### code chunk number 45: mle_2012.Rnw:822-827
###################################################
## this is not used
m1L3 <- mle2(pos~dbinom(prob=cfun(1-exp(b*(size-size0))),
                        size=tot),start=list(size0=90,b=0.001),data=sizetab)
pp3 <- pp1
pp3$pred <- predict(m1L3,pp3)


###################################################
### code chunk number 46: mle_2012.Rnw:850-851
###################################################
anova(m1L2,m1)


###################################################
### code chunk number 47: mle_2012.Rnw:884-885
###################################################
AICtab(m1,m1L,m1L2)


###################################################
### code chunk number 48: exercise-fit3 (eval = FALSE)
###################################################
## 
## m1L3 <- mle2(pos~dbinom(prob=cfun(1-exp(b*(size-a))),
##                        size=tot),start=list(a=50,b=-0.001),data=sizetab)
## 
## with(sizetab,plot(size,pos/tot,cex=tot,
##                  xlab="size",ylab="status/fraction seropositive"))
## abline(h=coef(m1),lty=2)                 
## lines(pp1$size,predict(m1L,pp1))                 
## lines(pp1$size,predict(m1L2,pp1),col=2)                 
## lines(pp1$size,predict(m1L3,pp1),col=4)  
## AICtab(m1,m1L,m1L2,m1L3)
## 


###################################################
### code chunk number 49: mle_2012.Rnw:935-936
###################################################
dev.off()


