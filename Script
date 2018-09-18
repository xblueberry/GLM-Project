#Generalized linear models final project
#June 2017
#Authors: Daniel Izquierdo Juncas, Nozomi Takemura

#libraries
library(GGally)
library(ggplot2)   
library(xtable)
library(MASS)
library(car)
library(lmtest)
library(sandwich)
library(coda)
library(MCMCpack)
library(zoo)
library(gpairs)
library(xtable)
library(MASS)
library(car)
library(lmtest)
library(sandwich)
library(MCMCpack)


Roadkills0<-read.table(choose.files(),header=T)
Roadkills<- subset(Roadkills0, select=c("TOT.N", "OPEN.L","MONT.S","POLIC","D.PARK",
                                        "SHRUB","WAT.RES","L.WAT.C","L.P.ROAD","D.WAT.COUR"))
View(Roadkills)

attach(Roadkills)
n<-nrow(Roadkills)

####Exploratory analysis####

str(Roadkills)
summary(Roadkills)
mean(TOT.N); var(TOT.N)
#stat.desc(Roadkills) #require(pastecs)
#xtable(stat.desc(Roadkills))
#mean(y)<<Var(y)

pairs(Roadkills)
gpairs(Roadkills)

####online code--- for scatter plot---####
 # Function to return points and geom_smooth
# allow for the method to be changed
my_fn <- function(data, mapping, method="loess", ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=method, ...)
  p
}

# Default loess curve    
#ggpairs(Roadkills[1:10], lower = list(continuous = my_fn))
#chart.Correlation(Roadkills, method="spearman",histogram=TRUE) #require(PerformanceAnalytics)

#Boxplots of response and predictor variables
par(mfrow=c(2,5))
boxplot(TOT.N,main="Number of amphibian fatalities")
boxplot(OPEN.L,main="Open lands (ha)")
boxplot(MONT.S,main="Montado with shrubs (ha)")
boxplot(POLIC,main="Policulture (ha)")
boxplot(D.PARK,main="Distance to Natural Park (m)")
boxplot(SHRUB,main="Shrubs (ha)")
boxplot(WAT.RES,main="Water reservoirs")
boxplot(L.WAT.C,main="Length of water courses (km)")
boxplot(L.P.ROAD,main="Paved road length (km)")
boxplot(D.WAT.COUR,main="Distance to water courses")
par(mfrow=c(1,1))

#Histogram of response
hist(TOT.N,nclas=10,col="light blue",xlab="Amphibian fatalities per segment",ylab="Frequency",main="")
#resembles a Poisson distribution


#####Multicollinearity diagnostic####
Roadkills.cor<-cor(Roadkills)
corx <- Roadkills.cor[-1,-1]
corx
Roadkills.VIF <- diag(solve(corx))#inverse, diagonal.
Roadkills.VIF
# It can be shown that VIFj equals the jth diagonal element of the inverse correlation
# matrix of X.
# VIFj = (R 
#         ???1
#         XX)jj 

xtable(as.data.frame(t(Roadkills.VIF)))
#VIFs small, no multicollinearity issues


####Poisson regression: model selection####

#Check poisson distribution
hist(TOT.N, nclas=10,col="light blue",prob=T,ylim=c(0,0.08),
     xlab="Amphibian fatalities per segment",ylab="Frequency",main="")
poismodel0 <- glm(TOT.N ~ 1,family = poisson,data = Roadkills) 
summary(poismodel0)
poissonmu <- exp(summary(poismodel0)$coefficients[,1])
lines(0:105,dpois(0:105,poissonmu),col="blue",lwd=3) #doesn't fit good the data
#Fit null model
#poismodel1 <- glm(TOT.N ~ 1,family = poisson(link="log"),data = Roadkills) 
#summary(poismodel1)


#Fit full (linear) model
poismodel1 <- glm(TOT.N ~ .,family = poisson(link="log"),data = Roadkills) 
summary(poismodel1)
xtable(poismodel1)

#Model selection
drop1(poismodel1,test="Chisq") #(Chisq=LRT)
drop1(poismodel1,test="Rao") #score

poismodel2 <- glm(TOT.N ~ .- D.WAT.COUR,family = poisson(link="log"),data = Roadkills) 
drop1(poismodel2,test="Chisq")
drop1(poismodel2,test="Rao")

#stepwise with AIC
step(poismodel0, scope = list(lower=poismodel0,upper=poismodel1),direction="both",data=Roadkills)

summary(poismodel2)
xtable(poismodel2)


#Fit---interaction--
poismodel.int.full <- glm(TOT.N~.*.,family = poisson(link="log"),data = Roadkills)
step(poismodel0, scope = list(lower=poismodel0,upper=poismodel.int.full),direction="both",data=Roadkills)
poismodel.int.from.null <- glm(formula = TOT.N ~ D.PARK + L.WAT.C + MONT.S + SHRUB + OPEN.L + 
                                 POLIC + WAT.RES + L.P.ROAD + L.WAT.C:MONT.S + MONT.S:OPEN.L + 
                                 L.WAT.C:L.P.ROAD + SHRUB:L.P.ROAD + OPEN.L:WAT.RES + SHRUB:POLIC + 
                                 D.PARK:L.P.ROAD + OPEN.L:L.P.ROAD + D.PARK:MONT.S, family = poisson, 
                               data = Roadkills)
step(poismodel.int.full, scope = list(lower=poismodel0,upper=poismodel.int.full),direction="both",data=Roadkills)
poismodel.int.from.full <- glm(formula = TOT.N ~ OPEN.L + MONT.S + POLIC + D.PARK + SHRUB + 
                                 WAT.RES + L.WAT.C + L.P.ROAD + D.WAT.COUR + OPEN.L:MONT.S + 
                                 OPEN.L:POLIC + OPEN.L:SHRUB + OPEN.L:L.P.ROAD + MONT.S:L.WAT.C + 
                                 MONT.S:L.P.ROAD + POLIC:SHRUB + POLIC:WAT.RES + POLIC:D.WAT.COUR + 
                                 D.PARK:SHRUB + D.PARK:L.P.ROAD + SHRUB:L.WAT.C + SHRUB:L.P.ROAD + 
                                 WAT.RES:D.WAT.COUR + L.WAT.C:D.WAT.COUR + L.P.ROAD:D.WAT.COUR, 
                               family = poisson(link = "log"), data = Roadkills)
summary(poismodel.int.from.full)
int.sum <- summary(poismodel.int.from.null)
xtable(int.sum$coefficients,digits = 4)
#Sqrt link
poismodel.sqrt0 <- glm(TOT.N ~ 1,family = poisson(link="sqrt"),data = Roadkills)
poismodel.sqrt1 <- glm(TOT.N ~ .,family = poisson(link="sqrt"),data = Roadkills)
step(poismodel.sqrt0, scope = list(lower=poismodel.sqrt0,upper=poismodel1),direction="both",data=Roadkills)
poismodel.sqrt2 <- glm(TOT.N ~ .-D.WAT.COUR-POLIC,family = poisson(link="sqrt"),data = Roadkills)
summary(poismodel.sqrt2)


####Poisson regression: diagnostics####

#Deviance
dev.pois<-summary(poismodel2)$deviance
df.pois<-summary(poismodel2)$df.residual
dev.pois/df.pois
1-pchisq(dev.pois,df.pois)
#residual deviance>>df, there's overdispersion
with(summary(poismodel2), cbind(res.deviance = deviance, 
              df = df.residual,p = pchisq(deviance, df.residual, lower.tail=FALSE)))
# res.deviance df            p
# [1,]     273.1239 43 3.101038e-35

par(mfrow=c(2,2))
# par(mfrow=c(1,1))
plot(poismodel2)

#Deviance residual plot (INCLUDE?)
par(mfrow=c(1,1))
r.pois <- residuals(poismodel2, type = "deviance")
#OPEN.L
plot(OPEN.L,r.pois,xlab="OPEN.L",ylab="Deviance residual",
     cex.lab=1.5,cex.axis=1.3)
loess.pois <- loess(r.pois~OPEN.L)
lo.pred.pois <- predict(loess.pois, se=T)
orderOPEN.L <- order(OPEN.L)
lines(OPEN.L[orderOPEN.L ],lo.pred.pois$fit[orderOPEN.L ],col="blue",lwd=3)
lines(OPEN.L[orderOPEN.L ],lo.pred.pois$fit[orderOPEN.L ]+2*lo.pred.pois$s[orderOPEN.L ], lty=2,col="red")
lines(OPEN.L[orderOPEN.L ],lo.pred.pois$fit[orderOPEN.L ]-2*lo.pred.pois$s[orderOPEN.L ], lty=2,col="red")
#L.WAT.C
plot(L.WAT.C,r.pois,xlab="L.WAT.C",ylab="Deviance residual",
     cex.lab=1.5,cex.axis=1.3, main = "Deviance residual plot (L.WAT.C)")
loess.pois <- loess(r.pois~L.WAT.C)
lo.pred.pois <- predict(loess.pois, se=T)
orderL.WAT.C <- order(L.WAT.C)
lines(L.WAT.C[orderL.WAT.C ],lo.pred.pois$fit[orderL.WAT.C ],col="blue",lwd=3)
lines(L.WAT.C[orderL.WAT.C ],lo.pred.pois$fit[orderL.WAT.C ]+2*lo.pred.pois$s[orderL.WAT.C ], lty=2,col="red")
lines(L.WAT.C[orderL.WAT.C ],lo.pred.pois$fit[orderL.WAT.C ]-2*lo.pred.pois$s[orderL.WAT.C ], lty=2,col="red")


#Deviance residual plot for L.WAT.C



#Influence plots
par(mfrow=c(2,2))
id <- 1:n
plot(hatvalues(poismodel2),rstudent(poismodel2))
plot(id,dffits(poismodel2),type="l")
plot(id,covratio(poismodel2),type="l")
plot(id,cooks.distance(poismodel2),type="l")


####Prediction--Comparison with poisson with interaction vs without--####
#for the poisson model with interaction
minage <- min(L.WAT.C,na.rm=TRUE)
maxage <- max(L.WAT.C,na.rm=TRUE)
grid <- seq(minage,maxage,0.1)
model.diet1.p <- data.frame(D.PARK=mean(D.PARK),
                            MONT.S=mean(MONT.S),
                            SHRUB=mean(SHRUB),
                            OPEN.L=mean(OPEN.L),
                            POLIC=mean(POLIC),
                            WAT.RES=mean(WAT.RES),
                            L.P.ROAD=mean(L.P.ROAD),L.WAT.C=grid)
p1 <- predict (poismodel.int.from.null,model.diet1.p,se=T,type = "response")
b <- p1$fit
LL <- p1$fit - 1.96 * p1$se.fit
UL <- p1$fit + 1.96 * p1$se.fit
a <- cbind(b,LL,UL)
#for poisson regression(no interaction) part
model <- glm(formula = TOT.N ~ OPEN.L+MONT.S+POLIC+D.PARK+SHRUB+WAT.RES+L.WAT.C+L.P.ROAD,
             family = poisson,data = test2)
p1 <- predict (model,model.diet1.p,se=T,type = "response")
b <-  p1$fit
LL <- p1$fit - 1.96 * p1$se.fit
UL <- p1$fit + 1.96 * p1$se.fit
c <- cbind(b,LL,UL)
#plot both in the one figure
par(mfrow=c(1,1))
par(pty="m")
par(mar = c(5,5,4,2)+0.1)
matplot(grid,a,lty=c(1,2,2),type="l",col=c(2,2,2),lwd=3,
        xlab="Length of water courses(km)",
        ylab="Predicted total number of amphibian fatalities",
        # ylab="Log of predicted amphibian fatalities",
        cex.lab=1.5,cex.axis=1.3 )
rug(Roadkills$L.WAT.C)
polygon(c(rev(grid), grid), c(rev(c[ ,3]), c[ ,2]), col = 'grey80', border = NA,density = c(50, 50) )
lines(grid, c[,1], lty=1, lwd=3,col=1)
lines(grid,c[,2], lty=4,lwd=3, col=1)
lines(grid,c[,3], lty=4,lwd=3, col=1)


####Negative binomial model: model selection####

#Visual inspection
hist(TOT.N, nclas=10,col="light blue",prob=T,ylim=c(0,0.08),
     xlab="Amphibian fatalities per segment",ylab="Frequency",main="")
lines(0:105,dpois(0:105,poissonmu),col="blue",lwd=2) #Poisson

ngbinmodel0 <- glm.nb(TOT.N ~ 1,data = Roadkills)
summary(ngbinmodel0)
nbmu <- exp(summary(ngbinmodel0)$coefficients[,1])
nbsize <- summary(ngbinmodel0)$theta
nbp <- 1-nbmu/(nbmu+nbsize)
ngbindmft <- dnbinom(0:105, prob=nbp, size=nbsize, log = FALSE)
lines(0:105,ngbindmft,col="red",lwd=2) #Negative binomial

legend("topright",lty = c(1,1),col = c("blue","red"),legend=c("Poisson","Negative binomial"))


#Fit full (linear) model
ngbinmodel1 <- glm.nb(TOT.N ~ .,data = Roadkills)
summary(ngbinmodel1)
xtable(ngbinmodel1)

#Drop with LRT test
drop1(ngbinmodel1,test="Chisq")
ngbinmodel2<-update(ngbinmodel1, . ~ . -POLIC)
drop1(ngbinmodel2,test="Chisq")
ngbinmodel3<-update(ngbinmodel2, . ~ . -D.WAT.COUR)
drop1(ngbinmodel3,test="Chisq")
ngbinmodel4 <- update(ngbinmodel3, . ~ . -WAT.RES)
drop1(ngbinmodel4,test="Chisq")
ngbinmodel5 <- update(ngbinmodel4, . ~ . -SHRUB)
drop1(ngbinmodel5,test="Chisq")
ngbinmodel6 <- update(ngbinmodel5, . ~ . -MONT.S)
drop1(ngbinmodel6,test="Chisq")
ngbinmodel7 <- update(ngbinmodel6, . ~ . -L.P.ROAD)
drop1(ngbinmodel7,test="Chisq")
drop1(ngbinmodel7,test="Rao")
#TOT.N ~ OPEN.L + D.PARK + L.WAT.C

#stepwise by AIC and BIC
step(ngbinmodel0, scope = list(upper=ngbinmodel1),direction="both",data=Roadkills) #AIC
#TOT.N ~ OPEN.L + D.PARK + L.WAT.C + L.P.ROAD
step(ngbinmodel0, scope = list(upper=ngbinmodel1),direction="both",k=log(n),data=Roadkills) #BIC
#TOT.N ~ OPEN.L + D.PARK + L.WAT.C

summary(ngbinmodel7)
xtable(ngbinmodel7)

#Interactions??
interaction<-glm.nb(TOT.N ~OPEN.L*D.PARK*L.WAT.C-OPEN.L : D.PARK : L.WAT.C ,data = Roadkills)
summary(interaction)
step(ngbinmodel0, scope = list(upper=interaction),direction="both",data=Roadkills) #BIC
interaction2<-glm.nb(TOT.N ~OPEN.L*D.PARK,data = Roadkills)

#Try all subsets
#require(MuMIn) #not included in class :(
#msnegbin<-dredge(ngbinmodel1,extra=alist("AIC","BIC"))
#msnegbin
#same results

####Negative binomial regression: diagnostic####

#Deviance
dev.negbin<-summary(ngbinmodel7)$deviance
df.negbin<-summary(ngbinmodel7)$df.residual
dev.negbin/df.negbin
1-pchisq(dev.negbin,df.negbin)
#no overdispersion

par(mfrow=c(2,2))
plot(ngbinmodel7)

#Residual plots
par(mfrow=c(1,3))
r.dev.ngb <- residuals(ngbinmodel7, type = "deviance")
#OPEN.L
plot(OPEN.L,r.dev.ngb,xlab="OPEN.L",ylab="Deviance residual",
     cex.lab=1.5,cex.axis=1.3)
loess.dev <- loess(r.dev.ngb~OPEN.L)
lo.pred <- predict(loess.dev, se=T)
orderOPEN.L <- order(OPEN.L)
lines(OPEN.L[orderOPEN.L ],lo.pred$fit[orderOPEN.L ],col="blue",lwd=3)
lines(OPEN.L[orderOPEN.L ],lo.pred$fit[orderOPEN.L ]+2*lo.pred$s[orderOPEN.L ], lty=2,col="red")
lines(OPEN.L[orderOPEN.L ],lo.pred$fit[orderOPEN.L ]-2*lo.pred$s[orderOPEN.L ], lty=2,col="red")

#D.PARK
plot(D.PARK,r.dev.ngb,xlab="D.PARK",ylab="Deviance residual",
     cex.lab=1.5,cex.axis=1.3)
loess.dev <- loess(r.dev.ngb~D.PARK)
lo.pred <- predict(loess.dev, se=T)
orderD.PARK <- order(D.PARK)
lines(D.PARK[orderD.PARK ],lo.pred$fit[orderD.PARK ],col="blue",lwd=3)
lines(D.PARK[orderD.PARK ],lo.pred$fit[orderD.PARK ]+2*lo.pred$s[orderD.PARK ], lty=2,col="red")
lines(D.PARK[orderD.PARK ],lo.pred$fit[orderD.PARK ]-2*lo.pred$s[orderD.PARK ], lty=2,col="red")

#L.WAT.C
#par(mfrow=c(1,1))
plot(L.WAT.C,r.dev.ngb,xlab="L.WAT.C",ylab="Deviance residual",
     cex.lab=1.5,cex.axis=1.3)
loess.dev <- loess(r.dev.ngb~L.WAT.C)
lo.pred <- predict(loess.dev, se=T)
orderL.WAT.C <- order(L.WAT.C)
lines(L.WAT.C[orderL.WAT.C ],lo.pred$fit[orderL.WAT.C ],col="blue",lwd=3)
lines(L.WAT.C[orderL.WAT.C ],lo.pred$fit[orderL.WAT.C ]+2*lo.pred$s[orderL.WAT.C ], lty=2,col="red")
lines(L.WAT.C[orderL.WAT.C ],lo.pred$fit[orderL.WAT.C ]-2*lo.pred$s[orderL.WAT.C ], lty=2,col="red")


#Influence plots
par(mfrow=c(2,2))
plot(hatvalues(ngbinmodel7),rstudent(ngbinmodel7))
plot(id,dffits(ngbinmodel7),type="l")
plot(id,covratio(ngbinmodel7),type="l")
plot(id,cooks.distance(ngbinmodel7),type="l")


####Negative binomial regression: Prediction####

Roadkills.p<-Roadkills
Roadkills.p$predy <- predict(ngbinmodel7,type="response")
par(mfrow=c(1,1))
plot(TOT.N ~ predy, data = Roadkills.p,pch = 16,xlab="Predicted response",ylab="Actual response")
predict.glm(ngbinmodel7,type="response",se.fit=T)



####Prediction Plot####
# NB regression
minage <- min(L.WAT.C,na.rm=TRUE)
maxage <- max(L.WAT.C,na.rm=TRUE)
grid <- seq(minage,maxage,0.1)
model.diet1.p <- data.frame(OPEN.L=mean(OPEN.L),
                            D.PARK=mean(D.PARK),L.WAT.C=grid)
p1 <- predict (ngbinmodel7,model.diet1.p,se=T,type = "response")
b <- p1$fit
LL <- p1$fit - 1.96 * p1$se.fit
UL <- p1$fit + 1.96 * p1$se.fit
a <- cbind(b,LL,UL)

par(mfrow=c(1,1))
par(pty="m")
par(mar = c(5,5,4,2)+0.1)
matplot(grid,a,lty=c(1,2,2),type="l",lwd=3,
        xlab="Length of water courses(km)",
        ylab="Predicted total number of amphibian fatalitis",
        cex.lab=1.5,cex.axis=1.3)
rug(Roadkills$L.WAT.C)

####Quasi-Poisson####

qpoismodel1 <- glm(TOT.N ~ .,family = quasipoisson(link="log"),data = Roadkills) 
summary(qpoismodel1)
quasi.full <- summary(qpoismodel1)
xtable(quasi.full$coefficients,digits = 4)
coeftest(qpoismodel1)
c <- as.data.frame(coeftest(qpoismodel1,vcov=sandwich))

xtable(c)
#Drop with LRT test
drop1(qpoismodel1,test="Chisq")
qpoismodel2<-update(qpoismodel1, . ~ . -D.WAT.COUR)
drop1(qpoismodel2,test="Chisq")
qpoismodel3<-update(qpoismodel2, . ~ . -POLIC)
drop1(qpoismodel3,test="Chisq")
qpoismodel4<-update(qpoismodel3, . ~ . -L.P.ROAD)
drop1(qpoismodel4,test="Chisq")
qpoismodel5<-update(qpoismodel4, . ~ . -WAT.RES)
drop1(qpoismodel5,test="Chisq")
qpoismodel6<-update(qpoismodel5, . ~ . -OPEN.L)
drop1(qpoismodel6,test="Chisq")
drop1(qpoismodel6,test="Rao")
qpoismodel7<-update(qpoismodel6, . ~ . -SHRUB)
drop1(qpoismodel7,test="Chisq")
drop1(qpoismodel7,test="Rao")

quasi.sum <- summary(qpoismodel7)
xtable(quasi.sum$coefficients,digits = 4)

dev.quasi<-quasi.sum$deviance
df.quasi<-quasi.sum$df.residual
dev.quasi/df.quasi
1-pchisq(dev.quasi,df.quasi)

####Influencial plot####
par(mfrow=c(2,2))
id <- 1:n
plot(hatvalues(qpoismodel7),rstudent(qpoismodel7))
plot(id,dffits(qpoismodel7),type="l")
plot(id,covratio(qpoismodel7),type="l")
plot(id,cooks.distance(qpoismodel7),type="l")

####quasi-poisson regression: Prediction####

Roadkills.p<-Roadkills
Roadkills.p$predy <- predict(qpoismodel7,type="response")
par(mfrow=c(1,1))
plot(TOT.N ~ predy, data = Roadkills.p,pch = 16,xlab="Predicted response",ylab="Actual response")
abline(a = 0,b = 1,col="red")
# predict.glm(qpoismodel7,type="response",se.fit=T)

####quasi-poisson regression: Prediction Plot (CI)####
#For quasi
minage <- min( L.WAT.C,na.rm=TRUE)
maxage <- max( L.WAT.C,na.rm=TRUE)
grid <- seq(minage,maxage,0.1)
model.qua <- data.frame(MONT.S=mean(MONT.S),
                        D.PARK=mean(D.PARK),L.WAT.C=grid)
p1.qua <- predict (qpoismodel7,model.qua,se=T,type = "response")
b.qua <- p1.qua$fit
LL.qua <- p1.qua$fit - 1.96 * p1.qua$se.fit
UL.qua <- p1.qua$fit + 1.96 * p1.qua$se.fit
a.qua <- cbind(b.qua,LL.qua,UL.qua)

par(mfrow=c(1,1))
par(pty="m")
par(mar = c(5,5,4,2)+0.1)
matplot(grid,a.qua,lty=c(1,2,2),type="l",lwd=3,
        xlab="Length of water courses(km)",
        ylab="Predicted total number of amphibian fatalitis",
        cex.lab=1.5,cex.axis=1.3)
# title("Quasi-Poisson regression")

rug(Roadkills$L.WAT.C)




#####Comparison of quasi model with the NB model##### 
poismodel3<-glm(TOT.N~MONT.S+D.PARK+L.WAT.C,family=poisson(link="log"),data=Roadkills)
summary(poismodel3)
par(mfrow=c(1,1))
plot(qpoismodel1,which = 5,main = "(Quasi-Poisson)")
plot(qpoismodel7,which = 5,main = "(Quasi-Poisson)")
plot(qpoismodel7,which = 5)
#plot(poismodel3,which = 4,main = "(Poisson)")
par(mfrow=c(1,1))
plot(poismodel2,which=5)
plot(ngbinmodel7,which=5)


####Prediction NB vs quasi over Length of water courses####
#For quasi
minage <- min( L.WAT.C,na.rm=TRUE)
maxage <- max( L.WAT.C,na.rm=TRUE)
grid <- seq(minage,maxage,0.1)
model.qua <- data.frame(MONT.S=mean(MONT.S),
                        D.PARK=mean(D.PARK),L.WAT.C=grid)
p1.qua <- predict (qpoismodel7,model.qua,se=T,type="response")
b.qua <- p1.qua$fit
LL.qua <- p1.qua$fit - 1.96 * p1.qua$se.fit
UL.qua <- p1.qua$fit + 1.96 * p1.qua$se.fit
a.qua <- cbind(b.qua,LL.qua,UL.qua)
#for NB part
model.diet1.nb <- data.frame(OPEN.L=mean(OPEN.L),
                             D.PARK=mean(D.PARK),L.WAT.C=grid)
p1.nb <- predict (ngbinmodel7,model.diet1.nb,se=T,type="response")
b.nb <- p1.nb$fit
LL.nb <- p1.nb$fit - 1.96 * p1$se.fit
UL.nb <- p1.nb$fit + 1.96 * p1$se.fit
c.nb <- cbind(b.nb,LL.nb,UL.nb)
#plot both in one fig 
par(mfrow=c(1,1))
par(pty="m")
par(mar = c(5,5,4,2)+0.1)
matplot(grid,a.qua,lty=c(1,2,2),type="l",col=c(2,2,2),lwd=3,
        xlab="Length of water courses(km)",
        ylab="Predicted total number of amphibian fatalities",
        # ylab="Log of predicted amphibian fatalities",
        cex.lab=1.5,cex.axis=1.3 )
polygon(c(rev(grid), grid), c(rev(c.nb[ ,3]), c.nb[ ,2]), col = 'grey80', border = NA,density = c(50, 50) )
lines(grid, c.nb[,1], lty=1, lwd=3,col=1)
lines(grid,c.nb[,2], lty=4,lwd=3, col=1)
lines(grid,c.nb[,3], lty=4,lwd=3, col=1)
rug(Roadkills$L.WAT.C)



####Bayesian Poisson regression####

poismodel.bayes2 <- MCMCpoisson(TOT.N ~.-D.WAT.COUR,data = Roadkills, burnin = 10000, mcmc=1000000)
# poismodel.bayes2 <- MCMCpoisson(TOT.N ~ OPEN.L+MONT.S+POLIC+D.PARK+SHRUB+WAT.RES+L.WAT.C+L.P.ROAD,data = Roadkills, burnin = 10000, mcmc=1000000)
# quasipo.bayes <- MCMCpoisson(TOT.N ~ MONT.S+ D.PARK + L.WAT.C,
#                              family=poisson, data = Roadkills,burnin = 1000, mcmc = 10000)
a <- summary(poismodel.bayes2)
xtable(a$statistics,digits = 4)
# xtable(a$quantiles)
summary(poismodel.bayes)
#trace plot
plot(poismodel.bayes2)
# 100(1?????)% equal tail CI [a,b]:
xtable(a$quantiles, digits = 4)


detach(Roadkills)
