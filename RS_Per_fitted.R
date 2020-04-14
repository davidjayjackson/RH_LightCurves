
setwd("C:/Users/Howe/Desktop/")
WD <- getwd()
source('Rprogs.R')
###   Specify the name and address of the data file,   ####
###   read it in, save it as RData

infile <- paste(WD, "RS_Per_aavso.csv", sep="/")
rstar <- data.frame(read.csv(infile, header=TRUE))
outfile <- paste(WD, "rstar.RData", sep="/")
save(rstar, file=outfile, ascii=FALSE)
summary(rstar)
###   Set as time series   #################################

#load("rstar.RData")
CLDdata <- ts(rstar)
t = CLDdata[,1]
x = CLDdata[,3]
plot(t,x,pch=".",cex=3,col="blue",lwd=3,
    xlab="JD ",
    ylab=" RS_Per",
    main="RS_Per",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()
#########################################
# simplest model: linear function of time
#########################################
xfit = lm(x ~ t)
summary(xfit)
##########################
# superimpose linear model
##########################
#lines(t,xfit$fitted.values,col="black",lwd=3)
plot(t,xfit$fitted.values,col="black",lwd=3,
    xlab="JD ",
    ylab=" RS_Per fitted values",
    main=" linear function",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()

################
# plot residuals
################
plot(t,xfit$residuals,pch=".",cex=3,col="blue",lwd=3,
    xlab="JD ",
    ylab="Residuals ",
    main=" linear model",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()
##############################
# quadratic model of residuals
##############################
tsquablack = t^2
xfit2 = lm(xfit$res ~ t+tsquablack)
summary(xfit2)
###########################
# superimpose model on data
###########################
#lines(t,xfit2$fit,col="black",lwd=3)
plot(t,xfit2$fit,col="black",lwd=3,
    xlab="JD ",
    ylab="Residuals fitted values (Quadratic) ",
    main=" linear model",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()


#############################
# quadratic model of raw data
#############################
xfit3 = lm(x ~ t+tsquablack)
summary(xfit3)
#################################
# plot data and superimpose model
#################################
plot(t,x,pch=".",cex=3,
    xlab="JD ",
    ylab=" ",
    main=" Quadratic",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()
#lines(t,xfit3$fit,col="black",lwd=3)
plot(t,xfit3$fit,col="black",lwd=3,
    xlab="JD ",
    ylab=" RS_Per fitted values",
    main=" Quadratic",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()
################
# plot residuals
################
plot(t,xfit3$resid,pch=".",cex=3,col="blue",lwd=3,
    xlab="JD ",
    ylab="Residuals ",
    main=" Quadratic",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()


############################
# quadratic + sinusoid model
############################
c1 = cos(2*pi*t)
s1 = sin(2*pi*t)
xfit4 = lm(x ~ t+tsquablack+c1+s1)
plot(t,x,pch=".",cex=3,
    xlab="JD ",
    ylab=" ",
    main=" sinusoid model",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()
#lines(t,xfit4$fit,col="black",lwd=2)
plot(t,xfit4$fit,col="black",lwd=2,
    xlab="JD ",
    ylab=" RS_Per fitted values",
    main=" sinusoid model",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()
##################################
# quadratic + sinusoid w/harmonics
##################################
c2 = cos(4*pi*t)
s2 = sin(4*pi*t)
c3 = cos(6*pi*t)
s3 = sin(6*pi*t)
xfit5 = lm(x ~ t+tsquablack+c1+s1+c2+s2+c3+s3)
plot(t,x,pch=".",cex=3,col="black",lwd=2,
    xlab="JD ",
    ylab=" ",
    main=" harmonics",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()
#lines(t,xfit5$fit,col="black",lwd=2)
plot(t,xfit5$fit,col="black",lwd=2,
    xlab="JD ",
    ylab=" RS_Per fitted values",
    main=" harmonics",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)
x11()
###########
# residuals
###########
plot(t,xfit5$resid,type="o",pch=".",cex=3,col="black",lwd=2,
    xlab="JD ",
    ylab="Residuals ",
    main=" harmonics",
    cex.axis=1.5,cex.lab=1.5,cex.main=2)




