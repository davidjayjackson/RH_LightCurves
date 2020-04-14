# Initialization


setwd("c:/users/howe/desktop/")
WD <- getwd()
source('Rprogs.R')

infile <- paste(WD, "RS_Per_AAVSO.csv ", sep="/")
alvestad_sw <- data.frame(read.csv(infile, header=TRUE))
outfile <- paste(WD, "RS_Per.RData", sep="/")
save(alvestad_sw, file=outfile, ascii=FALSE)

RS_Per <- ts(alvestad_sw)
#RS_Per = read.table("RS_Per.10d")

summary(RS_Per)
sd(RS_Per)


t = RS_Per[,1]
x = RS_Per[,3]
###########################
# just the data since day 1
###########################
n1 = findstart(t,1)
n2 = length(t)
t = t[n1:n2]
x = x[n1:n2]
plot(t,x)
x11()
plot(t,x,ylim=c(120,2))
x11()
plot(t,x,ylim=c(120,2),pch=20)
x11()
plot(t,x,pch=20,
    ylim=c(max(x),min(x)),
    xlab="Time",ylab="Sunspots",main="RS_Per",
    cex.lab=1.5,cex.axis=1.5,cex.main=2.5)
x11()
plot(t,x,type="o",pch=20,
    ylim=c(max(x),min(x)),
    xlab="Time",ylab="Sunspots",main="RS_Per",
    cex.lab=1.5,cex.axis=1.5,cex.main=2.5)
x11()
plot(t,x,type="l",
    ylim=c(max(x),min(x)),
    xlab="Time",ylab="Sunspots",main="RS_Per",
    cex.lab=1.5,cex.axis=1.5,cex.main=2.5)
x11()
#################
# AoV periodogram
#################
an1 = anova1(x,1.0)

ap = aovper(t,x)
x11()
########################
# identify all the peaks
########################
pp = peaks(ap)
#######################
# display the 1st three
#######################
pp[1:3,]
##################################
# identify just the strongest peak
##################################
pk1 = peak1(ap); pk1
period = pk1$per; period
###################
# DCDFT periodogram
###################
#dc = dcdft(t,x,.03,.05)
dc = dcdft(t,x,.03,.05)

###################################
# identify all peaks, display top 3
###################################
pp = peaks(dc)
pp[1:3,]
################
# strongest peak
################
peak1(dc)
#x11()
#####################################
# strongest peak with standard errors
# for frequency and amplitude
#####################################
peak1(dc,t,x)
#x11()
############################
# DCDFT closeup of main peak
############################
dc2 = dcdft(t,x,.03,.05)
#############################
# DCDFT closeup of main peak
# with 5x enhanced resolution
#############################
dc2 = dcdft(t,x,.03,.05,resmag=2)
pk1 = peak1(dc2,t,x); pk1
period = pk1$per
xfit = fourfit(t,x,period)
plot(t,x,type="o",pch=20,
    xlim=c(min(t),max(t)),ylim=c(min(x),max(x)),
    xlab="Days",ylab="Sunspots",main="RS_Per",
    cex.lab=1.5,cex.axis=1.5,cex.main=2.5)
lines(t,xfit$fit,col="red",lwd=3)
abline(v=0)
grid()
x11()

#############################
# folded light curve for RS_Per
#############################
fold = foldit(t,x,period)
fold = foldit(t,x,period,epoch=56400)


