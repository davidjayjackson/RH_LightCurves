########################################
# Grant Foster functions to be available for R
########################################
#############
#############
# anova1
# 1-way ANOVA
#############
#############
anova1 = function(x,fac,sort=F,progress=F){
    n = length(x)
    q = as.factor(fac)
    nlev = length(levels(q))
    #########################
    # compute counts and sums
    #########################
    xnum = NULL
    xsum = NULL
    xsum2 = NULL
    xlev = NULL
    for (j in 1:nlev){
        qlev = levels(q)[j]
        xx = x[q==qlev]
        xx = xx[is.finite(xx)]
        if (length(xx)>0){
            xnum = c(xnum,length(xx))
            xsum = c(xsum,sum(xx))
            xsum2 = c(xsum2,sum(xx^2))
            xlev = c(xlev,qlev)
        }
        if (progress){
            plot(0,0)
            title(paste(j,"/",nlev))
            x11()
        }
    }
    nlev = length(xlev)
    n = sum(xnum)
    xave = xsum/xnum                # average by level
    xtotave = sum(xsum)/sum(xnum)   # grand average
    xvary = xsum2 - xnum*xave^2     # variation
    xtotvary = sum(xsum2)-n*xtotave^2   # total variation
    xvar = xvary/(xnum-1)           # variance by level
    xsig = sqrt(xvar)               # std.dev. by level
    xtotvar = xtotvary/(n-1)        # total variance
    ###############################
    # standard error and confidence
    # limits for averages by level
    ###############################
    se = sqrt(xvar/xnum)
    tval = (xave-xtotave)/se
    for (j in 1:nlev){
        if (xnum[j] < 2){
            se[j] = 0
            tval[j] = 0
        }
    }
    lower.ci = xave-2*se
    upper.ci = xave+2*se
    #########################
    # plot averages by level
    # with 2-sigma error bars
    #########################
    t = 1:nlev
    ylo = min(lower.ci); yhi = max(upper.ci)
    plot(xave,ylim=c(ylo,yhi),
        xlab="",ylab="Average",cex.axis=1.5,cex.lab=1.5)
    arrows(t,lower.ci,t,upper.ci,length=.05,angle=90,code=3)
    abline(h=xtotave,lty="dashed")
    ############################
    # data for individual levels
    ############################
    zframe = data.frame(level=xlev,count=xnum,
        ave=xave,std.dev=xsig,se=se,t.val=round(tval,2))
    if (sort){
        zhold = zframe
        ord = order(-zhold$count)
        for (j in 1:length(ord)){
            zframe[j,] = zhold[ord[j],]
        }
    }
    ##################
    # F-test for ANOVA
    ##################
    Vwithin = sum(xvary)
    Vtotal = xtotvary
    Vbetween = Vtotal-Vwithin
    df.between = nlev - 1
    df.within = n-nlev
    Fstat = Vbetween*df.within/Vwithin/df.between
    pval = 1-pf(Fstat,df.between,df.within)
    Ftest = data.frame(Fstat,df.between,df.within,p=pval)
    Ftest = round(Ftest,6)
    zout = list(bylevel=zframe,grandave=xtotave,Ftest=Ftest)
    zout
}
#######################################
#                aovper
# Period search by Analysis of Variance
#######################################
aovper = function(t,x,
        lowfreq=0,hifreq=0,
        bins=4,resmag=1,
        plot=T,outfile=NULL){
    ndata = length(t)
    tspan = t[ndata]-t[1]
    tspan = tspan*ndata/(ndata-1)
    if (lowfreq <= 0){
        lowfreq = 1/tspan
    }
    if (hifreq <= 0){
        hifreq = (ndata-1)/tspan/2
    }
    if (resmag <=0){
        resmag = 1
    }
    nustep = .25/tspan/resmag
    xbar = mean(x)
    x2sum = sum(x^2)
    szero = x2sum - ndata*xbar^2
    freq = NULL
    per = NULL
    power = NULL
    nu = lowfreq
    ###########################
    # scan the test frequencies
    ###########################
    while (nu <= hifreq){
        freq = c(freq,nu)
        per = c(per,1/nu)
        ####################
        # phase-bin the data
        ####################
        phase = nu*t
        phase = phase - floor(phase)
        nbin = floor(bins*phase) + 1
        bincount = rep(0,bins)
        binave = bincount
        for (j in 1:bins){
            bindat = x[nbin==j]
            bincount[j] = length(bindat)
            binave[j] = sum(bindat)
            if (bincount[j]>0){
                binave[j] = binave[j]/bincount[j]
            }
        }
        s1 = sum(bincount*(binave-xbar)^2)
        s2 = szero - s1
        pow = s1/s2
        power = c(power,pow)
        nu = nu + nustep
    }
    power = power*(ndata-bins)/(bins-1)
    aovp = data.frame(freq,power,per)
    aovp = round(aovp,6)
    if (plot){
        plot(freq,power,type="l")
    }
    if (length(outfile)==1){
        write.table(aovp,file=outfile,row.names=F)
    }
    aovp
}
#############################################
#############################################
#                   dcdft
# date-compensated discrete Fourier transform
#############################################
#############################################
dcdft = function(t,x,lowfreq=0,hifreq=0,resmag=1,plot=T,outfile=NULL){
    ndata = length(t)
    xbar = mean(x)
    xvar = var(x)*(ndata-1)/ndata
    tspan = t[ndata]-t[1]
    tspan = tspan*ndata/(ndata-1)
    if (lowfreq <= 0){
        lowfreq = 1/tspan
    }
    if (hifreq <= 0){
        hifreq = (ndata-1)/tspan/2
    }
    if (resmag <=0){
        resmag = 1
    }
    nustep = .25/tspan/resmag
    freq = NULL
    per = NULL
    power = NULL
    amp = NULL
    nu = lowfreq
    while (nu <= hifreq){
        freq = c(freq,nu)
        per = c(per,1/nu)
        c1 = cos(2*pi*nu*t)
        s1 = sin(2*pi*nu*t)
        xfit = lm(x~c1+s1)
        pow = xfit$coef[2]^2 + xfit$coef[3]^2
        amp = c(amp,sqrt(pow))
        pow = var(xfit$fit)
        power = c(power,pow)
        nu = nu + nustep
    }
    power = power*ndata/xvar/2
    dcd = data.frame(freq,power,amp,per)
    dcd = round(dcd,8)
    if (plot){
        plot(freq,power,type="l",
        xlim=c(max(freq),min(freq)),ylim=c(max(power),min(power)),
        xlab="Frequency",ylab="Power",main="mag",
        cex.axis=1.5,cex.lab=1.5)
        x11()
        plot(freq,power,type="l")
        x11()
    }
    if (length(outfile)==1){
        write.table(dcd,file=outfile,row.names=F)
    }
    dcd
}
###########################################
# find starting index of a particular value
###########################################
findstart = function(X,start){
    n = length(X)
    nstart = 0
    for (j in 1:n){
        if (X[j] >= start) break
    }
    j
}
########################################
#               foldit
# fold a light curve with a trial period
########################################
foldit = function(t,x,period,epoch=0,plot=T,bin=.05){
    n = length(t)
    phase = (t-epoch)/period
    phase = phase - floor(phase)
    ord = order(phase)
    newt = numeric(n)
    newx = numeric(n)
    newfaz = numeric(n)
    for (j in 1:n){
        newfaz[j] = phase[ord[j]]
        newx[j] = x[ord[j]]
        newt[j] = t[ord[j]]
    }
    newfaz = c(newfaz-1,newfaz)
    newx = c(newx,newx)
    newt = c(newt,newt)
    folded = data.frame(phase=newfaz,x=newx,t=newt)
    ############################################
    # plot 
    ############################################
    ylo = min(newx)
    yhi = max(newx)
    if (plot){
        q = timave(newfaz,newx,bin,plot=F)
        plot(newfaz,newx,pch=".",cex=3,ylim=c(ylo,yhi),
            xlab="Phase (cycles)",ylab="Values",main="magnitude",
            cex.axis=1.5,cex.lab=1.5)
        #x11()
        lines(q$t,q$ave,type="o",pch=20,cex=2,lwd=2,col="red")
        arrows(q$t,q$lolim,q$t,q$uplim,
            length=.04,angle=90,code=3,col="red",lwd=2)
        x11()
    }
    folded = round(folded,6)
    folded
}
######################
# fourfit:
# Fit a Fourier series
######################
fourfit = function(t,x,p){
    w = 2*pi/p    # frequencies from periods
    ndata = length(t)
    nfreq = length(w)
    cs = matrix(nrow=ndata,ncol=2*nfreq)
    for (j in 1:nfreq){
        cs[,(2*j-1)] = cos(w[j]*t)
        cs[,(2*j)] = sin(w[j]*t)
    }
    if (nfreq==1){
        xfit = lm(x~cs[,1]+cs[,2])
    }
    if (nfreq==2){
        xfit = lm(x~cs[,1]+cs[,2]+cs[,3]+cs[,4])
    }
    if (nfreq==3){
        xfit = lm(x~cs[,1]+cs[,2]+cs[,3]+cs[,4]+cs[,5]+cs[,6])
    }
    if (nfreq==4){
        xfit = lm(x~cs[,1]+cs[,2]+cs[,3]+cs[,4]+cs[,5]+cs[,6]+cs[,7]+cs[,8])
    }
    if (nfreq==5){
        xfit = lm(x~cs[,1]+cs[,2]+cs[,3]+cs[,4]+cs[,5]+cs[,6]+cs[,7]+cs[,8]
        +cs[,9]+cs[,10])
    }
    if (nfreq==6){
        xfit = lm(x~cs[,1]+cs[,2]+cs[,3]+cs[,4]+cs[,5]+cs[,6]+cs[,7]+cs[,8]
        +cs[,9]+cs[,10]+cs[,11]+cs[,12])
    }
    if (nfreq>6){
        xfit = lm(x~cs[,1]+cs[,2]+cs[,3]+cs[,4]+cs[,5]+cs[,6]+cs[,7]+cs[,8]
        +cs[,9]+cs[,10]+cs[,11]+cs[,12]+cs[,13]+cs[,14])
    }
    xout = data.frame(time=t,fit=xfit$fitted.values,res=xfit$resid,data=x)
    xout = round(xout,6)
    xout
}
##########################
# compute a moving average
##########################
movave = function(X,L=12){
    xave = numeric(length(X)+1-L)
    for (j in 1:(length(X)+1-L)){
        sum = 0
        for (k in j:(j+L-1)){
            sum = sum + X[k]
        }
        xave[j] = sum/L
    }
    xave
}
##############################################
##############################################
# find the main peak (only) in DCDFT or AoVper
##############################################
##############################################
peak1 = function(dcd,t=NULL,x=NULL){
    peakpow = max(dcd$pow)          # max power level
    nmax = order(-dcd$pow)[1]       # index of maximum
    peakline = dcd[nmax,]
    peakfreq = dcd$freq[nmax]
    peakper = dcd$per[nmax]
    peakamp = dcd$amp[nmax]
    ############################
    # find FWHM lo and up limits
    ############################
    fwhmlo = peakfreq
    fwhmup = peakfreq
    for (j in 1:100){
        ndx = nmax-j
        if (ndx < 1 || ndx > length(dcd$pow)) break
        ndone = F
        if (dcd$pow[ndx] >= peakpow/2){
            fwhmlo = dcd$freq[ndx]
        } else{
            ndone = T
        }
        if (ndone) break
    }
    for (j in 1:100){
        ndx = nmax+j
        if (ndx < 1 || ndx > length(dcd$pow)) break
        ndone = F
        if (dcd$pow[ndx] >= peakpow/2){
            fwhmup = dcd$freq[ndx]
        } else{
            ndone = T
        }
        if (ndone) break
    }
    ################################
    # if given t,x then find period,
    # amplitude std.err based on
    # theoretical formulae
    ################################
    freq.se = NA
    amp.se = NA
    if (length(t) > 10){
        if (length(x) == length(t)){
            if (is.finite(peakamp)){
                phi = 2*pi*peakfreq*t
                c1 = cos(phi)
                s1 = sin(phi)
                xfit = lm(x~c1+s1)
                svar = var(xfit$res)
                numdat = length(t)
                totalt = t[numdat]-t[1]
                freq.se = 6*svar/numdat/(pi*peakamp*totalt)^2
                freq.se = sqrt(freq.se)
                amp.se = 2*svar/numdat
                amp.se = sqrt(amp.se)
            }
        }
    }
    ################
    # return results
    ################
    data.frame(peakline,freq.se,amp.se,
        fwhm.lo=fwhmlo,fwhm.up=fwhmup,per.lo=1/fwhmup,per.up=1/fwhmlo)
}
#####################################
# find the peaks in a DCDFT or AOVPER
#####################################
peaks = function(dcd,lofre=0,hifre=0,maxpeak=0,plot=T){
    nfreq = length(dcd$freq)
    if (hifre <= 0) { hifre = max(dcd$freq) }
    peaked = NULL
    peakpow = NULL
    for (j in 2:(nfreq-1)){
        if (dcd$freq[j] >= lofre && dcd$freq[j] <= hifre){
            if (dcd$power[j] > dcd$power[j-1] && dcd$power[j] > dcd$power[j+1]){
                peaked = rbind(peaked,dcd[j,])
                peakpow = c(peakpow,dcd$power[j])
            }
        }
    }
    npeak = length(peaked[,1])
    nstats = dim(peaked)[2]
    dummy = rep(0,npeak)
    peakout = cbind(peaked,quantile3=dummy)
    q = order(-peakpow)
    for (j in 1:npeak){
        for (k in 1:nstats){
            peakout[j,k] = peaked[q[j],k]
        }
        zunif = 1-j/(npeak+1)
        q3 = qchisq(zunif,3)
        peakout[j,nstats+1] = q3
    }
    if (maxpeak > 0 && npeak > maxpeak){
        peakout = peakout[1:maxpeak,]
    }
    if (plot){
        plot(peakout$quantile3,2*peakout$pow,pch=20,
        xlab="theoretical chi-square",ylab="observed chi-square",
        cex.axis=1.5,cex.lab=1.5)
        x11()
    }
    peakout
}
########################
########################
# time-based averaging
########################
########################
timave = function(t,x,t2ave=-1,t.off=0,n.min=2,wide=2,plot=T,lines=T,tit=NULL,big=F){
    N = length(x)
    if (t2ave<=0){
        t2ave = 1
    }
    n = 0
    tsum = 0
    sum = 0
    sum2 = 0
    time = NULL
    ave = NULL
    se = NULL
    num = NULL
    lasttblock = floor((t[1]-t.off)/t2ave)
    for (j in 1:N){
        tblock = floor((t[j]-t.off)/t2ave)
        if (tblock > lasttblock){
            if (n >= n.min){
                q = tsum/n
                time = c(time,q)
                q = sum/n
                ave = c(ave,q)
                if (n>1){
                    q = (sum2 - n*q^2)/(n-1)
                    if (q > 0){
                        q = sqrt(q/n)
                    } else{
                        q = 0
                    }
                }
                else{
                    q = 0
                }
                se = c(se,q)
                num = c(num,n)
                #########################
                # reset sums and counters
                #########################
                n = 0
                tsum = 0
                sum = 0
                sum2 = 0
            }
            lasttblock = tblock
        }
        if (is.finite(x[j])){
            n = n + 1
            tsum = tsum + t[j]
            sum = sum + x[j]
            sum2 = sum2 + x[j]^2
        }
    }
    ####################
    # do the final point
    ####################
    if (n>0){
        q = tsum/n
        time = c(time,q)
        q = sum/n
        ave = c(ave,q)
        if (n>1){
                    q = (sum2 - n*q^2)/(n-1)
                    if (q > 0){
                        q = sqrt(q/n)
                    } else{
                        q = 0
                    }
        }
        else{
            q = 0
        }
        se = c(se,q)
        num = c(num,n)
    }
    mu = mean(ave)
    for (j in 1:length(se)){
        if (se[j]<=0){
            q = abs(ave[j] - mu)+wide*sd(ave)
            q = q/wide
            se[j] = max(wide*sd(ave),q)
        }
    }
    lolim = ave - wide*se
    uplim = ave + wide*se
    dev = (ave - mu)/se
    norm = (ave - mu)/sd(ave)
    xout = data.frame(time,ave,se,num,lolim,uplim,dev,norm)
    xout = round(xout,6)
    if (plot){
        loplot=min(lolim)
        hiplot=max(uplim)
        if (lines){
            plot(time,ave,type="o",pch=".",cex=5,ylim=c(loplot,hiplot),xlab="",ylab="",main=tit)
           
        }
        else{
            plot(time,ave,pch=".",cex=5,ylim=c(loplot,hiplot),xlab="",ylab="",main=tit)
          
        }
        arrows(time,lolim,time,uplim,length=.04,angle=90,code=3)
        abline(mean(ave),0)
        abline(mean(ave)+wide*sd(ave),0,col="red",lty="dashed")
        abline(mean(ave)-wide*sd(ave),0,col="red",lty="dashed")
        if (big){
            bigt = time[abs(dev)>1.25*wide]
            biga = ave[abs(dev)>1.25*wide]
            points(bigt,biga,col="red",cex=2,lwd=2)
            bigt = time[abs(norm)>wide]
            biga = ave[abs(norm)>wide]
            points(bigt,biga,col="blue",cex=3,lwd=2)
        }
    }
    xout
}
