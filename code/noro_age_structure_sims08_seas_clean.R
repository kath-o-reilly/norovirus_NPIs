# use the code in King and Wearing to create SIR model with age structure, and use polymod data
# adapt for use applied to norovirus
# we have model parameters...
# at specific time points change the contact matrix

library(deSolve)
library(ggplot2)   
library(reshape2)
library(janitor)
library(qs) # comix data
library(RColorBrewer)

rm(list = ls())
setwd("~/Dropbox/Dropbox/Norovirus/Rcode/DeterministicModel")

# model attributes;
# - multiple age groups (30 age classes, first 20 years each have their own, further years of 5 yrs each)
# - bitrths to equal deaths
# - noro specific
births <- 1250

#ages <- c(c(1,5,10),seq(20,80,by=10)) # upper end of age classes (for Amar data)
ages <- c(4,14,24,34,44,54,64,80) # upper end of age classes (for Amar data)
da <- diff(c(0,ages))
length(ages)
# need to check age dist (in the absence of infection)
yinit <- c( S=c(rep(11500,7),21000),
            E=c(rep(0,length(ages))),
            I=c(rep(0,length(ages))),
            A=c(rep(0,length(ages))),
            R=c(rep(0,length(ages))),
            G=c(rep(0,length(ages))),
            cc=c(rep(0,length(ages))),
            aa=c(rep(0,length(ages)))
)
sum(yinit)
length(yinit)/8

sindex <- 1:length(ages) 
eindex <- length(ages)+(1:length(ages)) # 31:60 # 
iindex <- (length(ages)*2)+(1:length(ages))  # 61:90 
aindex <- (length(ages)*3)+(1:length(ages))  #91:120 
rindex <- (length(ages)*4)+(1:length(ages))  #121:150
gindex <- (length(ages)*5)+(1:length(ages))  #151:150
juvies <- 1
adults <- 2:length(ages)

aging <- diag(-1/da) 
aging[row(aging)-col(aging)==1] <- 1/head(da,-1)

#use polymod data
age.categories <- as.factor(ages)
# we want a symmetric contact matrix
library(socialmixr)
socm <- contact_matrix(polymod, countries = "United Kingdom", 
                       age.limits = c(0,as.numeric(as.character(age.categories[1:(length(ages)-1)]))), 
                       symmetric = TRUE)
dim(socm$matrix)   # units = "mean number of contacts"
poly <- data.frame(expand.grid(contactee=age.categories,contact=age.categories))
poly$vals <- as.vector(socm$matrix)
tapply(poly$vals,poly$contactee,sum)

# add in 'poly80' for 2021 onwards...
poly$vals80 <- poly$vals
# reduce adults by 80%
poly$adults <- 1
poly$adults[(poly$contact=="4"|poly$contact=="14") & (poly$contactee=="4"|poly$contactee=="14")] <- 0
poly$vals80[poly$adults==1] <- poly$vals80[poly$adults==1]*0.8

round(mean(tapply(poly$vals80,poly$contactee,sum)),2)

# Can treat the last X years of sol_old as years up to Dec 2019.
# Matrices have been adjusted to fit the age ranges here
# Have different scenarios and stitch the solutions together 
# (and follow on initial conditions)

comix <- read.csv("~/Dropbox/Dropbox/Norovirus/Data/comix_adjusted_08May2021.csv") %>% clean_names()
comix.period <- names(table(comix$period))
head(comix)
names(comix) <- c("contactee","contact","period","vals")  # so can be used in the function

comix$vals_bs <- 0

# 1
tmp <- qread("~/Dropbox/Dropbox/Norovirus/Data/share_matrices/bs1000_ngrps8_cap50_nwksx_sr1_Lockdown 1_scms.qs")
tmp1 <- apply(tmp,1,median)
comix[comix$period=="1. Lockdown 1",]$vals_bs <- tmp1
# 2
tmp <- qread("~/Dropbox/Dropbox/Norovirus/Data/share_matrices/bs1000_ngrps8_cap50_nwks8_sr11_Lockdown 1 easing_scms.qs")
tmp1 <- apply(tmp,1,median)
comix[comix$period=="2. Lockdown 1 easing",]$vals_bs <- tmp1
# 3
tmp <- qread("~/Dropbox/Dropbox/Norovirus/Data/share_matrices/bs1000_ngrps8_cap50_nwks5_sr19_Reduced restrictions_scms.qs")
tmp1 <- apply(tmp,1,median)
comix[comix$period=="3. Relaxed restrictions",]$vals_bs <- tmp1
# 4
tmp <- qread("~/Dropbox/Dropbox/Norovirus/Data/share_matrices/bs1000_ngrps8_cap50_nwks8_sr24_Schools open_scms.qs")
tmp1 <- apply(tmp,1,median)
comix[comix$period=="4. School reopening",]$vals_bs <- tmp1
# 5
tmp <- qread("~/Dropbox/Dropbox/Norovirus/Data/share_matrices/bs1000_ngrps8_cap50_nwks5_sr33_Lockdown 2_scms.qs")
tmp1 <- apply(tmp,1,median)
comix[comix$period=="5. Lockdown 2",]$vals_bs <- tmp1
# 6
tmp <- qread("~/Dropbox/Dropbox/Norovirus/Data/share_matrices/bs1000_ngrps8_cap50_nwks2_sr37_Lockdown 2 easing_scms.qs")
tmp1 <- apply(tmp,1,median)
comix[comix$period=="6. Lockdown 2 easing",]$vals_bs <- tmp1
# 7
tmp <- qread("~/Dropbox/Dropbox/Norovirus/Data/share_matrices/bs1000_ngrps8_cap50_nwks3_sr39_Christmas_scms.qs")
tmp1 <- apply(tmp,1,median)
comix[comix$period=="7. Christmas",]$vals_bs <- tmp1
# 8
tmp <- qread("~/Dropbox/Dropbox/Norovirus/Data/share_matrices/bs1000_ngrps8_cap50_nwks10_sr41_Lockdown 3_scms.qs")
tmp1 <- apply(tmp,1,median)
comix[comix$period=="8. Lockdown 3",]$vals_bs <- tmp1
# 9
tmp <- qread("~/Dropbox/Dropbox/Norovirus/Data/share_matrices/bs1000_ngrps8_cap50_nwks4_sr50_Lockdown 3 (with schools open)_scms.qs")
tmp1 <- apply(tmp,1,median)
comix[comix$period=="9. Lockdown 3 + schools",]$vals_bs <- tmp1

# check
#ggplot(comix[comix$period=="1. Lockdown 1",],aes(x=as.factor(contactee),y=as.factor(contact),fill=vals)) + geom_tile()
#ggplot(comix[comix$period=="4. School reopening",],aes(x=as.factor(contactee),y=as.factor(contact),fill=vals)) + geom_tile()
#max(comix[comix$period=="4. School reopening",]$vals)
#ggplot(poly,aes(x=as.factor(contactee),y=as.factor(contact),fill=vals)) + geom_tile()
#max(poly$vals)

comix2 <- comix
# check names
names(comix2) <- c("contactee","contact","period","vals_old","vals")

# Specify the dates ...
# see Chris' paper and 
# here https://www.instituteforgovernment.org.uk/sites/default/files/timeline-lockdown-web.pdf
# 1. Lockdown 1 = 23/03/2020
# 2. Lockdown 1 easing = 04/06/2020
# 3. Relaxed restrictions = 30/07/2020 ??
# 4. School reopening = 04/09/2020 ??
# 5. Lockdown 2 = 05/11/2020
# 6. Lockdown 2 easing = 03/12/2020
# 7. Christmas = 20/12/2020
# 8. Lockdown 3 = 05/01/2021
# 9. Lockdown 3 + schools = 09/03/2021
# 10. School reopening & relaxed restrictions (matrix 4) = 17/05/2021
# 11. Back to polymod = 17/07/2021 
# 12. End simulation = 31/12/2023 (with seasonality so don't need summer holidays)

# load in prevalence data from Harris
data <- read.csv("~/Dropbox/Dropbox/Norovirus/Data/harris_noro_incidence_withCI.csv") %>% clean_names()
head(data)

# fixed and estimated parameters...
mu_s <- 1/1  # rate of pre-inf to infectious
mu_a <- 1/2  # rate of infetious to asymptomatic
rho <- 1/15  # rate of asymtomatic to recovery

ccindex <- (length(ages)*6)+(1:length(ages))
aaindex <- (length(ages)*7)+(1:length(ages))

source("noro_functions_mcmc_clean.r")

parms <- list()
parms$mu_s <- mu_s
parms$mu_a <- mu_a
parms$rho <- 1/20#rho
parms$theta <- 5.1*365 #3181.696 #(~7 yrs)
parms$alpha <- 1 #0.5#alpha
parms$q <- 0.1824529
parms$w1 <- 0.1 # amplitude#5#94# 0.094 # amplitude (if 0 then not seasonal)
parms$w2 <- 0/12 #365/2 # shift (multiplied by pi)
parms$timestep <- 0.5
parms$intro <- 30   # frequency of infectious introduction
scalar <- list()
rr <- 2
for(j in 1:rr){
  scalar[[j]] <- 0.1 # 0.1
}
parms$scalar <- scalar
beta <- poly_fune(poly,parms,juvies,adults,post2021="same")
parms$beta <- beta
# and aging (which is assumed to be annual - so needs to be converted to days...)
parms$aging <- aging/(365)
parms$births <- (births)/(365)
parms$resist <- 0.20
# beta[adults,adults] <- 0

# run without any infections (to get a good pop size)
solch <- ode( y=yinit,times=seq(0,365*200,by=parms$timestep), func=ja.multistage.model.ii.seas,parms=parms )
# what does this look like? 
# 1. (check that aging happens sensibly)
# check that we maintain a pop size
suscets <- solch[,1+sindex] 
time <- solch[,1]
# and the age distribution?
equil <- drop(tail(solch[,1+sindex],1))[-1]
equil.pct <- equil/sum(equil)

# now add some infecteds - specify initial conditions
yinit <- drop(tail(solch[,],1))[-1]
yinit[iindex[2]] <- 1
yinit[sindex[2]] <- yinit[sindex[2]]-1

yrs <- 50
sol_old <- ode( y=yinit,times=seq(0,365*yrs,by=parms$timestep), func=ja.multistage.model.ii.seas,parms=parms )
sol <- sol_old

# *** check the incidence and the seasonality ***
like_old <- like_fun_harris(sol_old,data,sindex,eindex,iindex,aindex,rindex,gindex)
print(like_old)

# include a plot fo the model against the data
data$model <- cases_mod
data$age[2] <- c("05-14")

colors <- c("Model" = "grey60", "Data" = "blue")
m1 <- ggplot(data,aes(x=age,group=1)) + 
  geom_bar(aes(x=age,y=model,color="Model"),stat="identity",fill="grey60") + 
  geom_point(aes(x=age,y=obrien_avg,color="Data"),size=3) +
  geom_errorbar(aes(ymin = obrien_lwr, ymax = obrien_upr,color="Data"),width = 0.1) +
  labs(y = "Incidence (per 1,000 person-years)",
       x = "Age",
       color = "Legend") + 
  scale_color_manual(values = colors) +
  theme_bw()
  
# pdf("model_fit_harris_15June21.pdf",height=4,width=5)
# m1
# dev.off()

# timeline takes up the 3 years 2020 to 2023..
timeline <- c("01/01/2020","23/03/2020","04/06/2020","30/07/2020","04/09/2020","05/11/2020","03/12/2020","20/12/2020",
              "05/01/2021","09/03/2021","17/05/2021","17/07/2021","01/07/2023")
tstep <- as.integer(diff(as.Date(timeline,"%d/%m/%Y")))
tvals <- as.numeric(as.Date(timeline,"%d/%m/%Y")-as.Date("01/01/2020","%d/%m/%Y"))+365*yrs
# and now specific when we add infecteds
ivals <- seq(min(tvals),max(tvals),60) # new infection every 60 days
tcomm <- data.frame(x=c(tvals,ivals),y=c(rep("matrix",length(tvals)),rep("intro",length(ivals))))
# remove the intros where matrix exists
tcomm$z <- 0
for(i in 1:dim(tcomm)[1]){
  if(tcomm$y[i]=="intro"){
    oo <- which(tcomm$x == tcomm$x[i])
    if(length(oo)>1){
      tcomm$z[oo] <- 1
    }
  }
}
# remove the repeats
aa <- which(tcomm$z==1 & tcomm$y=="intro")
tcomm <- tcomm[-aa,]
tcomm2 <- tcomm[order(tcomm$x),]
# so use tcomm2 to work through the commands...

tunits <- length(tstep)
# need to refer to the right matrix...
#comix.period.adjust <- c(comix.period,"4. School reopening","4. School reopening")
comix.period.adjust <- c(comix.period,"4. School reopening","Polym")
#comix.period.adjust[4] <- "3. Relaxed restrictions"
#comix.period.adjust[11] <- "4. School reopening"
#<error checking>comix.period.adjust[4] <- "3. Relaxed restrictions"  # just to see what is it about school reopening...
# add comix.period to tcomm2
tcomm2$mat <- NA
ii <- 1
for(i in 2:dim(tcomm2)[1]){
  if(tcomm2$y[i]=="matrix"){
    tcomm2$mat[i] <- comix.period.adjust[ii]
    ii <- ii+1
  }
}

# check we have the right transmission coeficients
beta <- poly_fune(poly,parms,juvies,adults,post2021 = "same")
parms$beta <- beta
sol_old <- ode( y=yinit,times=seq(0,tcomm2$x[2],by=parms$timestep), func=ja.multistage.model.ii.seas,parms=parms )
# first round of sims
hold <- c(dim(sol_old)[1]-(5*365*(1/parms$timestep)),dim(sol_old)[1]) # hold the last 5 years
sol <- sol_old[hold[1]:hold[2],]

summary(sol[,1])/365  # check the time window for this
yinit_tt <- c(drop(tail(sol,1))[-1])  # initial conditions for future sims

casesrp <- sol_old[,1+ccindex] 
tmp <- (apply(casesrp,1,sum))  # sum (across ages) to get cases at each timestep
tt <- seq(0,dim(sol_old)[1],365*2)  # index of when a year goes by...
#diff(tmp[tt])/nall[tt[3:length(tt)]] # annual cases - turnover is too high by 10:1
#plot(diff(tmp[tt]))
# if we were to look at monthly incidence - look at the first year
ttmth <- round(seq(tt[49],tt[50],length.out = 13),0)
mthincid <- diff(tmp[ttmth])  # incidence for each month

# ok...so we can just churn through this loop...
post2021 <- "reduction"  # "same" #  **this refers to assumptions around contacts after 17th july**
for(tt in 2:(length(tcomm2$x)-1)){
  if(tcomm2$y[tt]=="intro"){
    # adjust the initial conditions but not matrice
    if(yinit_tt[sindex[2]]>1){
      yinit_tt[sindex[2]] <- yinit_tt[sindex[2]]-1
      yinit_tt[iindex[2]] <- yinit_tt[iindex[2]]+1
    }
  }
  if(tcomm2$y[tt]=="matrix" & tcomm2$mat[tt]!="Polym"){
    # use a different matrix
    cmat <- comix2[comix2$period==tcomm2$mat[tt],]
    print(table(cmat$period))
    cmat <- subset(cmat,select=c("contactee","contact","vals"))
    beta <- poly_fune(poly=cmat,parms,juvies,adults)  # cmat vs poly
    parms$beta <- beta
  }
  if(tcomm2$y[tt]=="matrix" & tcomm2$mat[tt]=="Polym" & post2021=="same"){
    # use a different matrix
    #cmat <- comix2[comix$period==tcomm2$mat[tt],]
    #print(table(cmat$period))
    #cmat <- subset(cmat,select=c("contactee","contact","vals"))
    beta <- poly_fune(poly=poly,parms,juvies,adults,post2021)
    parms$beta <- beta
  }
  if(tcomm2$y[tt]=="matrix" & tcomm2$mat[tt]=="Polym" & post2021=="reduction"){
    # use a different matrix
    #cmat <- comix2[comix$period==tcomm2$mat[tt],]
    #print(table(cmat$period))
    #cmat <- subset(cmat,select=c("contactee","contact","vals"))
    beta <- poly_fune(poly=poly,parms,juvies,adults,post2021=post2021)
    parms$beta <- beta
  }
  # with the updated params...run the model
  #sol_tt <- ode( y=yinit_tt,times=seq(0,tstep[tt],by=parms$timestep), func=ja.multistage.model.ii.seas,parms=parms )
  sol_tt <- ode( y=yinit_tt,times=seq(tcomm2$x[tt],tcomm2$x[tt+1],by=parms$timestep), func=ja.multistage.model.ii.seas,parms=parms )
  
  sol <- rbind(sol,sol_tt[2:dim(sol_tt)[1],])  # add to the bottom
  yinit_tt <- c(drop(tail(sol_tt,1))[-1])      # reset the starting conditions for the next round
  if(any(yinit_tt[sindex[1]:rindex[8]] < -0.5)){
    browser()
  }
  print(tt)
}

sol_df <- data.frame(sol)
names(sol_df)
max(sol_df$time)

# save the simulation
#save(sol_df,file="noro_sim_sc_20UP_20rho_A20_80_01Jul21.Rdat")

casesrp <- sol_df[,1+ccindex] 
incid <- apply(casesrp,1,sum)  # all cases

# look at daily cases (2020-100=1920)
tkeep <- seq(1,dim(sol_df)[1],round(1/parms$timestep))
length(tkeep)/365
#diff(incid[tvals])
plot((sol_df$time[tkeep[2:length(tkeep)]]/365)+1970,diff(incid[tkeep]),type="l",lwd=2,
     #ylim=c(0,60),
     xlim=c(2018,2023),
     ylab="Daily incidence (per 100,000)",xlab=" ")
abline(v=(cumsum(tstep)/365)+2020,lty=2,col="grey80")
abline(v=c(2015:2022),lty=1,col="red",lwd=0.5)
text(x=(cumsum(tstep[1:11])/365)+2020,y=50,labels=comix.period.adjust,pos=3,srt=90,cex=0.7)

# look at the age distribution of incidence...
# casesrp returns the cum incidence per age group
# we want for 3 time periods: 2019, 2020, 2021
summary((sol_df$time/365)+1920)
sol_df$time2 <- (sol_df$time/365)+1920
#
p1 <- which(sol_df$time2==2019)
p2 <- which(sol_df$time2==2020)
p3 <- which(sol_df$time2==2021)
p4 <- which(sol_df$time2==2022)
p5 <- which(sol_df$time2==2023)
calc_incid <- function(ivals){
  out <- rep(0,dim(ivals)[2])
  for(i in 1:dim(ivals)[2]){
    out[i] <- sum(diff(ivals[,i]))
  }
  return(out)
}

tmp <- casesrp[p1:p2,]
i2019 <- calc_incid(ivals=tmp)
tmp <- casesrp[p2:p3,]
i2020 <- calc_incid(ivals=tmp)
tmp <- casesrp[p3:p4,]
i2021 <- calc_incid(ivals=tmp)
tmp <- casesrp[p4:p5,]
i2022 <- calc_incid(ivals=tmp)
dat <- data.frame(expand.grid(age=ages,year=c(2019:2022)))
dat$incid <- c(i2019,i2020,i2021,i2022)
ggplot(dat,aes(x=as.factor(age),y=incid,group=year,fill=as.factor(age))) + 
  geom_bar(stat="identity") + facet_wrap(~year) +
  xlab("Age") + ylab("Incidence ('100,000 pop)") + labs(fill = "Age")
# total in each group...

# compare 2020 to to 2019
mean(100-(100*(dat$incid[dat$year==2020]/dat$incid[dat$year==2019])))
mean(100-(100*(dat$incid[dat$year==2021]/dat$incid[dat$year==2019])))
mean(100-(100*(dat$incid[dat$year==2022]/dat$incid[dat$year==2019])))

# proportion susceptible at end of each year?
suscets <- sol[,1+sindex] 
infects <- sol[,1+iindex] 
exposds <- sol[,1+eindex] 
recovds <- sol[,1+rindex] 
asympts <- sol[,1+aindex] 
nall <- suscets+infects+exposds+recovds+asympts
s2019 <- suscets[p2,]/nall[p2,]
s2020 <- suscets[p3,]/nall[p3,]
s2021 <- suscets[p4,]/nall[p4,]
s2022 <- suscets[p5,]/nall[p5,]
a2019 <- asympts[p2,]/nall[p2,]
dat$sprev <- c(s2019,s2020,s2021,s2022)
r2019 <- recovds[p2,]/nall[p2,]
r2020 <- recovds[p3,]/nall[p3,]
r2021 <- recovds[p4,]/nall[p4,]
r2022 <- recovds[p5,]/nall[p5,]
dat$rprev <- c(r2019,r2020,r2021,r2022)

head(dat)
ggplot(dat,aes(x=as.factor(age),y=sprev,group=year,fill=as.factor(age))) + 
  geom_bar(stat="identity") + facet_wrap(~year) + ylab("Proportion Susceptible") +
  scale_fill_brewer(palette="Reds") 
ggplot(dat,aes(x=as.factor(age),y=rprev,group=year,fill=as.factor(age))) + 
  geom_bar(stat="identity") + facet_wrap(~year) + ylab("Proportion Recovered (Sero +)") +
  scale_fill_brewer(palette="Greens") 







# end
