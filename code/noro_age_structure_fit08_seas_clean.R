# use the code in King and Wearing to create SIR model with age structure, and use polymod data
# adapt for use applied to norovirus

library(deSolve)
library(ggplot2)   
library(reshape2)
library(janitor)

rm(list = ls())
setwd("~/Dropbox/Dropbox/Norovirus/Rcode/DeterministicModel")

# model attributes;
# - multiple age groups (30 age classes, first 20 years each have their own, further years of 5 yrs each)
# - bitrths to equal deaths
# - noro specific
births <- 1250

ages <- c(4,14,24,34,44,54,64,80) # upper end of age classes (for Amar data)
da <- diff(c(0,ages))
length(ages)
yinit <- c( S=c(14000,19050,rep(14500,6)),
            E=c(rep(0,length(ages))),
            I=c(rep(0,1),1,rep(0,6)),
            A=c(rep(0,length(ages))),
            R=c(rep(0,length(ages))),
            G=c(rep(0,length(ages)))
)
sum(yinit)
length(yinit)/6

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
# try the comixr package data
# we want a symmetric contact matrix
library(socialmixr)
socm <- contact_matrix(polymod, countries = "United Kingdom", 
                       age.limits = c(0,as.numeric(as.character(age.categories[1:(length(ages)-1)]))), 
                       symmetric = TRUE)
dim(socm$matrix)   # units = "mean number of contacts"
poly <- data.frame(expand.grid(contactee=age.categories,contact=age.categories))
poly$vals <- as.vector(socm$matrix)
#ggplot(poly,aes(y=contacts,x=participant,fill=vals)) + geom_tile() + 
#  scale_fill_gradient(low="white",high="red")
round(mean(tapply(poly$vals,poly$contactee,sum)),2)

# load in prevalence data from Harris
data <- read.csv("~/Dropbox/Dropbox/Norovirus/Data/harris_noro_incidence.csv") %>% clean_names()
head(data)

mu_s <- 1/1  # rate of pre-inf to infectious
mu_a <- 1/2  # rate of infetious to asymptomatic
rho <- 1/15  # rate of asymtomatic to recovery
theta <- (5.1*365)  # rate of a loss of immunity
alpha <- 1#0.5

ccindex <- (length(ages)*6)+(1:length(ages))
aaindex <- (length(ages)*7)+(1:length(ages))

source("noro_functions_mcmc_clean.r")

# use a list to contain the parameter values
parms <- list()
parms$beta <- beta
parms$mu_s <- mu_s
parms$mu_a <- mu_a
parms$rho <- 1/20#rho
parms$theta <- 5.1*365 #
parms$alpha <- 1 #0.5#alpha
parms$q <- 0.23
parms$w1 <- 0.15# (if 0 then not seasonal)
parms$w2 <- 2/12 # # shift (multiplied by pi in the ODEs)
parms$timestep <- 1
scalar <- list()
rr <- 2
for(j in 1:rr){
  scalar[[j]] <- 0.1 # adult to child ratio
}
parms$scalar <- scalar
beta <- poly_fune(poly,parms,juvies,adults)
parms$beta <- beta
# and aging (which is assumed to be annual - so needs to be converted to days...)
parms$aging <- aging/(365)
parms$births <- (births)/(365)
parms$resist <- 0.20

# nitial conditions for model
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

solch <- ode( y=yinit,times=seq(0,365*200,by=parms$timestep), func=ja.multistage.model.ii.seas.B,parms=parms )
# what does this look like? 
#
# check that we maintain a pop size
suscets <- solch[,1+sindex] 
time <- solch[,1]

# and the age distribution?
equil <- drop(tail(solch[,1+sindex],1))#[-1]
equil.pct <- equil/sum(equil)

# now add some infecteds - specify initial conditions
yinit <- drop(tail(solch[,],1))[-1]
yinit[iindex[2]] <- 1
yinit[sindex[2]] <- yinit[sindex[2]]-1

sol_old <- ode( y=yinit,times=seq(0,365*50,by=parms$timestep), func=ja.multistage.model.ii.seas.B,parms=parms )
sol <- sol_old

# compare data to the model output...
like_old <- like_fun_harris(sol_old,data,sindex,eindex,iindex,aindex,rindex,gindex)
print(like_old)
out <- like_fun_harris2(sol_old,data,sindex,eindex,iindex,aindex,rindex,gindex)
print(out)

plot(1:12,out[[7]])

# step 1 - using mcmc - assuming each value is independent
iter <- 2000
select <- 1#+rr # q+rr
sd_prop <- 0.015 # update sd for scalar
sd_propq <- 0.003  # update sd for q
sd_propt <- 0.5  # update sd for theta (0.1 if theta much smaler if w1)
yrs <- 50 # years to run the model
parm_names <- c("q","alpha","theta","w1",paste0(rep("scalar",rr),1:rr))
runs <- data.frame(matrix(rep(NA,(length(parm_names)+4)*iter),ncol=(length(parm_names)+4)))
dim(runs)
names(runs) <- c(parm_names,"like","post","accept","select")
head(runs)
runs[1,] <- c(parms$q,parms$alpha,parms$theta,parms$w1,unlist(parms$scalar),like_old,0,0,0)

# run through the iterations
for(i in 1:iter){
  
  out <- move_param(parms,sol_old,sd_prop,sd_propq,sd_propt,select,yrs)
  print(paste0(i," - like: ",round(out$like_keep,1)," , ",round(out$parms$q,4)))
  # returns a list
  runs[i,] <- c(out$parms_old$q,out$parms_old$alpha,parms$theta,parms$w1,unlist(out$parms_old$scalar),
                out$like_keep,  #like
                out$like_keep+out$priors_keep,  # post
                out$accept,out$ss)
  
  # next iter input is the output of this one
  sol_old <- out$sol_old
  parms <- out$parms_old
  rm(out)
  #print(paste0(i," - like: ",round(out$like_keep,1)))
}

# check what has happened here
oo <- c(2:i)
plot(runs$post[oo])
table(runs$select,runs$accept)
table(runs$select[oo],runs$accept[oo])
runs[i-1,]

aa <- which(runs$post==max(runs$post,na.rm=T))
runs[aa,]

sol <- sol_old

write.csv(tmp,"ddump.csv",row.names = F)
save(runs,file="model_run_best_noseas_19May21.rdata")

#load(file="model_run_harris_noreinf_imm8_2poly_27April21.rdata")

# plot the best one

parms <- list()
parms$beta <- beta
parms$mu_s <- mu_s
parms$mu_a <- mu_a
parms$rho <- 1/15#rho
parms$theta <- 5.1*365 #3181.696 #(~7 yrs)
parms$alpha <- 1 #0.5#alpha
parms$q <- 0.1909186
parms$w1 <- 0.1# 0.094 # amplitude#5#94# 0.094 # amplitude (if 0 then not seasonal)
parms$w2 <- 11/12 #365/2 # shift (multiplied by pi)
parms$timestep <- 1
scalar <- list()
rr <- 2
for(j in 1:rr){
  scalar[[j]] <- 0.1 #0.1
}
parms$scalar <- scalar
beta <- poly_fune(poly,parms,juvies,adults)
parms$beta <- beta
# and aging (which is assumed to be annual - so needs to be converted to days...)
parms$aging <- aging/(365)
parms$births <- (births)/(365)
parms$resist <- 0.20

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

solch <- ode( y=yinit,times=seq(0,365*200,by=parms$timestep), func=ja.multistage.model.ii.seas.B,parms=parms )

# now add some infecteds - specify initial conditions
yinit <- drop(tail(solch[,],1))[-1]
yinit[iindex[2]] <- 1
yinit[sindex[2]] <- yinit[sindex[2]]-1


sol <- ode( y=yinit,times=seq(0,365*50,by=parms$timestep), func=ja.multistage.model.ii.seas.B,parms=parms )

# check everything worked ok
out <- like_fun_harris2(sol,data,sindex,eindex,iindex,aindex,rindex,gindex)
print(out)
out <- like_fun_harris(sol,data,sindex,eindex,iindex,aindex,rindex,gindex)
print(out)

#sol <- sol_old

suscets <- sol[,1+sindex] 
infects <- sol[,1+iindex] 
exposds <- sol[,1+eindex] 
recovds <- sol[,1+rindex] 
asympts <- sol[,1+aindex] 
casesrp <- sol[,1+ccindex] 
asymprp <- sol[,1+aaindex] 

njuv <- suscets[,juvies]+infects[,juvies]+exposds[,juvies]+asympts[,juvies]+recovds[,juvies]
nadu <- apply(suscets[,adults]+infects[,adults]+exposds[,adults]+asympts[,adults]+recovds[,adults],1,sum)
nall <- apply(suscets+infects+exposds+asympts+recovds,1,sum)
time <- sol[,1]
equil <- drop(tail(sol,1))[-1]

tmp <- (apply(casesrp,1,sum))  # sum (across ages) to get cases at each timestep
#plot(time[1:length(time)],tmp)
#length(time)
tt <- seq(0,365*50*(1/parms$timestep),365)  # index of when a year goes by...
diff(tmp[tt])/nall[tt[3:length(tt)]] # annual cases - turnover is too high by 10:1
plot(diff(tmp[tt]))

# if we were to look at monthly incidence - look at the last year
ttmth <- round(seq(tt[48],tt[50],length.out = 13*2),0)
mthincid <- diff(tmp[ttmth])  # incidence for each month
sum(mthincid[c(1,11,12)])/sum(mthincid[1:12])
sum(mthincid)
barplot(mthincid)

sol <- sol_old
ages_names <- data$age

#tmp <- cases_mod
pdf("Noro_model_asymptom_50infectious_fit.pdf",height=4,width=5)
par(mfrow=c(1,1))
x1 <- barplot(height=cases_mod,width=da,
              names=ages_names,las=2,col="grey",
              #main="Norovirus incidence by age",
              ylab="Incidence (per 1000-pyrs)",
              ylim=c(0,250))
points(x1,data$obrien_avg,col="blue",pch=19,cex=2)
for(i in 1:4){
  lines(x=c(x1[i],x1[i]),y=c(data$obrien_lwr[i],data$obrien_upr[i]),col="blue",lwd=2)
}
legend("topright",c("Model","Data"),pch=c(15,19),col=c("grey","blue"))
dev.off()

names(mthincid) <- month.abb

pdf("Noro_model_fitq01_seas_14May21.pdf",height=5,width=8)
par(mfrow=c(1,2))
x1 <- barplot(height=cases_mod,width=da,xaxt = "n",
              #names=ages_names,las=2,
              col="grey",
              #main="Norovirus incidence by age",
              ylab="Incidence (per 1000-pyrs)",
              ylim=c(0,250))
points(x1,data$obrien_avg,col="blue",pch=19,cex=2)
for(i in 1:4){
  lines(x=c(x1[i],x1[i]),y=c(data$obrien_lwr[i],data$obrien_upr[i]),col="blue",lwd=2)
}
text(x = x1,
     y = par("usr")[3]-20 ,
     labels = ages_names,
     xpd = NA,     ## Rotate the labels by 35 degrees.
     srt = 45,
     cex = 1.0)
legend("topright",c("Model","Data"),pch=c(15,19),col=c("grey","blue"))
x2 <- barplot(mthincid/100,ylab="Monthly incidence ('1000 total pop)",xaxt = "n")
text(x = x2,
     y = par("usr")[3]-0.6 ,
     labels = names(mthincid),
     xpd = NA,
     srt = 45,## Rotate the labels by 35 degrees.
     cex = 1.0)
dev.off()

# end
