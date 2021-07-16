# plot output of the 4 models for norovirus

library(deSolve)
library(ggplot2)   
library(reshape2)
library(janitor)
library(tidyverse)
library(qs) # comix data
library(RColorBrewer)
library(cowplot)

rm(list = ls())
setwd("~/Dropbox/Dropbox/Norovirus/Rcode/DeterministicModel")

load("noro_sim_sc_0UP_15rho_01Jul21.Rdat")
moda <- sol_df
load("noro_sim_sc_20UP_15rho_01Jul21.Rdat")
modb <- sol_df
load("noro_sim_sc_0UP_20rho_01Jul21.Rdat")
modc <- sol_df
load("noro_sim_sc_20UP_20rho_01Jul21.Rdat")
modd <- sol_df
timestep <- diff(moda$time[1:2])
ages <- c(4,14,24,34,44,54,64,80)

sindex <- 1:length(ages) 
eindex <- length(ages)+(1:length(ages)) # 31:60 # 
iindex <- (length(ages)*2)+(1:length(ages))  # 61:90 
aindex <- (length(ages)*3)+(1:length(ages))  #91:120 
rindex <- (length(ages)*4)+(1:length(ages))  #121:150
gindex <- (length(ages)*5)+(1:length(ages))  #151:150
ccindex <- (length(ages)*6)+(1:length(ages))

# check stuff in sims
casesrp <- moda[,1+ccindex] 
tmp <- (apply(casesrp,1,sum))  # sum (across ages) to get cases at each timestep

tt <- seq(0,dim(moda)[1],365*2)  # index of when a year goes by...
diff(tmp[tt])/nall[tt[3:length(tt)]] # annual cases - turnover is too high by 10:1
plot(diff(tmp[tt]))
# if we were to look at monthly incidence - look at the first year
ttmth <- round(seq(tt[1],tt[2],length.out = 13),0)
mthincid <- diff(tmp[ttmth])  # incidence for each month
sum(mthincid)
min(mthincid)
barplot(mthincid)

# timeline takes up the 2 years 2020 to 2021..
timeline <- c("01/01/2020","23/03/2020","04/06/2020","30/07/2020","04/09/2020","05/11/2020","03/12/2020","20/12/2020",
              "05/01/2021","09/03/2021","17/05/2021","19/07/2021","1/07/2023")
tstep <- as.integer(diff(as.Date(timeline,"%d/%m/%Y")))
ldlines <- (cumsum(tstep)/365)+2020

comix <- read.csv("~/Dropbox/Dropbox/Norovirus/Data/comix_adjusted_08May2021.csv") %>% clean_names()
comix.period <- names(table(comix$period))
comix.period.adjust <- c(comix.period,"10. School reopening","11. Polym")

# summary values
tmp <- comix[comix$period==comix.period.adjust[10],]
round(mean(tapply(tmp$estcon,tmp$contact,sum)),2)

# look at weekly cases (2020-100=1920) - will align better with official surveillance
tkeep <- seq(1,dim(moda)[1],7*round(1/timestep))
length(tkeep)

# plot together with dates for lockdown
# create a df that ggplot will like...
suscets <- moda[,1+sindex] 
infects <- moda[,1+iindex] 
exposds <- moda[,1+eindex] 
recovds <- moda[,1+rindex] 
asympts <- moda[,1+aindex] 
genetds <- moda[,1+gindex]
casesrp <- moda[tkeep,1+ccindex] 
incid <- apply(casesrp,1,sum)  # all cases
#tmp <- moda[tkeep,1+sindex]/nall
nall <- suscets+infects+exposds+asympts+recovds+genetds
suscept <- apply(suscets[tkeep,]/nall[tkeep,],1,mean) # average across all ages
length(suscept)
length(tkeep)
df <- data.frame(time=moda$time[tkeep][-1],incid=diff(incid),suscept=suscept[-1],sim=rep("A",length(tkeep)-1))
# repeat for the others
for(i in 1:3){
  if(i == 1){
    casesrp <- modb[tkeep,1+ccindex] 
    suscets <- modb[,1+sindex] 
    infects <- modb[,1+iindex] 
    exposds <- modb[,1+eindex] 
    recovds <- modb[,1+rindex] 
    asympts <- modb[,1+aindex] 
    genetds <- modb[,1+gindex]
    nall <- suscets+infects+exposds+asympts+recovds+genetds
    suscept <- apply(suscets[tkeep,]/nall[tkeep,],1,mean) # average across all ages
  }
  if(i == 2){
    casesrp <- modc[tkeep,1+ccindex] 
    suscets <- modc[,1+sindex] 
    infects <- modc[,1+iindex] 
    exposds <- modc[,1+eindex] 
    recovds <- modc[,1+rindex] 
    asympts <- modc[,1+aindex] 
    genetds <- modc[,1+gindex]
    nall <- suscets+infects+exposds+asympts+recovds+genetds
    suscept <- apply(suscets[tkeep,]/nall[tkeep,],1,mean) # average across all ages
  }
  if(i == 3){
    casesrp <- modd[tkeep,1+ccindex] 
    suscets <- modd[,1+sindex] 
    infects <- modd[,1+iindex] 
    exposds <- modd[,1+eindex] 
    recovds <- modd[,1+rindex] 
    asympts <- modd[,1+aindex] 
    genetds <- modd[,1+gindex]
    nall <- suscets+infects+exposds+asympts+recovds+genetds
    suscept <- apply(suscets[tkeep,]/nall[tkeep,],1,mean) # average across all ages
  }
  # if(i == 4){
  #   casesrp <- mode[tkeep,1+ccindex] 
  # }
  # if(i == 5){
  #   casesrp <- modf[tkeep,1+ccindex] 
  # }
  incid <- apply(casesrp,1,sum)  # all cases
  tmp <- data.frame(time=moda$time[tkeep][-1],incid=diff(incid),suscept=suscept[-1],sim=rep(toupper(letters[i+1]),length(tkeep)-1))
  df <- rbind(df,tmp)
}

df$time2 <- (df$time/365)+1970
summary(df$time2)
df$UP <- "0"
df$UP[df$sim=="B"|df$sim=="D"] <- "20"
df$shed <- "15"
df$shed[df$sim=="C"|df$sim=="D"] <- "20"

head(df)

RE <- 0.729

yxt <- 400
m1 <- ggplot(df,aes(x=(time/365)+1970,y=incid*RE,group=sim,col=sim)) + 
  geom_line(aes(linetype=UP, color=shed),lwd=1.5) +
  scale_color_brewer(palette="Reds") +
  xlab(" ") + ylab("Weekly incidence (per 100,000)") +
  xlim(2019,2023.6) +
  geom_vline(xintercept = ldlines[1:(length(ldlines)-1)], linetype="dotted", color = "grey50", size=0.7) +
  geom_text(x=2020.225,y=yxt,label=comix.period.adjust[1],col="black",angle = 90,hjust = 1) +
  geom_text(x=ldlines[1:(length(ldlines)-1)][2],y=yxt,label=comix.period.adjust[2],col="black",angle = 90,hjust = 1) +
  geom_text(x=ldlines[1:(length(ldlines)-1)][3],y=yxt,label=comix.period.adjust[3],col="black",angle = 90,hjust = 1) +
  geom_text(x=ldlines[1:(length(ldlines)-1)][4],y=yxt,label=comix.period.adjust[4],col="black",angle = 90,hjust = 1) +
  geom_text(x=ldlines[1:(length(ldlines)-1)][5],y=yxt,label=comix.period.adjust[5],col="black",angle = 90,hjust = 1) +
  geom_text(x=ldlines[1:(length(ldlines)-1)][6],y=yxt,label=comix.period.adjust[6],col="black",angle = 90,hjust = 1) +
  geom_text(x=ldlines[1:(length(ldlines)-1)][9],y=yxt,label=comix.period.adjust[9],col="black",angle = 90,hjust = 1) +
  geom_text(x=ldlines[1:(length(ldlines)-1)][10],y=yxt,label=comix.period.adjust[10],col="black",angle = 90,hjust = 1) +
  geom_text(x=ldlines[1:(length(ldlines)-1)][11],y=yxt,label=comix.period.adjust[11],col="black",angle = 90,hjust = 1) +
  theme_bw()
m2 <- ggplot(df,aes(x=(time/365)+1970,y=suscept,group=sim,col=sim)) +
  scale_color_brewer(palette="Reds") +
  geom_line(aes(linetype=UP, color=shed),lwd=1.5) + 
  xlab(" ") + ylab("Proportion susceptible to symptomatic infection ") +
  xlim(2019,2023.6) +
  geom_vline(xintercept = ldlines[1:(length(ldlines)-1)], linetype="dotted", color = "grey50", size=0.7) +
  theme_bw()

# incidence alone
# basic plot - Figure 2
pdf("Simulations_2018-2022_scUP_rho_01Jul21.pdf",height=5,width=10)
m1
dev.off()

# add in susceptibles  - Figure 2
pdf("Simulations_sus_2018-2022_scUP_rho_101Jul21.pdf",height=8,width=10)
plot_grid(m1,m2,labels = c('A', 'B'),ncol=1,align = 'v', axis = 'tblr' )
dev.off()

# add have barcharts of incidence by age...
# it's best to have by noro year, starting in wk 27...
# using first simulation
sol_df <- modd
  
sol_df$time2 <- (sol_df$time/365)+1970
casesrp <- sol_df[,1+ccindex]   # for entire timeseries
suscets <- sol_df[,1+sindex] 
p1 <- which((sol_df$time2-(2018+(27/52)))^2 == min((sol_df$time2-(2018+(27/52)))^2)) # jul 18
p2 <- which((sol_df$time2-(2019+(27/52)))^2 == min((sol_df$time2-(2019+(27/52)))^2)) # jul 19
p3 <- which((sol_df$time2-(2020+(27/52)))^2 == min((sol_df$time2-(2020+(27/52)))^2)) # jul 20
p4 <- which((sol_df$time2-(2021+(27/52)))^2 == min((sol_df$time2-(2021+(27/52)))^2)) # # jul 21
p5 <- which((sol_df$time2-(2022+(27/52)))^2 == min((sol_df$time2-(2022+(27/52)))^2))# jul 22
p6 <- which((sol_df$time2-(2023+(27/52)))^2 == min((sol_df$time2-(2023+(27/52)))^2)) # jul 23
calc_incid <- function(ivals){
  # ivals should be c incidence by age group
  out <- rep(0,dim(ivals)[2])
  for(i in 1:dim(ivals)[2]){
    out[i] <- sum(diff(ivals[,i]))
  }
  return(out)
}

tmp <- casesrp[p1:p2,]
i2018 <- calc_incid(ivals=tmp)
tmp <- casesrp[p2:p3,]
i2019 <- calc_incid(ivals=tmp)
tmp <- casesrp[p3:p4,]
i2020 <- calc_incid(ivals=tmp)
tmp <- casesrp[p4:p5,]
i2021 <- calc_incid(ivals=tmp)
tmp <- casesrp[p5:p6,]
i2022 <- calc_incid(ivals=tmp)
dat <- data.frame(expand.grid(sim=c("D"),age=ages,year=c(2018:2022)))
dat$incid <- c(i2018,i2019,i2020,i2021,i2022)

# for each age group report % difference compared to 2019
dat$i2018 <- rep(dat$incid[dat$year==2018],times=5)
dat$diff <- round((100*(dat$incid/dat$i2018)),1)
#View(dat)
write.csv(dat,file="~/Dropbox/Dropbox/Norovirus/Data/ddump.csv",row.names=F)

m3 <- ggplot(dat[dat$year!=2018,],aes(x=as.factor(age),y=incid*RE,group=year,fill=as.factor(age))) + 
  geom_bar(stat="identity") + facet_wrap(~year) +
  scale_fill_brewer(palette="Reds") +
  xlab("Age group") + ylab("Annual Incidence ('100,000 pop)") + labs(fill = "Age group") + theme_bw()


# save to file  ~/Dropbox/Dropbox/Norovirus/Data/
write.csv(dat,file="~/Dropbox/Dropbox/Norovirus/Data/ddump.csv",row.names=F)

# and we want to report the difference in susceptibility
# we want at each time point what the % sus is.
nall <- moda[,sindex+1]+moda[,iindex+1]+moda[,eindex+1]+moda[,aindex+1]+moda[,rindex+1]+moda[,gindex+1]  # total in each age
moda$time2 <- (moda$time/365)+1970
min(moda$time2)
p1 <- which((moda$time2-(2018+(27/52)))^2 == min((moda$time2-(2018+(27/52)))^2)) # jul 18
p2 <- which((moda$time2-(2019+(27/52)))^2 == min((moda$time2-(2019+(27/52)))^2)) # jul 19
p3 <- which((moda$time2-(2020+(27/52)))^2 == min((moda$time2-(2020+(27/52)))^2)) # jul 20
p4 <- which((moda$time2-(2021+(27/52)))^2 == min((moda$time2-(2021+(27/52)))^2)) # # jul 21
p5 <- which((moda$time2-(2022+(27/52)))^2 == min((moda$time2-(2022+(27/52)))^2))# jul 22
p6 <- which((moda$time2-(2023+(27/52)))^2 == min((moda$time2-(2023+(27/52)))^2)) # jul 23

i2018 <- moda[p1,1+sindex]/nall[p1,]  # at 2018.5 this is the Prop Sus 
i2019 <- moda[p2,1+sindex]/nall[p2,]
i2020 <- moda[p3,1+sindex]/nall[p3,]
i2021 <- moda[p4,1+sindex]/nall[p4,]
i2022 <- moda[p5,1+sindex]/nall[p5,]
dat <- data.frame(expand.grid(age=ages,year=c(2018:2022)))
dat$sus_pct <- unlist(c(i2018,i2019,i2020,i2021,i2022))
head(dat)
dat$i2018 <- rep(dat$sus_pct[dat$year==2018],times=5)
dat$diff <- round((100*(dat$sus_pct/dat$i2018)),1)

write.csv(dat,file="~/Dropbox/Dropbox/Norovirus/Data/ddump.csv",row.names=F)



# total in each group...

# what would the corresponding national surveillance look like?
sol_df <- df[df$sim=="A",]
summary(sol_df$incid) # so in a population of 100,000 these are the likely incidence
tmp <- 56290000/100000 # scalar (sims were 100,000) THIS IS FOR England 
sol_df$incid_scale <- sol_df$incid*tmp   # in the country
summary(sol_df$incid_scale)
p1 <- which((sol_df$time2-(2018+(27/52)))^2 == min((sol_df$time2-(2018+(27/52)))^2)) # jul 18
p2 <- which((sol_df$time2-(2019+(27/52)))^2 == min((sol_df$time2-(2019+(27/52)))^2)) # jul 19
p3 <- which((sol_df$time2-(2020+(27/52)))^2 == min((sol_df$time2-(2020+(27/52)))^2)) # jul 20
p4 <- which((sol_df$time2-(2021+(27/52)))^2 == min((sol_df$time2-(2021+(27/52)))^2)) # # jul 21
p5 <- which((sol_df$time2-(2022+(27/52)))^2 == min((sol_df$time2-(2022+(27/52)))^2))# jul 22
p6 <- which((sol_df$time2-(2023+(27/52)))^2 == min((sol_df$time2-(2023+(27/52)))^2)) # jul 23
sum(sol_df$incid_scale[p1:p2])  # cases (perfect surveillance) in the country in a year
(sum(sol_df$incid_scale[p1:p2])/56290000)*1000 # incidence per 1,000 person-years

sum(sol_df$incid[p1:p2])  # cases in a population of 100,000
sum(sol_df$incid[p1:p2])/100  # incidence per 1,000 person-years (just checking...)
sol_df$surv <- 0 
sol_df$surv_lwr <- 0
sol_df$surv_upr <- 0
prob <- 1/287.6      # Pr(case reported to surveillance) 
prob_lwr <- 1/346
prob_upr <- 1/239.1
for(i in 1:dim(sol_df)[1]){
  sol_df$surv[i] <- round(sol_df$incid_scale[i])*prob # not going to sample
  # confidence intervals
  # np+_sqrt((p/(1-p))/n)
  sol_df$surv_lwr[i] <- round(sol_df$incid_scale[i])*prob_lwr-sqrt((prob_lwr/(1-prob_lwr))/round(sol_df$incid_scale[i]))
  sol_df$surv_upr[i] <- round(sol_df$incid_scale[i])*prob_upr+sqrt((prob_upr/(1-prob_upr))/round(sol_df$incid_scale[i]))
}
summary(sol_df$surv)
summary(sol_df$surv_lwr)
summary(sol_df$surv_upr)

# combine in a plot...
sol_df$scenario <- "Typical Year"
sol_df$scenario[p2:p3] <- "2019/20"
sol_df$scenario[p3:p4] <- "2020/21"
sol_df$scenario[p4:p5] <- "2021/22"
sols <- sol_df[p1:p5,]
sols$week <- round(((sols$time2 - floor(sols$time2))*52),0)+1# 
sols$week_noro <- sols$week-26
sols$week_noro[sols$week_noro<0] <- sols$week[sols$week_noro<0]+26
sols$week_noro[sols$week==26] <- 52
#View(sols)

RE <- 0.729

head(sols)
m4 <- ggplot(sols[sols$scenario!="2019/20" & sols$scenario!="2020/21",],aes(x=week_noro,y=surv*RE,group=scenario,col=scenario)) + 
  geom_point() +
  theme_bw() + ylab("Estimated cases submitted to SGSS") + xlab("Norovirus Week Number (from late June") + xlim(0,52)
m4

tmp <- sols[sols$scenario=="2021",]

#View(tmp)
# what week is the lowest?
oo <- which(sols$surv[sols$scenario=="Typical Year"]==min(sols$surv[sols$scenario=="Typical Year"]))
sols[oo,] # week 31

dates <- as.Date(c("01/08/2021","01/09/2021","01/10/2021","01/11/2021","01/12/2021","01/01/2022","01/02/2022","01/03/2022","01/04/2022","01/05/2022","01/06/2022"),"%d/%m/%Y")
dates_wk <- as.numeric(format(dates,"%V"))-27
dates_wk[dates_wk<1] <- as.numeric(format(dates[dates_wk<1],"%V"))+25
dates_vals <- c(month.name[8:12],month.name[1:6])
tmp <- data.frame(wk=dates_wk,label=dates_vals)

# should add in July 19th and today...
s1 <- as.numeric(format(as.Date("2021-06-16"),"%V")) - 27
s2 <- as.numeric(format(as.Date("2021-07-19"),"%V")) - 27
m5 <- ggplot(sols[sols$scenario!="2019/20" & sols$scenario!="2020/21",],aes(x=week_noro,y=surv*RE,group=scenario,col=scenario)) + 
  geom_point() + 
  ylab("Estimated cases submitted to SGSS") + 
  xlab("Norovirus Week Number (from Week 27)") + 
  geom_line(data=sols[sols$scenario!="2019/20" & sols$scenario!="2020/21",],aes(x=week_noro,y=surv_lwr*RE)) +
  geom_line(data=sols[sols$scenario!="2019/20" & sols$scenario!="2020/21",],aes(x=week_noro,y=surv_upr*RE)) +
  scale_color_manual(values=c("#9970ab","#5ab4ac")) +
  geom_vline(xintercept=dates_wk,col="grey80",linetype="dashed") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlim(0,51) + ylim(0,1250)
  #geom_vline(xintercept=c(s1,s2),linetype="dotted", color = "grey50", size=0.7)
#m5 +  geom_label(data=tmp,aes(x=wk,y=370,label=label),angle = 90,vjust = "inward",col="grey50")

pdf("predict_surveillance_01Jul21.pdf",height=4,width=5)
m5
dev.off()

# add in the annual surveillance - figure 3
pdf("predict_surveillance_age_16June21.pdf",height=4,width=10)
plot_grid(m3,m5,labels = c('A', 'B'),ncol=2)
dev.off()

# add in SGSS data
sdat <- read.csv("~/Dropbox/Dropbox/Norovirus/Data/PHEdata/Noro_weekly_counts_02Jul21.csv") %>% clean_names()
head(sdat)

# average of the reports
table(sdat$year)
sdat$time <- round(sdat$year+(sdat$reporting_week-1)/52,2)  # calendar time
sdat$time_noro <- sdat$time - round((27-1)/52,2)   # noro weeks, starting in wk 27
sdat$year_noro <- as.numeric(substr(sdat$time_noro,1,4))
sdat$week_noro <- round((sdat$time_noro - sdat$year_noro)*52,0)+1
# remove noro year 2013... and include only up to 2018
sdat2 <- sdat[sdat$year_noro>2013 & sdat$year_noro<=2018,]
table(sdat2$year_noro)  # 5 years of data
# weekly summary
wdat <- sdat2 %>% group_by(week_noro) %>% summarize(avg=mean(count),min=min(count),max=max(count),n=n())
head(wdat)
tdat <- sols[sols$scenario=="Typical Year",]
head(tdat)
table(tdat$week_noro)
# merge...
wdat <- left_join(wdat,tdat,by=c("week_noro"))

ggplot(wdat,aes(x=week_noro,y=avg)) + geom_line() +
  geom_line(aes(x=week_noro,y=min),col="blue") +
  geom_line(aes(x=week_noro,y=max),col="blue") +  # ok great
  geom_line(data=tdat,aes(x=week_noro,y=surv),col="red")

# so from this data what is RE?
print(mean(wdat$avg/wdat$surv))
RE <- round(mean(wdat$avg/wdat$surv),3)  # 0.729

m6 <- ggplot(wdat,aes(x=week_noro,y=avg)) + geom_line(col="#a6611a",lwd=1.5) +
  geom_vline(xintercept=dates_wk,col="grey80",linetype="dashed") +
  geom_text(data=tmp,aes(x=wk,y=370,label=label),angle = 90,vjust = "inward",col="grey50") +
  geom_line(aes(x=week_noro,y=min),col="#dfc27d") +
  geom_line(aes(x=week_noro,y=max),col="#dfc27d") +  # ok great
  geom_line(data=tdat,aes(x=week_noro,y=surv*RE),col="#5ab4ac",lwd=1.5) +
  geom_line(data=tdat,aes(x=week_noro,y=surv_lwr*RE),col="#80cdc1") +
  geom_line(data=tdat,aes(x=week_noro,y=surv_upr*RE),col="#80cdc1") +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("Estimated cases submitted to SGSS") + xlab("Norovirus Week (from Week 27)") + 
  xlim(1,51) + ylim(0,400)

# just validation alone - Figure 2
pdf("Model_SGSS_compare_02Jul21.pdf",height=4,width=5)
m6
dev.off()

# How does 2020/21 data compare to the model?
# actually, include the year before too...
sdat3 <- sdat %>% filter(year_noro==2020 | year_noro==2019)
sols3 <- sols %>% filter(scenario=="2020/21" | scenario=="2019/20")
# need calendar time

label <- c(month.name[8:12],month.name[1:6])
m7 <- ggplot(sdat3,aes(x=week_noro,y=count)) + geom_line(col="#a6611a",lwd=1.5) +
  geom_vline(xintercept=dates_wk,col="grey80",linetype="dashed") +
  geom_text(data=tmp,aes(x=wk,y=290,label=label),angle = 90,vjust = "inward",col="grey50") +  
  geom_line(data=sols3,aes(x=week_noro,y=surv*RE),col="#5ab4ac",lwd=1.5) +
  geom_line(data=sols3,aes(x=week_noro,y=surv_lwr*RE),col="#80cdc1") +
  geom_line(data=sols3,aes(x=week_noro,y=surv_upr*RE),col="#80cdc1") + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("Cases submitted to SGSS (2020/21)") + xlab("Norovirus Week (from Week 27)") +
  ylim(0,320) #+ xlim(1,50)

tmp2 <- rbind(tmp,c(51,"July"),tmp)
tmp2$wk <- as.numeric(tmp2$wk)
tmp2$time <- 2019.5+round((tmp2$wk-1)/52,2)
tmp2$time[13:23] <- tmp2$time[13:23]+1

m7 <- ggplot(sdat3,aes(x=time,y=count)) + geom_line(col="#a6611a",lwd=1.5) +
  geom_vline(xintercept=tmp2$time,col="grey80",linetype="dashed") +
  geom_text(data=tmp2,aes(x=time,y=310,label=label),angle = 90,vjust = "inward",col="grey50") +  
  geom_line(data=sols3,aes(x=time2,y=surv*RE),col="#5ab4ac",lwd=1.5) +
  geom_line(data=sols3,aes(x=time2,y=surv_lwr*RE),col="#80cdc1") +
  geom_line(data=sols3,aes(x=time2,y=surv_upr*RE),col="#80cdc1") + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("Cases submitted to SGSS (2019-21)") + xlab("Time (from Norovirus Week 27)") +
  ylim(0,350) + xlim(2019.5,2021.4)

pdf("Model_SGSS_compare2019-21_02Jul21.pdf",height=4,width=8)
m7
dev.off()







