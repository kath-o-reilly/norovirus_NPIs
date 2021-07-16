
poly_fune <- function(poly,parms,juvies,adults,post2021 = "same"){
  #browser()
  if(sum(names(poly)=="vals")==0){
    warning("check names of poly")
  }
  if(post2021=="same"){
    x1 <- with( poly,
                tapply(vals,list(contact,contactee),unique)
    )
  }
  if(post2021=="reduction"){
    x1 <- with( poly,
                tapply(vals80,list(contact,contactee),unique)
    )
  }
  xsym <- ((x1+t(x1))/2)  # make contacts smetrical
  # convert to beta with transmission coefficinet
  beta <- parms$q*(xsym)  # default (matrix(rep(1,dim(poly)[1]),nrow=sqrt(dim(poly)[1]))#)
  q1 <- parms$q*parms$scalar[[1]]  # adult transmission
  q2 <- parms$q  # 
  beta[adults,adults] <- xsym[adults,adults]*q1
  
  return(beta)
}

ja.multistage.model.ii.seas.B <- function (t, x, parms) {
  s <- x[sindex]   # susceptibles
  e <- x[eindex]   # exposed
  i <- x[iindex]   # infecteds
  a <- x[aindex]   # asymptomatics
  r <- x[rindex]    # recovereds
  g <- x[gindex]
  n <- s+e+i+a+r+g
  Z <- 1 + parms$w1*cos((2*pi*t)/364 + parms$w2)
  beta_t <- parms$beta*Z
  lambda <- beta_t%*%((i+0.05*a+0.05*e)/n)
  #lambda <- beta_t%*%(i/n) 
  #lambda <- parms$beta%*%(i/n)
  #lambda <- parms$beta%*%((i+0.05*a+0.05*e)/n)   # force of infection (keep it simple for now)
  #browser()
  dsdt <- -lambda*s +(1/parms$theta)*r +parms$aging%*%s   # 
  dedt <- +lambda*s -parms$mu_s*e +parms$aging%*%e
  didt <- +parms$mu_s*e -parms$mu_a*i +parms$aging%*%i
  dadt <- +parms$mu_a*i -parms$rho*a  +parms$aging%*%a +parms$alpha*lambda*r
  drdt <-  +parms$rho*a -(1/parms$theta)*r  +parms$aging%*%r -parms$alpha*lambda*r
  dgdt <-  +parms$aging%*%g  # they just age
  # births - now 3-ways
  dgdt[1] <- dgdt[1]+parms$births*(parms$resist)  # genetically resistant don't contribute to transmission
  dsdt[1] <- dsdt[1]+parms$births*0.5*(1-parms$resist)  
  drdt[1] <- drdt[1]+parms$births*0.5*(1-parms$resist) 
  cc <- lambda*s# parms$mu_s*e
  aa <- parms$alpha*lambda*r
  
  list(
    c( dsdt,
       dedt,
       didt, 
       dadt,
       drdt,
       dgdt,
       cc,
       aa))
}

ja.multistage.model.ii.seas <- function (t, x, parms) {
  s <- x[sindex]   # susceptibles
  e <- x[eindex]   # exposed
  i <- x[iindex]   # infecteds
  a <- x[aindex]   # asymptomatics
  r <- x[rindex]    # recovereds
  g <- x[gindex]
  n <- s+e+i+a+r+g
  Z <- 1 + parms$w1*cos((2*pi*t)/364 + parms$w2*pi)
  beta_t <- parms$beta*Z
  #lambda <- beta_t%*%((i+0.05*a+0.05*e)/n)
  lambda <- beta_t%*%((i+0.05*a+0.05*e)/n)
  #lambda <- beta_t%*%(i/n) 
  #lambda <- parms$beta%*%(i/n)
  #lambda <- parms$beta%*%((i+0.05*a+0.05*e)/n)   # force of infection (keep it simple for now)
  #browser()
  dsdt <- -lambda*s +(1/parms$theta)*r +parms$aging%*%s   # 
  dedt <- +lambda*s -parms$mu_s*e +parms$aging%*%e
  didt <- +parms$mu_s*e -parms$mu_a*i +parms$aging%*%i
  dadt <- +parms$mu_a*i -parms$rho*a  +parms$aging%*%a +parms$alpha*lambda*r
  drdt <-  +parms$rho*a -(1/parms$theta)*r  +parms$aging%*%r -parms$alpha*lambda*r
  dgdt <-  +parms$aging%*%g  # they just age
  dsdt[1] <- dsdt[1]+parms$births*(1-parms$resist)
  dgdt[1] <- dgdt[1]+parms$births*(parms$resist)  # genetically resistant don't contribute to transmission
  cc <- lambda*s# parms$mu_s*e
  aa <- parms$alpha*lambda*r
  
  list(
    c( dsdt,
       dedt,
       didt, 
       dadt,
       drdt,
       dgdt,
       cc,
       aa))
}

like_fun_harris2 <- function(sol,data,sindex,eindex,iindex,aindex,rindex,gindex){
  # outputs
  suscets <- sol[,1+sindex] 
  casesrp <- sol[,1+ccindex] 
  asymprp <- sol[,1+aaindex] 
  asympts <- sol[,1+aindex] 
  recovds <- sol[,1+rindex] 
  genetds <- sol[,1+gindex]
  equil <- drop(tail(sol,1))[-1]
  n <- equil[sindex]+equil[eindex]+equil[iindex]+equil[aindex]+equil[rindex]+equil[gindex]  # total in each age class.
  
  yrs <- floor(dim(sol)[1]/((1/parms$timestep)*365))
  # incidence in the last year
  casesrp_ann <- casesrp[((yrs-1)*(365*(1/parms$timestep))):(yrs*(365*(1/parms$timestep))),] #casesrp[3650:7301,]
  cases_by_age <- ((casesrp_ann[dim(casesrp_ann)[1],]-casesrp_ann[1,])/n)*1000
  
  #              0-4, 5-14         15-65 (-24 to -64)       65+
  cases_mod <- c(cases_by_age[1:2],sum(cases_by_age[3:7])/5,cases_by_age[8])
  sprev <- 100*(sum(equil[rindex])/sum(n))  # as a pct from 0-100...
  aprev <- 100*(sum(equil[c(eindex,aindex,iindex)])/sum(n)) # a and i really but add in e for good measure
  rest <- 1/(sum(equil[sindex])/sum(n))
  
  # seasonality of cases
  incid <- apply(casesrp,1,sum)  # all cases
  ttmth <- round(seq((yrs-1)*(1/parms$timestep)*365,yrs*(1/parms$timestep)*365,length.out = 13),0)
  mthincid <- diff(incid[ttmth])  # incidence for each month
  # % of cases that are in months Nov-Jan (should be about 75%)
  prop_win <- sum(mthincid[c(1,11,12)])/sum(mthincid)
  
  # likelihood..we have 2 incidences per 1,000 but based on different data sizes

  ll1 <- sum(dnorm(data$obrien_avg,cases_mod,log=T))
  # and can we also fit to seroprevalence?
  ll2 <- dnorm(50,sprev,log=T)
  # output both separately
  out <- list()
  out[[1]] <- ll1
  out[[2]] <- ll2
  out[[3]] <- mean(cases_mod)
  out[[4]] <- sprev
  out[[5]] <- aprev
  out[[6]] <- rest
  out[[7]] <- mthincid
  out[[8]] <- prop_win
  return(out)
}


like_fun_harris <- function(sol,data,sindex,eindex,iindex,aindex,rindex,gindex){
  # outputs
  casesrp <- sol[,1+ccindex] 
  asymprp <- sol[,1+aaindex] 
  recovds <- sol[,1+rindex] 
  equil <- drop(tail(sol,1))[-1]
  n <- equil[sindex]+equil[eindex]+equil[iindex]+equil[aindex]+equil[rindex]+equil[gindex]  # total in each age class.
  
  yrs <- floor(dim(sol)[1]/((1/parms$timestep)*365))
  # incidence in the last year
  casesrp_ann <- casesrp[((yrs-1)*(365*(1/parms$timestep))):(yrs*(365*(1/parms$timestep))),] #casesrp[3650:7301,]
  cases_by_age <- ((casesrp_ann[dim(casesrp_ann)[1],]-casesrp_ann[1,])/n)*1000
  #              0-4, 5-14         15-65 (-24 to -64)       65+
  cases_mod <- c(cases_by_age[1:2],sum(cases_by_age[3:7])/5,cases_by_age[8])

  # likelihood..we have 2 incidences per 1,000 but based on different data sizes
  ll1 <- sum(dnorm(data$obrien_avg*1,cases_mod,log=T))
  out <- ll1
  return(out)
}

move_param <- function(parms,sol_old,sd_prop,sd_propq,sd_propt,select,yrs){
  # each iteration - select a parm to change
  ss <- round(runif(1,1,select),0)
  parms_old <- parms_new <- parms
  # priors for the full posterior
  p_1 <- dlnorm(parms_old$q,log(0.04),0.1,log=T) # hist(rlnorm(1000,log(0.2),0.5))
  p_2 <- dlnorm(parms_old$scalar[[1]],log(5),1,log=T) # hist(rlnorm(1000,log(5),1))
  p_3 <- 0#dlnorm(parms_old$scalar[[2]],log(1),1,log=T)
  p_4 <- dlnorm(parms_old$theta,log(5*365),0.5,log=T)  # hist(rlnorm(1000,log(5*365),0.5)/365)
  #p_4 <- 0#dlnorm(parms_old$w1,log(0.094),0.5,log=T)  # hist(rlnorm(1000,log(0.094),0.5))
  
  priors_old <- p_1 + p_2 + p_3 + p_4 
  
  if(ss == 1){ 
    # q (lognormal)
    new_val <- rlnorm(1,mean=log(parms$q), sd=sd_propq) # hist(rlnorm(1000,log(parms$q),sd_propq))
    # calculates probability of acceptance
    parms_new$q <- new_val  #
    beta <- poly_fune(poly,parms_new,juvies,adults)
    parms_new$beta <- beta 
    sol_new <- ode( y=yinit,times=seq(0,365*yrs,by=parms$timestep), func=ja.multistage.model.ii.seas.B,parms_new )  # incidence
  
    like_new <- like_fun_harris(sol_new,data,sindex,eindex,iindex,aindex,rindex,gindex)
    like_old <- like_fun_harris(sol_old,data,sindex,eindex,iindex,aindex,rindex,gindex)
    p_new <- dlnorm(parms_new$q,log(0.04),0.1,log=T)  
    priors_new <- p_new + p_2 + p_3 + p_4 
    # p_1 was here
    ratio_post <- (like_new+priors_new) - (like_old+priors_old)
    #correction <- 0 # no correction as move is symmetrical here (and 0 not 1 as on log scale)
    correction <- log(parms_new$q) - log(parms_old$q)  # lof scale
  }
  if(ss == 2){ 
    # theta (lognormal)
    new_val <- rlnorm(1,mean=log(parms$theta), sd=sd_propt)
    # calculates probability of acceptance
    parms_new$theta <- new_val  #
    #beta <- poly_fune(poly,parms_new,juvies,adults)
    #parms_new$beta <- beta 
    sol_new <- ode( y=yinit,times=seq(0,365*yrs,by=parms$timestep), func=ja.multistage.model.ii.seas,parms_new )  # incidence
    
    like_new <- like_fun_harris(sol_new,data,sindex,eindex,iindex,aindex,rindex,gindex)
    like_old <- like_fun_harris(sol_old,data,sindex,eindex,iindex,aindex,rindex,gindex)
    p_new <- dlnorm(parms_new$theta,log(5*365),0.5,log=T)   # hist(rlnorm(1000,log(1),0.5))
    priors_new <- p_1 + p_2 + p_3 + p_new 
    # p_1 was here
    ratio_post <- (like_new+priors_new) - (like_old+priors_old)
    #correction <- 0 # no correction as move is symmetrical here (and 0 not 1 as on log scale)
    correction <- log(parms_new$theta) - log(parms_old$theta)  # lof scale
  }
  
  # consistent for all parameters...
  p_accept <- ratio_post + correction # things are additive here as on log scale
  if(is.na(p_accept)){browser()}
  if(p_accept>0){
    p_accept <- 0
  }
  # accept/reject step
  tmp <- log(runif(1))
  if(tmp<p_accept){# accepting with a certain probability
    parms <- parms_new       #updated.val <- new.val
    like_keep <- like_new 
    #post_keep <- like_new+priors_new
    priors_keep <- priors_new
    accept <- 1
    sol_old <- sol_new
  }else # reject
  {#updated.val <- curr.val
    like_keep <- like_old
    priors_keep <- priors_old
    #post_keep <- like_old+priors_old
    accept <- 0
    # keep parms
    # keep the solution
    # keep the priors
  }	
  # return a list with rlevant info
  out <- list()
  out$sol_old <- sol_old
  out$parms_old <- parms
  out$priors_keep <- priors_keep
  out$like_keep <- like_keep
  out$ss <- ss
  out$accept <- accept
  #browser()
  #rm(like_keep,like_new,like_old)
  
  #out <- c(unlist(parms),like_keep,like_keep+priors_keep,accept,ss)
  return(out)
}

# # end