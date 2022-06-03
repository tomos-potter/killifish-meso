
# set length centering value
length.cen <- 6.5

#These functions are required in the likelihood function. They create the interaction surface.
myFunctionLPLP <- nimbleFunction(
  run = function(H = double(0), LP = double(0), philength.LP = double(0), rholength.LP = double(0), density.LP = double(0), le.comp.LP = double(1), le.focal.LPLP = double(1)) {
    numNumber <- LP * (sum( ((exp(philength.LP * le.comp.LP - philength.LP * le.focal.LPLP) ) *
                               exp((-(rholength.LP * le.comp.LP - rholength.LP * le.focal.LPLP)^2) / (2*H^2) ))[1:density.LP]
    ) ) / 1.4
    
    returnType(double(0))
    return(numNumber)
  }
)

myFunctionKOKO <- nimbleFunction(
  run = function(H = double(0), KO = double(0), philength.KO = double(0), rholength.KO = double(0), density.KO = double(0), le.comp.KO = double(1), le.focal.KOKO = double(1)) {
    numNumber <- KO * (sum( ((exp(philength.KO * le.comp.KO - philength.KO * le.focal.KOKO) ) *
                               exp((-(rholength.KO * le.comp.KO - rholength.KO * le.focal.KOKO)^2) / (2*H^2) ))[1:density.KO]
    ) ) / 1.4
    
    returnType(double(0))
    return(numNumber)
  }
)

#This is the Nimble model
code1 <- nimbleCode({
  #=====================================================================================
  #Priors for betas
  #=====================================================================================
  for(j in 1:8) {
    g.beta[j] ~ dnorm(0, sd = 100)
  }
  for(j in 1:6) {
    d.beta[j] ~ dnorm(0, sd = 100)
  }
  for(j in 1:6){
    f.beta[j] ~ dnorm(0, sd = 100)
  }
  
  g.sigma ~ dunif(0, 100)
  d.sigma ~ dunif(0, 100)
  
  
  #=====================================================================================
  #Priors for interaction surface
  #=====================================================================================
  philength.LP <- 0#~ dnorm(0.03, 0.010) #From Analysis with Interaction
  philength.KO <- 0#~ dnorm(0.03, 0.017) #From Analysis with Interaction
  rholength.LP <- 0#~ T(dnorm(0.03, 0.005), 0,1) #From Analysis with Interaction
  rholength.KO <- 0# ~ T(dnorm(0.03, 0.004), 0,1) #From Analysis with Interaction
  H <- 1 #~ dnorm(1, 100)
  
  #=====================================================================================
  #Priors for random effects
  #=====================================================================================
  
  ## random intercepts
  # random intercept of channel on growth
  g.sigma.int.c ~ dunif(0, 100)
  g.tau.int.c <- 1/(g.sigma.int.c*g.sigma.int.c)
  for(ii in 1:n_chan) {
    g.alpha.channel[ii] ~ dnorm(0, g.tau.int.c)
  }
  
  # These are for random intercepts of channel on fecundity and offspring size
  # f.sigma.int.c ~ dunif(0, 100)
  # f.tau.int.c <- 1/(f.sigma.int.c*f.sigma.int.c)
  # for(ii in 1:n_chan) {
  #   f.alpha.channel[ii] ~ dnorm(0, f.tau.int.c)
  # }
  
  # g.sigma.int.d ~ dunif(0, 100)
  # g.tau.int.d <- 1/(g.sigma.int.d*g.sigma.int.d)
  # for(ii in 1:n_drain) {
  #   g.alpha.drainage[ii] ~ dnorm(0, g.tau.int.d)
  # }
  
  
  #=====================================================================================
  #Likelihood
  #=====================================================================================
  
  # THE COMPETITION INTERACTION SURFACE
  
  for(i in 1:n_obs) {
    intsurf.LPLP[i] <-  myFunctionLPLP(1, LP[i], philength.LP, rholength.LP, density.LP[i], le.comp.LP[, chan.num[i]],le.focal.LPLP[i,])
    
    intsurf.KOKO[i] <- myFunctionKOKO(1, KO[i], philength.KO, rholength.KO, density.KO[i], le.comp.KO[, chan.num[i]],le.focal.KOKO[i,])
    
    
    ## GROWTH
    
    
    g.mu[i] <- KO[i] * (g.beta[1] + g.beta[2]*initial.length[i] + g.beta[3]*initial.length.squ[i] + g.beta[4]*intsurf.KOKO[i]) +
      LP[i] * (g.beta[5] + g.beta[6]*initial.length[i] + g.beta[7]*initial.length.squ[i] + g.beta[8]*intsurf.LPLP[i]) +
      g.alpha.channel[chan.num[i]] 
    
    
    growth[i] ~ dnorm(g.mu[i],sd=g.sigma)
    
    
    ## FECUNDITY
    
    log(f.mu[i]) <- KO[i] * (f.beta[1] + f.beta[2]*initial.length[i]  + f.beta[3]*intsurf.KOKO[i] ) +
      LP[i] * (f.beta[4] + f.beta[5]*initial.length[i]  + f.beta[6]*intsurf.LPLP[i] )# +
    #f.alpha.channel[chan.num[i]] #+ disp[i]
    
    
    fec[i] ~ dpois(f.mu[i])
    
    ## OFFSPRING SIZE
    
    d.mu[i] <-  KO[i] * (d.beta[1] + d.beta[2]*initial.length[i] + d.beta[3]*intsurf.KOKO[i]) +
      LP[i] * (d.beta[4] + d.beta[5]*initial.length[i] + d.beta[6]*intsurf.LPLP[i]) 
    
    
    off.length[i] ~ dnorm(d.mu[i],sd=d.sigma)
    
    
  }
})


##########################################################
#These are required to make the interaction surface
##########################################################
# pre-allocate a list and fill it with a loop
lengths.list <- vector("list", length(unique(gdata$channel.num)))
lengths.list <- setNames(lengths.list,unique(gdata$channel.num))
for (i in unique(gdata$channel.num)) {
  lengths.list[[i]] <- as.numeric(gdata[which(gdata$channel.num==i),'initial.total.length'])
}


#Need to cross list these with LP and KO individuals  
le.comp.LP <- array(0,c(24,length(unique(gdata$channel.num))))
for (i in 1:length(unique(gdata$channel.num))){
  for (j in 1:length(as.numeric(gdata[which(gdata$channel.num==i),'initial.total.length']))){
    le.comp.LP[j,i] <- as.numeric(gdata[which(gdata$channel.num==i & gdata$LP==1),'initial.total.length'])[j]-length.cen 
  }
}
le.comp.KO <- array(0,c(24,length(unique(gdata$channel.num))))
for (i in 1:length(unique(gdata$channel.num))){
  for (j in 1:length(as.numeric(gdata[which(gdata$channel.num==i),'initial.total.length']))){
    le.comp.KO[j,i] <- as.numeric(gdata[which(gdata$channel.num==i & gdata$KO==1),'initial.total.length'])[j]-length.cen 
  }
}
le.comp.LP[is.na(le.comp.LP)==T] <- 0
le.comp.KO[is.na(le.comp.KO)==T] <- 0

le.focal.LPLP <- array(0,c(dim(gdata)[1],24))
for (k in 1:dim(gdata)[1]){
  i <- gdata[k,'channel.num']
  # for (j in 1:length(as.numeric(gdata[which(gdata$channel.num==i),'initial.total.length']))){
  if (as.numeric(gdata[k,'LP'])==1){
    le.focal.LPLP[k,c(1:length(as.numeric(gdata[which(gdata$channel.num==i & gdata$LP==1),'initial.total.length'])))] <- as.numeric(gdata[k,'initial.total.length']) -length.cen
  }
  # }
}
le.focal.KOLP <- array(0,c(dim(gdata)[1],24))
for (k in 1:dim(gdata)[1]){
  i <- gdata[k,'channel.num']
  # for (j in 1:length(as.numeric(gdata[which(gdata$channel.num==i),'initial.total.length']))){
  if (as.numeric(gdata[k,'LP'])==1){
    le.focal.KOLP[k,c(1:length(as.numeric(gdata[which(gdata$channel.num==i & gdata$KO==1),'initial.total.length'])))] <- as.numeric(gdata[k,'initial.total.length']) -length.cen
  }
  # }
}

le.focal.KOKO <- array(0,c(dim(gdata)[1],24))
for (k in 1:dim(gdata)[1]){
  i <- gdata[k,'channel.num']
  # for (j in 1:length(as.numeric(gdata[which(gdata$channel.num==i),'initial.total.length']))){
  if (as.numeric(gdata[k,'KO'])==1){
    le.focal.KOKO[k,c(1:length(as.numeric(gdata[which(gdata$channel.num==i & gdata$KO==1),'initial.total.length'])))] <- as.numeric(gdata[k,'initial.total.length']) -length.cen
  }
  # }
}
le.focal.LPKO <- array(0,c(dim(gdata)[1],24))
for (k in 1:dim(gdata)[1]){
  i <- gdata[k,'channel.num']
  # for (j in 1:length(as.numeric(gdata[which(gdata$channel.num==i),'initial.total.length']))){
  if (as.numeric(gdata[k,'KO'])==1){
    le.focal.LPKO[k,c(1:length(as.numeric(gdata[which(gdata$channel.num==i & gdata$LP==1),'initial.total.length'])))] <- as.numeric(gdata[k,'initial.total.length']) -length.cen
  }
  # }
}


growth = as.numeric(gdata[,'growth']) 
fec = as.integer(gdata[,'num.embryos'])
fec[which(is.na(fec)==TRUE & gdata$final.sex=='F')] <- 0
off.length = as.numeric(gdata[,'off.length']) 
surv = as.numeric(gdata[,'survival']) 
prob.repro = (gdata[,'pregnant'])
n_block=max(as.numeric(gdata[,'block']))
n_chan=max(as.numeric(gdata[,'channel.num']))
n_drain=max(as.numeric(gdata[,'drain.num']))
n_obs = as.numeric(dim(gdata)[1])
chan.num=as.numeric(gdata[,'channel.num']) 
block.num=as.numeric(gdata[,'block']) 
drain.num=as.numeric(gdata[,'drain.num'])
initial.length = as.numeric(gdata[,'initial.total.length'] - length.cen)
initial.length.squ = (as.numeric(gdata[,'initial.total.length'] - length.cen))^2
stage = as.numeric(gdata[,'stage'])
le.comp.LP = le.comp.LP
le.comp.KO = le.comp.KO
le.focal.LPLP = le.focal.LPLP
le.focal.KOLP = le.focal.KOLP
le.focal.LPKO = le.focal.LPKO
le.focal.KOKO = le.focal.KOKO
LP = as.numeric(gdata[,'LP']) 
KO = as.numeric(gdata[,'KO'])
density.LP = as.numeric(gdata[,'density.LP'])
density.LP[which(density.LP==0)] <- 1
density.KO = as.numeric(gdata[,'density.KO'])
density.KO[which(density.KO==0)] <- 1
both.pheno = as.numeric(gdata[,'both.pheno'])

# constants
constants <- list(
  chan.num = chan.num,
  block.num = block.num,
  drain.num = drain.num,
  n_obs = n_obs,
  initial.length = initial.length,
  initial.length.squ = initial.length.squ,
  stage = stage,
  KO = KO,
  LP = LP,
  both.pheno = both.pheno,
  density.LP = density.LP,
  density.KO = density.KO,
  le.comp.LP = le.comp.LP,
  le.comp.KO = le.comp.KO,
  le.focal.LPLP = le.focal.LPLP,
  le.focal.KOLP = le.focal.KOLP,
  le.focal.LPKO = le.focal.LPKO,
  le.focal.KOKO = le.focal.KOKO
)

data <- list(
  growth = growth,
  fec = fec,
  off.length = off.length
)

inits <- list(
  g.beta = rep(0, 8),
  f.beta = rep(0, 6),
  d.beta = rep(0, 6),
  # philength.LP = 0,#rnorm(1, 0, 0.25),
  # philength.KO = 0,#rnorm(1, 0, 0.25),
  # rholength.LP = 0,#rnorm(1, 0, 0.25),
  # rholength.KO = 0,#rnorm(1, 0, 0.25),
  # H=1,
  g.sigma = 1,
  d.sigma = 1,
  
  #Random effects
  g.sigma.int.c = 1,
  g.alpha.channel = rep(0, n_chan)
  # f.sigma.int.c = 1,
  # f.alpha.channel = rep(0, n_chan)
  # g.sigma.int.d = 1,
  # g.alpha.drain = rep(0, n_drain)
  
)


