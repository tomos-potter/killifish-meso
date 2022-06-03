#=====================
# Killifish vital rates 
#    MAIN SCRIPT
#=====================
# Uploads data, fits Bayesian vital rate models simultaneously using NIMBLE
# model diagnostics, model selection, etc...
#=====================
#Load necessary libraries and functions
library(Rcpp)
library(devtools)
require(nimble)
library(plyr)
library(MCMCvis)

source('./upload_killifish_data.R')

#===============================================================================================
#Model without any phis, rhos, or H's.
#===============================================================================================

# MCMC parameters (yields 2000 samples of the posterior)
ni <- 70000 # total iterations per chain
nb <-  20000 # burn-in period
nt <- 100 # thinning interval (to save memory)
nc <- 4 # number of chains to run

params <- c('g.beta','f.beta','d.beta',
            'g.sigma','d.sigma', 'g.sigma.int.c')

# NARANJO
gdata <- upload_killifish_data("naranjo") # upload data
source('./nimble_killifish_model_null.R') # upload model code

# build the model
Rmodel <- nimbleModel(code1, constants, data, inits)
conf <- configureMCMC(Rmodel, monitors = params, enableWAIC = TRUE)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel,showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel,showCompilerOutput = TRUE)

# fit the model
set.seed(0)
samples_null_naranjo <- runMCMC(Cmcmc, niter = ni, nburnin = nb, thin = nt,
                                nchains = nc, WAIC = TRUE)

# plots of posteriors and trace to check convergence etc
MCMCtrace(samples_null_naranjo$samples, 
          params = params, pdf=F, ind = T, Rhat = T, n.eff = T)

# summary of parameter estimates
samplesSummary(samples_null_naranjo$samples[[1]])

# get WAIC
samples_null_naranjo$WAIC

save(samples_null_naranjo, file="samples_null_naranjo.Rda")

# GUANAPO
gdata <- upload_killifish_data("guanapo") # upload data
source('./nimble_killifish_model_null.R') # upload model code

# build the model
Rmodel <- nimbleModel(code1, constants, data, inits)
conf <- configureMCMC(Rmodel, monitors = params, enableWAIC = TRUE)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel,showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel,showCompilerOutput = TRUE)

# fit the model
set.seed(0)
samples_null_guanapo <- runMCMC(Cmcmc, niter = ni, nburnin = nb, thin = nt,
                                nchains = nc, WAIC = TRUE)

# plots of posteriors and trace to check convergence etc
MCMCtrace(samples_null_guanapo$samples, 
          params = params, pdf=F, ind = T, Rhat = T, n.eff = T)

# summary of parameter estimates
samplesSummary(samples_null_guanapo$samples[[1]])

# get WAIC
samples_null_guanapo$WAIC

save(samples_null_guanapo, file="samples_null_guanapo.Rda")

#===============================================================================================
#GLMs for comparison (just checking for Guanapo)
# this is just to see that we get similar results in a frequentist framework
#===============================================================================================
gdata$initial.total.length.cen <- gdata$initial.total.length - length.cen
gdata$initial.total.length.cen.squ <- (gdata$initial.total.length.cen)^2
gdata$density <- (gdata$density.KO + gdata$density.LP) / 1.4
gdata[which(gdata$KO==1),'pheno'] <- 'KO'
gdata[which(gdata$LP==1),'pheno'] <- 'LP'

gdata$fec <- gdata$num.embryos
gdata[which(is.na(gdata$fec)==TRUE & gdata$final.sex=='F'),'fec'] <- 0

#check growth model
mod <- glmmTMB(growth ~ 0 + pheno + pheno:initial.total.length.cen + 
             pheno:initial.total.length.cen.squ + pheno:density + (1|channel.num),
           data=gdata)

summary(mod) # compare with parameter estimates from NIMBLE model for growth

# check fecundity model
mod <- glmmTMB(fec ~ 0 + pheno + pheno:initial.total.length.cen + pheno:density,
           family='poisson',
           data=gdata)

summary(mod) # compare with parameter estimates from NIMBLE model for fecundity
#===============================================================================================
# Full model with population specific values of phi and rho
#===============================================================================================

params <- c('g.beta','f.beta','d.beta',
            'philength.LP','philength.KO' ,
            'rholength.LP','rholength.KO',
            'g.sigma','d.sigma', 'g.sigma.int.c')

ni <- 1200000
nb <-  200000
nt <- 1000
nc <- 4

 # NARANJO
gdata <- upload_killifish_data("naranjo") # upload data

source('./nimble_killifish_model_full.R')

Rmodel <- nimbleModel(code1, constants, data, inits)
conf <- configureMCMC(Rmodel, monitors = params, enableWAIC = TRUE)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel,showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel,showCompilerOutput = TRUE)

set.seed(0)
samples_full_naranjo <- runMCMC(Cmcmc, niter = ni, nburnin = nb, thin = nt, nchains = nc)

# plots of posteriors and trace to check convergence etc
MCMCtrace(samples_full_naranjo$samples, 
          params = params, pdf=F, ind = T, Rhat = T, n.eff = T)

# summary of parameter estimates
samplesSummary(samples_full_naranjo$samples[[1]])

# get WAIC
samples_full_naranjo$WAIC

save(samples_full_naranjo, file="samples_full_naranjo.Rda")

# GUANAPO
gdata <- upload_killifish_data("guanapo") # upload data

source('./nimble_killifish_model_full.R')

Rmodel <- nimbleModel(code1, constants, data, inits)
conf <- configureMCMC(Rmodel, monitors = params, enableWAIC = TRUE)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel,showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel,showCompilerOutput = TRUE)

set.seed(0)
samples_full_guanapo <- runMCMC(Cmcmc, niter = ni, nburnin = nb, thin = nt, nchains = nc)

# plots of posteriors and trace to check convergence etc
MCMCtrace(samples_full_guanapo$samples, 
          params = params, pdf=F, ind = T, Rhat = T, n.eff = T)

# summary of parameter estimates
samplesSummary(samples_full_guanapo$samples[[1]])

# get WAIC
samples_full_guanapo$WAIC

save(samples_full_guanapo, file="samples_full_guanapo.Rda")


#===============================================================================================
#GLMs
#===============================================================================================


myFunctionLPLP <- function(H, LP, philength.LP, rholength.LP,  le.comp.LP, le.focal.LPLP) {
    numNumber <- LP * (sum( ((exp(philength.LP * le.comp.LP - philength.LP * le.focal.LPLP) ) *
                               exp((-(rholength.LP * le.comp.LP - rholength.LP * le.focal.LPLP)^2) / (2*H^2) ))#[1:density.LP]
    ) ) / 1.4
    

    return(numNumber)
  }


myFunctionKOKO <- function(H , KO , philength.KO , rholength.KO , le.comp.KO , le.focal.KOKO ) {
    numNumber <- KO * (sum( ((exp(philength.KO * le.comp.KO - philength.KO * le.focal.KOKO) ) *
                               exp((-(rholength.KO * le.comp.KO - rholength.KO * le.focal.KOKO)^2) / (2*H^2) ))#[1:density.KO]
    ) ) / 1.4
    

    return(numNumber)
  }



for (i in 1:dim(gdata)[1]){
  if(gdata[i,'LP']==1){
    gdata[i,'eq.dens'] <- myFunctionLPLP(H=mean(samples2[[chain]][,'H']), gdata[i,'LP'], philength.LP=mean(samples2[[chain]][,'philength.LP']), rholength.LP=mean(samples2[[chain]][,'rholength.LP']), gdata[which(gdata$channel.num==gdata[i,'channel.num']),'initial.total.length.cen'], gdata[i,'initial.total.length.cen'])
  }
  if(gdata[i,'LP']==0){
    gdata[i,'eq.dens'] <- myFunctionKOKO(H=mean(samples2[[chain]][,'H']), gdata[i,'KO'], philength.KO=mean(samples2[[chain]][,'philength.KO']), rholength.KO=mean(samples2[[chain]][,'rholength.KO']), gdata[which(gdata$channel.num==gdata[i,'channel.num']),'initial.total.length.cen'], gdata[i,'initial.total.length.cen'])
  }
}

#===============================================================================================
#GLMs for comparison
#===============================================================================================

gdata[which(is.na(gdata$fec)==TRUE & gdata$final.sex=='F'),'fec'] <- 0


mod <- glmmTMB(growth ~ 0 + pheno + pheno:initial.total.length.cen + 
                 pheno:initial.total.length.cen.squ + pheno:eq.dens + 
                 (1|channel.num),data=gdata)
summary(mod)


mod <- glmmTMB(fec ~ 0 + pheno + pheno:initial.total.length.cen + pheno:eq.dens,family='poisson',data=gdata)
summary(mod)

#===============================================================================================
# Model with common values of phi and rho
#===============================================================================================

params <- c('g.beta','f.beta','d.beta',
            'philength','rholength',
            'g.sigma','d.sigma', 'g.sigma.int.c')

ni <- 1200000
nb <-  200000
nt <- 1000
nc <- 4

# NARANJO
gdata <- upload_killifish_data("naranjo") # upload data

source('./nimble_killifish_model_common_surface.R')

Rmodel <- nimbleModel(code1, constants, data, inits)
conf <- configureMCMC(Rmodel, monitors = params, enableWAIC = TRUE)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel,showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel,showCompilerOutput = TRUE)

set.seed(0)
samples_common_naranjo <- runMCMC(Cmcmc, niter = ni, nburnin = nb, thin = nt, nchains = nc)

# plots of posteriors and trace to check convergence etc
MCMCtrace(samples_common_naranjo$samples, 
          params = params, pdf=F, ind = T, Rhat = T, n.eff = T)

# summary of parameter estimates
samplesSummary(samples_common_naranjo$samples[[1]])

# get WAIC
samples_common_naranjo$WAIC

save(samples_common_naranjo, file="samples_common_naranjo.Rda")

# GUANAPO
gdata <- upload_killifish_data("guanapo") # upload data

source('./nimble_killifish_model_common_surface.R')

Rmodel <- nimbleModel(code1, constants, data, inits)
conf <- configureMCMC(Rmodel, monitors = params, enableWAIC = TRUE)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel,showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel,showCompilerOutput = TRUE)

set.seed(0)
samples_common_guanapo <- runMCMC(Cmcmc, niter = ni, nburnin = nb, thin = nt, nchains = nc)

# plots of posteriors and trace to check convergence etc
MCMCtrace(samples_common_guanapo$samples, 
          params = params, pdf=F, ind = T, Rhat = T, n.eff = T)

# summary of parameter estimates
samplesSummary(samples_common_guanapo$samples[[1]])

# get WAIC
samples_common_guanapo$WAIC

save(samples_common_guanapo, file="samples_common_guanapo.Rda")







