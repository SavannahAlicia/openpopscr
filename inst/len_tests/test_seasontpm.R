#### Jolly-Seber example 
library(openpopscr)
RcppParallel::setThreadOptions(numThreads = 10)

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 6000, lambda0 = 1.0, sigma = 30, phi = 0.7, beta = 0.9)

# make detectors array
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "proximity")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 10

# set temporary emigration for seasons (consider 2 alternating)
temporary_emigration <- list(matrix(c(1, 0, 0.4, 0.6), nrow = 2, byrow = TRUE),
                             matrix(c(.84, .16, 0, 1), nrow = 2, byrow = TRUE)
                            )
temporary_emigration <- temporary_emigration[c(rep(c(1,2), 4), 1)]

# set primaries and secondaries
primary <- c(sort(rep(1:10, 3)))

# simulate ScrData 
scrdat <- simulate_js_openscr(par = true_par, primary = primary,
                              detectors = detectors, mesh = mesh, seed = 19592,
                              temporary_emigration = temporary_emigration)

scrdat$add_covariate("season", as.factor(rep(c(1,2), 5)), "p")

# fit model ---------------------------------------------------------------

# create formulae 
par <- list(lambda0 ~ avail, 
            sigma ~ avail, 
            beta ~ 1, 
            phi ~ 1, 
            D ~ 1)

# get start values 
start <- #get_start_values(scrdat, model = "JsModel")
  list(lambda0 = .9, 
       sigma = 20, 
       beta = .8, 
       phi =.8, 
       D =1500)
  
statemod <- StateModel$new(data = scrdat, 
                           names = c("available", "unavailable"), 
                           structure = matrix(c(".", "~season", 
                                                "~season", "."), nr = 2, nc = 2, byrow = T), 
                           start = list(delta = c(0.5, 0.5), 
                                        tpm = matrix(c(0.8, 0.2,
                                                       0.2, 0.8), nr = 2, nc = 2, byrow = T)), 
                           cov = data.frame(avail = c(1, NA)))


# create model object 
oo <- JsModel$new(par, scrdat, start, statemod = statemod)
# compute initial likelihood 
oo$calc_llk()

# fit model 
nlm_args <- list(iterlim = 150, hessian = FALSE, gradtol = 1e-3)
oo$fit(nlm_args)

# see results 
oo
prod(exp(oo$estimates()$par[1:2,1]))
true_par
oo$get_par("D")
oo$get_par("lambda0", k = 1, j = 1)
oo$get_par("sigma", k = 1, j = 1)
oo$get_par("phi", k = 1, m = 1)
oo$get_par("beta", m = 1)

#appears to be some label switching here?
temporary_emigration[1:2]
oostate$tpm(k=1)
oostate$tpm(k=2)

oostate <- oo$state()
oostate()$trm()
oostate()$delta()
oostate()$par()
oostate()$nstates()
oostate()$groups()
#state parameter estimates 1 through 5
boot::inv.logit(oostate$estimates()[1,1]) - oostate$delta()[2] < 1e-10 #first parameter estimated is delta (initial state 2 membership)
exp(oostate$estimates()[2,1]) - oostate$trm(k=1)[2,1] < 1e-10 #transition from state 2 to 1 in season 1
exp(oostate$estimates()[4,1]) - oostate$trm(k=1)[1,2] < 1e-10 #transition from state 1 to 2 in season 1
exp(oostate$estimates()[2,1] + oostate$estimates()[3,1]) - oostate$trm(k=2)[2,1] < 1e-10 #transition from state 2 to 1 in season 2
exp(oostate$estimates()[4,1] + oostate$estimates()[5,1]) - oostate$trm(k=2)[1,2] < 1e-10 #transition from state 1 to 2 in season 2


#just to check if the delta we simulated matches (can add specify delta to simulation)
m_temporary_emigration <- matrix(c(mean(unlist(lapply(temporary_emigration, function(x){ x[1,1]}))),
                                   mean(unlist(lapply(temporary_emigration, function(x){ x[1,2]}))),
                                   mean(unlist(lapply(temporary_emigration, function(x){ x[2,1]}))),
                                   mean(unlist(lapply(temporary_emigration, function(x){ x[2,2]})))), nrow = 2, byrow = TRUE) 

stat_dist <- eigen(t(m_temporary_emigration))$vector[, 1]
stat_dist <- stat_dist / sum(stat_dist)

#again, label switching. Labels seem consistent throughout model output, but not between simulation and model output
oostate$delta()
stat_dist

#-------------------- check order of estimates wihtout NA (this model doesn't make sense, but just to see if estimate order changes)
#use factor for avail instead of NA
statemod2 <- StateModel$new(data = scrdat, 
                           names = c("available", "unavailable"), 
                           structure = matrix(c(".", "~season", 
                                                "~season", "."), nr = 2, nc = 2, byrow = T), 
                           start = list(delta = c(0.5, 0.5), 
                                        tpm = matrix(c(0.8, 0.2,
                                                       0.2, 0.8), nr = 2, nc = 2, byrow = T)), 
                           cov = data.frame(avail = as.factor(c(1, 2))))


# create model object 
oo2 <- JsModel$new(par, scrdat, start, statemod = statemod2)
oo2$fit(nlm_args)
oo2
prod(exp(oo2$estimates()$par[1:2,1]))
true_par
oo2$get_par("D")
oo2$get_par("lambda0", k = 1, j = 1)
oo2$get_par("sigma", k = 1, j = 1)
oo2$get_par("phi", k = 1, m = 1)
oo2$get_par("beta", m = 1)

oo2state$trm()

temporary_emigration[1:2]
oo2state$tpm(k=1)
oo2state$tpm(k=2)

stat_dist
oo2state$delta() #hold up, now labels aren't switched from simulation...

oo2state$delta()[2] -  boot::inv.logit(oo2state$par()[1]) < 1e-10 #first parameter is initial membership of state 2 (delta 2)
exp(oo2state$par()[2]) - oo2state$trm(k=1)[2,1] < 1e-10 #second parameter is still transition rate from state 2 to state 1 in season 1
exp(oo2state$par()[4]) - oo2state$trm(k=1)[1,2] < 1e-10 #third parameter is tr from 1 to 2 in season 1


