#### Jolly-Seber example 
library(openpopscr)
RcppParallel::setThreadOptions(numThreads = 10)
source("inst/len_tests/simulate_SR.R")

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 2000, lambda0 = 1.0, sigma = 30, phi = 0.7, beta = 0.9)

# make detectors array
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "proximity")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 10

# simulate ScrData 
scrdat <- simulate_js_openscr(par = true_par, n_occasions = n_occasions, 
                              n_sec_occasions = n_occasions,
                              detectors = detectors, mesh = mesh, seed = 19592,
                              temporary_emigration = matrix(c(0.9, 0.1, 0.6, 0.4), nrow = 2, byrow = TRUE))



# fit model ---------------------------------------------------------------

# create formulae 
par <- list(lambda0 ~ avail, 
            sigma ~ avail, 
            beta ~ 1, 
            phi ~ 1, 
            D ~ 1)
par2 <- list(lambda0 ~ as.factor(avail), 
            sigma ~ as.factor(avail), 
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
                           structure = matrix(c(".", "~1", 
                                                "~1", "."), nr = 2, nc = 2, byrow = T), 
                           start = list(delta = c(0.5, 0.5), 
                                        tpm = matrix(c(0.8, 0.2,
                                                       0.2, 0.8), nr = 2, nc = 2, byrow = T)), 
                           cov = data.frame(avail = c(1, NA)))

statemod2 <- StateModel$new(data = scrdat, 
                           names = c("available", "available2"), 
                           structure = matrix(c(".", "~1", 
                                                "~1", "."), nr = 2, nc = 2, byrow = T), 
                           start = list(delta = c(0.5, 0.5), 
                                        tpm = matrix(c(0.8, 0.2,
                                                       0.2, 0.8), nr = 2, nc = 2, byrow = T)), 
                           cov = data.frame(avail = c(1, 2)))


# create model object 
oo <- JsModel$new(par, scrdat, start, statemod = statemod)
oo2 <- JsModel$new(par2, scrdat, start, statemod = statemod2)
# compute initial likelihood 
oo$calc_llk()

# fit model 
nlm_args <- list(iterlim = 150, hessian = FALSE, gradtol = 1e-3)
oo$fit(nlm_args)
oo2$fit(nlm_args)
# see results 
oo
prod(exp(oo$estimates()$par[1:2,1]))
true_par
oo$get_par("D")
oo$get_par("lambda0", k = 1, j = 1)
oo$get_par("sigma", k = 1, j = 1)
oo$get_par("phi", k = 1, m = 1)
oo$get_par("beta", m = 1)

oo$state()$tpm()
oo$state()$delta()
oo$state()$par()
oo$state()$nstates()
oo$state()$groups()


oo2
prod(exp(oo2$estimates()$par[1:2,1]))
true_par
oo2$get_par("D")
oo2$get_par("lambda0", k = 1, j = 1)
oo2$get_par("sigma", k = 1, j = 1)
oo2$get_par("phi", k = 1, m = 1)
oo2$get_par("beta", m = 1)


#------------------------------------------------------------------------

statemodseason <- StateModel$new(data = scrdat, 
                                 names = c("available", "unavailable"), 
                                 structure = matrix(c(".", "~season", 
                                                      "~season", "."), nr = 2, nc = 2, byrow = T), 
                                 start = list(delta = c(0.5, 0.5), 
                                              tpm = matrix(c(0.8, 0.2,
                                                             0.2, 0.8), nr = 2, nc = 2, byrow = T)), 
                                 cov = data.frame(avail = c(1, NA)))
scrdat$add_covariate("season", factor(c(1,2,3,4,1,2,3,4,1,2)), "k")


oo_season <- JsModel$new(par, scrdat, start, statemod = statemodseason)
oo_season$fit(nlm_args)