#### Jolly-Seber example 
library(openpopscr)
library(ggplot2)
RcppParallel::setThreadOptions(numThreads = 10)
source("inst/len_tests_simulate_LT.R")

# simulate data -----------------------------------------------------------

# set truth 
true_par <- list(D = 4000, lambda0 = 1.0, sigma = 30, phi = 0.7, beta = 0.5)

# make detectors array
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "proximity")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 10

# simulate ScrData 
scrdat <- simulate_js_openscr(par = true_par, n_occasions = n_occasions, 
                              n_sec_occasions = n_occasions,
                              detectors = detectors, mesh = mesh, seed = 19592)

scrdat$add_covariate("firstprim", factor(c(1,rep(0, 9))), "p")
scrdat$add_covariate("realtime", scrdat$time(), "p")

# fit model ---------------------------------------------------------------

# create formulae 
par <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            beta ~ realtime, 
            phi ~ 1, 
            D ~ 1)
par2 <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            beta ~ realtime + firstprim, 
            phi ~ 1, 
            D ~ 1)
par3 <- list(lambda0 ~ 1, 
             sigma ~ 1, 
             beta ~ s(realtime, k = 3) , 
             phi ~ 1, 
             D ~ 1)
par4 <- list(lambda0 ~ 1, 
             sigma ~ 1, 
             beta ~ s(realtime, k = 3) + firstprim, 
             phi ~ 1, 
             D ~ 1)
par5 <- list(lambda0 ~ 1, 
             sigma ~ 1, 
             beta ~ s(realtime, k = 6) , 
             phi ~ 1, 
             D ~ 1)
par6 <- list(lambda0 ~ 1, 
             sigma ~ 1, 
             beta ~ s(realtime, k = 6) + firstprim, 
             phi ~ 1, 
             D ~ 1)

# get start values 
start <- #get_start_values(scrdat, model = "JsModel")
  list(lambda0 = .9, 
       sigma = 20, 
       beta = .7, 
       phi =.8, 
       D =1500)
  

# create model object 
oo <- JsModel$new(par, scrdat, start)
# compute initial likelihood 
oo$calc_llk()

# fit model 
nlm_args <- list(iterlim = 150, hessian = FALSE)
oo$fit(nlm_args)

oo2 <- JsModel$new(par2, scrdat, start)
oo2$fit(nlm_args)
oo3 <- JsModel$new(par3, scrdat, start)
oo3$fit(nlm_args)
oo4 <- JsModel$new(par4, scrdat, start)
oo4$fit(nlm_args)
oo5 <- JsModel$new(par5, scrdat, start)
oo5$fit(nlm_args)
oo6 <- JsModel$new(par6, scrdat, start)
oo6$fit(nlm_args)
# see results 

plotdat <- data.frame(beta = c(oo$get_par("beta", m = 1)[-1],
                               oo2$get_par("beta", m = 1)[-1],
                               oo3$get_par("beta", m = 1)[-1],
                               oo4$get_par("beta", m = 1)[-1],
                               #oo5$get_par("beta", m = 1)[-1],
                               #oo6$get_par("beta", m = 1)[-1],
                               c(rep((1 - true_par$beta) / (scrdat$n_occasions() - 1), scrdat$n_occasions() - 1))
                               ),
                      primary = rep(2:10, 5#7
                                    ),
                      model = c(rep("linear", 9),
                                rep("linear_wfirstprim", 9),
                                rep("smooth3df", 9),
                                rep("smooth3df_wfirstprim", 9),
                                #rep("smooth6df", 9),
                                #rep("smooth6df_wfirstprim", 9),
                                rep("true beta", 9)
                                ))

ggplot() +
  geom_line(data = plotdat, aes(x = primary, y = beta, col = model))

