#' Simulate Spatial Capture-Recapture data
#'
#' @param par named vector of D, lambda0, sigma parameters and sd, if movement is desired 
#' @param n_occasions number of occasions to simulate 
#' @param detectors detectors object made using SECR package make.traps 
#' @param mesh mesh object made using SECR make.mask
#' @param ihp relative density at each mesh point for inhomogeneous (so density at point = D*ihp)
#' @param move if TRUE, then activity centres move by Brownian motion and par$sd must be specified
#' @param time time of each occasions (e.g. 1 to 10, or 2002, 2003, 2004)
#' @param seed if supplied, random number generator is given this seed
#' @param print if TRUE then useful output is printed 
#'
#' @return a ScrData object 
#' @export
simulate_scr <- function(par, n_occasions, detectors, mesh, ihp = NULL, move = FALSE, time = NULL, seed = NULL, print = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(time)) time <- 1:n_occasions
  num.meshpts <- nrow(mesh)
  box <- attr(mesh, "boundingbox")
  D <- par$D
  if (!is.null(ihp)) {
    D <- D * ihp 
    model2D <- "IHP"
  } else {
    model2D <- "poisson"
  }
  D <- D / 100
  # simulate population
  if (print) cat("Simulating population and activity centres.......")
  pop <- sim.popn(D = D, core = mesh, model2D = model2D, Ndist = "poisson", buffertype = "rect")
  if (print) cat("done\n")
  dt <- rep(1, n_occasions - 1)
  if (!is.null(time))  dt <- diff(time)
  # generate capture histories
  lambda0 <- par$lambda0
  sigma <- par$sigma
  # move activity centres 
  nocc <- n_occasions
  nsess <- 1
  popn <- pop 
  trapn <- detectors
  if (move) {
    if (print) cat("Simulating moving activity centres......")
    poplist <- vector(mode = "list", length = n_occasions)
    poplist[[1]] <- pop
    traplist <- vector(mode = "list", length = n_occasions)
    traplist[[1]] <- detectors
    if (is.null(usage(detectors))) usage(detectors) <- matrix(1, nr = nrow(detectors), nc = n_occasions)
    usage(traplist[[1]]) <- matrix(usage(detectors)[,1], nr = nrow(detectors), nc = 1)
    for (k in 2:n_occasions) {
      for (i in 1:nrow(pop)) {
        dist <- pop[i,1] - mesh[,1]
        pr <- dnorm(dist, 0, par$sd * sqrt(dt[k-1]))
        dist2 <- pop[i,2] - mesh[,2]
        pr2 <- dnorm(dist2, 0, par$sd * sqrt(dt[k-1]))
        pr <- pr*pr2 / sum(pr*pr2)
        pop[i,] <- mesh[sample(1:nrow(mesh), size = 1, prob = pr),]    
      }
      poplist[[k]] <- pop
      traplist[[k]] <- detectors 
      usage(traplist[[k]]) <- matrix(usage(detectors)[,k], nr = nrow(detectors), nc = 1)
    }
    popn <- poplist
    trapn <- traplist
    nocc <- 1 
    nsess <- n_occasions
    cat("done\n")
  } 
  if (print) cat("Simulating capture histories......")
  if (move) {
    capturehistories <- vector(mode = "list", length = nsess)
    for (k in 1:nsess) {
      capturehistories[[k]] <- sim.capthist(trapn[[k]],  
                                            popn = popn[[k]], 
                                            detectfn = "HHN", 
                                            detectpar = list(lambda0 = lambda0, 
                                                             sigma = sigma), 
                                            noccasions = nocc, 
                                            renumber = FALSE)
    }
    capture_history <- join(capturehistories)
  } else {
    capture_history <- sim.capthist(trapn,  
                                    popn = popn, 
                                    detectfn = "HHN", 
                                    detectpar = list(lambda0 = lambda0, 
                                                     sigma = sigma), 
                                    noccasions = nocc,
                                    nsession = nsess, 
                                    renumber = FALSE)
  }
  if (print) cat("done\n")
  if (print) cat("Creating ScrData object.......")
  simdat <- ScrData$new(capture_history, mesh, time)
  if (print) cat("done\n")
  return(simdat)
}


#' Simulate SCR open population Cormack-Jolly-Seber survey
#'
#' @param par true parameters, named list of lambda0, sigma, phi
#' @param N number of individuals 
#' @param n_occasions total number of occasions in survey 
#' @param detectors secr trap object
#' @param mesh secr mesh object
#' @param move if TRUE then activity centres move and true_par$sd must be specified
#' @param time vector with time units between occasions 
#' @param primary index of primary periods that each occasion belongs to (default is 1 primary period)
#' @param seed seed to set before simulating 
#' @param print if TRUE then useful output is printed 
#'
#' @return ScrData object 
#' @export
simulate_cjs_openscr <- function(par, N, n_occasions, detectors, mesh,  move = FALSE, time = NULL, primary = NULL, seed = NULL, print = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(time)) {
    if (is.null(primary)) {
      time <- 1:n_occasions
    } else {
      time <- 1:max(primary)
    }
  }
  num_meshpts <- nrow(mesh)
  phi <- par$phi
  if (length(phi) == 1) phi <- rep(phi, n_occasions - 1)
  # simulate population
  area <- nrow(mesh) * attr(mesh, "area")
  D <- 1.5 * N / area
  pop <- matrix(1, nr = 1, nc = 1)
  if (print) cat("Simulating activity centres.......")
  while (nrow(pop) < N) pop <- sim.popn(D = D, core = mesh, Ndist = "fixed", buffertype = "convex")
  pop <- pop[1:N,]
  covariates(pop) <- covariates(pop)[1:N,]
  if (print) cat("done\n")
  life <- matrix(0, nr = nrow(pop), ncol = n_occasions) 
  life[, 1] <- 1
  diffprim <- diff(primary)
  if (all(diffprim == 0)) {
    diffprim <- rep(1, n_occasions - 1)
  } else {
    diffprim <- diff(primary)
  }
  for (k in 2:n_occasions) {
    if(diffprim[k - 1] > 0.5) {
      life[,k] <- rbinom(N, 1, phi[k - 1]) * life[, k - 1]
    } else {
      life[,k] <- life[,k-1]
    }
  }
  dt <- rep(1, n_occasions - 1)
  if (!is.null(time))  dt <- diff(time)
  lambda0 <- par$lambda0
  sigma <- par$sigma
  nocc <- n_occasions
  nsess <- 1
  popn <- pop 
  trapn <- detectors
  if (move) {
    if (print) cat("Simulating moving activity centres......")
    poplist <- vector(mode = "list", length = n_occasions)
    poplist[[1]] <- pop
    traplist <- vector(mode = "list", length = n_occasions)
    traplist[[1]] <- detectors
    if (is.null(usage(detectors))) usage(detectors) <- matrix(1, nr = nrow(detectors), nc = n_occasions)
    usage(traplist[[1]]) <- matrix(usage(detectors)[,1], nr = nrow(detectors), nc = 1)
    for (k in 2:n_occasions) {
      for (i in 1:nrow(pop)) {
        dist <- pop[i,1] - mesh[,1]
        pr <- dnorm(dist, 0, par$sd * sqrt(dt[k-1]))
        dist2 <- pop[i,2] - mesh[,2]
        pr2 <- dnorm(dist2, 0, par$sd * sqrt(dt[k-1]))
        pr <- pr*pr2 / sum(pr*pr2)
        pop[i,] <- mesh[sample(1:nrow(mesh), size = 1, prob = pr),]    
      }
      poplist[[k]] <- pop
      traplist[[k]] <- detectors 
      usage(traplist[[k]]) <- matrix(usage(detectors)[,k], nr = nrow(detectors), nc = 1)
    }
    popn <- poplist
    trapn <- traplist
    nocc <- 1 
    nsess <- n_occasions
    if (print) cat("done\n")
  } 
  if (print) cat("Simulating capture histories......")
  if (move) {
    capturehistories <- vector(mode = "list", length = nsess)
    for (k in 1:nsess) {
      capturehistories[[k]] <- sim.capthist(trapn[[k]],  
                                            popn = popn[[k]], 
                                            detectfn = "HHN", 
                                            detectpar = list(lambda0 = lambda0, 
                                                             sigma = sigma), 
                                            noccasions = nocc, 
                                            renumber = FALSE)
    }
    capture_history <- join(capturehistories)
  } else {
    capture_history <- sim.capthist(trapn,  
                                    popn = popn, 
                                    detectfn = "HHN", 
                                    detectpar = list(lambda0 = lambda0, 
                                                     sigma = sigma), 
                                    noccasions = nocc,
                                    nsession = nsess, 
                                    renumber = FALSE)
  }
  if (print) cat("done\n")
  # thin capture history by survival
  if (print) cat("Simulating life histories.......")
  n <- dim(capture_history)[1]
  ids <- as.numeric(rownames(capture_history))
  life <- life[ids,]
  seen <- rep(TRUE, n)
  full_cap <- capture_history 
  for (i in seq(n)) {
    life[i, 1:n_occasions] <- cumprod(life[i, 1:n_occasions])
    capture_history[i, ,] <- diag(life[i,]) %*% capture_history[i, ,]
    if (sum(capture_history[i, ,]) == 0) seen[i] <- FALSE
  }
  if (!is.null(attr(capture_history, "detectedXY"))) { 
    xy <- attr(capture_history, "detectedXY")
    inc <- rep(FALSE, nrow(xy))
    r <- 1 
    for (i in 1:length(full_cap)) {
      if (capture_history[i] > 0) inc[seq(r, r + capture_history[i] - 1)] <- TRUE
      r <- r + full_cap[i]
    }  
    xy <- xy[inc,]
    attributes(capture_history)$detectedXY <- xy 
  }
  capture_history <- subset(capture_history, subset = (1:n)[seen])
  A <- nrow(mesh) * attr(mesh, "area")
  if (print) cat("done\n")
  if (print) cat("Creating ScrData object........")
  simdat <- ScrData$new(capture_history, mesh, time, primary = primary) 
  if (print) cat("done\n")
  return(simdat)
}

#' Simulate SCR open population Jolly-Seber survey 
#'
#' @param par true parameters, named list of lambda0, sigma, phi, beta
#' @param n_occasions total number of primary occasions in survey 
#' @param n_sec_occasions number of secondary occasions in survey
#' @param detectors secr trap object
#' @param mesh secr mesh object
#' @param ihp relative density at each mesh point for inhomogeneous (so density at point = D*ihp)
#' @param move if TRUE, activity centres move and par$sd must be specified 
#' @param time vector with time units between occasions 
#' @param primary index of which primary each occasion is allocated to 
#' @param ne_trans user created transition layer for noneuclidean distances
#' @param seed seed to set before simulating 
#' @param print if TRUE, useful output is printed 
#'
#' @return ScrData object 
#' @export
simulate_js_openscr <- function(par, n_occasions, n_sec_occasions, detectors, mesh, ihp = NULL, move = FALSE, time = NULL, primary = NULL, ne_trans = NULL, seed = NULL, print = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(time)) {
    if (is.null(primary)) {
      time <- 1:n_occasions
    } else {
      time <- 1:max(primary)
    }
  }
  if(!is.null(ne_trans)){
    covariates(mesh) <- data.frame(noneuc = factor(rep("noneuc", nrow(mesh))))
    userdfn1 <- function (xy1, xy2, mesh) {
      if (missing(xy1)) return("noneuc")
      require(gdistance)
      costDistance(ne_trans, as.matrix(xy1), as.matrix(xy2))
    }
  } else {
    userdfn1 <- NULL
  }
  #note states only implemented for detection parameters
  boolstatemod <- length(par$lambda0) > 1 | length(par$sigma) > 1
  if (boolstatemod){
      warning("Mixture models only implemented for detection parameters")
  }
  num_meshpts <- nrow(mesh)
  D <- par$D #density is hrc per km^2
  if (!is.null(ihp)) {
    D <- D * ihp 
    model2D <- "IHP"
  } else {
    model2D <- "poisson"
  }
  D <- D / 100 #convert from hrc per km^2 to hrc per hectare
  phi <- par$phi
  if (length(phi) == 1) phi <- rep(phi, n_occasions - 1)
  beta <- par$beta
  if (length(beta) == 1) beta <- c(beta, rep((1 - beta) / (n_occasions - 1), n_occasions - 1))
  # simulate population
  if (print) cat("Simulating population and activity centres.......")
  pop <- sim.popn(D = D, core = mesh, model2D = model2D, 
                  Ndist = "poisson", buffertype = "rect")
  rownames(pop) <- NULL
  if (print) cat("done\n")
  birth_time <- sample(1:n_occasions, size = nrow(pop), prob = beta, replace = TRUE)
  dt <- diff(time)
  life <- matrix(0, nr = nrow(pop), ncol = n_occasions) 
  life[birth_time == 1, 1] <- 1
  for (k in 2:n_occasions){
    alives <- life[,k-1] == 1
    life[alives, k] <- rbinom(sum(alives), 1, phi[k - 1]^dt[k-1]) 
    life[birth_time == k, k] <- 1
  }
  
  lambda0 <- par$lambda0
  sigma <- par$sigma
  nocc <- n_sec_occasions
  nsess <- 1
  popn <- pop 
  trapn <- detectors
  if (move) {
    #note mixture model not yet implemented for moving hrcs
    if(!is.null(ne_trans)){
      if (print) cat("Calculating noneuclidean distance matrix......")
      meshbymesh <- userdfn1(mesh, mesh, mesh)
      trapoffset <- max(apply(as.array(1:nrow(detectors)), 1,
                              FUN = function(x){
                                max(c(min(abs(mesh$x - detectors[x,"x"])),
                                      min(abs(mesh$y - detectors[x,"y"]))))
                              }))
      if(trapoffset > 1) {
        warning(paste("Traps are offset ", trapoffset, "m from mesh pts", sep = ""))
      }
      trapismesh <- apply(as.array(1:nrow(detectors)), 1, 
                          FUN = function(x){which(
                            abs(mesh$x - detectors[x,"x"]) <= trapoffset  & 
                              abs(mesh$y - detectors[x,"y"]) <= trapoffset )})
      if (print) cat("done\n")
      
      userdistn <- userdfn1(mesh, pop, mesh)
      userdistlist <- vector(mode = "list", length = n_occasions)
      userdistlist[[1]] <- userdistn
    } else {
      userdistn <- NULL
      trapismesh <- NULL
    }
    if (print) cat("Simulating moving activity centres......")
    poplist <- vector(mode = "list", length = n_occasions)
    poplist[[1]] <- pop
    traplist <- vector(mode = "list", length = n_occasions)
    traplist[[1]] <- detectors
    if (is.null(usage(detectors))){ usage(detectors) <- matrix(1, nr = nrow(detectors), nc = nocc)}
    usecols <- which(primary == 1)
    usage(traplist[[1]]) <- matrix(usage(detectors)[,usecols], nr = nrow(detectors), nc = length(usecols))
    for (k in 2:n_occasions) {
      if(!is.null(ne_trans)){
        userdistlist[[k]] <- userdistlist[[k-1]]
      }
      for (i in 1:nrow(pop)) {
        #note this is euclidean distance
        if(is.null(userdistlist)){
          dist <- pop[i,1] - mesh[,1]
          dist2 <- pop[i,2] - mesh[,2]
          pr <- dnorm(dist, 0, par$sd * sqrt(dt[k-1]))
          pr2 <- dnorm(dist2, 0, par$sd * sqrt(dt[k-1]))
          pr <- pr*pr2 / sum(pr*pr2)
        } else {
          dist <- userdistlist[[k]][,i]
          pr <- dnorm(dist, 0, par$sd * sqrt(dt[k-1]))
          pr <- pr/sum(pr)
        }
        if (!all(is.na(pr))){
          #if pop i is down too narrow of a stream, don't move it
          whichmesh <- sample(1:nrow(mesh), size = 1, prob = pr)
          pop[i,] <- mesh[whichmesh,]   
          userdistlist[[k]][,i] <- meshbymesh[whichmesh,]
        }
      }
      poplist[[k]] <- pop
      traplist[[k]] <- detectors 
      usecols <- which(primary == k)
      usage(traplist[[k]]) <- matrix(usage(detectors)[,usecols], nr = nrow(detectors), nc = length(usecols))
    }
    popn <- poplist
    trapn <- traplist
    nsess <- n_occasions
    if(!is.null(ne_trans)){
      userdistn <- userdistlist
    }
    if (print) cat("done\n")
  } 
  if (print) cat("Simulating capture histories......")
  if (move) {
    capturehistories <- vector(mode = "list", length = nsess)
    for (k in 1:nsess) {
      capturehistories[[k]] <- sim.capthist(trapn[[k]],  
                                            popn = popn[[k]], 
                                            detectfn = "HHN", 
                                            detectpar = list(lambda0 = lambda0, 
                                                             sigma = sigma), 
                                            userdist = userdistn[[k]][trapismesh,],
                                            noccasions = length(which(primary == k)), 
                                            renumber = FALSE)
    }
    capture_history <- join(capturehistories)
  } else {
    if(boolstatemod){
      #assign state membership (NOTE only valid for lambda and sigma mixtures)
      nstates <- max(length(lambda0), length(sigma))
      stateprobs <- par$delta
      lambda0_state <- lambda0
      if (length(lambda0_state) == 1){
        lambda0_state <- rep(lambda0_state, nstates)
      }
      sigma_state <- sigma
      if (length(sigma_state) == 1){
        sigma_state <- rep(sigma_state, nstates)
      }
      istate <- sample(1:nstates, size = nrow(popn), prob = stateprobs, replace = TRUE)
      ch_list <- list()
      capture_history <- c()
      for (state in 1:nstates){
        ch_list[[state]] <- sim.capthist(trapn,  
                                         popn = popn[istate == state,], 
                                         detectfn = "HHN", 
                                         detectpar = list(lambda0 = lambda0_state[state], 
                                                          sigma = sigma_state[state]), 
                                         userdist = userdfn1,
                                         noccasions = n_sec_occasions,
                                         nsession = nsess, 
                                         renumber = FALSE)
      }
      capture_history <- do.call(rbind, ch_list)
    } else {
    capture_history <- sim.capthist(trapn,  
                                    popn = popn, 
                                    detectfn = "HHN", 
                                    detectpar = list(lambda0 = lambda0, 
                                                     sigma = sigma), 
                                    userdist = userdfn1,
                                    noccasions = n_sec_occasions,
                                    nsession = nsess, 
                                    renumber = FALSE)
    }
  }
  if (print) cat("done\n")
  # thin capture history by survival
  if (print) cat("Simulating life histories.......")
  n <- dim(capture_history)[1]
  ids <- as.numeric(rownames(capture_history))
  life <- life[ids,]
  birth_time <- birth_time[ids]
  seen <- rep(TRUE, n)
  for (i in seq(n)) {
    #life[i, birth_time[i]:n_occasions] <- cumprod(life[i, birth_time[i]:n_occasions])
    capture_history[i, life[i,][primary] ==0,] <- 0
    #capture_history[i, ,] <- diag(life[i,][primary]) %*% capture_history[i, ,]
    if (sum(capture_history[i, ,]) == 0) seen[i] <- FALSE
  } 
  if (!is.null(attr(capture_history, "detectedXY"))) { 
    xy <- attr(capture_history, "detectedXY")
    inc <- rep(FALSE, nrow(xy))
    r <- 1 
    for (i in 1:length(full_cap)) {
      if (capture_history[i] > 0) inc[seq(r, r + capture_history[i] - 1)] <- TRUE
      r <- r + full_cap[i]
    }  
    xy <- xy[inc,]
    attributes(capture_history)$detectedXY <- xy 
  }
  capture_history <- subset(capture_history, (1:n)[seen])
  A <- nrow(mesh) * attr(mesh, "area")
  if (print) cat("done\n")
  if (print) cat("Creating ScrData object......")
  simdat <- ScrData$new(capture_history, mesh, time, primary = primary) 
  if (print) cat("done\n")
  return(simdat)
}



