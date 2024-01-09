#' Noneuclidean distances for simulated population and 
#' 
#' @param trans_poly SpatialPolygons object that defines water
#' @param popgrid dataframe with x and y columns that is possible hrc grid  
#'  
#' @return transition layer
#' @export     
pop_netrans <- function(trans_poly, pop_grid){
  pop_pts <- st_as_sf(x = pop_grid, coords = c("x","y"), crs = crs(trans_poly))
  #then remaining traps
  connects <- nngeo::st_connect(poppts, st_as_sf(trans_poly))
  connects <- connects[ which(as.numeric(st_length(connects)) > 0.001)]
  fs <- seq(1:length(connects))*2 -1
  boxify <- function(f, connects, poly) {
    first <- st_coordinates(connects)[f,1:2]
    last <- st_coordinates(connects)[(f+1),1:2]
    #first make sure no pts are equal
    if(any(first == last)){
      print(paste("Pts directly in a line or identical ", f))
    }
    slope <- (first[2]-last[2])/(first[1]-last[1])
    invslope <- -1/slope
    b <- first[2] - slope*first[1]
    offset = 100
    
    #if first pt is to the left of last pt
    if (first[1] < last[1]){
      left <- first
      right <- last
    } else {
      left <- last
      right <- first
    }
    #finds two points that are offset away from pt on a line of slope 
    pt_away <- function(pt, offset, slope){
      newpt1 <- c(X = pt[1] + sqrt((offset^2)/(slope^2 + 1)), 
                  Y = sum(slope * sqrt((offset^2)/(slope^2 +1)), pt[2], na.rm = T))
      newpt2 <- c(X = pt[1] - sqrt((offset^2)/(slope^2 + 1)), 
                  Y = sum(-(slope * sqrt((offset^2)/(slope^2 +1))), pt[2], na.rm = T))
      rbind(newpt1, newpt2)
    }
    left_2 <- pt_away(left, offset, slope)
    #if identical (vertical line)
    if (all(left_2[1,] == left_2[2,])){
      #just double check vertical line
      if(left[1] == right[1]){
        bottom <- c(X = left[1], Y = min(left[2],right[2]))
        bottom_2 <- c(X = bottom[1], Y = bottom[2] - offset)
        Ss <- pt_away(bottom_2, offset, invslope)
        SW <- Ss[which(Ss[,1] == min(Ss[,1])),]
        SE <- Ss[which(Ss[,1] == max(Ss[,1])),]
        top <- c(X = left[1], Y = max(left[2], right[2]))
        top_2 <- c(X = top[1], Y = top[2] + offset)
        Ns <- pt_away(top_2, offset, invslope)
        NW <- Ns[which(Ns[,1] == min(Ns[,1])),]
        NE <- Ns[which(Ns[,1] == max(Ns[,1])),]
      }
    } else {
      #keep whichever is smaller x (more left)
      left_2 <- left_2[which(left_2[,1] == min(left_2[,1])),]
      Ws <- pt_away(left_2, offset, invslope)
      NW <- Ws[which(Ws[,2] == max(Ws[,2])),]
      SW <- Ws[which(Ws[,2] == min(Ws[,2])),]
      right_2 <- pt_away(right, offset, slope)
      #keep whichever is larger x (more right)
      right_2 <- right_2[which(right_2[,1] == max(right_2[,1])),]
      Es <- pt_away(right_2, offset, invslope)
      NE <- Es[which(Es[,2] == max(Es[,2])),]
      SE <- Es[which(Es[,2] == min(Es[,2])),]
    }
    box<- as_Spatial(st_sfc(st_polygon(list(rbind(
      NW,
      NE,
      SE,
      SW,
      NW
    ))),
    crs = st_crs(st_as_sf(poly))))
    return(box)
  }
  conn_polys_ls <- lapply(X = as.list(fs), FUN = boxify,
                          connects = connects, poly = trans_poly)
  connect_polys <- do.call(bind, conn_polys_ls)
  out_poly <- bind(connect_polys, new_poly)
  r <- raster(ncol = 1000, nrow = 1000)
  extent(r) <- extent(out_poly)
  rp <- rasterize(out_poly, r)
  values(rp)[!is.na(values(rp))] = 1
  rp_df <- as.data.frame(as(rp, "SpatialPixelsDataFrame"))
  colnames(rp_df) <- c("value", "x", "y")
  trans <- transition(rp, mean, directions = 16)
  return(trans)
}
  
#' #' Get distance
#' #' 
#' #' @param m row of mesh
#' #' @param j row of traps
#' #' @param trans transition layer
#' #' @param poly SpatialPolygon to reference crs
#' #' @param mesh mesh data.frame
#' #' @param traps traps data.frame
#' try_dist <- function(m, j, trans, poly, mesh, traps){
#'   get_dist <- function(m, j, trans, poly, mesh, traps){
#'     if (dist(rbind(mesh[m,1:2], traps[j,])) < 2) {
#'       distance = 0
#'     } else {
#'       path <- shortestPath(x = trans, 
#'                            origin = as.matrix(mesh.[m, c("x","y")]), 
#'                            goal = as.matrix(traps.[j, c("x","y")]),
#'                            output = "SpatialLines")
#'       crs(path) <- crs(poly)
#'       distance <- gLength(path)
#'     }
#'     return(distance)
#'   }
#'   out <- tryCatch(
#'     {
#'       get_dist(m, j, trans, poly, mesh, traps)
#'     },
#'     error=function(e){
#'       return(NA)
#'     },
#'     warning=function(w){
#'       return(NA)
#'     }
#'   )
#'   return(out)
#' }

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
#' @param poly bounding polygon for secr sim.popn
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
simulate_js_openscr <- function(par, n_occasions, n_sec_occasions, detectors, mesh, poly = NULL, ihp = NULL, move = FALSE, time = NULL, primary = NULL, ne_trans = NULL, seed = NULL, print = TRUE) {
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
  pop <- sim.popn(D = D, core = mesh, poly = poly, model2D = model2D, Ndist = "poisson", buffertype = "rect")
  rownames(pop) <- NULL
  if (print) cat("done\n")
  birth_time <- sample(1:n_occasions, size = nrow(pop), prob = beta, replace = TRUE)
  dt <- diff(time)
  life <- matrix(0, nr = nrow(pop), ncol = n_occasions) 
  life[birth_time == 1, 1] <- 1
  for (k in 2:n_occasions){
    alives <- life[,k-1] == 1
    life[alives, k] <- rbinom(sum(alives), 1, phi[k - 1]) 
    life[birth_time == k, k] <- 1
  }
  
  lambda0 <- par$lambda0
  sigma <- par$sigma
  nocc <- n_sec_occasions
  nsess <- 1
  popn <- pop 
  trapn <- detectors
  if(!is.null(ne_trans)){
    userdistn <- userdfn1(mesh, pop, mesh)
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
    
  } else {
    userdistn <- NULL
    trapismesh <- NULL
  }
    
  if (move) {
    if (print) cat("Simulating moving activity centres......")
    poplist <- vector(mode = "list", length = n_occasions)
    poplist[[1]] <- pop
    traplist <- vector(mode = "list", length = n_occasions)
    traplist[[1]] <- detectors
    if(!is.null(ne_trans)){
      userdistlist <- vector(mode = "list", length = n_occasions)
      userdistlist[[1]] <- userdistn
    }
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
    capture_history <- sim.capthist(trapn,  
                                    popn = popn, 
                                    detectfn = "HHN", 
                                    detectpar = list(lambda0 = lambda0, 
                                                     sigma = sigma), 
                                    userdist = userdistn[trapismesh,],
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
  birth_time <- birth_time[ids]
  seen <- rep(TRUE, n)
  for (i in seq(n)) {
    life[i, birth_time[i]:n_occasions] <- cumprod(life[i, birth_time[i]:n_occasions])
    capture_history[i, ,] <- diag(life[i,][primary]) %*% capture_history[i, ,]
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



