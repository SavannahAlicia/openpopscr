# Copyright (c) 2018 Richard Glennie
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#' SCR with transience model class 
#' 
#' @description Spatial capture-recapture model with transience: fits model, formats inference, and 
#' simulates from fitted model. This model inherits all the functions of the 
#' JsModel class. Functions that are different or additional are documented here. 
#' \itemize{
#'   \item form: a named list of formulae for each parameter (~1 for constant)
#'   \item scr_data: a ScrData object 
#'   \item start: a named list of starting values 
#'   \item print (default = TRUE): if TRUE then helpful output is printed
#' }
#' 
#' Methods include all those in ScrModel with following re-defined: 
#' \itemize{
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead)
#'  \item calc_D_llk(): computes the likelihood of the D parameter
#'   \item calc_Dpdet(): compute the integral int D(x)p(x) dx <- the overall intensity of detected individuals in the 
#'        survey area
#'  \item calc_pdet(): compute probability of being detected at least once during the survey
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#' }
#' 
ScrTransientModel <- R6Class("ScrTransientModel",
                             inherit = ScrModel, 
  public = list(
    
    initialize = function(form, data, start, detectfn = NULL, statemod = NULL, print = TRUE) {
      private$check_input(form, data, start, detectfn, print)
      
      # If user distance matrix, also need mesh distance matrix for moving AC
      if(!is.null(data$userdistmat())){
        if(is.null(data$usermeshdistmat())){
          stop("For moving AC, you must include a non-Euclidean distance matrix for the mesh as well as traps to mesh.")
        } else {
          private$noneuclidean_ <- TRUE
        }
      } else {
        private$noneuclidean_ <- FALSE
      }
      private$data_ <- data
      private$start_ <- start 
      private$dx_ <- attr(data$mesh(), "spacing")
      
      if(private$noneuclidean_){
        closest_dx <- max(apply(private$data_$usermeshdistmat(), 1, function(mdr){min(mdr[mdr!=0])}))
        # Need to check that there aren't pts that can't trans to all others
        if(length(which(apply(private$data_$usermeshdistmat(), 1, function(mdr){sum(mdr <= closest_dx)}) ==2)) > 1){
          closest_dx <- max(apply(private$data_$usermeshdistmat(), 1, function(mdr){min(mdr[mdr > closest_dx])}))
        } # Non Euclidean distance requires different "inside" indexing
        ncol_inside <- max(apply(private$data_$usermeshdistmat(),1, function(mdr){sum(mdr<= closest_dx)}))
      } else {
        closest_dx <- (1 + 1e-6) * private$dx_
        ncol_inside <- 4 #at most 4 cells can be adjacent
      }
      
      if(private$noneuclidean_){
        private$meshdistmat_ <- private$data_$usermeshdistmat()
      } else {
        private$meshdistmat_ <- apply(private$data_$mesh(), 1, function(meshx){
          sqrt((meshx[1] - private$data_$mesh()$x)^2 + 
                 (meshx[2]- private$data_$mesh()$y)^2)})
      }
      
      private$inside_ <- matrix(-1, nr = data$n_meshpts(), nc = ncol_inside)
      for (m in 1:data$n_meshpts()) {
        dis <- private$meshdistmat_[m,]
        # Which mesh pts are adjacent(e) or within radius such that all pts are connected (ne)
        wh <- which(dis <= closest_dx & dis > 1e-16) - 1 #zero based indexing for Rcpp
        if(length(wh) > 0) private$inside_[m, 1:length(wh)] <- as.numeric(wh) 
      }
      box <- attributes(data$mesh())$boundingbox
      region <- c(diff(box[1:2, 1]), diff(box[c(1, 3), 2]))
      private$num_cells_ <- numeric(3)
      private$num_cells_[1] <- data$n_meshpts()
      private$num_cells_[2] <- floor(region[1] / private$dx_)
      private$num_cells_[3] <- data$n_meshpts() / private$num_cells_[2]
      if (print) cat("Reading formulae.......")
      order <- c("sd", "D")
      private$read_formula(form, detectfn, statemod, order)
      private$par_type_[private$detfn_$npars() + 1] <- "k1ms"
      private$par_type_[private$detfn_$npars() + 2] <- "m"
      names(private$form_) <- c(private$detfn_$pars(), "sd", "D")
      # make parameter list 
      private$make_par() 
      # set link functions
      private$link2response_ <- c(private$detfn_$link2response(), list("exp"), list("exp"))
      names(private$link2response_) <- c(private$detfn_$pars(), "sd", "D")
      private$response2link_ <- c(private$detfn_$response2link(), list("log"), list("log"))
      names(private$response2link_) <- c(private$detfn_$pars(), "sd", "D")
      if (print) cat("done\n")
      if (print) cat("Initialising parameters.......")
      private$initialise_par(start)
      private$read_states() 
      if (print) cat("done\n")
      private$print_ = print 
    },
    
    # calc_initial_distribution = function() {
    #   n_mesh <- private$data_$n_meshpts()
    #   nstates <- private$state_$nstates()
    #   delta <- private$state_$delta() 
    #   pr0 <- matrix(delta, nrow = n_mesh, ncol = nstates, byrow = TRUE)
    #   a <- private$data_$cell_area() 
    #   D <- self$get_par("D", m = 1:n_mesh) * a 
    #   pr0 <- pr0 * D
    #   return(pr0)
    # },
    
    calc_Dpdet = function() {
      # compute probability of zero capture history 
      enc_rate <- self$encrate() 
      nstates <- self$state()$nstates()
      trap_usage <- usage(private$data_$traps())
      pr_empty <- list()
      for (j in 1:private$data_$n_occasions()) {
        pr_empty[[j]] <- matrix(1, nr = private$data_$n_meshpts(), nc = nstates)
        for (g in 1:nstates) {
          pr_empty[[j]][, g] <- exp(-t(trap_usage[, j]) %*% t(enc_rate[[g]][,,j]) )
        }
      }
      pr0 <- self$calc_initial_distribution()
      tpms <- self$calc_tpms()
      dt <- diff(self$data()$time())
      sd <- as.matrix(self$get_par("sd", s = 1:self$state()$nstates(), m = 1))
      sd[is.na(sd)] <- -10
      Dpdet <- C_calc_move_pdet(private$data_$n_occasions(), 
                               pr0, 
                               pr_empty, 
                               tpms, 
                               private$num_cells_,
                               private$inside_, 
                               private$dx, # was private$dx_, #change
                               dt, 
                               sd,
                               nstates, 
                               0, 
                               0,
                               private$noneuclidean_,
                               private$meshdistmat_); 
      a <- private$data_$cell_area() 
      D <- self$get_par("D", m = 1:private$data_$n_meshpts()) * a
      Dpdet <- sum(D) - Dpdet 
      return(Dpdet)
    },
    
    calc_llk = function(param = NULL, names = NULL) {
      if (!is.null(names)) names(param) <- names 
      if (!is.null(param)) {
        slen <- length(self$state()$par())
        param2 <- param 
        if (slen > 0) {
          ind <- seq(length(param) - slen + 1, length(param))
          self$state()$set_par(param[ind])
          param2 <- param[-ind]
        }
        self$set_par(private$convert_vec2par(param2));
      }
      # initial distribution 
      pr0 <- self$calc_initial_distribution()
      # compute probability of capture histories 
      # across all individuals, occasions and traps 
      pr_capture <- self$calc_pr_capture()
      n <- private$data_$n()
      n_occasions <- private$data_$n_occasions()
      n_meshpts <- private$data_$n_meshpts() 
      # get tpms for state model 
      nstates <- self$state()$nstates() 
      tpms <- self$calc_tpms()
      # get movement 
      dt <- diff(self$data()$time())
      sd <- as.matrix(self$get_par("sd", s = 1:self$state()$nstates(), m = 1))
      sd[is.na(sd)] <- -10
      # compute likelihood for each individual
      llk <- C_calc_move_llk(n, 
                             n_occasions,
                             pr0, 
                             pr_capture, 
                             tpms,
                             private$num_cells_, 
                             private$inside_, 
                             private$dx_, # was private$dx_, #change
                             dt, 
                             sd, 
                             nstates, 
                             0, 
                             0,
                             private$noneuclidean_,
                             private$meshdistmat_,
                             rep(0, private$data_$n()))
      # compute log-likelihood
      llk <- llk - n * log(self$calc_Dpdet())
      llk <- llk + self$calc_D_llk()
      if (private$print_) cat("llk:", llk, "\n")
      return(llk)
    }
  
),
                   
  private = list(
    dx_ = NULL, 
    inside_ = NULL,
    num_cells_ = NULL,
    noneuclidean_ = NULL,
    meshdistmat_ = NULL,
    
    initialise_par = function(start) {
      n_det_par <- private$detfn_$npars()
      names <- private$detfn_$pars()
      for (i in 1:n_det_par) {
        private$par_[[names[i]]][1] <- do.call(private$response2link_[[names[i]]], 
                                               list(start[[names[i]]]))
      }
      names(private$par_) <- c(names, "sd", "D")
      private$par_$sd[1] <-do.call(private$response2link_$sd, 
                                   list(start$sd))
      private$par_$D[1] <- do.call(private$response2link_$D, 
                                list(start$D))
      # compute initial parameters for each jkm
      private$compute_par()
      return(invisible())
    }
  )                 
)



