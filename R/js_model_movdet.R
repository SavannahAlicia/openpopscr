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
################################################################################

#' Jolly-Seber moving-detector model class 
#' 
#' @description Jolly-Seber moving detector model: fits model, formats inference, and 
#' simulates from fitted model. 
#' \itemize{
#'   \item form: a named list of formulae for each parameter (~1 for constant)
#'   \item data: a ScrData object for moving detectors
#'   \item start: a named list of starting values 
#'   \item print (default = TRUE): if TRUE then helpful output is printed to the screen
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item get_par(name, j, k, m): returns value of parameter "name" for detector j 
#'   on occasion k (if j, k omitted then returns value(s) for all)
#'  \item set_par(par): can change the parameter the model uses. Note, the model will simulate 
#'    data using this parameter, but will only present inference based on the maximum likelihood
#'    estimates. 
#'  \item set_mle(mle, V, llk): set maximum likleihood for this model with parameters mle, 
#'  covariance matrix V, and maximum likelihood value llk
#'  \item calc_D_llk(): computes the likelihood of the D parameter
#'  \item calc_initial_distribution(): computes initial distribution over life states (unborn, alive, dead)
#'  \item calc_pr_entry(): computes vector with entry j equal to probability of individual unborn up to occasion j
#'  being born just after occasion j
#'  \item calc_tpms(): returns list of transition probability matrix for each occasion 
#'  \item calc_pr_capture(): returns array where (i,k,m) is probability of capture record 
#'  on occasion k for individual i given activity centre at mesh point m
#'  \item calc_pdet(): compute probability of being detected at least once during the survey
#'  \item calc_llk(): compute log-likelihood at current parameter values 
#'  \item fit(): fit the model by obtaining the maximum likelihood estimates. Estimates of
#'        density are obtained from parametric boostrap with nsim resamples. 
#'  \item par(): return current parameter of the model 
#'  \item mle(): return maximum likelihood estimates for the fitted model 
#'  \item data(): return ScrData that the model is fit to 
#'  \item estimates(): return estimates in a easy to extract list 
#'  \item cov_matrix(): return variance-covariance matrix from fitted model (on working scale)
#'  \item mle_llk(): return log-likelihood value of maximum likelihood estimates 
#' }
#' 
JsModel <- R6Class("JsModelMovDet",
                   inherit = JsModel,
                   public = list(
                     
                     calc_pr_capture = function() {
                       n_occasions <- private$data_$n_occasions("all")
                       n_primary <- private$data_$n_primary()
                       nstates <- self$state()$nstates()
                       kstates <- private$known_states_
                       S <- private$data_$n_secondary() 
                       if (n_primary == 1) {
                         n_primary <- n_occasions
                         S <- rep(1, n_occasions)
                       }
                       enc_rate0 <- self$encrate()
                       trap_usage <- usage(private$data_$traps())
                       n <- private$data_$n()
                       n_meshpts <- private$data_$n_meshpts() 
                       n_traps <- private$data_$n_traps()
                       capthist <- private$data_$capthist()
                       imesh <- private$data_$imesh()
                       induse <- private$data_$induse()
                       prob <- C_calc_pr_capture_movdet(n, 
                                                 n_occasions, 
                                                 n_traps, 
                                                 n_meshpts, 
                                                 capthist, 
                                                 enc_rate0, 
                                                 trap_usage, 
                                                 induse,
                                                 nstates, 
                                                 1, 
                                                 1, 
                                                 kstates, 
                                                 self$data()$detector_type(),
                                                 n_primary, 
                                                 S,
                                                 rep(0, n), 
                                                 imesh, 
                                                 private$data_$capij())
                       return(prob)
                     }
                     
                   ),
                   
                   private = list(
                     Dk_ = NULL, 
                     var_ = NULL, 
                     confint_ = NULL
                     
                   )                 
)



