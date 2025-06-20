# Copyright (c) 2017 Richard Glennie
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
################################################################################
# openpopscr project: open population spatial capture-recapture 
#
# file description: scr data class based on secr package objects 
#
################################################################################

#' Spatial capture-recapture data class 
#' 
#' @description Encapsulates data from spatial capture-recapture survey into a single 
#' object. Object can be created using $new with arguments: 
#' \itemize{
#'   \item capthist: capture history object from secr package with trap information included
#'   \item mesh: mask object from secr package
#'   \item time: optional vector of numeric time each occasion took place at 
#'     (used for irregularly spaced occasions). Default is vector 1:n_occasions.
#'   \item primary: vector with index for each occasion in capthist that pools occasions into primary occasions
#'   \item userdistmat: user defined non
#'        
#' }
#' 
#' Methods include: 
#' \itemize{
#'  \item capthist(): return capture history object 
#'  \item traps(): return traps object 
#'  \item mesh(): return mesh object
#'  \item time(): return time vector 
#'  \item plot_mesh(cov, nbreaks = 10, ...): if cov not NULL then plots covariate value over mesh, nbreaks 
#'        is number of breaks in covariate value plotted, ... is passed to geom_point calls
#'  \item detector_type(): return detector type index (1 = count, 2 = proximity, 
#'        3 = multi/single)
#'  \item get_cov_list(): return list of covariates
#'  \item covs(j, k, m): return covariate values at detector j 
#'        on occasion k at mesh point m, if any of j,k are NULL then returns for all values; if m is NULL returns spatial mean 
#'  \item n(): return number of individuals seen over the entire survey
#'  \item n_occasions(): return number of primary occasions in the survey, n_occasions("all") returns #'        total number of secondary occasions, n_occasions(i) returns number of secondary occasions #'        in primary occasion i 
#'  \item n_traps(): return number of traps used at some time in the survey
#'  \item n_meshpts(): return number of mesh points 
#'  \item n_primary(): return number of primary occasions 
#'  \item n_secondary(): return vector of number of secondary occasions per primary
#'  \item primary(): return vector showing primary occasion each secondary occasion is assigned to
#'  \item encrate(each = FALSE): return empirical estimate of encounter rate, each = TRUE returns for each primary occasion
#'  \item encrange(k = NULL, each = FALSE): return empirical estimate of encounter range, each = TRUE returns for each primary occasion, k can be a vector/numeric specifying which occasions to compute over
#'  \item unique: return number of unique captures per occasion 
#'  \item area(): return total area of the mesh 
#'  \item cell_area(): return area of a single grid cell on mesh
#'  \item replace_mesh(newmesh): replace stored mesh with newmesh
#'  \item distances(): return matrix with (i,j)th entry the distance from trap i to mesh point j
#'  \item add_covariate(covariate_name, covariate_vector/matrix, covariate_type): add covariate to object
#'  \item remove_covariate(covariate_name): remove covariate from data object
#' }
#' 
ScrData <- R6Class("ScrDataMovDet", 
                   inherit = ScrModel,
  public = list( 
    initialize = function(capthist, 
                          mesh, 
                          time = NULL, 
                          primary = NULL,
                          userdm = NULL,
                          usermeshdm = NULL,
                          induse = NULL) {
      private$check_input(capthist, mesh, time, primary, userdm, usermeshdm, induse) 
      ## detectors
      private$detector_type_ <- switch(attr(traps(capthist), "detector")[1], 
                                       count = 1, 
                                       proximity = 2, 
                                       multi = 3,
                                       single = 3, 
                                       transect = 4, 
                                       transectX = 5, 
                                       polygon = 6, 
                                       polygonX = 7)
      if (attr(traps(capthist), "detector")[1] == "single") warning("Single-catch detectors treated as multi-catch detectors.")
      private$capthist_ <- capthist 
      if (is.null(secr::usage(secr::traps(capthist)))) {
        if (private$detector_type_ %in% 4:7) {
          nids <- length(unique(attr(traps(capthist), "polyID")))
          secr::usage(secr::traps(private$capthist_)) <- matrix(1, nr = nids, nc = dim(capthist)[2])
        } else {
          secr::usage(secr::traps(private$capthist_)) <- matrix(1, nr = dim(capthist)[3], nc = dim(capthist)[2])
        }
      }
      if (private$detector_type_ %in% 4:7) {
        outputdetector <- ifelse(private$detector_type_ %in% c(4, 6), "count", "proximity")
        spacing <- attr(mesh, "polygon_trap_spacing")
        if (is.null(spacing)) spacing <- attr(mesh, "spacing")
        private$capthist_ <- secr::discretize(capthist, 
                                              spacing = spacing, 
                                              outputdetector = outputdetector, 
                                              cell.overlap = TRUE)
        if (outputdetector == "proximity") attributes(traps(private$capthist_))$detector <- "multi"
        private$detector_type_ <- ifelse(outputdetector == "count", 1, 3)
      }
      if (private$detector_type_ == 3) {
        private$capij_ <- matrix(0, nr = dim(private$capthist_)[1], nc = dim(private$capthist_)[2])
        for (j in 1:dim(private$capthist_)[2]) {
          #if zero detections that occasion, return NAs
          #otherwise, return which are detections
          whichtrap <- as.numeric(apply(private$capthist_[,j,], 1, FUN = function(x) {which(x > 0)}))
          if(length(whichtrap) > 0){
            private$capij_[,j] <- whichtrap
          } else {
            private$capij_[,j] <- rep(NA, dim(private$capthist_)[1])
          }
        }
        private$capij_[is.na(private$capij_)] <- -9
        private$capij_ <- private$capij_ - 1 
      } else {
        private$capij_ <- matrix(-10, nr = 1, nc = 1)
      }
      ## split capthist into primary occasions 
      if (is.null(primary)) primary <- seq(1, dim(capthist)[2])
      private$primary_ <- primary 
      private$n_primary_ <- max(primary)
      private$n_occasions_ <- length(unique(primary))
      private$mesh_ <- mesh
      if (is.null(time)) {
        private$time_ <- seq(1, self$n_occasions("all"))
      } else {
        private$time_ <- time
      }
      ## add built-in covariates
      private$cov_$t <- (1:self$n_occasions("all") - 1)
      private$cov_$T <- as.factor((1:self$n_occasions("all") - 1))
      private$cov_type_ <- c("k", "k")     
      private$cov_$primary <- (1:self$n_occasions() - 1)
      private$cov_$Primary <- as.factor((1:self$n_occasions() - 1))
      private$cov_type_ <- c(private$cov_type_, "p", "p")   
      private$cov_$x <- scale(private$mesh_[,1])[,1]
      private$cov_$y <- scale(private$mesh_[,2])[,1]
      private$cov_type_ <- c(private$cov_type_, "m")
      private$cov_type_ <- c(private$cov_type_, "m")
      if (!(private$detector_type_ %in% 1:7)) stop("Detector type not implemented.")

      if (is.null(userdm)){
        ## compute distance trap-to-mesh
        self$calc_distances()
        private$userdistmat_ <- NULL
      } else {
        private$distances_ <- userdm
        private$userdistmat_ <- userdm
      }
      private$usermeshdistmat_ <- usermeshdm
      private$induse_ <- induse
       

      ## compute distance centroid-to-mesh
      private$ibuf_ <- attributes(mesh)$ibuffer
      private$make_imesh()
    },
    
    #### OUTPUT FUNCTIONS 
    
    #### ACCESSORS 
 
    induse = function(){return(private$induse_)}
    
    #### SUMMARY STATISTICS 

    #### FUNCTIONS 
    
  ), 
  
  private = list(
    capthist_ = NULL, # capture history object 
    capij_ = NULL, 
    mesh_ = NULL, # mesh object 
    time_ = NULL, # vector of occasion start times 
    cov_ = NULL, # list of covariates 
    cov_type_ = NULL, # type of each covariate in cov_ list 
    distances_ = NULL, # matrix of distances from trap to mesh 
    userdistmat_ = NULL, # user defined distance matrix from traps to mesh
    usermeshdistmat_ = NULL, # user defined distance matrix from mesh to mesh
    detector_type_ = NULL, # type of detectors as integer (see initialize)
    primary_ = NULL, # vector of indices for each occasion indexing what primary it belongs to
    n_primary_ = NULL, # number of primary occasions
    n_occasions_ = NULL, # number of occasions
    ibuf_ = NULL, # buffer around individual centroids
    imesh_ = NULL, # mesh points within ibuf_ to each individual's centroid 
    induse_ = NULL, # individual trap usage from moving detector
    
    #### FUNCTIONS
    
    # check input into intialize 
    check_input = function(capthist, mesh, time, primary, userdm, usermeshdm, induse) {
      if (!("capthist" %in% class(capthist))) stop("Invalid capture history object.")
      if (!("mask" %in% class(mesh))) stop("Invalid mesh object.")
      if (!is.null(time)) {
        if (!is.null(primary)) {
          if (length(time) != length(unique(primary))) stop("Length of time vector not equal to number of primary occasions.")
        } else if (length(time) != dim(capthist)[2]) {
          stop("Length of time vector not equal to number of occasions.")
        }
        if (!is.numeric(time) | !is.vector(time)) stop("Time is not a numeric vector.")
        if (max(abs(sort(time) - time)) > 1e-10) stop("Time must be an increasing vector of numbers.")
      }
      if (!is.null(primary)) {
        if (!is.numeric(primary) | !is.vector(primary)) stop("Primary is not a numeric vector.")
        if (length(primary) != dim(capthist)[2]) stop("Length of primary vector not equal to number of occasions.")
        if (max(abs(round(primary,0) - primary)) > 1e-10) stop("Primary must be integer labels.")
        if (max(abs(sort(primary) - primary)) > 1e-10) stop("Primary must be an increasing vector of numbers.")
        nprim <- max(primary)
        testprim <- sort(unique(primary))
        if (length(testprim) != nprim) stop("Primary must contain integers 1:maximum primary.")
        if (max(abs(testprim - 1:nprim)) > 0.5) stop("Primary must contain integers 1:maximum primary.")
      }
      if (!is.null(userdm)){
        if(dim(userdm)[1] != dim(capthist)[3]) stop("Distance matrix dimension 1 must be number of traps.")
        if(dim(userdm)[2] != nrow(mesh)) stop("Distance matrix dimention 2 must be number of mesh points.")
        if(!is.numeric(userdm)) stop("Distance matrix must be numeric.")
      }
      if (!is.null(usermeshdm)){
        if(is.null(userdm)) stop("If using a user defined distance matrix for the mesh, you must supply a distance matrix for traps by mesh.")
        if(dim(usermeshdm)[1] != dim(usermeshdm)[2]) stop("Mesh distance matrix must be a square matrix.")
        if(!is.numeric(usermeshdm)) stop("Mesh distance matrix must be numeric.")
        if(!all(diag(meshdist)==0)) stop("Diagonals of mesh distance matrix should be zeros.")
      }
      if(!is.null(induse)){
        if(length(dim(induse))!=3) stop("Individual use matrix should be i x j x k dimension.")
        if(dim(induse)[1] != dim(capthist)[1]) stop("Individual use must be specified for each individual in capture history (first dimension).")
        if(dim(induse)[2] != dim(capthist)[3]) stop("Individual use must be specified for each trap (second dimension).")
        if(dim(induse)[3] != dim(capthist)[2]) stop("Individual use must be specified for each secondary occasion (third dimesion).")
        if(is.null(secr::usage(secr::traps(capthist)))) stop("Trap usage must be specified if also using individual usage matrix.")
       }
      return(0)
    }
  )
)
