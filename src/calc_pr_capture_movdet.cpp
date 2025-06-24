// Copyright (c) 2019 Richard Glennie, University of St Andrews
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files, to deal in the software
// without restriction, including without limitation the right to use, copy,
// modify, publish, distribute, sublicense, and/or sell copies of the software,
// and to permit persons to whom the software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS  OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE THE USE OR OTHER DEALINGS IN
// THE SOFTWARE
//
////////////////////////////////////////////////////////////////////////////////
// openpopscr project: open population spatial capture-recapture modelling
//
// moveds.cpp: C++ functions to compute capture history probabilities 
//
////////////////////////////////////////////////////////////////////////////////
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel; 

struct PrCaptureCalculator : public Worker {
  
  // input 
  const int n; //number of individuals
  const int K; // number of occasions 
  const int J; // number of traps 
  const int M; // number of mesh points 
  const int alive_col; // column that contains alive state 
  const arma::cube& capthist; // capthist history records: individuals x occasion x trap 
  const Rcpp::List enc; // encounter rate: occasion x mesh point x trap
  const arma::mat& usage; // usage of traps: trap x occasion 
  const arma::cube induse; // individual usage: individuals x trap x occasion
  const int num_states; // number of hidden states in life history model 
  const arma::cube& known_state; // -1 if known not to be in that state 
  const int minstate; // states before alive state 
  const int maxstate; // states after alives states
  const int detector_type; // type of detector, see initialize in ScrData  
  const int n_prim; // number of primary occasions 
  const arma::vec& S; // number of secondary occasion per primary 
  const arma::vec& entry; // occasion each individual enters survey 
  const arma::field<arma::vec>& imesh; // mesh points within range of each detector 
  const arma::mat& capik; // efficient storage for multi-detectors 
  
  // derived variables 
  std::vector<arma::cube> enc0; //vector of cubes
  std::vector<arma::mat> total_enc; 
  std::vector<arma::cube> total_enc_ind;  //total enc*induse up to detection
  std::vector<arma::mat> log_total_enc; 
  std::vector<arma::mat> log_total_penc; 
  std::vector<arma::cube> logenc0; 
  std::vector<arma::cube> log_penc; 
  std::vector<arma::cube> log_penc_trap; //log(1-exp(-enc*induse)) for single trap

  // output 
  // capture probability for record of individual x occasion x mesh point
  arma::field<arma::cube>& probfield;
  
  // initialise
  PrCaptureCalculator(const int n,
                      const int K, 
                        const int J, 
                        const int M, 
                        const int alive_col, 
                        const arma::cube& capthist, 
                        const Rcpp::List enc, 
                        const arma::mat& usage, 
                        const arma::cube induse,
                        const int num_states,
                        const int minstate, 
                        const int maxstate, 
                        const arma::cube& known_state, 
                        const int detector_type,
                        const int n_prim, 
                        const arma::vec& S, 
                        const arma::vec& entry,
                        const arma::field<arma::vec>& imesh, 
                        const arma::mat& capik, 
                        arma::field<arma::cube>& probfield) : n(n), K(K), J(J), M(M), 
                        alive_col(alive_col), 
                        capthist(capthist), 
                        enc(enc), 
                        usage(usage), 
                        induse(induse),
                        num_states(num_states), 
                        minstate(minstate), 
                        maxstate(maxstate), 
                        known_state(known_state), 
                        detector_type(detector_type),
                        n_prim(n_prim), 
                        S(S), 
                        entry(entry), 
                        imesh(imesh), 
                        capik(capik), 
                        probfield(probfield) {
    
    //// compute dervied quantities 
    // encounter rate per state 
    enc0.resize(num_states); 
    for (int g = 0; g < num_states; ++g) {
      enc0[g] = Rcpp::as<arma::cube>(enc[g]); //mesh x trap x occasion for each state
    }
    if (detector_type == 3) {
      total_enc.resize(num_states);
      total_enc_ind.resize(num_states);
      log_total_enc.resize(num_states);
      log_total_penc.resize(num_states);
      log_penc_trap.resize(num_states);
      auto idx4D = [=](int n, int m, int j, int k) {
        return n + N * (m + M * (j + J * k));
      }; //for indexing log_penc_trap
    } else if (detector_type == 2) {
      stop("Error: Moving detector likelihood has only been formulated for multi-catch type traps.")
    }
    for (int g = 0; g < num_states; ++g) {
      if (detector_type == 3) { 
        // multi detector: 
        //  log_total_enc: log total hazard over detectors for each mesh point x occasion 
        //  log_total_penc: log total probability of detection for each mesh point x occasion 
        total_enc[g] = arma::zeros<arma::mat>(M, K); 
        total_enc_ind[g] = arma::zeros<arma::cube>(M, n, K); 
        log_total_enc[g] = arma::zeros<arma::mat>(M, K); 
        log_total_penc[g] = arma::zeros<arma::mat>(M, K); 
        log_penc_trap[g] = arma::zeros<arma::vec>(n * M * J * K); // log(1-exp(-enc*induse)) for single trap
        int k = -1; //steps through all secondaries
        for (int prim = 0; prim < n_prim; ++prim) {
          for (int s = 0; s < S(prim); ++s) {
            ++k; 
            total_enc[g].col(k) = enc0[g].slice(k) * usage.col(k); //matrix multiplication
            total_enc_ind[g].slice(k) = enc0[g].slice(k) * t(induse.slice(k));
            for (int i = 0; i < n; ++n){
              for (int j = 0; j < J; ++j) {
                for (int m = 0; m < M; ++m){
                  log_penc_trap[g](idx4D(n, m, j, k)) = enc0(m, j, k) * induse(n, j, k);
                }
              }
            }
          }
        }
        log_penc_trap[g] = 1.0 - exp(-log_penc_trap[g]);
        log_penc_trap[g] = log(log_penc_trap[g] + 1e-16);
        log_total_enc[g] = log(total_enc[g]); 
        log_total_penc[g] = 1.0 - exp(-total_enc[g]); 
        log_total_penc[g] = log(log_total_penc[g] + 1e-16); 
      } else if (detector_type == 2) {
        // proximity detector: 
       stop("Error: Moving detector likelihood has only been formulated for multi-catch type traps.")
      }
      // log encounter rate occasion x mesh point x trap
      logenc0.resize(num_states);
      logenc0[g] = log(enc0[g]); 
    }
  } 
  
  void operator()(std::size_t begin, std::size_t end) { 
    // loop over individuals 
    for (int i = begin; i < end; ++i) {
      arma::vec savedenc(M); 
      double sumcap; 
      int g = 0; // current state processed 
      int k = -1; // current occasion processed 
      // loop over primary occasions
      for (int prim = 0; prim < n_prim; ++prim) {
        bool unseen = true; // individual seen or not in this primary?
        if (prim > entry(i) - 1) {
          // loop over secondary occasions 
          for (int s = 0; s < S(prim); ++s) {
            ++k; // increment occasion number 
            // loop over detectable states
            for (int gp = minstate; gp < minstate + num_states; ++gp) {
              g = gp - minstate; // current state 
              // if state known to be impossible, leave probfield zero 
              if (known_state(i, k, gp) < 0) {
                continue; 
              } else {
                // state is possible 
                // independent detectors 
                if (detector_type != 3) {
                 stop("Error: Moving detector likelihood has only been formulated for multi-catch type traps.")
                }
                // dependent detectors 
                if (detector_type == 3) {
                  if(capik(i,k) > -1) unseen = false; 
                  sumcap = capik(i, k) > -1 ? 1 : 0; //1 if seen, 0 if not
                  if (capik(i, k) > -1) { //if seen
                    //logenc vector length mesh, induse is single value
                    int base_idx = i + n * (0 + M * (capik(i,k) + J * k)); //indexing for 4D array
         
                    for (int m = 0; m < imesh(i).size(); ++m) probfield(i)(imesh(i)(m), gp, prim) {
                      savedenc(imesh(i)(m)) = log_penc_trap[g](base_idx + n * imesh(i)(m));
                      += savedenc(imesh(i)(m)) - sumcap * total_enc_ind[g](imesh(i)(m), i, k); // log(exp(-sum(enc*induse))) + log(1-exp(-enc*induse))
                    }
                    }                                                                               //if detected, this term is 0                   
                  for (int m = 0; m < imesh(i).size(); ++m) probfield(i)(imesh(i)(m), gp, prim) += -(1.0 - sumcap) * total_enc[g](imesh(i)(m), k); 
                }
              }
            }
          }
          for (int gp = minstate; gp < minstate + num_states; ++gp) {
            if (known_state(i, k, gp) < 0) continue; 
            for (int m = 0; m < imesh(i).size(); ++m) probfield(i)(imesh(i)(m), gp, prim) = exp(probfield(i)(imesh(i)(m), gp, prim));
          }
        } else {
          // if not entered in this primary yet 
          k = k + S(prim); // move occasion number to end of primary 
        }
        if (unseen) {
          for (int gpi = 0; gpi < minstate; ++gpi) probfield(i).slice(prim).col(gpi).ones(); 
          for (int gpi = minstate + num_states; gpi < minstate + num_states + maxstate; ++gpi) probfield(i).slice(prim).col(gpi).ones(); 
        }
      }
    }
  }
};  
  
//' Computes probability of each capture record
//'
//' @param n number of individuals 
//' @param K total number of occasions 
//' @param J total number of traps ever used  
//' @param M total number of mesh points
//' @param capthist capthist array 
//' @param enc0 encounter rate array, see calc_pr_capture() in JsModel
//' @param usage matrix with K x J where (k,j) entry is usage of trap j in occasion k
//' @param num_states number of alive states 
//' @param minstate number of states before alive (Scr,Cjs = 0, JS = 1)
//' @param maxstate number of states after alive (Scr = 0, Cjs/js = 1)
//' @param detector_type 1 = count, 2 = proximity/binary, 3 = multi-catch, 4 = transect 
//' @param n_prim number of primary occasions 
//' @param S number of secondary occasions per primary occasion 
//' @param entry occasion each individual entered survey 
//'
//' @return  Array with (i,k,m) entry the probability of capture record for individual i in occasion k given activity centre at mesh point m  
//' 
// [[Rcpp::export]]
arma::field<arma::cube> C_calc_pr_capture(const int n, const int K, const int J, const int M, 
                             const arma::cube& capthist, 
                             const Rcpp::List enc0,
                             const arma::mat usage, 
                             const int num_states,
                             const int minstate, 
                             const int maxstate,
                             const arma::cube& known_state, 
                             const int detector_type, 
                             const int n_prim, 
                             const arma::vec S, 
                             const arma::vec entry, 
                             const arma::field<arma::vec>& imesh, 
                             const arma::mat& capik) {
  int alive_col = 1; 
  if (num_states < 3) alive_col = 0; 
  arma::field<arma::cube> probfield(n);
  for (int i = 0; i < n; ++i) probfield(i) = arma::zeros<arma::cube>(M, num_states + minstate + maxstate, n_prim); 
  PrCaptureCalculator pr_capture_calc(n, K, J, M, alive_col, capthist, enc0, usage, num_states, minstate, maxstate, known_state, detector_type, n_prim, S, entry, imesh, capik, probfield); 
  parallelFor(0, n, pr_capture_calc, 10); 
  return(probfield);
}
