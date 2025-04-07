// Copyright (c) 2018 Richard Glennie, University of St Andrews
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
// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>
#include "forward_declare.h"

using namespace RcppParallel; 

//'
 //' @param num_cells vector length 3, first entry is how many mesh 
 //' @param sd movement parameter for each occasion
 //' @param dx mesh spacing
 //' @param inside matrix of which cells are closeby 
 //' @param meshdistmat distances between mesh  
 //'
 //' @return transition rate matrix 
 //' 
 // [[Rcpp::export]]
arma::sp_mat CalcTrm(const arma::vec num_cells, 
                     const double sd, 
                     const double dx, 
                     const arma::mat inside,
                     const arma::mat meshdistmat) {
  arma::sp_mat tpr = arma::zeros<arma::sp_mat>(num_cells(0), num_cells(0)); //sparse square matrix, dim number of mesh cells
  int icol = inside.n_cols;
  //arma::mat meshdistmat(num_cells(0), num_cells(0), arma::fill::value(dx));
  double sum; 
  for (int s = 0; s < num_cells(0); ++s) {
    sum = 0; 
    for (int i = 0; i < icol; ++i) {
      if (inside(s, i) > -0.5) {
        //rate
        tpr(s, inside(s, i)) = sd * sd / (2 * meshdistmat(s, inside(s, i)) * meshdistmat(s, inside(s, i)));
        sum += tpr(s, inside(s, i)); 
      }
    }
    tpr(s, s) = -sum; 
  }
  return tpr;
}

arma::vec ExpG(const arma::vec& v_in, //prob ch & m (for kp and s), length [m]
                  const arma::sp_mat& a, //trm (for kp and s), dim [m x m]
                  const double& t, //dt for kp
                  const int& krylov_dim, //dimension of Krylov subspace for approx
                  const double& tol) { //error tolerance for approx
  arma::rowvec v(v_in.t()); //transposes v_in
  double m = fmin(a.n_rows, krylov_dim); //Krylov subspace dim
  double anorm = arma::norm(a, "Inf"); //infinity norm (max abs rowSum), should be +rate leaving m
  double mxrej = 10;
  double mx;
  double btol = 1.0e-7;
  double gamma = 0.9;
  double mb = m;
  int nstep = 0;
  double t_now = 0;
  double t_step;
  double delta = 1.2;
  double t_out = fabs(t); //absolute value
  double s_error = 0;
  double rndoff = anorm * 1e-16;
  
  int k1 = 1;
  double xm = 1 / m;
  double normv = norm(v); //Euclidean norm (if v is almost 0, don't bother)
  if (normv < 1e-16) return v_in; 
  double avnorm;
  double beta = normv;
  double fact = std::pow((m + 1) / std::exp(1), m + 1) * std::sqrt(2 * M_PI * (m + 1));
  double t_new = (1.0 / anorm) * std::pow((fact * tol) / (4 * beta * anorm), xm);
  double s = std::pow(10, std::floor(std::log10(t_new)) - 1);
  t_new = std::ceil(t_new / s) * s;
  double sgn = t > 0 ? 1 : -1;
  int ireject;
  double err_loc;
  double phi1;
  double phi2;
  
  arma::rowvec w = v;
  double hump = normv;
  arma::mat vmat = arma::zeros<arma::mat>(m + 1, a.n_rows);
  arma::mat hmat = arma::zeros<arma::mat>(m + 2, m + 2);
  arma::mat fmat;
  arma::rowvec p;
  while (t_now < t_out) {
    ++nstep;
    t_step = fmin(t_out - t_now, t_new);
    vmat.zeros();
    hmat.zeros();
    vmat.row(0) = (1 / beta) * w;
    for (int j = 0; j < m; ++j) {
      p = vmat.row(j) * a;
      for (int i = 0; i <= j; ++i) {
        hmat(j, i) = arma::dot(vmat.row(i), p);
        p -= vmat.row(i) * hmat(j, i);
      }
      s = norm(p);
      if (s < btol) {
        k1 = 0;
        mb = j;
        t_step = t_out - t_now;
        break;
      }
      hmat(j, j + 1) = s;
      vmat.row(j + 1) = (1 / s) * p;
    }
    if (k1 != 0) {
      hmat(m, m + 1) = 1;
      avnorm = arma::norm(vmat.row(m) * a);
    }
    ireject = 0;
    while (ireject <= mxrej) {
      mx = mb + k1;
      fmat = arma::expmat(sgn * t_step * hmat.submat(0, 0, mx, mx));
      if (k1 == 0) {
        err_loc = btol;
        break;
      }
      else {
        phi1 = fabs(beta * fmat(0, m));
        phi2 = fabs(beta * fmat(0, m + 1) * avnorm);
        if (phi1 > 10 * phi2) {
          err_loc = phi2;
          xm = 1 / m;
        }
        else if (phi1 > phi2) {
          err_loc = (phi1 * phi2) / (phi1 - phi2);
          xm = 1 / m;
        }
        else {
          err_loc = phi1;
          xm = 1 / (m - 1);
        }
      }
      if (err_loc <= delta * t_step * tol) break;
      else {
        t_step = gamma * t_step * std::pow(t_step * tol / err_loc, xm);
        s = std::pow(10, std::floor(std::log10(t_step)) - 1);
        t_step = std::ceil(t_step / s) * s;
        if (ireject == mxrej) {
          //std::cout << "error: requested tolerance too high for Krylov approximation" << std::endl;
        }
        ++ireject;
      }
    }
    mx = mb + fmax(0, k1 - 1);
    w = beta * fmat.row(0).cols(0, mx) * vmat.rows(0, mx);
    beta = arma::norm(w);
    hump = fmax(hump, beta);
    
    t_now = t_now + t_step;
    t_new = gamma * t_step * std::pow(t_step * tol / err_loc, xm);
    s = std::pow(10, std::floor(std::log10(t_new) - 1));
    t_new = std::ceil(t_new / s) * s;
    
    err_loc = fmax(err_loc, rndoff);
    s_error += err_loc;
  }
  double err = s_error;
  hump = hump / normv;
  arma::vec v_out(w.t()); 
  return v_out; 
}

struct MoveLlkCalculator : public Worker { //inherits parallelization
  
  // input 
  const int n; 
  const int Kp;
  const arma::mat pr0; 
  const Rcpp::List pr_capture; 
  const Rcpp::List tpms;
  const arma::vec num_cells; 
  const arma::mat inside; 
  const double dx; 
  const arma::vec dt; 
  const arma::mat sd; 
  const int num_states;
  const int minstate; 
  const int maxstate; 
  const arma::mat meshdistmat;
  const arma::vec entry; 
  
  // transform 
  std::vector<arma::mat> tpm; 
  std::vector<arma::cube> pr_cap; 
  std::vector<arma::sp_mat> trm; 
  int alivestates; 
  
  // output 
  arma::vec& illk; 
  
  // initialiser
  MoveLlkCalculator(const int n, 
                    const int Kp, //number of primary occasions
                const arma::mat pr0, 
                const Rcpp::List pr_capture, 
                const Rcpp::List tpms,
                const arma::vec num_cells, 
                const arma::mat inside, 
                const double dx, 
                const arma::vec dt, 
                const arma::mat sd,
                const int num_states,
                const int minstate, 
                const int maxstate, 
                const arma::mat meshdistmat,
                const arma::vec entry,
                arma::vec& illk) : n(n), Kp(Kp), pr0(pr0), pr_capture(pr_capture), tpms(tpms), num_cells(num_cells), inside(inside), dx(dx), dt(dt), sd(sd), num_states(num_states), minstate(minstate), maxstate(maxstate), meshdistmat(meshdistmat), entry(entry), illk(illk) {
    if (num_states > 1) {
      tpm.resize(Kp);  //length of vector (of matrices)
      for (int kp = 0; kp < Kp - 1; ++kp) tpm[kp] = Rcpp::as<arma::mat>(tpms[kp]); 
    }
    alivestates = num_states - minstate - maxstate; 
    trm.resize(Kp * alivestates);  //different transitions for each alive state
    for (int kp = 0; kp < Kp - 1; ++kp) {
      for (int g = minstate; g < minstate + alivestates; ++g) {
        if (sd(kp, g - minstate) < 0) continue; 
        //trm is vector of matrices
        trm[g - minstate + kp * alivestates] = CalcTrm(num_cells, sd(kp, g - minstate), dx, inside, meshdistmat);  //returns m x m sparse rate matrix
      }
    }
    pr_cap.resize(n); //vector of cubes
    for (int i = 0; i < n; ++i) {
      Rcpp::NumericVector pr_capvec(pr_capture[i]);
      arma::cube pr_icap(pr_capvec.begin(), num_cells(0), num_states, Kp, false);
      pr_cap[i] = pr_icap; //shallow copy, restructure from list to vector
    }
  }
  
  void operator()(std::size_t begin, std::size_t end) { //set up parallel
    for (int i = begin; i < end; ++i) {
      double llk = 0;
      double sum_pr;
      arma::mat pr = pr0; //initial distribution m x s
      arma::cube prcap; //vector (length n) of prob capthist m x s x kp
      for (int kp = entry(i); kp < Kp - 1; ++kp) { //entry is vector of 0's
        pr %= pr_cap[i].slice(kp); //initial probability * pr capthist [m x s]
        if (num_states > 1) {
          pr *= tpm[kp]; //times transition probs (if multiple alive states)
        }
        for (int g = minstate; g < minstate + alivestates; ++g) {
          if (sd(kp, g - minstate) < 0) continue; 
          try {
            // for kp g, pr.col(g) is prob of ch & state [m] of s cols
            pr.col(g) = ExpG(pr.col(g), trm[g - minstate + kp * alivestates], dt(kp));
          } catch(...) {
            llk = -arma::datum::inf;
          }
        }
        sum_pr = accu(pr);
        llk += log(sum_pr);
        pr /= sum_pr;
      }
      pr %= pr_cap[i].slice(Kp - 1);
      llk += log(accu(pr));
      illk(i) = llk;
    }
  }
};


//' Computes log-likelihood 
//'
//' @param n number of individuals 
//' @param Kp total number of occasions 
//' @param pr0 initial distribution over life states
//' @param pr_capture output of calc_pr_capture() in JsModel
//' @param tpms output of calc_tpms() in JsModel
//' @param num_cells number of cells in x,y,total 
//' @param inside 0 if meshpt outside survey region, 1 otherwise 
//' @param dx mesh spacing 
//' @param dt time between occasions 
//' @param sd movement parameter for each occasion 
//' @param num_states 2 = CJS model, 3 = JS model 
//' @param entry time each individual entered survey 
//'
//' @return log-likelihood value 
//' 
// [[Rcpp::export]]
double C_calc_move_llk(const int n, const int Kp,
                       const arma::mat pr0, 
                       const Rcpp::List pr_capture, 
                       const Rcpp::List tpms,
                       const arma::vec num_cells, 
                       const arma::mat inside, 
                       const double dx, 
                       const arma::vec dt, 
                       const arma::mat sd, 
                       const int num_states,
                       const int minstate, 
                       const int maxstate, 
                       const arma::mat meshdistmat,
                       const arma::vec entry) {
  
  arma::vec illk(n);
  MoveLlkCalculator move_llk_calulator(n, Kp, pr0, pr_capture, tpms, num_cells, inside, dx, dt, sd, num_states, minstate, maxstate, meshdistmat, entry, illk); 
  parallelFor(0, n, move_llk_calulator, 10); 
  return(arma::accu(illk)); 
}

//' Computes detection probability (seen at least once) for Jolly-Seber model 
//'
//' @param Kp total number of occasions 
//' @param pr0 initial distribution over life states
//' @param pr_captures list of empty capture histories, see calc_pdet() in JsModel
//' @param tpms output of calc_tpms() in JsModel
//' @param num_cells number of cells in x,y,total 
//' @param inside 0 if meshpt outside survey region, 1 otherwise 
//' @param dx mesh spacing 
//' @param dt time between occasions 
//' @param sd movement parameter for each occasion 
//' @param num_states 2 = CJS model, 3 = JS model 
//'
//' @return pdet = probability seen at some time on the survey 
//' 
// [[Rcpp::export]]
double C_calc_move_pdet(const int Kp, 
                   arma::mat pr0, 
                   Rcpp::List pr_captures,
                   Rcpp::List tpms,
                   const arma::vec num_cells, 
                   const arma::mat inside, 
                   const double dx, 
                   const arma::vec dt,
                   const arma::mat sd, 
                   const int num_states, 
                   const int minstate, 
                   const int maxstate,
                   const arma::mat meshdistmat) {
  
  double pdet = 0; 
  arma::mat pr(pr0);
  arma::sp_mat trm(num_cells(0), num_cells(0)); 
  double sum_pr;
  arma::mat tpm;
  arma::mat pr_capture;
  int alive_col = 0; 
  int alivestates = num_states - minstate - maxstate; 
  if (num_states > 2) alive_col = 1; 
  for (int kp = 0; kp < Kp - 1; ++kp) {
    pr_capture = Rcpp::as<arma::mat>(pr_captures[kp]);
    pr %= pr_capture;
    if (num_states > 1) {
      tpm = Rcpp::as<arma::mat>(tpms[kp]); 
      pr *= tpm; 
    }
    for (int g = minstate; g < minstate + alivestates; ++g) {
      if (sd(kp, g - minstate) < 0) continue; 
      trm = CalcTrm(num_cells, sd(kp, g - minstate), dx, inside, meshdistmat);
      try {
        pr.col(g) = ExpG(pr.col(g), trm, dt(kp)); 
      } catch(...) {
        return -arma::datum::inf; 
      }
    }
    sum_pr = accu(pr); 
    pdet += log(sum_pr); 
    pr /= sum_pr; 
  }
  pr_capture = Rcpp::as<arma::mat>(pr_captures[Kp - 1]);
  pr %= pr_capture;
  pdet += log(accu(pr)); 
  pdet = exp(pdet); 
  return(pdet); 
}



