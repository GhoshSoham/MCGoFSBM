#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// A - adjacency matrix of a graph
// p_mle - k*k matrix of MLE table
// C - numeric vector of size n of block assignment; from 1 to k
// n - no of nodes in the graph
// k - no of blocks in the graph
// [[Rcpp::export]]
double graphchi_cpp(const arma::mat& A, const arma::mat& p_mle, const arma::vec& C, int n, int k){
  // Initialize degree sequence and expected nodes matrix with all zeros
  arma::mat degseq(n, k); degseq.zeros();
  arma::mat exp_mat(n, k); exp_mat.zeros();

  for(int inode=0; inode<n; inode++){
    arma::vec indexRow = arma::vectorise(A.row(inode));
    for(int iblock=0; iblock<k; iblock++){
      arma::uvec indexCol= arma::find(C == iblock + 1);
      // Observed no of neighbors of a node in each block
      degseq(inode, iblock) = arma::accu(indexRow(indexCol));
      int ni = indexCol.size(); // no of nodes in ith block
      // Expected no of neighbors of a node in each block
      exp_mat(inode, iblock) = p_mle(iblock, C(inode)-1) * ni;
    }
  }

  // Calculate the value of the chi-sq statistic
  double output = arma::accu(arma::pow(degseq.elem(arma::find(exp_mat!=0)) - exp_mat.elem(arma::find(exp_mat!=0)), 2) / exp_mat.elem(arma::find(exp_mat!=0)));

  // Return output
  // Output:
  // the value of the chi-sq statistic
  return output;
}
