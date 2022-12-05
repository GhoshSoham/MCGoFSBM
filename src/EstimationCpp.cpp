#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// The true model is ERSBM
// Estimate the corresponding parameter values i.e,
// q_i,j (the probability of edges between block i and j)

// A - adjacency matrix of a graph
// C - numeric vector of size n of block assignment; from 1 to k
// k - no of blocks in the graph
// [[Rcpp::export]]
arma::mat get_mle_cpp(const arma::mat& A, const arma::vec& C, const int k){

  // Initialize block_card, a vector of length of k with all zeros
  arma::colvec block_card(k); block_card.zeros();


  // Initialize table_obs, a k x k matrix with all zeros
  arma::mat table_obs(k, k); table_obs.zeros();

  for (int row = 0; row < k; row ++) {
    arma::uvec indexRow = arma::find(C == row + 1);
    // block_card will contain no of elements in each block
    block_card(row) = indexRow.size();
    for (int col = 0; col < k; col ++) {
      arma::uvec indexCol = arma::find(C == col + 1);
      // table_obs will contain total observed edges between block and within block
      table_obs(row, col) = arma::accu(A(indexRow, indexCol));
    }
  }

  // Calculate total possible edges between and within block based on nodes information
  arma::mat tot_edge = block_card * block_card.t() - arma::diagmat(block_card);


  // Return output
  // Output:
  // the estimated q_i,j matrix (the probability of edges between block i and j)
  return arma::mat (table_obs / tot_edge);
}
