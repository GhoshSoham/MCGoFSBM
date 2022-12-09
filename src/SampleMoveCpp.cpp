#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Getting edge list in matrix
// adj - adjacency matrix of a graph
// E - no of edges in the graph
// n - no of nodes in the graph
// [[Rcpp::export]]
arma::umat get_edgelist_cpp(const arma::mat& adj, const arma::vec& C, const int E, const int n){
  // Initialize edge_mat with all zeros
  arma::umat edge_mat(E, 4); edge_mat.zeros();

  // If there is an edge between two nodes, it will add those two nodes in edge_mat
  int e = 0;
  for (int row = 0; row < n - 1; row ++) {
    for (int col = row + 1; col < n; col ++) {
      if (adj(row, col) == 1){
        edge_mat(e, 0) = row + 1;
        edge_mat(e, 1) = col + 1;
        edge_mat(e, 2) = C(row);
        edge_mat(e, 3) = C(col);
        e = e + 1;
      }
    }
  }

  // Return output
  // edge_mat - a E x 2 matrix containing the edges of the graph
  return edge_mat;
}

// Getting edge list of a complement graph in matrix
// adj - adjacency matrix of a graph
// E - no of edges in the complement graph
// n - no of nodes in the graph
// [[Rcpp::export]]
arma::umat get_edgelist_comp_cpp(const arma::mat& adj, const arma::vec& C, const int E, const int n){
  // Initialize edge_mat with all zeros
  arma::umat edge_mat(E, 4); edge_mat.zeros();

  // If there is no edge between two nodes i.e., there is an edge in the complement graph
  // it will add those two nodes in edge_mat
  int e = 0;
  for (int row = 0; row < n - 1; row ++) {
    for (int col = row + 1; col < n; col ++) {
      if (adj(row, col) == 0){
        edge_mat(e, 0) = row + 1;
        edge_mat(e, 1) = col + 1;
        edge_mat(e, 2) = C(row);
        edge_mat(e, 3) = C(col);
        e = e + 1;
      }
    }
  }

  // Return output
  // edge_mat - a E x 2 matrix containing the edges of the complement graph
  return edge_mat;
}

// A - adjacency matrix of a graph
// p_mle - k*k matrix of MLE table
// C - numeric vector of size n of block assignment; from 1 to k
// n - no of nodes in the graph
// k - no of blocks in the graph
// [[Rcpp::export]]
double graphchi_cpp(const arma::mat& A, const arma::mat& p_mle, const arma::vec& C, const int n, const int k){
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
  double output = arma::accu(arma::pow(degseq.elem(arma::find(exp_mat != 0)) - exp_mat.elem(arma::find(exp_mat != 0)), 2) / exp_mat.elem(arma::find(exp_mat != 0)));

  // Return output
  // Output:
  // the value of the chi-sq statistic
  return output;
}

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
