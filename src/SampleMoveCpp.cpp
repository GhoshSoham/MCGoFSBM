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


// Estimate the parameter of the model
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


// Calculate the chi square test statistic using the estimated value of the parameters
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


// Sample a graph in the same fiber
// The true model is ERSBM
// adj - adjacency matrix of a graph
// C - numeric vector of size n of block assignment; from 1 to k
// n - no of nodes in the graph
// k - no of blocks in the graph
// [[Rcpp::export]]
arma::mat sample_a_move_cpp(const arma::mat& adj, const arma::vec& C, const int n, const int k){
  // Assigning adj to a new variable for ease of use
  arma::mat A = adj;

  // Getting edge information of the graph
  int num_edges = arma::accu(A) / 2;
  arma::umat all_edges = get_edgelist_cpp(A, C, num_edges, n);

  // Getting edge information of the complement graph
  int num_edges_comp = (n * (n - 1)) / 2 - num_edges ;
  arma::umat comp_edges = get_edgelist_comp_cpp(A, C, num_edges_comp, n);

  // Determine whether the move is inter-block or intra-block
  int type = Rcpp::sample(2, 1)(0);

  if (type == 1) {
    // Adding and deleting edges from same block
    // Sample a block
    int s = Rcpp::sample(k, 1)(0);

    //Find edges within the sampled block for both the graph and complement graph
    arma::uvec index_all = arma::find((all_edges.col(2) == s) && (all_edges.col(3) == s));
    arma::umat to_delete = all_edges.rows(index_all);
    arma::uvec index_comp = arma::find((comp_edges.col(2) == s) && (comp_edges.col(3) == s));
    arma::umat to_add = comp_edges.rows(index_comp);

    // Getting dimension information
    int to_delete_l = to_delete.n_rows;
    int to_add_l = to_add.n_rows;

    // Check whether the sampled block has at least one edge
    if ((to_delete_l > 0) && (to_add_l > 0)) {
      // Sample an edge to add from complement graph and delete from the graph
      int delete_edge = Rcpp::sample(to_delete_l, 1)(0);
      int add_edge = Rcpp::sample(to_add_l, 1)(0);

      // If the sampled block has only one edge in the complement graph, then that will be added,
      // otherwise it will add one edge by sampling randomly from complement graph
      if (to_add_l == 1) {
        A(to_add(0) - 1, to_add(1) - 1) = 1;
        A(to_add(1) - 1, to_add(0) - 1) = 1;
      } else {
        A(to_add(add_edge - 1, 0) - 1, to_add(add_edge - 1, 1) - 1) = 1;
        A(to_add(add_edge - 1, 1) - 1, to_add(add_edge - 1, 0) - 1) = 1;
      }

      // # If the sampled block has only one edge in the graph, then that will be deleted
      // # otherwise it will delete one edge by sampling randomly from graph
      if (to_delete_l == 1) {
        A(to_delete(0) - 1, to_delete(1) - 1) = 0;
        A(to_delete(1) - 1, to_delete(0) - 1) = 0;
      } else {
        A(to_delete(delete_edge - 1, 0) - 1, to_delete(delete_edge - 1, 1) - 1) = 0;
        A(to_delete(delete_edge - 1, 1) - 1, to_delete(delete_edge - 1, 0) - 1) = 0;
      }

    }
  } else if (type == 2) {
    // Adding and deleting edges between two different blocks
    // Sample two different blocks
    IntegerVector temp = Rcpp::sample(k, 2, false);
    int s = temp(0);
    int t = temp(1);

    // Find edges between two fixed blocks for both the graph and complement graph
    arma::uvec index_all = arma::find(((all_edges.col(2) == s) && (all_edges.col(3) == t)) || ((all_edges.col(2) == t) && (all_edges.col(3) == s)));
    arma::umat inter = all_edges.rows(index_all);
    arma::uvec index_comp = arma::find(((comp_edges.col(2) == s) && (comp_edges.col(3) == t)) || ((comp_edges.col(2) == t) && (comp_edges.col(3) == s)));
    arma::umat comp_inter = comp_edges.rows(index_comp);

    // Getting dimension information
    int inter_l = inter.n_rows;
    int comp_inter_l = comp_inter.n_rows;

    // Check whether the sampled blocks have at least one edge
    if ((inter_l > 0) && (comp_inter_l > 0)) {
      // Sample an edge to add from complement graph and delete from the graph
      int delete_edge = Rcpp::sample(inter_l, 1)(0);
      int add_edge = Rcpp::sample(comp_inter_l, 1)(0);

      // If the sampled blocks have only one between edge in the complement graph, then that will be added,
      // otherwise it will add one between edge by sampling randomly from complement graph
      if (comp_inter_l == 1) {
        A(comp_inter(0) - 1, comp_inter(1) - 1) = 1;
        A(comp_inter(1) - 1, comp_inter(0) - 1) = 1;
      } else {
        A(comp_inter(add_edge - 1, 0) - 1, comp_inter(add_edge - 1, 1) - 1) = 1;
        A(comp_inter(add_edge - 1, 1) - 1, comp_inter(add_edge - 1, 0) - 1) = 1;
      }

      // If the sampled blocks have only one between edge in the graph, then that will be added,
      // otherwise it will add one between edge by sampling randomly from graph
      if (inter_l == 1) {
        A(inter(0) - 1, inter(1) - 1) = 0;
        A(inter(1) - 1, inter(0) - 1) = 0;
      } else {
        A(inter(delete_edge - 1, 0) - 1, inter(delete_edge - 1, 1) - 1) = 0;
        A(inter(delete_edge - 1, 1) - 1, inter(delete_edge - 1, 0) - 1) = 0;
      }
    }
  }
  // Output:
  // the adjacency matrix of the graph after one random move
  return A;
}


// Computing chi-square GoF statistic several times for each of the sampled graph on
// same fiber after estimating the parameter
// A - adjacency matrix of a graph
// C - numeric vector of size n of block assignment; from 1 to k
// n - no of nodes in the graph
// k - no of blocks in the graph
// numGraphs - no of times chain will move to produce a new graph
// [[Rcpp::export]]
arma::vec graph_chain_on_fiber(const arma::mat& A, const arma::vec& C, const int n, const int k, const int numGraphs){
  // Defining new variables to store the test statistics value
  arma::vec chi_seq_cpp(numGraphs); chi_seq_cpp.zeros();
  arma::mat A_new = A;
  arma::mat A_current = A_new;

  // Estimating the parameters from observed graph
  // which will remain same after generating a new graph on same fiber
  arma::mat p_mle = get_mle_cpp(A, C, k);

  // Storing the GoF statistic on observed graph
  chi_seq_cpp(0) = graphchi_cpp(A, p_mle, C, n, k);
  for(int iter = 1; iter<numGraphs; iter++){
    // Sampling a new graph on same fiber
    A_current = sample_a_move_cpp(A_new, C, n, k);

    // Computing chi-square GoF test statistic on new sampled graph
    chi_seq_cpp(iter) = graphchi_cpp(A_current, p_mle, C, n, k);
    A_new = A_current;
  }

  // Output
  // chi_seq_cpp - a vector containing value of the statistic of each sample
  return chi_seq_cpp;
}
