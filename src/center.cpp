#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
RcppExport SEXP center(SEXP Uin, SEXP ID) {
  NumericMatrix U_tmp(Uin);
  mat U = as<mat>(U_tmp);

  IntegerVector id_tmp(ID);
  uvec id = as<uvec>(id_tmp);
  uvec uid = unique(id);

  int n_obs = U.n_rows;
  int n_col = U.n_cols;
  // int n_clust = uid.n_elem;
  // int n_clust(nclust);

  // mat Uc_means(size(U));
  mat Uout(U);

  // Initialize the cluster identifier flag
  unsigned int clust_id = id[0];
  // Initialize the cluster index flags
  unsigned int clust_start_idx = 0;
  unsigned int clust_end_idx = 0;

  // double c_sum = 0;
  // Temporary vector to store the sums
  // for each column in the matrix Uin
  // for each cluster
  vec u_sums(n_col);
  u_sums.fill(0.0);
  // int c_size = 0; 

  // Loop over the rows of the matrix Uin
  for(int k = 0; k < n_obs; ++k){

    // If id[k] is a new cluster
    if( id[k] != clust_id ){

      // Update the out matrix
      clust_end_idx = k - 1;

      // Loop over the columns of U for the cluster id[k]
      // to create centered elements in Uout
      for(int l = 0; l < n_col; ++l){
	Uout(span(clust_start_idx, clust_end_idx),l) -= u_sums(l) / (clust_end_idx - clust_start_idx + 1);
	// Uout(clust_start_idx,l) = u_sums(l) / (clust_end_idx - clust_start_idx + 1);
	// Uout(clust_end_idx,l) = u_sums(l) / (clust_end_idx - clust_start_idx + 1);
	// Uout(clust_start_idx,l) = u_sums(l) ;
	// Uout(clust_end_idx,l) = u_sums(l) ;
	// Uout(clust_end_idx,l) = k;
      }

      // Update the cluster_id to the new one
      clust_id = id[k];
      // Update the cluster_start_idx to the new one
      clust_start_idx = k;
      // Update the cluster sums
      u_sums = U.row(k).t();
      // // Reset the cluster sums
      // u_sums.fill(0.0);

    } else {

      // Update the cluster sums
      u_sums += U.row(k).t();

    }

  }

  // Center the last cluster
  clust_end_idx = n_obs - 1;

  for(int l = 0; l < n_col; ++l){
    Uout(span(clust_start_idx, clust_end_idx),l) -= u_sums(l) / (clust_end_idx - clust_start_idx + 1);
    // Uout(clust_start_idx,l) -= u_sums(l) / (clust_end_idx - clust_start_idx + 1);
    // Uout(clust_end_idx,l) -= u_sums(l);

  }

  // Uout.rows(clust_start_idx, clust_end_idx).each_col() -= u_sums / (clust_end_idx - clust_start_idx + 1);

  return wrap(Uout);

}
