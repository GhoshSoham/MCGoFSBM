// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_edgelist
arma::mat get_edgelist(const arma::mat& adj, const int E, const int n);
RcppExport SEXP _MCGoFSBM_get_edgelist(SEXP adjSEXP, SEXP ESEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type adj(adjSEXP);
    Rcpp::traits::input_parameter< const int >::type E(ESEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(get_edgelist(adj, E, n));
    return rcpp_result_gen;
END_RCPP
}
// get_edgelist_comp
arma::mat get_edgelist_comp(const arma::mat& adj, const int E, const int n);
RcppExport SEXP _MCGoFSBM_get_edgelist_comp(SEXP adjSEXP, SEXP ESEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type adj(adjSEXP);
    Rcpp::traits::input_parameter< const int >::type E(ESEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(get_edgelist_comp(adj, E, n));
    return rcpp_result_gen;
END_RCPP
}
// get_mle_cpp
arma::mat get_mle_cpp(const arma::mat& A, const arma::vec& C, const int k);
RcppExport SEXP _MCGoFSBM_get_mle_cpp(SEXP ASEXP, SEXP CSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(get_mle_cpp(A, C, k));
    return rcpp_result_gen;
END_RCPP
}
// graphchi_cpp
double graphchi_cpp(const arma::mat& A, const arma::mat& p_mle, const arma::vec& C, const int n, const int k);
RcppExport SEXP _MCGoFSBM_graphchi_cpp(SEXP ASEXP, SEXP p_mleSEXP, SEXP CSEXP, SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type p_mle(p_mleSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(graphchi_cpp(A, p_mle, C, n, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MCGoFSBM_get_edgelist", (DL_FUNC) &_MCGoFSBM_get_edgelist, 3},
    {"_MCGoFSBM_get_edgelist_comp", (DL_FUNC) &_MCGoFSBM_get_edgelist_comp, 3},
    {"_MCGoFSBM_get_mle_cpp", (DL_FUNC) &_MCGoFSBM_get_mle_cpp, 3},
    {"_MCGoFSBM_graphchi_cpp", (DL_FUNC) &_MCGoFSBM_graphchi_cpp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_MCGoFSBM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
