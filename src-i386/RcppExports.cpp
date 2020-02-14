// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Generate_Modelprediction_rnd
Rcpp::S4 Generate_Modelprediction_rnd(Rcpp::S4 DDModel_, long trials_);
RcppExport SEXP _DDModeling_Generate_Modelprediction_rnd(SEXP DDModel_SEXP, SEXP trials_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type DDModel_(DDModel_SEXP);
    Rcpp::traits::input_parameter< long >::type trials_(trials_SEXP);
    rcpp_result_gen = Rcpp::wrap(Generate_Modelprediction_rnd(DDModel_, trials_));
    return rcpp_result_gen;
END_RCPP
}
// Fit_observed_data_rnd
Rcpp::List Fit_observed_data_rnd(Rcpp::S4 DDModel_, Rcpp::S4 DDRep_);
RcppExport SEXP _DDModeling_Fit_observed_data_rnd(SEXP DDModel_SEXP, SEXP DDRep_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type DDModel_(DDModel_SEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type DDRep_(DDRep_SEXP);
    rcpp_result_gen = Rcpp::wrap(Fit_observed_data_rnd(DDModel_, DDRep_));
    return rcpp_result_gen;
END_RCPP
}
// Fit_observed_data_grid
Rcpp::S4 Fit_observed_data_grid(Rcpp::List calc_cluster);
RcppExport SEXP _DDModeling_Fit_observed_data_grid(SEXP calc_clusterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type calc_cluster(calc_clusterSEXP);
    rcpp_result_gen = Rcpp::wrap(Fit_observed_data_grid(calc_cluster));
    return rcpp_result_gen;
END_RCPP
}
// Calculate_Parameter_Combinations
void Calculate_Parameter_Combinations(Rcpp::S4 DDModel_, std::string wd, std::string name, std::vector<int> steps, int nSplit);
RcppExport SEXP _DDModeling_Calculate_Parameter_Combinations(SEXP DDModel_SEXP, SEXP wdSEXP, SEXP nameSEXP, SEXP stepsSEXP, SEXP nSplitSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type DDModel_(DDModel_SEXP);
    Rcpp::traits::input_parameter< std::string >::type wd(wdSEXP);
    Rcpp::traits::input_parameter< std::string >::type name(nameSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< int >::type nSplit(nSplitSEXP);
    Calculate_Parameter_Combinations(DDModel_, wd, name, steps, nSplit);
    return R_NilValue;
END_RCPP
}
// Grid_calc
void Grid_calc(Rcpp::List calc_cluster);
RcppExport SEXP _DDModeling_Grid_calc(SEXP calc_clusterSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type calc_cluster(calc_clusterSEXP);
    Grid_calc(calc_cluster);
    return R_NilValue;
END_RCPP
}
// Generate_DDRep
Rcpp::S4 Generate_DDRep(Rcpp::S4 DDModel_, Rcpp::List RAW_);
RcppExport SEXP _DDModeling_Generate_DDRep(SEXP DDModel_SEXP, SEXP RAW_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type DDModel_(DDModel_SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type RAW_(RAW_SEXP);
    rcpp_result_gen = Rcpp::wrap(Generate_DDRep(DDModel_, RAW_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DDModeling_Generate_Modelprediction_rnd", (DL_FUNC) &_DDModeling_Generate_Modelprediction_rnd, 2},
    {"_DDModeling_Fit_observed_data_rnd", (DL_FUNC) &_DDModeling_Fit_observed_data_rnd, 2},
    {"_DDModeling_Fit_observed_data_grid", (DL_FUNC) &_DDModeling_Fit_observed_data_grid, 1},
    {"_DDModeling_Calculate_Parameter_Combinations", (DL_FUNC) &_DDModeling_Calculate_Parameter_Combinations, 5},
    {"_DDModeling_Grid_calc", (DL_FUNC) &_DDModeling_Grid_calc, 1},
    {"_DDModeling_Generate_DDRep", (DL_FUNC) &_DDModeling_Generate_DDRep, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_DDModeling(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
