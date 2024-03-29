// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// elements
DataFrame elements();
RcppExport SEXP _mzAnnotation_elements() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(elements());
    return rcpp_result_gen;
END_RCPP
}
// suitableElementRanges
List suitableElementRanges(float mass);
RcppExport SEXP _mzAnnotation_suitableElementRanges(SEXP massSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type mass(massSEXP);
    rcpp_result_gen = Rcpp::wrap(suitableElementRanges(mass));
    return rcpp_result_gen;
END_RCPP
}
// generate
DataFrame generate(double measured_mass, double ppm, int charge, List element_ranges);
RcppExport SEXP _mzAnnotation_generate(SEXP measured_massSEXP, SEXP ppmSEXP, SEXP chargeSEXP, SEXP element_rangesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type measured_mass(measured_massSEXP);
    Rcpp::traits::input_parameter< double >::type ppm(ppmSEXP);
    Rcpp::traits::input_parameter< int >::type charge(chargeSEXP);
    Rcpp::traits::input_parameter< List >::type element_ranges(element_rangesSEXP);
    rcpp_result_gen = Rcpp::wrap(generate(measured_mass, ppm, charge, element_ranges));
    return rcpp_result_gen;
END_RCPP
}
// ppmRange
List ppmRange(double mz, double ppm);
RcppExport SEXP _mzAnnotation_ppmRange(SEXP mzSEXP, SEXP ppmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mz(mzSEXP);
    Rcpp::traits::input_parameter< double >::type ppm(ppmSEXP);
    rcpp_result_gen = Rcpp::wrap(ppmRange(mz, ppm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mzAnnotation_elements", (DL_FUNC) &_mzAnnotation_elements, 0},
    {"_mzAnnotation_suitableElementRanges", (DL_FUNC) &_mzAnnotation_suitableElementRanges, 1},
    {"_mzAnnotation_generate", (DL_FUNC) &_mzAnnotation_generate, 4},
    {"_mzAnnotation_ppmRange", (DL_FUNC) &_mzAnnotation_ppmRange, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_mzAnnotation(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
