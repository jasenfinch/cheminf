// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cnvrt
std::string cnvrt(std::string input, const char* inputType, const char* outputType);
RcppExport SEXP _cheminf_cnvrt(SEXP inputSEXP, SEXP inputTypeSEXP, SEXP outputTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type input(inputSEXP);
    Rcpp::traits::input_parameter< const char* >::type inputType(inputTypeSEXP);
    Rcpp::traits::input_parameter< const char* >::type outputType(outputTypeSEXP);
    rcpp_result_gen = Rcpp::wrap(cnvrt(input, inputType, outputType));
    return rcpp_result_gen;
END_RCPP
}
// smilesToMF
std::string smilesToMF(std::string smile);
RcppExport SEXP _cheminf_smilesToMF(SEXP smileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type smile(smileSEXP);
    rcpp_result_gen = Rcpp::wrap(smilesToMF(smile));
    return rcpp_result_gen;
END_RCPP
}
// smilesToAccurateMass
double smilesToAccurateMass(std::string smile);
RcppExport SEXP _cheminf_smilesToAccurateMass(SEXP smileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type smile(smileSEXP);
    rcpp_result_gen = Rcpp::wrap(smilesToAccurateMass(smile));
    return rcpp_result_gen;
END_RCPP
}
// smartsSearch
int smartsSearch(std::string smile, std::string smart);
RcppExport SEXP _cheminf_smartsSearch(SEXP smileSEXP, SEXP smartSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type smile(smileSEXP);
    Rcpp::traits::input_parameter< std::string >::type smart(smartSEXP);
    rcpp_result_gen = Rcpp::wrap(smartsSearch(smile, smart));
    return rcpp_result_gen;
END_RCPP
}
// descriptor
double descriptor(std::string smile, const char* desc);
RcppExport SEXP _cheminf_descriptor(SEXP smileSEXP, SEXP descSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type smile(smileSEXP);
    Rcpp::traits::input_parameter< const char* >::type desc(descSEXP);
    rcpp_result_gen = Rcpp::wrap(descriptor(smile, desc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cheminf_cnvrt", (DL_FUNC) &_cheminf_cnvrt, 3},
    {"_cheminf_smilesToMF", (DL_FUNC) &_cheminf_smilesToMF, 1},
    {"_cheminf_smilesToAccurateMass", (DL_FUNC) &_cheminf_smilesToAccurateMass, 1},
    {"_cheminf_smartsSearch", (DL_FUNC) &_cheminf_smartsSearch, 2},
    {"_cheminf_descriptor", (DL_FUNC) &_cheminf_descriptor, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_cheminf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}