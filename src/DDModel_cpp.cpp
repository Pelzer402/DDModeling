#include "DDModel_cpp.h"

Rcpp::S4 DDModel_cpp::Convert_to_S4(){
  Rcpp::S4 M("DDModel");
  M.slot("ID") = ID;
  M.slot("MM") = MM;
  M.slot("DM") = DM;
  M.slot("SP") = SP;
  M.slot("RF") = RF;
  return(M);
}
