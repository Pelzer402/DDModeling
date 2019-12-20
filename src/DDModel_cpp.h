#pragma once
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include "PRNG.h"
#include <vector>

class DDModel_cpp
{
public:
  // constructor
  DDModel_cpp(Rcpp::S4 DDModel_)
  {
    ID   = Rcpp::as<std::string>(DDModel_.slot("ID"));
    DM   = Rcpp::as<Rcpp::NumericMatrix>(DDModel_.slot("DM"));
    SP   = Rcpp::as<Rcpp::NumericMatrix>(DDModel_.slot("SP"));
    RF   = Rcpp::as<Rcpp::List>(DDModel_.slot("RF"));
    for (int i = 0; i<Rcpp::as<Rcpp::List>(DDModel_.slot("MM")).length();++i)
    {
      MM.push_back(Rcpp::as<Rcpp::NumericMatrix>(Rcpp::as<Rcpp::List>(DDModel_.slot("MM"))[i]));
    }
    MM.names() = Rcpp::as<Rcpp::List>(DDModel_.slot("MM")).names();
    Conditions = MM.names();
    Parameter = Rcpp::colnames(DM);
    ModelParameter = Rcpp::rownames(MM[0]);
    dt = SP(0,0);
    sigma = SP(0,1);
  }
  // variables
  std::string ID;
  Rcpp::ListMatrix MM;
  Rcpp::NumericMatrix DM;
  Rcpp::NumericMatrix SP;
  Rcpp::List RF;
  Rcpp::CharacterVector Parameter;
  Rcpp::CharacterVector ModelParameter;
  Rcpp::CharacterVector Conditions;
  double dt;
  double sigma;
  // functions
  Rcpp::S4 Convert_to_S4();
};
