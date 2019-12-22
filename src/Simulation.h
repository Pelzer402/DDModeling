#pragma once
#include "DDModel_cpp.h"
#include "DDRep_cpp.h"
#include <fstream>
#include <iostream>

class EVAL_format
{
public:
  EVAL_format(){}
  EVAL_format(DDModel_cpp DDModel_cpp_)
  {
    Parameter.resize(DDModel_cpp_.Parameter.length());
    Rep = DDRep_cpp(DDModel_cpp_);
  }
  std::vector<double> Parameter;
  double              Fit = 0.0;
  DDRep_cpp             Rep;
  bool operator<(const EVAL_format& rhs) const
  {
    return Fit < rhs.Fit;
  }
};

class Simulation
{
public:
  Simulation(DDModel_cpp DDModel_cpp_) :Model{DDModel_cpp_}
  {
    PAR_Model.resize(Model.ModelParameter.length());
    trials=1000;
    EVAL.resize(Model.Parameter.length()+4); // +4 ==  4 additional buffer points for SIMPLEX
    for (int i = 0; i<EVAL.size();++i)
    {
      EVAL[i] = EVAL_format(DDModel_cpp_);
    }
  }
  Simulation(DDModel_cpp DDModel_cpp_, long trials_) :Model{DDModel_cpp_}
  {
    PAR_Model.resize(Model.ModelParameter.length());
    EVAL.resize(Model.Parameter.length()+4); // +4 ==  4 additional buffer points for SIMPLEX
    for (int i = 0; i<EVAL.size();++i)
    {
      EVAL[i] = EVAL_format(DDModel_cpp_);
    }
    trials = trials_;
  }
  Simulation(DDModel_cpp DDModel_cpp_,Rcpp::S4 DDRep_) :Model{DDModel_cpp_}
  {
    PAR_Model.resize(Model.ModelParameter.length());
    EVAL.resize(Model.Parameter.length()+4); // +4 ==  4 additional buffer points for SIMPLEX
    for (int i = 0; i<EVAL.size();++i)
    {
      EVAL[i] = EVAL_format(DDModel_cpp_);
    }
    TBF = EVAL_format(DDModel_cpp_);
    TBF.Rep = DDRep_cpp(DDRep_);
    trials = 5000;
  }
  // variables
  std::string                       Dir;
  DDModel_cpp                       Model;
  long                              trials;       // Anzahl von Trials
  RAW_format                        response;     // Eine Antwort [0]=Zeit, [1] = resp einer Trial
  std::vector<double>               PAR_Model;    // Ueberfuehrte PAR zu MPAR
  std::vector<EVAL_format>          EVAL;         //Matrix von PAR [Set][Parameter]
  EVAL_format                       TBF;
  std::vector<EVAL_format>          RESULT;
  // functions
  void Simulate(int Set);
  void Simulate_and_Fit(int Set);
  void FitCrit_Get(int Set);
  void PAR_Model_Get(int Set, int cond);
  void PAR_Init_Rnd();
  void PAR_Init_Man(std::vector<double> PAR_);
  void PAR_Init_from_Result(int ind);
  void response_DSTP();
  void REP_Get(int Set);
  void SIMPLEX();
  void SIMPLEX_TransformSimplex(int ihi, double fac);
  void SIMPLEX_TransformSimplex_star2_beta(int ihi, double fac);
  void SIMPLEX_TransformSimplex_star2_gamma(int ihi, double fac);
  void GRID_Get_ParComb(std::vector<int> maxes);
  void GRID_Split(int nS,std::string name);
  void GRID_Get(int depth, std::vector<std::vector<double>> & SEARCH, std::vector<double> & INIT, std::vector<int> & maxes);
  void GRID_IN(std::ifstream &grid);
  void GRID_Read(std::vector<std::string> grid_parts);
  void GRID_Read_ParComb(std::ifstream &pc);
  void GRID_Simulate_ParComb(std::ofstream& outstream);
  Rcpp::S4 Get_DDFit_EVAL(int Set);
  Rcpp::S4 Get_DDFit_RESULT(int Set);
};
