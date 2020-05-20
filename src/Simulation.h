#pragma once
#include "DDModel_cpp.h"
#include "DDRep_cpp.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <chrono>
#include <iomanip>

// class to store a complete evaluation of a given "Parameter" Set with "Fit" und Representation "Rep"
class EVAL_format
{
public:
  EVAL_format(){}
  EVAL_format(DDModel_cpp DDModel_cpp_)
  {
    Rep = DDRep_cpp(DDModel_cpp_);
  }
  double                Fit = 0.0;
  DDRep_cpp             Rep;
  bool operator<(const EVAL_format& rhs) const
  {
    return Fit < rhs.Fit;
  }
};

// class for a simulation
class Simulation
{
public:
  Simulation(DDModel_cpp DDModel_cpp_) :Model{DDModel_cpp_}
  {
    PAR_Model.resize(Model.ModelParameter.length());
    EVAL.resize(Model.Parameter.length()+4); // +4 ==  4 additional buffer points for SIMPLEX
    for (std::size_t i = 0; i<EVAL.size();++i)
    {
      EVAL[i] = EVAL_format(DDModel_cpp_);
    }
    trials=10000;
    S_Sampling = false;
    n_s_sample = 1;
  }
  Simulation(DDModel_cpp DDModel_cpp_, long trials_) :Model{DDModel_cpp_}
  {
    PAR_Model.resize(Model.ModelParameter.length());
    EVAL.resize(Model.Parameter.length()+4); // +4 ==  4 additional buffer points for SIMPLEX
    for (std::size_t i = 0; i<EVAL.size();++i)
    {
      EVAL[i] = EVAL_format(DDModel_cpp_);
    }
    trials = trials_;
    S_Sampling = false;
    n_s_sample = 1;
  }
  Simulation(DDModel_cpp DDModel_cpp_,Rcpp::S4 DDRep_) :Model{DDModel_cpp_}
  {
    PAR_Model.resize(Model.ModelParameter.length());
    EVAL.resize(Model.Parameter.length()+4); // +4 ==  4 additional buffer points for SIMPLEX
    for (std::size_t i = 0; i<EVAL.size();++i)
    {
      EVAL[i] = EVAL_format(DDModel_cpp_);
    }
    TBF = EVAL_format(DDModel_cpp_);
    TBF.Rep = DDRep_cpp(DDRep_);
    trials = 10000;
    n_s_sample = 1;
    S_Sampling = false;
  }
  // variables
  std::string                       Dir;          // working directory
  DDModel_cpp                       Model;        // Model used
  // Simulation variables
  long                              trials;       // number of trials in simulation
  bool                              S_Sampling;   // logical for application of super sampling
  int                               n_s_sample;   // Number of super samples
  long                              n_s_trials;   // Trials inside a super sample
  RAW_format                        response;     // Single response of da diffusion proces [0]=time[ms], [1] = response [1==correct,0==incorrect]
  std::vector<double>               PAR_Model;    // Buffer for the true modelparameters (after transformation with modelmatrix)
  // Fit variables
  EVAL_format                       TBF;          // Data that should be fitted
  std::vector<EVAL_format>          EVAL;         // Evaluationsbuffer used in SIMPLEX
  std::vector<EVAL_format>          RESULT;       // Buffer of possible results
  std::string                       fit_method;   // Fit method used
  std::string                       start_method; // Start value method used
  std::string                       perf_ana;
  // Grid variables
  int                               n_GRID;       // Number of values that are to be saved from grid
  // SIMPLEX variables
  std::vector<int>                  SIMPLEX_struc; // Structure of simplexes used for fitting
  double                            SIMPLEX_alpha; // Reflection coef
  double                            SIMPLEX_beta;  // Contraction coef
  double                            SIMPLEX_gamma; // Expansion coef
  double                            SIMPLEX_sigma; // Shrink coef
  double                            SIMPLEX_rtoll; // Simplex tollerance
  int                               SIMPLEX_nfunc; // maximum number of iterations
  int                               SIMPLEX_nshrink; // maximum number of shrinks
  // Misc vars (perfomrmance analysis)
  double                            self_fit;


  // Functions for Simulation/Parameterinitialisation/Fit-evalutation
  void Simulate(int Set);                             // Simulates a given parameterset in EVAL[Set] and writes the corresponding representation
  void Simulate_and_Fit(int Set);                     // Simulates and fits a given parameterset in EVAL[Set] to data in TBF
  void FitCrit_Get(int Set);                          // Calculates the Fit of given Data in EVAL[Set] to "TBF"
  void PAR_Model_Get(int Set, int cond);              // Calculates the true modelparameters "PAR_Model" for a given EVAL[Set] and condition "cond"
  void PAR_Init_Rnd();                                // Initialize EVAL[0] with random parameters in the domain specified in "Model"
  void PAR_Init_Man(std::vector<double> PAR_);        // Initialize EVAL[0] manually with a parameterset
  void PAR_Init_from_Result(int ind);                 // Initialize EVAL[0] with parameters from the result buffer
  void response_DSTP();                               // Generate a response from the DSTP
  void response_DMC();                                // Generate a response from the DMC
  void response_SSP();                                // Generate a response from the SSP
  void response_DDM_classic();                                // Generate a response from the classic DDM
  void REP_Get(int Set);                              // Generates a Representation for a given EVAL[Set]
  void SIM_Init_SS(int m_trials);                     // Initializes Super Sampling

  // Functions for the Simplex
  void SIMPLEX();                                     // Runs a Simplex on a given EVAL
  void SIMPLEX_TransformSimplex(int ihi);               // Subroutine for Simplex-manipulation
  void SIMPLEX_TransformSimplex_star2_beta(int ihi);    // Subroutine for Simplex-manipulation
  void SIMPLEX_TransformSimplex_star2_gamma(int ihi);   // Subroutine for Simplex-manipulatio
  void SIMPLEX_Init_coef(Rcpp::List lcoef);
  void Run_SIMPLEX_struc(bool rnd);                                 // Runs Simplex structure
  // Functions for the Grid-calculation
  void GRID_Get_ParComb(std::vector<int> maxes);                    // Calculates parametercombinations with stepsizes equal to "maxes" and writes them in RESULT
  void GRID_Split(int nS,std::string name);                         // Exports parametercombinations (generated by GRID_Get_ParComb) to "nS" numbered files with the basename "name" into a subdirecotry
  void GRID_Get(int depth, std::vector<std::vector<double>> & SEARCH, std::vector<double> & INIT, std::vector<int> & maxes);  // Recursive function to generate the parametercombination of a given Grid
  void GRID_IN(std::ifstream &grid, bool eval);                                   // Reads one line (i.e. one evalutaion point) of a given Grid
  void GRID_Read(std::vector<std::string> grid_parts,bool eval);                  // Reeas a Grid, build from "grid_parts" files
  void GRID_Simulate_ParComb(std::ofstream& outstream, std::ifstream &instream);  // Simulates parametercombination found in instream and writes the evaluation in the ofstream

  // Misc
  Rcpp::S4 Get_DDFit(EVAL_format &EF);                                            // Returns a DDFit object for a given EVAL_format object
  void Performance_Analysis(EVAL_format &EF);


};
