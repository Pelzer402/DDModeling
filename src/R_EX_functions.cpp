#include "Simulation.h"
#include "PRNG.h"

// Here functions for the export in the R-environment are listed


// TESST
// [[Rcpp::export(.TEST)]]
std::vector<double> TEST(double lb, double ub, int n){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  std::vector<double> OUT;
  for (int i = 0; i<n; ++i)
  {
    OUT.push_back(urand(lb,ub));
  }
  return(OUT);
}

// TESST
// [[Rcpp::export(.TEST2)]]
std::vector<double> TEST2(double mu, double s, int n){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  std::vector<double> OUT;
  for (int i = 0; i<n; ++i)
  {
    OUT.push_back(mu+nrand()*s);
  }
  return(OUT);
}

// Functions for Simulation
// [[Rcpp::export(.Sim_DDModel_rnd)]]
Rcpp::S4 Generate_Modelprediction_rnd(Rcpp::S4 DDModel_,long trials_){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  DDModel_cpp M(DDModel_);
  Simulation S(M,trials_);
  S.PAR_Init_Rnd();
  S.Simulate(0);
  return(S.EVAL[0].Rep.Convert_to_S4());
}


// [[Rcpp::export(.Sim_DDModel_par)]]
Rcpp::S4 Generate_Modelprediction_par(Rcpp::S4 DDModel_,long trials_,std::vector<double> param_){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  DDModel_cpp M(DDModel_);
  Simulation S(M,trials_);
  S.PAR_Init_Man(param_);
  S.Simulate(0);
  return(S.EVAL[0].Rep.Convert_to_S4());
}

// Functions for fitting
// [[Rcpp::export(.Fit_DDModel_rnd)]]
Rcpp::S4 Fit_observed_data_rnd(Rcpp::List calc_cluster){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  Rcpp::S4 DDModel_ = Rcpp::as<Rcpp::S4>(calc_cluster[0]);
  Rcpp::S4 DDRep_   = Rcpp::as<Rcpp::S4>(calc_cluster[1]);
  bool S_Sampling = Rcpp::as<bool>(calc_cluster[2]);
  int m_trials = Rcpp::as<int>(calc_cluster[3]);
  DDModel_cpp M(DDModel_);
  Simulation S(M,DDRep_);
  S.start_method = "Random";
  S.SIMPLEX_struc = Rcpp::as<std::vector<int>>(calc_cluster[4]);
  Rcpp::List scoef = Rcpp::as<Rcpp::List>(calc_cluster[5]);
  S.SIMPLEX_Init_coef(scoef);
  S.fit_method = "";
  for (std::size_t i = 0; i<S.SIMPLEX_struc.size(); ++ i)
  {
    S.fit_method = S.fit_method + std::to_string(S.SIMPLEX_struc[i]) + " Simplex";
    if (i<S.SIMPLEX_struc.size()-1)
    {
      S.fit_method =  S.fit_method  + " -> ";
    }
  }
  S.S_Sampling = S_Sampling;
  S.SIM_Init_SS(m_trials);
  S.Run_SIMPLEX_struc(true);
  return(S.Get_DDFit(S.RESULT[0]));
}

// [[Rcpp::export(.Fit_DDModel_grid)]]
Rcpp::S4 Fit_observed_data_grid(Rcpp::List calc_cluster){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  Rcpp::S4 DDModel_ = Rcpp::as<Rcpp::S4>(calc_cluster[0]);
  Rcpp::S4 DDRep_   = Rcpp::as<Rcpp::S4>(calc_cluster[1]);
  std::vector<std::string> grid_parts =Rcpp::as<std::vector<std::string>>(calc_cluster[2]);
  bool S_Sampling = Rcpp::as<bool>(calc_cluster[3]);
  int m_trials = Rcpp::as<int>(calc_cluster[4]);
  std::string grid_path = Rcpp::as<std::string>(calc_cluster[5]);
  DDModel_cpp M(DDModel_);
  Simulation S(M,DDRep_);
  S.start_method = "Grid: " + grid_path;
  S.SIMPLEX_struc = Rcpp::as<std::vector<int>>(calc_cluster[6]);
  Rcpp::List scoef = Rcpp::as<Rcpp::List>(calc_cluster[7]);
  S.SIMPLEX_Init_coef(scoef);
  if (S.SIMPLEX_struc[0] == 0)
  {
    S.n_GRID = 1;
  }
  else
  {
    S.n_GRID = S.SIMPLEX_struc[0];
  }
  S.fit_method = "";
  for (std::size_t i = 0; i<S.SIMPLEX_struc.size(); ++ i)
  {
    S.fit_method = S.fit_method + std::to_string(S.SIMPLEX_struc[i]) + " Simplex";
    if (i<S.SIMPLEX_struc.size()-1)
    {
      S.fit_method =  S.fit_method  + " -> ";
    }
  }
  for ( std::size_t i = 0; i<grid_parts.size();++i)
  {
    std::replace(grid_parts[i].begin(),grid_parts[i].end(),'/','\\');
  }
  S.GRID_Read(grid_parts,true);
  S.S_Sampling = S_Sampling;
  S.SIM_Init_SS(m_trials);
  S.Run_SIMPLEX_struc(false);
  return(S.Get_DDFit(S.RESULT[0]));
}

// [[Rcpp::export(.Fit_DDModel_DL)]]
Rcpp::S4 Fit_observed_data_DL(Rcpp::List calc_cluster){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  Rcpp::S4 DDModel_ = Rcpp::as<Rcpp::S4>(calc_cluster[0]);
  Rcpp::S4 DDRep_   = Rcpp::as<Rcpp::S4>(calc_cluster[1]);
  bool S_Sampling = Rcpp::as<bool>(calc_cluster[2]);
  int m_trials = Rcpp::as<int>(calc_cluster[3]);
  DDModel_cpp M(DDModel_);
  Simulation S(M,DDRep_);
  S.start_method = "DL Prediction";
  S.SIMPLEX_struc = Rcpp::as<std::vector<int>>(calc_cluster[4]);
  std::vector<double> PRE = Rcpp::as<std::vector<double>>(calc_cluster[5]);
  Rcpp::List scoef = Rcpp::as<Rcpp::List>(calc_cluster[6]);
  S.SIMPLEX_Init_coef(scoef);
  S.fit_method = "";
  for (std::size_t i = 0; i<S.SIMPLEX_struc.size(); ++ i)
  {
    S.fit_method = S.fit_method + std::to_string(S.SIMPLEX_struc[i]) + " Simplex";
    if (i<S.SIMPLEX_struc.size()-1)
    {
      S.fit_method =  S.fit_method  + " -> ";
    }
  }
  S.S_Sampling = S_Sampling;
  S.SIM_Init_SS(m_trials);
  S.PAR_Init_Man(PRE);
  for (int i = 0; i <S.Model.Parameter.size()+1;++i)
  {
    S.Simulate_and_Fit(i);
  }
  S.RESULT.push_back(S.EVAL[0]);
  if (S.SIMPLEX_struc[0] == 0)
  {
    return(S.Get_DDFit(S.RESULT[0]));
  }
  else
  {
    S.SIMPLEX();
    std::sort(S.RESULT.begin(), S.RESULT.end());
    return(S.Get_DDFit(S.RESULT[0]));
  }
}


// Functions for Grid Calculation:
// [[Rcpp::export(.Get_ParComb_cpp)]]
void Calculate_Parameter_Combinations(Rcpp::S4 DDModel_,std::string wd,std::string name,std::vector<int> steps,int nSplit)
{
  DDModel_cpp M(DDModel_);
  Simulation S(M);
  S.Dir = wd;
  std::replace(S.Dir.begin(), S.Dir.end(), '/', '\\');
  S.GRID_Get_ParComb(steps);
  S.GRID_Split(nSplit,name);
}
// [[Rcpp::export(.Get_GRID_cpp)]]
void Grid_calc(Rcpp::List calc_cluster)
{
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  Rcpp::S4 DDModel_ = Rcpp::as<Rcpp::S4>(calc_cluster[0]);
  std::string pc_path = Rcpp::as<std::string>(calc_cluster[1]);
  std::replace(pc_path.begin(), pc_path.end(), '/', '\\');
  std::string out_path = Rcpp::as<std::string>(calc_cluster[2]);
  std::replace(out_path.begin(), out_path.end(), '/', '\\');
  DDModel_cpp M(DDModel_);
  Simulation S(M);
  S.trials = 10000;
  std::ifstream ParComb_in(pc_path.c_str());
  std::ofstream out(out_path.c_str());
  S.GRID_Simulate_ParComb(out,ParComb_in);
}

// [[Rcpp::export(.GRID_to_DDRep)]]
Rcpp::List GRID_to_DDRep(Rcpp::List calc_cluster)
{
  Rcpp::S4 DDModel_ = Rcpp::as<Rcpp::S4>(calc_cluster[0]);
  std::vector<std::string> grid_parts =Rcpp::as<std::vector<std::string>>(calc_cluster[1]);
  for ( std::size_t i = 0; i<grid_parts.size();++i)
  {
    std::replace(grid_parts[i].begin(),grid_parts[i].end(),'/','\\');
  }
  DDModel_cpp M(DDModel_);
  Simulation S(M);
  S.n_GRID = Rcpp::as<int>(calc_cluster[2]);
  S.GRID_Read(grid_parts,false);
  Rcpp::List OUT;
  for (std::size_t i = 0; i<S.RESULT.size();++i)
  {
    OUT.push_back(S.RESULT[i].Rep.Convert_to_S4());
  }
  return(OUT);
}


// Functions for DDRep calculation
// [[Rcpp::export(.DDRep_cpp)]]
Rcpp::S4 Generate_DDRep(Rcpp::S4 DDModel_,Rcpp::List RAW_){
  DDModel_cpp M(DDModel_);
  Simulation S(M);
  DDRep_cpp REP_buff(M,RAW_);
  S.EVAL[0].Rep = REP_buff;
  S.REP_Get(0);
  return(S.EVAL[0].Rep.Convert_to_S4());
}

// [[Rcpp::export(.Reshape_DDRep_cpp)]]
Rcpp::S4 Reshape_DDRep(Rcpp::S4 DDModel_,Rcpp::S4 DDRep_){
  DDModel_cpp M(DDModel_);
  Simulation S(M);
  DDRep_cpp REP_buff(M);
  DDRep_cpp REP_IN(DDRep_);
  REP_buff.RAW =  REP_IN.RAW;
  REP_buff.PAR_n  =  REP_IN.PAR_n;
  REP_buff.PAR_v  =  REP_IN.PAR_v;
  S.EVAL[0].Rep = REP_buff;
  S.REP_Get(0);
  return(S.EVAL[0].Rep.Convert_to_S4());
}

