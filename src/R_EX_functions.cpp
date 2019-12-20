#include "Simulation.h"
#include "PRNG.h"

// [[Rcpp::export(.Sim_DDModel)]]
Rcpp::S4 Generate_Modelprediction_rnd(Rcpp::S4 DDModel_,long trials_){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  DDModel_cpp M(DDModel_);
  Simulation S(M,trials_);
  S.PAR_Init_Rnd();
  S.Simulate(0);
  return(S.EVAL[0].Rep.Convert_to_S4(0));
}

// [[Rcpp::export(.Fit_DDModel_rnd)]]
Rcpp::List Fit_observed_data_rnd(Rcpp::S4 DDModel_,Rcpp::S4 DDRep_){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  DDModel_cpp M(DDModel_);
  Simulation S(M,DDRep_);
  for ( int j = 0;j<1;++j)
  {
    S.PAR_Init_Rnd();
    for (int i = 0; i <S.Model.Parameter.size()+1;++i)
    {
      S.Simulate_and_Fit(i);
    }
    S.SIMPLEX();
  }
  std::sort(S.RESULT.begin(), S.RESULT.end());
  Rcpp::List RES;
  for ( int i = 0;i<S.RESULT.size();++i)
  {
    RES.push_back(S.Get_DDFit_RESULT(i));
  }
  return(RES);
}

// [[Rcpp::export(.Fit_DDModel_grid)]]
Rcpp::List Fit_observed_data_grid(Rcpp::S4 DDModel_,Rcpp::S4 DDRep_, std::string gridpath){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  DDModel_cpp M(DDModel_);
  Simulation S(M,DDRep_);
  // TO DO:
  // (1): Implement READ GRID
  // (2): Implement EVAL GRID
  // (3): Implement SORT RESULT BUFFER
  for ( int j = 0;j<1;++j)
  {
    S.PAR_Init_Rnd();
    for (int i = 0; i <S.Model.Parameter.size()+1;++i)
    {
      S.Simulate_and_Fit(i);
    }
    S.SIMPLEX();
  }
  std::sort(S.RESULT.begin(), S.RESULT.end());
  Rcpp::List RES;
  for ( int i = 0;i<S.RESULT.size();++i)
  {
    RES.push_back(S.Get_DDFit_RESULT(i));
  }
  return(RES);
}

// [[Rcpp::export(.Get_ParComb)]]
void Calculate_Parameter_Combinations(Rcpp::S4 DDModel_,std::string wd,std::vector<int> steps,int nSplit)
{
  DDModel_cpp M(DDModel_);
  Simulation S(M);
  S.Dir = wd;
  std::replace(S.Dir.begin(), S.Dir.end(), '/', '\\');
  S.GRID_Get_ParComb(steps);
  S.GRID_Split(nSplit);
}

// [[Rcpp::export(.Get_GRID)]]
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
  std::ifstream ParComb_in(pc_path.c_str());
  S.GRID_Read_ParComb(ParComb_in);
  std::ofstream out(out_path.c_str());
  S.GRID_Simulate_ParComb(out);
}

// [[Rcpp::export(.DDRep_cpp)]]
Rcpp::S4 Generate_DDRep(Rcpp::S4 DDModel_,Rcpp::List RAW_){
  DDModel_cpp M(DDModel_);
  Simulation S(M);
  DDRep_cpp REP_buff(M,RAW_);
  S.EVAL[0].Rep = REP_buff;
  S.REP_Get(0);
  return(S.EVAL[0].Rep.Convert_to_S4(1));
}

