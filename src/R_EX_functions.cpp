#include "Simulation.h"
#include "PRNG.h"

// Here functions for the export in the R-environment are listed

// [[Rcpp::export(.Sim_DDModel)]]
Rcpp::S4 Generate_Modelprediction_rnd(Rcpp::S4 DDModel_,long trials_){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  DDModel_cpp M(DDModel_);
  Simulation S(M,trials_);
  S.PAR_Init_Rnd();
  S.Simulate(0);
  return(S.EVAL[0].Rep.Convert_to_S4());
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
    RES.push_back(S.Get_DDFit(S.RESULT[i]));
  }
  return(RES);
}

// [[Rcpp::export(.Fit_DDModel_grid)]]
Rcpp::S4 Fit_observed_data_grid(Rcpp::List calc_cluster){
  seed_nrand(std::chrono::system_clock::now().time_since_epoch().count());
  Rcpp::S4 DDModel_ = Rcpp::as<Rcpp::S4>(calc_cluster[0]);
  Rcpp::S4 DDRep_   = Rcpp::as<Rcpp::S4>(calc_cluster[1]);
  std::vector<std::string> grid_parts =Rcpp::as<std::vector<std::string>>(calc_cluster[2]);
  bool S_Sampling = Rcpp::as<bool>(calc_cluster[3]);
  int m_trials = Rcpp::as<int>(calc_cluster[4]);
  DDModel_cpp M(DDModel_);
  Simulation S(M,DDRep_);
  for ( int i = 0; i<grid_parts.size();++i)
  {
    std::replace(grid_parts[i].begin(),grid_parts[i].end(),'/','\\');
  }
  //auto start = std::chrono::system_clock::now();//DEL
  S.GRID_Read(grid_parts);
  //auto end = std::chrono::system_clock::now();//DEL
  //std::chrono::duration<double> elapsed_seconds = end-start; //DEL
  //Rcpp::Rcout << "Finished grid read. Time: " << elapsed_seconds.count() << "s\n";//DEL
  S.S_Sampling = S_Sampling;
  if (S_Sampling == false)
  {
    S.trials = m_trials;
  }
  else
  {
    S.n_s_trials = m_trials;
    S.trials = 0;
    for (int p = 0; p<S.TBF.Rep.CAF[0].size();++p)
    {
      S.trials += S.TBF.Rep.CAF[0][p].N_A + S.TBF.Rep.CAF[0][p].N_B;
    }
    S.n_s_sample = S.n_s_trials/S.trials;
  }
  //start = std::chrono::system_clock::now(); //DEL
  for ( int j = 0;j<20;++j)
  {
    S.PAR_Init_from_Result(j);
    for (int i = 0; i <S.Model.Parameter.size()+1;++i)
    {
      S.Simulate_and_Fit(i);
    }
    S.SIMPLEX();
  }
  std::sort(S.RESULT.begin(), S.RESULT.end());
  for (int j = 0; j<2;++j)
  {
    S.PAR_Init_from_Result(j);
    for (int i = 0; i <S.Model.Parameter.size()+1;++i)
    {
      S.Simulate_and_Fit(i);
    }
    S.SIMPLEX();
  }
  std::sort(S.RESULT.begin(), S.RESULT.end());
  //end = std::chrono::system_clock::now(); //DEL
  //elapsed_seconds = end-start;//DEL
  //Rcpp::Rcout << "Finished Simplex. Time: " << elapsed_seconds.count() << "s\n";//DEL
  S.RESULT.push_back(S.TBF);
  return(S.Get_DDFit(S.RESULT[0]));
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



