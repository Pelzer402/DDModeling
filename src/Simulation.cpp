#include "Simulation.h"
#include "Distributions.h"
void Simulation::Simulate(int Set){
  for(int c = 0; c<Model.Conditions.length();++c)
  {
   EVAL[Set].Rep.RAW[c].clear();
   PAR_Model_Get(Set,c);
    response.cond = Model.Conditions[c];
    for (int t=0; t<trials;++t)
    {
      if (Model.ID == "DSTP")
      {
        response_DSTP();
        EVAL[Set].Rep.RAW[c].push_back(response);
      }
      else if (Model.ID == "DMC")
      {
        response_DMC();
        EVAL[Set].Rep.RAW[c].push_back(response);
      }
      else if (Model.ID == "SSP")
      {
        response_SSP();
        EVAL[Set].Rep.RAW[c].push_back(response);
      }
    }
    std::sort(EVAL[Set].Rep.RAW[c].begin(),EVAL[Set].Rep.RAW[c].end());
  }
  REP_Get(Set);
}

void Simulation::Simulate_and_Fit(int Set){
  if (S_Sampling)
  {
    std::vector<EVAL_format> E_buff;
    for (int i = 0; i<n_s_sample;++i)
    {
      Simulate(Set);
      FitCrit_Get(Set);
      E_buff.push_back(EVAL[Set]);
    }
    std::sort(E_buff.begin(),E_buff.end());
    EVAL[Set] = E_buff[0];
  }
  else
  {
    Simulate(Set);
    FitCrit_Get(Set);
  }
}

void Simulation::response_DSTP(){
  int ready_RS1 = 0;
  int ready_RS2 = 0;
  long t = 0;
  double I_RS = 0; // x0 RS
  double I_SS = 0; // x0 SS
  double a          = PAR_Model[1];
  double A          = a/2;
  double B          = -A;
  double c          = PAR_Model[2];
  double C          = c/2;
  double D          = -c/2;
  double mu_RS1     = PAR_Model[3];
  double mu_RS2_C   = PAR_Model[4];
  double mu_RS2_D   = PAR_Model[5];
  double mu_SS      = PAR_Model[6];
  double dt         = Model.dt;
  double sqrt_dt    = std::sqrt(dt);
  double sigma      = Model.sigma;
  long   Ter        = (long)(PAR_Model[0]/dt);  //Parameter Ter is in seconds

  while (!ready_RS1)
  {
    I_SS += mu_SS*dt + sigma*sqrt_dt*nrand();
    I_RS += mu_RS1*dt + sigma*sqrt_dt*nrand();
    t++;
    if ((I_RS <= B) || (I_RS >= A))
    {
      ready_RS1 = 1;
    }
    else
    {
      if ((I_SS <= D) || (I_SS >= C))
      {
        ready_RS1 = 2;
      }
    }
  }
  if (ready_RS1 == 1)
  {
    if (I_RS <= B)
    {
      response.resp = 0;
    }
    else
    {
      response.resp = 1;
    }
  }
  else if (ready_RS1 == 2)
  {
    if (I_SS >= C)
    {
      while(!ready_RS2)
      {
        I_RS += mu_RS2_C*dt+sigma*sqrt_dt*nrand();
        t++;
        if (I_RS <=B)
        {
          response.resp = 0;
          ready_RS2 = 1;
        }
        else if (I_RS >= A)
        {
          response.resp = 1;
          ready_RS2 = 1;
        }
      }
    }
    else if (I_SS <= D)
    {
      while(!ready_RS2)
      {
        I_RS += mu_RS2_D*dt+sigma*sqrt_dt*nrand();
        t++;
        if (I_RS <=B)
        {
          response.resp = 0;
          ready_RS2 = 1;
        }
        else if (I_RS >= A)
        {
          response.resp = 1;
          ready_RS2 = 1;
        }
      }
    }
  }
  response.time = t + Ter;
}

void Simulation::response_DMC()
{
  int ready;											// Pseudo bool for loop exit
  long t = 0;
  double I = 0; // x0
  double a          = PAR_Model[1];
  double A          = I + a/2.0;
  double B          = I - a/2.0;
  double zeta       = PAR_Model[2];
  double alpha      = PAR_Model[3];
  double mu_c       = PAR_Model[4];
  double tau        = PAR_Model[5];
  double dt         = Model.dt*1000; // conversion!!
  double sqrt_dt    = std::sqrt(dt);
  double sigma      = Model.sigma;
  long   Ter        = (long)(PAR_Model[0]);  //Parameter Ter is in milliseconds
  ready = 0;
  while (!ready)
  {
    I += zeta * (std::exp(-(t + 1) / tau)* std::pow(std::exp(1)*(t + 1) / ((alpha - 1)*tau), alpha - 1)*((alpha - 1) / (t + 1) - 1 / tau)) * dt + mu_c * dt + sqrt_dt *sigma*nrand();												// Iterate Process
    t++;												// Add 1 time quant to the overall time

    // Evaluation of the exit condition
    if (I >= A)
    {
      ready = 1;
      response.resp = 1;
    }
    else if (I <= B)
    {
      ready = 1;
      response.resp = 0;
    }
  }
  response.time = t + Ter;
}
void Simulation::response_SSP()
{
  int ready;											// Pseudo bool for loop exit
  long t = 0;
  double I = 0; // x0
  double a          = PAR_Model[1];
  double A          = I + a/2;
  double B          = I - a/2;
  double P          = PAR_Model[2];
  double sda        = PAR_Model[3];
  double rd         = PAR_Model[4];
  double dt         = Model.dt;
  double sqrt_dt    = std::sqrt(dt);
  double sigma      = Model.sigma;
  long   Ter        = (long)(PAR_Model[0]/dt);  //Parameter Ter is in milliseconds
  double tar = 2.5;									// Location of Target
  double inner, outer, cent;

  ready = 0;
  while (!ready)
  {
    if ((sda - rd * (t + 1)) < 0.001)
    {
      NormalDistribution Gauss1(tar, 0.001);		// Make a normal distribution (mu = tar, sigma = sda-rd*(i+1))
      inner = 2.0*(Gauss1.cdf(2.0) - Gauss1.cdf(1.0));
      cent = Gauss1.cdf(3.0) - Gauss1.cdf(2.0);
      outer = 1 - inner - cent;
      I += (inner*P + outer * P + cent * std::abs(P))*dt + sigma*sqrt_dt*nrand();
      t++;												// Add 1 time quant to the overall time
      // Evaluation of the exit condition
      if (I >= A)
      {
        ready = 1;
        response.resp = 1;
      }
      else if (I <= B)
      {
        ready = 1;
        response.resp = 0;
      }
     // Rcpp::Rcout << "Drift_m: " << P << " VALUE: "<<inner*P + outer * P + cent * std::abs(P) << std::endl;
    }
    else
    {
      NormalDistribution Gauss1(tar, sda - rd * (t + 1));		// Make a normal distribution (mu = tar, sigma = sda-rd*(i+1))
      inner = 2.0*(Gauss1.cdf(2.0) - Gauss1.cdf(1.0));
      cent = Gauss1.cdf(3.0) - Gauss1.cdf(2.0);
      outer = 1 - inner - cent;
      I += (inner*P + outer * P + cent * std::abs(P))*dt + sigma*sqrt_dt*nrand();
      t++;												// Add 1 time quant to the overall time
      // Evaluation of the exit condition
      if (I >= A)
      {
        ready = 1;
        response.resp = 1;
      }
      else if (I <= B)
      {
        ready = 1;
        response.resp = 0;
      }
      //Rcpp::Rcout << "Drift: " << P << " VALUE: "<<inner*P + outer*P  + cent*std::abs(P)  << std::endl;
    }
  }
  response.time = t + Ter;
}



void Simulation::PAR_Model_Get(int Set, int cond){
  double buff;
  Rcpp::NumericMatrix MM = Model.MM[cond];
  for (int mp = 0; mp<Model.ModelParameter.length();++mp)
  {
    buff = 0.0;
    for (int cp = 0; cp<Model.Parameter.length();++cp)
    {
      buff +=MM(mp,cp)*EVAL[Set].Rep.PAR_v[cp];
    }
    PAR_Model[mp]= buff;
  }
}

void Simulation::Run_SIMPLEX_struc(bool rnd){
  if (rnd)
  {
    for (int SORTING = 0; SORTING<SIMPLEX_struc.size(); ++SORTING)
    {
      for ( int j = 0;j<SIMPLEX_struc[SORTING];++j)
      {
        if (SORTING == 0)
        {
          PAR_Init_Rnd();
        }
        else
        {
          PAR_Init_from_Result(j);
        }
        for (int i = 0; i <Model.Parameter.size()+1;++i)
        {
          Simulate_and_Fit(i);
        }
        SIMPLEX();
      }
      std::sort(RESULT.begin(), RESULT.end());
    }
  }
  else
  {
    for (int SORTING = 0; SORTING<SIMPLEX_struc.size(); ++SORTING)
    {
      for ( int j = 0;j<SIMPLEX_struc[SORTING];++j)
      {
        PAR_Init_from_Result(j);
        for (int i = 0; i <Model.Parameter.size()+1;++i)
        {
          Simulate_and_Fit(i);
        }
        SIMPLEX();
      }
      std::sort(RESULT.begin(), RESULT.end());
    }
  }
}

void Simulation::SIM_Init_SS(int m_trials){
  if (S_Sampling == false)
  {
    trials = m_trials;
    n_s_sample = 1;
    n_s_trials = trials;
  }
  else
  {
    n_s_trials = m_trials;
    trials = 0;
    for (int p = 0; p<TBF.Rep.CAF[0].size();++p)
    {
      trials += TBF.Rep.CAF[0][p].N_A + TBF.Rep.CAF[0][p].N_B;
    }
    n_s_sample = n_s_trials/trials;
  }
}

void Simulation::PAR_Init_Rnd(){
  for (int cp = 0; cp<Model.Parameter.length();++cp)
  {
    EVAL[0].Rep.PAR_v[cp] = urand(Model.DM(1,cp),Model.DM(0,cp));
  }
  std::vector<double> t;
  t.resize(Model.Parameter.length());
  t[0] = 0.05;
  for (int sp = 1; sp < Model.Parameter.length()+1; ++sp)
  {
    for (int cp = 0; cp < Model.Parameter.length(); ++cp)
    {
      EVAL[sp].Rep.PAR_v[cp] = EVAL[0].Rep.PAR_v[cp]*(1.0+t[cp]);
    }
    t.insert(t.begin(), t[t.size() - 1]);
    t.erase(t.end() - 1);
  }
}
void Simulation::PAR_Init_Man(std::vector<double> PAR_)
{
  for (int i = 0; i < Model.Parameter.length(); ++i)
  {
    EVAL[0].Rep.PAR_v[i] = PAR_[i];
  }
  std::vector<double> t;
  t.resize(Model.Parameter.length());
  t[0] = 0.05;
  for (int sp = 1; sp < Model.Parameter.length()+1; ++sp)
  {
    for (int cp = 0; cp < Model.Parameter.length(); ++cp)
    {
      EVAL[sp].Rep.PAR_v[cp] = EVAL[0].Rep.PAR_v[cp]*(1.0+t[cp]);
    }
    t.insert(t.begin(), t[t.size() - 1]);
    t.erase(t.end() - 1);
  }
}

void Simulation::PAR_Init_from_Result(int ind)
{
  for (int i = 0; i < Model.Parameter.length(); ++i)
  {
    EVAL[0].Rep.PAR_v[i] = RESULT[ind].Rep.PAR_v[i];
  }
  std::vector<double> t;
  t.resize(Model.Parameter.length());
  t[0] = 0.05;
  for (int sp = 1; sp < Model.Parameter.length()+1; ++sp)
  {
    for (int cp = 0; cp < Model.Parameter.length(); ++cp)
    {
      EVAL[sp].Rep.PAR_v[cp] = EVAL[0].Rep.PAR_v[cp]*(1.0+t[cp]);
    }
    t.insert(t.begin(), t[t.size() - 1]);
    t.erase(t.end() - 1);
  }
}


void Simulation::REP_Get(int Set){

  for (int c = 0; c<Model.Conditions.length();++c)
  {
    std::sort(EVAL[Set].Rep.RAW[c].begin(),EVAL[Set].Rep.RAW[c].end());
    //split Data in correct and incorrect
    std::vector<RAW_format> RAW_c;
    std::vector<RAW_format> RAW_i;
    for (int i = 0; i<EVAL[Set].Rep.RAW[c].size();++i)
    {
      if(EVAL[Set].Rep.RAW[c][i].resp == 1)
      {
        RAW_c.push_back(EVAL[Set].Rep.RAW[c][i]);
      }
      else
      {
        RAW_i.push_back(EVAL[Set].Rep.RAW[c][i]);
      }
    }
    // CDF korrekter Antworten
    if (RAW_c.size()<EVAL[Set].Rep.CDF[c].size())
    {
      for (int i = 0; i<EVAL[Set].Rep.CDF[c].size();++i)
      {
        EVAL[Set].Rep.CDF[c][i].time = 0.0;
        EVAL[Set].Rep.CDF[c][i].N = 0;
      }
    }
    else
    {
      for (int i = 0; i<EVAL[Set].Rep.CDF[c].size();++i)
      {
        EVAL[Set].Rep.CDF[c][i].time = RAW_c[(long)(RAW_c.size()*EVAL[Set].Rep.CDF[c][i].perc)].time;
        EVAL[Set].Rep.CDF[c][i].N = (RAW_c.size()*EVAL[Set].Rep.CDF[c][i].perc);
      }
    }
    // CAF
    std::vector<RAW_format> tmp_RAW;
    for (int i = 0; i<EVAL[Set].Rep.RAW[c].size();++i)
    {
      for (int j = 0; j<EVAL[Set].Rep.CAF[c].size();++j)
      {
        tmp_RAW.push_back(EVAL[Set].Rep.RAW[c][i]);
      }
    }
    double f_start = 0.0;
    double f_end = 0.0;
    long space = 0;
    double buff_time;
    long buff_N_A;
    for (int b = 0; b < EVAL[Set].Rep.CAF[c].size(); ++b)
    {
      buff_time = 0.0;
      buff_N_A = 0;
      f_start = EVAL[Set].Rep.RF[1][b];
      f_end = EVAL[Set].Rep.RF[1][b+1];
      for (long ub = (long)(f_start*tmp_RAW.size()); ub < (long)(f_end*tmp_RAW.size()); ++ub)
      {
        buff_time += tmp_RAW[ub].time;
        buff_N_A += tmp_RAW[ub].resp;
      }
      space = (long)(f_end*tmp_RAW.size()) - (long)(f_start*tmp_RAW.size());
      EVAL[Set].Rep.CAF[c][b].time = std::round(buff_time / space);
      EVAL[Set].Rep.CAF[c][b].N_B  = space - buff_N_A ;
      EVAL[Set].Rep.CAF[c][b].N_B = std::round(EVAL[Set].Rep.CAF[c][b].N_B *(f_end-f_start));
      EVAL[Set].Rep.CAF[c][b].N_A = std::round(buff_N_A  *(f_end-f_start));
      EVAL[Set].Rep.CAF[c][b].acc = EVAL[Set].Rep.CAF[c][b].N_A / (double)(EVAL[Set].Rep.CAF[c][b].N_A+EVAL[Set].Rep.CAF[c][b].N_B);
      EVAL[Set].Rep.CAF[c][b].perc = (f_end+f_start)/2.0;
    }
  }
}

void Simulation::FitCrit_Get(int Set)
{
  double tmp_FitCrit;																	// Final return Criterium
  double tmp;																			// Tmp buffer
  double prop_sim, prop_sim2;															// proportions of the simulation
  double n_cdf_tbf = 0.0;
  double n_cdf_eval = 0.0;
  tmp = 0.0;
  tmp_FitCrit = 0.0;
  for (int c = 0; c < Model.Conditions.length(); ++c)
  {
    tmp = 0.0;
    n_cdf_tbf= 0.0;
    for (int i = 0; i<TBF.Rep.CAF[c].size();++i )
    {
      n_cdf_tbf += TBF.Rep.CAF[c][i].N_A;
    }
    if (n_cdf_tbf < TBF.Rep.CDF[c].size())
    {
    }
    else
    {
      n_cdf_eval = 0.0;
      for (int i = 0; i<EVAL[Set].Rep.CAF[c].size();++i )
      {
        n_cdf_eval += EVAL[Set].Rep.CAF[c][i].N_A;
      }
      for (int p = 0; p < TBF.Rep.CDF[c].size(); ++p)
      {
        prop_sim = (double)EVAL[Set].Rep.CDF[c][p].time/ (double)TBF.Rep.CDF[c][p].time;
        tmp += std::pow((1 - prop_sim), 2.0);
      }
      prop_sim2 = ((double)EVAL[Set].Rep.CDF[c][0].N /(double)n_cdf_eval)/((double)TBF.Rep.CDF[c][0].N/(double)n_cdf_tbf);
      if (prop_sim2 == 0.0)
      {
        prop_sim2 = 1.0;
      }
      tmp += (TBF.Rep.CDF[c].size() * std::pow((1 - prop_sim2), 2.0));
    }
    tmp_FitCrit += tmp;
    tmp = 0.0;
    for (int p = 0; p < TBF.Rep.CAF[c].size(); ++p)
    {
      if (TBF.Rep.CAF[c][p].acc == 0.0)
      {
        prop_sim = (double)EVAL[Set].Rep.CAF[c][p].time/(double)TBF.Rep.CAF[c][p].time;
        prop_sim2 = (double)EVAL[Set].Rep.CAF[c][p].acc/ 0.001;
        if (prop_sim2 == 0.0)
        {
          prop_sim2 = 1.0;
        }
      }
      else
      {
        prop_sim = (double)EVAL[Set].Rep.CAF[c][p].time/(double)TBF.Rep.CAF[c][p].time;
        prop_sim2 = (double)EVAL[Set].Rep.CAF[c][p].acc/(double)TBF.Rep.CAF[c][p].acc;
      }
      tmp += std::pow((1 - prop_sim), 2.0) + std::pow((1 - prop_sim2), 2.0);
    }
    tmp_FitCrit += tmp;
  }
  EVAL[Set].Fit = tmp_FitCrit;
}

void Simulation::Performance_Analysis(EVAL_format &EF){
  int flag = 0;
  for (int i = 0; i<Model.Parameter.size();++i)
  {
    if (TBF.Rep.PAR_v[i] == 0.0)
    {
      flag += 1;
    }
  }
  if (flag == 0)
  {
    PAR_Init_Man(TBF.Rep.PAR_v);
    Simulate_and_Fit(0);
    self_fit = EVAL[0].Fit;
    bool eta = true;
    double eta_buff = 0.0;
    for (int i = 0; i<Model.Parameter.size();++i)
    {
      eta_buff = std::abs(TBF.Rep.PAR_v[i] - EF.Rep.PAR_v[i])/(Model.DM(0,i)-Model.DM(1,i));
      if (eta_buff >0.05)
      {
        eta = false;
      }
    }
    if (EF.Fit < self_fit)
    {
      if (eta)
      {
       perf_ana = "S";
      }
      else
      {
       perf_ana = "ME";
      }
    }
    else
    {
      if (eta)
      {
        perf_ana = "S";
      }
      else
      {
        perf_ana = "FE";
      }
    }
  }
  else
  {
    self_fit = 9999;
    perf_ana = "NA";
  }
}


Rcpp::S4 Simulation::Get_DDFit(EVAL_format &EF){
  Performance_Analysis(EF);
  Rcpp::S4 DDFit_buff("DDFit");
  DDFit_buff.slot("INP_REP") = TBF.Rep.Convert_to_S4();
  DDFit_buff.slot("FIT_REP") = EF.Rep.Convert_to_S4();
  Rcpp::S4 DDFitPar_buff("DDFitPar");
  DDFitPar_buff.slot("FIT_V") = EF.Fit;
  DDFitPar_buff.slot("FIT_N") = TBF.Rep.CDF.size()*TBF.Rep.CDF[0].size() + TBF.Rep.CAF.size()*TBF.Rep.CAF[0].size() * 2;
  DDFitPar_buff.slot("S_SAMPLING") = S_Sampling;
  DDFitPar_buff.slot("TRIALS_TOTAL") = n_s_sample*trials;
  DDFitPar_buff.slot("TRIALS_SAMPLE") = trials;
  DDFitPar_buff.slot("METHOD") = fit_method;
  DDFitPar_buff.slot("START_VALUE") = start_method;
  DDFitPar_buff.slot("FIT_V_self") = self_fit;
  DDFitPar_buff.slot("PERF") = perf_ana;
  DDFit_buff.slot("MODEL") = Model.Convert_to_S4();
  DDFit_buff.slot("FIT") = DDFitPar_buff;
  return(DDFit_buff);
}


void Simulation::GRID_Get_ParComb(std::vector<int> maxes)
{
  double buffer;	// a double buffer for loops
  std::vector<std::vector<double>> search;								// [Parameter][Step] Array of all parameter values in a given grid search
  std::vector<double> INIT;												// loop array for a parameter point

  search.resize(Model.Parameter.length());
  INIT.resize(Model.Parameter.length());

  for (int i = 0; i < Model.Parameter.length(); ++i)
  {
    buffer = Model.DM(0,i) - Model.DM(1,i);						// Calculate the range of each parameter space
    buffer /= (maxes[i] + 1);													// Get a step size
    for (int j = 0; j < maxes[i]; ++j)
    {
      search[i].push_back(Model.DM(1,i) + buffer * (j + 1));			// save a given step size value for the grid search
    }
  }
  GRID_Get(Model.Parameter.length(), search, INIT, maxes);
}

void Simulation::GRID_Get(int depth, std::vector<std::vector<double>> & SEARCH, std::vector<double> & INIT, std::vector<int> & maxes)
{
  if (depth>0)
  {
    for (int i = 0; i<maxes[depth - 1]; i++)
    {
      INIT[depth - 1] = SEARCH[depth - 1][i];
      GRID_Get(depth - 1, SEARCH, INIT, maxes);
    }
  }
  else
  {
    PAR_Init_Man(INIT);
    EVAL_format tmp;
    tmp.Rep.PAR_v = EVAL[0].Rep.PAR_v;
    RESULT.push_back(tmp);
  }
}

void Simulation::GRID_Split(int nS,std::string name)
{
  int range = (int)RESULT.size();
  int steps = range / nS;
  int mod = range % nS;
  int uB;
  int lB;
  for (int t = 0; t < nS; ++t)
  {
    std::ofstream of_Grid((Dir + "\\" + name + "_" + std::to_string(t+1) +".ParComb").c_str());
    if (t == nS - 1)
    {
      uB = steps * (t + 1) + mod;
      lB = steps * t;
    }
    else
    {
      uB = steps * (t + 1);
      lB = steps * t;
    }
    of_Grid << uB-lB << std::endl; // number of parameter combinatione
    for (int it = lB; it < uB; ++it)
    {
      for (int i = 0; i < Model.Parameter.length(); ++i)
      {
        of_Grid << RESULT[it].Rep.PAR_v[i];
        if (i < Model.Parameter.length() - 1)
        {
          of_Grid << " ";
        }
      }
      of_Grid << std::endl;
    }
  }
}

void Simulation::GRID_IN(std::ifstream &grid, bool eval)
{
  char ichar[64];
  for (int cond = 0; cond < EVAL[0].Rep.CDF.size(); ++cond)
  {
    for (int bin = 0; bin < EVAL[0].Rep.CDF[cond].size(); ++bin)
    {
      grid >> ichar;
      EVAL[0].Rep.CDF[cond][bin].time = atof(ichar);
      // >> ichar;
      //EVAL[0].Rep.CDF[cond][bin].perc = atof(ichar);
      grid >> ichar;
      EVAL[0].Rep.CDF[cond][bin].N = atof(ichar);
    }
  }
  for (int cond = 0; cond < EVAL[0].Rep.CAF.size(); ++cond)
  {
    for (int bin = 0; bin < EVAL[0].Rep.CAF[cond].size(); ++bin)
    {
      grid >> ichar;
      EVAL[0].Rep.CAF[cond][bin].time= atof(ichar);
      //grid >> ichar;
      //EVAL[0].Rep.CAF[cond][bin].perc= atof(ichar);
      grid >> ichar;
      EVAL[0].Rep.CAF[cond][bin].acc = atof(ichar);
      grid >> ichar;
      EVAL[0].Rep.CAF[cond][bin].N_A = atol(ichar);
      grid >> ichar;
      EVAL[0].Rep.CAF[cond][bin].N_B= atol(ichar);
    }
  }
  /*
  for (int cond = 0; cond < EVAL[0].Rep.CAF.size(); ++cond)
  {
    grid >> ichar;
    grid >> ichar;
  }
  */
  for (int par = 0; par < Model.Parameter.length(); ++par)
  {
    grid >> ichar;
    EVAL[0].Rep.PAR_v[par] = atof(ichar);
  }
  if(eval)
  {
    FitCrit_Get(0);
  }
  RESULT.push_back(EVAL[0]);
  if (RESULT.size()>n_GRID)
  {
    std::sort(RESULT.begin(), RESULT.end());
    RESULT.pop_back();
  }
}

void Simulation::GRID_Simulate_ParComb(std::ofstream &outstream, std::ifstream &instream)
{
  long n_pc = 0;
  instream >> n_pc;
  outstream << n_pc;
  outstream << std::endl;
  for (int i = 0; i<n_pc;++i)
  {
    for (int ip = 0; ip < Model.Parameter.length(); ++ip)
    {
      instream >> EVAL[0].Rep.PAR_v[ip];
    }
    Simulate(0);
    for (int cond = 0; cond < Model.Conditions.length() ; ++cond)
    {
        for (int bin = 0; bin < EVAL[0].Rep.CDF[cond].size(); ++bin)
        {
          outstream << EVAL[0].Rep.CDF[cond][bin].time << " ";
          //outstream << EVAL[0].Rep.CDF[cond][bin].perc << " ";
          outstream << EVAL[0].Rep.CDF[cond][bin].N	   << " ";
        }
    }
    for (int cond = 0; cond < Model.Conditions.length(); ++cond)
    {
      for (int bin = 0; bin < EVAL[0].Rep.CAF[cond].size(); ++bin)
      {
        outstream << EVAL[0].Rep.CAF[cond][bin].time << " ";
        //outstream << EVAL[0].Rep.CAF[cond][bin].perc << " ";
        outstream << EVAL[0].Rep.CAF[cond][bin].acc << " ";
        outstream << EVAL[0].Rep.CAF[cond][bin].N_A << " ";
        outstream << EVAL[0].Rep.CAF[cond][bin].N_B << " ";
      }
    }
    for (int cond = 0; cond < Model.Conditions.length(); ++cond)
    {
      long NA_buff = 0;
      long NB_buff = 0;
      for (int bin = 0; bin < EVAL[0].Rep.CAF[cond].size(); ++bin)
      {
        NA_buff += EVAL[0].Rep.CAF[cond][bin].N_A;
        NB_buff += EVAL[0].Rep.CAF[cond][bin].N_B;
      }
      //outstream << NA_buff << " " << NB_buff << " ";
    }
    for (int par = 0; par < Model.Parameter.length(); ++par)
    {
      outstream << EVAL[0].Rep.PAR_v[par] << " ";
    }
    outstream << std::endl;
  }
}

void Simulation::GRID_Read(std::vector<std::string> grid_parts, bool eval)
{
  for ( int i = 0; i<grid_parts.size(); ++i)
  {
    std::ifstream grid_in(grid_parts[i].c_str());
    int n_grid = 0;
    grid_in >> n_grid;
    for ( int k = 0; k<n_grid;++k)
    {
      GRID_IN(grid_in, eval);
    }
    grid_in.close();
  }
  if (eval)
  {
    std::sort(RESULT.begin(), RESULT.end());
  }
}
