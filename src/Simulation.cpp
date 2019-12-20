#include "Simulation.h"

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

      }
      else if (Model.ID == "SSP")
      {

      }
    }
    std::sort(EVAL[Set].Rep.RAW[c].begin(),EVAL[Set].Rep.RAW[c].end());
  }
  REP_Get(Set);
}

void Simulation::Simulate_and_Fit(int Set){
  Simulate(Set);
  FitCrit_Get(Set);
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

void Simulation::PAR_Model_Get(int Set, int cond){
  double buff;
  Rcpp::NumericMatrix MM = Model.MM[cond];
  for (int mp = 0; mp<Model.ModelParameter.length();++mp)
  {
    buff = 0.0;
    for (int cp = 0; cp<Model.Parameter.length();++cp)
    {
      buff +=MM(mp,cp)*EVAL[Set].Parameter[cp];
    }
    PAR_Model[mp]= buff;
  }
}

void Simulation::PAR_Init_Rnd(){
  for (int cp = 0; cp<Model.Parameter.length();++cp)
  {
    EVAL[0].Parameter[cp] = urand(Model.DM(1,cp),Model.DM(0,cp));
  }
  std::vector<double> t;
  t.resize(Model.Parameter.length());
  t[0] = 0.05;
  for (int sp = 1; sp < Model.Parameter.length()+1; ++sp)
  {
    for (int cp = 0; cp < Model.Parameter.length(); ++cp)
    {
      EVAL[sp].Parameter[cp] = EVAL[0].Parameter[cp]*(1.0+t[cp]);
    }
    t.insert(t.begin(), t[t.size() - 1]);
    t.erase(t.end() - 1);
  }
}
void Simulation::PAR_Init_Man(std::vector<double> PAR_)
{
  for (int i = 0; i < Model.Parameter.length(); ++i)
  {
    EVAL[0].Parameter[i] = PAR_[i];
  }
}

void Simulation::PAR_Init_from_Result(int ind)
{
  for (int i = 0; i < Model.Parameter.length(); ++i)
  {
    EVAL[0].Parameter[i] = RESULT[ind].Parameter[i];
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
  double prop_data;																	// proportion of the input
  double prop_sim, prop_sim2;															// proportions of the simulation
  double n_cdf_tbf = 0.0;
  double n_cdf_eval = 0.0;
  tmp = 0.0;
  tmp_FitCrit = 0.0;
  for (int c = 0; c < Model.Conditions.length(); ++c)
  {
    tmp = 0.0;
    n_cdf_tbf= 0.0;
    for (int i = 0; i<tbf_REP.CAF[c].size();++i )
    {
      n_cdf_tbf += tbf_REP.CAF[c][i].N_A;
    }
    if (n_cdf_tbf < tbf_REP.CDF[c].size())
    {
    }
    else
    {
      n_cdf_eval = 0.0;
      for (int i = 0; i<EVAL[Set].Rep.CAF[c].size();++i )
      {
        n_cdf_eval += EVAL[Set].Rep.CAF[c][i].N_A;
      }
      for (int p = 0; p < tbf_REP.CDF[c].size(); ++p)
      {
        prop_sim = (double)EVAL[Set].Rep.CDF[c][p].time/ (double)tbf_REP.CDF[c][p].time;
        tmp += pow((1 - prop_sim), 2.0);
      }
      prop_sim2 = ((double)EVAL[Set].Rep.CDF[c][0].N /(double)n_cdf_eval)/((double)tbf_REP.CDF[c][0].N/(double)n_cdf_tbf);
      if (prop_sim2 == 0.0)
      {
        prop_sim2 = 1.0;
      }
      tmp += (tbf_REP.CDF[c].size() * pow((1 - prop_sim2), 2.0));
    }
    tmp_FitCrit += tmp;
    tmp = 0.0;
    for (int p = 0; p < tbf_REP.CAF[c].size(); ++p)
    {
      if (tbf_REP.CAF[c][p].acc == 0.0)
      {
        prop_sim = (double)EVAL[Set].Rep.CAF[c][p].time/(double)tbf_REP.CAF[c][p].time;
        prop_sim2 = (double)EVAL[Set].Rep.CAF[c][p].acc/ 0.001;
        if (prop_sim2 == 0.0)
        {
          prop_sim2 = 1.0;
        }
      }
      else
      {
        prop_sim = (double)EVAL[Set].Rep.CAF[c][p].time/(double)tbf_REP.CAF[c][p].time;
        prop_sim2 = (double)EVAL[Set].Rep.CAF[c][p].acc/(double)tbf_REP.CAF[c][p].acc;
      }
      tmp += pow((1 - prop_sim), 2.0) + pow((1 - prop_sim2), 2.0);
    }
    tmp_FitCrit += tmp;
  }
  EVAL[Set].Fit = tmp_FitCrit;
}

Rcpp::S4 Simulation::Get_DDFit_EVAL(int Set){
  Rcpp::S4 DDFit_buff("DDFit");
  DDFit_buff.slot("Input_Rep") = tbf_REP.Convert_to_S4(0);
  DDFit_buff.slot("Fit_Rep") = EVAL[Set].Rep.Convert_to_S4(0);
  Rcpp::S4 DDFitPar_buff("DDFitPar");
  Rcpp::DataFrame Parameter_frame =   Rcpp::DataFrame::create();
  for ( int i = 0; i<EVAL[Set].Parameter.size();++i)
  {
    Parameter_frame.push_back(EVAL[Set].Parameter[i]);
  }
  Parameter_frame.names() = Model.Parameter;
  Parameter_frame = Rcpp::as<Rcpp::DataFrame>(Parameter_frame);
  DDFitPar_buff.slot("Parameter") =Parameter_frame;
  DDFitPar_buff.slot("Fit") = EVAL[Set].Fit;
  DDFitPar_buff.slot("nFit") = 40;
  DDFit_buff.slot("Model") = Model.Convert_to_S4();
  DDFit_buff.slot("Fit") = DDFitPar_buff;
  return(DDFit_buff);
}

Rcpp::S4 Simulation::Get_DDFit_RESULT(int Set){
  Rcpp::S4 DDFit_buff("DDFit");
  DDFit_buff.slot("Input_Rep") = tbf_REP.Convert_to_S4(0);
  DDFit_buff.slot("Fit_Rep") = RESULT[Set].Rep.Convert_to_S4(0);
  Rcpp::S4 DDFitPar_buff("DDFitPar");
  Rcpp::DataFrame Parameter_frame =   Rcpp::DataFrame::create();
  for ( int i = 0; i<RESULT[Set].Parameter.size();++i)
  {
    Parameter_frame.push_back(RESULT[Set].Parameter[i]);
  }
  Parameter_frame.names() = Model.Parameter;
  Parameter_frame = Rcpp::as<Rcpp::DataFrame>(Parameter_frame);
  DDFitPar_buff.slot("Parameter") =Parameter_frame;
  DDFitPar_buff.slot("Fit") = RESULT[Set].Fit;
  DDFitPar_buff.slot("nFit") = 40;
  DDFit_buff.slot("Model") = Model.Convert_to_S4();
  DDFit_buff.slot("Fit") = DDFitPar_buff;
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
    tmp.Parameter = EVAL[0].Parameter;
    RESULT.push_back(tmp);
  }
}

void Simulation::GRID_Split(int nS)
{
  int range = (int)RESULT.size();
  int steps = range / nS;
  int mod = range % nS;
  int uB;
  int lB;
  for (int t = 0; t < nS; ++t)
  {
    std::ofstream of_Grid((Dir + "\\" + "GRID_" + std::to_string(t+1) +".ParComb").c_str());
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
        of_Grid << RESULT[it].Parameter[i];
        if (i < Model.Parameter.length() - 1)
        {
          of_Grid << " ";
        }
      }
      of_Grid << std::endl;
    }
  }
}
/*
void Simulation::GRID_Read(std::ifstream &grid)
{
  for (int cond = 0; cond < Set.PAR_sim.nCond; ++cond)
  {
    for (int type = 0; type < Set.PAR_sim.nRespType; ++type)
    {
      for (int bin = 0; bin < Set.PAR_sim.nCDF; ++bin)
      {
        grid >> MODEL.CDF_t[cond][type][bin];			// [cond][type][bin] (type=0=correct, type=1=incrorrect) time of bin
        grid >> MODEL.CDF_perc[cond][type][bin];			// [cond][type][bin] (type=0=correct, type=1=incrorrect) percentil information
        grid >> MODEL.CDF_n[cond][type][bin];			// [cond][type][bin] (type=0=correct, type=1=incrorrect) number of responses in bin
      }

    }
  }
  for (int cond = 0; cond < Set.PAR_sim.nCond; ++cond)
  {
    for (int bin = 0; bin < Set.PAR_sim.nCAF; ++bin)
    {
      grid >> MODEL.CAF_t[cond][bin];								// [cond][bin] time of bin
      grid >> MODEL.CAF_p[cond][bin];								// [cond][bin] percentage of correct responses
      for (int type = 0; type < Set.PAR_sim.nRespType; ++type)
      {
        grid >> MODEL.CAF_n[cond][type][bin];					// [cond][type][bin] number of (type=0 correct, type=1incorrect) responses in bin
      }
    }
  }
  for (int cond = 0; cond < Set.PAR_sim.nCond; ++cond)
  {
    for (int type = 0; type < Set.PAR_sim.nRespType; ++type)
    {
      grid >> MODEL.CDF_nR[cond][type];							// [cond][type] number of total responses
    }
  }

  for (int par = 0; par < Set.PAR_sim.nPar; ++par)
  {
    grid >> Set.PAR_MAT[0][par];
  }
}
*/
void Simulation::GRID_Read_ParComb(std::ifstream &pc)
{
  long n_pc = 0;
  pc >> n_pc;
  for (int it = 0; it < n_pc; ++it)
  {
    EVAL_format tmp(Model);
    for (int ip = 0; ip < Model.Parameter.length(); ++ip)
    {
      pc >> tmp.Parameter[ip];
    }
    RESULT.push_back(tmp);
  }
}

void Simulation::GRID_Simulate_ParComb(std::ofstream &outstream)
{
  for (int i = 0; i<RESULT.size();++i)
  {
    PAR_Init_from_Result(i);
    Simulate(0);
    for (int cond = 0; cond < Model.Conditions.length() ; ++cond)
    {
        for (int bin = 0; bin < EVAL[0].Rep.CDF[cond].size(); ++bin)
        {
          outstream << EVAL[0].Rep.CDF[cond][bin].time << " ";
          outstream << EVAL[0].Rep.CDF[cond][bin].perc << " ";
          outstream << EVAL[0].Rep.CDF[cond][bin].N	   << " ";
        }
    }
    for (int cond = 0; cond < Model.Conditions.length(); ++cond)
    {
      for (int bin = 0; bin < EVAL[0].Rep.CAF[cond].size(); ++bin)
      {
        outstream << EVAL[0].Rep.CAF[cond][bin].time << " ";
        outstream << EVAL[0].Rep.CAF[cond][bin].perc << " ";
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
      outstream << NA_buff << " " << NB_buff << " ";
    }
    for (int par = 0; par < Model.Parameter.length(); ++par)
    {
      outstream << EVAL[0].Parameter[par] << " ";
    }
    outstream << std::endl;
  }
}
