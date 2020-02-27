#include "DDRep_cpp.h"

Rcpp::S4 DDRep_cpp::Convert_to_S4(){
  Rcpp::List RAW_List;
  Rcpp::List RF_List;
  Rcpp::List CDF_List;
  Rcpp::List CAF_List;
  Rcpp::List REP_List;
   for (int c = 0; c<CDF.size();++c)
   {
   Rcpp::NumericVector time;
   Rcpp::NumericVector resp;
   Rcpp::CharacterVector cond;
     for (int i = 0; i<RAW[c].size();++i)
     {
     time.push_back(RAW[c][i].time);
     resp.push_back(RAW[c][i].resp);
     cond.push_back(RAW[c][i].cond);
     }
   RAW_List.push_back(Rcpp::DataFrame::create( Rcpp::Named("cond") = cond,Rcpp::Named("resp") = resp,Rcpp::Named("time")=time));
   }
  Rcpp::CharacterVector cond_names_buff;
  for (int c = 0; c< CDF.size();++c)
  {
    Rcpp::NumericVector time;
    Rcpp::NumericVector N;
    Rcpp::CharacterVector cond;
    Rcpp::NumericVector perc;
    for (int i = 0; i<CDF[c].size();++i)
    {
      time.push_back(CDF[c][i].time);
      N.push_back(CDF[c][i].N);
      cond.push_back(CDF[c][i].cond);
      perc.push_back(CDF[c][i].perc);
    }
    cond_names_buff.push_back(CDF[c][0].cond);
    Rcpp::DataFrame CDF_frame = Rcpp::DataFrame::create( Rcpp::Named("cond") = cond,Rcpp::Named("perc") = perc,Rcpp::Named("time")=time,Rcpp::Named("N")=N);
    CDF_List.push_back(CDF_frame);
  }
  CDF_List.names() = cond_names_buff;
  for (int c = 0; c<CAF.size();++c)
  {
    Rcpp::NumericVector time_caf;
    Rcpp::NumericVector N_A;
    Rcpp::NumericVector N_B;
    Rcpp::CharacterVector cond_caf;
    Rcpp::NumericVector perc_caf;
    Rcpp::NumericVector acc;
    for (int i = 0; i<CAF[c].size();++i)
    {
      time_caf.push_back(CAF[c][i].time);
      N_A.push_back(CAF[c][i].N_A);
      N_B.push_back(CAF[c][i].N_B);
      cond_caf.push_back(CAF[c][i].cond);
      perc_caf.push_back(CAF[c][i].perc);
      acc.push_back(CAF[c][i].acc);
    }
    Rcpp::DataFrame CAF_frame = Rcpp::DataFrame::create( Rcpp::Named("cond") = cond_caf,Rcpp::Named("perc") = perc_caf,Rcpp::Named("time")=time_caf,Rcpp::Named("acc")=acc,Rcpp::Named("N_A")=N_A,Rcpp::Named("N_B")=N_B);
    CAF_List.push_back(CAF_frame);
  }
  CAF_List.names()=cond_names_buff;
  Rcpp::NumericVector CDF_Form;
  Rcpp::NumericVector CAF_Form;
  for (int i = 0; i<RF[0].size();++i)
  {
    CDF_Form.push_back(RF[0][i]);
  }
  for (int i = 0; i<RF[1].size();++i)
  {
    CAF_Form.push_back(RF[1][i]);
  }
  RF_List.push_back(CDF_Form);
  RF_List.push_back(CAF_Form);
  RF_List.names() = Rcpp::CharacterVector({"CDF","CAF"});
  REP_List.push_back(CDF_List);
  REP_List.push_back(CAF_List);
  REP_List.names() = Rcpp::CharacterVector({"CDF","CAF"});
  RAW_List.names() = cond_names_buff;
  Rcpp::DataFrame PAR_frame = Rcpp::DataFrame::create();
  for (int i = 0; i<PAR_v.size();++i)
  {
    PAR_frame.push_back(PAR_v[i]);
  }
  PAR_frame.names() = PAR_n;
  Rcpp::S4 DDREP_buff("DDRep");
  DDREP_buff.slot("RAW") = RAW_List;
  DDREP_buff.slot("REP") = REP_List;
  DDREP_buff.slot("RF") = RF_List;
  DDREP_buff.slot("PAR") = Rcpp::as<Rcpp::DataFrame>(PAR_frame);
  return(DDREP_buff);
}
Rcpp::List DDRep_cpp::Convert_FORM_to_List(){
  Rcpp::List FORM_vec;
  Rcpp::NumericVector CDF_Form;
  Rcpp::NumericVector CAF_Form;
  for (int i = 0; i<RF[0].size();++i)
  {
    CDF_Form.push_back(RF[0][i]);
  }
  for (int i = 0; i<RF[1].size();++i)
  {
    CAF_Form.push_back(RF[1][i]);
  }
  FORM_vec.push_back(CDF_Form);
  FORM_vec.push_back(CAF_Form);
  FORM_vec.names() = Rcpp::CharacterVector({"CDF","CAF"});
  return(FORM_vec);
}
