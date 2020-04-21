#pragma once
#include <Rcpp.h>
#include "DDModel_cpp.h"
#include "PRNG.h"
#include <vector>


struct RAW_format
{
  double time;
  int  resp;
  std::string  cond;
  bool operator<(const RAW_format& rhs) const
  {
    return time < rhs.time;
  }
};

struct CDF_perc
{
  std::string cond;
  double perc =0.0;
  double time = 0.0;
  long N =0;
};

struct CAF_perc
{
  std::string cond;
  double perc = 0.0;
  double  time = 0.0;
  double acc =0.0;
  long  N_A =0;
  long N_B=0;
};

class DDRep_cpp
{
public:
  DDRep_cpp(){}
  DDRep_cpp(DDModel_cpp DDModel_cpp_)
  {
    RF.resize(DDModel_cpp_.RF.length());
    Rcpp::NumericVector cdf_perc_in = Rcpp::as<Rcpp::NumericVector>(DDModel_cpp_.RF["CDF"]);
    Rcpp::NumericVector caf_perc_in = Rcpp::as<Rcpp::NumericVector>(DDModel_cpp_.RF["CAF"]);
    for (int i = 0; i< DDModel_cpp_.Parameter.length();++i)
    {
      PAR_v.push_back(0.0);
      PAR_n.push_back(Rcpp::as<std::string>(DDModel_cpp_.Parameter[i]));
    }
    for (int i = 0; i<cdf_perc_in.size();++i)
    {
      RF[0].push_back(cdf_perc_in[i]);
    }
    for (int i = 0; i<caf_perc_in.size();++i)
    {
      RF[1].push_back(caf_perc_in[i]);
    }
    std::vector<double> cdf_buff;
    std::vector<double> caf_buff;
    for (std::size_t p = 0; p< RF[0].size(); ++p)
    {
      cdf_buff.push_back(RF[0][p]);
    }

    for (std::size_t p = 0; p< (RF[1].size()-1); ++p)
    {
      caf_buff.push_back((RF[1][p+1]+RF[1][p])/2);
    }
    CAF.resize(DDModel_cpp_.Conditions.length());
    CDF.resize(DDModel_cpp_.Conditions.length());
    for (int c = 0; c<DDModel_cpp_.Conditions.length();++c)
    {
      CDF[c].resize(cdf_buff.size());
      CAF[c].resize(caf_buff.size());
      for (std::size_t i = 0; i<CDF[c].size();++i)
      {
        CDF[c][i].perc = cdf_buff[i];
        CDF[c][i].cond = DDModel_cpp_.Conditions[c];
      }
      for (std::size_t i = 0; i<CAF[c].size();++i)
      {
        CAF[c][i].perc = caf_buff[i];
        CAF[c][i].cond = DDModel_cpp_.Conditions[c];
      }
    }
    RAW.resize(DDModel_cpp_.Conditions.length());
  }
  DDRep_cpp(DDModel_cpp DDModel_cpp_,Rcpp::List RAW_)
  {
    RF.resize(DDModel_cpp_.RF.length());
    for (int i = 0; i< DDModel_cpp_.Parameter.length();++i)
    {
      PAR_v.push_back(0.0);
      PAR_n.push_back(Rcpp::as<std::string>(DDModel_cpp_.Parameter[i]));
    }
    Rcpp::NumericVector cdf_perc_in = Rcpp::as<Rcpp::NumericVector>(DDModel_cpp_.RF["CDF"]);
    Rcpp::NumericVector caf_perc_in = Rcpp::as<Rcpp::NumericVector>(DDModel_cpp_.RF["CAF"]);
    for (int i = 0; i<cdf_perc_in.size();++i)
    {
      RF[0].push_back(cdf_perc_in[i]);
    }
    for (int i = 0; i<caf_perc_in.size();++i)
    {
      RF[1].push_back(caf_perc_in[i]);
    }
    std::vector<double> cdf_buff;
    std::vector<double> caf_buff;
    for (std::size_t p = 0; p< RF[0].size(); ++p)
    {
      cdf_buff.push_back(RF[0][p]);
    }

    for (std::size_t p = 0; p< (RF[1].size()-1); ++p)
    {
      caf_buff.push_back((RF[1][p+1]+RF[1][p])/2);
    }
    CAF.resize(DDModel_cpp_.Conditions.length());
    CDF.resize(DDModel_cpp_.Conditions.length());
    RAW.resize(DDModel_cpp_.Conditions.length());
    RAW_format RAW_single_buff;
    for (int c = 0; c<DDModel_cpp_.Conditions.length();++c)
    {
      Rcpp::DataFrame RAW_frame_buff = Rcpp::as<Rcpp::DataFrame>(RAW_[c]);
      Rcpp::CharacterVector RAW_cond_buff = RAW_frame_buff["cond"];
      Rcpp::NumericVector RAW_time_buff = RAW_frame_buff["time"];
      Rcpp::NumericVector RAW_resp_buff = RAW_frame_buff["resp"];
      CDF[c].resize(cdf_buff.size());
      CAF[c].resize(caf_buff.size());
      for (std::size_t i = 0; i<CDF[c].size();++i)
      {
        CDF[c][i].perc = cdf_buff[i];
        CDF[c][i].cond = DDModel_cpp_.Conditions[c];
      }
      for (std::size_t i = 0; i<CAF[c].size();++i)
      {
        CAF[c][i].perc = caf_buff[i];
        CAF[c][i].cond = DDModel_cpp_.Conditions[c];
      }
      for (int i = 0; i<RAW_frame_buff.nrow();++i)
      {
        RAW_single_buff.cond = RAW_cond_buff[i];
        RAW_single_buff.resp = RAW_resp_buff[i];
        RAW_single_buff.time = RAW_time_buff[i];
        RAW[c].push_back(RAW_single_buff);
      }
    }
  }
  DDRep_cpp(Rcpp::S4 DDRep_)
  {
    Rcpp::List RAW_buff = Rcpp::as<Rcpp::List>(DDRep_.slot("RAW"));
    Rcpp::List CDF_buff = Rcpp::as<Rcpp::List>(DDRep_.slot("REP"))["CDF"];
    Rcpp::List CAF_buff = Rcpp::as<Rcpp::List>(DDRep_.slot("REP"))["CAF"];
    Rcpp::List RF_buff = Rcpp::as<Rcpp::List>(DDRep_.slot("RF"));
    Rcpp::DataFrame PAR_buff = Rcpp::as<Rcpp::DataFrame>(DDRep_.slot("PAR"));
    Rcpp::CharacterVector PAR_names = PAR_buff.names();
    for (int i = 0; i<PAR_buff.length();++i)
    {
      PAR_v.push_back(Rcpp::as<double>(PAR_buff[i]));
      PAR_n.push_back(Rcpp::as<std::string>(PAR_names[i]));
    }
    Rcpp::CharacterVector cond_name_buff = CDF_buff.names();
    RF.resize(RF_buff.length());
    Rcpp::NumericVector cdf_perc_in = Rcpp::as<Rcpp::NumericVector>(RF_buff["CDF"]);
    Rcpp::NumericVector caf_perc_in = Rcpp::as<Rcpp::NumericVector>(RF_buff["CAF"]);
    for (int i = 0; i<cdf_perc_in.size();++i)
    {
      RF[0].push_back(cdf_perc_in[i]);
    }
    for (int i = 0; i<caf_perc_in.size();++i)
    {
      RF[1].push_back(caf_perc_in[i]);
    }
    std::vector<double> cdf_buff;
    std::vector<double> caf_buff;
    for (std::size_t p = 0; p< RF[0].size(); ++p)
    {
      cdf_buff.push_back(RF[0][p]);
    }

    for (std::size_t p = 0; p< (RF[1].size()-1); ++p)
    {
      caf_buff.push_back((RF[1][p+1]+RF[1][p])/2);
    }
    CAF.resize(cond_name_buff.length());
    CDF.resize(cond_name_buff.length());
    RAW.resize(cond_name_buff.length());
    RAW_format RAW_single_buff;
    for (int c = 0; c<cond_name_buff.length();++c)
    {
      CAF[c].resize(caf_buff.size());
      CDF[c].resize(cdf_buff.size());
      Rcpp::DataFrame CDF_frame_buff = Rcpp::as<Rcpp::DataFrame>(CDF_buff[c]);
      Rcpp::DataFrame CAF_frame_buff = Rcpp::as<Rcpp::DataFrame>(CAF_buff[c]);
      Rcpp::DataFrame RAW_frame_buff = Rcpp::as<Rcpp::DataFrame>(RAW_buff[c]);
      Rcpp::NumericVector cdf_time_buff = CDF_frame_buff["time"];
      Rcpp::NumericVector cdf_N_buff = CDF_frame_buff["N"];
      Rcpp::NumericVector caf_time_buff = CAF_frame_buff["time"];
      Rcpp::NumericVector caf_acc_buff = CAF_frame_buff["acc"];
      Rcpp::NumericVector caf_N_A_buff = CAF_frame_buff["N_A"];
      Rcpp::NumericVector caf_N_B_buff = CAF_frame_buff["N_B"];
      Rcpp::CharacterVector RAW_cond_buff = RAW_frame_buff["cond"];
      Rcpp::NumericVector RAW_time_buff = RAW_frame_buff["time"];
      Rcpp::NumericVector RAW_resp_buff = RAW_frame_buff["resp"];
      for (std::size_t i = 0; i<cdf_buff.size();++i)
      {
        CDF[c][i].perc = cdf_buff[i];
        CDF[c][i].cond = cond_name_buff[c];
        CDF[c][i].time = cdf_time_buff[i];
        CDF[c][i].N = cdf_N_buff[i];
      }
      for (std::size_t i = 0; i<caf_buff.size();++i)
      {
        CAF[c][i].perc = caf_buff[i];
        CAF[c][i].cond = cond_name_buff[c];
        CAF[c][i].time = caf_time_buff[i];
        CAF[c][i].acc = caf_acc_buff[i];
        CAF[c][i].N_A = caf_N_A_buff[i];
        CAF[c][i].N_B = caf_N_B_buff[i];
      }
      for (int i = 0; i<RAW_frame_buff.nrow();++i)
      {
        RAW_single_buff.cond = RAW_cond_buff[i];
        RAW_single_buff.resp = RAW_resp_buff[i];
        RAW_single_buff.time = RAW_time_buff[i];
        RAW[c].push_back(RAW_single_buff);
      }
    }
  }
  //variables
  std::vector<std::vector<CDF_perc>> CDF; //[condition][percentile]
  std::vector<std::vector<CAF_perc>> CAF; //[condition][percentile]
  std::vector<std::vector<RAW_format>> RAW; //[condition][data]
  std::vector<std::vector<double>> RF;  // 0: CDF, 1: CAF
  std::vector<double> PAR_v;// Parameter values for the representation
  std::vector<std::string> PAR_n;// Parameter names for the representation
  //functions
  Rcpp::S4 Convert_to_S4();
  Rcpp::List Convert_FORM_to_List();
};

