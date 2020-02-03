#pragma once
const double pi = 3.14159265;
// class for a normal distribution with mean _mu und standarddeviation _sigma
class NormalDistribution
{
public:
  NormalDistribution(double _mu, double _sigma) : mu(_mu), sigma(_sigma) {}
  inline double pdf(double x);
  double cdf(double x);
private:
  double mu;
  double sigma;
};

inline double NormalDistribution::pdf(double x)
{  //Standard gaussian function
  return exp(-1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * pi));
}
double NormalDistribution::cdf(double x)
{
  // cdf of standard gaussian function
  return 0.5 * (1 + std::erf((x - mu) / (sigma * sqrt(2.))));
}
