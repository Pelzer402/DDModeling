#include "Simulation.h"
#include <math.h>
#include <algorithm>

using namespace std;

//Simplex:
void Simulation::SIMPLEX()
{
  int i,j;				// loop var
  int ihi;				// highes(worst) point index
  int ilo;				// lowest(best) point index
  int inhi;				// second highest(worst) point index
  double rtol;			// Break Critria
  double sum;				// tmp var for summation
  int   nfunk;			// Number of function calls
  int  Sr_count = 0;		//Nr. auf shrinks performed

  nfunk = 0;
  for (;;)
  {

    ilo = 0;

    // First we must determine which point is the highest(worst), next - highest, and lowest (best), by looping over the points in the simplex.
    if (EVAL[0].Fit > EVAL[1].Fit)
    {
      ihi = 0;
      inhi = 1;
    }
    else
    {
      ihi = 1;
      inhi = 0;
    }
    for (i = 0; i < Model.Parameter.size()+1; ++i)
    {
      if (EVAL[i].Fit <= EVAL[ilo].Fit)
      {
        ilo = i;
      }
      if (EVAL[i].Fit  > EVAL[ihi].Fit )
      {
        inhi = ihi;
        ihi = i;
      }
      else if (EVAL[i].Fit  > EVAL[inhi].Fit  && i != ihi)
      {
        inhi = i;
      }
    }

    // calculate Sum over parameter and chi2 (except i==ihi)
    for (j = 0; j < Model.Parameter.size(); ++j)
    {
      sum = 0.0;
      for (i = 0; i < Model.Parameter.size()+1; ++i)
      {
        if (i != ihi)
        {
          sum += EVAL[i].Parameter[j];
        }
      }
      EVAL[Model.Parameter.size()+1].Parameter[j] = sum;
    }

    EVAL[Model.Parameter.size()+1].Fit = 0.0;
    for (j = 0; j < Model.Parameter.size()+1; ++j)
    {
      if (i != ihi)
      {
        EVAL[Model.Parameter.size()+1].Fit += EVAL[j].Fit;
      }

    }
    //Compute the fractional range from highest to lowest
    rtol = 0.0;
    for (i = 0; i < Model.Parameter.size()+1; ++i)
    {
      if (i != ihi)
      {
        rtol += pow((EVAL[i].Fit - (EVAL[Model.Parameter.size()+1].Fit / Model.Parameter.size())), 2.0) / Model.Parameter.size();
      }
    }

    // Test display ??

   // Rcpp::Rcout << std::endl;
   // Rcpp::Rcout << "Funktion Call:  " << nfunk << ", Fit:  " << rtol << endl << endl;
   // std::setprecision(5);
    int wd = 15;
    for (int ok = 0; ok < Model.Parameter.size()+1; ++ok)
    {
     // Rcpp::Rcout<< std::setw(5) << std::left << ok << std::setw(wd) << std::left << EVAL[ok].Fit;
      for (int it = 0; it < Model.Parameter.size(); ++it)
      {
       // Rcpp::Rcout << std::setw(wd) << std::left << EVAL[ok].Parameter[it];
      }
      //Rcpp::Rcout << std::endl;
    }
   // Rcpp::Rcout << endl << "---------- DETAILS ----------" << endl << endl;
    //Rcpp::Rcout << "Nr. of Space Shrink:  " << Sr_count << endl << endl;
    //Rcpp::Rcout << "B Chi2:  " << ilo << "   " <<EVAL[ilo].Fit << endl << endl;
    //Rcpp::Rcout << "W Chi2:  " << ihi << "   " <<EVAL[ihi].Fit << endl << endl;
   // Rcpp::Rcout << "------------------------------------------------------------" << endl << endl;

    // Test Display ??

    if (rtol < 0.000000000001)
    {
      //std::cout << " Simplex terminated!" << endl;
      //If returning, put best point and value in slot 0.
      std::swap(EVAL[0].Fit, EVAL[ilo].Fit);
      std::swap(EVAL[0].Parameter, EVAL[ilo].Parameter);
      //std::swap(SIMPLEX_FINAL[0], SIMPLEX_FINAL[ilo]);
      //std::swap(SIMPLEX_FINAL[0], SIMPLEX_FINAL[ilo]);
      break;
    }

    if (nfunk >= 120)
    {
      //std::cout << "Number of function calls exceeded!" << endl;
      std::swap(EVAL[0].Fit, EVAL[ilo].Fit);
      std::swap(EVAL[0].Parameter, EVAL[ilo].Parameter);
      //std::swap(SIMPLEX_FINAL[0], SIMPLEX_FINAL[ilo]);
      //std::swap(SIMPLEX_FINAL[0], SIMPLEX_FINAL[ilo]);
      break;
    }
    if (Sr_count >3 )
    {
      //std::cout << "Shrink count exceeded!" << endl;
      std::swap(EVAL[0].Fit, EVAL[ilo].Fit);
      std::swap(EVAL[0].Parameter, EVAL[ilo].Parameter);
      //std::swap(SIMPLEX_FINAL[0], SIMPLEX_FINAL[ilo]);
      //std::swap(SIMPLEX_FINAL[0], SIMPLEX_FINAL[ilo]);
      break;
    }
    //Begin a new iteration.


    //First reflexion by a factor alpha=1 through the face of the simplex across from the high point, i.e. reflect the simplex from the high point.
    SIMPLEX_TransformSimplex( ihi, 1.0);

    if (EVAL[Model.Parameter.size()+2].Fit <= EVAL[ilo].Fit)
    {
      // Gives a result better than the best point, so try an additional Expansion by a factor gamma=2
      SIMPLEX_TransformSimplex_star2_gamma(Model.Parameter.size()+2, 2.0);

      if (EVAL[Model.Parameter.size()+3].Fit <= EVAL[ilo].Fit)
      {
        EVAL[ihi].Fit = EVAL[Model.Parameter.size()+3].Fit;
        EVAL[ihi].Parameter = EVAL[Model.Parameter.size()+3].Parameter;
      }
      else
      {
        EVAL[ihi].Fit = EVAL[Model.Parameter.size()+2].Fit;
        EVAL[ihi].Parameter = EVAL[Model.Parameter.size()+2].Parameter;
      }
    }
    else if (EVAL[Model.Parameter.size()+2].Fit >= EVAL[inhi].Fit)
    {
      // The reflected point is worse than the second - highest, so look for an intermediate lower point, i.e., do a one - dimensional contraction.
      if (EVAL[Model.Parameter.size()+2].Fit > EVAL[ihi].Fit)
      {
        SIMPLEX_TransformSimplex_star2_beta(ihi, 0.5);

        if (EVAL[Model.Parameter.size()+3].Fit > EVAL[ihi].Fit)
        {
          for (i = 0; i < Model.Parameter.size()+1; ++i)
          {
            if (i != ilo)
            {
              for (j = 0; j < Model.Parameter.size(); ++j)
              {
                EVAL[i].Parameter[j] = (double)(0.5*(EVAL[i].Parameter[j] + EVAL[ilo].Parameter[j]));
              }
              Simulate_and_Fit(i);
            }
          }
          Sr_count++;
        }
        else
        {
          EVAL[ihi].Fit = EVAL[Model.Parameter.size()+3].Fit;
          EVAL[ihi].Parameter = EVAL[Model.Parameter.size()+3].Parameter;
        }
      }
      else
      {
        EVAL[ihi].Fit = EVAL[Model.Parameter.size()+2].Fit;
        EVAL[ihi].Parameter = EVAL[Model.Parameter.size()+2].Parameter;

        SIMPLEX_TransformSimplex_star2_beta(ihi, 0.5);

        if (EVAL[Model.Parameter.size()+3].Fit > EVAL[ihi].Fit)
        {
          for (i = 0; i < Model.Parameter.size()+1; ++i)
          {
            if (i != ilo)
            {
              for (j = 0; j < Model.Parameter.size(); ++j)
              {
                EVAL[i].Parameter[j] = (double)(0.5*(EVAL[i].Parameter[j] + EVAL[ilo].Parameter[j]));
              }
              Simulate_and_Fit(i);
            }
          }
          Sr_count++;
        }
        else
        {
          EVAL[ihi].Fit = EVAL[Model.Parameter.size()+3].Fit;
          EVAL[ihi].Parameter = EVAL[Model.Parameter.size()+3].Parameter;
        }
      }
    }
    else
    {
      EVAL[ihi].Fit = EVAL[Model.Parameter.size()+2].Fit;
      EVAL[ihi].Parameter = EVAL[Model.Parameter.size()+2].Parameter;
    }

    nfunk++;
  }// Go back for the test of doneness , iteration.
  //Set.E_BUFFER.push_back(SIMPLEX_FINAL[0]);
  //std::cout << " Funktion Calls:  " << nfunk << ", Final Chi2:  " << Set.FitCrit[0] << endl;
  //EVAL_BUFFER tmp;
  //tmp.FIT = Set.FitCrit[0];
  //tmp.PAR = Set.PAR_MAT[0];
  //Set.E_BUFFER.push_back(tmp);
  RESULT.push_back(EVAL[0]);
}

//TransformSimplex:
/*
 Input:
 ihi		-> Index of the highest value in point cluster
 fac		-> Factor of the transformation
 mod     -> Model Modifikator
 */
void Simulation::SIMPLEX_TransformSimplex(int ihi, double fac)
{
  int j;
  double fac1;		// factor for mean transformation
  double fac2;		// facotr for *   transformation

  fac1 = (1.0 + fac);
  fac2 = -fac;
  for (j = 0; j < Model.Parameter.size(); ++j)
  {
    EVAL[Model.Parameter.size()+2].Parameter[j] = (EVAL[Model.Parameter.size()+1].Parameter[j]/Model.Parameter.size()) * fac1 + EVAL[ihi].Parameter[j] * fac2;
    if (EVAL[Model.Parameter.size()+2].Parameter[j] < Model.DM(1,j))
    {
      EVAL[Model.Parameter.size()+2].Parameter[j] =  Model.DM(1,j);
    }
    if (EVAL[Model.Parameter.size()+2].Parameter[j] > Model.DM(0,j))
    {
      EVAL[Model.Parameter.size()+2].Parameter[j] = Model.DM(0,j);
    }
  }

  Simulate_and_Fit(Model.Parameter.size()+2);    //Evaluate the function at the trial point.
}
void Simulation::SIMPLEX_TransformSimplex_star2_beta(int ihi, double fac)
{
  int j;
  double fac1;		// factor for mean transformation
  double fac2;		// facotr for *   transformation

  fac1 = (1.0 - fac);
  fac2 = fac;
  for (j = 0; j < Model.Parameter.size(); ++j)
  {
    EVAL[Model.Parameter.size()+3].Parameter[j] = (EVAL[Model.Parameter.size()+1].Parameter[j] / Model.Parameter.size()) * fac1 + EVAL[ihi].Parameter[j] * fac2;
    if (EVAL[Model.Parameter.size()+3].Parameter[j] <  Model.DM(1,j))
    {
      EVAL[Model.Parameter.size()+3].Parameter[j] =  Model.DM(1,j);
    }
    if (EVAL[Model.Parameter.size()+3].Parameter[j] > Model.DM(0,j))
    {
      EVAL[Model.Parameter.size()+3].Parameter[j] = Model.DM(0,j);
    }
  }

  Simulate_and_Fit(Model.Parameter.size()+3);    //Evaluate the function at the trial point.

}
void Simulation::SIMPLEX_TransformSimplex_star2_gamma(int ihi, double fac)
{
  int j;
  double fac1;		// factor for mean transformation
  double fac2;		// facotr for *   transformation

  fac1 = -fac;
  fac2 = fac+1.0;
  for (j = 0; j < Model.Parameter.size(); ++j)
  {
    EVAL[Model.Parameter.size()+3].Parameter[j] = (EVAL[Model.Parameter.size()+1].Parameter[j] / Model.Parameter.size()) * fac1 + EVAL[ihi].Parameter[j] * fac2;
    if (EVAL[Model.Parameter.size()+3].Parameter[j] <  Model.DM(1,j))
    {
      EVAL[Model.Parameter.size()+3].Parameter[j] =  Model.DM(1,j);
    }
    if (EVAL[Model.Parameter.size()+3].Parameter[j] >  Model.DM(0,j))
    {
      EVAL[Model.Parameter.size()+3].Parameter[j] = Model.DM(0,j);
    }
  }

  Simulate_and_Fit(Model.Parameter.size()+3);		    //Evaluate the function at the trial point.

}
