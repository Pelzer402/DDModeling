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
          sum += EVAL[i].Rep.PAR_v[j];
        }
      }
      EVAL[Model.Parameter.size()+1].Rep.PAR_v[j] = sum;
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
        rtol += std::pow((EVAL[i].Fit - (EVAL[Model.Parameter.size()+1].Fit / Model.Parameter.size())), 2.0) / Model.Parameter.size();
      }
    }
    if (rtol < 0.000000000001)
    {
      // Simplex terminated
      std::swap(EVAL[0], EVAL[ilo]);
      break;
    }

    if (nfunk >= 120)
    {
      // Number of function calls exceeded
      std::swap(EVAL[0], EVAL[ilo]);
      break;
    }
    if (Sr_count >3 )
    {
      // Shrink count exceeded
      std::swap(EVAL[0], EVAL[ilo]);
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
        EVAL[ihi].Rep.PAR_v = EVAL[Model.Parameter.size()+3].Rep.PAR_v;
      }
      else
      {
        EVAL[ihi].Fit = EVAL[Model.Parameter.size()+2].Fit;
        EVAL[ihi].Rep.PAR_v = EVAL[Model.Parameter.size()+2].Rep.PAR_v;
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
                EVAL[i].Rep.PAR_v[j] = (double)(0.5*(EVAL[i].Rep.PAR_v[j] + EVAL[ilo].Rep.PAR_v[j]));
              }
              Simulate_and_Fit(i);
            }
          }
          Sr_count++;
        }
        else
        {
          EVAL[ihi].Fit = EVAL[Model.Parameter.size()+3].Fit;
          EVAL[ihi].Rep.PAR_v = EVAL[Model.Parameter.size()+3].Rep.PAR_v;
        }
      }
      else
      {
        EVAL[ihi].Fit = EVAL[Model.Parameter.size()+2].Fit;
        EVAL[ihi].Rep.PAR_v = EVAL[Model.Parameter.size()+2].Rep.PAR_v;

        SIMPLEX_TransformSimplex_star2_beta(ihi, 0.5);

        if (EVAL[Model.Parameter.size()+3].Fit > EVAL[ihi].Fit)
        {
          for (i = 0; i < Model.Parameter.size()+1; ++i)
          {
            if (i != ilo)
            {
              for (j = 0; j < Model.Parameter.size(); ++j)
              {
                EVAL[i].Rep.PAR_v[j] = (double)(0.5*(EVAL[i].Rep.PAR_v[j] + EVAL[ilo].Rep.PAR_v[j]));
              }
              Simulate_and_Fit(i);
            }
          }
          Sr_count++;
        }
        else
        {
          EVAL[ihi].Fit = EVAL[Model.Parameter.size()+3].Fit;
          EVAL[ihi].Rep.PAR_v = EVAL[Model.Parameter.size()+3].Rep.PAR_v;
        }
      }
    }
    else
    {
      EVAL[ihi].Fit = EVAL[Model.Parameter.size()+2].Fit;
      EVAL[ihi].Rep.PAR_v = EVAL[Model.Parameter.size()+2].Rep.PAR_v;
    }

    nfunk++;
  }
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
    EVAL[Model.Parameter.size()+2].Rep.PAR_v[j] = (EVAL[Model.Parameter.size()+1].Rep.PAR_v[j]/Model.Parameter.size()) * fac1 + EVAL[ihi].Rep.PAR_v[j] * fac2;
    /*
    if (EVAL[Model.Parameter.size()+2].Rep.PAR_v[j] < Model.DM(1,j))
    {
      EVAL[Model.Parameter.size()+2].Rep.PAR_v[j] =  Model.DM(1,j);
    }
    if (EVAL[Model.Parameter.size()+2].Rep.PAR_v[j] > Model.DM(0,j))
    {
      EVAL[Model.Parameter.size()+2].Rep.PAR_v[j] = Model.DM(0,j);
    }
    */
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
    EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] = (EVAL[Model.Parameter.size()+1].Rep.PAR_v[j] / Model.Parameter.size()) * fac1 + EVAL[ihi].Rep.PAR_v[j] * fac2;
    /*
    if (EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] <  Model.DM(1,j))
    {
      EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] =  Model.DM(1,j);
    }
    if (EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] > Model.DM(0,j))
    {
      EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] = Model.DM(0,j);
    }
    */
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
    EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] = (EVAL[Model.Parameter.size()+1].Rep.PAR_v[j] / Model.Parameter.size()) * fac1 + EVAL[ihi].Rep.PAR_v[j] * fac2;
    /*
    if (EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] <  Model.DM(1,j))
    {
      EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] =  Model.DM(1,j);
    }
    if (EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] >  Model.DM(0,j))
    {
      EVAL[Model.Parameter.size()+3].Rep.PAR_v[j] = Model.DM(0,j);
    }
    */
  }

  Simulate_and_Fit(Model.Parameter.size()+3);		    //Evaluate the function at the trial point.

}
