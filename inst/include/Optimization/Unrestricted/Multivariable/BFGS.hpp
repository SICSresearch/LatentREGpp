// BFGS Algorithm
// Adapted from: Compact Numerical Methods for Computers. J.C. Nash
// 
// acctol     acceptable point of tolerance
// reltest    to check equality of parameter
// accpoint   to indicate an acceptable point
// B          an aproximation of inverse matrix
// c          to store last gradient
// count      to check for parameter equality
// D1, D2     temporaly working storage
// f_eval     temporary function value
// funcount   count of function evaluations
// gradient   to hold the gradient
// gradproj   gradient proyection on search vector
// ilast      record last step at which B was initialized to a unit matrix
// s          inner product accumulator
// steplength linear search steplength
// t          to store working vector for line search

#ifndef BFGS_HPP
#define BFGS_HPP

#include <Optimization/Optim.hpp>
#include <cmath>
#include <iostream>

namespace spgo
{

enum INDEX{FUNCTION, GRADIENT, HESSIAN};

template <typename Function, typename Parameter, typename Info>
class BFGS : public Optim<Function, Parameter, Info>
{
  public:
    BFGS()
    {
      setControlConstants(1e30, 0.02, 0.0001, 10.0, 0.00001, 1e-8);

      maxiter = 15;
      size    = -1;
    }

    BFGS(double INF, double stepredn, double acctol, double reltest, double abstol, double reltol)
    {
      setControlConstants(INF, stepredn, acctol, reltest, abstol, reltol);
      maxiter = 15;
    }

    void setControlConstants(double INF, double stepredn, double acctol, double reltest, double abstol, double reltol)
    {
      this->INF      = INF;
      this->stepredn = stepredn;
      this->acctol   = acctol;
      this->reltest  = reltest;
      this->abstol   = abstol;
      this->reltol   = reltol;
      iterations     = 0;
    }

    void setParameterSize(int size) { this->size = size; }

    void allocateMemory(unsigned int n)
    {
      t = new double[n];
      X = new double[n];
      c = new double[n];
      B = new double*[n];

      for(unsigned int i = 0; i < n; i++) { B[i] = new double[n]; }

      this->n = n;
    }

    // f[0] = Function to optimize
    // f[1] = Gradient of the function
    // p = Parameters to change
    // info = Static parameters, other things that function may use.
    int run(std::vector<Function> f, Parameter &p, Info info)
    {
      allocateMemory(getSize(p));

      this->f = f[FUNCTION];
      this->g = f[GRADIENT];
      f_eval  = this->f(p, info)[0];
    
      if (f_eval >= INF) return (2);

      fmin     = f_eval;
      funcount = gradcount = 1;
      gradient = this->g(p, info);
      ilast    = gradcount;

      iterations++;

      iterate(p, info);

      for(unsigned int i = 0; i < n; i++) { delete [] B[i]; }

      delete[] B;
      delete[] t;
      delete[] c;
      delete[] X;

      if (iterations < maxiter) { return (0); }

      return (3);
    }

    ~BFGS(){};

  private:
    Parameter gradient;
    Function  f,
              g;
    double  abstol, reltol, acctol,
            INF, reltest, stepredn,
            steplength,
            gradproj, f_eval, fmin, s, D1, D2,
          * t, *c, *X, **B;
    int iterations, ilast, n,
        i, j, count, gradcount,
        funcount, maxiter, size;
    bool enough, accpoint;

    int getSize(std::vector<double> vect)
    {
      return vect.size();
    }

    int getSize(double * vect)
    {
      return this->size;
    }

    void initializeMatrix()
    {
      for (i = 0; i < n; i++)
      {
        for (j = 0; j < i; j++) { B[i][j] = 0.0; }
        B[i][i] = 1.0;
      }
    }

    void gradientProjection()
    {
      // To save the tT*g inner product
      gradproj = 0.0;

      for (i = 0; i < n; i++)
      {
        // To acumulate element of B*g
        s = 0.0;

        for (j = 0; j <= i; j++) { s -= B[i][j] * gradient[j]; }
        for (j = i + 1; j < n; j++) { s -= B[j][i] * gradient[j]; }

        t[i]      = s;
        gradproj += s * gradient[i];
      }
    }

    void linearSearch(Parameter &p, Info info)
    {
      // Always try full step first
      steplength = 1.0;
      // Don't have a good point yet
      accpoint = false;

      do
      {
        // To count unchanged parameters
        count = 0;

        for (i = 0; i < n; i++)
        {
          p[i] = X[i] + steplength * t[i];

          if (reltest + X[i] == reltest + p[i])
          {
            count++;
          }
        }

        if (count < n)
        {
          // Function calculation
          f_eval = f(p, info)[0];
          funcount++;

          accpoint = (f_eval < INF) && (f_eval <= fmin + gradproj * steplength * acctol);
          if (!accpoint) steplength *= stepredn;
        }
      } while (!(count == n || accpoint));
    }

    void progress(Parameter &p, Info info)
    {
      if (count < n)
      {
        // Save function value
        fmin = f_eval;

        gradient = g(p, info);

        gradcount++;
        iterations++;

        D1 = 0.0;

        for (i = 0; i < n; i++)
        {
          // To compute vector y
          t[i] = steplength * t[i];
          c[i] = gradient[i] - c[i];
          // To compute inner product t*y
          D1   += t[i] * c[i];
        }

        if (D1 > 0)
        {
          // Test if update is possible
          D2 = 0.0;
          for (i = 0; i < n; i++)
          {
            s = 0.0;
            for (j = 0; j <= i; j++) { s += B[i][j] * c[j]; }
            for (j = i + 1; j < n; j++) { s += B[j][i] * c[j]; }

            X[i] = s;
            D2 += s * c[i];
          }

          D2 = 1.0 + D2 / D1;

          for (i = 0; i < n; i++)
          {
            for (j = 0; j <= i; j++)
              {
                B[i][j] -= (X[i]*t[j] + t[i]*X[j] - D2*t[i]*t[j]) / D1;
              }
          }
        } else ilast = gradcount;
      }
      else
      {
        if (ilast < gradcount)
        {
          count = 0;
          ilast = gradcount;
        }
      }
    }

    void iterate(Parameter &p, Info info)
    {
      do
      {
        if (ilast == gradcount) initializeMatrix();

        for (i = 0; i < n; i++)
        {
          // Save the best parameters
          X[i] = p[i];
          // Save gradient
          c[i] = gradient[i];
        }

        gradientProjection();

        // Test for descent direction
        if (gradproj < 0.0)
        {
          linearSearch(p, info);

          enough = (f_eval > abstol) && fabs(f_eval - fmin) > reltol * (fabs(fmin) + reltol);

          if (!enough)
          {
            count = n;
            fmin  = f_eval;
          }

          progress(p, info);
        }
        else
        {
          count = 0;

          if (ilast == gradcount) { count = n; }
          else { ilast = gradcount; }
        }

        if (iterations >= maxiter) { break; }

        if (gradcount - ilast > 2 * n) { ilast = gradcount; }
      } while (count != n || ilast != gradcount);
    }
};

}

#endif
