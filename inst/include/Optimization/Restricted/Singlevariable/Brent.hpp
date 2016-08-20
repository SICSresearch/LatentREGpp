// Brent golden section Algorithm
// Adapted from:  Brent, R.P. 1973, Algorithms for Minimization without Derivatives (Englewood Cliffs, NJ: Prentice-Hall)
// 

#ifndef BRENT_HPP
#define BRENT_HPP

#include <Optimization/Optim.hpp>
#include <cmath>
#include <float.h>

using namespace std;

namespace spgo
{

  template <typename Function, typename Parameter, typename Info>
  class Brent : public Optim<Function, Parameter, Info>
  {
    public:

      Brent() { setControlConstants(0.0001220703, sqrt(DBL_EPSILON), DBL_EPSILON + 1.); }

      Brent(double tol, double eps, double tol1) { setControlConstants(tol, eps, tol1); }

      ~Brent() {}

      void setControlConstants(double tol, double eps, double tol1)
      {
        this->tol1 = tol1;
        this->eps  = eps;
        this->tol  = tol;
      }

      void setBounds(double lower, double upper)
      {
        this->a = lower;
        this->b = upper;
      }

      // f[0] = Function to optimize
      // p = Parameter to change
      int run(vector<Function> f, Parameter &p, Info info)
      {
        double c;

        c    = (3. - sqrt(5.)) * .5;
        v    = a + c * (b - a);
        w    = v;
        p    = v;
        d    = 0.;
        e    = 0.;
        fx   = f[0](p, info);
        fv   = fx;
        fw   = fx;
        tol3 = tol / 3.;

        for (;;)
        {
          xm   = (a + b) * .5;
          tol1 = eps * fabs(p) + tol3;
          t2   = tol1 * 2.;

          // check stopping criterion
          if (fabs(p - xm) <= t2 - (b - a) * .5) { break; }
                
          p1 = 0.;
          q  = 0.;
          r  = 0.;
                
          if (fabs(e) > tol1)
          {
            // fit parabola
            r  = (p - w) * (fx - fv);
            q  = (p - v) * (fx - fw);
            p1 = (p - v) * q - (p - w) * r;
            q  = (q - r) * 2.;

            if (q > 0.) { p1 = -p1; }
            else { q = -q; }

            r = e;
            e = d;
          }
          if (fabs(p1) >= fabs(q * .5 * r) || p1 <= q * (a - p) || p1 >= q * (b - p))
          {
            // a golden-section step
            if (p < xm) { e = b - p; }
            else { e = a - p; }
            d = c * e;
          }
          else
          {
            // a parabolic-interpolation step
            d = p1 / q;
            u = p + d;

            // f must not be evaluated too close to ax or bx
            if (u - a < t2 || b - u < t2)
            {
              d = tol1;
              if (p >= xm) { d = -d; }
            }
          }

          // f must not be evaluated too close to x
          if (fabs(d) >= tol1) { u = p + d; }
          else if (d > 0.) { u = p + tol1; }
          else { u = p - tol1; }

          fu = f[0](u, info);

          // update a, b, v, w, and x
          if (fu <= fx)
          {
            if (u < p) { b = p; }
            else { a = p; }

            v  = w;
            w  = p;
            p  = u;
            fv = fw;
            fw = fx;
            fx = fu;
          }
          else
          {
            if (u < p) { a = u; }
            else { b = u; }
            if (fu <= fw || w == p)
            {
              v = w;
              fv = fw;
              w = u;
              fw = fu;
            }
            else if (fu <= fv || v == p || v == w)
            {
              v = u;
              fv = fu;
            }
          }
        }

        return 0;
      }

    private:
      /* Local variables */
      double a, b, d, e, p1, q, r, u, v, w;
      double t2, fu, fv, fw, fx, xm, eps, tol, tol1, tol3;
  };

}

#endif