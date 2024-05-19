#include <stdio.h>
#include <math.h>
#include "parm.h"

double swsh (double , double, double);
double fact (double);

double fact (double n){
double c;
double fac = 1.0; 
for (c = 1.0; c <= n; c++)
  {  fac = fac * c; }
return fac;
}

double swsh (double ell, double ml, double the){
double tmp, c1, c2, t, norm;

tmp = 0.0;

for (t = 0; t <= ell; t++){
norm = fact(ell)/fact(t)/fact(ell-t)*fact(ell+t)/fact(t)/fact(ell);
tmp = tmp + norm*pow((the-1.0)/2.0,t);
}
tmp = tmp * (2.0*ell+1.0)/2.0;

/*
if ((int)ell == 0)
tmp = 1.0;

if ((int)ell == 1)
tmp = the;

if ((int)ell == 2)
tmp = 1.5*the*the-0.5;

if ((int)ell == 3)
tmp = 2.5*the*the*the-1.5*the;

if ((int)ell == 4)
tmp = (35.0*pow(the,4)-30.0*the*the+3.0)/8.0;
*/

return tmp;
}

double
project (double *qq, double l, double m)
{

  int i, j;
  double stheta, ctheta, dtheta, hlm, th;
  dtheta = -2.0 / (M-1);
  hlm = 0.0;

  for (i = 0; i < (M-1) / 8 ; i++)
    {

      th = 1.0 + 8.0 * dtheta * (double) i;
      hlm += 989.0 * qq[8*i] * dtheta * swsh(l, m, th);

      th += dtheta;
      hlm += 5888.0 * qq[8*i+1] * dtheta * swsh(l, m, th);

      th += dtheta;
      hlm += (-928.0) * qq[8*i+2] * dtheta * swsh(l, m, th);

      th += dtheta;
      hlm += 10496.0 * qq[8*i+3] * dtheta * swsh(l, m, th);

      th += dtheta;
      hlm += (-4540.0) * qq[8*i+4] * dtheta * swsh(l, m, th);

      th += dtheta;
      hlm +=  10496.0 * qq[8*i+5] * dtheta * swsh(l, m, th);

      th += dtheta;
      hlm += (-928.0) * qq[8*i+6] * dtheta * swsh(l, m, th);

      th += dtheta;
      hlm += 5888.0 * qq[8*i+7] * dtheta * swsh(l, m, th);

      th += dtheta;
      hlm += 989.0 * qq[8*i+8] * dtheta * swsh(l, m, th);

    }

  hlm = 4.0 * hlm / 14175.0;
  return hlm;

}
