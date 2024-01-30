#include "FormFactor.h"
#include <cstring>



void Helm::Setsval(double s) {

  sval = s;

}

double Helm::Getsval() { return sval;}


double Helm::FFval(double Q)
{

  double ff = 1;
  //double R = 1.14*pow(A,1./3.);
  double R = 1.2*pow(A,1./3.);

  // This is for varying Rn but keeping sval fixed (will distort the shape
  // but in practice very small difference )

  // double Rnorig = sqrt(3./5.*pow(R,2)+3*sval*sval);
  //double Rmod = sqrt(5./3.*(pow(Rnorig*Rfac,2)-3*sval*sval));
  //double qR = Q*Rmod;

  // This scales the radius
  Q *= Rfac;

   double qR = Q*R;

  ff= (3*(sin(qR)/(qR*qR)-cos(qR)/qR)/(qR))*exp(-1.*Q*Q*sval*sval/2.);
  //  ff2= pow(3*(sin(qR)/(qR*qR)-cos(qR)/qR)/(qR),2)*exp(-1.*Q*Q*sval*sval);
  if (isnan(ff)) {ff=1.;}

  return ff;

}

Helm::Helm(double s, double Rf) {

  sval = s;
  Rfac = Rf;
  strcpy(fftype,"helm");

}
