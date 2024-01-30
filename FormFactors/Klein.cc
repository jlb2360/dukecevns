#include "FormFactor.h"
#include <cstring>


Klein::Klein(double ak, double sf, double Rf) {

  akval = ak;
  skinfac = sf;
  Rfac = Rf;
  strcpy(fftype,"klein");

}

void Klein::Setakval(double ak) {

  akval = ak;

}

double Klein::Getskinfac() { return skinfac;}

void Klein::Setskinfac(double sf) {

  skinfac = sf;

}

double Klein::Getakval() { return akval;}



double Klein::FFval(double Q)
{

  double ff = 1;

  // Gutlein... probably wrong
  // double R2 = 1.14*pow(A,1./3.);

  // Incorrect way to add skin
  //double R2 = 1.2*pow(A,1./3.)+skinfac*1.01*(double(A)-2.*Z)/double(A);

  // Adding a skin
  double skindelta = skinfac*1.01*(double(A)-2.*Z)/double(A);
  double R2 = 1.2*pow(A,1./3.);
  if (skindelta != 0) {

    R2 = sqrt(R2*R2+ 2*sqrt(15)/3*sqrt(R2*R2+10*akval*akval)*skindelta + 5*skindelta*skindelta/3);
  }

  //double Ravg = sqrt(3*R2*R2/5+6*akval*akval);
  //  std::cout << "A, Z, R2, skindelta, Ravg: "<<A<<" "<<Z<<" "<<R2<<" "<<skindelta<<" "<<Ravg<<std::endl;

      // This scales the radius by Rfac
  Q *= Rfac;

  double qR = Q*R2;


  ff= (3*(sin(qR)/(qR*qR)-cos(qR)/qR)/(qR))*(1./(1+akval*akval*Q*Q));
  //  ff2= pow(3*(sin(qR)/(qR*qR)-cos(qR)/qR)/(qR),2)*pow(1./(1+akval*akval*Q*Q),2);


  if (isnan(ff)) {ff=1.;}

  return ff;

}

/////////
