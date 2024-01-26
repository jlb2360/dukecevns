#include <cstdlib>
#include <cstring>
#include <iostream>
#include <math.h>
#include "NuFlux.h"



double Reactor::fluxval(double Enu, int flavor, double ebinsize)
{

 // Polynomials from Mueller 2011.  Gives flux in per MeV per fission
  // Output from this function is for the specified ebinsize

  // parent array gives the relative contributions from each fissioning parent, which will actually vary with time

  double alpha235U[6] = {3.217, -3.111, 1.395, -3.690e-1, 4.445e-2, -2.053e-3};
  double alpha238U[6] = {4.833e-1, 1.927e-1, -1.283e-1, -6.762e-3, 2.233e-3, -1.536e-4};

  double alpha239Pu[6] = {6.413,-7.432,3.535,-0.8820,0.1025,-.004550};
  double alpha241Pu[6] = {3.251,-3.204,1.428,-.3675,0.04254,-0.001896};

  double nuen = Enu;

  // The polynomials are only good down to about 2 MeV.  The following is an extrapolation at low energy, which is probably better than using the anomalous polynomial values there, but for which information should be improved.

  // These are extrapolated exponential parameters below 1.8 MeV
  double a_235U = -0.7256;
  double b_235U = 1.7112;
  double a_238U = -0.3167;
  double b_238U = 0.9961;
  double a_239Pu = -1.0380;
  double b_239Pu = 2.2021;
  double a_241Pu = -1.0380;
  double b_241Pu = 2.2021;

  // nuebar only
  if (flavor != -1) {
        return 0;
  }


  double flux235 = 0.;
  double flux238 = 0.;
  double flux239 = 0.;
  double flux241 = 0.;

  if (nuen<1.8) {


    flux235 = exp(b_235U)*exp(a_235U*nuen)*ebinsize;
    flux238 = exp(b_238U)*exp(a_238U*nuen)*ebinsize;
    flux239 = exp(b_239Pu)*exp(a_239Pu*nuen)*ebinsize;
    flux241 = exp(b_241Pu)*exp(a_241Pu*nuen)*ebinsize;

  } else {
    int numterms = 6;
    int p;

    for (p=1; p<=numterms; p++) {
      flux235 += alpha235U[p-1]*pow(nuen,p-1);
      flux238 += alpha238U[p-1]*pow(nuen,p-1);
      flux239 += alpha239Pu[p-1]*pow(nuen,p-1);
      flux241 += alpha241Pu[p-1]*pow(nuen,p-1);
    }

    flux235 = exp(flux235)*ebinsize;
    flux238 = exp(flux238)*ebinsize;
    flux239 = exp(flux239)*ebinsize;
    flux241 = exp(flux241)*ebinsize;


  }



  double fluxtot = flux235*parentfrac[0]+flux238*parentfrac[1]+flux239*parentfrac[2]+flux241*parentfrac[3];


  return fluxtot*norm;


}


double Reactor::maxEnu()
{

  // Return the maximum energy in MeV

  double maxEnu = 8.;

  return maxEnu;

}

void Reactor::SetParentFrac(double* frac) {

  int i;
  for (i=0;i<numparent;i++) {
    parentfrac[i] = frac[i];
  }

}

double* Reactor::GetParentFrac() {
  return parentfrac;
}

/////////////



double Monochromatic::fluxval(double Enu, int flavor, double ebinsize)
{

  double flux;
  if (flavor == monoflavor) {
    flux = 1.;
  } else {
    flux = 0.;
  }

  if (fabs(Enu-monoenergy)<ebinsize/2.) {
    flux *= norm;
  } else {
    flux = 0.;
  }

  return flux;

}


double Monochromatic::maxEnu()
{

  // Return the maximum energy in MeV

  double maxEnu = monoenergy;

  return maxEnu;

}

void Monochromatic::SetFlavor(int flavor) {

  monoflavor = flavor;

}

 int Monochromatic::GetFlavor() {
  return monoflavor;
}



void Monochromatic::SetEnergy(double energy) {

  monoenergy = energy;

}

double Monochromatic::GetEnergy() {
  return monoenergy;
}


/////////

void NumericalFlux::ReadFluxFile()
{

  // Input file has energies in GeV
  double enu,nue,numu,nutau,nuebar,numubar,nutaubar;
  std::ifstream fluxfile;
  std::string fluxfilename = filename;
  fluxfile.open(fluxfilename.c_str());
  if (!fluxfile) {
    std::cout << "File "<<fluxfilename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    while(! fluxfile.eof() )
      {
        fluxfile >> enu >> nue>>numu>>nutau>>nuebar>>numubar>>nutaubar;
        if (! fluxfile.eof()) {
	  _nuefluxmap[enu*1000.] = nue;
	  _numufluxmap[enu*1000.] = numu;
	  _nutaufluxmap[enu*1000.] = nutau;
	  _nuebarfluxmap[enu*1000.] = nuebar;
	  _numubarfluxmap[enu*1000.] = numubar;
	  _nutaubarfluxmap[enu*1000.] = nutaubar;

	}
      }
  }

  fluxfile.close();
}


double NumericalFlux::fluxval(double enu,int flavor, double ebinsize)
{
  double flux = 1;

  // Output is flux per ebinsize, for ebinsize in MeV.  enu is in MeV.

    //http://www.bnikolic.co.uk/blog/cpp-map-interp.html
    // Interpolate from the map.  Must have been initalized for output to make sense

  typedef std::map<double, double>::const_iterator i_t;

  // Q is scaled by Rfac (see email from Chuck, Oct 17, 2017)

  std::map<double, double> _fluxmap;

  switch (flavor) {
  case 1: _fluxmap = _nuefluxmap;
    break;
  case 2: _fluxmap = _numufluxmap;
    break;
  case 3: _fluxmap = _nutaufluxmap;
    break;
  case -1: _fluxmap = _nuebarfluxmap;
    break;
  case -2: _fluxmap = _numubarfluxmap;
    break;
  case -3: _fluxmap = _nutaubarfluxmap;
    break;
    std::cout<< "Wrong flavor "<<std::endl;
  }



  i_t i=_fluxmap.upper_bound(enu);
  if(i==_fluxmap.end())
    {
      return (--i)->second;
    }
  if (i==_fluxmap.begin())
    {
      return i->second;
    }
  i_t l=i; --l;

  const double delta=(enu- l->first)/(i->first - l->first);
  flux= delta*i->second +(1-delta)*l->second;

  if (isnan(flux)) {flux=0.;}

  return flux*ebinsize*norm;

}

void NumericalFlux::SetFluxFilename(const char * fname) {
  strcpy(filename, fname);
}

const char * NumericalFlux::GetFluxFilename() {
  return filename;
}

double NumericalFlux::maxEnu()
{

  // Return the maximum energy in MeV.  Map is such that
  //  this should be the same for all flavors

  typedef std::map<double, double>::const_reverse_iterator i_t;
  i_t it = _nuefluxmap.rbegin();
  double maxEnu = it->first;

  return maxEnu;

}


////

double PinchedThermal::fluxval(double Enu, int flavor, double ebinsize)
{

  //  flavor index goes nue, numu, nutau, nuebar, numubar, nutaubar
  // -3, -2, -1, 1, 2, 3
  // Energy should be in MeV

  int j;

  switch (flavor) {
  case 1: j=0;
    break;
  case 2: j=1;
    break;
  case 3: j=2;
    break;
  case -1: j=3;
    break;
  case -2: j=4;
    break;
  case -3: j=5;
    break;
  default: std::cout<< "Incorrect flavor "<<std::endl;
    exit(-1);
  }


  // Conversions for flux at 10 kpc

  double fluxtot= 0.;
  if (alpha[j]>0 && avgen[j]>0  && luminosity[j]>0) {
    const double dist=3.08568025e22; // [dist]=cm, 10 kpc

    double N=pow((alpha[j]+1.),(alpha[j]+1.))/(avgen[j]*tgamma(alpha[j]+1.));
    double phi=N*pow((Enu/avgen[j]),alpha[j])*exp((-1.)*(alpha[j]+1.)*Enu/avgen[j]);

    fluxtot = 1./(4*M_PI*dist*dist)*luminosity[j]/avgen[j]*phi*ebinsize;

  }
  return fluxtot*norm;

}


double PinchedThermal::maxEnu()
{

  // Return the maximum energy in MeV

  double maxEnu = 50.;

  return maxEnu;

}

void PinchedThermal::SetLuminosity(double* lumi) {

  // Input luminosity is erg/s, store as MeV/s
  const double mevpererg = 624150.;

  int i;
  for (i=0;i<6;i++) {
    luminosity[i] = lumi[i]*mevpererg;
  }

}

double* PinchedThermal::GetLuminosity() {
  return luminosity;
}

void PinchedThermal::SetAvgEn(double* en) {

  // In MeV
  int i;
  for (i=0;i<6;i++) {
    avgen[i] = en[i];
  }

}

double* PinchedThermal::GetAvgEn() {
  return avgen;
}

void PinchedThermal::SetAlpha(double* a) {

  int i;
  for (i=0;i<6;i++) {
    alpha[i] = a[i];
  }

}

double* PinchedThermal::GetAlpha() {
  return alpha;
}


////////

NuFlux::NuFlux(){
  norm = 1.;
}

NuFlux::NuFlux(const char * type)
{
   strcpy(fluxtype,type);
   norm = 1.;
}


void NuFlux::Setfluxtype(const char * type) {
  strcpy(fluxtype, type);
}

const char * NuFlux::Getfluxtype() {
  return fluxtype;
}

void NuFlux::SetNorm(double normval) {
  norm = normval;
}

double NuFlux::GetNorm() {
  return norm;
}


void NuFlux::SetOscParam(double* ua4, double m, double b) {
  doosc = 1;
  // Unitarity constraint... user must ensure not violated
  double us4_2 = 1.-pow(ua4[0],2)-pow(ua4[1],2)-pow(ua4[2],2);
  sin22thes = 4.*pow(ua4[0],2)*us4_2;
  sin22thmus = 4.*pow(ua4[1],2)*us4_2;
  sin22thtaus = 4.*pow(ua4[2],2)*us4_2;
  dm2 = m; // ev^2
  baseline = b; // cm
}

void NuFlux::GetOscParam(double* oscparam) {

  oscparam[0]=sin22thes;
  oscparam[1]=sin22thmus;
  oscparam[2]=sin22thtaus;
  oscparam[3]= dm2;

}

////
