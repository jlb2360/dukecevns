#include "FormFactor.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <math.h>


void Horowitz::ReadFFfile()
{

  double q,ff;
  std::ifstream fffile;
  std::string fffilename = filename;
  fffile.open(fffilename.c_str());
  if (!fffile) {
    std::cout << "File "<<fffilename<<" does not exist!" <<std::endl;
    exit(-1);
  } else {
    while(! fffile.eof() )
      {
        fffile >> q >> ff;
	//	std::cout << q << " "<<ff << std::endl;
        _ffmap[q] = ff;
      }
  }

  fffile.close();
}

double Horowitz::FFval(double Q)
{
  double ff = 1;

    //http://www.bnikolic.co.uk/blog/cpp-map-interp.html
    // Interpolate from the map.  Must have been initalize for output to make sense

  typedef std::map<double, double>::const_iterator i_t;

  // Q is scaled by Rfac (see email from Chuck, Oct 17, 2017)


  // Scale the radius
  Q *=Rfac;

  i_t i=_ffmap.upper_bound(Q);
  if(i==_ffmap.end())
    {
      return (--i)->second;
    }
  if (i==_ffmap.begin())
    {
      return i->second;
    }
  i_t l=i; --l;

  const double delta=(Q- l->first)/(i->first - l->first);
  ff= delta*i->second +(1-delta)*l->second;

  if (isnan(ff)) {ff=1.;}

  // Note not squared in the file
  return ff;

}


void Horowitz::SetFFfilename(const char * fname) {
  strcpy(filename, fname);
}

const char * Horowitz::GetFFfilename() {
  return filename;
}


////////

FormFactor::FormFactor()
{
  Rfac=1.;
}

FormFactor::FormFactor(const char * type)
{
  strcpy(fftype,type);
  Rfac=1.;

}


void FormFactor::SetA(int Aval) {

  A = Aval;

}

int FormFactor::GetA() { return A;}


void FormFactor::SetZ(int Zval) {

  Z = Zval;

}

int FormFactor::GetZ() { return Z;}



void FormFactor::SetRfac(double Rfacval) {

  Rfac = Rfacval;

}

double FormFactor::GetRfac() { return Rfac;}

void FormFactor::Setfftype(const char * type) {
  strcpy(fftype, type);
}

const char * FormFactor::Getfftype() {
  return fftype;
}
