#ifndef _DetectorResponse_
#define _DetectorResponse_

#include <map>
#include <vector>
#include <fstream>
#include <math.h>
#include <string>

// Energies in MeV
class DetectorResponse
{

 protected:

  char detectortype[80];

  // For QF in numerical format
  std::map<double,double> _qfmap;
  char qffilename[80];

  // For QF in polynomial format
  char qfpolyfilename[80];
  std::vector<double> qfpolycoeff;
  double qfpolyrange[2]; // Range of validity for polynomial


  // For Gaussian smearing in polynomial format
  char gspolyfilename[80];
  std::vector<double> gspolycoeff;
  double gspolyrange[2]; // Range of validity for polynomial


  // For a step-function threshold in MeVr

  double step_thresh=0.;

  // For efficiency in numerical format
  std::map<double,double> _efficmap;
  char efficfilename[80];

 public: 
  DetectorResponse();
  DetectorResponse(const char *);
  ~DetectorResponse(){};

  void Setdetectortype(const char *);
  const char * Getdetectortype();

  // For QF in numerical format
  void SetQFFilename(const char * qffilename);
  const char * GetQFFilename();
  void ReadQFFile();
  double qfnum(double);
  double maxErec();

  // For QF in polynomial format
  void SetQFPolyRange(double*);
  double* GetQFPolyRange();

  void SetQFPolyFilename(const char * qfpolyfilename);
  const char * GetQFPolyFilename();
  void ReadQFPolyFile();
  double qfpoly(double);
  double qfpolyderiv(double);


  // For Gaussian smearing in polynomial formats

  int gstype;
  void SetGSType(int);
  int GetGSType();

  void SetGSPolyRange(double*);
  double* GetGSPolyRange();

  void SetGSPolyFilename(const char * gspolyfilename);
  const char * GetGSPolyFilename();
  void ReadGSPolyFile();
  double gspoly(double);
  std::map<double,double> Smear(std::map<double,double>);

  int NEeeBin;
  void SetNEeeBin(int);
  int GetNEeeBin();


  double maxEee;
  void SetMaxEee(double);
  double GetMaxEee();
  double** SmearingMatrix;
  // Not bothering to clean this up with a delete method, I'm a bad person
  void SetGaussSmearingMatrix();

  // For efficiency as a function of Erec, file in numerical format

  void SetEfficFilename(const char * efficfilename);
  const char * GetEfficFilename();
  void ReadEfficFile();
  double efficnum(double);
  double maxEfficErec();

  void SetStepThresh(double);
  double GetStepThresh();

};



#endif
