#include "NuFlux.h"
#include <cstring>
#include <nlohmann/json.hpp>
#include <cstdlib>
#include <math.h>
#include "NuFlux.h"
#include <iostream>

void get_flavor_weight(int, double,double,double*,double*, double*, double*);

PiDAR::PiDAR(nlohmann::json j){
    const char* ft = "PiDAR";
    std::strcpy(fluxtype, ft);
    norm = 1.;

    wnumu=1.;
    wnumubar=1.;
    wnue=1.;

    convolved=0;
    perpot = 0;
    usernorm = 1.;

    // Can't seem to do this with a nested key easily; make convolved not nested
    if (j.find("convolved") != j.end()) {
          convolved = j["convolved"];
    }
    std::cout << "Convolved "<<convolved<<std::endl;


    if (j.find("perpot") != j.end()) {
          perpot = j["perpot"];
    }

    // A generic normalization factor for the recoil spectrum

    if (j.find("usernorm") != j.end()) {
        usernorm = j["usernorm"];
        std::cout << "User normalization: "<<usernorm<<std::endl;
    }


    // Don't use flavor weights if using snsflux numerical flux; it that should take care of the weighting
    //get_flavor_weight(1400.,7400.,&wnumu,&wnumubar,&wnue);

    tw1 = j["timewindow"]["start"];
    tw2 = j["timewindow"]["end"];


    if (convolved>0) {
        double teffic_params[3]={0,0,0};

        if (j.find("teffic") != j.end()) {
            teffic_params[0] = j["teffic"]["offset"];
            teffic_params[1] = j["teffic"]["a"];
            teffic_params[2] = j["teffic"]["b"];
        }

        get_flavor_weight(convolved,tw1,tw2,teffic_params,&wnumu,&wnumubar,&wnue);
    }

    // Flavor weighting from file if chosen. Do individually so as not to override convolved

    if (j.find("wnumu") != j.end()) {
        wnumu = j["wnumu"];
    }
    if (j.find("wnumubar") != j.end()) {
        wnumubar = j["wnumubar"];
    }
    if (j.find("wnue") != j.end()) {
        wnue = j["wnue"];
    }


    std::cout << "Flavor weights: "<< wnumu<<" "<<wnumubar<<" "<<wnue<<std::endl;

    double mevperproton = j["flux"]["mevperproton"];
    double jperproton= mevperproton*1.e6*1.6021e-19;
    double beampower = j["flux"]["power"];
    beampower*=1.e6; // in Joules/s
    protonspersec = beampower/jperproton;
    double nusperprotonperflavor = j["flux"]["nusperprotonperflavor"];
    nuspersecperflavor = nusperprotonperflavor*protonspersec;
    double dist = j["distance"];
    std::cout << "Nus per sec per flavor "<<nuspersecperflavor<<" "<<dist<<" flux "<<nuspersecperflavor/(4*M_PI*dist*dist)<<std::endl;

    SetNorm(nuspersecperflavor/(4*M_PI*dist*dist));

    if (j.find("doosc") != j.end()) {
        int doosc = j["doosc"];
        if (doosc == 1) {
            double ua4[3];
            ua4[0] = j["osc"]["ue4"];
            ua4[1] = j["osc"]["umu4"];
            ua4[2] = j["osc"]["utau4"];
            double dm2 = j["osc"]["dm2"];
            SetOscParam(ua4,dm2,dist);
        }
    }

    double hours =j["flux"]["hours"];
    exposure = 3600.*hours;
}


double PiDAR::fluxval(double Enu, int flavor, double ebinsize)
{


  // Energies in MeV
  // 1 = e, 2 = mu
  //const double mmu = 105.6;
  const double mmu = 105.66837;
  //  const double Enumu = 29.9;
  const double Enumu = 29.792;

  const double a= 2/mmu;
  double flux = 0.;

  if (Enu>mmu/2.) {
    flux = 0.;
    return flux;
  }

  if (flavor == 1) {
    flux = norm*12*pow(a*Enu,2)*(1-a*Enu)*a*ebinsize;
   } else if (flavor == 2) {

    if (fabs(Enu-Enumu)<ebinsize/2.) {
      flux = norm;
    }

  } else if (flavor == -2) {
     flux = norm*2*pow(a*Enu,2)*(3-2*a*Enu)*a*ebinsize;
  } else {

    flux = 0.;
  }

  // Oscillate if requested

  //  std::cout << "1: "<< flux << std::endl;
  double sin22th;
  if (doosc==1) {

    if (abs(flavor) == 1) {
      sin22th = sin22thes;
    } else if (abs(flavor) == 2) {
      sin22th = sin22thmus;
    }  else if (abs(flavor) == 3) {
      sin22th = sin22thtaus;
    } else {
      sin22th = 0.;
    }

    // Simple sterile disappearance,  baseline in cm (converted to m in expression), Enu in MeV
    flux *= (1-sin22th*pow(sin(1.27*dm2*(baseline/100.)/Enu),2));


  }

  //  std::cout << 1.-sin22th<<" "<<dm2<<" "<<Enu<<" "<<" "<<baseline<<" "<<pow(sin(1.27*dm2*baseline/Enu),2)<<std::endl;
  //std::cout << flux << std::endl;

  return flux;

}

double PiDAR::maxEnu()
{

  // Return the maximum energy in MeV

  // To compare with old
  //  double maxEnu = 105.6/2.;
  double maxEnu = 105.66837/2.;
  return maxEnu;

}

////
