#ifndef _xscns_
#define _xscns_

#include <iostream>
#include <nlohmann/json.hpp>
#include <queue>
#include "../fluxes/NuFlux.h"
class Coupling {
    public:
	    double GV_wff_e;
	    double GV_wff_ebar;
	    double GV_wff_mu;
	    double GV_wff_mubar;
	    double GV_wff_tau;
	    double GV_wff_taubar;
	    double GA_wff;
	    double GA_bar_wff;

	    double munu_e=0.;
	    double munu_ebar=0.;
	    double munu_mu=0.;
	    double munu_mubar=0.;
	    double munu_tau=0.;
	    double munu_taubar=0.;

        Coupling(nlohmann::json, double ff[4],int ,int ,int ,int, double);
        Coupling();
        void Set_Couplings(nlohmann::json j, double, double, double, double, int Z, int Nn, int Ndiff, int Zdiff, double Q);

};

class Xscn{

    public:
	    double drate_e=0;
	    double drate_ebar=0;
	    double drate_mu=0;
	    double drate_mubar=0;
	    double drate_tau=0;
	    double drate_taubar=0;
        double knuMin=0;
        double knuMax=0;
        double knuStep=0;

        Xscn();

        void Set_0() {
            drate_e=0;
            drate_ebar=0;
            drate_mu=0;
            drate_mubar=0;
            drate_tau=0;
            drate_taubar=0;}

        void Set_knu(double min, double max, double step);


        virtual void calc_drate(double mass, double erec, NuFlux *flux)=0;//needs to be overloaded

};


class VectorXscn : public Xscn {
    public:
        VectorXscn();
        void calc_drate(double mass, double erec,NuFlux *flux);
};

class AxialXscn : public Xscn {
    public:
        AxialXscn();
        void calc_drate(double mass, double erec,NuFlux *flux);
};

class InterfXscn : public Xscn {
    public:
        InterfXscn();
        void calc_drate(double mass, double erec,NuFlux *flux);
};

class MagXscn : public Xscn {
    public:
        MagXscn();
        void calc_drate(double mass, double erec,NuFlux *flux);
};



double diffxscnvec(double, double, double);
double diffxscnaxial(double, double, double);
double diffxscninterf(double, double, double);
double diffxscnmag(double, double);

void sm_vector_couplings(int, double*);
void sm_axial_couplings(int, int, double*);

double GV_SM(int,int, int);
double GA_SM(int,int, int, int, int, int);


double GV_nsi_nonuniv(int, int, int, double, double, double, double, double, double, double, double, double, double);

double GV_nsi_fc2(int, int, int, double, double, double, double, double, double, double, double, double, double);

double chgradcorr(int,int);
double chgradcorr_tomalak(double,int);

        void nsi_vector_couplings(double* , double , double );
double diffnuelectronxscn(int,double, double);
double diffnuelectronxscn2(int,double, double, double, double, int);

double diffangdist(double, double);

double mufactor(double);
double taufactor(double);

std::queue<Xscn*> XscnQueue(nlohmann::json);
#endif
