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
        std::string name;

	    double drate_e=0;
	    double drate_ebar=0;
	    double drate_mu=0;
	    double drate_mubar=0;
	    double drate_tau=0;
	    double drate_taubar=0;
        double knuMin=0;
        double knuMax=0;
        double knuStep=0;


        // Contributions for each component
        double diffrate_e[15] = {0.};
        double diffrate_ebar[15] = {0.};
        double diffrate_mu[15] = {0.};
        double diffrate_mubar[15] = {0.};
        double diffrate_tau[15] = {0.};
        double diffrate_taubar[15] = {0.};

        double sum_diffrate_e=0;
        double sum_diffrate_ebar=0;
        double sum_diffrate_mu=0;
        double sum_diffrate_mubar=0;
        double sum_diffrate_tau=0;
        double sum_diffrate_taubar=0;

        int Z=0;

        Xscn();

        void Set_0() {
            drate_e=0;
            drate_ebar=0;
            drate_mu=0;
            drate_mubar=0;
            drate_tau=0;
            drate_taubar=0;}

        void Set_Z(double z) {Z=z;}

        void Set_diffrate_0() {
            for (int i=0; i<15; i++) {
                diffrate_e[i]=0;
                diffrate_ebar[i]=0;
                diffrate_mu[i]=0;
                diffrate_mubar[i]=0;
                diffrate_tau[i]=0;
                diffrate_taubar[i]=0;
            }

            sum_diffrate_e=0;
            sum_diffrate_ebar=0;
            sum_diffrate_mu=0;
            sum_diffrate_mubar=0;
            sum_diffrate_tau=0;
            sum_diffrate_taubar=0;
        }

        void Set_knu(double min, double max, double step);


        virtual void calc_drate(double mass, double erec, NuFlux *flux)=0;//needs to be overloaded

        virtual void calc_diffrate(int is, double mufact, double ntfac, double mass_fraction[15], NuFlux *flux, Coupling *couplings)=0;

};


class VectorXscn : public Xscn {
    public:
        VectorXscn();
        void calc_drate(double mass, double erec,NuFlux *flux);
        void calc_diffrate(int is, double mufact, double ntfac, double mass_fraction[15], NuFlux *flux, Coupling *couplings);
};

class AxialXscn : public Xscn {
    public:
        AxialXscn();
        void calc_drate(double mass, double erec,NuFlux *flux);
        void calc_diffrate(int is, double mufact, double ntfac, double mass_fraction[15], NuFlux *flux, Coupling *couplings);
};

class InterfXscn : public Xscn {
    public:
        InterfXscn();
        void calc_drate(double mass, double erec,NuFlux *flux);
        void calc_diffrate(int is, double mufact, double ntfac, double mass_fraction[15], NuFlux *flux, Coupling *couplings);
};

class MagXscn : public Xscn {
    public:
        MagXscn();
        void calc_drate(double mass, double erec,NuFlux *flux);
        void calc_diffrate(int is, double mufact, double ntfac, double mass_fraction[15], NuFlux *flux, Coupling *couplings);
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
