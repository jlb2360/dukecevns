#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

#include "DetectorResponse/DetectorResponse.h"
#include "FormFactors/FormFactor.h"
#include "fluxes/NuFlux.h"


#include <algorithm>
#include <cstdlib>
#include <string>
#include <math.h>
#include <map>

#include "CrossSections/xscns.h"


void get_flavor_weight(int, double,double,double*,double*, double*, double*);


int main(int argc, char * argv[] )
{
    //These need to be in main so that it will compile with types.
    #include "DetectorResponse/isomaps.h"
    #include "DetectorResponse/mixtures.h"

  if (argc<2) {
    std::cout << "Usage:  ./sns_rates [jsonfile]"<<std::endl;
    exit(0);
  }

  const char * jsonfile = argv[1];

  std::string jsonfilename = "jsonfiles/"+std::string(jsonfile)+".json";

  // creating variables
  double M;
  double Delta;
  int Nn,Z,A;
  int Zdiff, Ndiff;

  double Ntargets = 0; // Total targets
  int is=0;


// Read a JSON file with the parameters
    std::ifstream i(jsonfilename);
    json j;
    i >> j;

    // print values
    std::cout << j << '\n';

    std::cout << j["flux"]["nusperprotonperflavor"]<<std::endl;


  // Set up the form factor

  std::string ffname = j["formfactor"]["type"];

  int noff = 0.;
  if (ffname == "unity") {
    noff = 1.;
  }


  // Array of pointers to form factors, for protons and neutrons separately, axial and vector separately
  // (although small differences)

  // max components is 15, found in mixtures.h
  FormFactor** ffpv;
  ffpv = new FormFactor*[max_components];
  FormFactor** ffpa;
  ffpa = new FormFactor*[max_components];

  FormFactor** ffnv;
  ffnv = new FormFactor*[max_components];
  FormFactor** ffna;
  ffna = new FormFactor*[max_components];


  // Set up the neutrino flux
    PiDAR* snsflux = new PiDAR(j);

  // Set up the detector response-- this is an overall detector response
  DetectorResponse* detresp = new DetectorResponse(j);

  Coupling* couplings = new Coupling();

  std::queue<Xscn*> xscnQ = XscnQueue(j);

  double kmax = snsflux->maxEnu();
  std::cout << "kmax "<<kmax << std::endl;


  std::ofstream outfile;
  std::string outfilename;
  outfilename = "out/sns_diff_rates-"+std::string(jsonfile)+"-"+detresp->material+"-"+ffname+".out";
  outfile.open(outfilename);
  std::cout << outfilename <<std::endl;

  std::vector<std::string>::iterator v = detresp -> isotope_component.begin();
  std::string isotope;

  // First get the total mass.  Get also the maximum recoil values

  double Mtot = 0;

  DetectorResponse** qffunc;
  qffunc = new DetectorResponse*[max_components];

  // First loop over materials, to set up material-specific arrays
  double minM = 1.e10;
  while( v != detresp->isotope_component.end()) {

    isotope = *v;
    std::cout << "isotope"<< isotope << std::endl;
    std::string isoname = std::string(isotope);

    Z = Zs[std::string(isotope)];
    Nn = Ns[std::string(isotope)];
    Delta = Deltas[std::string(isotope)];
    M = (Z+Nn)*amu - Z*me + Delta;
    if (M<minM) {minM=M;}
    Mtot += M*detresp->fraction[is];


    // Set up the form factor for this isotope

    double nvrfact = j["formfactor"]["nvrfact"];
    double narfact = j["formfactor"]["narfact"];
    double pvrfact = j["formfactor"]["pvrfact"];
    double parfact = j["formfactor"]["parfact"];


    if (ffname == "helm") {

      ffnv[is] = new Helm(j["formfactor"]["nvsfact"],j["formfactor"]["nvrfact"]);

      ffna[is] = new Helm(j["formfactor"]["nasfact"],j["formfactor"]["narfact"]);

      ffpv[is] = new Helm(j["formfactor"]["pvsfact"],j["formfactor"]["pvrfact"]);

      ffpa[is] = new Helm(j["formfactor"]["pasfact"],j["formfactor"]["parfact"]);

    }
    else if (ffname == "klein") {

      ffnv[is] = new Klein(j["formfactor"]["nvak"],j["formfactor"]["nvskin"],j["formfactor"]["nvrfact"]);

      ffna[is] = new Klein(j["formfactor"]["naak"],j["formfactor"]["naskin"],j["formfactor"]["narfact"]);

      ffpv[is] = new Klein(j["formfactor"]["pvak"],j["formfactor"]["pvskin"],j["formfactor"]["pvrfact"]);

      ffpa[is] = new Klein(j["formfactor"]["paak"],j["formfactor"]["paskin"],j["formfactor"]["parfact"]);


    }
    else if  (ffname =="horowitz"){
      Horowitz* horowitzffnv = new Horowitz();
      ffnv[is] = horowitzffnv;

      Horowitz* horowitzffna = new Horowitz();
      ffna[is] = horowitzffna;

      Horowitz* horowitzffpv = new Horowitz();
      ffpv[is] = horowitzffpv;

      Horowitz* horowitzffpa = new Horowitz();
      ffpa[is] = horowitzffpa;


      std::transform(isoname.begin(), isoname.end(),isoname.begin(), ::toupper);
      std::string horowitz_filename = isoname+".FF";
      std::cout << horowitz_filename << std::endl;
      horowitzffnv->SetFFfilename(horowitz_filename.c_str());
      horowitzffnv->ReadFFfile();
      horowitzffnv->SetRfac(nvrfact);

      horowitzffna->SetFFfilename(horowitz_filename.c_str());
      horowitzffna->ReadFFfile();
      horowitzffna->SetRfac(narfact);

      //      std::cout << horowitz_filename << std::endl;
      // Not really appropriate for protons, but using the structure
      horowitzffpv->SetFFfilename(horowitz_filename.c_str());
      horowitzffpv->ReadFFfile();
      horowitzffpv->SetRfac(pvrfact);


      horowitzffpa->SetFFfilename(horowitz_filename.c_str());
      horowitzffpa->ReadFFfile();
      horowitzffpa->SetRfac(parfact);


    }

    A = Nn + Z;
    if (noff == 0) {
      ffnv[is]->SetA(A);
      ffna[is]->SetA(A);
      ffpv[is]->SetA(A);
      ffpa[is]->SetA(A);

      ffnv[is]->SetZ(Z);
      ffna[is]->SetZ(Z);
      ffpv[is]->SetZ(Z);
      ffpa[is]->SetZ(Z);

    }

  // Set up detector quenching factors for each component

    if (detresp->qftype == "poly") {
      std::string qffilename;
      qffilename = "qf/"+std::string(detresp->qfname)+"_"+isoname+"_qf.txt";
      DetectorResponse* qf = new DetectorResponse();

      std::cout << "Quenching factor: "<<qffilename<<std::endl;
      qffunc[is] = qf;
      qffunc[is]->SetQFPolyFilename(qffilename.c_str());
      qffunc[is]->ReadQFPolyFile();
    } else if (detresp->qftype == "numerical") {
      std::string qffilename;
      qffilename = "qf/"+std::string(detresp->qfname)+"_"+isoname+"_qfnum.txt";
      DetectorResponse* qf = new DetectorResponse();

      std::cout << "Quenching factor: "<<qffilename<<std::endl;
      qffunc[is] = qf;
      qffunc[is]->SetQFFilename(qffilename.c_str());
      qffunc[is]->ReadQFFile();
    }
    v++; is++;
  }

    // Overall norm factor


  double norm_factor = detresp->detector_mass*snsflux->exposure*snsflux->usernorm;

  if (snsflux->perpot == 1) {
    double pot = snsflux->protonspersec*snsflux->exposure;
    norm_factor/= pot;
    std::cout << "Protons on target: "<<pot<<" ; output per pot "<<std::endl;
  }


  // Set up arrays for quenched total rates... could make this a stl vec
   // but this is probably more efficient



  // Use the mass of the lightest component
  double erecmaxall = 2*kmax*kmax/(minM+2*kmax);

  //double recoilthresh = 0.013; //MeVr
  double erecstart = detresp->recoilthresh;
  double erecend = detresp->recoilupperthresh > detresp->recoilthresh ?
                   std::min(erecmaxall, detresp->recoilupperthresh) :
                   erecmaxall;

  double erecstep = 0.0001;
  if (j.find("erecstep") != j.end()) {
    erecstep = j["erecstep"];
  }

  // Set up the recoil energy arrays

  int maxiq;
  maxiq = int((erecend-erecstart)/erecstep)+1;

  double* Er = new double[maxiq];
  double** Eee = new double*[max_components];
  double** dNdEee = new double*[max_components];
  double** dNdEr = new double*[max_components];
  for (int iee = 0;iee < max_components;iee++) {
    Eee[iee] = new double[maxiq];
    dNdEee[iee] = new double[maxiq];
    dNdEr[iee] = new double[maxiq];
  }

  double* dNdErall = new double[maxiq];

  // Now compute the differential recoil spectra

    double Erec;
    double knu;

    std::cout << "erecmaxall "<<erecmaxall<<std::endl;

    double knustep = 0.0001;
    if (j.find("knustep") != j.end()) {
      knustep = j["knustep"];
    }

   // The totals
   double toterecoil = 0.;
   double totevents = 0.;
   double toteventsnue = 0.;
   double toteventsnuebar = 0.;
   double toteventsnumu = 0.;
   double toteventsnumubar = 0.;
   double toteventsnutau= 0.;
   double toteventsnutaubar = 0.;


   int iq=0;
   // Loop over recoil energy
   for (Erec=erecstart+erecstep;Erec<=erecend; Erec+=erecstep) {

     Er[iq] = Erec;

     // Contributions for each component
     double diffrate_e_vec[max_components]={0.};
     double diffrate_ebar_vec[max_components]={0.};
     double diffrate_mu_vec[max_components]={0.};
     double diffrate_mubar_vec[max_components]={0.};
     double diffrate_tau_vec[max_components]={0.};
     double diffrate_taubar_vec[max_components]={0.};

     double diffrate_e_axial[max_components]={0.};
     double diffrate_ebar_axial[max_components]={0.};
     double diffrate_mu_axial[max_components]={0.};
     double diffrate_mubar_axial[max_components]={0.};
     double diffrate_tau_axial[max_components]={0.};
     double diffrate_taubar_axial[max_components]={0.};

     double diffrate_e_interf[max_components]={0.};
     double diffrate_ebar_interf[max_components]={0.};
     double diffrate_mu_interf[max_components]={0.};
     double diffrate_mubar_interf[max_components]={0.};
     double diffrate_tau_interf[max_components]={0.};
     double diffrate_taubar_interf[max_components]={0.};

     double diffrate_e_mag[max_components]={0.};
     double diffrate_ebar_mag[max_components]={0.};
     double diffrate_mu_mag[max_components]={0.};
     double diffrate_mubar_mag[max_components]={0.};
     double diffrate_tau_mag[max_components]={0.};
     double diffrate_taubar_mag[max_components]={0.};



     // Sum for each component,  not quenched

     double sum_diffrate_e_vec=0;
     double sum_diffrate_ebar_vec=0;
     double sum_diffrate_mu_vec=0;
     double sum_diffrate_mubar_vec=0;
     double sum_diffrate_tau_vec=0;
     double sum_diffrate_taubar_vec=0;

     double sum_diffrate_e_axial=0;
     double sum_diffrate_ebar_axial=0;
     double sum_diffrate_mu_axial=0;
     double sum_diffrate_mubar_axial=0;
     double sum_diffrate_tau_axial=0;
     double sum_diffrate_taubar_axial=0;

     double sum_diffrate_e_interf=0;
     double sum_diffrate_ebar_interf=0;
     double sum_diffrate_mu_interf=0;
     double sum_diffrate_mubar_interf=0;
     double sum_diffrate_tau_interf=0;
     double sum_diffrate_taubar_interf=0;

     double sum_diffrate_e_mag=0;
     double sum_diffrate_ebar_mag=0;
     double sum_diffrate_mu_mag=0;
     double sum_diffrate_mubar_mag=0;
     double sum_diffrate_tau_mag=0;
     double sum_diffrate_taubar_mag=0;



     // With efficiency, which is a function of Erec in MeV in this formuation

     double recoil_eff_factor = 1;
     if (detresp->effname != "none" && detresp->eff_type == "erecoil") {
       recoil_eff_factor = detresp->efficnum(Erec);
     }

     //     double eff_factor = 1.;
     //     std::cout << "recoil eff factor "<<recoil_eff_factor<<std::endl;
     // Skip if too small contribution
     if (recoil_eff_factor>0) {


       v = detresp->isotope_component.begin();
       // Now loop over components
       is=0;
       while( v != detresp->isotope_component.end()) {

	        isotope = *v;
	        //	  std::cout << "isotope"<< isotope << std::endl;

	        Z = Zs[std::string(isotope)];
	        Nn = Ns[std::string(isotope)];
	        Delta = Deltas[std::string(isotope)];
	        M = (Z+Nn)*amu - Z*me + Delta;

	        Zdiff = Zdiffs[std::string(isotope)];
	        Ndiff = Ndiffs[std::string(isotope)];

	        mass_fraction[is] = M/Mtot*detresp->fraction[is];

	        A = Nn + Z;
	        //	 std::cout << " Z "<<Z<<" N "<<Nn<<" A "<<A<<" M "<<M << " "<<mass_fraction[is]<<std::endl;

	        // Loop over neutrino energy contributions

	        // Minimum neutrino energy contributing to a given recoil energy

	        double knumin = 0.5*(Erec+sqrt(Erec*Erec+2*M*Erec));
	        double hbarc = 197.327; // MeV-fm, convert for Q in MeV for ff
	        double Q = sqrt(2*M*Erec+Erec*Erec); // MeV
	        //	 double Q = sqrt(2*M*Erec); // MeV

	        double qq = Q/hbarc;
	        //    double ff2 = helmff->FFval(qq);
	        double ffnvval;
	        double ffnaval;
	        double ffpvval;
	        double ffpaval;
	        if (noff == 0) {
	            ffnvval = ffnv[is]->FFval(qq);
	            ffnaval = ffna[is]->FFval(qq);
	            ffpvval = ffpv[is]->FFval(qq);
	            ffpaval = ffpa[is]->FFval(qq);
	        } else {
	            ffnvval = 1.;
	            ffnaval = 1.;
	            ffpvval = 1.;
	            ffpaval = 1.;
	        }


	 // SM Couplings
     couplings->Set_Couplings(j,ffpvval,ffnvval,ffpaval,ffnaval,Z,Nn,Ndiff,Zdiff,Q);



	 // Targets for one ton of material
	 // Will be weighted by mass fraction

	 double Nt = 1.e6/(M/amu)*6.0221409e23;


	 //	 std::cout << "Number of targets "<<Nt<<std::endl;
	 // A2: G^2/(2Pi) * hbarcinmeters^-4
	 double ntfac = Nt;

	  // Quenching factor for this component and Eee for this Erec

	 double qfderiv=1;
	 if (detresp->qftype == "poly") {
	  Eee[is][iq] = qffunc[is]->qfpoly(Erec)*Erec;
	  qfderiv = abs(qffunc[is]->qfpolyderiv(Erec));
	 } else if (detresp->qftype == "numerical") {
	  Eee[is][iq] = qffunc[is]->qfnum(Erec)*Erec;
	  qfderiv = abs(qffunc[is]->qfnumderiv(Erec));
	 }


     queue<Xscn*> tempQ;

     Xscn* xscn;

     while (!xscnQ.empty()) {
         xscn = xscnQ.front();
         xscnQ.pop();
         tempQ.push(xscn);

         xscn->Set_0();
         xscn->Set_knu(knumin, kmax, knustep);

         xscn->calc_drate(M, Erec, snsflux);
     }

	 // Now multiply by target-dependent factors and add up this recoil energy bin

     double mufact = 1;


     xscn = tempQ.front();
     tempQ.pop();
     xscnQ.push(xscn);
	 diffrate_e_vec[is] += ntfac*pow(couplings->GV_wff_e,2)*mass_fraction[is]*xscn->drate_e*snsflux->wnue;
	 diffrate_ebar_vec[is] += ntfac*pow(couplings->GV_wff_ebar,2)*mass_fraction[is]*xscn->drate_ebar;
	 diffrate_mu_vec[is] += ntfac*pow(couplings->GV_wff_mu,2)*mass_fraction[is]*xscn->drate_mu*mufact*snsflux->wnumu;
	 diffrate_mubar_vec[is] += ntfac*pow(couplings->GV_wff_mubar,2)*mass_fraction[is]*xscn->drate_mubar*mufact*snsflux->wnumubar;
	 diffrate_tau_vec[is] +=  ntfac*pow(couplings->GV_wff_tau,2)*mass_fraction[is]*xscn->drate_tau;
	 diffrate_taubar_vec[is] += ntfac*pow(couplings->GV_wff_taubar,2)*mass_fraction[is]*xscn->drate_taubar;

     xscn = tempQ.front();
     tempQ.pop();
     xscnQ.push(xscn);

	 diffrate_e_axial[is] += ntfac*pow(couplings->GA_wff,2)*mass_fraction[is]*xscn->drate_e*snsflux->wnue;
	 diffrate_ebar_axial[is] += ntfac*pow(couplings->GA_bar_wff,2)*mass_fraction[is]*xscn->drate_ebar;
	 diffrate_mu_axial[is] += ntfac*pow(couplings->GA_wff,2)*mass_fraction[is]*xscn->drate_mu*mufact*snsflux->wnumu;
	 diffrate_mubar_axial[is] += ntfac*pow(couplings->GA_bar_wff,2)*mass_fraction[is]*xscn->drate_mubar*mufact*snsflux->wnumubar;
	 diffrate_tau_axial[is] +=  ntfac*pow(couplings->GA_wff,2)*mass_fraction[is]*xscn->drate_tau;
	 diffrate_taubar_axial[is] +=  ntfac*pow(couplings->GA_bar_wff,2)*mass_fraction[is]*xscn->drate_taubar;

     xscn = tempQ.front();
     tempQ.pop();
     xscnQ.push(xscn);

	 diffrate_e_interf[is] += ntfac*couplings->GV_wff_e*couplings->GA_wff*mass_fraction[is]*xscn->drate_e*snsflux->wnue;
	 diffrate_ebar_interf[is] += ntfac*couplings->GV_wff_ebar*couplings->GA_bar_wff*mass_fraction[is]*xscn->drate_ebar;
	 diffrate_mu_interf[is] += ntfac*couplings->GV_wff_mu*couplings->GA_wff*mass_fraction[is]*xscn->drate_mu*mufact*snsflux->wnumu;
	 diffrate_mubar_interf[is] += ntfac*couplings->GV_wff_mubar*couplings->GA_bar_wff*mass_fraction[is]*xscn->drate_mubar*mufact*snsflux->wnumubar;
	 diffrate_tau_interf[is] +=  ntfac*couplings->GV_wff_tau*couplings->GA_wff*mass_fraction[is]*xscn->drate_tau;
	 diffrate_taubar_interf[is] +=  ntfac*couplings->GV_wff_taubar*couplings->GA_bar_wff*mass_fraction[is]*xscn->drate_taubar;

     if (j.find("magmon") != j.end()) {
        xscn = tempQ.front();
        tempQ.pop();
        xscnQ.push(xscn);
	    diffrate_e_mag[is] += ntfac*pow(couplings->munu_e,2)*pow(Z,2)*mass_fraction[is]*xscn->drate_e*snsflux->wnue;
	    diffrate_ebar_mag[is] += ntfac*pow(couplings->munu_ebar,2)*pow(Z,2)*mass_fraction[is]*xscn->drate_ebar;
	    diffrate_mu_mag[is] += ntfac*pow(couplings->munu_mu,2)*pow(Z,2)*mass_fraction[is]*xscn->drate_mu*snsflux->wnumu;
	    diffrate_mubar_mag[is] += ntfac*pow(couplings->munu_mubar,2)*pow(Z,2)*mass_fraction[is]*xscn->drate_mubar*snsflux->wnumubar;
	    diffrate_tau_mag[is] +=  ntfac*pow(couplings->munu_tau,2)*pow(Z,2)*mass_fraction[is]*xscn->drate_tau;
	    diffrate_taubar_mag[is] +=  ntfac*pow(couplings->munu_taubar,2)*pow(Z,2)*mass_fraction[is]*xscn->drate_taubar;
     }

	  // Now add the contribution from this isotope to the sum


	 sum_diffrate_e_vec += diffrate_e_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_ebar_vec += diffrate_ebar_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mu_vec += diffrate_mu_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mubar_vec += diffrate_mubar_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_tau_vec += diffrate_tau_vec[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_taubar_vec += diffrate_taubar_vec[is]*norm_factor*recoil_eff_factor;

	 sum_diffrate_e_axial += diffrate_e_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_ebar_axial += diffrate_ebar_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mu_axial += diffrate_mu_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mubar_axial += diffrate_mubar_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_tau_axial += diffrate_tau_axial[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_taubar_axial += diffrate_taubar_axial[is]*norm_factor*recoil_eff_factor;

	 sum_diffrate_e_interf += diffrate_e_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_ebar_interf += diffrate_ebar_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mu_interf += diffrate_mu_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mubar_interf += diffrate_mubar_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_tau_interf += diffrate_tau_interf[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_taubar_interf += diffrate_taubar_interf[is]*norm_factor*recoil_eff_factor;


	 sum_diffrate_e_mag += diffrate_e_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_ebar_mag += diffrate_ebar_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mu_mag += diffrate_mu_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_mubar_mag += diffrate_mubar_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_tau_mag += diffrate_tau_mag[is]*norm_factor*recoil_eff_factor;
	 sum_diffrate_taubar_mag += diffrate_taubar_mag[is]*norm_factor*recoil_eff_factor;



	  // Sum for this Erec and isotope
	  double sum_events_iso = 0;
	  sum_events_iso = diffrate_e_vec[is] + diffrate_ebar_vec[is] + diffrate_mu_vec[is]+ diffrate_mubar_vec[is]+ diffrate_tau_vec[is] + diffrate_taubar_vec[is];
	  sum_events_iso += diffrate_e_axial[is] + diffrate_ebar_axial[is] + diffrate_mu_axial[is]+ diffrate_mubar_axial[is]+ diffrate_tau_axial[is] + diffrate_taubar_axial[is];

	  sum_events_iso+= diffrate_e_interf[is] + diffrate_ebar_interf[is] + diffrate_mu_interf[is]+ diffrate_mubar_interf[is]+ diffrate_tau_interf[is] + diffrate_taubar_interf[is];


	  sum_events_iso+= diffrate_e_mag[is] + diffrate_ebar_mag[is] + diffrate_mu_mag[is]+ diffrate_mubar_mag[is]+ diffrate_tau_mag[is] + diffrate_taubar_mag[is];



	  sum_events_iso *= norm_factor*recoil_eff_factor;

	  // Now apply the quenching for this Ee and isotope component
	// sum_events_iso is dNderec

	  dNdEr[is][iq] = sum_events_iso;
	  dNdErall[iq] += sum_events_iso;

	  if (qfderiv>0) {
	    dNdEee[is][iq] = sum_events_iso/qfderiv;
	  } else {
	    dNdEee[is][iq] = 0.;
	  }


	 v++;is++;


       } // End of loop over material components

     } // End of efficiency factor check


     // Only want diff values in scientific format
     std::cout.unsetf(ios::fixed | ios::scientific);


     outfile << Erec<<scientific<<" "<<sum_diffrate_e_vec<<" "<<sum_diffrate_ebar_vec<<" "<<sum_diffrate_mu_vec<<" "<<sum_diffrate_mubar_vec<<" "<<sum_diffrate_tau_vec<<" "<<sum_diffrate_taubar_vec<<" "<<sum_diffrate_e_axial<<" "<<sum_diffrate_ebar_axial<<" "<<sum_diffrate_mu_axial<<" "<<sum_diffrate_mubar_axial<<" "<<sum_diffrate_tau_axial<<" "<<sum_diffrate_taubar_axial<<" "<<sum_diffrate_e_interf<<" "<<sum_diffrate_ebar_interf<<" "<<sum_diffrate_mu_interf<<" "<<sum_diffrate_mubar_interf<<" "<<sum_diffrate_tau_interf<<" "<<sum_diffrate_taubar_interf <<" "<<sum_diffrate_e_mag<<" "<<sum_diffrate_ebar_mag<<" "<<sum_diffrate_mu_mag<<" "<<sum_diffrate_mubar_mag<<" "<<sum_diffrate_tau_mag<<" "<<sum_diffrate_taubar_mag<<std::endl;
	// Reset the format
     std::cout.unsetf(ios::fixed | ios::scientific);

	double events=0;
	events = sum_diffrate_e_vec + sum_diffrate_ebar_vec + sum_diffrate_mu_vec+ sum_diffrate_mubar_vec+ sum_diffrate_tau_vec + sum_diffrate_taubar_vec;
	events += sum_diffrate_e_axial + sum_diffrate_ebar_axial + sum_diffrate_mu_axial+ sum_diffrate_mubar_axial+ sum_diffrate_tau_axial + sum_diffrate_taubar_axial;
        events += sum_diffrate_e_interf + sum_diffrate_ebar_interf + sum_diffrate_mu_interf+ sum_diffrate_mubar_interf+ sum_diffrate_tau_interf + sum_diffrate_taubar_interf;

        events += sum_diffrate_e_mag + sum_diffrate_ebar_mag + sum_diffrate_mu_mag+ sum_diffrate_mubar_mag+ sum_diffrate_tau_mag + sum_diffrate_taubar_mag;

	toteventsnue+= erecstep*(sum_diffrate_e_vec+sum_diffrate_e_axial+sum_diffrate_e_interf+sum_diffrate_e_mag);

	toteventsnumu+= erecstep*(sum_diffrate_mu_vec+sum_diffrate_mu_axial+sum_diffrate_mu_interf+sum_diffrate_mu_mag);

	toteventsnumubar+= erecstep*(sum_diffrate_mubar_vec+sum_diffrate_mubar_axial+sum_diffrate_mubar_interf+sum_diffrate_mubar_mag);


	totevents+=events*erecstep;

	toterecoil += events*Erec*erecstep;

	// Increment bin for quenching

	iq++;


  } // End of loop over Erec

   if (detresp->recoilupperthresh > detresp->recoilthresh) {
      std::cout << "Total events over "<< detresp->recoilthresh*1000.<<" keVr and under "<<detresp->recoilupperthresh*1000<<" keVr: "<<totevents<< std::endl;
   }
   else {
      std::cout << "Total events over "<< detresp->recoilthresh*1000.<<" keVr: "<<totevents<< std::endl;
   }
   std::cout << "Total recoil energy deposited:  "<< toterecoil<< std::endl;

   outfile.close();


  std::ofstream integraloutfile;
  outfilename = "out/sns_diff_rates-"+std::string(jsonfile)+"-"+detresp->material+"-"+ffname+"-integral.out";

  std::cout << outfilename << std::endl;
  integraloutfile.open(outfilename);


  integraloutfile << j << '\n';
  if (detresp->recoilupperthresh > detresp->recoilthresh) {
    integraloutfile << "Total events over "<< detresp->recoilthresh*1000.<<" keVr and under "<<detresp->recoilupperthresh*1000<<" keVr: "<<totevents<< std::endl;
  }
  else {
    integraloutfile << "Total events over "<< detresp->recoilthresh*1000.<<" keVr: "<<totevents<< std::endl;
  }
  integraloutfile.close();

  // Output by isotope, integrated over flavor.
  //Can also output quenched stuff here
   // Loop over isotopes


  v = detresp->isotope_component.begin();
  // Now loop over components
  is=0;
  while( v != detresp->isotope_component.end()) {

    isotope = *v;
    //	  std::cout << "isotope"<< isotope << std::endl;
    std::string isoname = std::string(isotope);

    std::ofstream isooutfile;
    outfilename = "out/sns_diff_rates-"+std::string(jsonfile)+"-"+detresp->material+"-"+ffname+"-"+isoname+".out";

    std::cout << outfilename << std::endl;
    isooutfile.open(outfilename);

    int ie;
    for (ie=0;ie<iq;ie++) {

      isooutfile << Er[ie]<< "  "<<dNdEr[is][ie]<<endl;


    }
    isooutfile.close();
    v++;is++;
  }

  // Integrated over flavor and isotope

    std::ofstream allisooutfile;
    outfilename = "out/sns_diff_rates_alliso-"+std::string(jsonfile)+"-"+detresp->material+"-"+ffname+".out";

    std::cout << outfilename << std::endl;
    allisooutfile.open(outfilename);

    int ie;
    for (ie=0;ie<iq;ie++) {

      allisooutfile << Er[ie]<< "  "<<dNdErall[ie]<<endl;

    }
    allisooutfile.close();

    ////////////  QUENCHED  RESPONSE ///////////////

    // Now dump the quenched output, by isotope.  This has non-uniform Eee energy bins.  At the same time fill some maps to be interpolated for a sum, and get the maximum quenched energy for use for that/

    // Don't have this broken down by flavor and interaction... need to do that

    if (detresp->qfname != "none") {

      double maxeee = 0;
      // One of these per component
      std::map<double, double> _quenchedmap[max_components];

      // The total response
      std::map<double, double> _quenchedtot;

      v = detresp->isotope_component.begin();
      // Now loop over components
      is=0;

      while( v != detresp->isotope_component.end()) {

	        isotope = *v;
	        //	  std::cout << "isotope"<< isotope << std::endl;
	        std::string isoname = std::string(isotope);

	        std::ofstream qisooutfile;
	        outfilename = "out/sns_diff_rates_quenched-"+std::string(jsonfile)+"-"+detresp->material+"-"+ffname+"-"+isoname+".out";

	        std::cout << outfilename << std::endl;
	        qisooutfile.open(outfilename);

	        int ie;
	        for (ie=0;ie<iq;ie++) {

                if (Eee[is][ie]>maxeee) {maxeee = Eee[is][ie];}

	            qisooutfile << Eee[is][ie]<< "  "<<dNdEee[is][ie]<<std::endl;
	            _quenchedmap[is][Eee[is][ie]] = dNdEee[is][ie];
	        }
	        qisooutfile.close();
	        v++;is++;
      }

      // Retrieve the smearing function, which should have Gaussian sigma as a function of Eee
      std::string gsname = j["detectorresponse"]["gsname"];

      // Now interpolated rates for quenched, summed over components

      std::ofstream qoutfile;
      outfilename = "out/sns_diff_rates_quenched-alliso-"+std::string(jsonfile)+"-"+detresp->material+"-"+ffname+".out";

      std::cout << outfilename << std::endl;
      qoutfile.open(outfilename);

      // Fixed number of steps in quenched energy
      int  ieee = 0;

      //      int neeebin = j["detectorresponse"]["neeebin"];
      int neeebin = iq;

      double eeestep = maxeee/neeebin;
      double eee = 0;
      double dndeee=0.;

      double nquenchedtot=0;

      for (ieee=0;ieee<neeebin;ieee++) {

	eee += eeestep;
	_quenchedtot[eee] = 0.;

	// Now loop over components
	v = detresp->isotope_component.begin();

	is=0;
	while( v != detresp->isotope_component.end()) {

	  isotope = *v;
	  //	  std::cout << "isotope"<< isotope << std::endl;
	  std::string isoname = std::string(isotope);

	  // Interpolate dNdEee value for isotope is, at this eee
	  // Should encapsulate this in an interpolation routine

	  typedef std::map<double, double>::const_iterator i_t;

	  i_t i=_quenchedmap[is].upper_bound(eee);
	  if(i==_quenchedmap[is].end())
	    {
	      dndeee = (--i)->second;
	    } else if (i==_quenchedmap[is].begin())
	    {
	      dndeee =  i->second;
	    } else {
	    i_t l=i; --l;

	    const double delta=(eee- l->first)/(i->first - l->first);
	    dndeee= delta*i->second +(1-delta)*l->second;
	  }
	  if (::isnan(dndeee)) {dndeee=0.;}

	  _quenchedtot[eee] += dndeee;

	  v++;is++;

	} // End of loop over components


	  // Apply the Eee efficiency here, if requested and not smeared
	double eee_eff_factor = 0.;
	if (gsname == "none"  && detresp->eff_type == "eee"){

	  if (eee>=detresp->eethresh &&
              (detresp->eeupperthresh > detresp->eethresh ? eee <= detresp->eeupperthresh : true)) {
	    if (detresp->effname != "none") {
	      eee_eff_factor = detresp->efficnum(eee);
	    } else {
	      eee_eff_factor = 1;
	    }
	  }

	} else {
	  eee_eff_factor = 1;
	}


      // Now output the total quenched output, per MeVee
	qoutfile <<eee<<" "<<_quenchedtot[eee]*eee_eff_factor<<std::endl;
	std::cout.unsetf(ios::fixed | ios::scientific);

	nquenchedtot += _quenchedtot[eee]*eee_eff_factor*eeestep;



      } // End of loop over Eee

      qoutfile.close();


      std::cout << "Total quenched: "<<nquenchedtot<<std::endl;


      // Now pad the end of the quenched array to allow for smearing

      for (ieee=neeebin;ieee<neeebin*2;ieee++) {
	eee += eeestep;
	_quenchedtot[eee] = 0.;
      }


    // Now do Gaussian smearing, if requested.  Quenching must be requested also if this is to be invoked.  Apply also efficiency here


      if (gsname != "none") {


	// Read the smearing parameters from the file and set them
	std::string gsfilename;
	gsfilename = "gs/"+std::string(gsname)+"_gs.txt";

	DetectorResponse* gs = new DetectorResponse();


	gs->SetGSPolyFilename(gsfilename.c_str());
	gs->ReadGSPolyFile();

	gs->SetMaxSmearEn(maxeee*2);
	gs->SetNSmearBin(neeebin*2);
	gs->SetGaussSmearingMatrix();

	// Do the smearing
	std::map<double,double> _smearedmap= gs->Smear(_quenchedtot);

      // Output the smeared output file

	std::ofstream smoutfile;
	outfilename = "out/sns_diff_rates_smeared-alliso-"+std::string(jsonfile)+"-"+detresp->material+"-"+ffname+".out";

	std::cout << outfilename << std::endl;
	smoutfile.open(outfilename);

	double eeei=0.;
	double totweff = 0.;
	for (ieee=0;ieee<neeebin*2;ieee++) {

	  // Apply the Eee efficiency here, if requested
	  double eee_eff_factor = 1.;
	  if (eeei>=detresp->eethresh &&
              (detresp->eeupperthresh > detresp->eethresh ? eeei <= detresp->eeupperthresh : true)) {
	    if (detresp->effname != "none" && detresp->eff_type == "eee"){
	      eee_eff_factor = detresp->efficnum(eeei);
	    }
	    eeei += eeestep;
	    smoutfile << eeei<<" "<<_smearedmap[eeei]*eee_eff_factor<<std::endl;
	    totweff += _smearedmap[eeei]*eee_eff_factor*eeestep;
	  } else {
	    eeei += eeestep;
	  }
	}

	smoutfile.close();
	std::cout << "Total smeared with efficiency: "<<totweff<<std::endl;

      } // End of do-smearing case

    ////////////  QC-LEVEL RESPONSE ///////////////
    // For this, require quenching also
      // Collected charge, either pe or ADC, etc.


      // In principle this can be a nonlinear response function... linear for now
      double qcperkeVee = j["detectorresponse"]["qcperkeVee"];

      // The total qc distribution (all components... could break it out by isotope)
      std::map<double, double> _qcmapall;
      //Eee in MeVee

      double qcperMeVee = qcperkeVee*1000.;
      double maxmeanqc = maxeee*qcperMeVee;
      double maxqc = maxeee*qcperMeVee*2; // for smearing matrix

      //      std::cout << "maxqc "<<maxqc<<" maxmeanqc "<<maxmeanqc<<std::endl;
      // Loop over qc's, as integers

      // Do the Poisson qc smearing if requested

      std::string qcsmearing = j["detectorresponse"]["qcsmearing"];

      if (qcsmearing != "none") {

	int iqc;
	double qc;

	double totinqc = 0.;
	int qcbinning = j["detectorresponse"]["qcbinning"];

	for (iqc=0;iqc<=int(maxqc);iqc+=qcbinning) {

	  // Interpolate dNdqc from the quenchedmap

	  qc = double(iqc);
	  if (qc<=maxmeanqc+0.01) {


	    typedef std::map<double, double>::const_iterator i_t;

	    //	  double mevee = qc/qcperMeVee;

	  // Do a more fine-grained interpolation and integrate over qc bin,
	  // to reduce binned integration error
	  // Not always really necessary

	    double fracqc;
	    double qcstep = 0.1*qcbinning;

	    double dndqcinbin=0.;
	    double dndqcinterp;
	    double dndqc;

	    // if <qc+0.5 includes last point
	    for (fracqc=qc-0.5*qcbinning;fracqc<qc+0.49*qcbinning;fracqc+=qcstep) {
	      if (fracqc>0) {

	      double mevee2 = fracqc/qcperMeVee;

	      i_t i=_quenchedtot.upper_bound(mevee2);

	      if(i==_quenchedtot.end())
		{
		  dndqc = (--i)->second;
		}
	      else if (i==_quenchedtot.begin())
		{
		  dndqc =  i->second;

		} else {

		i_t l=i; --l;

		const double delta=(mevee2- l->first)/(i->first - l->first);
		dndqc = delta*i->second +(1-delta)*l->second;
	      }
	      dndqcinbin += dndqc*qcstep;

	      //	      std::cout <<qc <<" "<<fracqc<<" "<<mevee2<<" "<<dndqc*qcstep<<" "<<dndqcinbin<<std::endl;

	      } // End of >0 qc case

	  } // End of loop over fractional qc integration
	    //	std::cout <<qc <<" "<<mevee<<" "<<dndqcinbin<<std::endl;


	    dndqcinterp = dndqcinbin/qcperMeVee;

	    if (::isnan(dndqcinterp)) {dndqcinterp=0.;}

	  //	if (qc>0) {
	    totinqc += dndqcinterp;
	    //	}
	    _qcmapall[qc] = dndqcinterp;

	    //	    std::cout << qc<<" "<<dndqcinterp<<" "<<totinqc<<std::endl;

	  } else {

	    // Pad the end with zeroes to allow for smearing
	    _qcmapall[qc] = 0.;

	  } // end of qc<=maxmeanqc check


	}

	std::cout << "Integral of qc dist, including zero bin: "<<totinqc<<std::endl;


	// Poisson smear includes the zero bin
	DetectorResponse* qcsmear = new DetectorResponse();

	qcsmear->SetQCBinning(qcbinning);
	qcsmear->SetNSmearBin(int(maxqc/qcbinning)+1);
	qcsmear->SetMaxSmearEn(double(int(maxqc)+1));

	if (qcsmearing == "poisson") {
	  qcsmear->SetPoissonSmearingMatrix();

	}  else if (qcsmearing == "gamma") {
	    double gammasmearpars[2];
	    gammasmearpars[0] = j["detectorresponse"]["aparam"];
	    gammasmearpars[1] = j["detectorresponse"]["bparam"];
	    qcsmear->SetGammaSmearPars(gammasmearpars);
	    qcsmear->SetGammaSmearingMatrix();
	}  else {

	    std::string qcgsname = j["detectorresponse"]["qcgsname"];

	    std::string qcsmearfilename;
	    qcsmearfilename = "gs/"+std::string(qcgsname)+"_qcsmear.txt";
	    qcsmear->SetGSPolyFilename(qcsmearfilename.c_str());
	    qcsmear->ReadGSPolyFile();
	    qcsmear->SetGaussSmearingMatrix();
	}


	std::cout << "do qc smearing" <<std::endl;
	// Do the smearing
	std::map<double,double> _smearedqcmap = qcsmear->Smear(_qcmapall);


	// Output the qc distribution, applying efficiency if requested
	// (should not also have recoil or quenched efficiency)

	std::ofstream qcoutfile;
	outfilename = "out/sns_diff_rates_qc-alliso-"+std::string(jsonfile)+"-"+detresp->material+"-"+ffname+".out";

	std::cout << outfilename << std::endl;
	qcoutfile.open(outfilename);

	double totev = 0.;
	double totevunsmeared = 0.;

	for (iqc=0;iqc<=int(maxqc);iqc+=qcbinning) {

	// Apply the qc efficiency here, if requested

	  qc = double(iqc);
	  double qc_eff_factor = 1.;

	  if (iqc>=detresp->qcthresh &&
              (detresp->qcupperthresh > detresp->qcthresh ? iqc <= detresp->qcupperthresh : true)) {
	    if (detresp->effname != "none" && detresp->eff_type == "qc"){
	      qc_eff_factor = detresp->efficnum(qc);
	    }

	    qcoutfile << iqc <<" "<<_smearedqcmap[qc]<<" "<<_smearedqcmap[qc]*qc_eff_factor<<" "<<_qcmapall[qc]<<" "<<_qcmapall[qc]*qc_eff_factor<<std::endl;

	    //	    std::cout << iqc <<" "<<_smearedqcmap[qc]<<" "<<_smearedqcmap[qc]*qc_eff_factor<<" "<<_qcmapall[qc]<<" qc "<<qc<<" eff "<<qc_eff_factor<<" "<<_qcmapall[qc]*qc_eff_factor<<std::endl;

	    // It's events per qc bin
	    totev += _smearedqcmap[qc]*qc_eff_factor;
	    totevunsmeared += _qcmapall[qc]*qc_eff_factor;
	  }
	}


	qcoutfile.close();


	cout << "Total qc events: "<<totev<<" unsmeared "<<totevunsmeared<<endl;

      } // End of do qc smearing case

    }  // End of do-quenching case


    // Total flux-averaged xscn

    double dist = j["distance"];
    Ntargets = 1.e6/(Mtot/amu)*6.0221409e23*detresp->detector_mass;
    double totalnus = snsflux->nuspersecperflavor*(snsflux->wnue+snsflux->wnumu+snsflux->wnumubar)/(4*M_PI*dist*dist)*snsflux->exposure;
    double totalnue = snsflux->nuspersecperflavor*snsflux->wnue/(4*M_PI*dist*dist)*snsflux->exposure;
    double totalnumu = snsflux->nuspersecperflavor*snsflux->wnumu/(4*M_PI*dist*dist)*snsflux->exposure;
    double totalnumubar = snsflux->nuspersecperflavor*snsflux->wnumubar/(4*M_PI*dist*dist)*snsflux->exposure;


    std::cout << "Total flux-averaged cross section: "<< totevents/Ntargets/totalnus*1e40<<" x 10-40 cm^2"<<std::endl;
    std::cout << "Total flux-averaged cross section, nue: "<< toteventsnue/Ntargets/totalnue*1e40<<" x 10-40 cm^2"<<std::endl;
    std::cout << "Total flux-averaged cross section, numu: "<< toteventsnumu/Ntargets/totalnumu*1e40<<" x 10-40 cm^2"<<std::endl;
    std::cout << "Total flux-averaged cross section, numubar: "<< toteventsnumubar/Ntargets/totalnumubar*1e40<<" x 10-40 cm^2"<<std::endl;
    std::cout << "Total flux-averaged cross section, numu+numubar: "<< (toteventsnumu+toteventsnumubar)/Ntargets/(totalnumu+totalnumubar)*1e40<<" x 10-40 cm^2"<<std::endl;

    delete[] Er;
    delete[] Eee;
    delete[] dNdEee;
    delete[] dNdEr;
    delete[] dNdErall;

  return 0;

}
