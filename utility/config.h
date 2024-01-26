
#include <string>
#include <nlohmann/json.hpp>

#ifndef CONFIG_H
#define CONFIG_H

class Config{
    public:
        int convolved;
        int pot;
        double usernorm;
        double tw1;
        double tw2;
        double teffic_params[3];
        double wnumu;
        double wnumubar;
        double wnue;
        std::string ffname;
        double nuspersecperflavor;
        double dist;
        double ua4[3];
        double dm2;
        std::string eff_type;
        std::string effname;
        double lowerThresh;
        double upperThresh;
        std::string material;
        std::string qftype;
        std::string qfname;
        double nvrfact;
        double narfact;
        double pvrfact;
        double parfact;

        // For Helm
        double nvsfact;
        double nasfact;
        double pvsfact;
        double pasfact;

        // For Klein
        double nvak;
        double naak;
        double pvak;
        double paak;
        double nvskin;
        double naskin;
        double pvskin;
        double paskin;

        // for detector
        double detector_mass;
        double hours;


        double erecstep;
        double knustep;


	    int chgcorrtype;


        std::string gsname;


        Config(nlohmann::json config);






};
#endif // !CONFIG_H
