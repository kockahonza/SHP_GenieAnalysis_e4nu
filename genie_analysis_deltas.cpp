#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <tuple>

#include "GenieAnalysis/GenieAnalysisDeltaStudies.h"
#include "GenieAnalysis/misc.h"

int main(int argc, char *argv[]) {
    string input_file, output_file;
    if (argc == 2) {
        std::string arg{argv[1]};
        if (arg == "local") {
            input_file = "/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root";
            output_file = "output_local.root";
        } else if (arg == "full") {
            input_file = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/"
                         "apapadop_SuSav2_C12_2261GeV_master.root";
            output_file = "output_full.root";
        }
    } else {
        std::cout << "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code"
                  << std::endl;
        return -1;
    }

    GenieAnalysis1PionStaged ga{input_file.c_str(),
                                output_file.c_str(),
                                GenieAnalysisDeltaStudies::Target::C12,
                                GenieAnalysisDeltaStudies::BeamEnergy::MeV_2261,
                                {"nocut", "π+", "π-"},
                                {"W", "wght", "el_phi", "el_cos_theta", "el_p", "el_E", "el_acceptance", "pi_phi",
                                 "pi_cos_theta", "pi_p", "pi_E", "pi_acceptance"},
                                {"ALL", "QE", "RES_ALL", "DELTA1232", "DIS"}};

    ga.m_do_precuts = false;
    ga.m_do_electron_fiducials = false;
    ga.m_do_pion_fiducials = false;
    ga.m_do_photon_fiducials = false;

    ga.m_p_pion_momentum_threshold = 0;

    ga.runAnalysis();

    return 0;
}
