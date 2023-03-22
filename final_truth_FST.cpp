#include <iostream>
#include <limits>
#include <optional>
#include <tuple>
#include <utility>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

#include "GenieAnalysis/GACLAS6FinalFST.h"

int main(int argc, char *argv[]) {
    string input_file, output_file;
    if (argc == 2) {
        std::string arg{argv[1]};
        if (arg == "local") {
            input_file = "/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root";
            output_file = "local";
        } else if (arg == "full") {
            input_file = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/"
                         "apapadop_SuSav2_C12_2261GeV_master.root";
            output_file = "full";
        }
    } else {
        std::cout << "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code"
                  << std::endl;
        return -1;
    }

    /* Final1Pion1NucleonTruth ga{input_file.c_str(), */
    /*                            ("final_test_" + output_file + ".root").c_str(), */
    /*                            Final1Pion1NucleonTruth::PionType::Minus, */
    /*                            Final1Pion1NucleonTruth::NucleonType::Proton, */
    /*                            Final1Pion1NucleonTruth::RunType::PrimaryState, */
    /*                            {}, */
    /*                            {}, */
    /*                            {"ALL", "QE", "DELTA1232", "RES_OTHER", "DIS"}, */
    /*                            { */
    /*                                {"pi_phi", "reco_pi_phi", {}}, */
    /*                                {"pi_ct", "reco_pi_ct", {}}, */
    /*                                {"pi_p", "reco_pi_p", {}}, */
    /*                                {"pi_E", "reco_pi_E", {}}, */
    /*                                {"pi_resc", "nuc_resc", {}}, */
    /*                            }}; */
    /* ga.runAnalysis(); */

    vector<tuple<string, Final1Pion1NucleonTruth::PionType, Final1Pion1NucleonTruth::NucleonType>> variants{
        {"_0pip1pim1p0n", Final1Pion1NucleonTruth::PionType::Minus, Final1Pion1NucleonTruth::NucleonType::Proton},
        {"_1pip0pim0p1n", Final1Pion1NucleonTruth::PionType::Plus, Final1Pion1NucleonTruth::NucleonType::Neutron},

        {"_1pip0pim1p0n", Final1Pion1NucleonTruth::PionType::Plus, Final1Pion1NucleonTruth::NucleonType::Proton},
        {"_0pip1pim0p1n", Final1Pion1NucleonTruth::PionType::Minus, Final1Pion1NucleonTruth::NucleonType::Neutron},
    };

    for (auto const &[ext, pi_t, nuc_t] : variants) {
        Final1Pion1NucleonTruth ga{input_file.c_str(),
                                   ("final_ps_" + output_file + ext + ".root").c_str(),
                                   pi_t,
                                   nuc_t,
                                   Final1Pion1NucleonTruth::RunType::PrimaryState,
                                   {},
                                   {},
                                   {"ALL", "QE", "DELTA1232", "RES_OTHER", "DIS"},
                                   {
                                       {"pi_phi", "reco_pi_phi", {}},
                                       {"pi_ct", "reco_pi_ct", {}},
                                       {"pi_p", "reco_pi_p", {}},
                                       {"pi_E", "reco_pi_E", {}},
                                       {"pi_resc", "nuc_resc", {}},
                                   }};
        ga.runAnalysis();
    }

    for (auto const &[ext, pi_t, nuc_t] : variants) {
        Final1Pion1NucleonTruth ga{input_file.c_str(),
                                   ("final_fsr_" + output_file + ext + ".root").c_str(),
                                   pi_t,
                                   nuc_t,
                                   Final1Pion1NucleonTruth::RunType::FinalStateResc1,
                                   {},
                                   {},
                                   {"ALL", "QE", "DELTA1232", "RES_OTHER", "DIS"},
                                   {
                                       {"pi_phi", "reco_pi_phi", {}},
                                       {"pi_ct", "reco_pi_ct", {}},
                                       {"pi_p", "reco_pi_p", {}},
                                       {"pi_E", "reco_pi_E", {}},
                                   }};
        ga.runAnalysis();
    }

    for (auto const &[ext, pi_t, nuc_t] : variants) {
        Final1Pion1NucleonTruth ga{input_file.c_str(),
                                   ("final_fs_" + output_file + ext + ".root").c_str(),
                                   pi_t,
                                   nuc_t,
                                   Final1Pion1NucleonTruth::RunType::FinalState,
                                   {},
                                   {},
                                   {"ALL", "QE", "DELTA1232", "RES_OTHER", "DIS"},
                                   {
                                       {"pi_phi", "reco_pi_phi", {}},
                                       {"pi_ct", "reco_pi_ct", {}},
                                       {"pi_p", "reco_pi_p", {}},
                                       {"pi_E", "reco_pi_E", {}},
                                       {"resc_same", "reco_pi_p_test", {}},
                                       {"resc_same", "p_change_mag", {}},
                                       {"resc_same", "E_change", {}},
                                   }};
        ga.runAnalysis();
    }
}
