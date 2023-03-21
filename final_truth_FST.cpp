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
            output_file = "output_local";
        } else if (arg == "full") {
            input_file = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/"
                         "apapadop_SuSav2_C12_2261GeV_master.root";
            output_file = "output_full";
        }
    } else {
        std::cout << "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code"
                  << std::endl;
        return -1;
    }

    Final1Pion1NucleonTruth gamp{input_file.c_str(),
                                 ("final_" + output_file + "_0pip1pim1p0n.root").c_str(),
                                 Final1Pion1NucleonTruth::PionType::Minus,
                                 Final1Pion1NucleonTruth::NucleonType::Proton,
                                 {},
                                 {"W", "el_phi", "kaka", "pi_phi"},
                                 {},
                                 {{"pi_resc", "nuc_resc", {}}}};
    gamp.runAnalysis();
}
