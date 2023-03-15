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

    GenieAnalysisDeltaStudies ga1{input_file.c_str(),
                                  (output_file + "_FSt_new.root").c_str(),
                                  {},
                                  {},
                                  {},
                                  {
                                      {"W", "Ws", {}},
                                      {"W", "el_p", {}},
                                      {"W", "bjorken_x", {}},
                                      {"Ws", "bjorken_x", {}},
                                      {"passed_num_pip", "passed_num_pim", {}},
                                      {"fs_num_pim", "passed_num_pim", {}},
                                      {"fs_num_pip", "passed_num_pip", {}},
                                      {"fs_num_pip", "fs_num_neutrons", {}},
                                      {"fs_num_pim", "fs_num_neutrons", {}},
                                      {"fs_num_pip", "fs_num_protons", {}},
                                      {"fs_num_pim", "fs_num_protons", {}},
                                  }};
    ga1.runAnalysis();

    return 0;
}
