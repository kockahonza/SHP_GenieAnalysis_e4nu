#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <tuple>

#include "GenieAnalysis/GACLA6Data.h"
#include "GenieAnalysis/misc.h"

int main(int argc, char *argv[]) {
    string input_file, output_file;

    if (argc == 2) {
        std::string arg{argv[1]};
        if (arg == "local") {
            input_file = "/home/honza/Sync/University/CurrentCourses/SHP/data/Skimmed_1M_e4vWorkshop_C12_2261.root";
            output_file = "jlab_local";
        } else if (arg == "full") {
            /* input_file = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/" */
            /*              "apapadop_SuSav2_C12_2261GeV_master.root"; */
            /* output_file = "jlab_full"; */
        }
    } else {
        std::cout << "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code"
                  << std::endl;
        return -1;
    }

    GACLAS6Data ga1{
        input_file.c_str(), (output_file + "_FSt_resc.root").c_str(), {}, {}, {"ALL"}, {{"W", "Ws", {}}},
    };
    ga1.runAnalysis();

    return 0;
}