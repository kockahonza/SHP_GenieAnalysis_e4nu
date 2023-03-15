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

    GenieAnalysisPiNucleonCounts gapp{input_file.c_str(), (output_file + "_1pip0pim1p0n.root").c_str(), 1, 0, 1, 0};
    gapp.runAnalysis();

    GenieAnalysisPiNucleonCounts gapn{input_file.c_str(), (output_file + "_1pip0pim0p1n.root").c_str(), 1, 0, 0, 1};
    gapn.runAnalysis();

    GenieAnalysisPiNucleonCounts gamp{input_file.c_str(), (output_file + "_0pip1pim1p0n.root").c_str(), 0, 1, 1, 0};
    gamp.runAnalysis();

    GenieAnalysisPiNucleonCounts gamn{input_file.c_str(), (output_file + "_0pip1pim0p1n.root").c_str(), 0, 1, 0, 1};
    gamn.runAnalysis();

    return 0;
}
