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

            GenieAnalysis1Pion ga{
                input_file.c_str(), "output_local_FSt_pip.root", GenieAnalysis1Pion::PionType::Plus, {"nocut"}};
            ga.runAnalysis();

            GenieAnalysis1Pion ga2{
                input_file.c_str(), "output_local_FSt_pim.root", GenieAnalysis1Pion::PionType::Minus, {"nocut"}};
            ga2.runAnalysis();

            GenieAnalysisDeltaStudies ga3{input_file.c_str(), "output_local_FSt_e.root"};
            ga3.runAnalysis();

            GenieAnalysisDeltaStudies ga4{input_file.c_str(), "output_local_FSt_e_noa.root"};
            ga4.m_do_electron_acceptance = false;
            ga4.m_do_pion_acceptance = false;
            ga4.runAnalysis();
        } else if (arg == "full") {
            input_file = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/"
                         "apapadop_SuSav2_C12_2261GeV_master.root";

            /* GenieAnalysis1Pion gap{ */
            /*     input_file.c_str(), "output_full_pip.root", {}, properties, types,
             * GenieAnalysis1Pion::PionType::Plus}; */
            /* gap.runAnalysis(); */

            /* GenieAnalysis1Pion gam{ */
            /*     input_file.c_str(), "output_full_pim.root", {}, properties, types,
             * GenieAnalysis1Pion::PionType::Minus}; */
            /* gam.runAnalysis(); */

            /* GenieAnalysis1Pion gae{input_file.c_str(), */
            /*                        "output_full_pie.root", */
            /*                        {}, */
            /*                        properties, */
            /*                        types, */
            /*                        GenieAnalysis1Pion::PionType::Either}; */
            /* gae.runAnalysis(); */
        }
    } else {
        std::cout << "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code"
                  << std::endl;
        return -1;
    }

    /* ga.m_do_precuts = false; */
    /* ga.m_do_electron_fiducials = false; */
    /* ga.m_do_pion_fiducials = false; */
    /* ga.m_do_photon_fiducials = false; */

    /* ga.m_p_pion_momentum_threshold = 0; */

    return 0;
}
