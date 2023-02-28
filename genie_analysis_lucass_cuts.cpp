#include <iostream>

#include "GenieAnalysis/GenieAnalysisOriginalCuts.h"

int main(int argc, char *argv[]) {
    GenieAnalysisOriginalCuts ga{"/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root",
                                 "output_lucass_cuts_allcuts.root",
                                 GenieAnalysisOriginalCuts::Target::C12,
                                 GenieAnalysisOriginalCuts::BeamEnergy::MeV_2261,
                                 {"nocut"},
                                 {"W", "el_smeared_mag", "el_smeared_phi", "el_smeared_cos_theta", "el_acceptance"},
                                 {"ALL", "QE", "RES", "DIS"}};

    ga.runAnalysis();

    return 0;
}
