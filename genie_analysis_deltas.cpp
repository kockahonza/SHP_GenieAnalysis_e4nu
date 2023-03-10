#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <tuple>

#include "GenieAnalysis/GenieAnalysisDeltaStudies.h"
#include "GenieAnalysis/misc.h"


int main(int argc, char *argv[]) {
    GenieAnalysis1Pion ga{gst_path_jan,
                          "output_jan.root",
                          GenieAnalysisDeltaStudies::Target::C12,
                          GenieAnalysisDeltaStudies::BeamEnergy::MeV_2261,
                          {"nocut", "π+", "π-"},
                          {"W", "wght", "el_phi", "el_cos_theta", "el_p", "el_E", "el_acceptance", "pi_phi",
                           "pi_cos_theta", "pi_p", "pi_E", "pi_acceptance"},
                          {"ALL", "QE", "RES_ALL", "DELTA1232", "DIS"}};
    /* "PIP", "QE_PIP", "RES_ALL_PIP", "DELTA1232_PIP", "DIS_PIP", */
    /* "PIM", "QE_PIM", "RES_ALL_PIM", "DELTA1232_PIM", "DIS_PIM"}, true, true, false}; */

    /* GenieAnalysisDeltaStudies ga{gst_path_jan, */
    /*                              "output_jan.root", */
    /*                              GenieAnalysisDeltaStudies::Target::C12, */
    /*                              GenieAnalysisDeltaStudies::BeamEnergy::MeV_2261, */
    /*                              {"nocut", "p_gdoc", "p_efid"}, */
    /*                              {"W", "el_phi", "el_cos_theta", "el_p", "el_E", "el_acceptance"}, */
    /*                              {"ALL", "QE", "RES_ALL", "DELTA1232", "DIS"}, false, true, false}; */

    /* GenieAnalysis1Pion ga{genie_machine_path, */
    /*                       "output_jan_full.root", */
    /*                       GenieAnalysisDeltaStudies::Target::C12, */
    /*                       GenieAnalysisDeltaStudies::BeamEnergy::MeV_2261, */
    /*                       {"PIP", "PIM"}, */
    /*                       {"W", "wght", "el_phi", "el_cos_theta", "el_p", "el_E", "el_acceptance", "pi_phi", */
    /*                        "pi_cos_theta", "pi_p", "pi_E", "pi_acceptance"}, */
    /*                       {"ALL", "QE", "RES_ALL", "DELTA1232", "DIS"}}; */

    ga.m_do_precuts = false;
    ga.m_do_electron_fiducials = false;
    ga.m_do_pion_fiducials = false;
    ga.m_do_photon_fiducials = false;

    ga.m_p_pion_momentum_threshold = 0;

    ga.runAnalysis();

    return 0;
}
