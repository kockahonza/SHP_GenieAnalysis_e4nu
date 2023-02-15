#include <iostream>

#include <TRandom.h>
#include <TLorentzVector.h>

#include "GenieAnalysis/GenieAnalysisAuto.h"
#include "GenieAnalysis/misc.h"


class GenieAnalysisLucassCuts : public GenieAnalysisAutoTH1Fs {
  private:
    // Parameters
    double smearing_reso_p{0.01};   // smearing for the proton
    double smearing_reso_el{0.005};  // smearing for the electrons
    double smearing_reso_pi{0.007}; // smearing for pions, executive decision by Larry (28.08.19)

    // Properties of the loaded event to be set and used in passesCuts and useEntry
    TVector3 V3_el;
    TLorentzVector V4_el;


  public:
    using GenieAnalysisAutoTH1Fs::GenieAnalysisAutoTH1Fs;

    bool passesCuts() {
        double SmearedPe{gRandom->Gaus(m_ge.pl, smearing_reso_el * m_ge.pl)};
        double SmearedEe{sqrt(SmearedPe * SmearedPe + e_mass * e_mass)};

        V3_el.SetXYZ(SmearedPe / m_ge.pl * m_ge.pxl, SmearedPe / m_ge.pl * m_ge.pyl, SmearedPe / m_ge.pl * m_ge.pzl);
        V4_el.SetPxPyPzE(V3_el.X(), V3_el.Y(), V3_el.Z(), SmearedEe);
        
    }
};


int main(int argc, char *argv[]) {

    GenieAnalysisLucassCuts ga{"/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root",
                                  "output_lucass_cuts.root",
                                  {"W"},
                                  {"QE", "RES", "DIS"}};

    ga.runAnalysis();

    return 0;
}
