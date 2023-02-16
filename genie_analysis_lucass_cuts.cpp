#include <TVector3.h>
#include <iostream>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>

#include "GenieAnalysis/Fiducial.h"
#include "GenieAnalysis/GenieAnalysisAuto.h"
#include "GenieAnalysis/misc.h"

class FiducialWrapper {
  private:
    Fiducial m_fiducial;

  public:
    FiducialWrapper() : m_fiducial{} {
        // Set up fiducial for 2.261Gev and carbon 12 target
        m_fiducial.InitPiMinusFit("2261");
        m_fiducial.InitEClimits();

        m_fiducial.SetConstants(2250, "12C", {{"1161", 1.161}, {"2261", 2.261}, {"4461", 4.461}});
        m_fiducial.SetFiducialCutParameters("2261");
    }

    bool electron_cut(TVector3 electron_3V) { return m_fiducial.EFiducialCut("2261", electron_3V); }
};

// So far I'm only focusing on exactly the stuff used with the command line arguments in the README,
// specifically this means 2.261GeV beam on C12 target and more.
class GenieAnalysisLucassCuts : public GenieAnalysisAutoTH1Fs {
  private:
    // Paths and other system stuff
    // TODO: Move these to the constructor in some nice way, should not just be hardcoded
    const char *m_el_acceptance_filename{"original/e2a_maps/e2a_maps_12C_E_2_261.root"};
    const std::unique_ptr<TFile> m_el_acceptance_file{TFile::Open(m_el_acceptance_filename, "READ")};
    const std::unique_ptr<TH3D> m_el_acceptance_acc{(TH3D *)m_el_acceptance_file->Get("Accepted Particles")};
    const std::unique_ptr<TH3D> m_el_acceptance_gen{(TH3D *)m_el_acceptance_file->Get("Generated Particles")};

    FiducialWrapper m_fiducials;

    // Parameters
    static constexpr double m_smearing_reso_p{0.01};   // smearing for the proton
    static constexpr double m_smearing_reso_el{0.005}; // smearing for the electrons
    static constexpr double m_smearing_reso_pi{0.007}; // smearing for pions, executive decision by Larry (28.08.19)

    // Properties of the loaded event to be set and used in passesCuts and useEntry
    TVector3 m_smeared_el_V3;
    TLorentzVector m_smeared_el_V4;

  public:
    using GenieAnalysisAutoTH1Fs::GenieAnalysisAutoTH1Fs;
    GenieAnalysisLucassCuts(const char *filename, const char *output_filename, const vector<string> &properties,
                            const vector<string> &types, const char *gst_ttree_name = "gst")
        : GenieAnalysisAutoTH1Fs(filename, output_filename, properties, types, gst_ttree_name), m_fiducials{} {}

    bool passesCuts() {
        const double smeared_pl{gRandom->Gaus(m_ge.pl, m_smearing_reso_el * m_ge.pl)};
        const double smeared_El{sqrt(smeared_pl * smeared_pl + e_mass * e_mass)};

        m_smeared_el_V3.SetXYZ(smeared_pl / m_ge.pl * m_ge.pxl, smeared_pl / m_ge.pl * m_ge.pyl,
                               smeared_pl / m_ge.pl * m_ge.pzl);
        m_smeared_el_V4.SetPxPyPzE(m_smeared_el_V3.X(), m_smeared_el_V3.Y(), m_smeared_el_V3.Z(), smeared_El);

        // GENIE coordinate system flipped with respect to CLAS -- blindly taken from original
        m_smeared_el_V3.SetPhi(m_smeared_el_V3.Phi() + TMath::Pi());

        // Electron theta and momentum fiducial (I think) cut, the values are specifically for C12 ~2GeV set by
        // inspecting
        // https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing and
        // previous values
        const double theta_min{TMath::DegToRad() * (16 + 10.5 / m_smeared_el_V3.Mag())};
        if ((m_smeared_el_V3.Theta() < theta_min) || (m_smeared_el_V3.Theta() < 0) || (m_smeared_el_V3.Theta() > 60)) {
            return false;
        }

        // This was originally later on but I think it makes more sense here
        if (!m_fiducials.electron_cut(m_smeared_el_V3)) {
            return false;
        }

        // Filter some specific sectors, this was enabled and probably makes some sense, essentially only
        // use sectors 0, 1 and 5
        // Sectors are a bit weird, the first goes from -30 to 30 and the rest go as expected to 330
        double temp{(m_smeared_el_V3.Phi() + TMath::Pi() / 6)};
        if (temp < 0) {
            temp += 2 * TMath::Pi();
        }
        const int ElectronSector = temp / (TMath::Pi() / 3);
        if ((ElectronSector >= 2) && (ElectronSector <= 4)) {
            return false;
        }

        // ## There was a cut on electron phi here but unused given the command line options
        // ## Cut on electron sectors here, but unused given our command line options

        // Here there's a longish part calculating weights but the "Mott_cross_section" is hardcoded
        // to 1 and wghts are I think all 1 so essentially it's just the e acceptance weight

        // Calculation of kinematic quantities (nu, Q2, x bjorken, q and W) -- mostly taken from original though
        // specific 2Gev beam
        const TLorentzVector el_change{m_smeared_el_V4 - TLorentzVector{0, 0, 2.261, 2.261}};
        const double reco_Q2 = -el_change.Mag2();

        const double nu = -el_change.E();
        const double x_bjk = reco_Q2 / (2 * m_prot * nu);

        // Get the electron acceptance weight from the e2a map
        const double electron_acceptance_weight{
            electronAcceptance(smeared_pl, m_smeared_el_V3.CosTheta(),
                               m_smeared_el_V3.Phi() + TMath::Pi())}; // could be issue here - angles and CoTheta

        return true;
    }

    double electronAcceptance(const double &p, const double &cos_theta, double phi) {
        // map 330 till 360 to [-30:0] for the acceptance map histogram
        if (phi > (2 * TMath::Pi() - TMath::Pi() / 6.)) {
            phi -= 2 * TMath::Pi();
        }

        int redef = 0; // or -30 I think, it was this way in the original

        // Find number of generated events
        double pbin_gen = m_el_acceptance_gen->GetXaxis()->FindBin(p);
        double tbin_gen = m_el_acceptance_gen->GetYaxis()->FindBin(cos_theta);
        double phibin_gen = m_el_acceptance_gen->GetZaxis()->FindBin(phi * 180 / TMath::Pi() + redef);
        double num_gen = m_el_acceptance_gen->GetBinContent(pbin_gen, tbin_gen, phibin_gen);

        // Find number of accepted events
        double pbin_acc = m_el_acceptance_acc->GetXaxis()->FindBin(p);
        double tbin_acc = m_el_acceptance_acc->GetYaxis()->FindBin(cos_theta);
        double phibin_acc = m_el_acceptance_acc->GetZaxis()->FindBin(phi * 180 / TMath::Pi() + redef);
        double num_acc = m_el_acceptance_acc->GetBinContent(pbin_acc, tbin_acc, phibin_acc);

        double acc_ratio = num_acc / num_gen;
        double acc_err = sqrt(acc_ratio * (1 - acc_ratio)) / num_gen;

        return acc_ratio;
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
