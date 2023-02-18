#include <iostream>

#include <TH3D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <memory>

#include "GenieAnalysis/Fiducial.h"
#include "GenieAnalysis/GenieAnalysisAuto.h"
#include "GenieAnalysis/misc.h"

class FiducialWrapper {
  public:
    enum class PiPhotonId : int { Minus = -1, Photon = 0, Plus = 1 };

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

    bool electronCut(TVector3 momentum_V3) { return m_fiducial.EFiducialCut("2261", momentum_V3); }

    bool piAndPhotonCuts(TVector3 momentum_V3, PiPhotonId which_particle) {
        return m_fiducial.Pi_phot_fid_united("2261", momentum_V3, static_cast<int>(which_particle));
    }
};

// So far I'm only focusing on exactly the stuff used with the command line arguments in the README,
// specifically this means 2.261GeV beam on C12 target and more.
class GenieAnalysisLucassCuts : public GenieAnalysisAutoTH1Fs {
  private:
    // Paths and other system stuff
    // TODO: Move these to the constructor in some nice way, should not just be hardcoded
    const std::unique_ptr<TFile> m_el_acceptance_file{
        TFile::Open("original/e2a_maps/e2a_maps_12C_E_2_261.root", "READ")};
    const std::unique_ptr<TFile> m_p_acceptance_file{
        TFile::Open("original/e2a_maps/e2a_maps_12C_E_2_261_p.root", "READ")};
    const std::unique_ptr<TFile> m_pip_acceptance_file{
        TFile::Open("original/e2a_maps/e2a_maps_12C_E_2_261_pip.root", "READ")};
    const std::unique_ptr<TFile> m_pim_acceptance_file{
        TFile::Open("original/e2a_maps/e2a_maps_12C_E_2_261_pim.root", "READ")};

    // It seems to be necessary for I suppose performance reason? that these are found beforehand and not at each call
    // to acceptanceJoined
    const std::unique_ptr<TH3D> m_acc_el_gen, m_acc_el_acc, m_acc_p_gen, m_acc_p_acc, m_acc_pip_gen, m_acc_pip_acc,
        m_acc_pim_gen, m_acc_pim_acc;

    FiducialWrapper m_fiducials;

    // Parameters
    static constexpr double m_smearing_reso_p{0.01};   // smearing for the proton
    static constexpr double m_smearing_reso_el{0.005}; // smearing for the electrons
    static constexpr double m_smearing_reso_pi{0.007}; // smearing for pions, executive decision by Larry (28.08.19)

    // Properties of the loaded event to be set and used in passesCuts and new properties
    // electron
    TVector3 m_smeared_el_V3; // Smeared and rotated by pi
    TLorentzVector m_smeared_el_V4;

    // hadrons -- all of these only contain information on particles passing relevant cuts (pions need momentum above
    // 0.15 for example)
    vector<tuple<TLorentzVector, TVector3, double>>
        m_passed_pi_plus; // tuple has the smeared 4 momentum, smeared 3 momentum and the pions calculated acceptance
    vector<tuple<TLorentzVector, TVector3, double>>
        m_passed_pi_minus; // tuple has the smeared 4 momentum, smeared 3 momentum and the pions calculated acceptance
    vector<tuple<TLorentzVector, TVector3, double>>
        m_passed_photons; // tuple has the smeared 4 momentum, smeared 3 momentum and the pions calculated acceptance

    double electron_acceptance_weight;

  protected:
    // Extensions to automatic TH1Fs
    map<string, AutoProperty> m_new_known_properties{
        {"el_smeared_phi",
         {{720, -30, 330},
          [this]() {
              double phi_deg{m_smeared_el_V3.Phi() * TMath::RadToDeg()};
              return (phi_deg < -30) ? phi_deg + 360 : phi_deg;
          }}},
        {"el_smeared_cos_theta", {{720, -1, 1}, [this]() { return m_smeared_el_V3.CosTheta(); }}},
        {"el_smeared_mag", {{720, 0, 3}, [this]() { return m_smeared_el_V3.Mag(); }}},
        {"el_acceptance", {{100, 0, 1}, [this]() { return electron_acceptance_weight; }}}};

  public:
    GenieAnalysisLucassCuts(const char *filename, const char *output_filename, const vector<string> &stages,
                            const vector<string> &properties, const vector<string> &types,
                            const char *gst_ttree_name = "gst")
        : GenieAnalysisAutoTH1Fs(filename, output_filename, stages, properties, types, gst_ttree_name), m_fiducials{},
          m_acc_el_gen{(TH3D *)m_el_acceptance_file->Get("Generated Particles")},
          m_acc_el_acc{(TH3D *)m_el_acceptance_file->Get("Accepted Particles")},
          m_acc_p_gen{(TH3D *)m_p_acceptance_file->Get("Generated Particles")},
          m_acc_p_acc{(TH3D *)m_p_acceptance_file->Get("Accepted Particles")},
          m_acc_pip_gen{(TH3D *)m_pip_acceptance_file->Get("Generated Particles")},
          m_acc_pip_acc{(TH3D *)m_pip_acceptance_file->Get("Accepted Particles")},
          m_acc_pim_gen{(TH3D *)m_pim_acceptance_file->Get("Generated Particles")},
          m_acc_pim_acc{(TH3D *)m_pim_acceptance_file->Get("Accepted Particles")} {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }

    bool passesCuts() {
        const double smeared_pl{gRandom->Gaus(m_ge.pl, m_smearing_reso_el * m_ge.pl)};
        const double smeared_El{sqrt(smeared_pl * smeared_pl + e_mass * e_mass)};

        m_smeared_el_V3.SetXYZ(smeared_pl / m_ge.pl * m_ge.pxl, smeared_pl / m_ge.pl * m_ge.pyl,
                               smeared_pl / m_ge.pl * m_ge.pzl);
        m_smeared_el_V4.SetPxPyPzE(m_smeared_el_V3.X(), m_smeared_el_V3.Y(), m_smeared_el_V3.Z(), smeared_El);

        // GENIE coordinate system flipped with respect to CLAS -- blindly taken from original
        m_smeared_el_V3.SetPhi(m_smeared_el_V3.Phi() + TMath::Pi());

        useEntryAtStage("nocut");

        // Electron theta and momentum fiducial (essentially I think) cut, the values are specifically for C12 2.261GeV
        // set by inspecting
        // https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing and
        // previous values
        const double theta_min{TMath::DegToRad() * (16 + 10.5 / m_smeared_el_V3.Mag())};
        if ((m_smeared_el_V3.Theta() < theta_min) || (m_smeared_el_V3.Theta() < 0) || (m_smeared_el_V3.Theta() > 60)) {
            return false;
        }

        // This was originally later on but I think it makes more sense here
        if (!m_fiducials.electronCut(m_smeared_el_V3)) {
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
        // specific to 2.261Gev beam
        const TLorentzVector el_change{m_smeared_el_V4 - TLorentzVector{0, 0, 2.261, 2.261}};
        const double reco_Q2 = -el_change.Mag2();

        const double nu = -el_change.E();
        const double x_bjk = reco_Q2 / (2 * m_prot * nu);

        // Get the electron acceptance weight from the e2a map
        electron_acceptance_weight =
            electronAcceptance(smeared_pl, m_smeared_el_V3.CosTheta(),
                               m_smeared_el_V3.Phi() + TMath::Pi()); // could be issue here - angles and CoTheta

        // Hadron loop -- need to add
        for (int i{0}; i < m_ge.nf; i++) {
            // pi-
            if (m_ge.pdgf[i] == -221) {
                // required momentum for detection
                if (m_ge.pf[i] < 0.15) {
                    continue;
                }
                const double smeared_p = gRandom->Gaus(m_ge.pf[i], m_smearing_reso_pi * m_ge.pf[i]);
                const double smeared_E = sqrt(smeared_p * smeared_p + m_pion * m_pion);

                TVector3 V3{smeared_p / m_ge.pf[i] * m_ge.pxf[i], smeared_p / m_ge.pf[i] * m_ge.pyf[i],
                            smeared_p / m_ge.pf[i] * m_ge.pzf[i]};
                V3.SetPhi(V3.Phi() + TMath::Pi());
                TLorentzVector V4{V3, smeared_E};

                if (!m_fiducials.piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Minus)) {
                    continue;
                }

                m_passed_pi_minus.push_back({V4, V3, piMinusAcceptance(V3.Mag(), V3.CosTheta(), V3.Phi())});
            }
            // pi+
            if (m_ge.pdgf[i] == 221) {
                // required momentum for detection
                if (m_ge.pf[i] < 0.15) {
                    continue;
                }
                const double smeared_p = gRandom->Gaus(m_ge.pf[i], m_smearing_reso_pi * m_ge.pf[i]);
                const double smeared_E = sqrt(smeared_p * smeared_p + m_pion * m_pion);

                TVector3 V3{smeared_p / m_ge.pf[i] * m_ge.pxf[i], smeared_p / m_ge.pf[i] * m_ge.pyf[i],
                            smeared_p / m_ge.pf[i] * m_ge.pzf[i]};
                V3.SetPhi(V3.Phi() + TMath::Pi());
                TLorentzVector V4{V3, smeared_E};

                if (!m_fiducials.piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Plus)) {
                    continue;
                }

                m_passed_pi_plus.push_back({V4, V3, piPlusAcceptance(V3.Mag(), V3.CosTheta(), V3.Phi())});
            }
            // photons
            if (m_ge.pdgf[i] == 22) {
                // required momentum for detection
                if (m_ge.pf[i] < 0.3) {
                    continue;
                }

                // No photon smearing
                TVector3 V3{m_ge.pxf[i], m_ge.pyf[i], m_ge.pzf[i]};
                V3.SetPhi(V3.Phi() + TMath::Pi());
                TLorentzVector V4{V3, m_ge.Ef[i]};

                if (!m_fiducials.piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Photon)) {
                    continue;
                }

                // TODO: The following should cut on presence of radiation electrons but it's worth reviewing, a comment
                // in original says "within 40 degrees in theta and 30 in phi" but that's not what it did, even after
                // fixing a clear(ish) bug For the second condition there's a phi angle difference the +2pi is to bring
                // it to the 0 to 4pi range adn then mod 2pi and compare to
                /* if ((V3.Angle(m_smeared_el_V3) < 40) && (TMath::Abs(fmod(V3.Phi() - m_smeared_el_V3.Phi() + 2 *
                 * TMath::Pi(), 2 * TMath::Pi())) < 30)) { */
                const double positive_phi_difference{TMath::Abs(V3.Phi() - m_smeared_el_V3.Phi()) *
                                                     TMath::RadToDeg()}; // will be in the [0, 360) range]
                if ((V3.Angle(m_smeared_el_V3) < 40) &&
                    ((positive_phi_difference < 30) || (positive_phi_difference > 330))) {
                    return false;
                }

                m_passed_photons.push_back({V4, V3, photonAcceptance(V3.Mag(), V3.CosTheta(), V3.Phi())});
            }
        }

        return true;
    }

    double acceptanceJoined(const double &p, const double &cos_theta, double phi,
                            const std::unique_ptr<TH3D> &generated, const std::unique_ptr<TH3D> &accepted) {
        // map 330 till 360 to [-30:0] for the acceptance map histogram
        if (phi > (2 * TMath::Pi() - TMath::Pi() / 6.)) {
            phi -= 2 * TMath::Pi();
        }

        int redef = 0; // or -30 I think, it was this way in the original

        // Find number of generated events
        double pbin_gen = generated->GetXaxis()->FindBin(p);
        double tbin_gen = generated->GetYaxis()->FindBin(cos_theta);
        double phibin_gen = generated->GetZaxis()->FindBin(phi * 180 / TMath::Pi() + redef);
        double num_gen = generated->GetBinContent(pbin_gen, tbin_gen, phibin_gen);

        // Find number of accepted events
        double pbin_acc = accepted->GetXaxis()->FindBin(p);
        double tbin_acc = accepted->GetYaxis()->FindBin(cos_theta);
        double phibin_acc = accepted->GetZaxis()->FindBin(phi * 180 / TMath::Pi() + redef);
        double num_acc = accepted->GetBinContent(pbin_acc, tbin_acc, phibin_acc);

        double acc_ratio = num_acc / num_gen;
        double acc_err = sqrt(acc_ratio * (1 - acc_ratio)) / num_gen;

        return acc_ratio;
    }

    double electronAcceptance(const double &p, const double &cos_theta, const double &phi) {
        return acceptanceJoined(p, cos_theta, phi, m_acc_el_gen, m_acc_el_acc);
    }

    double photonAcceptance(const double &p, const double &cos_theta, const double &phi) {
        return acceptanceJoined(p, cos_theta, phi, m_acc_p_gen, m_acc_p_acc);
    }

    double piPlusAcceptance(const double &p, const double &cos_theta, const double &phi) {
        return acceptanceJoined(p, cos_theta, phi, m_acc_pip_gen, m_acc_pip_acc);
    }

    double piMinusAcceptance(const double &p, const double &cos_theta, const double &phi) {
        return acceptanceJoined(p, cos_theta, phi, m_acc_pim_gen, m_acc_pim_acc);
    }
};

int main(int argc, char *argv[]) {
    GenieAnalysisLucassCuts ga{"/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root",
                               "output_lucass_cuts_allcuts.root",
                               {"nocut"},
                               {"W", "el_smeared_mag", "el_smeared_phi", "el_smeared_cos_theta", "el_acceptance"},
                               {"ALL", "QE", "RES", "DIS"}};

    ga.runAnalysis();

    return 0;
}
