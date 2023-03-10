#include "GenieAnalysisDeltaStudies.h"
#include <TMath.h>

Double_t GenieAnalysisDeltaStudies::passesCuts() {
    useEntryAtStage("nocut", m_ge.wght);

    // Do electron smearing
    const double smeared_pl{gRandom->Gaus(m_ge.pl, m_smearing_reso_el * m_ge.pl)};
    const double smeared_El{sqrt(smeared_pl * smeared_pl + e_mass * e_mass)};
    m_smeared_el_V3.SetXYZ(smeared_pl / m_ge.pl * m_ge.pxl, smeared_pl / m_ge.pl * m_ge.pyl,
                           smeared_pl / m_ge.pl * m_ge.pzl);
    m_smeared_el_V4.SetPxPyPzE(m_smeared_el_V3.X(), m_smeared_el_V3.Y(), m_smeared_el_V3.Z(), smeared_El);

    // GENIE coordinate system flipped with respect to CLAS, all the 4momentas phis are added pi to -- blindly taken
    // from original
    m_smeared_el_V3.SetPhi(m_smeared_el_V3.Phi() + TMath::Pi());

    // Electron theta and momentum fiducial (essentially I think) cut, the values are specifically for C12 2.261GeV
    // set by inspecting
    // https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing and
    // previous values
    if (m_do_precuts) {
        const double theta_min{TMath::DegToRad() * (m_precut_parameter1 + m_precut_parameter2 / m_smeared_el_V3.Mag())};
        if ((m_smeared_el_V3.Theta() < theta_min) || (m_smeared_el_V3.Theta() < 0 * TMath::DegToRad()) ||
            (m_smeared_el_V3.Theta() > 80 * TMath::DegToRad())) {
            return 0;
        }
        useEntryAtStage("p_gdoc", m_ge.wght);
    }

    // This was originally later on but I think it makes more sense here
    if (m_do_electron_fiducials) {
        if (!m_fiducials->electronCut(m_smeared_el_V3)) {
            return 0;
        }
        useEntryAtStage("p_efid", m_ge.wght);
    }

    // Filter some specific sectors, this was enabled and probably makes some sense, essentially only
    // use sectors 0, 1 and 5, o this only if m_do_sectors is true
    // Sectors are a bit weird, the first goes from -30 to 30 and the rest go as expected to 330
    if (m_do_sectors) {
        double temp{(m_smeared_el_V3.Phi() + TMath::Pi() / 6)};
        if (temp < 0) {
            temp += 2 * TMath::Pi();
        }
        const int ElectronSector = temp / (TMath::Pi() / 3);
        if ((ElectronSector >= 2) && (ElectronSector <= 4)) {
            return 0;
        }
    }

    // Calculation of kinematic quantities (nu, Q2, x bjorken, q and W) -- mostly taken from original though
    // specific to 2.261Gev beam
    // TODO: Unusued now, delete or use
    /* const TLorentzVector el_change{m_smeared_el_V4 - TLorentzVector{0, 0, 2.261, 2.261}}; */
    /* const double reco_Q2 = -el_change.Mag2(); */

    /* const double nu = -el_change.E(); */
    /* const double x_bjk = reco_Q2 / (2 * m_prot * nu); */

    // Get the electron acceptance weight from the e2a map
    m_electron_acceptance_weight = electronAcceptance(smeared_pl, m_smeared_el_V3.CosTheta(), m_smeared_el_V3.Phi());
    if (m_electron_acceptance_weight != TMath::Abs(m_electron_acceptance_weight)) {
        throw "Electron acceptance not reasonable";
    }

    // Hadron loop, first clear the vectors though
    m_passed_pi_minus.clear();
    m_passed_pi_plus.clear();
    m_passed_photons.clear();

    // Temp variables for the loop, declared here for performance
    double smeared_p;
    double smeared_E;
    double pi_acceptance;
    TVector3 V3;
    double positive_phi_difference;
    for (int i{0}; i < m_ge.nf; i++) {
        // pi-
        if (m_ge.pdgf[i] == -211) {
            // required momentum for detection
            if (m_ge.pf[i] < m_p_pion_momentum_threshold) {
                continue;
            }
            smeared_p = gRandom->Gaus(m_ge.pf[i], m_smearing_reso_pi * m_ge.pf[i]);
            smeared_E = sqrt(smeared_p * smeared_p + m_pion * m_pion);

            V3.SetXYZ(smeared_p / m_ge.pf[i] * m_ge.pxf[i], smeared_p / m_ge.pf[i] * m_ge.pyf[i],
                      smeared_p / m_ge.pf[i] * m_ge.pzf[i]);
            V3.SetPhi(V3.Phi() + TMath::Pi());

            if (m_do_pion_fiducials) {
                if (!m_fiducials->piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Minus)) {
                    continue;
                }
            }

            // And only use pirns with reasonable acceptance, quite often they are negative or so on, this is kindof a
            // workaround
            pi_acceptance = piMinusAcceptance(V3.Mag(), V3.CosTheta(), V3.Phi());
            if (pi_acceptance == TMath::Abs(pi_acceptance)) {
                m_passed_pi_minus.push_back({{V3, smeared_E}, V3, pi_acceptance});
            }
        }
        // pi+
        if (m_ge.pdgf[i] == 211) {
            // required momentum for detection
            if (m_ge.pf[i] < m_p_pion_momentum_threshold) {
                continue;
            }
            smeared_p = gRandom->Gaus(m_ge.pf[i], m_smearing_reso_pi * m_ge.pf[i]);
            smeared_E = sqrt(smeared_p * smeared_p + m_pion * m_pion);

            V3.SetXYZ(smeared_p / m_ge.pf[i] * m_ge.pxf[i], smeared_p / m_ge.pf[i] * m_ge.pyf[i],
                      smeared_p / m_ge.pf[i] * m_ge.pzf[i]);
            V3.SetPhi(V3.Phi() + TMath::Pi());

            if (m_do_pion_fiducials) {
                if (!m_fiducials->piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Plus)) {
                    continue;
                }
            }

            // And only use pirns with reasonable acceptance, quite often they are negative or so on, this is kindof a
            // workaround
            pi_acceptance = piMinusAcceptance(V3.Mag(), V3.CosTheta(), V3.Phi());
            if (pi_acceptance == TMath::Abs(pi_acceptance)) {
                m_passed_pi_plus.push_back({{V3, smeared_E}, V3, piPlusAcceptance(V3.Mag(), V3.CosTheta(), V3.Phi())});
            }
        }
        // photons
        if (m_ge.pdgf[i] == 22) {
            // required momentum for detection
            if (m_ge.pf[i] < m_p_photon_momentum_threshold) {
                continue;
            }

            // No photon smearing
            V3.SetXYZ(m_ge.pxf[i], m_ge.pyf[i], m_ge.pzf[i]);
            V3.SetPhi(V3.Phi() + TMath::Pi());

            if (m_do_photon_fiducials) {
                if (!m_fiducials->piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Photon)) {
                    continue;
                }
            }

            // TODO: The following should cut on presence of radiation electrons but it's worth reviewing, a comment
            // in original says "within 40 degrees in theta and 30 in phi" but that's not what it did, even after
            // fixing a clear(ish) bug For the second condition there's a phi angle difference the +2pi is to bring
            // it to the 0 to 4pi range adn then mod 2pi and compare to
            if (m_do_radiation_check) {
                positive_phi_difference = TMath::Abs(V3.Phi() - m_smeared_el_V3.Phi()) * TMath::RadToDeg();
                if ((V3.Angle(m_smeared_el_V3) < m_p_radiation_photon_angle) &&
                    ((positive_phi_difference < m_p_radiation_photon_phi_diff) ||
                     (positive_phi_difference > (360 - m_p_radiation_photon_phi_diff)))) {
                    return 0;
                }
            }

            m_passed_photons.push_back({{V3, m_ge.Ef[i]}, V3, photonAcceptance(V3.Mag(), V3.CosTheta(), V3.Phi())});
        }
    }

    return m_electron_acceptance_weight * m_ge.wght; // I'm pretty sure all the events have wght=1
}

double GenieAnalysisDeltaStudies::acceptanceJoined(const double &p, const double &cos_theta, double phi,
                                                   const std::unique_ptr<TH3D> &generated,
                                                   const std::unique_ptr<TH3D> &accepted) {
    // first map -pi, pi to [0, 2pi
    if (phi < 0) {
        phi += 2 * TMath::Pi();
    }
    phi *= TMath::RadToDeg();
    // map 330 till 360 to [-30:0] but in radians for the acceptance map histogram
    if (phi > (330)) {
        phi -= 360;
    }

    int redef = 0; // or -30 I think, it was this way in the original

    // Find number of generated events
    double pbin_gen = generated->GetXaxis()->FindBin(p);
    double tbin_gen = generated->GetYaxis()->FindBin(cos_theta);
    double phibin_gen = generated->GetZaxis()->FindBin(phi + redef);
    double num_gen = generated->GetBinContent(pbin_gen, tbin_gen, phibin_gen);

    // Find number of accepted events
    double pbin_acc = accepted->GetXaxis()->FindBin(p);
    double tbin_acc = accepted->GetYaxis()->FindBin(cos_theta);
    double phibin_acc = accepted->GetZaxis()->FindBin(phi + redef);
    double num_acc = accepted->GetBinContent(pbin_acc, tbin_acc, phibin_acc);

    double acc_ratio = num_acc / num_gen;
    // TODO: Unused, use or delete soon
    /* double acc_err = sqrt(acc_ratio * (1 - acc_ratio)) / num_gen; */

    return acc_ratio;
}
