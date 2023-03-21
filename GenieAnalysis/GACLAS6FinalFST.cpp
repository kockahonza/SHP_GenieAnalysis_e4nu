#include "GACLAS6FinalFST.h"

Double_t ElectronFiducials::passesCuts() {
    // Electron 4 momentum, rotate to be the same as CLAS6
    m_el_V4.SetPxPyPzE(m_ge.pxl, m_ge.pyl, m_ge.pzl, m_ge.El);
    m_el_V4.SetPhi(m_el_V4.Phi() + TMath::Pi());

    // Electron theta and momentum fiducial (essentially I think) cut, the values are specifically for C12 2.261GeV
    // set by inspecting
    // https://docs.google.com/presentation/d/1ghG08JfCYXRXh6O8hcXKrhJOFxkAs_9i5ZfoIkiiEHU/edit?usp=sharing and
    // previous values
    const double theta_min{TMath::DegToRad() * (m_el_precut_parameter1 + m_el_precut_parameter2 / m_el_V4.P())};
    if ((m_el_V4.Theta() < theta_min) || (m_el_V4.Theta() > 80 * TMath::DegToRad())) {
        return 0;
    }
    // Electron fiducials
    if (!m_fiducials->electronCut(m_el_V4.Vect())) {
        return 0;
    }

    // Calculation of kinematic quantities (nu, Q2, x bjorken, q and W) -- literally taken from original though
    const TLorentzVector el_change{m_el_V4 - TLorentzVector{0, 0, m_beam_energy_val, m_beam_energy_val}};
    m_reco_Q2 = -el_change.Mag2();
    const double nu = -el_change.E();
    m_reco_x = m_reco_Q2 / (2 * mass_proton * nu);
    m_reco_W = TMath::Sqrt((mass_proton + nu) * (mass_proton + nu) - el_change.Vect().Mag2());

    if (m_Wcut) {
        if ((m_reco_W < m_Wcut->first) || (m_reco_W > m_Wcut->second)) {
            std::cout << "Wcut trigerred" << std::endl;
            return 0;
        }
    }

    return m_ge.wght;
}

Double_t Final1Pion1NucleonTruth::passesCuts() {
    Double_t weight{ElectronFiducials::passesCuts()};
    if (weight == 0) {
        return 0;
    }

    m_found_pi = false;
    m_found_nuc = false;
    // Temp variables for the primary state hadrons(and photon) loop, declared here for performance, these are the
    // most interesting hadrons for comparison to experiment
    TVector3 V3;
    for (int i{0}; i < m_ge.ni; i++) {
        if (m_ge.pdgi[i] == 2212) { // proton
            V3.SetXYZ(m_ge.pxi[i], m_ge.pyi[i], m_ge.pzi[i]);
            V3.SetPhi(V3.Phi() + TMath::Pi());

            if (V3.Mag() < m_p_proton_momentum_threshold) {
                continue;
            }

            if (!m_fiducials->protonCut(V3)) {
                continue;
            }

            // If we got here, it means it qualifies as signal, if we are still looking for a proton, make this
            // found, if we are looking for neutrons or have already found a proton, the event is invalid so return
            // 0
            if ((m_nuc_type == NucleonType::Proton) && (!m_found_nuc)) {
                m_found_nuc = true;
                m_nuc_V4 = TLorentzVector(V3, m_ge.Ei[i]);
                m_ps_nuc_resc = m_ge.resc[i];
            } else {
                return 0;
            }

        } else if (m_ge.pdgi[i] == 2112) { // neutron
            V3.SetXYZ(m_ge.pxi[i], m_ge.pyi[i], m_ge.pzi[i]);
            V3.SetPhi(V3.Phi() + TMath::Pi());

            // If we got here, it means it qualifies as signal, if we are still looking for a neutron, make this
            // found, if we are looking for protons or have already found a nucleon, the event is invalid so return
            // 0
            if ((m_nuc_type == NucleonType::Neutron) && (!m_found_nuc)) {
                m_found_nuc = true;
                m_nuc_V4 = TLorentzVector(V3, m_ge.Ei[i]);
                m_ps_nuc_resc = m_ge.resc[i];
            } else {
                return 0;
            }

        } else if (m_ge.pdgi[i] == 211) { // pi+
            V3.SetXYZ(m_ge.pxi[i], m_ge.pyi[i], m_ge.pzi[i]);
            V3.SetPhi(V3.Phi() + TMath::Pi());

            if (V3.Mag() < m_p_pion_momentum_threshold) {
                continue;
            }

            if (!m_fiducials->piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Plus)) {
                continue;
            }

            // If we got here, it means it qualifies as signal, if we are still looking for a pi+, make this found,
            // if we are looking for pi- or have already found a suitable pion, the event is invalid so return 0
            if ((m_pi_type == PionType::Plus) && (!m_found_pi)) {
                m_found_pi = true;
                m_pi_V4 = TLorentzVector(V3, m_ge.Ei[i]);
                m_ps_pi_resc = m_ge.resc[i];
            } else {
                return 0;
            }

        } else if (m_ge.pdgi[i] == -211) { // pi-
            V3.SetXYZ(m_ge.pxi[i], m_ge.pyi[i], m_ge.pzi[i]);
            V3.SetPhi(V3.Phi() + TMath::Pi());

            if (V3.Mag() < m_p_pion_momentum_threshold) {
                continue;
            }

            if (!m_fiducials->piAndPhotonCuts(V3, FiducialWrapper::PiPhotonId::Minus)) {
                continue;
            }

            // If we got here, it means it qualifies as signal, if we are still looking for a pi-, make this found,
            // if we are looking for pi+ or have already found a suitable pion, the event is invalid so return 0
            if ((m_pi_type == PionType::Minus) && (!m_found_pi)) {
                m_found_pi = true;
                m_pi_V4 = TLorentzVector(V3, m_ge.Ei[i]);
                m_ps_pi_resc = m_ge.resc[i];
            } else {
                return 0;
            }

        } else if (m_ge.pdgi[i] == 22) { // photon
            // Nothing now
        } else {
            /* std::cout << "extra FS particle -- " << m_ge.pdgf[i] << std::endl; */
        }
    }

    if ((!m_found_pi) || (!m_found_nuc)) {
        return 0;
    }

    return m_ge.wght;
}