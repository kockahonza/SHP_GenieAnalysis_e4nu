#ifndef GACLAS6COMMON_H
#define GACLAS6COMMON_H

#include "Fiducial.h"

using std::unique_ptr, std::optional, std::string;

namespace GACLAS6Common {
// Configuration options for the major target/energy runs
enum class Target { C12, Fe56 };
// If the ~1GeV energy is to be added, smearing resolutions were tripled!, do't forget to add that
enum class BeamEnergy { MeV_1161, MeV_2261, MeV_4461 };

/**
 * A very simple wrapper around the Fiducial class, which is taken exactly as in original, so that it can be referred
 * to in a more modern way. No need to pas strings to specify what is to be done in each call.
 */
class FiducialWrapper {
  public:
    enum class PiPhotonId : int { Minus = -1, Photon = 0, Plus = 1 };

    const GACLAS6Common::Target m_target;
    const GACLAS6Common::BeamEnergy m_beam_energy;
    const string m_target_str;
    const string m_beam_energy_str;

  private:
    Fiducial m_fiducial;

  public:
    FiducialWrapper(const GACLAS6Common::Target &target, const GACLAS6Common::BeamEnergy &beam_energy)
        : m_target(target), m_beam_energy(beam_energy), m_target_str{targetStr(m_target)},
          m_beam_energy_str{beamEnergyStr(m_beam_energy)}, m_fiducial{} {

        m_fiducial.InitPiMinusFit(m_beam_energy_str);
        m_fiducial.InitEClimits();

        // The first value is torusCurrent, values taken from original
        m_fiducial.SetConstants(m_beam_energy == BeamEnergy::MeV_1161 ? 750 : 2250, m_target_str,
                                {{"1161", 1.161}, {"2261", 2.261}, {"4461", 4.461}});
        m_fiducial.SetFiducialCutParameters(m_beam_energy_str);
    }

    bool electronCut(const TVector3 &momentum_V3) { return m_fiducial.EFiducialCut(m_beam_energy_str, momentum_V3); }

    bool protonCut(const TVector3 &momentum_V3) { return m_fiducial.PFiducialCut(m_beam_energy_str, momentum_V3); }

    bool piAndPhotonCuts(const TVector3 &momentum_V3, const PiPhotonId &which_particle) {
        return m_fiducial.Pi_phot_fid_united(m_beam_energy_str, momentum_V3, static_cast<int>(which_particle));
    }

    static string targetStr(const Target &target) {
        if (target == Target::C12) {
            return "12C";
        } else if (target == Target::Fe56) {
            return "12C"; // There's no dedicated Fe file and original used 12C for anything except He
        }
        throw "This should not happen";
    }

    static string beamEnergyStr(const BeamEnergy &beam_energy) {
        if (beam_energy == BeamEnergy::MeV_1161) {
            return "1161";
        } else if (beam_energy == BeamEnergy::MeV_2261) {
            return "2261";
        } else if (beam_energy == BeamEnergy::MeV_4461) {
            return "4461";
        }
        throw "This should not happen";
    }
};
} // namespace GACLAS6Common

#endif
