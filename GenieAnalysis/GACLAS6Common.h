#ifndef GACLAS6COMMON_H
#define GACLAS6COMMON_H

#include <optional>
#include <stdexcept>

#include "Fiducial.h"

using std::unique_ptr, std::optional, std::string;

namespace GACLAS6Common {
// Configuration options for the major target/energy runs
enum class Target { C12, Fe56 };
// If the ~1GeV energy is to be added, smearing resolutions were tripled!, do't forget to add that
enum class BeamEnergy { MeV_1161, MeV_2261, MeV_4461 };

string beam_energy_str_fiducials(const BeamEnergy &beam_energy);

/**
 * A very simple wrapper around the Fiducial class, which is taken exactly as in original, so that it can be referred
 * to in a more modern way. No need to pas strings to specify what is to be done in each call.
 */
class FiducialWrapper {
  public:
    enum class PiPhotonId : int { Minus = -1, Photon = 0, Plus = 1 };

    const GACLAS6Common::BeamEnergy m_beam_energy;
    const string m_beam_energy_str;

  private:
    Fiducial m_fiducial;

  public:
    FiducialWrapper(const GACLAS6Common::BeamEnergy &beam_energy)
        : m_beam_energy(beam_energy), m_beam_energy_str{beam_energy_str_fiducials(m_beam_energy)}, m_fiducial{} {

        m_fiducial.InitPiMinusFit(m_beam_energy_str);
        m_fiducial.InitEClimits();

        // The first value is torusCurrent, values taken from original
        m_fiducial.SetConstants(m_beam_energy == BeamEnergy::MeV_1161 ? 750 : 2250, "",
                                {{"1161", 1.161}, {"2261", 2.261}, {"4461", 4.461}});
        m_fiducial.SetFiducialCutParameters(m_beam_energy_str);
    }

    bool electronCut(const TVector3 &momentum_V3) { return m_fiducial.EFiducialCut(m_beam_energy_str, momentum_V3); }

    bool protonCut(const TVector3 &momentum_V3) { return m_fiducial.PFiducialCut(m_beam_energy_str, momentum_V3); }

    bool piAndPhotonCuts(const TVector3 &momentum_V3, const PiPhotonId &which_particle) {
        return m_fiducial.Pi_phot_fid_united(m_beam_energy_str, momentum_V3, static_cast<int>(which_particle));
    }
};

string target_str_filename(const Target &target);
string get_data_filename(GACLAS6Common::Target t, GACLAS6Common::BeamEnergy be, string arg);

} // namespace GACLAS6Common

#endif
