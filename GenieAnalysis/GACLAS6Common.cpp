#include "GACLAS6Common.h"

using namespace GACLAS6Common;

string GACLAS6Common::targetStr(const Target &target) {
    if (target == Target::C12) {
        return "12C";
    } else if (target == Target::Fe56) {
        return "12C"; // There's no dedicated Fe file and original used 12C for anything except He
    }
    throw std::runtime_error("This should not happen");
}

string GACLAS6Common::beamEnergyStr(const BeamEnergy &beam_energy) {
    if (beam_energy == BeamEnergy::MeV_1161) {
        return "1161";
    } else if (beam_energy == BeamEnergy::MeV_2261) {
        return "2261";
    } else if (beam_energy == BeamEnergy::MeV_4461) {
        return "4461";
    }
    throw std::runtime_error("This should not happen");
}

string GACLAS6Common::get_data_filename(GACLAS6Common::Target t, GACLAS6Common::BeamEnergy be, string arg) {
    // this is specific to my (Jan's) computer
    if (arg == "local") {
        if ((t != Target::C12) || (be != BeamEnergy::MeV_2261)) {
            throw std::runtime_error("Only C12 ~2GeV is available locally (this is specific to my (Jan's) computer)");
        } else {
            return "/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root";
        }
    } else if (arg == "full") {
        return "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/" + targetStr(t) + "_" +
               beamEnergyStr(be) + "GeV/apapadop_SuSav2_C12_2261GeV_master.root";
    } else {
        throw std::runtime_error("Needs an argument to specify which way to run (should be \"local\" or \"full\")");
    }
}
