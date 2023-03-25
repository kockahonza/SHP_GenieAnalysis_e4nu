#include "GACLAS6Common.h"

using namespace GACLAS6Common;

string GACLAS6Common::beam_energy_str_fiducials(const BeamEnergy &beam_energy) {
    if (beam_energy == BeamEnergy::MeV_1161) {
        return "1161";
    } else if (beam_energy == BeamEnergy::MeV_2261) {
        return "2261";
    } else if (beam_energy == BeamEnergy::MeV_4461) {
        return "4461";
    }
    throw std::runtime_error("This should not happen");
}

string GACLAS6Common::target_str_filename(const Target &target) {
    if (target == Target::C12) {
        return "C12";
    } else if (target == Target::Fe56) {
        return "56Fe";
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
        return "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/" + target_str_filename(t) + "_" +
               beam_energy_str_fiducials(be) + "GeV/apapadop_SuSav2_" + target_str_filename(t) + "_" +
               beam_energy_str_fiducials(be) + "GeV_master.root";
    } else {
        throw std::runtime_error("Needs an argument to specify which way to run (should be \"local\" or \"full\")");
    }
}
