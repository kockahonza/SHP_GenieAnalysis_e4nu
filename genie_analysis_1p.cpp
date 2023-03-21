#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <tuple>

#include "GenieAnalysis/GACLAS6MC.h"
#include "GenieAnalysis/misc.h"

class GA1Proton : public GACLAS6MC {
  public:
    GA1Proton(const char *filename, const char *output_filename,
              // Specify the analysis - which stages, properties and types to do histograms for
              const vector<string> &stages = {}, const vector<string> &properties = {},
              const vector<string> &types = {}, const vector<GAAutoHistograms::AutoVsPlot> &vs_property_plots = {},
              // Select run
              const GACLAS6Common::Target &target = GACLAS6Common::Target::C12,
              const GACLAS6Common::BeamEnergy &beam_energy = GACLAS6Common::BeamEnergy::MeV_2261,
              // Pass this to GenieAnalysis
              const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          GACLAS6MC(filename, output_filename, stages, properties, types, vs_property_plots, target, beam_energy) {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }

    Double_t passesCuts() override {
        Double_t weight{GACLAS6MC::passesCuts()};

        if (m_fs_number_of_protons != 1) {
            return 0;
        }

        return weight;
    }
};

int main(int argc, char *argv[]) {
    string input_file, output_file;
    if (argc == 2) {
        std::string arg{argv[1]};
        if (arg == "local") {
            input_file = "/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root";
            output_file = "output_local";
        } else if (arg == "full") {
            input_file = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/"
                         "apapadop_SuSav2_C12_2261GeV_master.root";
            output_file = "output_full";
        }
    } else {
        std::cout << "Needs an argument to specify which way to run (should be \"local\" or \"full\"), check the code"
                  << std::endl;
    }

    GA1Proton ga1{input_file.c_str(),
                  (output_file + "_1p.root").c_str(),
                  {},
                  {
                      "W",
                      "Ws",
                      "resc",
                      "reco_W",
                      "bjorken_x",
                      "el_p",
                      "passed_num_pip",
                      "passed_num_pim",
                      "ps_num_pip",
                      "ps_num_pim",
                      "ps_num_protons",
                      "ps_num_neutrons",
                      "fs_num_pip",
                      "fs_num_pim",
                      "fs_num_protons",
                      "fs_num_neutrons",
                  },
                  {},
                  {
                      {"W", "Ws", {}},
                      {"W", "el_p", {}},
                      {"W", "bjorken_x", {}},
                      {"Ws", "bjorken_x", {}},
                      {"passed_num_pip", "passed_num_pim", {}},
                      {"fs_num_pim", "passed_num_pim", {}},
                      {"fs_num_pip", "passed_num_pip", {}},
                  }};
    ga1.runAnalysis();

    return 0;
}
