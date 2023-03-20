#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <tuple>

#include "GenieAnalysis/GACLAS6MC.h"
#include "GenieAnalysis/misc.h"

class GAFinalStateTransparency : public GACLAS6MC {
  protected:
    bool m_resc_all_same;

    bool m_resc_all_mone;
    bool m_resc_all_zero;
    bool m_resc_all_one;

    map<string, AutoProperty> m_new_known_properties{
        {"resc_same", {"Are all resc values the same", {2, -0.5, 1.5}, [this]() { return m_resc_all_same; }}},
        {"resc_all_mone", {"Are all resc values -1", {2, -0.5, 1.5}, [this]() { return m_resc_all_mone; }}},
        {"resc_all_zero", {"Are all resc values 0", {2, -0.5, 1.5}, [this]() { return m_resc_all_zero; }}},
        {"resc_all_one", {"Are all resc values 1", {2, -0.5, 1.5}, [this]() { return m_resc_all_one; }}},
    };

  public:
    GAFinalStateTransparency(const char *filename, const char *output_filename, const vector<string> &stages = {},
                             const vector<string> &properties = {}, const vector<string> &types = {},
                             const vector<GAAutoHistograms::AutoVsPlot> &vs_property_plots = {},
                             const Target &target = Target::C12, const BeamEnergy &beam_energy = BeamEnergy::MeV_2261,
                             const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), GACLAS6MC{filename, output_filename,   stages, properties,
                                                             types,    vs_property_plots, target, beam_energy} {
        m_known_properties.insert(m_new_known_properties.begin(), m_new_known_properties.end());
    }

    Double_t passesCuts() override {
        Double_t weight{GACLAS6MC::passesCuts()};

        if (weight == 0) {
            return 0;
        }

        m_resc_all_same = true;
        for (int i{1}; i < m_ge.ni; i++) {
            if (m_ge.resc[i] != m_ge.resc[0]) {
                m_resc_all_same = false;
            }
        }

        m_resc_all_mone = false;
        m_resc_all_zero = false;
        m_resc_all_one = false;

        if (m_resc_all_same) {
            if (m_ge.resc[0] == -1) {
                m_resc_all_mone = true;
            } else if (m_ge.resc[0] == 0) {
                m_resc_all_zero = true;
            } else if (m_ge.resc[0] == 1) {
                m_resc_all_one = true;
            }
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
        return -1;
    }

    /* GACLAS6MC ga1{input_file.c_str(), */
    /*                               (output_file + "_FSt_new.root").c_str(), */
    /*                               {}, */
    /*                               {}, */
    /*                               {}, */
    /*                               { */
    /*                                   {"W", "Ws", {}}, */
    /*                                   {"W", "el_p", {}}, */
    /*                                   {"W", "bjorken_x", {}}, */
    /*                                   {"Ws", "bjorken_x", {}}, */
    /*                                   {"ps_num_pim", "fs_num_pim", {}}, */
    /*                                   {"ps_num_pip", "fs_num_pip", {}}, */
    /*                                   {"fs_num_pim", "passed_num_pim", {}}, */
    /*                                   {"fs_num_pip", "passed_num_pip", {}}, */
    /*                                   {"passed_num_pip", "passed_num_pim", {}}, */
    /*                                   /1* {"fs_num_pim", "passed_num_pim", {}}, *1/ */
    /*                                   /1* {"fs_num_pip", "passed_num_pip", {}}, *1/ */
    /*                                   /1* {"fs_num_pip", "fs_num_neutrons", {}}, *1/ */
    /*                                   /1* {"fs_num_pim", "fs_num_neutrons", {}}, *1/ */
    /*                                   /1* {"fs_num_pip", "fs_num_protons", {}}, *1/ */
    /*                                   /1* {"fs_num_pim", "fs_num_protons", {}}, *1/ */
    /*                               }}; */
    /* ga1.runAnalysis(); */

    GAFinalStateTransparency ga1{input_file.c_str(),
                                 (output_file + "_FSt_resc.root").c_str(),
                                 {},
                                 {},
                                 {},
                                 {{"W", "Ws", {}},
                                  {"fs_num_pip", "nfpip", {}},
                                  {"fs_num_protons", "nfp", {}},
                                  {"resc", "passed_num_pip", {}},
                                  {"resc", "passed_num_pim", {}},
                                  {"resc", "passed_num_protons", {}},
                                  {"resc", "fs_num_pip", {}},
                                  {"resc", "fs_num_pim", {}},
                                  {"resc", "fs_num_protons", {}}}};
    ga1.runAnalysis();

    return 0;
}
