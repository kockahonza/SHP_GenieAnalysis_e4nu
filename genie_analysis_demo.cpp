#include <iostream>

#include <TH1F.h>

#include "GenieAnalysis/GenieAnalysis.h"
#include "GenieAnalysis/GAAutoHistograms.h"
#include "GenieAnalysis/misc.h"

constexpr Int_t W_hist_bins{1000};
constexpr Double_t W_hist_min{0};
constexpr Double_t W_hist_max{4};

class GenieAnalysisDemoW : public GenieAnalysis {
  private:
    const std::unique_ptr<TFile> m_output_file;

  public:
    TH1F m_hist_W_qe;
    TH1F m_hist_W_res;
    TH1F m_hist_W_dis;
    TH1F m_hist_W_total;

  public:
    GenieAnalysisDemoW(const char *filename, const char *output_filename, const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name), m_output_file(TFile::Open(output_filename, "RECREATE")),
          m_hist_W_qe("W_qe", "W of Quasielastic events", W_hist_bins, W_hist_min, W_hist_max),
          m_hist_W_res("W_res", "W of Resonant events", W_hist_bins, W_hist_min, W_hist_max),
          m_hist_W_dis("W_dis", "W of Deep Inelastic events", W_hist_bins, W_hist_min, W_hist_max),
          m_hist_W_total("W_total", "W of all events", W_hist_bins, W_hist_min, W_hist_max) {}

    void useEntry(const Double_t &weight) override {
        if (m_ge.qel) {
            m_hist_W_qe.Fill(m_ge.W);
        }
        if (m_ge.res) {
            m_hist_W_res.Fill(m_ge.W);
        }
        if (m_ge.dis) {
            m_hist_W_dis.Fill(m_ge.W);
        }
        m_hist_W_total.Fill(m_ge.W);
    }

    ~GenieAnalysisDemoW() {
        std::cout << "Destroying DemoW" << std::endl;
        m_hist_W_qe.Write();
        m_hist_W_res.Write();
        m_hist_W_dis.Write();
        m_hist_W_total.Write();
    }
};

int main(int argc, char *argv[]) {

    GenieAnalysisDemoW gademo{"/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root",

                              "genie_analysis_demo_output.root"};

    gademo.runAnalysis();

    return 0;
}
