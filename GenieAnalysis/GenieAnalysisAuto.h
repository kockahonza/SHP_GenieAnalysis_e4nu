#ifndef GENIE_ANALYSIS_AUTO_H
#define GENIE_ANALYSIS_AUTO_H

#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TVector3.h>

#include "GenieAnalysis.h"
#include "misc.h"

using std::vector, std::string, std::map, std::tuple, std::function;

struct AutoProperty {
    tuple<Int_t, Int_t, Int_t> bin_params; // How to initialize the TH1F - nbinsx, xlow, xup in order
    function<double()> get_property;
};

class GenieAnalysisAutoTH1Fs : public GenieAnalysis {
  private:
    const std::unique_ptr<TFile> m_output_file;

    const vector<string> m_properties;
    const vector<string> m_types;

    bool m_hists_initialized{false};
    map<string, map<string, TH1F>> m_hists;

  protected:
    map<string, AutoProperty> m_known_properties{{"W", {{1000, 0, 4}, [this]() { return m_ge.W; }}},
                                                 {"wght", {{100, 0, 2}, [this]() { return m_ge.wght; }}}};

    map<string, function<bool()>> m_known_types{
        {"ALL", [this]() { return true; }},
        {"QE", [this]() { return m_ge.qel; }},
        {"RES", [this]() { return m_ge.res; }},
        {"DIS", [this]() { return m_ge.dis; }},
    };

  public:
    GenieAnalysisAutoTH1Fs(const char *filename, const char *output_filename, const vector<string> &properties,
                           const vector<string> &types, const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          m_output_file(TFile::Open(output_filename, "RECREATE")), m_properties{properties}, m_types{types} {}

    void createTH1Fs() {
        string name_and_title;
        Int_t nbinsx, xlow, xup;

        m_output_file->cd();
        for (string property : m_properties) {
            for (string type : m_types) {
                name_and_title = makeName(property, type);
                std::tie(nbinsx, xlow, xup) = m_known_properties[property].bin_params;
                m_hists[property][type] = TH1F(name_and_title.c_str(), name_and_title.c_str(), nbinsx, xlow, xup);
            }
        }
    }

    void runPreAnalysis() override {
        if (!m_hists_initialized) {
            createTH1Fs();
        }
    }

    void runPostAnalysis() override {
        for (string property : m_properties) {
            for (string type : m_types) {
                m_hists[property][type].Write();
            }
        }
    }

    void useEntry() override {
        for (string property : m_properties) {
            for (string type : m_types) {
                if (m_known_types[type]()) {
                    m_hists[property][type].Fill(m_known_properties[property].get_property());
                }
            }
        }
    }

    virtual const string makeName(const string &property, const string &type) { return property + "_" + type; }
};

#endif
