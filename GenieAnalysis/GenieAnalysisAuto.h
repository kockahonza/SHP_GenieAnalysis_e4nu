#ifndef GENIE_ANALYSIS_AUTO_H
#define GENIE_ANALYSIS_AUTO_H

#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

#include "GenieAnalysis.h"

using std::vector, std::string, std::map, std::tuple;

class GenieAnalysisAutoTH1Fs : public GenieAnalysis {
  private:
    const std::unique_ptr<TFile> m_output_file;

    const vector<string> m_properties;
    const vector<string> m_observables;

    bool m_hists_initialized{false};
    map<string, map<string, TH1F>> m_hists;

  public:
    // Specifies how to initialize the histograms for each property
    map<string, tuple<Int_t, Int_t, Int_t>> m_bin_params{{"W", {1000, 0, 4}}};

  public:
    GenieAnalysisAutoTH1Fs(const char *filename, const char *output_filename, const vector<string> &properties,
                           const vector<string> &observables, const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          m_output_file(TFile::Open(output_filename, "RECREATE")), m_properties{properties}, m_observables{
                                                                                                 observables} {
        m_output_file->cd();
    }

    void createTH1Fs() {
        string name_and_title;
        Int_t nbinsx, xlow, xup;

        for (string property : m_properties) {
            for (string observable : m_observables) {
                name_and_title = makeName(property, observable);
                std::tie(nbinsx, xlow, xup) = m_bin_params[property];
                m_hists[property][observable] = TH1F(name_and_title.c_str(), name_and_title.c_str(), nbinsx, xlow, xup);
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
            for (string observable : m_observables) {
                m_hists[property][observable].Write();
            }
        }
    }

    void useEntry() override {
        for (string property : m_properties) {
            for (string observable : m_observables) {
                if (isObservable(observable, m_loaded_event)) {
                    m_hists[property][observable].Fill(getProperty(property, m_loaded_event));
                }
            }
        }
    }

    virtual const string makeName(const string &property, const string &observable) {
        return property + "_" + observable;
    }

    static bool isObservable(const string &observable, const GenieEvent &ge);
    static Double_t getProperty(const string &property, const GenieEvent &ge);
};

#endif
