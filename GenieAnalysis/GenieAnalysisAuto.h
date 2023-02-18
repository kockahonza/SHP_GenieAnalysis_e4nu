#ifndef GENIE_ANALYSIS_AUTO_H
#define GENIE_ANALYSIS_AUTO_H

#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
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

    const vector<string> m_stages;
    const vector<string> m_properties;
    const vector<string> m_types;

    bool m_hists_initialized{false};
    map<string, map<string, TH1F>> m_hists;
    map<string, map<string, map<string, TH1F>>> m_staged_hists;

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
    GenieAnalysisAutoTH1Fs(const char *filename, const char *output_filename, const vector<string> &stages,
                           const vector<string> &properties, const vector<string> &types,
                           const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          m_output_file(TFile::Open(output_filename, "RECREATE")), m_stages{stages},
          m_properties{properties}, m_types{types} {}

    void createTH1Fs() {
        string name_and_title;
        Int_t nbinsx, xlow, xup;

        m_output_file->cd();
        // Create "final" stage (post all cuts) hists
        for (string property : m_properties) {
            for (string type : m_types) {
                name_and_title = makeName(property, type);
                std::tie(nbinsx, xlow, xup) = m_known_properties[property].bin_params;
                m_hists[property][type] = TH1F(name_and_title.c_str(), name_and_title.c_str(), nbinsx, xlow, xup);
            }
        }

        // Create hists for possible other stages
        for (string stage : m_stages) {
            for (string property : m_properties) {
                for (string type : m_types) {
                    name_and_title = makeName(stage, property, type);
                    std::tie(nbinsx, xlow, xup) = m_known_properties[property].bin_params;
                    m_staged_hists[stage][property][type] =
                        TH1F(name_and_title.c_str(), name_and_title.c_str(), nbinsx, xlow, xup);
                }
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
            THStack th_p{property.c_str(), (property + " stack").c_str()};

            Color_t color{0};
            for (string type : m_types) {
                // Using colors 1 to 49 right now - this is a specific part of the default ROOT color table
                m_hists[property][type].SetLineColor(1 + color);
                color = (color + 1) % 10;
                m_hists[property][type].Write();
                th_p.Add(&m_hists[property][type]);
            }

            th_p.Write();
        }

        for (string stage : m_stages) {
            for (string property : m_properties) {
                THStack th_p{(stage + "_" + property).c_str(), (stage + "_" + property + " stack").c_str()};

                Color_t color{0};
                for (string type : m_types) {
                    m_staged_hists[stage][property][type].SetLineColor(1 + color);
                    color = (color + 1) % 10;
                    m_staged_hists[stage][property][type].Write();
                    th_p.Add(&m_staged_hists[stage][property][type]);
                }

                th_p.Write();
            }
        }
    }

    void useEntryAtStage(string stage) {
        for (string property : m_properties) {
            for (string type : m_types) {
                if (m_known_types[type]()) {
                    m_staged_hists[stage][property][type].Fill(m_known_properties[property].get_property());
                }
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
    virtual const string makeName(const string &stage, const string &property, const string &type) {
        return stage + "_" + property + "_" + type;
    }
};

#endif
