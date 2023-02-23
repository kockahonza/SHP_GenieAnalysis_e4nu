#ifndef GENIE_ANALYSIS_AUTO_H
#define GENIE_ANALYSIS_AUTO_H

#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TTree.h>
#include <TVector3.h>

#include "GenieAnalysis.h"
#include "misc.h"

using std::vector, std::string, std::map, std::tuple, std::function;

struct AutoProperty {
    string title;
    tuple<Int_t, Int_t, Int_t> bin_params; // How to initialize the TH1F - nbinsx, xlow, xup in order
    function<double()> get_property;
};

struct AutoType {
    string title;
    function<bool()> is_type;
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
    map<string, AutoProperty> m_known_properties{
        {"W", {"W from gst", {1000, 0, 4}, [this]() { return m_ge.W; }}},
        {"wght", {"wght from gst", {100, 0, 2}, [this]() { return m_ge.wght; }}}};

    map<string, AutoType> m_known_types{
        {"ALL", {"All events", [this]() { return true; }}},
        {"QE", {"Quasi-Elastic events", [this]() { return m_ge.qel; }}},
        {"RES_ALL", {"Resonant events", [this]() { return m_ge.res; }}},
        {"DELTA1232", {"Resonant events with a Delta1232", [this]() { return (m_ge.res && (m_ge.resid == 0)); }}},
        {"DIS", {"Deep-inelastic events", [this]() { return m_ge.dis; }}},
    };

  public:
    GenieAnalysisAutoTH1Fs(const char *filename, const char *output_filename, const vector<string> &stages,
                           const vector<string> &properties, const vector<string> &types,
                           const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          m_output_file(TFile::Open(output_filename, "RECREATE")), m_stages{stages},
          m_properties{properties}, m_types{types} {}

    void createTH1Fs() {
        string name;
        string title;
        Int_t nbinsx, xlow, xup;

        m_output_file->cd();
        // Create "final" stage (post all cuts) hists
        for (string property : m_properties) {
            for (string type : m_types) {
                name = makeName(property, type);
                title = makeTitle(property, type);
                std::tie(nbinsx, xlow, xup) = m_known_properties[property].bin_params;
                m_hists[property][type] = TH1F(name.c_str(), title.c_str(), nbinsx, xlow, xup);
            }
        }

        // Create hists for possible other stages
        for (string stage : m_stages) {
            for (string property : m_properties) {
                for (string type : m_types) {
                    name = makeName(stage, property, type);
                    title = makeTitle(stage, property, type);
                    std::tie(nbinsx, xlow, xup) = m_known_properties[property].bin_params;
                    m_staged_hists[stage][property][type] = TH1F(name.c_str(), title.c_str(), nbinsx, xlow, xup);
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
            THStack hist_stack{property.c_str(), m_known_properties[property].title.c_str()};
            TLegend stack_legend{33, 16, "By interaction types"};
            stack_legend.SetName((property + "_legend").c_str());

            Color_t color{0};
            for (string type : m_types) {
                // Using colors 1 to 49 right now - this is a specific part of the default ROOT color table
                m_hists[property][type].SetLineColor(1 + color);
                color = (color + 1) % 10;
                m_hists[property][type].Write();
                hist_stack.Add(&m_hists[property][type]);
                stack_legend.AddEntry(&m_hists[property][type], m_known_types[type].title.c_str());
            }

            hist_stack.Write();
            stack_legend.Write();
        }

        for (string stage : m_stages) {
            for (string property : m_properties) {
                THStack hist_stack{(stage + "_" + property).c_str(),
                                   (m_known_properties[property].title + " - at " + stage).c_str()};
                TLegend stack_legend{33, 16, "By interaction types"};
                stack_legend.SetName((stage + "_" + property + "_legend").c_str());

                Color_t color{0};
                for (string type : m_types) {
                    m_staged_hists[stage][property][type].SetLineColor(1 + color);
                    color = (color + 1) % 10;
                    m_staged_hists[stage][property][type].Write();
                    hist_stack.Add(&m_staged_hists[stage][property][type]);
                    stack_legend.AddEntry(&m_hists[property][type], m_known_types[type].title.c_str());
                }

                hist_stack.Write();
                stack_legend.Write();
            }
        }
    }

    void useEntryAtStage(string stage, Double_t weight = 1) {
        for (string property : m_properties) {
            for (string type : m_types) {
                if (m_known_types[type].is_type()) {
                    m_staged_hists[stage][property][type].Fill(m_known_properties[property].get_property(), weight);
                }
            }
        }
    }

    void useEntry(Double_t weight) override {
        for (string property : m_properties) {
            for (string type : m_types) {
                if (m_known_types[type].is_type()) {
                    m_hists[property][type].Fill(m_known_properties[property].get_property(), weight);
                }
            }
        }
    }

    virtual const string makeName(const string &property, const string &type) { return property + "_" + type; }
    virtual const string makeName(const string &stage, const string &property, const string &type) {
        return stage + "_" + makeName(property, type);
    }
    virtual const string makeTitle(const string &property, const string &type) {
        return m_known_properties[property].title + " - " + m_known_types[type].title;
    }
    virtual const string makeTitle(const string &stage, const string &property, const string &type) {
        return makeTitle(property, type) + " - at " + stage;
    }
};

#endif
