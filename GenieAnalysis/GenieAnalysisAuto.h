#ifndef GENIE_ANALYSIS_AUTO_H
#define GENIE_ANALYSIS_AUTO_H

#include <array>
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

using std::vector, std::string, std::map, std::tuple, std::pair, std::function;

struct AutoProperty {
    string title;
    tuple<Int_t, Int_t, Int_t> bin_params; // How to initialize the TH1F - nbinsx, xlow, xup in order
    function<double()> get_property;
};

struct AutoType {
    string title;
    function<bool()> is_type;
};

/**
 * This implements some of the most important functionoality of the package, it takes care of the automatic histogram
 * creation and filling. There are types which designate a type of interaction like QE, RES, DIS and so on, there are
 * essentially defined by a string name and a function that returns a bool. Next there are properties and these describe
 * the interesting properties of an event like the W, or the electron momentum (implemented later on).
 *
 * For each there is a map `m_known_types` and `m_known_properties` which each map a string name of a given type or
 * property to a tuple of values that describe that type, property (different for each, see the code). These can also be
 * added in inheriting classes. Then a user can specifies the lists `m_types` and `m_properties` of an instance to
 * specify which of the known types and properties to use. Then TH1F histograms are created and filled for each
 * combination of these types and properties showing a histogram of all passing events.
 *
 * There is also support for pairing up properties against each other to form 2D histograms of TH2F types, this is done
 * through
 *
 *
 * As an additional feature mainly useful for debugging as it is not foolproof are stages. In case, for each stage
 * string passed to the constructor another set of all the TH1Fs is created and filled whenever `useEntryAtStage` is
 * called with that particular stage string passed as an argument.
 */
class GenieAnalysisAutoHistograms : public GenieAnalysis {
  private:
    const std::unique_ptr<TFile> m_output_file;

    vector<string> m_stages;

    // Which properties and types to make TH1Fs for, if they are empty (default) all known types or properties are used
    vector<string> m_properties;
    vector<string> m_types;

    // Specify what pairs of properties to make TH2Fs for, the first 2 strings specify the properties and the vector of
    // string specified for which types to make the TH2F for, if empty do for all
    vector<tuple<string, string, vector<string>>> m_vs_property_plots;

    bool m_hists_initialized{false};
    map<string, map<string, TH1F>> m_TH1Fs; // The 1D histograms for each property, type combination
    map<string, map<string, TH1F>> m_TH2Fs; // The 2D histograms showing pairs of properties vs each other
    map<string, map<string, map<string, TH1F>>> m_staged_hists; // same as m_TH1Fs but for the stages

  protected:
    map<string, AutoProperty> m_known_properties{
        {"Ws", {"Ws (generated in GENIE) from gst [GeV]", {1000, 0, 4}, [this]() { return m_ge.Ws; }}},
        {"W", {"W (computed in GENIE) from gst [GeV]", {1000, 0, 4}, [this]() { return m_ge.W; }}},
        {"fspl", {"Final state lepton PDG from gst", {51, 0, 50}, [this]() { return m_ge.fspl; }}},
        {"wght", {"wght from gst", {100, 0, 2}, [this]() { return m_ge.wght; }}}};

    map<string, AutoType> m_known_types{
        {"ALL", {"All events", [this]() { return true; }}},
        {"QE", {"Quasi-Elastic events", [this]() { return m_ge.qel; }}},
        {"RES_ALL", {"Resonant events", [this]() { return m_ge.res; }}},
        {"DELTA1232", {"Resonant events with a Delta1232", [this]() { return (m_ge.res && (m_ge.resid == 0)); }}},
        {"DIS", {"Deep-inelastic events", [this]() { return m_ge.dis; }}},
    };

  public:
    constexpr static size_t m_number_colors{7};
    std::array<Color_t, m_number_colors> m_colors{kBlack, kRed, kGreen, kBlue, kMagenta, kCyan, kOrange};

    GenieAnalysisAutoHistograms(const char *filename, const char *output_filename, const vector<string> &stages = {},
                                const vector<string> &properties = {}, const vector<string> &types = {},
                                const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          m_output_file(TFile::Open(output_filename, "RECREATE")), m_stages{stages},
          m_properties{properties}, m_types{types} {}

    void prepareAutoHists() {
        string name;
        string title;
        Int_t nbinsx, xlow, xup;

        if (m_properties.empty()) {
            for (const auto &known_property_keyval : m_known_properties) {
                m_properties.push_back(known_property_keyval.first);
            }
        }

        if (m_types.empty()) {
            for (const auto &known_type_keyval : m_known_types) {
                m_types.push_back(known_type_keyval.first);
            }
        }

        m_output_file->cd();
        // Create "final" stage (post all cuts) hists
        for (const string &property : m_properties) {
            for (const string &type : m_types) {
                name = makeTH1FName(property, type);
                title = makeTH1FTitle(property, type);
                std::tie(nbinsx, xlow, xup) = m_known_properties[property].bin_params;
                m_TH1Fs[property][type] = TH1F(name.c_str(), title.c_str(), nbinsx, xlow, xup);
            }
        }

        // Create hists for possible other stages
        for (const string &stage : m_stages) {
            for (const string &property : m_properties) {
                for (const string &type : m_types) {
                    name = makeTH1FName(stage, property, type);
                    title = makeTH1FTitle(stage, property, type);
                    std::tie(nbinsx, xlow, xup) = m_known_properties[property].bin_params;
                    m_staged_hists[stage][property][type] = TH1F(name.c_str(), title.c_str(), nbinsx, xlow, xup);
                }
            }
        }
    }

    void runPreAnalysis() override {
        if (!m_hists_initialized) {
            prepareAutoHists();
        }
    }

    void runPostAnalysis() override {
        int color_i;
        for (const string &property : m_properties) {
            THStack hist_stack{property.c_str(), m_known_properties[property].title.c_str()};
            TLegend stack_legend{33, 16, "By interaction types"};
            stack_legend.SetName((property + "_legend").c_str());

            color_i = 0;
            for (const string &type : m_types) {
                m_TH1Fs[property][type].SetLineColor(m_colors[color_i++ % m_number_colors]);
                m_TH1Fs[property][type].Write();
                hist_stack.Add(&m_TH1Fs[property][type]);
                stack_legend.AddEntry(&m_TH1Fs[property][type], m_known_types[type].title.c_str());
            }

            hist_stack.Write();
            stack_legend.Write();
        }

        for (const string &stage : m_stages) {
            for (const string &property : m_properties) {
                THStack hist_stack{(stage + "_" + property).c_str(),
                                   (m_known_properties[property].title + " - " + stage).c_str()};
                TLegend stack_legend{33, 16, "By interaction types"};
                stack_legend.SetName((stage + "_" + property + "_legend").c_str());

                color_i = 0;
                for (const string &type : m_types) {
                    m_staged_hists[stage][property][type].SetLineColor(m_colors[color_i++ % m_number_colors]);
                    m_staged_hists[stage][property][type].Write();
                    hist_stack.Add(&m_staged_hists[stage][property][type]);
                    stack_legend.AddEntry(&m_TH1Fs[property][type], m_known_types[type].title.c_str());
                }

                hist_stack.Write();
                stack_legend.Write();
            }
        }
    }

    void useEntryAtStage(const string &stage, const Double_t &weight = 1) {
        if (auto stage_hists{m_staged_hists.find(stage)}; stage_hists != m_staged_hists.end()) {
            for (const string &property : m_properties) {
                for (const string &type : m_types) {
                    if (m_known_types[type].is_type()) {
                        (stage_hists->second)[property][type].Fill(m_known_properties[property].get_property(), weight);
                    }
                }
            }
        }
    }

    void useEntry(const Double_t &weight) override {
        for (const string &property : m_properties) {
            for (const string &type : m_types) {
                if (m_known_types[type].is_type()) {
                    m_TH1Fs[property][type].Fill(m_known_properties[property].get_property(), weight);
                }
            }
        }
    }

    virtual const string makeTH1FName(const string &property, const string &type) { return property + "_" + type; }
    virtual const string makeTH1FName(const string &stage, const string &property, const string &type) {
        return stage + "_" + makeTH1FName(property, type);
    }
    virtual const string makeTH1FTitle(const string &property, const string &type) {
        return m_known_properties[property].title + " - " + m_known_types[type].title;
    }
    virtual const string makeTH1FTitle(const string &stage, const string &property, const string &type) {
        return makeTH1FTitle(property, type) + " - at " + stage;
    }
};

#endif
