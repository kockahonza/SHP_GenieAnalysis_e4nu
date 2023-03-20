#ifndef GENIE_ANALYSIS_AUTO_H
#define GENIE_ANALYSIS_AUTO_H

#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TTree.h>
#include <TVector3.h>

#include "GenieAnalysis.h"
#include "misc.h"

using std::vector, std::string, std::map, std::tuple, std::pair, std::function;

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
class GAAutoHistograms : public virtual GenieAnalysis {
  public:
    struct AutoProperty {
        string title;
        tuple<Int_t, Double_t, Double_t> bin_params; // How to initialize the TH1F - nbinsx, xlow, xup in order
        function<double()> get_property;
    };

    struct AutoType {
        string title;
        function<bool()> is_type;
    };

    struct AutoVsPlot {
        string property1;
        string property2;
        vector<string> types;
    };

  private:
    const std::unique_ptr<TFile> m_output_file;

    vector<string> m_stages;

    // Which properties and types to make simple property histograms for, if they are empty (default) all known types or
    // properties are used
    vector<string> m_properties;
    vector<string> m_types;

    // Specify what pairs of properties to make TH2Fs for, the first 2 strings specify the properties and the vector of
    // string specified for which types to make the TH2F for, if empty do for all
    vector<AutoVsPlot> m_vs_property_plots;

    bool m_hists_initialized{false};
    map<string, map<string, TH1F>> m_simple_property_hists; // The TH1F objects, map keys are property, type in order
    map<string, map<string, map<string, TH1F>>> m_staged_hists; // same as m_simple_property_hists but for the stages,
                                                                // map keys are stage, property and type in order
    map<string, map<string, map<string, TH2F>>> m_vs_property_hists; // The TH2F object to vs property plots, the map
                                                                     // keys are property1, property2 and type in order
    static constexpr Double_t m_type_legend_w{28};
    static constexpr Double_t m_type_legend_h{16};
    static constexpr const char *m_type_legend_label{""};

  protected:
    map<string, AutoProperty> m_known_properties{
        {"Ws", {"Ws (generated in GENIE) from gst [GeV]", {1000, 0, 4}, [this]() { return m_ge.Ws; }}},
        {"W", {"W (computed in GENIE) from gst [GeV]", {1000, 0, 4}, [this]() { return m_ge.W; }}},
        {"fspl", {"Final state lepton PDG from gst", {51, -0.5, 50.5}, [this]() { return m_ge.fspl; }}},
        {"resc", {"resc code from gst", {301, -150.5, 150.5}, [this]() { return m_ge.resc[0]; }}},
        {"nfpip", {"nfpip from gst (number of final state pi plus)", {6, -0.5, 5.5}, [this]() { return m_ge.nfpip; }}},
        {"nfpim", {"nfpim from gst (number of final state pi minus)", {6, -0.5, 5.5}, [this]() { return m_ge.nfpim; }}},
        {"nfp", {"nfp from gst (number of final state (anti)protons)", {6, -0.5, 5.5}, [this]() { return m_ge.nfp; }}},
        {"nfn", {"nfn from gst (number of final state (anti)neutrons)", {6, -0.5, 5.5}, [this]() { return m_ge.nfn; }}},
        {"wght", {"wght from gst", {100, -2, 2}, [this]() { return m_ge.wght; }}}};

    map<string, AutoType> m_known_types{
        {"ALL", {"All", [this]() { return true; }}},
        {"QE", {"Quasi-Elastic", [this]() { return m_ge.qel; }}},
        {"MEC", {"Meson-exchange", [this]() { return m_ge.mec; }}},
        {"RES_ALL", {"Resonant", [this]() { return m_ge.res; }}},
        {"DELTA1232", {"Delta1232", [this]() { return (m_ge.res && (m_ge.resid == 0)); }}},
        {"DIS", {"Deep-inelastic", [this]() { return m_ge.dis; }}},
    };

  public:
    constexpr static size_t m_number_colors{7};
    std::array<Color_t, m_number_colors> m_colors{kBlack, kRed, kGreen, kBlue, kMagenta, kCyan, kOrange};

    GAAutoHistograms(const char *filename, const char *output_filename, const vector<string> &stages = {},
                                const vector<string> &properties = {}, const vector<string> &types = {},
                                const vector<AutoVsPlot> vs_property_plots = {}, const char *gst_ttree_name = "gst")
        : GenieAnalysis(filename, gst_ttree_name),
          m_output_file(TFile::Open(output_filename, "RECREATE")), m_stages{stages},
          m_properties{properties}, m_types{types}, m_vs_property_plots{vs_property_plots} {}

    void prepareAutoHists();

    void runPreAnalysis() override {
        if (!m_hists_initialized) {
            prepareAutoHists();
        }
    }

    void runPostAnalysis(const Long64_t &) override;

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
                    m_simple_property_hists[property][type].Fill(m_known_properties[property].get_property(), weight);
                }
            }
        }

        for (auto const &[property1, property2, types] : m_vs_property_plots) {
            for (auto const &type : types) {
                if (m_known_types[type].is_type()) {
                    m_vs_property_hists[property1][property2][type].Fill(m_known_properties[property1].get_property(),
                                                                         m_known_properties[property2].get_property(),
                                                                         weight);
                }
            }
        }
    }

    virtual const string makeSimplePropertyHistName(const string &property, const string &type) {
        return property + "_" + type;
    }
    virtual const string makeSimplePropertyHistName(const string &stage, const string &property, const string &type) {
        return stage + "_" + makeSimplePropertyHistName(property, type);
    }
    virtual const string makeSimplePropertyHistTitle(const string &property, const string &type) {
        return m_known_properties[property].title + " - " + m_known_types[type].title;
    }
    virtual const string makeSimplePropertyHistTitle(const string &stage, const string &property, const string &type) {
        return makeSimplePropertyHistTitle(property, type) + " - at " + stage;
    }

    virtual const string makeVsPlotName(const string &property1, const string &property2) {
        return "vs_" + property1 + "_" + property2;
    }
    virtual const string makeVsPlotTitle(const string &property1, const string &property2) {
        return m_known_properties[property1].title + " vs " + m_known_properties[property2].title;
    }
    virtual const string makeVsPlotName(const string &property1, const string &property2, const string &type) {
        return makeVsPlotName(property1, property2) + "_" + type;
    }
    virtual const string makeVsPlotTitle(const string &property1, const string &property2, const string &type) {
        return makeVsPlotTitle(property1, property2) + " - " + m_known_types[type].title;
    }
};

#endif
