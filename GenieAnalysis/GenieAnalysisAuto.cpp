#include "GenieAnalysisAuto.h"

void GenieAnalysisAutoHistograms::prepareAutoHists() {
    string name;
    string title;
    Int_t nbinsx;
    Double_t xlow, xup;

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
            name = makeSimplePropertyHistName(property, type);
            title = makeSimplePropertyHistTitle(property, type);
            std::tie(nbinsx, xlow, xup) = m_known_properties[property].bin_params;
            m_simple_property_hists[property][type] = TH1F(name.c_str(), title.c_str(), nbinsx, xlow, xup);
        }
    }

    // Create hists for possible other stages
    for (const string &stage : m_stages) {
        for (const string &property : m_properties) {
            for (const string &type : m_types) {
                name = makeSimplePropertyHistName(stage, property, type);
                title = makeSimplePropertyHistTitle(stage, property, type);
                std::tie(nbinsx, xlow, xup) = m_known_properties[property].bin_params;
                m_staged_hists[stage][property][type] = TH1F(name.c_str(), title.c_str(), nbinsx, xlow, xup);
            }
        }
    }

    Int_t nbinsy;
    Double_t ylow, yup;
    for (auto &[property1, property2, types] : m_vs_property_plots) {
        if (types.empty()) {
            for (const auto &known_type_keyval : m_known_types) {
                types.push_back(known_type_keyval.first);
            }
        }
        for (auto const &type : types) {
            name = makeVsPlotName(property1, property2, type);
            std::tie(nbinsx, xlow, xup) = m_known_properties[property1].bin_params;
            std::tie(nbinsy, ylow, yup) = m_known_properties[property2].bin_params;
            m_vs_property_hists[property1][property2][type] =
                TH2F(name.c_str(), m_known_types[type].title.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup);
        }
    }
}

void GenieAnalysisAutoHistograms::runPostAnalysis(const Long64_t &number_of_entries) {
    int color_i;
    for (const string &property : m_properties) {
        THStack hist_stack{property.c_str(), m_known_properties[property].title.c_str()};
        TLegend stack_legend{m_type_legend_w, m_type_legend_h, m_type_legend_label};
        stack_legend.SetName((property + "_legend").c_str());

        color_i = 0;
        for (const string &type : m_types) {
            m_simple_property_hists[property][type].SetLineColor(m_colors[color_i++ % m_number_colors]);
            m_simple_property_hists[property][type].Scale(1.0 / number_of_entries);
            m_simple_property_hists[property][type].Write();

            m_simple_property_hists[property][type].SetTitle(m_known_types[type].title.c_str());
            hist_stack.Add(&m_simple_property_hists[property][type]);
            stack_legend.AddEntry(&m_simple_property_hists[property][type], m_known_types[type].title.c_str());
        }

        hist_stack.Write();
        stack_legend.Write();
    }

    for (const string &stage : m_stages) {
        for (const string &property : m_properties) {
            THStack hist_stack{(stage + "_" + property).c_str(),
                               (m_known_properties[property].title + " - " + stage).c_str()};
            TLegend stack_legend{m_type_legend_w, m_type_legend_h, m_type_legend_label};
            stack_legend.SetName((stage + "_" + property + "_legend").c_str());

            color_i = 0;
            for (const string &type : m_types) {
                m_staged_hists[stage][property][type].SetLineColor(m_colors[color_i++ % m_number_colors]);
                m_staged_hists[stage][property][type].Scale(1.0 / number_of_entries);
                m_staged_hists[stage][property][type].Write();

                m_staged_hists[stage][property][type].SetTitle(m_known_types[type].title.c_str());
                hist_stack.Add(&m_staged_hists[stage][property][type]);
                stack_legend.AddEntry(&m_simple_property_hists[property][type], m_known_types[type].title.c_str());
            }

            hist_stack.Write();
            stack_legend.Write();
        }
    }

    for (auto const &[property1, property2, types] : m_vs_property_plots) {
        THStack hist_stack{makeVsPlotName(property1, property2).c_str(), makeVsPlotTitle(property1, property2).c_str()};
        TLegend stack_legend{m_type_legend_w, m_type_legend_h, m_type_legend_label};
        stack_legend.SetName((makeVsPlotName(property1, property2) + "_legend").c_str());

        color_i = 0;
        for (auto const &type : types) {
            m_vs_property_hists[property1][property2][type].SetMarkerColor(m_colors[color_i++ % m_number_colors]);
            m_vs_property_hists[property1][property2][type].GetXaxis()->SetTitle(
                m_known_properties[property1].title.c_str());
            m_vs_property_hists[property1][property2][type].GetYaxis()->SetTitle(
                m_known_properties[property2].title.c_str());
            m_vs_property_hists[property1][property2][type].Scale(1.0 / number_of_entries);
            m_vs_property_hists[property1][property2][type].Write();

            m_vs_property_hists[property1][property2][type].SetTitle(m_known_types[type].title.c_str());
            hist_stack.Add(&m_vs_property_hists[property1][property2][type]);
            stack_legend.AddEntry(&m_vs_property_hists[property1][property2][type], m_known_types[type].title.c_str());
        }
        hist_stack.Write();
        stack_legend.Write();
    }
}
