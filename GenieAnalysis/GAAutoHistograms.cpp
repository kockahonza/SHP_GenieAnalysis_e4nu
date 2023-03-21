#include "GAAutoHistograms.h"

void GAAutoHistograms::prepareAutoHists() {
    string name;
    string title;
    Int_t nbinsx;
    Double_t xlow, xup;

    if (m_properties.empty()) {
        for (const auto &known_property_keyval : m_known_properties) {
            m_properties.push_back(known_property_keyval.first);
        }
    } else {
        for (const string &property : m_properties) {
            if (m_known_properties.count(property) != 1) {
                throw std::runtime_error("There is an unknown property \"" + property + "\" given to GAAutoHistograms");
            }
        }
    }

    if (m_types.empty()) {
        for (const auto &known_type_keyval : m_known_types) {
            m_types.push_back(known_type_keyval.first);
        }
    } else {
        for (const string &type : m_types) {
            if (m_known_types.count(type) != 1) {
                throw std::runtime_error("There is an unknown type \"" + type + "\" given to GAAutoHistograms");
            }
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
        if (m_known_properties.count(property1) != 1) {
            throw std::runtime_error("There is an unknown property \"" + property1 +
                                     "\" given to GAAutoHistograms for a vs plot");
        }
        if (m_known_properties.count(property2) != 1) {
            throw std::runtime_error("There is an unknown property \"" + property2 +
                                     "\" given to GAAutoHistogram for a vs plots");
        }

        if (types.empty()) {
            for (const auto &known_type_keyval : m_known_types) {
                types.push_back(known_type_keyval.first);
            }
        } else {
            for (const string &type : types) {
                if (m_known_types.count(type) != 1) {
                    throw std::runtime_error("There is an unknown type \"" + type +
                                             "\" given to GAAutoHistograms for a vs plot");
                }
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

void GAAutoHistograms::runPostAnalysis(const Long64_t &number_of_entries) {
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
