# SHP_genie_data_analysis

A rework of the code from [https://github.com/ltaling/MSc-Thesis](https://github.com/ltaling/MSc-Thesis) which is based on [https://github.com/adishka/e4nu/](https://github.com/adishka/e4nu/) for the purposes of my Senior Honours Project at Edinburgh ( `original` has part of the original code to keep as a reference).
There I am looking at the delta1232 resonances so a substantial part of the code is geared towards that (the `GenieAnalysisDeltaStudies` classes).

The .cpp files in the root folder are executables that each do slightly different parts of the analysis and all use classes based on `GenieAnalysis`.
The likely most generalizeable part is the `GenieAnalysisAutoHistograms` class that has what I call types of events (specifically QE, RES, DELTA1232, MEC, ALL and others can easily be added) and properties (for example event properties from the gst).
Then the user can specify which types to separate and which properties to make which histograms for and they are taken care of, also other types and properties can be added very easily.
