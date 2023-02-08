#include <iostream>

#include "GenieAnalysis/GenieAnalysis.h"


class GenieAnalysisDemoW : public GenieAnalysis {
    // "Inherits" all the constructors as we don't need an changes, saves writing
    using GenieAnalysis::GenieAnalysis;
};


int main(int argc, char *argv[]) {

    GenieAnalysisDemoW ga{"/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root"};
    ga.runAnalysis();

    return 0;
}
