#include <iostream>

#include "genie_analysis.h"


int main(int argc, char *argv[]) {

    GenieAnalysis ga{"/home/honza/Sync/University/CurrentCourses/SHP/data/Genie_gst_2000000.root"};
    ga.runAnalysis();

    return 0;
}
