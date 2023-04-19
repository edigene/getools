#ifndef GETOOLS_GS2AS_H
#define GETOOLS_GS2AS_H

#include "util.h"
#include "common.h"

struct ctg_t{
    std::string refn;
    std::vector<std::string> seqns;
    std::vector<std::string> seqss;
    std::vector<int64_t> seqnn;
    std::vector<double> seqps;

    ctg_t(){};
    ~ctg_t(){};

    void cmp_seqnn(double edeff, int64_t totr);
    void out_seqs(std::ofstream& fw, bool details = false);
};

struct gs2as_opt_t{
    std::string igtntsv;
    std::string ogsntsv = "gs2as.simple.tsv";
    std::string ogsdtsv = "gs2as.details.tsv";
    double edeff = 0.8;
    int64_t totr = 10000;

    gs2as_opt_t(){};
    ~gs2as_opt_t(){};

    bool valid();
    void gs2as();
};

void gs2as_usage(gs2as_opt_t* opt, char* arg0);
int gs2as_main(int argc, char** argv);

#endif
