#ifndef AMP_IDX_H
#define AMP_IDX_H

#include <vector>
#include "mrm.h"
#include "util.h"
#include "common.h"
#include "amplicon.h"
#include "obwa.h"

struct amp_idx_t{
    bool sematch = false; // single end match
    int maxs = 0; // max contigs
    int minseed = 19; // min match needed
    int minscore = 19; // min score output
    obwa_t* obwa; // online bwa aligner
    std::vector<int> sidx; // ref tid to sample index

    amp_idx_t();
    ~amp_idx_t();

    void sets();
    void init(const AmpliconList& ampl, bool sem = false);
    bool match(krec1_t* r);
    bool match(krec1_t* r1, krec1_t* r2);
    void msecond(krec1_t* r, MatResutMap& mrs);
};

#endif
