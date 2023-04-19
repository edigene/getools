#ifndef SPLIT_REPORT_H
#define SPLIT_REPORT_H

#include "options.h"
#include "mergepairs.h"
#include "splitresult.h"

// reporter 
struct Reporter{
    Options* mOpt; // opt

    // constructor
    Reporter(Options* opt){
        mOpt = opt;
    }

    // destructor
    ~Reporter(){}

    // report json
    void reportJSON(SplitResult* sr, QcStat* qs1, QcStat* qs2, PairMerger* merge);
    // report tsv
    void reportTSV(SplitResult* sr, QcStat* qs1, QcStat* qs2, PairMerger* merge);
    // report
    void report(SplitResult* sr, QcStat* qs1, QcStat* qs2, PairMerger* merge);
};

#endif
