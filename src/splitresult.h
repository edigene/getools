#ifndef SPLITRESULT_H
#define SPLITRESULT_H

#include <string>
#include <vector>
#include <ostream>
#include <iomanip>
#include "krec.h"
#include "qcstat.h"
#include "options.h"
#include "htslib/kstring.h"

// split result
class SplitResult{
    public:
        Options *mOptions; // pointer to Options
        uint64_t mTotalReads; // total reads number feed to spliter
        uint64_t mDropReads; // dropped read number
        uint64_t mSplitFailedReads; // split failed reads
        uint64_t mQCFailedReads; // qc failed reads
        std::vector<uint64_t> mDropCount8Reasons; // split failure by reason counts
        std::vector<uint64_t> mBaseCount; // bases split into each sample counts
        std::vector<uint64_t> mSplitCount; // reads split into each sample counts
        std::vector<kstring_t*> mSplitRead1; // read1 split result of each sample 
        std::vector<kstring_t*> mSplitRead2; // read2 split result of each sample 
        std::vector<QcStat*> mQCRead1; // read1 split result QC
        std::vector<QcStat*> mQCRead2; // read2 split result QC
        std::vector<uint64_t> mFRCount; // forward/reverse reads count in input r1/2
        std::vector<uint64_t> mRFCount; // reverse/forward reads count in input r1/2
        std::vector<uint64_t> mMMCount; // mismatch count of each split read/s
        uint64_t mTotalFR; // total FR
        uint64_t mTotalRF; // total RF
        double mGotRate; // split got rate
        bool mIsPE; // is pe if true
        bool mSummarized; // summarized or not

    public:
        // SplitResult constructor
        SplitResult(Options *opt);

        // SplitResult destructor
        ~SplitResult();

        // init before each run
        void init();

        // increase mDropReads
        void addDropRead(krec1_t* r1, krec1_t* r2, bool isQCFail);

        /** add one read/pair successfuly split to a sample
         * @param r1 read1
         * @param r2 read2
         * @param mm read1/2 total mismatch
         */
        void addSplitRead(krec1_t* r1, krec1_t* r2, int mm);
        
        /** merge a vector of SplitResult obects
         * @param vresults reference of a vector of SplitResult
         * @return pointer to a merged SplitResult
         */
        static SplitResult* mergeResult(const std::vector<SplitResult*>& vresults);

        // summarize
        void summary();

        // report json
        void reportJSON(kstring_t* s, const char* dh, const char* dm);
        void reportHTML(kstring_t* s);

        // report tsv head
        void tsvHead(kstring_t* s, bool m=false);

        // report tsv body
        void tsvBody(kstring_t* s, bool m=false);
};

#endif
