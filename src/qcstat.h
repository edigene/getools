#ifndef QC_STAT_H
#define QC_STAT_H

#include <cstdint>
#include <string>
#include "options.h"
#include "htslib/kstring.h"
#include "util.h"
#include "krec.h"

// QC stat of one read library
struct QcStat{
    Options* mOpt; // options
    size_t mReads; // total raw reads
    size_t mBases; // total raw bases
    size_t mMinReadLen; // min read length
    size_t mMaxReadLen; // max read length
    size_t mMeanReadLen; // mean read length
    int mMinQual; // min qual
    int mMaxQual; // max qual
    float mMeanQual; // mean qual
    size_t mBuflen; // max buffer length
    size_t mCycle; // cycle length
    int mQualBase; // ASCII quality base
    char mQ20Chr; // Q20 ASCII chr
    char mQ30Chr; // Q30 ASCII chr
    char mLowQChr; // low quality chr
    size_t mLowQual; // low quality reads count
    size_t mNReads; // contain N reads
    size_t mQ20Total; // q gt 20 bases
    size_t mQ30Total; // q gt 30 bases
    size_t mLengthSum; // total reads length
    size_t mQualSum; // total quality
    size_t mQ20Bases[5]; // q gt 20 ACGTN count
    size_t mQ30Bases[5]; // q gt 30 ACGTN count
    size_t mBaseContents[5]; // ACGTN count
    double mGCPercent; // GC content
    size_t *mCycleQ20Bases[5];// q gt 20 ACGTN cycle count
    size_t *mCycleQ30Bases[5]; // q gt 20 ACGTN cycle count
    size_t *mCycleBaseContents[5]; // ACGTN cycle count
    size_t *mCycleBaseQuality[5]; // ACGTN cycle qual count
    size_t *mCycleTotalBase; // cycle total base
    size_t *mCycleTotalQuality; // cycle total qual
    std::map<std::string, double*> mContentCurves; // content in each cycle
    std::map<std::string, double*> mQualityCurves; // quality in each cycle
    bool mSummarized; // if summary or not
    bool mGoodRead; //  good read marker for statRead usage

    // constructor
    QcStat(Options* opt, size_t bufLen = 1000);

    // destructor
    ~QcStat();

    // allocate res
    void allocateRes();

    // reallocate res
    void reallocRes(size_t bufLen);

    // stat read, return false if failed QC
    bool statRead(const krec1_t* r, int off = 0);

    // merge bunch of QcStat
    static QcStat* merge(const std::vector<QcStat*>& qcs);

    // summary
    void summary();

    // generate cycle curve
    void cycleCurve();

    // generate json format result
    void reportJSON(kstring_t* s, const char* dh, const char* dm, int r);
    void reportHTML(kstring_t* s, int r);

    // generate tsv format head
    static void tsvHead(kstring_t* s);

    // generate tsv format content
    void tsvBody(kstring_t* s, QcStat* r2);
};

#endif
