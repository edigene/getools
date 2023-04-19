#ifndef PRE_IDX_H
#define PRE_IDX_H

#include <vector>
#include "mrm.h"
#include "prefix.h"
#include "common.h"
#include "obwa.h"

#define MIN_SEED_LEN 6

struct pre_idx_t{
    int maxs = 0; // max contigs
    size_t maxl = 0; // max primer length
    size_t minl = 0xffff; // min primer length
    int maxo = 2; // max off set
    bool sematch = false; // single end match
    bool droppre = false; // drop matched part
    int mins = 8; // min seed length
    int maxm = 2; // max mismatch
    bool fwdonly = true; // forward only
    // debug options
    int64_t nom = 0; // no match
    int64_t nop = 0; // no proper match
    char* inr1 = NULL; // input r1
    char* inr2 = NULL; // input r2
    char* cfgf = NULL; // configure file
    // debug options
    std::vector<int> sidx; // bwt primer ref tid to sample index, odd fwd, even rev
    obwa_t* obwaa = NULL; // online bwa to do fwd/rev match
    std::vector<obwa_t*> obwas; // online bwa of partner of each primer
    std::vector<int8_t> strands; // strands of each primer
    std::vector<int> lens; // length of each primer
    std::vector<int> boff; // barcode/index length to skip
    std::vector<bool> seidx; // this pair of primer is just single index

    pre_idx_t();
    ~pre_idx_t();

    void sets(obwa_t* obwa);
    void init(const PrefixList& pfl, int mm, int mo, bool sem = false, bool dpr = false);
    bool match(obwa_t* obwa, krec1_t* r, int& m, int& tid, int& scr);
    bool match(krec1_t* r, int& m);
    bool match(krec1_t* r1, krec1_t* r2, int& m);
    bool dmatch(krec1_t* r1, krec1_t* r2, int& m);
    bool dmatch(obwa_t* obwa, krec1_t* r, MatResutMap& mrs);
};

void pri_usage(pre_idx_t* opt, char* arg0);
int pri_main(int argc, char** argv);

#endif
