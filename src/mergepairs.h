/* reconstructed from VSEARCH authored by Torbjorn Rognes, Frederic Mahe and Tomas Flouri
 * the original source code can be found in https://github.com/torognes/vsearch
 * I have changed some behaviors and add some other functionality for better usage with this project
 */
#ifndef MERGE_PAIRS_H
#define MERGE_PAIRS_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <string>
#include <limits.h>
#include "util.h"
#include "kmerhash.h"
#include "htslib/kstring.h"
#include "common.h"
#include "options.h"
#include "htmlutil.h"
#include "krec.h"

#ifndef MERGE_REASON_CNT
#define MERGE_REASON_CNT 17
#endif

// merge fail/success reason enum
enum reason_enum{
    undefined, // default value, undefined
    okbykmers, // merged okay by kmers
    okbyoverlap, // merged okay by overlap
    minlen, // input seq too short
    maxlen, // input seq too long
    maxns, // too many Ns in input
    minovlen, // overlap too short
    maxdiffs, // too many differences(mismatch bases counts)
    maxdiffpct, // too high percentage of differences
    staggered, // staggered, rev seq has 3' overhang outside of overlap region
    indel, // indels in overlap region
    repeat, // potential repeats in overlap region / multiple overlaps
    minmergelen, // merged sequence too short
    maxmergelen, // merged sequence too long
    maxee, // expected error too high
    minscore, // alignment score too low, insignificant, potential indel
    nokmers // too few kmers on same diag found
};

// reason to str
static const char* reason_str[MERGE_REASON_CNT] = {
    "undefined",
    "OkByKmer",            "OkByOverlap",     "ReadsTooShort",     "ReadsTooLong",      
    "ReadsTooManyN",       "OverlapTooShort", "TooManyDiff",       "TooHighDiffPct",
    "ReadsStaggered",      "IndelInOverlap",  "RepeatRegion",      "TooShortMergeProduct", 
    "TooLongMergeProduct", "TooHighMismatch", "TooLowAlignScore",  "NoSeedKmers"
};

// type to store merge information
struct merge_data_t{
    char *fwd_header = NULL; // fwd name
    char *rev_header = NULL; // rev name
    char *fwd_sequence = NULL; // fwd seq
    char *rev_sequence = NULL; // rev seq
    char *fwd_quality = NULL; // fwd qual
    char *rev_quality = NULL; // rev qual
    int64_t seq_alloc = 0; // total len allocated to store seq
    int64_t fwd_length = 0; // fwd len
    int64_t rev_length = 0; // rev len
    int64_t fwd_trunc = 0; // fwd len after truncation by qual 
    int64_t rev_trunc = 0; // rev len after truncation by qual
    char *merged_sequence = NULL; // merged seq
    char *merged_quality = NULL; // merged qual
    int64_t merged_length = 0; // merged length
    int64_t merged_seq_alloc = 0; // total length alloated for merged sequence
    double ee_merged = 0; // merged sequence base error probablity sum
    double ee_fwd = 0; // forward sequence merged base error probablity sum
    double ee_rev = 0; // reverse sequence merged base error probablity sum
    int64_t fwd_errors = 0; // fwd base error in merged region(original vs merge generated)
    int64_t rev_errors = 0; // rev base error in merged region(original vs merge generated)
    int64_t offset = 0; // number of bases consumed when merge end from 3'end of rev read
    bool merged = 0; // 1 if merge successfuly, 0 if merge failed
    reason_enum reason = undefined; // reason of merge status

    merge_data_t(){}

    ~merge_data_t(){
        if(merged_sequence){ free(merged_sequence); merged_sequence = NULL; }
        if(merged_quality){ free(merged_quality); merged_quality = NULL; }
    }
};

// class to do pair merge
struct PairMerger{
    /* external options */
    Options* mOpt;
    /* scores in bits */
    int k                        = 5; // kmer length
    double merge_minscore        = .0; // min alignment score needed to merge pair
    double merge_dropmax         = 16.0; // max score drop allowed
    int merge_mindiagcount       = 4; // minimal consecutive kmers needed to merge pair
    int merge_minrepeatdiagcount = 12;// minimal consecutive kmers needed to rep a repeat unit (12+4)
    double merge_mismatchmax     = -4.0;
    /* merge threshold */
    int64_t opt_fastq_ascii = 33; // ASCII Base of fastq qual(33 or 64)
    int64_t opt_fastq_qmin = 0; // min fastq qual val
    int64_t opt_fastq_qmax = 41; // max fastq qual val
    int64_t opt_fastq_qminout = 0; // min fastq qual val output in merged fq
    int64_t opt_fastq_qmaxout = 41; // max fastq qual val output in merged fq
    float opt_fastq_maxee = DBL_MAX; // max error allowed in merged fq
    int64_t opt_fastq_maxdiffs = 10; // max different bases allowed in merge region
    float opt_fastq_maxdiffpct = 25.0; // max didderent base rate allowed in merged region
    int64_t opt_fastq_minovlen = 10; // min overlap length needed in merged region
    int64_t opt_fastq_minmergelen = 0; // minimum merged product length needed
    int64_t opt_fastq_maxmergelen  = 1000000; // max merged product length allowed
    int64_t opt_fastq_minlen = 1; // fastq with length shorter than this will not be merged
    int64_t opt_fastq_maxlen = LONG_MAX; // fastq with length longer than this will not be merged
    bool opt_trunc_by_qual = false; // truncate all 5' region of reads by qual if true
    int64_t opt_fastq_truncqual = LONG_MIN; // truncate all 5' region of reads with base qual less than this
    bool opt_filter_nbase = false; // mask N base qual as 0 and filter reads by max N bases
    int64_t opt_fastq_maxns = LONG_MAX; // total n base count allowed in each read
    bool opt_fastq_allowmergestagger = true; // allow merge reads with 5' part in merge region and 3' part overhang
    int64_t opt_olpm_minolplen = 15;
    float opt_olpm_maxdiffpct = 0.2;
    float opt_olpm_maxseqerr = 0.1;
    /* merge status */
    int64_t merged = 0; // total number of merged pairs
    int64_t notmerged = 0; // total number of pairs not merged
    int64_t total = 0; // total number of pairs processed
    int64_t peovhang = 0; // paired and 3' ends overhang far to 5' ends(usually adapters included)
    /* stat info */
    double sum_read_length = 0.0; // total read length processed
    double sum_squared_fragment_length = 0.0; // square sum of merge product length
    double sum_fragment_length = 0.0; // sum of merged product length
    double sum_ee_fwd = 0.0; // error sum of fowward reads in merged product
    double sum_ee_rev = 0.0; // error sum of reverse reads in merged product
    double sum_ee_merged = 0.0; // error sum of merged product
    uint64_t sum_errors_fwd = 0; // total mismatches in merged region on forward sequence
    uint64_t sum_errors_rev = 0; // total mismatches in merged region on reverse sequence 
    double errors_rate_fwd = 0; // mismatch rate in merged region on forward sequence
    double errors_rate_rev = 0; // mismatch rate in merged region on reverse sequence
    uint64_t merge_ops_cnt[MERGE_REASON_CNT] = {0}; // merge status region count
    uint64_t* merge_len_cnt = NULL; // merged length dist
    int64_t merge_len_min = 0; // minimal merged length
    int64_t merge_len_max = 0; // maxmal merged length
    uint64_t merge_len_mod = 0; // mod merged length
    uint64_t merge_len_quant[4] = {0}; // 25% step quant
    int64_t merge_len_alloc = 0; // merge length dist allocated
    bool summarized = false; // summarized if true
    /* qual operation table */
    char merge_qual_same[128][128]; // merged qual of same bases
    char merge_qual_diff[128][128]; // merged qual of diff bases,(make sure idx first > sec when fetch)
    double match_score[128][128]; // estimated match score of two bases with observed quals
    double mism_score[128][128]; // estimated mismatch score of two bases with ovserved quals
    double q2p[128]; // store qual_chr to err_prob array
    /* kmer hash */
    kh_handle_t* kmerhash = NULL;
    /* merge struct */
    merge_data_t *ip = NULL;

    PairMerger(){
        kmerhash = kh_init();
        ip = new merge_data_t();
        precompute_qual();
    }

    PairMerger(Options* opt){
        mOpt = opt;
        kmerhash = kh_init();
        ip = new merge_data_t();
        precompute_qual();
        opt_fastq_allowmergestagger = opt->mpe.allowMergeStagger;
        opt_fastq_maxdiffs = opt->mpe.maxdiffs;
        opt_fastq_maxdiffpct = opt->mpe.maxdiffpct;
        merge_len_alloc = 1;
        merge_len_cnt = (uint64_t*)calloc(merge_len_alloc, sizeof(uint64_t));
    }

    ~PairMerger(){
        if(kmerhash){
            kh_exit(kmerhash);
            kmerhash = NULL;
        }
        if(ip){
            delete ip;
            ip = NULL;
        }
        if(merge_len_cnt){
            free(merge_len_cnt);
            merge_len_cnt = NULL;
        }
    }

    // base qual to err prob
    double q_to_p(int q);
    // base qual chr to int
    int get_qual(char q);
    // fill merged qual array, score array
    void precompute_qual();
    // merge two overlap base/qual into one
    void merge_sym(char* sym, char* qual, char fwd_sym, char rev_sym, char fwd_qual, char rev_qual);
    // update stat info of successful merged pair
    void keep();
    // update stat info of merge failure pair
    void discard();
    // do merge work in computed range
    void merge();
    // get optimal merge range
    int64_t optimize();
    // do pair merge work
    void process();
    // init ip with krec1_t
    void parse_pair(krec1_t* fwd, krec1_t* rev);
    // keep or discard a pair
    void keep_or_discard();
    // merge two nucleotide krec1_t 
    void merge(krec1_t *fwd, krec1_t *rev);
    // merge two nucleotide strings
    void merge(const std::string& fwd, const std::string& rev);
    // merge two nucleotide cstrings
    void merge(const char* fwd, const char* rev, const char* qf, const char* qr);
    // merge by overlap analysis
    void merge_by_olp();
    // summary
    void summary();
    // generate json format result
    void reportJSON(kstring_t* s, const char* dh, const char* dm);
    void reportHTML(kstring_t* s);
    // generate tsv format head
    static void tsvHead(kstring_t* s);
    // generate tsv format content
    void tsvBody(kstring_t* s);
};

extern int mergepe_main(int argc, char** argv);
extern int mergefq_main(int argc, char** argv);

#endif
