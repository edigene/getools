#ifndef GET_SBC_H
#define GET_SBC_H

#include "util.h"
#include "writer.h"
#include "readgen.h"
#include "threadpool.h"
#ifdef GET_SBC_DEBUG
#include "editdistance.h"
#endif
#include "htslib/kstring.h"
#include <unordered_map>
#include <map>

#ifndef MISSIONBIO_BARCODE_COUNT
#define MISSIONBIO_BARCODE_COUNT 1536
#endif

#ifndef MISSIONBIO_BARCODE_CONST_MAX_OFFSET
#define MISSIONBIO_BARCODE_CONST_MAX_OFFSET 3
#endif

#ifndef MISSIONBIO_BARCODE_LENGTH
#define MISSIONBIO_BARCODE_LENGTH 9
#endif

#ifndef MISSIONBIO_BARCODE_INDEX_DEFAULT
#define MISSIONBIO_BARCODE_INDEX_DEFAULT 3
#endif

#ifndef MISSIONBIO_BARCODE_INDEX_MASK
#define MISSIONBIO_BARCODE_INDEX_MASK 0xFFF
#endif

#ifndef MISSIONBIO_BARCODE_COUNT_MASK
#define MISSIONBIO_BARCODE_COUNT_MASK 0xFFF
#endif

#ifndef MISSIONBIO_BARCODE_MISMATCH_MASK
#define MISSIONBIO_BARCODE_MISMATCH_MASK 0x3
#endif

#ifndef MISSIONBIO_INDEX_BCIDX_SHIFT
#define MISSIONBIO_INDEX_BCIDX_SHIFT 20
#endif

#ifndef MISSIONBIO_INDEX_BCCNT_SHIFT
#define MISSIONBIO_INDEX_BCCNT_SHIFT 8
#endif

#define set_barcode_index(x, i, c, m) (x=(i<<MISSIONBIO_INDEX_BCIDX_SHIFT)|(c<<MISSIONBIO_INDEX_BCCNT_SHIFT)|m)
#define get_bardoce_bcindex(x) ((x>>MISSIONBIO_INDEX_BCIDX_SHIFT)&MISSIONBIO_BARCODE_INDEX_MASK)
#define get_barcode_bccount(x) ((x>>MISSIONBIO_INDEX_BCCNT_SHIFT)&MISSIONBIO_BARCODE_COUNT_MASK)
#define get_barcode_mismatch(x) (x&MISSIONBIO_BARCODE_MISMATCH_MASK)

#ifndef BARCODE_MATCH_STATUS_ARR_LEN
#define BARCODE_MATCH_STATUS_ARR_LEN 5
#endif

#ifndef BARCODE_MISMATCH0_CNTARR_INDEX
#define BARCODE_MISMATCH0_CNTARR_INDEX 0
#endif

#ifndef BARCODE_MISMATCH1_CNTARR_INDEX
#define BARCODE_MISMATCH1_CNTARR_INDEX 1
#endif

#ifndef BARCODE_MISMATCH2_CNTARR_INDEX
#define BARCODE_MISMATCH2_CNTARR_INDEX 2
#endif

#ifndef BARCODE_MISMATCH3_CNTARR_INDEX
#define BARCODE_MISMATCH3_CNTARR_INDEX 3
#endif

#ifndef BARCODE_MISMATCH4_CNTARR_INDEX
#define BARCODE_MISMATCH4_CNTARR_INDEX 4
#endif

#ifndef BARCODE_INS1_CNTARR_INDEX
#define BARCODE_INS1_CNTARR_INDEX 3
#endif

#ifndef BARCODE_DEL1_CNTARR_INDEX
#define BARCODE_DEL1_CNTARR_INDEX 4
#endif

#ifndef BARCODE_DROP_STATUS_ARR_LEN
#define BARCODE_DROP_STATUS_ARR_LEN 7
#endif

#ifndef BARCODE_DROP_BY_CONSTSEQ1_NOT_FOUND
#define BARCODE_DROP_BY_CONSTSEQ1_NOT_FOUND 0
#endif

#ifndef BARCODE_DROP_BY_CONSTSEQ2_NOT_FOUND
#define BARCODE_DROP_BY_CONSTSEQ2_NOT_FOUND 1
#endif

#ifndef BARCODE_DROP_BY_NINEBP_MM_TOO_MANY
#define BARCODE_DROP_BY_NINEBP_MM_TOO_MANY 2
#endif

#ifndef BARCODE_DROP_BY_NINEBP_HIT_TOO_MANY
#define BARCODE_DROP_BY_NINEBP_HIT_TOO_MANY 3
#endif

#ifndef BARCODE_DROP_BY_INDEL_TOO_MANY
#define BARCODE_DROP_BY_INDEL_TOO_MANY 4
#endif

#ifndef BARCODE_DROP_BY_INDEL_HIT_TOO_MANY
#define BARCODE_DROP_BY_INDEL_HIT_TOO_MANY 5
#endif

#ifndef BARCODE_DROP_BY_READ_TOO_SHORT
#define BARCODE_DROP_BY_READ_TOO_SHORT 6
#endif

static const char* MISSIONBIO_BARCODE_DROP_REASON[7] = {
    "const_seq1_not_found", "const_seq2_not_found", 
    "nine_bp_mm_too_many",  "nine_bp_hit_too_many",
    "indel_too_many",       "indel_hit_too_many", 
    "read_too_short",
};

static const char* MISSIONBIO_BARCODE_CONST_MAX_OFFSTR = "CTG"; // rev seq of the original "GTC"

// split stat info
struct spl_res_t{
    int64_t totr = 0; // tot reads
    int64_t dropr = 0; // drop reads
    int64_t gotr = 0; // got reads
    int64_t b1mm[BARCODE_MATCH_STATUS_ARR_LEN] = {0}; // mm0, mm1, mm2, ins1, del1 for barcode1
    int64_t b2mm[BARCODE_MATCH_STATUS_ARR_LEN] = {0}; // mm0, mm1, mm2, ins1, del1 for barcode2
    int64_t bbmm[BARCODE_MATCH_STATUS_ARR_LEN] = {0}; // mm0, mm1, mm2, ins1, del1 total for barcode1/2
    int64_t r2drop[BARCODE_DROP_STATUS_ARR_LEN] = {0}; // no cont seq, 9bp mm too large; 9bp hit too many; indel too many, read too short
    int64_t** rcmap = NULL; // read counts of each cell barcodes
    int64_t offc[MISSIONBIO_BARCODE_CONST_MAX_OFFSET+1] = {0}; // offset count

    spl_res_t();
    ~spl_res_t();

    void res2jsn(kstring_t* s);

};

// this is for missionbio single cell barcode split
struct sbc_biom_t{
    int bc1l = MISSIONBIO_BARCODE_LENGTH; // cell barcode1 length
    int bc2l = MISSIONBIO_BARCODE_LENGTH; // cell barcode2 length
    int maxoff = MISSIONBIO_BARCODE_CONST_MAX_OFFSET; // offset between cell barcode1 and const sequence1
    const char* cs1s = "AGTACGTACGAGTC"; // const sequence1 seq
    const char* cs2s = "GTACTCGCAGTAGTC"; // const sequence2 seq
    int cs1l = 14; // const sequence1 length
    int cs2l = 15; // const sequence2 length
    int sbeg1 = 0; // search beg of cs1
    int srng1 = 0; // search len of cs1
    int sbeg2 = 0; // search beg of cs2
    int srng2 = 0; // search len of cs2
    int minl = 0; // min read length needed
    int minbc = 1000; // min barcode support of clean reads
    bool outdrop = false; // out dropped reads if true
    char* inr1 = NULL; // input read1
    char* inr2 = NULL; // input read2
    std::string outdir = "./"; // output directory
    std::string outr1 = "Raw.Fq.R1.gz"; // output read1
    std::string outr2 = "Raw.Fq.R2.gz"; // output read2
    std::string outc1 = "Clean.Fq.R1.gz"; // output read1 with more than N barcode support
    std::string outc2 = "Clean.Fq.R2.gz"; // output read2 with more than N barcode support
    std::string outd1 = "Drop.Fq.R1.gz"; // output dropped read1
    std::string outd2 = "Drop.Fq.R2.gz"; // output dropped read2
    std::string outtsv = "mbspl.tsv"; // output tsv file of stat of raw reads
    std::string outjsn = "mbspl.json"; // output json file of stat of split
    int thread = 4; // split threads
    int maxp = 3; // max packs of reads
    int maxr = 100000;// max reads in each pack
    ThreadPool* tpl = NULL; // thread pool
    Writer* wr1 = NULL; // writer of r1
    Writer* wr2 = NULL; // writer of r2
    Writer* dw1 = NULL; // writer of drop read1
    Writer* dw2 = NULL; // writer of drop read2
    uint64_t totr = 0; // total reads
    uint64_t dropr = 0; // dropped reads
    uint64_t gotr = 0; // valid reads
    std::mutex lock; // mutex lock for write
    std::vector<kstring_t*> br1; // read1 write buffer
    std::vector<kstring_t*> br2; // read2 write buffer
    std::vector<kstring_t*> dbr1; // read1 drop write buffer
    std::vector<kstring_t*> dbr2; // read2 drop write buffer
    std::vector<spl_res_t*> sprs; // split results
    kstring_t* jks = NULL; // json kstring tmp str
    // index for whitelist barcodes
    char* wlist = NULL; // barcode white list
    char* iwidx = NULL; // barcode white list index input
    char* owidx = NULL; // barcode white list index output
    std::vector<std::string> warray; // white list barcode sequence array
    std::vector<uint32_t> wmindex; // white list barcode index for mismatch and exact match
    std::vector<uint32_t> wiindex; // white list barcode index for 1bp insertion
    std::vector<uint32_t> wdindex; // white list barcode idnex for 1bp deletion
    int32_t** bcidx = NULL; // barcode index to NO

    sbc_biom_t();
    ~sbc_biom_t();

    bool valid4spl();
    bool valid4idx();
    
    void init();

    void bcs2idx();
    void bcs2mmidx();
    void bcs2ididx();

    void idx2file();
    void file2idx();

    void split_all();

    void split_rng(ReadPack** rp, int beg, int end, int t);
    
    void raw2clean_all();
    void raw2clean_rng(ReadPack** rp, int beg, int end, int t);
};

void sbc_biom_usage(sbc_biom_t* opt, char* arg0);
int sbc_biom_main(int argc, char** argv);

void sbc_bidx_usage(char* arg0);
int sbc_bidx_main(int argc, char** argv);

#endif
