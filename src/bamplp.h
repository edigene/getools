#ifndef BAM_PLP_H
#define BAM_PLP_H

#include "ksw4get.h"
#include "common.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>

#ifndef SC_SNP_CALLED
#define SC_SNP_CALLED 0
#endif

#ifndef SC_INS_CALLED
#define SC_INS_CALLED 1
#endif

#ifndef SC_DEL_CALLED
#define SC_DEL_CALLED 2
#endif

#ifndef SC_DIN_CALLED
#define SC_DIN_CALLED 3
#endif

struct bam_info_aux_t{
    int cc = 0;
    char* bi = NULL;
    int* ofs = NULL;
    int mfs = 0;
    int nfs = 0;
    int* cbc = NULL;
};

struct mplp_aux_t{
    samFile* fp = NULL; // file pointer to BAM
    bam_hdr_t* hdr = NULL; // bam header
    hts_idx_t* sidx = NULL; // sam index
    hts_itr_t* iter = NULL; // iterator of BAM
    std::vector<bam1_t*>* bams = NULL; // bam records
    size_t idx = 0; // bams index
    char* ref = NULL; // ref seq
    KSW_FTYPE dropmask = 0; // dropmask
};

struct mplp_pileup_t{
    int n;
    int *n_plp = NULL;
    const bam_pileup1_t** plp = NULL;
};

struct call_ret_t{
    int tid = -1;
    int pos = -1;
    int32_t depth = 0;
    int32_t totcnt = 0;
    int32_t snp[5] = {0}; // ACGTN
    float af[10] = {.0}; // dynamic usage
    int32_t ad[10] = {0}; // dynamic usage
    int32_t dp[2] = {0}; // dynamic usage
    int32_t tr[2] = {0}; // dynamic usage
    int ridx = 4;
    std::unordered_map<std::string, int64_t> ins;
    std::unordered_map<int64_t, int64_t> del;
    std::unordered_map<std::string, int64_t> din;
    uint8_t vm = 0;
    kstring_t s = {0, 0, 0};
    // single cell count
    uint8_t sc_called[4] = {0};
    std::vector<std::unordered_map<std::string, int>> sc_snp;
    std::unordered_map<std::string, std::unordered_map<std::string, int>> sc_ins;
    std::unordered_map<int64_t, std::unordered_map<std::string, int>> sc_del;
    std::unordered_map<std::string, std::unordered_map<std::string, int>> sc_din;
    kstring_t* ts = NULL;

    call_ret_t(){
        ts = (kstring_t*)calloc(1, sizeof(kstring_t));
    }
    
    ~call_ret_t(){
        if(s.s) free(s.s);
        if(ts){
            if(ts->s) free(ts->s);
            free(ts);
        }
    }

    void init(){
        memset(snp, 0, 5 * sizeof(int32_t));
        memset(sc_called, 0, 4 * sizeof(uint8_t));
        sc_snp.clear();
        sc_ins.clear();
        sc_del.clear();
        sc_din.clear();
        tid = -1;
        pos = -1;
        depth = 0;
        ridx = -1;
        ins.clear();
        del.clear();
        din.clear();
        vm = 0;
        s.l = 0;
        sc_ins.clear();
        sc_del.clear();
        sc_din.clear();
    };

};

#ifndef SNP_CALLED
#define SNP_CALLED  0x1
#endif

#ifndef INS_CALLED
#define INS_CALLED  0x2
#endif

#ifndef DEL_CALLED
#define DEL_CALLED  0x4
#endif

#ifndef DIN_CALLED
#define DIN_CALLED  0x8
#endif

#ifndef GEV_CALLED
#define GEV_CALLED  0x10
#endif

#ifndef VAR_CALLED
#define VAR_CALLED  0xf
#endif

void mplp_aux_destroy(mplp_aux_t* ma);
void mplp_pileup_init(mplp_pileup_t* mp, int _n);
void mplp_pileup_destroy(mplp_pileup_t* mp);
int mplp_func(void *data, bam1_t* b);
int pileup_constructor(void *data, const bam1_t *b, bam_pileup_cd *cd);
int pileup_destructor(void *data, const bam1_t *b, bam_pileup_cd *cd);
void plp_call(int tid, int pos, const bam_pileup1_t* plp, int n_plp, call_ret_t* cr, uint16_t vtag);
int call2bcf1(call_ret_t* cr, bcf_hdr_t* hdr, std::vector<bcf1_t*>& recs, char* ref);
int call2bcf2(call_ret_t* cr1, call_ret_t* cr2, bcf_hdr_t* hdr, std::vector<bcf1_t*>& recs, char* ref);

#endif
