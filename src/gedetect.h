#ifndef GEDETECT_H
#define GEDETECT_H

#include "cnt.h"
#include "ksw4get.h"
#include "krec.h"
#include "vars.h"
#include "kmeans.h"
#include "options.h"
#include "amplicon.h"
#include <vector>
#include <array>
#include <tuple>
#include "hapcnt.h"
#include "bamsorter.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"

#define TOTAL_DROP_REASON 6

// drop reason stat
struct drop_schem_t{
    int KSW_FHIERRLOWQ_IDX = 0;
    int KSW_FPRIMDIMER_IDX = 1;
    int KSW_FMISPRIMINN_IDX = 2;
    int KSW_FMISPRIMOUT_IDX = 3;
    int KSW_FLOWSURPT_IDX = 4;
    int KSW_FMANYVARS_IDX = 5;
    int tot_reason = TOTAL_DROP_REASON;
};

const static char* drop_reasons_str[TOTAL_DROP_REASON] = {
    "LowQualHighVar",  "PrimerDimer",  "ConsecutiveInDel",
    "MissPrime",       "LowSupport",   "TooManyVars",
};

// single cell result
struct sc_ge_t{
    int32_t edicnt = 0; // edited
    int32_t refcnt = 0; // reftype
    int32_t recpct = 0; // recombination(partial match, any hit)
    int32_t reccnt = 0; // recombination(partial match, all hit)
    int32_t recent = 0; // recombination(exact match)
    int32_t othcnt = 0; // others
    int32_t droped = 0; // dropped
};

// class to detect edit efficience of one ref
struct GEDetector{
    // option
    Options* opt;
    // hapcnt
    hap_opt_t* hap = NULL;
    // ref/sgrn
    char* ref = NULL;
    char* name = NULL;
    char* sgr = NULL;
    char* donorseq = NULL;
    int donorlen = 0;
    int donorbeg = 0;
    int donorenda = 0;
    int donorendd = 0;
    std::vector<std::string> donorins;
    std::vector<std::string> donordis;
    std::vector<int32_t> donordel;
    std::vector<int32_t> donorsnv;
    int donormutcnt = 0;
    int32_t rtid;
    int32_t sgrbeg;
    int32_t sgrend;
    int32_t clsbeg;
    int32_t clsend;
    int32_t cutbuffer;
    int32_t vvbeg; // valid var beg
    int32_t vvend; // valid var end
    int vvmask;
    int fplen = 0;
    int rvlen = 0;
    int mpmlen = 30;
    bool greedy = false;
    float minaf = 0.01;
    int mincnt = 2;
    int maxvar = 10;
    // score strategy
    int match = 5;
    int mismatch = 4;
    int gapopen = 25;
    int gapext = 0;
    // align variables
    int8_t* score_mat = NULL;
    uint8_t* donorints = NULL;
    uint8_t* refints = NULL;
    int rlen = 0;
    uint8_t flag = 1;
    uint16_t filters = 0;
    int32_t filterd = 0;
    int32_t masklen = 15;
    int maxciter = 10000;
    uint32_t mincpnt = 1000;
    bool cfworked = false;
    int maxdel = 0xffff;
    float seqerr = 0.01;
    // variants
    KSW_FTYPE dropmask = 0;
    KSW_FTYPE varmask = 0;
    KSW_FTYPE recmask = 0;
    std::vector<bam1_t*> vars;
    std::vector<bcf1_t*> bcfs;
    int vsorted = 0;
    // gene edit computation
    int clspos = 0; // cleavage site
    int64_t dropcnt = 0; // filtered out sequcnes count
    int64_t refcnt = 0; // reftype reads
    int64_t recpct = 0; // recombination(partial match, any hit)
    int64_t reccnt = 0; // recombination(partial match, all hit)
    int64_t recent = 0; // recombination(exact match)
    int64_t othcnt = 0; // other reads
    int64_t edicnt = 0; // edit event supporting event, limit to [sgrbeg-cutbuffer, sgrend+cutbuffer]
    int64_t totcnt = 0; // total valid reads(ref+alt+other)
    int64_t allcnt = 0; // total reads assigned to this amplican
    int64_t mutcnt = 0; // edit event supporting event, in all amplicon range
    int64_t ctrltt = 0; // control total valid sequences

    double edieff = .0; // edit efficience
    double recpef = .0; // recombination efficience(partial match, any hit)
    double receff = .0; // recombination efficience(partial match, all hit)
    double reeeff = .0; // recombination efficience(exact match);
    double muteff = .0; // mut rate across all amplicon range
    int64_t maxvc = 0; // max variant count, used in kcluster
    kstring_t* kmcjs = NULL; // kcluster html output
    // cnt
    var_num_t* tot_cnt_dist = NULL; // variant count dist
    var_num_t* ins_cnt_dist = NULL; // insertion count dist
    var_num_t* del_cnt_dist = NULL; // deletion count dist
    var_num_t* din_cnt_dist = NULL; // delins count dist
    var_num_t* snv_cnt_dist = NULL; // snv count dist
    var_num_t* scl_cnt_dist = NULL; // sc count dist
    // len
    var_num_t* ins_len_dist = NULL; // insertion length dist
    var_num_t* del_len_dist = NULL; // deletion length dist
    var_num_t* din_len_dist = NULL; // delins length dist
    var_num_t* scl_len_dist = NULL; // sc length dist
    // pos
    var_num_t* ins_pos_dist = NULL; // insertion position dist
    var_num_t* del_pos_dist = NULL; // deletion position dist
    var_num_t* din_pos_dist = NULL; // delins position dist
    var_num_t* snv_pos_dist = NULL; // snv position dist
    var_num_t* scl_pos_dist = NULL; // sc position dist
    // pos-wise avlen
    var_num_t* ins_lps_dist = NULL; // position-wise ins length sum
    var_num_t* del_lps_dist = NULL; // position-wise del length sum
    var_num_t* din_lps_dist = NULL; // psoition-wise delins length sum
    int64_t** nccnt = NULL; // ACGTN count on each cycle
    // reftype freq
    double* refrate = NULL;
    int64_t* refcov = NULL;
    int64_t* altcov = NULL;
    bool summarized = false;
    // single cell
    std::unordered_map<std::string, sc_ge_t*> scretm; // single cell results
    // drop reason count
    drop_schem_t dropschem;
    int64_t droprcnt[TOTAL_DROP_REASON] = {0};
    // ctrl var
    std::string cbam;
    SeqVarCntMap cins;
    SeqVarCntMap cdis;
    LenVarCntVector cdel;
    SNPVarCntVector csnv;
    int minctrlc;

    // constructor
    GEDetector(){};

    // constructor
    GEDetector(amplicon_t* a, Options* op){
        opt = op;
        for(int i = 0; i < TOTAL_DROP_REASON; ++i) droprcnt[i] = 0;
        match = opt->aln.match;
        dropmask = (KSW_FHIERRLOWQ | KSW_FPRIMDIMER | KSW_FMISPRIMINN | KSW_FMISPRIMOUT | KSW_FLOWSURPT | KSW_FMANYVARS);
        vvmask = opt->edo.vartypem;
        varmask = recmask = 0;
        if(vvmask & GEVAR_SNV) varmask |= KSW_FSNVINRNG;
        if(vvmask & GEVAR_INS) varmask |= KSW_FINSINRNG;
        if(vvmask & GEVAR_DEL) varmask |= KSW_FDELINRNG; 
        if(vvmask & GEVAR_DIN) varmask |= KSW_FDIINRNG;
        switch(opt->edo.rectypem){
            case REC_EXACT_IS_EDIT:
                recmask |= KSW_FRECEXACT;
                break;
            case REC_ALLHIT_IS_EDIT:
                recmask |= KSW_FRECALLHIT;
                break;
            case REC_ANYHIT_IS_EDIT:
                recmask |= KSW_FRECANYHIT;
                break;
        }
        if(a->donor.size()){
            donorseq = strdup(a->donor.c_str());
            donorlen = a->donor.size();
            donorints = kswge_seq2ints(donorseq, donorlen);
            donorbeg = a->dbeg;
            donorenda = a->denda;
            donorendd = a->dendd;
            donorins = a->dins;
            donordis = a->ddis;
            donordel = a->ddel;
            donorsnv = a->dsnv;
            donormutcnt = a->dmutcnt;
        }else{
            donorseq = NULL;
            donorints = NULL;
        }
        mismatch = opt->aln.mismatch;
        gapopen = opt->aln.gapopen;
        gapext = opt->aln.gapext;
        maxdel = a->maxdel;
        minaf = op->edo.minaf;
        mincnt = op->edo.minseqc;
        maxvar = opt->edo.maxvc;
        mpmlen = op->edo.mpmatch;
        greedy = op->edo.greedy;
        rvlen = a->rvlen;
        fplen = a->fplen;
        // beg update mpmlen
        if(mpmlen == 0){
            mpmlen = rvlen;
            if(mpmlen > fplen) mpmlen = fplen;
        }
        // end update mpmlen
        ref = strdup(a->aseq.c_str());
        name = strdup(a->aname.c_str());
        sgr = strdup(a->sgrseq.c_str());
        rlen = a->aseq.size();
        sgrbeg = a->sgrbeg;
        sgrend = a->sgrend;
        clsbeg = a->clsbeg;
        clsend = a->clsend;
        cutbuffer = opt->edo.cutbuflen;
        vvbeg = MAX(0, clsbeg - cutbuffer);
        vvend = MIN(clsend + cutbuffer, rlen-1);
        score_mat = kswge_gen_smat(match, mismatch);
        refints = kswge_seq2ints(ref, rlen);
        // hapcnt
        if(opt->edo.hapcnt && donorseq){
            hap = new hap_opt_t();
            hap->beg = sgrbeg+1;
            hap->end = sgrend+1;
            hap->far = strdup(ref);
            hap->faa = strdup(donorseq);
            hap->aps.resize(a->dsnv.size(), 4);
            for(size_t iii = 0; iii < a->dsnv.size(); ++iii){
                if(a->dsnv[iii] >= 0){
                    hap->aps[iii] = a->dsnv[iii];
                }
            }
            hap->gtn = 1000;
            hap->maxh = 14;
            hap->minr = 0;
            hap->otn = -1;
            hap->extsnv = opt->extsnv;
            hap->outdir = util::joinpath(opt->hapdir, a->aname);
            hap->jscdn = opt->hrjsn;
            hap->update();
            hap->dropmask = dropmask;
        }
        // count 
        tot_cnt_dist = new var_num_t(rlen);
        ins_cnt_dist = new var_num_t(rlen);
        del_cnt_dist = new var_num_t(rlen);
        din_cnt_dist = new var_num_t(rlen);
        snv_cnt_dist = new var_num_t(rlen);
        scl_cnt_dist = new var_num_t(rlen);
        // length
        ins_len_dist = new var_num_t(rlen);
        del_len_dist = new var_num_t(rlen);
        din_len_dist = new var_num_t(rlen);
        scl_len_dist = new var_num_t(rlen);
        // pos
        ins_pos_dist = new var_num_t(rlen);
        del_pos_dist = new var_num_t(rlen);
        din_pos_dist = new var_num_t(rlen);
        snv_pos_dist = new var_num_t(rlen);
        scl_pos_dist = new var_num_t(rlen);
        // lps
        ins_lps_dist = new var_num_t(rlen);
        del_lps_dist = new var_num_t(rlen);
        din_lps_dist = new var_num_t(rlen);
        // nccnt
        nccnt = (int64_t**)malloc(5 * sizeof(int64_t*));
        for(int i = 0; i < 5; ++i) nccnt[i] = (int64_t*)calloc(rlen, sizeof(int64_t));
        // refrate
        refrate = (double*)calloc(rlen, sizeof(double));
        refcov = (int64_t*)calloc(rlen, sizeof(int64_t));
        altcov = (int64_t*)calloc(rlen, sizeof(int64_t));
        // ctrl
        if(opt->ctrl.size()) cbam = opt->ctrl;
        minctrlc = opt->edo.ctfc;
        // kmeans cluster js
        kmcjs = (kstring_t*)calloc(1, sizeof(kstring_t));
    }

    // destructor
    ~GEDetector(){
        if(ref){ free(ref); ref = NULL; }
        if(name){ free(name); name = NULL; }
        if(sgr){ free(sgr); sgr = NULL; }
        if(donorseq){ free(donorseq); donorseq = NULL; }
        if(score_mat){ free(score_mat); score_mat = NULL; }
        if(refints){ free(refints); refints = NULL; }
        if(donorints){ free(donorints); donorints = NULL; }
        // count
        if(tot_cnt_dist){ delete tot_cnt_dist; tot_cnt_dist = NULL; }
        if(ins_cnt_dist){ delete ins_cnt_dist; ins_cnt_dist = NULL; }
        if(del_cnt_dist){ delete del_cnt_dist; del_cnt_dist = NULL; }
        if(din_cnt_dist){ delete din_cnt_dist; din_cnt_dist = NULL; }
        if(snv_cnt_dist){ delete snv_cnt_dist; snv_cnt_dist = NULL; }
        if(scl_cnt_dist){ delete scl_cnt_dist; scl_cnt_dist = NULL; }
        // length
        if(ins_len_dist){ delete ins_len_dist; ins_len_dist = NULL; }
        if(del_len_dist){ delete del_len_dist; del_len_dist = NULL; }
        if(din_len_dist){ delete din_len_dist; din_len_dist = NULL; }
        if(scl_len_dist){ delete scl_len_dist; scl_len_dist = NULL; }
        // pos
        if(ins_pos_dist){ delete ins_pos_dist; ins_pos_dist = NULL; }
        if(del_pos_dist){ delete del_pos_dist; del_pos_dist = NULL; }
        if(din_pos_dist){ delete din_pos_dist; din_pos_dist = NULL; }
        if(snv_pos_dist){ delete snv_pos_dist; snv_pos_dist = NULL; }
        if(scl_pos_dist){ delete scl_pos_dist; scl_pos_dist = NULL; }
        // lps
        if(ins_lps_dist){ delete ins_lps_dist; ins_lps_dist = NULL; }
        if(del_lps_dist){ delete del_lps_dist; del_lps_dist = NULL; }
        if(din_lps_dist){ delete din_lps_dist; din_lps_dist = NULL; }
        // nccnt
        if(nccnt){
            for(int i = 0; i < 5; ++i){ if(nccnt[i]){ free(nccnt[i]); nccnt[i] = NULL; }}
            free(nccnt); nccnt = NULL;
        }
        // refrate
        if(refrate){ free(refrate); refrate = NULL; }
        if(refcov){ free(refcov); refcov = NULL; }
        if(altcov){ free(altcov); altcov = NULL; }
        // empty sc result
        for(auto& p: scretm){
            if(p.second){
                free(p.second);
                p.second = NULL;
            }
        }
        scretm.clear();
        // kmeans cluster js
        if(kmcjs){
            if(kmcjs->s){ free(kmcjs->s); kmcjs->s = NULL; }
            free(kmcjs);
        }
    }

    // generate bam header
    bam_hdr_t* get_bam_hdr();
    
    // align to kswr_t
    kswr_t* align(char* qseq, int qlen);

    // ksw to bam
    bam1_t* ksw2bam(kswr_t* ret, char* qseq, int qlen, char* qname, int qnlen, int mask);

    // add variant alignment result
    void add_var(bam1_t* b);

    // compute cleavage site
    void cal_clspos();

    // k-means cluster and filter
    void kcluster();

    // mark var in range
    void mark_var();

    // get ins del snv count/ len of an bam
    void stat_var(bam1_t* b, bool toadd);

    // miss priming event detection by extra match
    void mark_misprimeout();

    // miss priming event detection by IMD, DMI
    void mark_misprimeinn();

    // miss priming marking 
    void mark_misprime(){
        mark_misprimeout();
        mark_misprimeinn();
    }

    // count variants in one alignment and mark range status
    void var_count(bam1_t* b);

    // fix pe info of KSW_FMISPRIMINN, KSW_FPRIMDIMER, KSW_FHIERRLOWQ
    void fixpe();

    // count variants in two alignments and mark range status
    void var_count(bam1_t* b1, bam1_t* b2);

    // compute edit efficiency
    void cal_edieff(); // bulk
    void cal_scedff(); // single cell level

    // hapcnt
    void hapcnt(){
        if(hap){
            hap->hapcnt8bv(vars, rlen);
            hap->hap2tsv();
            hap->bam2tsv();
            hap->aas2tsv();
            hap->hap2jsn();
            hap->stat2tsv();
            hap->hap2html();
            hap->snv2tsv();
            hap->aac2tsv();
            hap->snv2foc();
            hap->bias2tsv();
            delete hap;
            hap = NULL;
        }
    }

    // af filter, must be done after cal_edieff
    void af_filter();

    // summary
    void summary();

    // generate json format result
    void reportJSON(kstring_t* s, const char* dh, const char* dm);
    void reportHTML(kstring_t* s);

    // generate tsv format head
    static void tsvHead(kstring_t* s);

    // generate tsv format content
    void tsvBody(kstring_t* s);

    // ctrl bam to vars
    void cbam2vars();

    // mark ref from ctrl
    void markctrlref();

    // bam2bcf
    void bam2bcf();

    // getttr
    int64_t get_ttr(const char* bam);
};

#endif
