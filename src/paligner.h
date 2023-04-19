#ifndef PALIGNER_H
#define PALIGNER_H

#include "ksw4get.h"
#include "krec.h"
#include "options.h"
#include "amplicon.h"
#include <vector>
#include <array>
#include <tuple>
#include "bamsorter.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"


// class to detect edit efficience of one ref
struct PAligner{
    // option
    Options* opt;
    // ref/sgrn
    char* ref = NULL;
    char* name = NULL;
    int32_t rtid;
    // score strategy
    int match = 5;
    int mismatch = 4;
    int gapopen = 25;
    int gapext = 0;
    // align variables
    int8_t* score_mat = NULL;
    uint8_t* refints = NULL;
    int rlen = 0;
    uint8_t flag = 1;
    // bams
    std::vector<bam1_t*> bams;
    bool vsorted = false;
    // stats se
    int64_t* nmmdist_se = NULL;
    int64_t* ninsdist_se = NULL;
    int64_t* ndeldist_se = NULL;
    int64_t* nvardist_se = NULL;
    int64_t* mmdist_se = NULL;
    // stats pe
    int64_t* nmmdist_pe = NULL;
    int64_t* ninsdist_pe = NULL;
    int64_t* ndeldist_pe = NULL;
    int64_t* nvardist_pe = NULL;
    int64_t* mmdist_pe = NULL;
    int64_t mapped = 0;
    int64_t unmapped = 0;
    int maxdstl_se = 0; // max useful dist array length for se
    int maxdstl_pe = 0; // max useful dist array lengto for pe
    bool summarized = false;

    // constructor
    PAligner(){};

    // constructor
    PAligner(amplicon_t* a, Options* op){
        opt = op;
        match = opt->aln.match;
        mismatch = opt->aln.mismatch;
        gapopen = opt->aln.gapopen;
        gapext = opt->aln.gapext;
        ref = strdup(a->aseq.c_str());
        name = strdup(a->aname.c_str());
        rlen = a->aseq.size();
        score_mat = kswge_gen_smat(match, mismatch);
        refints = kswge_seq2ints(ref, rlen);
        nmmdist_se = (int64_t*)calloc(rlen, sizeof(int64_t));
        ninsdist_se = (int64_t*)calloc(rlen, sizeof(int64_t));
        ndeldist_se = (int64_t*)calloc(rlen, sizeof(int64_t));
        nvardist_se = (int64_t*)calloc(rlen, sizeof(int64_t));
        mmdist_se = (int64_t*)calloc(rlen, sizeof(int64_t));
        nmmdist_pe = (int64_t*)calloc(rlen, sizeof(int64_t));
        ninsdist_pe = (int64_t*)calloc(rlen, sizeof(int64_t));
        ndeldist_pe = (int64_t*)calloc(rlen, sizeof(int64_t));
        nvardist_pe = (int64_t*)calloc(rlen, sizeof(int64_t));
        mmdist_pe = (int64_t*)calloc(rlen, sizeof(int64_t));
    }

    // destructor
    ~PAligner(){
        if(ref){ free(ref); ref = NULL; }
        if(name){ free(name); name = NULL; }
        if(score_mat){ free(score_mat); score_mat = NULL; }
        if(refints){ free(refints); refints = NULL; }
        if(nmmdist_se){ free(nmmdist_se); nmmdist_se = NULL; }
        if(ninsdist_se){ free(ninsdist_se); ninsdist_se = NULL; }
        if(ndeldist_se){ free(ndeldist_se); ndeldist_se = NULL; }
        if(nvardist_se){ free(nvardist_se); nvardist_se = NULL; }
        if(mmdist_se){ free(mmdist_se); mmdist_se = NULL; }
        if(nmmdist_pe){ free(nmmdist_pe); nmmdist_pe = NULL; }
        if(ninsdist_pe){ free(ninsdist_pe); ninsdist_pe = NULL; }
        if(ndeldist_pe){ free(ndeldist_pe); ndeldist_pe = NULL; }
        if(nvardist_pe){ free(nvardist_pe); nvardist_pe = NULL; }
        if(mmdist_pe){ free(mmdist_pe); mmdist_pe = NULL; }
    }

    // generate bam header
    bam_hdr_t* get_bam_hdr();
    
    // align to kswr_t
    kswr_t* align(char* qseq, int qlen);

    // ksw to bam
    bam1_t* ksw2bam(kswr_t* ret, char* qseq, int qlen, char* qname, int qnlen, int mask);

    // add one bam to store
    void add_var(bam1_t* b);

    // stat muts
    void stat_var(kswr_t* r1, kswr_t* r2);

    // summary
    void summary();

    // generate json format result
    void reportJSON(kstring_t* s, const char* dh, const char* dm);

    // generate tsv format head
    static void tsvHead(kstring_t* s);

    // generate tsv format content
    void tsvBody(kstring_t* s);

    // generate tsv format head of diffx
    static void tsv4DiffxHead(kstring_t* s, int ml);

    // generate tsv format content of diffx
    void tsv4DiffxBodySE(kstring_t* s);
    void tsv4DiffxBodyPE(kstring_t* s);

};

#endif
