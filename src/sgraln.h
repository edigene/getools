#ifndef SGR_ALN_H
#define SGR_ALN_H

#include <stdlib.h>
#include "htslib/kstring.h"
#include "util.h"
#include "krec.h"
#include "ksw4get.h"

struct sgr_seq_t{
    char* sgrfwd = NULL;
    char* sgrrev = NULL;
    uint8_t* sgrintf = NULL;
    uint8_t* sgrintr = NULL;
    kswr_t* faln = NULL;
    kswr_t* raln = NULL;
    kswr_t* best = NULL;
    int sgrlen = 0;
    int maxscr = 0;
    bool maxisr = false;

    sgr_seq_t(){
    }

    sgr_seq_t(const char* seq){
        sgrfwd = strdup(seq);
        sgrlen = strlen(sgrfwd);
        sgrrev = util::revComp2NewSeq(sgrfwd, sgrlen);
        sgrintf = kswge_seq2ints(sgrfwd, sgrlen);
        sgrintr = kswge_seq2ints(sgrrev, sgrlen);
    }

    ~sgr_seq_t(){
        if(sgrfwd){ free(sgrfwd); sgrfwd = NULL; }
        if(sgrrev){ free(sgrrev); sgrrev = NULL; }
        if(sgrintf){ free(sgrintf); sgrintf = NULL; }
        if(sgrintr){ free(sgrintr); sgrintr = NULL; }
        if(faln){ kswr_destroy(faln); faln = NULL; }
        if(raln){ kswr_destroy(raln); raln = NULL; }
        best = NULL;
    }

    void clear(){
        if(faln){ kswr_destroy(faln); faln = NULL; }
        if(raln){ kswr_destroy(raln); raln = NULL; }
        best = NULL;
    }

};

struct sgr_aln_opt{
    // io
    char* inref = NULL; // references fa to find sgrna
    char* sgrfwd = NULL; // sgrna_pam sequence to find sgrna
    std::vector<sgr_seq_t*> sgrseqs; // sgrna seqs
    int8_t* mat = NULL;
    char* outref = NULL; // out reference with sgrna_pam uppercased
    char* outtsv = NULL; // out tsv file
    char* outaln = NULL; // out aln file
    int breakp = 0; // break point
    bool getbp = false; // break point passed
    int flk = 30; // flank length
    // score
    int match = 4; 
    int mismatch = 3;
    int gapo = 3;
    int gape = 3;
    int minscore = 0;

    sgr_aln_opt();
    ~sgr_aln_opt();

    bool valid();
    void init();
    void search();
};

void sgraln_usage(sgr_aln_opt* opt, char* arg0);
int sgraln_main(int argc, char** argv);

#endif
