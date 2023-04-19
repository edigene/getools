#ifndef GET_B2G_H
#define GET_B2G_H

#include <stdio.h>
#include <libgen.h>
#include <unistd.h>
#include <map>
#include "util.h"
#include "common.h"
#include "htmlopt.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

#ifndef GEVAR_SNV
#define GEVAR_SNV 0x1
#endif

#ifndef GEVAR_INS
#define GEVAR_INS 0x2
#endif

#ifndef GEVAR_DEL
#define GEVAR_DEL 0x4
#endif

#ifndef GEVAR_DIN
#define GEVAR_DIN 0x8
#endif

static const char* GEV_STR_ARR4_B2G[9] = {
    "REF", "SNV", "INS", "", "DEL",
    "-"  , "-"  , "-"  , "DIN"};
    
// variant type
struct b2g_var_t{
    std::string ref; // ref
    std::string alt; // alt
    int32_t start; // start
    int32_t end; // end
    double freq; // freq
    int width; // width
    std::string type; // type

    b2g_var_t(){
        ref = "REF";
        alt = "ALT";
        start = 0;
        end = 0;
        freq = .0;
        width = 0;
        type = "Others";
    }

    ~b2g_var_t(){}
};

// variant sorter
struct sort8freq{
    bool operator()(const b2g_var_t* v1, const b2g_var_t* v2) const {
        return v1->freq > v2->freq;
    }
};


// options 
struct b2g_opt_t{
    char* inbcf = NULL;
    char* infa = NULL;
    std::string outdir = "./";
    std::string outfa = "gtn.fa";
    std::string outtsv = "gtn.tsv";
    faidx_t* fai = NULL;
    int nctc = 3;
    int type = 6;
    int topn = 50;

    b2g_opt_t();
    ~b2g_opt_t();

    bool valid();
    void init();
    void gtn2fa();
};

void b2g_usage(b2g_opt_t* opt, char* arg0);

int b2g_main(int argc, char** argv);

#endif
