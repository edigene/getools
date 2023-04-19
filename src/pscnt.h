#ifndef PSCNT_H
#define PSCNT_H

#include "ksw4get.h"
#include "util.h"
#include "common.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include <unordered_map>

struct ps_t{
    std::string smp;
    std::string pat;
    std::string anno;
    int64_t cnt;
};

struct pscnt_opt_t{
    std::string inpaf; // input patname,patseq to annotate patstrs
    std::string inref; // indexed reference of bam
    std::string inbam; // input var.bam from getools
    std::string outdir = "./"; // output directory
    std::string outf = PSCNT_PSC_TSV; // output stat file
    std::string outs = PSCNT_RDC_TSV; // reads cnt
    std::string outx = PSCNT_NXC_TSV; // pos-wise nuc cnt
    int skr = -1; // use read1 if true
    int pos = 1; // 1based pos to look for pattern from there on
    int len = 10; // pattern length
    float ppid = 0.9; // percentage before match
    bool ngt = false; // bam not from getools caledit
    bool nidl = false; // skip reads with indel
    bool npidl = false; // skip pattern with indel
    std::unordered_map<std::string, std::string> pam; // pattern annotation map
    std::vector<std::unordered_map<std::string, int64_t>> pcnt; // pattern count of each contig
    std::vector<int64_t> ttvr; // total valid seq
    std::vector<int64_t> ttpr; // total pattern count seq
    faidx_t* fai = NULL;
    int64_t** nuccnt = NULL;
    
    pscnt_opt_t();
    ~pscnt_opt_t();

    void init();
    bool valid();
    void fill_pam();
    void ana();
};

void pscnt_usage(pscnt_opt_t* opt, char* arg0);

int pscnt_main(int argc, char** argv);

#endif
