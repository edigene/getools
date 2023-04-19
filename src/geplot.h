#ifndef GENE_EDIT_PLOT_H
#define GENE_EDIT_PLOT_H

#include "ksw4get.h"
#include "util.h"
#include "txtrs.h"
#include <zlib.h>
#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include "common.h"
#include "htslib/sam.h"
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

// html var record
struct html_var_t{
    kstring_t* seq = NULL; // sequence
    int fscnt = 0; // frame shift
    kstring_t* cigar = NULL; // alignment cigar
    int cc = 0; // cluster count
    double af = .0; // af
    int cn = 0; // cell count
    std::vector<kstring_t*> iseq; // ins seq
    std::vector<int32_t> ipos; // ins pos
};

// variants
struct var_rec_t{
    bam1_t* b = NULL; // bam record
    int cc = 0; // cluster count
    double af = .0; // af
    int cn = 0; // cell count
    KSW_FTYPE vs = 0; // variant status
    int ins = 0; // ins length spanning sgrna±flklen
    int del = 0; // del length spanning sgrna±flklen
    char* seq = NULL; // sequence
    int lseq = 0; // sequence length
    std::vector<int> ipos; // insertion positions
    int lpos = 0; // left most pos of seq relative to sgrna beg
};

// variants sorter
struct var_sort_t{
    inline bool operator()(const var_rec_t* v1, const var_rec_t* v2) const {
        return v1->cc > v2->cc;
    }
};


// gene edit plot options
struct gep_opt_t{
    const char* inbam; // input bam
    const char* config; // getools configure file
    const char* amplicon; // amplicon sequence
    const char* name; // smplicon name
    const char* outdir; // output parent dir
    char* outpdf; // output pdf
    char* outtsv; // output tsv
    char* outbam; // output top n BAM
    char* outrsc; // output R script
    bam_hdr_t* bamhdr; // bam header
    int rlen; // amplicon len
    int topn; // top n to plot
    int sgrbeg; // sgrna beg pos on amplicon(0based)
    int sgrend; // sgrna end pos on amplicon(0based)
    int clsbpos; // cleavage site beg pos, relative to sgRNA_PAM end(1based)
    int clsepos; // cleavage site end pos, relative to sgRNA_PAM end(1based)
    int cutlen; // cut buffer length around sgrna to call variant
    int clsbeg; // cleavage site end pos(0based)
    int clsend; // cleavage site beg pos(0based)
    int vbeg; // variant beg pos considered as edited event(0based)
    int vend; // variant end pos considered as edited event(0based)
    int flklen; // flank length to plot var
    int refbeg; // beg pos to plot on amplicon(0based)
    int refend; // end pos to plot on amplicon(0based)
    int refgot; // total length to plot (refend - refbeg + 1)
    int hrbeg; // html ref beg
    int hrend; // html ref end
    int hfigh; // html figure height
    int64_t allrs; // all rs
    int64_t fscnt; // frame shift count
    int64_t edcnt; // edited count
    int64_t macnt; // match count
    KSW_FTYPE topm; // topn mask
    KSW_FTYPE dropm; // drop mask
    bool ismt; // is mitochondrion
    bool empty; // is variants to plot empty
    bool genhtvar; // generate html plot vars
    std::vector<html_var_t*> hvars; // html plot usage

    gep_opt_t(){
        inbam = config = amplicon = name = NULL;
        outdir = outpdf = outtsv = outrsc = outbam = NULL;
        bamhdr = NULL;
        flklen = 25;
        cutlen = 5;
        clsbpos = -1;
        clsepos = -1;
        topn = 50;
        topm = (KSW_FMAYVARSEQ | KSW_FREFTYPE);
        dropm = (KSW_FHIERRLOWQ | KSW_FLOWFREQ | KSW_FPRIMDIMER | KSW_FMISPRIMINN | KSW_FMISPRIMOUT | KSW_FLOWSURPT | KSW_FMANYVARS);
        fscnt = edcnt = macnt = 0;
        ismt = false;
        empty = false;
        genhtvar = false;
    }

    ~gep_opt_t(){
        if(outpdf){ free(outpdf); outpdf = NULL; }
        if(outtsv){ free(outtsv); outtsv = NULL; }
        if(outrsc){ free(outrsc); outrsc = NULL; }
        if(outbam){ free(outbam); outbam = NULL; }
        if(bamhdr){ sam_hdr_destroy(bamhdr); bamhdr = NULL; }
        if(genhtvar){
            for(auto& e: hvars){
                if(e){
                    if(e->seq){ 
                        if(e->seq->s){ free(e->seq->s); e->seq->s = NULL; }
                        free(e->seq); e->seq = NULL;
                    }
                    if(e->cigar){
                        if(e->cigar->s){ free(e->cigar->s); e->cigar->s = NULL; }
                        free(e->cigar); e->cigar = NULL;
                    }
                    free(e);
                    e = NULL;
                }
            }
        }
    }

    bool valid_opt();
    void parse_cfg();
    void update_iof();
    void update_amp();
    void plot();
    void bam2hv(bam1_t* b, html_var_t* v);
    void plot(const std::vector<bam1_t*>& vbs);
    void plot(const std::vector<var_rec_t*>& vrs);
    void topn2html(kstring_t* s);
};

extern int geplot_main(int argc, char** argv);

#endif
