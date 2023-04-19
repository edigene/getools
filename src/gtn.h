#ifndef GET_GTNC_H
#define GET_GTNC_H

#include <stdio.h>
#include <libgen.h>
#include <unistd.h>
#include <map>
#include "util.h"
#include "common.h"
#include "htmlopt.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/kstring.h"
#include <unordered_set>


static const char* GEV_STR_ARR[10] = {
    "REF", "SNV", "INS", "-", "DEL",
    "-"  , "-"  , "-"  , "DIN", "Others"};
    
// variant type
struct variant_t{
    std::string ref; // ref
    std::string alt; // alt
    int32_t start; // start
    int32_t end; // end
    double freq; // freq
    double sfreq; // shared average freq across samples
    int sidx; // sample index
    int width; // width
    std::string type; // type
    int ntype; // type ints
    int count; // count
    int varid; // variant id, same variant share same id
    int isgtn; // is topN

    variant_t(){
        ref = "REF";
        alt = "ALT";
        start = 0;
        end = 0;
        freq = .0;
        sidx = -1;
        width = 0;
        type = GEV_STR_ARR[GEVAR_OTH];
        ntype = GEVAR_OTH;
        count = 0;
        varid = -1;
    }

    ~variant_t(){}

    bool same_var(const variant_t* v) const {
        return ref == v->ref && alt == v->alt && start == v->start && end == v->end && ntype == v->ntype;
    }

    bool dup_var(variant_t* v){
        return sidx == v->sidx && ref == v->ref && alt == v->alt && start == v->start && end == v->end && ntype == v->ntype;
    }

    static void out_head(FILE* f){
        fprintf(f, "Sample\tVarID\tFreq\tCount\tF\tStart\tEnd\tType\tREF\tALT\n");
    }

    void out_rec(FILE* f, const std::vector<std::string>& sps){
        fprintf(f, "%s\t%d\t%f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", sps[sidx].c_str(), varid, freq, count, width, start, end, type.c_str(), ref.c_str(), alt.c_str());
    }
};

// variant sorter
struct VarSortBySpTp{
    bool operator()(const variant_t* v1, const variant_t* v2) const {
        return (v1->sidx > v2->sidx) ||
               (v1->sidx == v2->sidx && v1->ref > v2->ref) ||
               (v1->sidx == v2->sidx && v1->ref == v2->ref && v1->alt > v2->alt) ||
               (v1->sidx == v2->sidx && v1->ref == v2->ref && v1->alt == v2->alt && v1->start > v2->start) ||
               (v1->sidx == v2->sidx && v1->ref == v2->ref && v1->alt == v2->alt && v1->start == v2->start && v1->end > v2->end) ||
               (v1->sidx == v2->sidx && v1->ref == v2->ref && v1->alt == v2->alt && v1->start == v2->start && v1->end == v2->end && v1->ntype > v2->ntype) ||
               (v1->sidx == v2->sidx && v1->ref == v2->ref && v1->alt == v2->alt && v1->start == v2->start && v1->end == v2->end && v1->ntype == v2->ntype && v1->freq > v2->freq);
    }
};

// variant sorter
struct VarSortByType{
    bool operator()(const variant_t* v1, const variant_t* v2) const {
        return (v1->ref > v2->ref) ||
               (v1->ref == v2->ref && v1->alt > v2->alt) ||
               (v1->ref == v2->ref && v1->alt == v2->alt && v1->start > v2->start) ||
               (v1->ref == v2->ref && v1->alt == v2->alt && v1->start == v2->start && v1->end > v2->end) ||
               (v1->ref == v2->ref && v1->alt == v2->alt && v1->start == v2->start && v1->end == v2->end && v1->ntype > v2->ntype) ||
               (v1->ref == v2->ref && v1->alt == v2->alt && v1->start == v2->start && v1->end == v2->end && v1->ntype == v2->ntype && v1->freq > v2->freq);
    }
};

// variant sorter
struct VarSortByFreq{
    bool operator()(const variant_t* v1, const variant_t* v2) const {
        return v1->freq > v2->freq;
    }
};

// variant sorter
struct VarSortByShareFreq{
    bool operator()(const variant_t* v1, const variant_t* v2) const {
        return v1->sfreq > v2->sfreq;
    }
};

struct VarSortByShareFreqAndType{
    bool operator()(const variant_t* v1, const variant_t* v2) const {
        return (v1->sfreq > v2->sfreq) ||
               (v1->sfreq == v2->sfreq && v1->varid > v2->varid);
    }
};

// getools caledit edit details schema
struct get_vt_schema{
    int freq = 0;
    int count = 1;
    int width = 2;
    int start = 3;
    int end = 4;
    int type = 5;
    int ref = 6;
    int alt = 7;
    int tmark = 8;
};

// options 
struct gtn_opt_t{
    char* inlist = NULL;
    char* outpre = NULL;
    char* colsys = NULL;
    std::string hrjsn;
    std::map<std::string, std::string> colmap;
    std::string colunknown;
    std::string colother;
    std::string unknownkey = "unknown";
    std::string otherkey = "other";
    std::string sortn; // sort by name
    // bcf input usage beg
    char* amp = NULL;
    int beg = -1;
    int end = -1;
    int nctc = 3;
    float mrv = 0.5; // min ref freq needed to consider as ref
    float mxv = 0.3; // min ref freq needed to consider as Wt/Mt
    // bcf input usage end
    int type = 6;
    int topn = 50;
    bool uniqa = false;
    bool calg = false;
    int width = 12;
    int height = 8;
    int hw = 800; // html width
    float lr = 0.8;
    float rr = 0.2;
    get_vt_schema gs;
    HtmlOpt hmo;
    // variables
    std::vector<std::string> samples; // sample names
    std::vector<std::string> grps; // group names
    std::set<std::string> sgrps; // group set
    std::vector<std::string> sid2grp; // sample index to group names
    std::vector<std::string> sidngrp; // sample index to group number
    std::vector<int> grns; // sample count in each group
    std::vector<variant_t*> vars; // all variants from bcf in ranges
    std::vector<int> sharen; // share count
    std::vector<double> sharef; // share vaf
    std::vector<double> sumaf; // summed vaf
    std::vector<double> otfreq; // other acc freq
    double maxaf = .0; // max acc af
    int lastid = 0; // id mark same var across samples
    // edeff
    std::vector<double> edieff, recpef, receff, reeeff, muteff;
    std::vector<int64_t> rawcnt, totcnt;
    std::vector<int> recstat, gttype;
    std::vector<std::vector<double>> gtnacc;
    int needrec = 0;
    int nsamples = 0; // pure samples [+ groups]
    int psamples = 0; // pure samples
    int nallele = 0;
    bool showt = false; // show sample name in DR
    int drmethod = 0; // dimension reduction method
    double** pa4s = NULL;
    std::unordered_map<int, std::string> s2gm; // sample idx to group map
    std::vector<std::string> uniqs; // uniq group


    gtn_opt_t(){}
    ~gtn_opt_t(){
        for(auto& v: vars){
            delete v;
        }
        if(pa4s){
            for(int i = 0; i < nsamples; ++i){
                free(pa4s[i]);
            }
            free(pa4s);
            pa4s = NULL;
        }
    }

    int cons_ana();
    int parse_bcfs();
    int call_gt(std::vector<variant_t*>& vs, double ef);
    void ana8R();
    void ana8js();
    void html2head(kstring_t* s);
    void html2body(kstring_t* s);
    void html2foot(kstring_t* s);
};

void gtn_usage(gtn_opt_t* opt, char* arg0);

int gtn_main(int argc, char** argv);

#endif
