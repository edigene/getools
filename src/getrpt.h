#ifndef GET_RPT_H
#define GET_RPT_H

#include "gtn.h"
#include "util.h"
#include "common.h"
#include "grpmut.h"
#include "htmlopt.h"
#include "htslib/kstring.h"

#ifndef MAX_MUTMM_DIST_LEN
#define MAX_MUTMM_DIST_LEN 200
#endif

struct rpt_dir_t{
    std::string consdir = "GE.Consistence";
    std::string qcdir = "QC.Pooled";
    std::string qcfile = "pooled.qc.html";
    std::string sumdir = "GE.Summary";
    std::string rptdir = "GE.Report";
    std::string hapdir = "GE.HapCNT";
    std::string snvdir = "GE.SNVCNT";
    // html
    std::string sbiasfile = "snvbias.summary.html";
    std::string hapfile = "hapcnt.summary.html";
    std::string snvfile = "snvcnt.summary.html";
    std::string aacfile = "aaccnt.summary.html";
    std::string snvfocf = "snvfoc.summary.html";
    std::string mutfreq = "mutfreq.summary.html";
    std::string mutn8rdist = "mutn8rdist.summary.html";
    std::string edfile = "pooled.ed.html";
    // tsv
    std::string tsv_sbiasfile = "snvbias.summary.tsv";
    std::string tsv_qcfile = "pooled.qc.tsv";
    std::string tsv_hapfile = "hapcnt.summary.tsv";
    std::string tsv_snvfile = "snvcnt.summary.tsv";
    std::string tsv_aacfile = "aaccnt.summary.tsv";
    std::string tsv_snvfocf = "snvfoc.summary.tsv";
    std::string tsv_mutfreq = "mutfreq.summary.tsv";
    std::string tsv_mutn8rdist = "mutn8rdist.summary.tsv";
    std::string tsv_edfile = "pooled.ed.tsv";
};

struct rpt_opt_t{
    char* indir = NULL; // getools analysis output dir of a series samples
    char* splist = NULL; // 1 column sample name list to generate report
    std::string spstr; // sapmle names, seperated by ','
    char* outdir = NULL; // output directory of report files
    char* title = NULL; // output title
    std::string hapindir; // hap indir
    bool haprpt = false; // hap report ?
    std::string ampgrpf; // amplicon group info tsv(col1:amp_name, col2:amp_grp)
    std::unordered_map<std::string, std::string> amp2grp; // amplicon to group name
    std::string extsnv; // extra snv to focus on
    GroupMutColType amnv; // all mnvs
    std::vector<std::string> nmnv; // mnv name
    std::vector<std::string> outmsnv; // output mnv file basename
    GroupMutColType afreq; // all freqs
    std::vector<std::string> nfreq; // mfreq name
    std::vector<std::string> outmfreq; // output freq file basename
    GroupMutColType arnct; // all dist
    std::vector<std::string> nrcnt; // rcnt name
    std::vector<std::string> outrcnt; // output rcnt dist file basename
    GroupMutColType abias; // all bias
    std::vector<std::string> nbias; // bias name
    std::vector<std::string> outbias; // output bias dist file basename
    size_t maxp4focsnv = 0; // max pos count of nuc
    size_t maxp4focfreq = 0; // max pos count of freq
    size_t maxp4focrcnt = 0; // max pos coutn of dist
    size_t maxp4bias = 25; // max pos countn of bias
    size_t maxmdistl = 10; // max nucleotide mismatch to cal from inframe region
    int nrow = 0; // table row number
    bool gen4web = false; // generate for index.html located at indir
    bool diffamps = false; // amplicon configure file are different for each samples if set
    char*  consamp = NULL; // consistence analysis amplicons 
    int na4cons = 50; // consistence analysis alleles
    int topnW = 800; // topN allele graph width
    bool uniqa = false; // consistence analysis for all unique alleles
    bool showt = false;
    bool calg = false;
    int drmethod = 0;
    int ncol = 6; // number of samples per row
    std::string jscfg;
    std::vector<std::string> samples; // samples
    kstring_t* ks = NULL; // html string
    kstring_t* ps = NULL; // temp string for path
    kstring_t* ss = NULL; // temp string for string 
    kstring_t* ts = NULL; // temp string for string of tsv
    HtmlOpt hmo;
    rpt_dir_t rptd;

    rpt_opt_t();
    ~rpt_opt_t();

    bool valid();
    void init();

    void reportHTMLHeader();
    void reportHTMLBody();
    void reportHTML();
    void writeHTMLTableHeader(kstring_t* s, std::string ttitle);
    void tsv2HTMLTable(char* f, kstring_t* s);

    void gen_consf();
    void gen_sumqcf();
    void gen_sumedf();
    void gen_qcf();
    void gen_hap();
    void gen_snv();
    void gen_aac();
    void gen_fsnv();
    void gen_fmutf();
    void gen_fmndist();
    void gen_biasf();
    void cp_rpf();
    void reportFiles();

    void report();
};

void rpt_usage(rpt_opt_t* opt, char* arg0);
int rpt_main(int argc, char** argv);

#endif
