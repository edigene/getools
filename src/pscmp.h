#ifndef PSCMP_H
#define PSCMP_H

#include "pscnt.h"
#include "htmlopt.h"


struct contig_pscnt_qc_t{
    int64_t total_reads;
    int64_t with_anno_pat;
    int64_t with_other_pat;
    int64_t no_pat_found;

    void html2tbody(kstring_t* s, const std::string& n){
        ksprintf(s, "<tr>");
        ksprintf(s, "<td class='col1'>%s</td>", n.c_str());
        ksprintf(s, "<td class='col2'>%lld</td>", total_reads);
        ksprintf(s, "<td class='col2'>%lld</td>", with_anno_pat);
        ksprintf(s, "<td class='col2'>%lld</td>", with_other_pat);
        ksprintf(s, "<td class='col2'>%lld</td>", no_pat_found);
        ksprintf(s, "</tr>");
    }

    static void html2thead(kstring_t* s){
        ksprintf(s, "<thead><tr>");
        ksprintf(s, "<th class='col1'>Contig</th>");
        ksprintf(s, "<th class='col1'>TotalReads</th>");
        ksprintf(s, "<th class='col1'>WithAnnoPat</th>");
        ksprintf(s, "<th class='col1'>WithOtherPat</th>");
        ksprintf(s, "<th class='col1'>NoPaFoundt</th>");
        ksprintf(s, "</tr></thead>");
    }

    void tsv2body(kstring_t* s, const std::string& n){
        ksprintf(s, "%s\t%lld\t%lld\t%lld\t%lld\n",
                     n.c_str(),
                     total_reads,
                     with_anno_pat,
                     with_other_pat,
                     no_pat_found);
    }

    static void tsv2head(kstring_t* s){
        ksprintf(s, "Contig\tTotalReads\tWithAnnoPat\tWithOtherPat\tNoPatFound\n");
    }
};
typedef std::unordered_map<std::string, contig_pscnt_qc_t*> ContigPsCntQCMap;

// pattern count of case/ctrl
struct pscnt_info_t{
    int64_t count[2] = {0, 0};
    double freq[2] = {0, 0};
    bool val4logo = false;
    std::string anno;

    static void tsv2head(kstring_t* s){
        ksprintf(s, "Contig\tPattern\tCount\tFreq\tAnno\tCountInCtrl\tFreqInCtrl\tFoldChange\tLog2FC\n");
    }

    void tsv2body(kstring_t* s, const std::string& n, const std::string& seq){
        ksprintf(s, "%s\t%s\t%lld\t%lf\t%s\t%lld\t%lf\t%lf\t%lf\n",
                    n.c_str(),
                    seq.c_str(),
                    count[0],
                    freq[0],
                    anno.c_str(),
                    count[1],
                    freq[1],
                    freq[1] > 0 ? freq[0]/freq[1] : .0,
                    freq[1] > 0 ? log2(freq[0]/freq[1]) : .0
                    );
    }
};

typedef std::unordered_map<std::string, pscnt_info_t*> PsCntMap;
struct contig_pscnt_t{
    PsCntMap pscnt;
    int64_t ttpr = 0;
    bool isnorm = false;
};
typedef std::unordered_map<std::string, contig_pscnt_t*> ContigPsCntMap;

struct contig_nuc_freq_t{
    int len = 0;
    double** cfreq = NULL;
};
typedef std::unordered_map<std::string, contig_nuc_freq_t*> ContigPsNucFreqMap;

struct pscmp_opt_t{
    std::string splist; // sample list to do analysis(1 column name) 
    std::string indir; // output directory of pscnt(contain subdir of sample in splist)
    std::string cnorm; // normal contig name to use as control if any
    std::string outdir = "./"; // output directory
    std::string hrjsn; // jscdn
    std::string outhtml = "index.html"; // output html
    std::string n_outpsqc = "psqc.tsv";
    std::string outpsqc; // output psqc
    std::string n_outpscnt = "pscnt.tsv";
    std::string outpscnt; // output pscnt
    std::string n_outlogps = "pslogo.tsv";
    std::string outlogps; // output seq used to plot logo
    int plen = 10; // pattern length
    float maxfcl = -2.5; // max log fold change
    float fcymax = 0.4; // fold change y max
    float fcymin = -0.4; // fold change y min
    std::vector<std::string> samples; // samples
    std::vector<std::string> contigs; // contigs
    ContigPsCntQCMap map4ctgqc; // qc of each contig
    ContigPsCntMap map4pscnt; // pscnt of each contig
    ContigPsNucFreqMap map4psnuc; // nuc freq 
    HtmlOpt hmo; // html options
    kstring_t* ts = NULL; // temp kstring for tsv
    kstring_t* ps = NULL; // temp kstring for path
    double** nfreq = NULL; // control freq dist
    int64_t nttr = 0; // control total reads
    int64_t** vnuc = NULL; // valid nucleotides
    contig_pscnt_t* pcpnorm = NULL; // pointer to normal samples

    pscmp_opt_t();
    ~pscmp_opt_t();

    bool valid();
    void init();
    void parse_samples();
    void parse_in();
    void parse_qc();
    void parse_pscnt();
    void parse_freq();

    void genrpt();
    void html2head(kstring_t* s);
    void html2body(kstring_t* s);
    void html2foot(kstring_t* s);

    void allqc2html(kstring_t* s);
    void allps2html(kstring_t* s);
    void allnuc2html(kstring_t* s);
    void oneqc2html(kstring_t* s, const std::string& n);
    void oneps2html(kstring_t* s, const std::string& n);
    void onenuc2html(kstring_t* s, const std::string& n);
    void nuc2freq(kstring_t* s, const std::string& n, contig_nuc_freq_t* f);
    void nuc2logo4exist2(kstring_t* s, const std::string& n);
    void nuc2logo4exist1(kstring_t* s, const std::string& n);
};

void pscmp_usage(pscmp_opt_t* opt, char* arg0);
int pscmp_main(int argc, char** argv);

#endif
