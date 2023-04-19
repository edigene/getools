#ifndef GETOOLS_HAPCNT_H
#define GETOOLS_HAPCNT_H

#include <string>
#include "util.h"
#include <unistd.h>
#include <libgen.h>
#include "common.h"
#include "htmlopt.h"
#include "txtrs.h"
#include "ksw4get.h"
#include <unordered_map>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

typedef uint32_t TYPE_AAS;

#ifndef AAS_TYPE_SYNONYM
#define AAS_TYPE_SYNONYM 0x1
#endif

#ifndef AAS_TYPE_NONSYNO
#define AAS_TYPE_NONSYNO 0x2
#endif

#ifndef AAS_TYPE_WITHALT
#define AAS_TYPE_WITHALT 0x4
#endif

#ifndef AAS_TYPE_ONLYALT
#define AAS_TYPE_ONLYALT 0x8
#endif

#ifndef NUC_WITHALT
#define NUC_WITHALT 0x10
#endif

typedef std::unordered_map<std::string, TYPE_AAS> AAS2TypeMap;

struct bam_rec_t{
    int32_t cnt;
    int32_t ma2a;
    int32_t mn2a;
    int32_t ma2r;
    int32_t mn2r;
    std::string muts;
    std::string refs;
    std::string seq;
    std::string aas;
    TYPE_AAS aat = 0;
    std::string oriaa;
    std::string mutaa;
    std::vector<int32_t> pos;
    std::vector<int32_t> apos;

    static void bam2head(std::ofstream& ofs){
        ofs << "ReadSeq\tReadCount\tRate2All\tRate2Muts\tRate2HapAll\tRate2HapMuts\tHapLength\tRefType\tHapType\tHapPos\n";
    }

    static void aas2head(std::ofstream& ofs){
        ofs << "AASeq\tReadCount\tRate2All\tRate2Muts\tRate2HapAll\tRate2HapMuts\tMutAACnt\tOrigType\tMutType\tMutPos\tAnno\n";
    }

    void bam2rec(std::ofstream& ofs, const int64_t& ttr, const int64_t& tar, const int64_t& tcr, const int64_t& tcm){
        ofs << seq << "\t";
        ofs << cnt << "\t";
        ofs << (double)cnt/(double)ttr << "\t";
        ofs << (double)cnt/(double)tar << "\t";
        ofs << (double)cnt/(double)tcr << "\t";
        ofs << (double)cnt/(double)tcm << "\t";
        ofs << pos.size() << "\t" << refs << "\t" << muts << "\t";
        for(auto& p: pos){
            ofs << p+1 << ",";
        }
        ofs << "\n";
    }

    void aas2rec(std::ofstream& ofs, const int64_t& ttr, const int64_t& tar, const int64_t& tcr, const int64_t& tcm){
        ofs << aas << "\t";
        ofs << cnt << "\t";
        ofs << (double)cnt/(double)ttr << "\t";
        ofs << (double)cnt/(double)tar << "\t";
        ofs << (double)cnt/(double)tcr << "\t";
        ofs << (double)cnt/(double)tcm << "\t";
        ofs << apos.size() << "\t" << oriaa << "\t" << mutaa << "\t";
        for(auto& p: apos){
            ofs << p+1 << ",";
        }
        ofs << "\t";
        if(aat & AAS_TYPE_SYNONYM){
            ofs << "Synonym,";
        }else if(aat & AAS_TYPE_ONLYALT){
            ofs << "OnlyDesiredMut,";
        }else if(aat & AAS_TYPE_WITHALT){
            ofs << "WithDesiredMut,";
        }else{
            ofs << "OtherMut,";
        }
        ofs << "\n";
    }
};

struct bam_rec_cnt_sorter{
    bool operator()(const bam_rec_t* r1, const bam_rec_t* r2) const {
        return (r1->cnt > r2->cnt) ||
               (r1->cnt == r2->cnt && r1->pos.size() > r2->pos.size());
    }
};

struct bam_rec_seq_sorter{
    bool operator()(const bam_rec_t* r1, const bam_rec_t* r2) const {
        return (r1->seq < r2->seq) ||
               (r1->seq == r2->seq  && r1->cnt > r2->cnt);
    }
};

struct bam_rec_aas_sorter{
    bool operator()(const bam_rec_t* r1, const bam_rec_t* r2) const {
        return (r1->aas < r2->aas) ||
               (r1->aas == r2->aas  && r1->cnt > r2->cnt);
    }
};

struct hap_rec_t{
    int32_t cnt;
    std::string muts;
    std::string refs;
    std::vector<int32_t> pos;

    static void hap2head(std::ofstream& ofs){
        ofs << "HapSupportReads\tHapRate2All\tHapRate2Muts\tHapRate2HapAll\tHapRate2HapMuts\tHapLength\tRefType\tHapType\tHapPos\n";
    }

    void hap2rec(std::ofstream& ofs, const int64_t& ttr, const int64_t& tmr, const int64_t& tcr, const int64_t& tcm){
        ofs << cnt << "\t";
        ofs << (double)cnt/(double)ttr << "\t";
        ofs << (double)cnt/(double)tmr << "\t";
        ofs << (double)cnt/(double)tcr << "\t";
        ofs << (double)cnt/(double)tcm << "\t";
        ofs << pos.size() << "\t" << refs << "\t" << muts << "\t";
        for(auto& p: pos){
            ofs << p+1 << ",";
        }
        ofs << "\n";
    }
};

struct hap_rec_sorter{
    bool operator()(const hap_rec_t* r1, const hap_rec_t* r2) const {
        return (r1->cnt > r2->cnt) ||
               (r1->cnt == r2->cnt && r1->pos.size() > r2->pos.size());
    }
};

typedef std::vector<hap_rec_t*> HapRecList;
typedef std::vector<std::vector<int>> RecombBitMask;
typedef std::vector<RecombBitMask> RecombBitMaskList;
typedef std::vector<RecombBitMaskList> RecombBitMaskSet;
typedef std::unordered_map<std::string, hap_rec_t*> HapMap;
typedef std::vector<bam_rec_t*> BamRecList;

void merge_and_sort_bam_rec8seq(BamRecList& brl);
void merge_and_sort_bam_rec8aas(BamRecList& brl);

struct hap_opt_t{
    std::string inbam; // input bam
    std::string contig; // contig
    std::string extsnv; // extra cnv to focus
    int tid; // tid
    int beg; // beg pos, 1 based inclusive
    int end; // end pos, 1 based inclusive
    int subl; // sub length to get
    int sual; // sub aa length
    std::string ref; // ref fasta file
    std::string alt; // alt fasta file
    int64_t ttr = 0; // total reads
    int64_t tmr = 0; // total reads with desired snp mutations
    int64_t tar = 0; // total reads with any mutations(different to reference type)
    int64_t tfr = 0; // total reads with indel/soft clip in [beg, end] query region
    int64_t tnr = 0; // total reads with non snp in [beg, end] query region
    int64_t txr = 0; // total reads with too many snp in [beg, end] query region
    int64_t trr = 0; // total reference type reads
    int64_t tcr = 0; // total cleaned reads(ttr-tfr-txr)
    int64_t tat = 0; // total animo acid seq types in inframe region
    bool ismt = false; // is mitochondrion
    bool nget = false; // bam not from getools
    int maxh = 10; // max linked snp hapcount
    int minr = 3; // min support reads to output
    int otn = 1000; // top n event to output to table
    int gtn = 100; // top n event to plot
    int nthn = 1000; // haha
    int refl = 0; // ref length
    bool l2s = false; // long hap accumulate short hap if true
    size_t umaxh = 10; // size_t type max hapcount
    KSW_FTYPE dropmask = 0; // drop mask
    std::string outdir = "./out"; // output dir
    std::string hapb2tsv = HAPBIAS2_OUTTSV; // 2base mute bias
    std::string hapb1tsv = HAPBIAS1_OUTTSV; // 1base mute bias
    std::string outtsv = HAPCNT_OUTTTSV; // output tsv
    std::string outhtml = "hapcnt.html"; // output html
    std::string outjsn = "hapcnt.json"; // output json
    std::string outstsv = HAPCNT_OUTSTSV; // stat info of all reads
    std::string outntsv = "hapseqn.tsv"; // output topN bam record
    std::string outatsv = "hapaasn.tsv"; // output topN aas record
    std::string outsnv = HAPCNT_OUTSNV; // output snv info in whole range
    std::string outaac = HAPCNT_OUTAAC; // output aac info in inframe range
    std::string outmfaac = HAPCNT_OUTMFAAC; // output aac mut freq in inframe range
    std::string outmfnuc = HAPCNT_OUTMFNUC; // output nuc mut freq in inframe range
    std::string outrcntnall2ref = HAPCNT_OUTRCNTNALL2REF; // all nuc mut2ref count dist by read in inframe range
    std::string outrcntnall2alt = HAPCNT_OUTRCNTNALL2ALT; // all nuc mut2alt count dist by read in inframe range
    std::string outrcntaall2ref = HAPCNT_OUTRCNTAALL2REF; // all acc mut2ref count dist by read in inframe range
    std::string outrcntaall2alt = HAPCNT_OUTRCNTAALL2ALT; // all acc mut2alt count dist by read in inframe range
    std::string outrcntnndef2ref = HAPCNT_OUTRCNTNNDEF2REF; // all nuc mut2ref count dist by read in inframe range(no desired mut)
    std::string outrcntnndef2alt = HAPCNT_OUTRCNTNNDEF2ALT; // all nuc mut2alt count dist by read in inframe range(no desired mut)
    std::string outrcntandef2ref = HAPCNT_OUTRCNTANDEF2REF; // all acc mut2ref count dist by read in inframe range(no desired mut)
    std::string outrcntandef2alt = HAPCNT_OUTRCNTANDEF2ALT; // all acc mut2alt count dist by read in inframe range(no desired mut)
    std::string outrcntnydef2ref = HAPCNT_OUTRCNTNYDEF2REF; // all nuc mut2ref count dist by read in inframe range(with desired mut)
    std::string outrcntnydef2alt = HAPCNT_OUTRCNTNYDEF2ALT; // all nuc mut2alt count dist by read in inframe range(with desired mut)
    std::string outrcntaydef2ref = HAPCNT_OUTRCNTAYDEF2REF; // all acc mut2ref count dist by read in inframe range(with desired mut)
    std::string outrcntaydef2alt = HAPCNT_OUTRCNTAYDEF2ALT; // all acc mut2alt count dist by read in inframe range(with desired mut)
    std::vector<std::string> outmsnv; // other snv info stat
    std::vector<kstring_t*> outmstr; // other snv info stat
    std::vector<kstring_t*> outpstr; // acgt pos str
    std::string jscdn; // jscdn
    HapMap hapm; // hap count map
    HapRecList hapl; // hap list
    HtmlOpt hmo; // html opt
    faidx_t* fai = NULL;
    char* far = NULL; // ref seq
    char* faa = NULL; // alt seq
    std::vector<int32_t> aps; // mut pos
    std::vector<uint8_t> mnv; // mut snv
    std::vector<uint8_t> mmv; // mut mark
    int mc8nuc = 0;
    int mc8aac = 0;
    BamRecList brl; 
    RecombBitMaskSet rbms;
    kstring_t* ks = NULL;
    kstring_t* ms = NULL;
    kstring_t* rs = NULL;
    char* bseq = NULL;
    std::vector<int32_t> bpos;
    std::vector<int32_t> tpos;
    AAS2TypeMap aas2type; // aas to type map
    std::string xaas, xras, refas, altas;
    std::vector<int> xapos;
    std::vector<char> xamark;
    std::vector<int> refais;
    int allaltc = 0;
    int64_t xasyn = 0;
    int64_t xaonlym = 0;
    int64_t xawithm = 0;
    int64_t xaother = 0;
    uint8_t* bnucref = NULL;
    int64_t** bnuccnt = NULL;
    int64_t** snuccnt = NULL;
    int64_t** bnucttt = NULL;
    int64_t** nuccnt = NULL;
    int64_t** nucttt = NULL;
    int64_t** aaccnt = NULL;
    int64_t* ttdepn = NULL;
    int64_t* ttdepa = NULL;
    int64_t* ttmutn = NULL;
    int64_t* ttmuta = NULL;
    int64_t* mutn2rcntall2ref = NULL;
    int64_t* muta2rcntall2ref = NULL;
    int64_t* mutn2rcntndef2ref = NULL;
    int64_t* muta2rcntndef2ref = NULL;
    int64_t* mutn2rcntydef2ref = NULL;
    int64_t* muta2rcntydef2ref = NULL;
    int64_t* mutn2rcntall2alt = NULL;
    int64_t* muta2rcntall2alt = NULL;
    int64_t* mutn2rcntndef2alt = NULL;
    int64_t* muta2rcntndef2alt = NULL;
    int64_t* mutn2rcntydef2alt = NULL;
    int64_t* muta2rcntydef2alt = NULL;
    uint8_t* refints = NULL;
    uint8_t* altints = NULL;
    int8_t cntdep = 0;

    hap_opt_t(){}
    ~hap_opt_t(){
        if(fai) { fai_destroy(fai); fai = NULL; }
        if(far) { free(far); far = NULL; }
        if(faa) { free(faa); faa = NULL; }
        if(bseq) { free(bseq); bseq = NULL; }
        if(nuccnt) {
            for(int i = 0; i < 5; ++i){
                if(nuccnt[i]){ free(nuccnt[i]); nuccnt[i] = NULL; }
            }
            free(nuccnt);
            nuccnt = NULL;
        }
        if(bnucref){
            free(bnucref); bnucref = NULL;
        }
        if(bnuccnt) {
            for(int i = 0; i < 37; ++i){
                if(bnuccnt[i]){ free(bnuccnt[i]); bnuccnt[i] = NULL; }
            }
            free(bnuccnt);
            bnuccnt = NULL;
        }
        if(snuccnt) {
            for(int i = 0; i < 5; ++i){
                if(snuccnt[i]){ free(snuccnt[i]); snuccnt[i] = NULL; }
            }
            free(snuccnt);
            snuccnt = NULL;
        }
        if(nucttt) {
            for(int i = 0; i < 5; ++i){
                if(nucttt[i]){ free(nucttt[i]); nucttt[i] = NULL; }
            }
            free(nucttt);
            nucttt = NULL;
        }
        if(bnucttt) {
            for(int i = 0; i < 37; ++i){
                if(bnucttt[i]){ free(bnucttt[i]); bnucttt[i] = NULL; }
            }
            free(bnucttt);
            bnucttt = NULL;
        }
        if(aaccnt){
            for(int i = 0; i < ALL_AA_CNT; ++i){
                if(aaccnt[i]){ free(aaccnt[i]); aaccnt[i] = NULL; }
            }
            free(aaccnt);
            aaccnt = NULL;
        }
        if(ttdepn){ free(ttdepn); ttdepn = NULL; }
        if(ttdepa){ free(ttdepa); ttdepa = NULL; }
        if(ttmutn){ free(ttmutn); ttmutn = NULL; }
        if(ttmuta){ free(ttmuta); ttmuta = NULL; }
        if(mutn2rcntall2ref){ free(mutn2rcntall2ref); mutn2rcntall2ref = NULL; }
        if(muta2rcntall2ref){ free(muta2rcntall2ref); muta2rcntall2ref = NULL; }
        if(mutn2rcntndef2ref){ free(mutn2rcntndef2ref); mutn2rcntndef2ref = NULL; }
        if(muta2rcntndef2ref){ free(muta2rcntndef2ref); muta2rcntndef2ref = NULL; }
        if(mutn2rcntydef2ref){ free(mutn2rcntydef2ref); mutn2rcntydef2ref = NULL; }
        if(muta2rcntydef2ref){ free(muta2rcntydef2ref); muta2rcntydef2ref = NULL; }
        if(mutn2rcntall2alt){ free(mutn2rcntall2alt); mutn2rcntall2alt = NULL; }
        if(muta2rcntall2alt){ free(muta2rcntall2alt); muta2rcntall2alt = NULL; }
        if(mutn2rcntndef2alt){ free(mutn2rcntndef2alt); mutn2rcntndef2alt = NULL; }
        if(muta2rcntndef2alt){ free(muta2rcntndef2alt); muta2rcntndef2alt = NULL; }
        if(mutn2rcntydef2alt){ free(mutn2rcntydef2alt); mutn2rcntydef2alt = NULL; }
        if(muta2rcntydef2alt){ free(muta2rcntydef2alt); muta2rcntydef2alt = NULL; }
        if(refints){ free(refints); refints = NULL; }
        if(altints){ free(altints); altints = NULL; }
        for(auto& e: outmstr){
            if(e){
                if(e->s) free(e->s);
                free(e);
                e = NULL;
            }
        }
        for(auto& e: outpstr){
            if(e){
                if(e->s) free(e->s);
                free(e);
                e = NULL;
            }
        }
        if(ks){
            if(ks->s) free(ks->s);
            free(ks);
            ks = NULL;
        }
        if(ms){
            if(ms->s) free(ms->s);
            free(ms);
            ms = NULL;
        }
        if(rs){
            if(rs->s) free(rs->s);
            free(rs);
            rs = NULL;
        }
        for(auto& e: hapl){
            delete e;
            e = NULL;
        }
        for(auto& e: hapm){
            e.second = NULL;
        }
        hapl.clear();
        hapm.clear();
    }

    bool valid();
    void update();
    void allocmem();
    void buildam();
    void ref2bnuc();
    void ref2altp();
    void naa2dep();
    void foc2mnv();
    void snv2foc();
    void hapcnt8bv(std::vector<bam1_t*>& inbl, int blen);
    void hapcnt8bf();
    void hapcnt1b(bam1_t* b);
    void annoaas();
    void hap2tsv();
    void stat2tsv();
    void bam2tsv();
    void aas2tsv();
    void hap2html();
    void hap2jsn();
    void snv2tsv();
    void aac2tsv();
    void bias2tsv();
    void bm2ex(bam1_t* b);
    void mc2tsv(int64_t* vec, int32_t len, int64_t sum, const std::string& outf); // len is lenof(vec)-1
    void printHeader(kstring_t* s);
    void printBody(kstring_t* s);
    void printFooter(kstring_t* s);
};

void hapc_usage(hap_opt_t* opt, char* arg0);
int hapc_main(int argc, char** argv);
void out_rbms(RecombBitMaskSet& rbms);

#endif
