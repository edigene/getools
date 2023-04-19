#ifndef AMPLICON_H
#define AMPLICON_H

#include <string>
#include <vector>

// amplicon info
struct amplicon_t{
    std::string aname; // amplicon name
    std::string aseq; // amplicon sequence
    std::string sgrseq; // sgr sequence [sgRNA_PAM forward strand sequence]
    std::string donor; // donor sequence used for recombination insertion
    int dbeg; // donor align beg pos [0 base inclusive]
    int denda; // donor align end pos on aseq[0 base inclusive]
    int dendd; // donor align end pos on donor[0 base inclusive]
    std::vector<std::string> dins; // donor ins
    std::vector<std::string> ddis; // donor delins
    std::vector<int32_t> ddel; // donor del
    std::vector<int32_t> dsnv; // donor snv
    int dmutcnt; // donor total mutations
    int sgrbeg; // sgr begin[0based, included]
    int sgrend; // sgr end[0based, included]
    int clsbeg; // cleavage site begin[0based, included]
    int clsend; // cleavage site end[0based, included]
    int maxdel; // max del allowed
    int fplen; // fwd primer length
    int rvlen; // rev primer length
    int maxrpl; // max single nucleotide repeat length

    amplicon_t(){}

    ~amplicon_t(){}
};

typedef std::vector<amplicon_t*> AmpliconList;

#endif
