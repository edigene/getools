#ifndef AMP_PREFIX_H
#define AMP_PREFIX_H

#include <map>
#include <string>
#include <unordered_map>
#include "util.h"

// prefix info used to split reads
struct prefix_t{
    std::string fseq; // forward seq
    std::string rseq; // reverse seq
    uint64_t fhash; // forward hash
    uint64_t rhash; // reverse hash
    bool fwd; // fseq is from forward strand of templates if true
    int sidx; // sample index this prefix belongs to
    bool seidx = false; // sample only constructed with single index
    int fbo = 0; // fseq barcode/index length to skip
    int rbo = 0; // rseq barcode/index length to skip

    // constructor1
    prefix_t(){}

    // constructor2
    prefix_t(const std::string& f, const std::string& r, char s, int si){
        fseq = f;
        rseq = r;
        fwd = s == '+' ? true : false;
        sidx = si;
        fhash = rhash = 0;
        if(fseq.length() < 32 && rseq.length() < 32){
            for(size_t i = 0; i < fseq.length(); ++i){
                fhash <<= 2UL;
                fhash |= nuc_to_2bit[(int)fseq[i]];
            }
            for(size_t i = 0; i < rseq.length(); ++i){
                rhash <<= 2UL;
                rhash |= nuc_to_2bit[(int)rseq[i]];
            }
        }
    }

    // destructor
    ~prefix_t(){}

};

// prefix_t sorter
struct prefix_sort_t{
    inline bool operator()(const prefix_t* p1, const prefix_t* p2){
        return p1->fseq.length() < p2->fseq.length();
    }
};

// a list of prefix_t
typedef std::vector<prefix_t*> PrefixList;

// a list of same len hash hash 
typedef std::vector<std::unordered_map<uint64_t, PrefixList>> PrefixHash;

// a list of same len string hash
typedef std::vector<std::unordered_map<std::string, PrefixList>> PreStrHash;

inline void file2prl(char* prf, PrefixList& prl){
    // parse prefix cfg
    std::ifstream fr(prf);
    std::string line;
    std::vector<std::string> vstr;
    std::map<std::string, int> smps;
    int count = 0, idx = 0;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        if(vstr.size() < 3){
            fprintf(stderr, "input configure file format wrong, it must be a 3-column TSV format\n");
            fprintf(stderr, "the 3 fields are: name,fwd_primer,rev_primer\n");
            fr.close();
            exit(EXIT_FAILURE);
        }
        util::nuc2upper(vstr[1]);
        util::nuc2upper(vstr[2]);
        auto iter = smps.find(vstr[0]);
        if(iter == smps.end()){
            idx = count;
            smps[vstr[0]] = count++;
        }else{
            idx = iter->second;
        }
        prefix_t* pr = new prefix_t(vstr[1], vstr[2], '+', idx);
        prl.push_back(pr);
    }
    fr.close();
}

#endif
