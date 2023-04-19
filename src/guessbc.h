#ifndef GUESS_BC_H
#define GUESS_BC_H

#include "util.h"
#include "krec.h"
#include "htslib/kstring.h"
#include <unordered_map>

struct bcs1_t{
    std::string s;
    int count;
};

struct bcs_cmp{
    bool operator()(const bcs1_t* bs1, const bcs1_t* bs2) const {
        return bs1->count > bs2->count;
    }
};

struct gbc_opt_t{
    char* inf1 = NULL;
    char* inf2 = NULL;
    char* outfa = NULL;
    char* outtsv = NULL;
    int len1 = 9;
    int len2 = 9;
    int off1 = 0;
    int off2 = 0;
    int mins = 3;
    int64_t maxc = 0;
    bool cat = false;
    int topn = 1536;
    int topr = 1000000;
    std::vector<bcs1_t*> bcs;

    gbc_opt_t();
    ~gbc_opt_t();

    void guess();
    void guess1(std::unordered_map<std::string, int64_t>& ret, char* infq, int len, int off);
    void guess2(std::unordered_map<std::string, int64_t>& ret, char* infq1, char* infq2);
};

void gbc_usage(gbc_opt_t* opt, char* arg0);
int gbc_main(int argc, char** argv);

#endif
