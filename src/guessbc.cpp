#include "guessbc.h"

gbc_opt_t::gbc_opt_t(){
}

gbc_opt_t::~gbc_opt_t(){
    for(auto& b: bcs){
        if(b){ free(b); b = NULL; }
    }
}

void gbc_opt_t::guess2(std::unordered_map<std::string, int64_t>& ret, char* infq1, char* infq2){
    if(!infq1 || !infq2) return;
    kstring_t key = {0, 0, 0};
    gzFile ifp1 = gzopen(infq1, "r");
    gzFile ifp2 = gzopen(infq2, "r");
    kseq1_t* ks1 = kseq1_init(ifp1);
    kseq1_t* ks2 = kseq1_init(ifp2);
    krec1_t* kr1 = krec1_init();
    krec1_t* kr2 = krec1_init();
    int ttr = 0;
    while(kseq1_read(ks1, kr1) >= 0 && kseq1_read(ks2, kr2) >= 0){
        if(kr1->seq.l < off1 + len1 || kr2->seq.l < off2 + len2) continue;
        ++ttr;
        key.l = 0;
        kputsn(kr1->seq.s+off1, len1, &key);
        kputsn(kr2->seq.s+off2, len2, &key);
        auto iter = ret.find(key.s);
        if(iter == ret.end()){
            if(maxc == 0) maxc = 1;
            ret[key.s] = 1;
        }else{
            ++iter->second;
            if(maxc < iter->second) maxc = iter->second;
        }
        if(ttr > topr) break;
    }
    kseq1_destroy(ks1);
    kseq1_destroy(ks2);
    krec1_destroy(kr1);
    krec1_destroy(kr2);
    gzclose(ifp1);
    gzclose(ifp2);
}

void gbc_opt_t::guess1(std::unordered_map<std::string, int64_t>& ret, char* infq, int len, int off){
    if(!infq) return;
    gzFile ifp = gzopen(infq, "r");
    kseq1_t* ks = kseq1_init(ifp);
    krec1_t* kr = krec1_init();
    int ttr = 0;
    while(kseq1_read(ks, kr) >= 0){
        if(kr->seq.l < off + len) continue;
        ++ttr;
        std::string key(kr->seq.s + off, len);
        auto iter = ret.find(key);
        if(iter == ret.end()){
            if(maxc == 0) maxc = 1;
            ret[key] = 1;
        }else{
            ++iter->second;
            if(maxc < iter->second) maxc = iter->second;
        }
        if(ttr > topr) break;
    }
    kseq1_destroy(ks);
    krec1_destroy(kr);
    gzclose(ifp);
}

void gbc_opt_t::guess(){
    std::unordered_map<std::string, int64_t> bcm;
    if(cat){
        guess2(bcm, inf1, inf2);
    }else{
        if(inf1) guess1(bcm, inf1, len1, off1);
        if(inf2) guess1(bcm, inf2, len2, off2);
    }
    if(bcm.empty()) return;
    std::vector<int64_t> cal(maxc+1, 0);
    for(auto& iter: bcm) ++cal[iter.second];
    int acc = cal[maxc];
    int i = maxc;
    while(acc < topn && i >= 1){
        acc += cal[i-1];
        --i;
    }
    int mi = MAX(i, mins);
    for(auto& iter: bcm){
        if(iter.second >= mi){
            bcs1_t* bc = (bcs1_t*)calloc(1, sizeof(bcs1_t));
            bc->count = iter.second;
            bc->s = iter.first;
            bcs.push_back(bc);
        }
    }
    std::sort(bcs.begin(), bcs.end(), bcs_cmp());
    fprintf(stderr, "total valid candidates type: %ld\n", bcs.size());
    int maxn = bcs.size();
    if(maxn < topn) topn = maxn;
    if(outfa){
        std::ofstream fw(outfa);
        for(int i = 0; i < topn; ++i){
            fw << ">seq" << i+1 << " " << bcs[i]->count << "\n";
            fw << bcs[i]->s << "\n";
        }
        fw.close();
    }
    if(outtsv){
        std::ofstream fw(outtsv);
        fw << "topN\tcount\tseq\n";
        for(int i = 0; i < topn; ++i){
            fw << i+1 << "\t" << bcs[i]->count << "\t" << bcs[i]->s << "\n";
        }
        fw.close();
    }
}

void gbc_usage(gbc_opt_t* opt, char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s [options]\n\n", arg0);
    fprintf(stderr, "Options: -i FILE input fa/q1\n");
    fprintf(stderr, "         -I FILE input fa/q2\n");
    fprintf(stderr, "         -l INT  prefix length r1 [%d]\n", opt->len1);
    fprintf(stderr, "         -L INT  prefix length r2 [%d]\n", opt->len2);
    fprintf(stderr, "         -b INT  beg pos of prefix r1 [%d]\n", opt->off1);
    fprintf(stderr, "         -B int  beg pos of prefix r2 [%d]\n", opt->off2);
    fprintf(stderr, "         -n INT  topN to get [%d]\n", opt->topn);
    fprintf(stderr, "         -t INT  total reads to use [%d]\n", opt->topr);
    fprintf(stderr, "         -s INT  min support [%d]\n", opt->mins);
    fprintf(stderr, "         -o FILE output fa\n");
    fprintf(stderr, "         -O FILE output tsv\n");
    fprintf(stderr, "         -c      cat p1/p2 if set\n");
    fprintf(stderr, "\n");
}

int gbc_main(int argc, char** argv){
    gbc_opt_t opt;
    if(argc == 1){
        gbc_usage(&opt, argv[0]);
        return 0;
    }

    int c = -1;
    while((c = getopt(argc, argv, "i:I:l:L:n:o:O:b:B:s:t:ch")) >= 0){
        switch(c){
            case 'i': opt.inf1 = optarg; break;
            case 'I': opt.inf2 = optarg; break;
            case 'l': opt.len1 = atoi(optarg); break;
            case 'L': opt.len2 = atoi(optarg); break;
            case 'b': opt.off1 = atoi(optarg); break;
            case 'B': opt.off2 = atoi(optarg); break;
            case 'n': opt.topn = atoi(optarg); break;
            case 't': opt.topr = atoi(optarg); break;
            case 'o': opt.outfa = optarg; break;
            case 'O': opt.outtsv = optarg; break;
            case 's': opt.mins = atoi(optarg); break;
            case 'c': opt.cat = true; break;
            case 'h': gbc_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.inf1 == NULL){
        fprintf(stderr, "must provede at least one fq/a file\n");
        return 1;
    }
    opt.guess();
    return 0;
}


#ifdef GBC_MAIN_FUN
int main(int argc, char** argv){
    return gbc_main(argc, argv);
}
#endif
