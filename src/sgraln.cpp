#include "sgraln.h"

sgr_aln_opt::sgr_aln_opt(){
}

sgr_aln_opt::~sgr_aln_opt(){
    if(mat){ free(mat); mat = NULL; }
    for(auto& e: sgrseqs){
        delete e;
        e = NULL;
    }
}

bool sgr_aln_opt::valid(){
    if(sgrfwd == NULL || inref == NULL){
        fprintf(stderr, "sgRNA_PAM sequence and reference file must be provided\n");
        return false;
    }
    if(outref == NULL && outaln == NULL && outtsv == NULL){
        fprintf(stderr, "output file must be provided at least one kine\n");
        return false;
    }
    if(match < 0 || mismatch < 0 || gape < 0 || gapo < 0){
        fprintf(stderr, "all score options must be positive integers\n");
        return false;
    }
    if(flk < 0){
        fprintf(stderr, "flank length must be positive integers\n");
        return false;
    }
    if(getbp){
        if(breakp < 0){
            fprintf(stderr, "breakpoint must be positive integer\n");
            return false;
        }
    }
    return true;
}

void sgr_aln_opt::init(){
    std::vector<std::string> vsgrs;
    util::split(sgrfwd, vsgrs, ",");
    for(auto& s : vsgrs){
        util::nuc2upper(s);
        sgr_seq_t* sgrs = new sgr_seq_t(s.c_str());
        sgrseqs.push_back(sgrs);
    }
    mat = kswge_gen_smat(match, mismatch);
    // N treat as match
    int i = 0, k = 0, j = 0;
    for(; i < 4; ++i){
        for(j = 0; j < 4; ++j) ++k;
        mat[k++] = match;
    }
    for(j = 0; j < 5; ++j) mat[k++] = match;
}

void sgr_aln_opt::search(){
    gzFile ifp = gzopen(inref, "r");
    kseq1_t* seq = kseq1_init(ifp);
    krec1_t* rec = krec1_init();
    char* ref = NULL;
    int rlen = 0;
    uint8_t* rints = NULL;
    kswr_t* best = NULL;
    sgr_seq_t* sgr = NULL;
    bool isrev = false;
    int qbeg = 0, qend = 0, rbeg = 0, rend = 0;
    int roff = 0, realb = 0, reale = 0;
    if(getbp) roff = MAX(breakp - flk, 0);
    kstring_t ofas = {0, 0, 0};
    kstring_t otbs = {0, 0, 0};
    FILE* ofaln = NULL;
    if(outaln) ofaln = fopen(outaln, "w");
    while(kseq1_read(seq, rec) >= 0){ 
        for(size_t i = 0; i < rec->seq.l; ++i) rec->seq.s[i] = nuc_to_lower[(int)rec->seq.s[i]];
        if(!getbp){
            breakp = 0.5 * rec->seq.l;
            roff = MAX(breakp - flk, 0);
        }
        ref = rec->seq.s + roff; rlen = MIN(2 * flk, (int)rec->seq.l - roff);
        rints = kswge_seq2ints(ref, rlen);
        for(auto& sgrs: sgrseqs){
            sgrs->faln = kswge_semi_global(sgrs->sgrlen, sgrs->sgrintf, rlen, rints, 5, mat, gapo, gape);
            sgrs->raln = kswge_semi_global(sgrs->sgrlen, sgrs->sgrintr, rlen, rints, 5, mat, gapo, gape);
            sgrs->best = sgrs->faln; sgrs->maxisr = false;
            if(sgrs->raln->score > sgrs->faln->score){
                sgrs->best = sgrs->raln;
                sgrs->maxisr = true;
            }
            if(best == NULL || best->score < sgrs->best->score){
                best = sgrs->best;
                isrev = sgrs->maxisr;
                sgr = sgrs;
            }
        }
        if(best->score < minscore){
            fprintf(stderr, "%s not found, score too low %d\n", rec->name.s, best->score);
            for(auto& e : sgrseqs) e->clear();
            best = NULL;
            sgr = NULL;
            continue;
        }
        qbeg = best->qb;
        qend = best->qe;
        rbeg = best->tb;
        rend = best->te;
        realb = roff + rbeg;
        reale = roff + rend;
        if(qbeg > 0) rbeg = MAX(0, rbeg - qbeg);
        if(qend < sgr->sgrlen - 1) rend = MIN(rlen-1, rend+(sgr->sgrlen-1-qend));
        // out aln
        if(ofaln){
            if(isrev) kswge_ret_output(ofaln, best, ref, sgr->sgrrev, kswge_nt2int, rec->name.s, 60);
            else kswge_ret_output(ofaln, best, ref, sgr->sgrfwd, kswge_nt2int, rec->name.s, 60);
        }
        // process fa
        for(int i = rbeg; i <= rend; ++i) ref[i] = nuc_to_upper[(int)ref[i]];
        // out fa
        ksprintf(&ofas, ">%s\n%s\n", rec->name.s, rec->seq.s);
        // out tsv
        ksprintf(&otbs, "%s\t%s\t", rec->name.s, rec->seq.s);
        if(isrev){
            for(int i = rbeg; i <= rend; ++i) kputc(nuc_to_cmp[(int)ref[rend-(i-rbeg)]], &otbs);
        }else{
            for(int i = rbeg; i <= rend; ++i) kputc(ref[i], &otbs);
        }
        ksprintf(&otbs, "\t%c\t%d\t%d\t%d\n", isrev ? '-' : '+', best->score, realb, reale);
        // clear
        for(auto& e: sgrseqs) e->clear();
        best = NULL;
        sgr = NULL;
    }
    if(ofaln) fclose(ofaln);
    if(outtsv){
        FILE* ofptsv = fopen(outtsv, "w");
        fwrite(otbs.s, sizeof(char), otbs.l, ofptsv);
        fclose(ofptsv);
    }
    if(outref){
        FILE* ofpfa = fopen(outref, "w");
        fwrite(ofas.s, sizeof(char), ofas.l, ofpfa);
        fclose(ofpfa);
    }
    if(otbs.s) free(otbs.s);
    if(ofas.s) free(ofas.s);
    kseq1_destroy(seq);
    krec1_destroy(rec);
    gzclose(ifp);
}

void sgraln_usage(sgr_aln_opt* opt, char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options: -i FILE input fa file to detect sgRNA\n");
    fprintf(stderr, "         -s FILE sgRNA_PAM sequences, seperated by comma\n");
    fprintf(stderr, "         -b INT  breakpoint posision\n");
    fprintf(stderr, "         -f INT  flank length around breakpoint to look for sgRNA\n");
    fprintf(stderr, "         -o FILE output fa format results\n");
    fprintf(stderr, "         -t FILe output tsv format results\n");
    fprintf(stderr, "         -a FILE output align format results\n");
    fprintf(stderr, "         -M INT  match score [%d]\n", opt->match);
    fprintf(stderr, "         -X INT  mismatch penalty [%d]\n", opt->mismatch);
    fprintf(stderr, "         -O INT  gap open penalty [%d]\n", opt->gapo);
    fprintf(stderr, "         -E INT  gap extend penalty [%d]\n", opt->gape);
    fprintf(stderr, "         -S INT  min match score [%d]\n", opt->minscore);
    fprintf(stderr, "\n");
}

int sgraln_main(int argc, char** argv){
    sgr_aln_opt opt;
    if(argc == 1){
        sgraln_usage(&opt, argv[0]);
        return 0;
    }
    int c = 0;
    while((c = getopt(argc, argv, "i:s:b:f:o:t:a:M:X:O:E:S:h")) >= 0){
        switch(c){
            case 'i': opt.inref = optarg; break;
            case 's': opt.sgrfwd = optarg; break;
            case 'b': opt.breakp = atoi(optarg); opt.getbp = true; break;
            case 'f': opt.flk = atoi(optarg); break;
            case 'o': opt.outref = optarg; break;
            case 't': opt.outtsv = optarg; break;
            case 'a': opt.outaln = optarg; break;
            case 'M': opt.match = atoi(optarg); break;
            case 'X': opt.mismatch = atoi(optarg); break;
            case 'O': opt.gapo = atoi(optarg); break;
            case 'E': opt.gape = atoi(optarg); break;
            case 'S': opt.minscore = atoi(optarg); break;
            case 'h': sgraln_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.valid()){
        opt.init();
        opt.search();
        return 0;
    }
    return 1;
}

#ifdef SGRALN_MAIN_FUN
int main(int argc, char** argv){
    return sgraln_main(argc, argv);
}
#endif
