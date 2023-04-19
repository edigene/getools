#include "b2g.h"

b2g_opt_t::b2g_opt_t(){
}

b2g_opt_t::~b2g_opt_t(){
    if(fai){
        fai_destroy(fai);
        fai = NULL;
    }
}

bool b2g_opt_t::valid(){
    if(!inbcf){
        fprintf(stderr, "input bcf file must be provided\n");
        return false;
    }
    if(!infa){
        fprintf(stderr, "input ref fasta file must be provided\n");
        return false;
    }
    if(outdir.empty()){
        fprintf(stderr, "output directory must be provided\n");
        return false;
    }
    return true;
}

void b2g_opt_t::init(){
    fai = fai_load(infa);
    if(!util::exists(outdir)) util::makedir(outdir);
    outfa = util::joinpath(outdir, outfa);
    outtsv = util::joinpath(outdir, outtsv);
}

void b2g_opt_t::gtn2fa(){
    std::ofstream fwfa(outfa);
    std::ofstream fwtsv(outtsv);
    fwtsv << "Contig\tName\tFreq\tSequence\n";
    kstring_t* tmps = (kstring_t*)calloc(1, sizeof(kstring_t));
    kstring_t* tmpn = (kstring_t*)calloc(1, sizeof(kstring_t));
    htsFile* fp = bcf_open(inbcf, "r");
    hts_idx_t* idx = bcf_index_load(inbcf);
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    bcf1_t* b = bcf_init1();
    for(int i = 0; i < faidx_nseq(fai); ++i){
        const char* amp = faidx_iseq(fai, i);
        int rlen = 0;
        char* fas = fai_fetch(fai, amp, &rlen);
        fwfa << ">" << amp << "\n" << fas << "\n";
        fwtsv << amp << "\t" << amp << "_ref" << "\t0\t" << fas << "\n";
        int tid = bcf_hdr_name2id(hdr, amp);
        int cid = bcf_hdr_id2int(hdr, BCF_DT_CTG, amp);
        bcf_hrec_t* hrec = bcf_hdr_id2hrec(hdr, BCF_DT_CTG, BCF_HL_CTG, cid);
        int hidx = bcf_hrec_find_key(hrec, "geb");
        int beg = atoi(hrec->vals[hidx]);
        hidx = bcf_hrec_find_key(hrec, "gee");
        int end = atoi(hrec->vals[hidx]);
        hts_itr_t* itr = bcf_itr_queryi(idx, tid, beg, end);
        int* vtv = NULL; float* afv = NULL; int* adv = NULL;
        int nvt = 0, naf = 0, nad = 0;
        int nsmp = bcf_hdr_nsamples(hdr);
        // get all vars
        std::vector<b2g_var_t*> vars;
        while(bcf_itr_next(fp, itr, b) >= 0){
            bcf_unpack(b, BCF_UN_INFO);
            if(bcf_get_info_int32(hdr, b, "VT", &vtv, &nvt) > 0){
                if(vtv[0] & type){
                    bcf_unpack(b, BCF_UN_ALL);
                    bcf_get_format_float(hdr, b, "AF", &afv, &naf);
                    bcf_get_format_int32(hdr, b, "AD", &adv, &nad);
                    for(int a = 0; a < b->n_allele-1; ++a){
                        if(nsmp == 1 ||
                           (adv[a] > 0 && adv[nsmp*b->n_allele+a] < nctc)){
                            b2g_var_t* v = new b2g_var_t();
                            v->ref = b->d.allele[0];
                            v->alt = b->d.allele[a+1];
                            v->freq = afv[a];
                            v->start = b->pos;
                            if(vtv[0] & (GEVAR_DEL | GEVAR_DIN)){
                                v->end = v->start + v->ref.size()-1;
                            }else if(vtv[0] & (GEVAR_SNV | GEVAR_INS)){
                                v->end = v->start;
                            }
                            v->width = (int)v->alt.size()-(int)v->ref.size();
                            v->type = GEV_STR_ARR4_B2G[vtv[0]];
                            vars.push_back(v);
                        }
                    }
                }
            }
        }
        if(afv){ free(afv); afv = NULL; }
        if(adv){ free(adv); adv = NULL; }
        if(vtv){ free(vtv); vtv = NULL; }
        hts_itr_destroy(itr);
        // sort vars by freq
        std::sort(vars.begin(), vars.end(), sort8freq());
        // output topN vars to fa
        size_t maxn = topn;
        if(maxn > vars.size()) maxn = vars.size();
        for(size_t i = 0; i < maxn; ++i){
            tmps->l = 0;
            kputsn(fas, vars[i]->start, tmps);
            kputs(vars[i]->alt.c_str(), tmps);
            kputs(fas+vars[i]->end+1, tmps);
            tmpn->l = 0;
            ksprintf(tmpn, "%s_t%ld_p%d_%s%d_%s_%s",
                           amp,
                           i+1,
                           vars[i]->start,
                           vars[i]->type.c_str(),
                           abs(vars[i]->width),
                           vars[i]->ref.c_str(),
                           vars[i]->alt.c_str());
            fwfa << ">" << tmpn->s << "\n" << tmps->s << "\n";
            fwtsv << amp << "\t" << tmpn->s << "\t" << vars[i]->freq << "\t" << tmps->s << "\n";
        }
        // clean up vars
        for(size_t i = 0; i < vars.size(); ++i){
            delete vars[i];
            vars[i] = NULL;
        }
        vars.clear();
        free(fas); fas = NULL;
    }
    // close bcf
    bcf_close(fp);
    hts_idx_destroy(idx);
    bcf_hdr_destroy(hdr);
    bcf_destroy(b);
    // close file
    fwfa.close();
    fwtsv.close();
    if(tmpn->s){ free(tmpn->s); free(tmpn); }
    if(tmps->s){ free(tmps->s); free(tmps); }
}

void b2g_usage(b2g_opt_t* opt, char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options:  -i FILE input bcf generated by getools(var.bcf)\n");
    fprintf(stderr, "          -r FILE input fasta ref generated by getools(ref.fa)\n");
    fprintf(stderr, "          -o STR  output directory path [%s]\n", opt->outdir.c_str());
    fprintf(stderr, "          -n INT  number of top variants to output [%d]\n", opt->topn);
    fprintf(stderr, "          -v INT  masks of variant type to output(1:SNV,2:INSERTION,4:DELETION,8:DELINS) [%d]\n", opt->type);
    fprintf(stderr, "\n");
}

int b2g_main(int argc, char** argv){
    b2g_opt_t opt;
    if(argc == 1){
        b2g_usage(&opt, argv[0]);
        return 0;
    }
    int c = 0;
    while((c = getopt(argc, argv, "i:r:o:n:v:h")) >= 0){
        switch(c){
            case 'i': opt.inbcf = optarg; break;
            case 'r': opt.infa = optarg; break;
            case 'o': opt.outdir = optarg; break;
            case 'n': opt.topn = atoi(optarg); break;
            case 'v': opt.type = atoi(optarg); break;
            case 'h': b2g_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.valid()){
        opt.init();
        opt.gtn2fa();
        return 0;
    }else{
        return 1;
    }
}
