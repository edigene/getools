#include "stats.h"

void stat_opt_t::stat(){
    // alloc dist store 
    // vartypecnt[i] = number of sequences with total i vartype
    std::vector<int64_t> snpcnt(buflen, 0); // snp
    std::vector<int64_t> inscnt(buflen, 0); // ins
    std::vector<int64_t> delcnt(buflen, 0); // del
    std::vector<int64_t> sclcnt(buflen, 0); // clip
    std::vector<int64_t> varcnt(buflen, 0); // snp+ins+del(total var)
    // vartypelen[i] = number of sequences with vartype length of i 
    std::vector<int64_t> snplen(buflen, 0); // snp
    std::vector<int64_t> inslen(buflen, 0); // ins
    std::vector<int64_t> dellen(buflen, 0); // del
    std::vector<int64_t> scllen(buflen, 0); // clip
    // vartype*pos[i] = number of sequences with vartype var at position i on ref or qry
    std::vector<int64_t> snprpos(buflen, 0); // snp on ref
    std::vector<int64_t> insrpos(buflen, 0); // ins on ref
    std::vector<int64_t> delrpos(buflen, 0); // del on ref
    std::vector<int64_t> sclrpos(buflen, 0); // scl on ref
    std::vector<int64_t> snpqpos(buflen, 0); // snp on qry
    std::vector<int64_t> insqpos(buflen, 0); // ins on qry
    std::vector<int64_t> delqpos(buflen, 0); // del on qry
    std::vector<int64_t> sclqpos(buflen, 0); // scl on qry
    // linkage stat vartype*lcnt[i][j] = number of sequences with j vartype var in [0,i] on ref or qry
    std::vector<std::vector<int64_t>> snvrlcnt(buflen, std::vector<int64_t>(buflen, 0)); // linkage snv on ref
    std::vector<std::vector<int64_t>> insrlcnt(buflen, std::vector<int64_t>(buflen, 0)); // linkage ins on ref
    std::vector<std::vector<int64_t>> delrlcnt(buflen, std::vector<int64_t>(buflen, 0)); // linkage del on ref
    std::vector<std::vector<int64_t>> snvqlcnt(buflen, std::vector<int64_t>(buflen, 0)); // linkage snv on qry
    std::vector<std::vector<int64_t>> insqlcnt(buflen, std::vector<int64_t>(buflen, 0)); // linkage ins on qry
    std::vector<std::vector<int64_t>> delqlcnt(buflen, std::vector<int64_t>(buflen, 0)); // linkage del on qry
    // do stats
    samFile* ifp = sam_open(inbam, "r");
    sam_hdr_t* hdr = sam_hdr_read(ifp);
    hts_idx_t* idx = sam_index_load(ifp, inbam);
    hts_itr_t* itr = NULL;
    if(name) itr = sam_itr_querys(idx, hdr, name);
    else itr = sam_itr_querys(idx, hdr, ".");
    bam1_t* b = bam_init1();
    uint64_t ttr = 0;
    uint8_t* vsd = NULL;
    uint8_t* ccd = NULL;
    int64_t vsv = 0, ccv = 0;
    uint32_t* cigar = NULL;
    int oplen = 0, opint = 0;
    int32_t rpos = 0, qpos = 0;
    int ic = 0, sc = 0, dc = 0, zc = 0;
    int lsnv = 0, lins = 0, lind = 0, msnv = 0, mins = 0, mdel = 0;
    while(sam_itr_next(ifp, itr, b) >= 0){
        if(b->core.flag & BAM_FUNMAP) continue;
        vsd = bam_aux_get(b, "VS");
        vsv = bam_aux2i(vsd);
        if(((vsv & imask) == imask) && (!(vsv & emask))){
            ic = sc = dc = zc = 0;
            lsnv = lins = lind = 0;
            rpos = b->core.pos;
            qpos = 0;
            ccd = bam_aux_get(b, "CC");
            ccv = bam_aux2i(ccd);
            ttr += ccv;
            cigar = bam_get_cigar(b);
            for(uint32_t c = 0; c < b->core.n_cigar; ++c){
                opint = bam_cigar_op(cigar[c]);
                oplen = bam_cigar_oplen(cigar[c]);
                switch(opint){
                    case BAM_CINS:
                        inslen[oplen] += ccv;
                        insrpos[rpos] += ccv;
                        insqpos[qpos] += ccv;
                        ++lins;
                        insqlcnt[qpos][lins] += ccv;
                        insrlcnt[rpos][lins] += ccv;
                        ++ic;
                        qpos += oplen;
                        break;
                    case BAM_CDEL:
                        dellen[oplen] += ccv;
                        delrpos[rpos] += ccv;
                        delqpos[qpos] += ccv;
                        ++lind;
                        delqlcnt[qpos][lind] += ccv;
                        delrlcnt[rpos][lind] += ccv;
                        ++dc;
                        rpos += oplen;
                        break;
                    case BAM_CDIFF:
                        snplen[oplen] += ccv;
                        snprpos[rpos] += ccv;
                        snpqpos[qpos] += ccv;
                        ++lsnv;
                        snvqlcnt[qpos][lsnv] += ccv;
                        snvrlcnt[rpos][lsnv] += ccv;
                        ++sc;
                        rpos += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CEQUAL:
                        rpos += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CSOFT_CLIP:
                        scllen[oplen] += ccv;
                        sclrpos[rpos] += ccv;
                        sclqpos[qpos] += ccv;
                        ++zc;
                        qpos += oplen;
                        break;
                    default: break;
                }
            }
            snpcnt[sc] += ccv;
            inscnt[ic] += ccv;
            delcnt[dc] += ccv;
            sclcnt[zc] += ccv;
            varcnt[ic+dc+sc] += ccv;
            maxl = MAX(qpos, maxl);
            maxl = MAX(rpos, maxl);
            msnv = MAX(msnv, sc);
            mins = MAX(mins, ic);
            mdel = MAX(mdel, dc);
        }
    }
    sam_close(ifp);
    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    // output
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    out_head(ks);
    out_rec1d(ks, "snpcnt", snpcnt);
    out_rec1d(ks, "inscnt", inscnt);
    out_rec1d(ks, "delcnt", delcnt);
    out_rec1d(ks, "sclcnt", sclcnt);
    out_rec1d(ks, "varcnt", varcnt);
    out_rec1d(ks, "snplen", snplen);
    out_rec1d(ks, "inslen", inslen);
    out_rec1d(ks, "dellen", dellen);
    out_rec1d(ks, "scllen", scllen);
    out_rec1d(ks, "snprpos", snprpos);
    out_rec1d(ks, "insrpos", insrpos);
    out_rec1d(ks, "delrpos", delrpos);
    out_rec1d(ks, "sclrpos", sclrpos);
    out_rec1d(ks, "snpqpos", snpqpos);
    out_rec1d(ks, "insqpos", insqpos);
    out_rec1d(ks, "delqpos", delqpos);
    out_rec1d(ks, "sclqpos", sclqpos);
    out_rec2d(ks, "linkrsnv", msnv, snvrlcnt);
    out_rec2d(ks, "linkrdel", mdel, delrlcnt);
    out_rec2d(ks, "linkrins", mins, insrlcnt);
    out_rec2d(ks, "linkqsnv", msnv, snvqlcnt);
    out_rec2d(ks, "linkqdel", mdel, delqlcnt);
    out_rec2d(ks, "linkqins", mins, insqlcnt);
    ksprintf(ks, "#totr\t%lld\n", ttr);
    FILE* ofp = fopen(outtsv, "w");
    fwrite(ks->s, sizeof(char), ks->l, ofp);
    fclose(ofp);
}

void stat_opt_t::out_head(kstring_t* ks){
    ksprintf(ks, "Index\t");
    for(int i = 0; i < maxl; ++i) ksprintf(ks, "%d\t", i);
    ks->s[ks->l-1] = '\n';
}

void stat_opt_t::out_rec1d(kstring_t* ks, const char*name, const std::vector<int64_t>& res){
    ksprintf(ks, "%s\t", name);
    for(int i = 0; i < maxl; ++i) ksprintf(ks, "%lld\t", res[i]);
    ks->s[ks->l-1] = '\n';
}

void stat_opt_t::out_rec2d(kstring_t* ks, const char* name, int mh, const std::vector<std::vector<int64_t>>& res){
    ksprintf(ks, "%s\t", name);
    for(int i = 0; i < maxl; ++i) ksprintf(ks, "%d\t", i);
    ks->s[ks->l-1] = '\n';
    for(int i = 1; i <= mh; ++i){
        ksprintf(ks, "___l%d\t", i);
        for(int j = 0; j < maxl; ++j){
            ksprintf(ks, "%lld\t", res[j][i]);
        }
        ks->s[ks->l-1] = '\n';
    }
}

void stats_usage(stat_opt_t* opt, char* arg0){
    fprintf(stderr, "\nUsage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options:  -i  input BAM file\n");
    fprintf(stderr, "          -l  max(rlen, qlen)[%d]\n", opt->rlen);
    fprintf(stderr, "          -n  contig name to stat\n");
    fprintf(stderr, "          -m  VS flag must be met[%d]\n", opt->imask);
    fprintf(stderr, "          -e  VS flag must not met[%d]\n", opt->emask);
    fprintf(stderr, "          -o  output tsv file\n");
    fprintf(stderr, "\n");
}

int stats_main(int argc, char** argv){
    stat_opt_t opt;
    if(argc == 1){
        stats_usage(&opt, argv[0]);
        return 0;
    }
    int c = 0;
    while((c = getopt(argc, argv, "i:l:n:m:e:o:h")) >= 0){
        switch(c){
            case 'i': opt.inbam = optarg; break;
            case 'l': opt.rlen = atoi(optarg); break;
            case 'n': opt.name = optarg; break;
            case 'm': opt.imask = atoi(optarg); break;
            case 'e': opt.emask = atoi(optarg); break;
            case 'o': opt.outtsv = optarg; break;
            case 'h': stats_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    opt.init();
    opt.stat();
    return 0;
}
