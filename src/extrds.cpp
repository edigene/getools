#include "extrds.h"

void ext_opt_t::ext(){
    samFile* ifp = sam_open(inbam, "r");
    bam_hdr_t* h = sam_hdr_read(ifp);
    hts_idx_t* idx = sam_index_load(ifp, inbam);
    hts_itr_t* itr = NULL;
    if(contig) itr = sam_itr_querys(idx, h, contig);
    else itr = sam_itr_querys(idx, h, ".");
    bam1_t* b = bam_init1();
    std::vector<bam1_t*> ubs;
    uint8_t* data = NULL;
    uint8_t* bcd = NULL;
    char* bcv = 0;
    bool mbc = false;
    int64_t val = 0;
    int64_t ttr = 0;
    while(sam_itr_next(ifp, itr, b) >= 0){
        data = bam_aux_get(b, "VS");
        if(data){
            if(cellbc){
                mbc = false;
                bcd = bam_aux_get(b, CELL_BARCODE_ID_TAG);
                if(bcd){
                    bcv = bam_aux2Z(bcd);
                    char* pcbv = strdup(bcv);
                    int smax = 0, *sofs = 0;
                    int snf = ksplit_core(pcbv, '|', &smax, &sofs);
                    for(int sni = 0; sni < snf; ++sni){
                        char* psstr = pcbv+sofs[sni];
                        char* pcstr = strchr(psstr, '-');
                        *pcstr = '\0';
                        if(strcmp(psstr, cellbc) == 0){
                            mbc = true;
                            break;
                        }
                    }
                    free(pcbv);
                    free(sofs);
                }
                if(!mbc) continue;
            }
            val = bam_aux2i(data);
            if(((val & imask) == imask) && (!(val & emask))){
                if(mergs < 0 || (mergs >= 0 && get_merge_status(val) == mergs)){
                    ttr += bam_aux2i(bam_aux_get(b, "CC"));
                    if(outbam){
                        ubs.push_back(b);
                        b = bam_init1();
                    }
                }
            }
        }
    }
    sam_close(ifp);
    bam_destroy1(b);
    hts_idx_destroy(idx);
    hts_itr_destroy(itr);
    if(outbam){
        if(sbyc) std::sort(ubs.begin(), ubs.end(), BamSortByCoord());
        else std::sort(ubs.begin(), ubs.end(), BamSortByQname());
        samFile* ofp = sam_open(outbam, "wb");
        assert(sam_hdr_write(ofp, h) >= 0);
        for(auto& e: ubs){
            assert(sam_write1(ofp, h, e) >= 0);
            bam_destroy1(e); e = NULL;
        }
        sam_close(ofp);
        if(sbyc) assert(sam_index_build(outbam, 0) >= 0);
    }
    fprintf(stderr, "Total CC(original reads count) got: %lld\n", ttr);
}

void extrds_usage(ext_opt_t* opt, char* arg0){
     fprintf(stderr, "\n");
     fprintf(stderr, "Usage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
     fprintf(stderr, "Options: -i FILE input bam generated by getools\n");
     fprintf(stderr, "         -o FILE output bam to store extract records\n");
     fprintf(stderr, "         -t STR  contig name to extract reads\n");
     fprintf(stderr, "         -b STR  cell barcode to extract reads\n");
     fprintf(stderr, "         -m INT  masks sequences should met [%d]\n", opt->imask);
     fprintf(stderr, "         -r INT  masks sequences should not met [%d]\n", opt->emask);
     fprintf(stderr, "         -s INT  merge status sequences should met [%d]\n", opt->mergs);
     fprintf(stderr, "         -c      sort output by coordinates\n");
     fprintf(stderr,  "\n");
}

int extrds_main(int argc, char** argv){
    ext_opt_t opt;
    int c = 0;
    while((c = getopt(argc, argv, "i:o:t:b:m:r:s:ch")) >= 0){
        switch(c){
            case 'i': opt.inbam = optarg; break;
            case 'o': opt.outbam = optarg; break;
            case 't': opt.contig = optarg; break;
            case 'b': opt.cellbc = optarg; break;
            case 'm': opt.imask = atoll(optarg); break;
            case 'r': opt.emask = atoll(optarg); break;
            case 's': opt.mergs = atoi(optarg); break;
            case 'c': opt.sbyc = true; break;
            case 'h': extrds_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(!opt.valid()){
        extrds_usage(&opt, argv[0]);
        return 0;
    }
    opt.ext();
    return 0;
}